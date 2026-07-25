#!/bin/bash
# ==============================================================================
# Integrated RNA-Seq Analysis Pipeline (v1.0.0)
# ==============================================================================
# Description:
#   A master script to run a complete RNA-Seq analysis pipeline.
#   - Execution flow is strictly controlled by flags.
#   - Default flow: FastQ setup -> FastQC -> Salmon -> EdgeR.
#   - Automatically supports Paired/Unpaired designs via CSV structure.
#   - Includes --start-at flag to resume pipeline from specific steps (e.g., edger).
#
# Usage:
#   rnsequant.sh <csv_file|directory> <species> <gencode_release> [OPTIONS]
# ==============================================================================

set -o pipefail

# --- Script Path Definitions ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
SUBSCRIPT_DIR="${SCRIPT_DIR}/rnsequant"

DOWNLOAD_SCRIPT="${SUBSCRIPT_DIR}/run_fastq_download.sh"
FASTQC_SCRIPT="${SUBSCRIPT_DIR}/run_fastqc.sh"
SALMON_SCRIPT="${SUBSCRIPT_DIR}/run_salmon.sh"
EDGER_SCRIPT="${SUBSCRIPT_DIR}/run_edgeR_core.sh"
DRIMSEQ_SCRIPT="${SUBSCRIPT_DIR}/drimseq.sh"
QAPA_SCRIPT="${SUBSCRIPT_DIR}/qapa.sh"

# --- Defaults ---
RETRY_COUNT=0
OUTPUT_BASE_DIR=""
EXEC_MODE="docker"
TRANSCRIPT_TYPE_FLAG=""
START_AT="download"

# Flow Control Flags
RUN_DRIMSEQ=false
RUN_QAPA=false
LOCAL_FASTQ_DIR=""

# Component args
DOWNLOAD_ARGS=""
FASTQC_ARGS=""
SALMON_ARGS=""
EDGER_ARGS=""
DRIMSEQ_ARGS=""
QAPA_ARGS=""
DRIMSEQ_QAPA_ARGS=""
CONTEXT=""

# --- Help Message ---
usage() {
  cat <<EOM
Usage: $(basename "$0") <csv_file|directory> <species> <gencode_release> [OPTIONS]

Description:
  Run a complete RNA-Seq analysis pipeline.
  
  [Default Flow]
  1. Download FastQ (or prepare local files).
  2. FastQC -> Salmon -> EdgeR.
  3. Stops after EdgeR.

REQUIRED ARGUMENTS:
  <csv_file|directory>    Path to a single .csv file OR a directory containing .csv files.
  <species>               "human" or "mouse".
  <gencode_release>       GENCODE release version (e.g., 46).

OPTIONS:
  -o, --output-dir <dir>  Base directory for outputs.
  --retry [num]           Enable retry on failure (default 2 if num omitted).
  
  --start-at <step>       Start pipeline from a specific step.
                          Valid: download, fastqc, salmon, edger, drimseq, qapa
                          (Useful to skip FastQ/Salmon steps if outputs already exist).
  
  --local-fastq <dir>     Use local FastQ files from <dir> instead of downloading.
                          - Skips download step.
  
  --run-drimseq           Execute DRIMSeq (DTU analysis) after EdgeR.
  --run-qapa              Execute QAPA (APA analysis) after EdgeR.
                          (Includes post-QAPA DRIMSeq step).
                          
  --local                 Use locally installed tools instead of Docker/Singularity.
  
  Transcript Type Options:
  (default)               Protein-coding transcripts (GENCODE)
  --all-transcripts       All transcripts (GENCODE)
  --basic-protein-coding  Basic protein-coding transcripts (GENCODE)
  --basic-all-transcripts Basic all transcripts (GENCODE)
  
  Pass-through Arguments:
  --args-for-download ...
  --args-for-fastqc ...
  --args-for-salmon ...
  --args-for-edger ...
  --args-for-drimseq ...
  --args-for-qapa ...

  -h, --help              Display this help message.
EOM
  exit 1
}

# --- Check Required Scripts ---
check_required_scripts() {
    local missing=0
    local required_scripts=("$FASTQC_SCRIPT" "$SALMON_SCRIPT" "$EDGER_SCRIPT")
    
    # downloadスクリプトは、local_fastqが未指定 かつ start-atがdownloadの時のみ必須とする
    if [[ -z "$LOCAL_FASTQ_DIR" && "$START_AT" == "download" ]]; then 
        required_scripts+=("$DOWNLOAD_SCRIPT")
    fi
    if [[ "$RUN_DRIMSEQ" == true || "$RUN_QAPA" == true ]]; then required_scripts+=("$DRIMSEQ_SCRIPT"); fi
    if [[ "$RUN_QAPA" == true ]]; then required_scripts+=("$QAPA_SCRIPT"); fi

    for script in "${required_scripts[@]}"; do
        if [[ ! -f "$script" ]]; then echo "ERROR: Required script not found: $script" >&2; missing=1; fi
    done
    if [[ $missing -eq 1 ]]; then exit 1; fi
}

# --- Argument Parsing ---
if [[ "$1" == "-h" || "$1" == "--help" ]]; then usage; fi
if [[ $# -lt 3 ]]; then echo "ERROR: Missing required arguments." >&2; usage; fi

INPUT_TARGET="$1"; shift
SPECIES="$1"; shift
GENCODE_RELEASE="$1"; shift

# Parse Options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output-dir) OUTPUT_BASE_DIR="$2"; shift 2 ;;
    
    --retry)
        if [[ -n "$2" && "$2" =~ ^[0-9]+$ ]]; then RETRY_COUNT="$2"; shift 2;
        else RETRY_COUNT=2; shift 1; fi ;;
        
    --start-at)
        START_AT="$2"
        if [[ ! "$START_AT" =~ ^(download|fastqc|salmon|edger|drimseq|qapa)$ ]]; then
            echo "ERROR: Invalid --start-at step. Valid: download, fastqc, salmon, edger, drimseq, qapa" >&2
            exit 1
        fi
        
        if [[ "$START_AT" == "drimseq" ]]; then
            RUN_DRIMSEQ=true
        elif [[ "$START_AT" == "qapa" ]]; then
            RUN_QAPA=true
        fi
        
        shift 2 ;;
    
    --local-fastq) LOCAL_FASTQ_DIR="$2"; shift 2 ;;
    --run-drimseq) RUN_DRIMSEQ=true; shift ;;
    --run-qapa)    RUN_QAPA=true; shift ;;
    
    --local) EXEC_MODE="local"; shift ;;
    --all-transcripts|--basic-protein-coding|--basic-all-transcripts)
        if [[ -n "$TRANSCRIPT_TYPE_FLAG" ]]; then echo "ERROR: Only one transcript type option allowed." >&2; exit 1; fi
        TRANSCRIPT_TYPE_FLAG="$1"; shift ;;
    
    --args-for-download) CONTEXT="download"; shift ;;
    --args-for-fastqc)   CONTEXT="fastqc"; shift ;;
    --args-for-salmon)   CONTEXT="salmon"; shift ;;
    --args-for-edger)    CONTEXT="edger"; shift ;;
    --args-for-drimseq)  CONTEXT="drimseq"; shift ;;
    --args-for-qapa)     CONTEXT="qapa"; shift ;;
    --args-for-drimseq_qapa) CONTEXT="drimseq_qapa"; shift ;;
    
    *)
      if [[ -n "$CONTEXT" ]]; then
        case "$CONTEXT" in
          download) DOWNLOAD_ARGS+="$1 " ;; 
          fastqc)   FASTQC_ARGS+="$1 "   ;; 
          salmon)   SALMON_ARGS+="$1 "   ;; 
          edger)    EDGER_ARGS+="$1 "    ;;
          drimseq)  DRIMSEQ_ARGS+="$1 "  ;; 
          qapa)     QAPA_ARGS+="$1 "     ;; 
          drimseq_qapa) DRIMSEQ_QAPA_ARGS+="$1 " ;;
        esac; shift
      else echo "ERROR: Unknown option: $1" >&2; usage; fi ;;
  esac
done

# --- QAPA Compatibility Check ---
if [[ "$RUN_QAPA" == true ]]; then
    if [[ "$SPECIES" == "mouse" && "$GENCODE_RELEASE" != "22" ]]; then
        echo "ERROR: QAPA analysis (--run-qapa) for mouse is currently restricted to GENCODE release '22'." >&2
        exit 1
    elif [[ "$SPECIES" == "human" && "$GENCODE_RELEASE" != "31" ]]; then
        echo "ERROR: QAPA analysis (--run-qapa) for human is currently restricted to GENCODE release '31'." >&2
        exit 1
    elif [[ "$SPECIES" != "mouse" && "$SPECIES" != "human" ]]; then
        echo "ERROR: QAPA analysis (--run-qapa) is only supported for 'human' and 'mouse'." >&2
        exit 1
    fi
fi

check_required_scripts

# --- Initial Logging ---
echo "========================================================================"
echo "INFO: RNSequant Pipeline Started (v1.0)"
echo "INFO: Target: $INPUT_TARGET"
echo "INFO: Species: $SPECIES, Release: $GENCODE_RELEASE"
if [[ -n "$LOCAL_FASTQ_DIR" ]]; then echo "INFO: Mode: Local FastQ Input ($LOCAL_FASTQ_DIR)"; fi
echo "INFO: Options: DRIMSeq=$RUN_DRIMSEQ, QAPA=$RUN_QAPA"
echo "INFO: Starting at step: $START_AT"
echo "========================================================================"

prepare_local_fastq() {
    local sample_id="$1"
    local source_dir="$2"
    local dest_dir="$3"
    
    echo "INFO: Searching for local FastQ files for sample: $sample_id"
    local files=($(find "$source_dir" -maxdepth 1 -type f -name "*${sample_id}*" | sort))
    local count=${#files[@]}
    
    if [[ $count -eq 0 ]]; then
        echo "ERROR: No FastQ files found for Sample ID '$sample_id' in $source_dir" >&2
        return 1
    elif [[ $count -eq 1 ]]; then
        echo "  -> Found 1 file (Single-End detected): ${files[0]}"
        ln -sf "${files[0]}" "$dest_dir/${sample_id}.fastq.gz"
    elif [[ $count -eq 2 ]]; then
        echo "  -> Found 2 files (Paired-End detected):"
        echo "     R1: ${files[0]}"
        echo "     R2: ${files[1]}"
        ln -sf "${files[0]}" "$dest_dir/${sample_id}_1.fastq.gz"
        ln -sf "${files[1]}" "$dest_dir/${sample_id}_2.fastq.gz"
    else
        echo "ERROR: Ambiguous file match. Found $count files for '$sample_id'." >&2
        return 1
    fi
    return 0
}

# ==============================================================================
# Core Pipeline Function
# ==============================================================================
run_pipeline_logic() {
    local samples_file="$1"
    
    local samples_file_abs="$(cd "$(dirname "$samples_file")" && pwd)/$(basename "$samples_file")"
    local prefix=$(basename "$samples_file" .csv)
    local ACTIVE_CSV="$(basename "$samples_file")"
    
    local output_dir=""
    
    if [[ -n "$OUTPUT_BASE_DIR" ]]; then
        output_dir="${OUTPUT_BASE_DIR}"
    else
        output_dir="$(pwd)"
    fi

    if [[ "$(basename "$output_dir")" != "$prefix" ]]; then
        output_dir="${output_dir}/${prefix}"
    fi
    
    mkdir -p "$output_dir"
    local output_dir_abs="$(cd "$output_dir" && pwd)"
    
    cp "$samples_file_abs" "$output_dir_abs/$ACTIVE_CSV"

    local local_flag=""; if [[ "$EXEC_MODE" == "local" ]]; then local_flag="--local"; fi

    local species_p=""; if [[ "$SPECIES" == "mouse" ]]; then species_p="M"; fi
    local idx_s="pc"
    if [[ "$TRANSCRIPT_TYPE_FLAG" == "--all-transcripts" ]]; then idx_s="all";
    elif [[ "$TRANSCRIPT_TYPE_FLAG" == "--basic-protein-coding" ]]; then idx_s="pc.basic";
    elif [[ "$TRANSCRIPT_TYPE_FLAG" == "--basic-all-transcripts" ]]; then idx_s="all.basic"; fi
    
    local version_prefix="v${species_p}${GENCODE_RELEASE}.${idx_s}"
    local qapa_version_prefix="v${species_p}${GENCODE_RELEASE}"

    cd "$output_dir_abs" || exit 1

    # Check for previously merged CSV (Silent Handoff for resumed runs)
    if [[ -f "merged_${ACTIVE_CSV}" ]]; then
        echo "INFO: Auto-detected previously merged CSV: merged_${ACTIVE_CSV}"
        ACTIVE_CSV="merged_${ACTIVE_CSV}"
    fi

    # --- Pre-Step: Handle Local FastQ ---
    # Only run local fastq symlinking if we are actually starting at an early step
    if [[ -n "$LOCAL_FASTQ_DIR" && "$START_AT" =~ ^(download|fastqc|salmon)$ ]]; then
        mkdir -p fastq
        while IFS=, read -r col1 rest; do
             if [[ -z "$col1" || "$col1" == "file_name" ]]; then continue; fi
             prepare_local_fastq "$col1" "$LOCAL_FASTQ_DIR" "$(pwd)/fastq" || exit 1
        done < "$ACTIVE_CSV"
    fi

    # --- Pipeline Steps ---
    local steps=("download" "fastqc" "salmon" "edger" "drimseq" "qapa" "drimseq_qapa")
    local step_active=false

    for step in "${steps[@]}"; do
        # Activate execution once we reach the START_AT step
        if [[ "$step" == "$START_AT" ]]; then
            step_active=true
        fi
        
        if [[ "$step_active" == false ]]; then
            echo "INFO: [Step: $step] Skipped (--start-at $START_AT)"
            continue
        fi

        # 1. Skip download if using local fastq
        if [[ "$step" == "download" && -n "$LOCAL_FASTQ_DIR" ]]; then continue; fi
        # 2. Skip DRIMSeq if not requested
        if [[ "$step" == "drimseq" && "$RUN_DRIMSEQ" == false ]]; then continue; fi
        # 3. Skip QAPA if not requested
        if [[ "$step" == "qapa" && "$RUN_QAPA" == false ]]; then continue; fi
        # 4. Skip Post-QAPA DRIMSeq if QAPA is not requested
        if [[ "$step" == "drimseq_qapa" && "$RUN_QAPA" == false ]]; then continue; fi

        echo "----------------------------------------------------"
        echo "INFO: [Step: $step] Processing $prefix"
        
        case "$step" in
            "download")
                bash "$DOWNLOAD_SCRIPT" "$ACTIVE_CSV" -o fastq $DOWNLOAD_ARGS $local_flag || exit 1
                if [[ -f "fastq/merged_${ACTIVE_CSV}" ]]; then
                    cp "fastq/merged_${ACTIVE_CSV}" "merged_${ACTIVE_CSV}"
                    ACTIVE_CSV="merged_${ACTIVE_CSV}"
                fi
                ;;
            "fastqc")
                bash "$FASTQC_SCRIPT" "$ACTIVE_CSV" -i "$(pwd)/fastq" -o trimed_fastq $FASTQC_ARGS $local_flag || exit 1
                ;;
            "salmon")
                bash "$SALMON_SCRIPT" "$ACTIVE_CSV" "$SPECIES" "$GENCODE_RELEASE" \
                    --input-fastq-dir "$(pwd)/trimed_fastq" -o . \
                    $TRANSCRIPT_TYPE_FLAG $SALMON_ARGS $local_flag || exit 1
                ;;
            "edger")
                local counts_file="${version_prefix}_salmon_quant/${prefix}_${version_prefix}_formatted_for_edgeR_counts.tsv"
                local edger_out="${version_prefix}_edger_results"
                local gene_map="${version_prefix}_gencode_files/geneid2symbol_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"
                mkdir -p "$edger_out"
                bash "$EDGER_SCRIPT" "$ACTIVE_CSV" "$counts_file" \
                    --prefix "$edger_out/$prefix" \
                    --gene-symbol-map "$gene_map" \
                    $local_flag $EDGER_ARGS || exit 1
                ;;
            "drimseq")
                local salmon_counts="${version_prefix}_salmon_quant/${prefix}_transcript_counts_matrix.tsv"
                local tx2gene="${version_prefix}_gencode_files/tx2gene_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"
                bash "$DRIMSEQ_SCRIPT" "$ACTIVE_CSV" "$SPECIES" "$GENCODE_RELEASE" \
                    "$salmon_counts" $TRANSCRIPT_TYPE_FLAG \
                    --output-dir "${version_prefix}_drimseq" --tx2gene "$tx2gene" $DRIMSEQ_ARGS $local_flag || exit 1
                ;;
            "qapa")
                bash "$QAPA_SCRIPT" "$ACTIVE_CSV" "$SPECIES" "$GENCODE_RELEASE" \
                    --input-dir "$(pwd)/trimed_fastq" $QAPA_ARGS $local_flag || exit 1
                ;;
            "drimseq_qapa")
                local qapa_counts="${qapa_version_prefix}_qapa/qapa_transcript_counts_matrix.tsv"
                local qapa_tx2gene="${qapa_version_prefix}_qapa_fasta/apaid2geneid.tsv"
                local qapa_gene_map="${qapa_version_prefix}_qapa_fasta/geneid2genesymbol.tsv"
                
                bash "$DRIMSEQ_SCRIPT" "$ACTIVE_CSV" "$SPECIES" "$GENCODE_RELEASE" \
                    "$qapa_counts" $TRANSCRIPT_TYPE_FLAG \
                    --tx2gene "$qapa_tx2gene" --add-gene-symbol "$qapa_gene_map" \
                    --output-dir "${qapa_version_prefix}_qapa_drimseq" \
                    $DRIMSEQ_QAPA_ARGS $local_flag || exit 1
                ;;
        esac
    done
}

# ==============================================================================
# Retry Wrapper & Main Execution
# ==============================================================================
run_with_retry() {
    local target_csv="$1"
    local attempt=1
    local max_attempts=$((RETRY_COUNT + 1))
    local csv_name=$(basename "$target_csv")
    
    while [[ $attempt -le $max_attempts ]]; do
        echo ""
        echo "========================================================================"
        echo "INFO: PROCESSING: $csv_name (Attempt $attempt/$max_attempts)"
        echo "========================================================================"
        
        if (run_pipeline_logic "$target_csv"); then
            echo "INFO: SUCCESS: Pipeline finished for $csv_name"
            return 0
        else
            echo "WARNING: Pipeline failed for $csv_name at attempt $attempt."
            if [[ $attempt -lt $max_attempts ]]; then
                echo "INFO: Retrying in 5 seconds..."; sleep 5; ((attempt++))
            else
                echo "ERROR: Max retries reached for $csv_name." >&2; return 1
            fi
        fi
    done
}

if [[ -f "$INPUT_TARGET" ]]; then
    echo "INFO: Single file mode detected."
    if ! run_with_retry "$INPUT_TARGET"; then echo "CRITICAL: Pipeline failed." >&2; exit 1; fi
elif [[ -d "$INPUT_TARGET" ]]; then
    echo "INFO: Batch directory mode detected."
    count=0; fail_count=0
    while IFS= read -r csv_file; do
        ((count++))
        if ! run_with_retry "$csv_file"; then echo "CRITICAL ERROR: Failed $csv_file" >&2; ((fail_count++)); fi
    done < <(find "$INPUT_TARGET" -maxdepth 1 -name "*.csv" | sort)
    
    if [[ $count -eq 0 ]]; then echo "ERROR: No CSV files found." >&2; exit 1; fi
    if [[ $fail_count -gt 0 ]]; then exit 1; else exit 0; fi
else
    echo "ERROR: Input '$INPUT_TARGET' is invalid." >&2; exit 1
fi

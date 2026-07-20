#!/bin/bash
# ==============================================================================
# Script Name: run_salmon.sh
# Version: 1.0
# Description:
#   Automates the RNA-Seq quantification workflow using Salmon.
#   - Auto-accepts Paired CSV design if 'subject' column is present.
#   - Downloads GENCODE annotations (Human/Mouse).
#   - Supports Comprehensive and Basic gene annotations (handling via awk).
#   - Quantifies expression (Gene & Transcript level).
#   - Aggregates results using tximport (TPM, ScaledTPM, Raw Counts).
#   - Generates QC reports using MultiQC.
#   - Formats count matrices for downstream analysis (e.g., edgeR).
# ==============================================================================

set -e
set -o pipefail
set -u

# --- Global Variables (Container Images) ---
SALMON_CONTAINER="quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5"
TXIMPORT_CONTAINER="ezojika7713/tximport:v1.22"
MULTIQC_CONTAINER="quay.io/biocontainers/multiqc:1.10.1--py_0"

# --- Help Message ---
usage() {
  cat <<EOM
Usage: $(basename "$0") [OPTIONS] <samples.csv> <species> <gencode_release>
   or: $(basename "$0") -h | -v

Description:
  Automates Salmon quantification and QC. All outputs (raw matrices, formatted
  files for edgeR, and logs) are saved into a version-specific salmon
  quantification directory.

Required Arguments:
  <samples.csv>     Sample information CSV file. Header must be exactly one of:
                      Unpaired: file_name,sample_name,condition
                      Paired:   file_name,sample_name,condition,subject
  <species>         Annotation species. Supported: "human" or "mouse".
  <gencode_release> GENCODE release version number (e.g., 46 for human, 36 for mouse).

Optional Flags:
  --input-fastq-dir <path>  Directory containing FASTQ files (default: current directory).
  -o, --output-dir <path>   Base directory for all output files (default: current directory).
  --threads <int>           Number of threads to use for Salmon (default: 8).
  --all-transcripts         Use all transcripts (including non-coding). Default is protein-coding only.
  --basic-protein-coding    Use GENCODE basic protein-coding transcripts (creates FASTA from basic GTF).
  --basic-all-transcripts   Use GENCODE basic all transcripts (creates FASTA from basic GTF).
  --local                   Use locally installed tools (salmon, R, multiqc) instead of Docker.
  -v                        Display tool versions and exit.
  -h                        Display this help message and exit.

Note: Only one of --all-transcripts, --basic-protein-coding, or --basic-all-transcripts can be specified.
EOM
  exit 1
}

# --- 1. Argument Parsing & Configuration ---
if [[ "$1" == "-h" || "$1" == "--help" ]]; then usage; fi
if [[ "$1" == "-v" || "$1" == "--version" ]]; then
    echo "--- Querying Tool Versions ---"
    echo "This script uses:"
    echo "  - Salmon (Docker: ${SALMON_CONTAINER})"
    echo "  - tximport (Docker: ${TXIMPORT_CONTAINER})"
    echo "  - MultiQC (Docker: ${MULTIQC_CONTAINER})"
    exit 0
fi

EXEC_MODE="docker"
ANNOTATION_TYPE="comprehensive"  # comprehensive, basic
TRANSCRIPT_TYPE="pc_"
TRANSCRIPT_TYPE_MSG="protein-coding (default)"
THREADS=8
FASTQ_DIR="."
OUTPUT_BASE_DIR="."
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local)
      EXEC_MODE="local"
      shift
      ;;
    --all-transcripts)
      if [[ "$TRANSCRIPT_TYPE" != "pc_" ]]; then
        echo "ERROR: Cannot specify multiple transcript type options." >&2
        exit 1
      fi
      TRANSCRIPT_TYPE=""
      TRANSCRIPT_TYPE_MSG="all transcripts"
      shift
      ;;
    --basic-protein-coding)
      if [[ "$TRANSCRIPT_TYPE" != "pc_" ]]; then
        echo "ERROR: Cannot specify multiple transcript type options." >&2
        exit 1
      fi
      ANNOTATION_TYPE="basic"
      TRANSCRIPT_TYPE="pc_"
      TRANSCRIPT_TYPE_MSG="GENCODE basic protein-coding"
      shift
      ;;
    --basic-all-transcripts)
      if [[ "$TRANSCRIPT_TYPE" != "pc_" ]]; then
        echo "ERROR: Cannot specify multiple transcript type options." >&2
        exit 1
      fi
      ANNOTATION_TYPE="basic"
      TRANSCRIPT_TYPE=""
      TRANSCRIPT_TYPE_MSG="GENCODE basic all transcripts"
      shift
      ;;
    --input-fastq-dir)
      if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
        echo "ERROR: --input-fastq-dir option requires a directory path." >&2
        exit 1
      fi
      if [[ ! -d "$2" ]]; then
        echo "ERROR: FASTQ directory '$2' does not exist." >&2
        exit 1
      fi
      FASTQ_DIR="$2"
      shift 2
      ;;
    -o|--output-dir)
      if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
        echo "ERROR: -o/--output-dir option requires a directory path." >&2
        exit 1
      fi
      OUTPUT_BASE_DIR="$2"
      shift 2
      ;;
    --threads)
      if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
        echo "ERROR: --threads option requires a number." >&2
        exit 1
      fi
      if ! [[ "$2" =~ ^[0-9]+$ ]] || [[ "$2" -eq 0 ]]; then
        echo "ERROR: --threads must be a positive integer." >&2
        exit 1
      fi
      THREADS="$2"
      shift 2
      ;;
    -*)
      echo "ERROR: Unknown option '$1' in $(basename "$0")" >&2
      usage
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}"

if [[ $# -ne 3 ]]; then
    echo "ERROR: Missing or invalid required arguments: <samples.csv> <species> <gencode_release>" >&2
    usage
fi

SAMPLES_FILE=$1
SPECIES=$2
GENCODE_RELEASE=$3
PREFIX=$(basename "$SAMPLES_FILE" .csv)

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "ERROR: Sample file '$SAMPLES_FILE' not found." >&2
    exit 1
fi

# Check local dependencies if requested
if [[ "$EXEC_MODE" == "local" ]]; then
    for cmd in salmon Rscript multiqc wget; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "ERROR: In --local mode, '$cmd' must be installed and in your PATH." >&2
            exit 1
        fi
    done
fi

# --- 2. Version Prefix Generation ---
SPECIES_PREFIX=""
if [[ "$SPECIES" == "mouse" ]]; then
    SPECIES_PREFIX="M"
elif [[ "$SPECIES" != "human" ]]; then
    echo "ERROR: Species must be 'human' or 'mouse', got '$SPECIES'" >&2
    exit 1
fi

INDEX_SUFFIX="pc"
if [[ -z "$TRANSCRIPT_TYPE" ]]; then
    INDEX_SUFFIX="all"
fi

if [[ "$ANNOTATION_TYPE" == "basic" ]]; then
    INDEX_SUFFIX="${INDEX_SUFFIX}.basic"
fi

VERSION_PREFIX="v${SPECIES_PREFIX}${GENCODE_RELEASE}.${INDEX_SUFFIX}"

# --- 3. Define Paths and Create Directories ---
mkdir -p "$OUTPUT_BASE_DIR"

GENCODE_DIR="${OUTPUT_BASE_DIR}/${VERSION_PREFIX}_gencode_files"
SALMON_OUT_DIR="${OUTPUT_BASE_DIR}/${VERSION_PREFIX}_salmon_quant"
MULTIQC_DIR="${OUTPUT_BASE_DIR}/${VERSION_PREFIX}_multiqc_salmon"
LOG_FILE="${SALMON_OUT_DIR}/run_pipeline_${PREFIX}_${VERSION_PREFIX}.log"

# Output files
TPM_FILE="${SALMON_OUT_DIR}/${PREFIX}_gene_tpm_matrix.tsv"
SCALEDTPM_FILE="${SALMON_OUT_DIR}/${PREFIX}_gene_scaledtpm_matrix.tsv"
COUNTS_FILE="${SALMON_OUT_DIR}/${PREFIX}_gene_counts_matrix.tsv"
TRANSCRIPT_COUNTS_FILE="${SALMON_OUT_DIR}/${PREFIX}_transcript_counts_matrix.tsv"
GENE_SYMBOL_MAP="${GENCODE_DIR}/geneid2symbol_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"

EDGE_PREFIX="${PREFIX}_${VERSION_PREFIX}_formatted_for_edgeR"
EDGE_COUNTS_FILE="${SALMON_OUT_DIR}/${EDGE_PREFIX}_counts.tsv"

# --- Initialize Logging ---
mkdir -p "$SALMON_OUT_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "========================================================================"
echo "INFO: Salmon Quantification & QC Pipeline Started (v1.0)"
echo "INFO: Date               : $(date)"
echo "INFO: Execution Mode     : $EXEC_MODE"
echo "INFO: Annotation Type    : $ANNOTATION_TYPE"
echo "INFO: Transcript Type    : $TRANSCRIPT_TYPE_MSG"
echo "INFO: Threads            : $THREADS"
echo "INFO: FASTQ Directory    : $FASTQ_DIR"
echo "INFO: Output Directory   : $OUTPUT_BASE_DIR"
echo "INFO: Sample Info        : $SAMPLES_FILE"
echo "INFO: Prefix             : $PREFIX"
echo "INFO: Version Prefix     : $VERSION_PREFIX"
echo "INFO: Species            : $SPECIES"
echo "INFO: GENCODE Release    : $GENCODE_RELEASE"
echo "========================================================================"

# --- Strict Header Validation (Paired/Unpaired Auto-Detect) ---
EXPECTED_HEADER_UNPAIRED="file_name,sample_name,condition"
EXPECTED_HEADER_PAIRED="file_name,sample_name,condition,subject"
ACTUAL_HEADER=$(head -1 "${SAMPLES_FILE}" | tr -d '\r')

if [[ "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_UNPAIRED}" && "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_PAIRED}" ]]; then
  echo "ERROR: Invalid CSV header format in ${SAMPLES_FILE}" >&2
  echo "       Expected:" >&2
  echo "         Unpaired: ${EXPECTED_HEADER_UNPAIRED}" >&2
  echo "         Paired:   ${EXPECTED_HEADER_PAIRED}" >&2
  echo "       Actual:   ${ACTUAL_HEADER}" >&2
  exit 1
fi
echo "INFO: CSV header validated: ${ACTUAL_HEADER}"

# --- 4. Tool Version Check ---
echo "----------------------------------------------------"
echo "INFO: Step 1/7 - Recording Tool Versions (Mode: $EXEC_MODE)"
if [[ "$EXEC_MODE" == "docker" ]]; then
    echo "INFO:   [Salmon]"
    docker run --rm "$SALMON_CONTAINER" salmon --version
    echo "INFO:   [MultiQC]"
    docker run --rm "$MULTIQC_CONTAINER" multiqc --version
    echo "INFO:   [R Packages in tximport container]"
    docker run --rm "$TXIMPORT_CONTAINER" Rscript -e 'cat("R version:", R.version.string, "\n\n"); suppressPackageStartupMessages({ libs <- c("tximport", "readr", "dplyr"); for (lib in libs) { if (require(lib, character.only=TRUE)) { cat(lib, ": ", as.character(packageVersion(lib)), "\n") } else { cat(lib, ": Not Installed\n") } } })'
else
    echo "INFO:   [Salmon]"
    salmon --version
    echo "INFO:   [MultiQC]"
    multiqc --version
    echo "INFO:   [R Packages]"
    Rscript -e 'cat("R version:", R.version.string, "\n\n"); suppressPackageStartupMessages({ libs <- c("tximport", "readr", "dplyr"); for (lib in libs) { if (require(lib, character.only=TRUE)) { cat(lib, ": ", as.character(packageVersion(lib)), "\n") } else { cat(lib, ": Not Installed\n") } } })'
fi

# --- 5. Pre-flight Check (FASTQ Existence) ---
echo "----------------------------------------------------"
echo "INFO: Step 2/7 - Pre-flight check for FASTQ files"
missing_files=()
while IFS=, read -r file_name sample_name condition rest; do
    [[ -z "$file_name" ]] && continue
    
    if [[ -f "${SALMON_OUT_DIR}/${file_name}/quant.sf" ]]; then
        echo "INFO:   - quant.sf found for '$file_name'. Skipping FASTQ check."
        continue
    fi

    found_count=$(find "$FASTQ_DIR" -maxdepth 1 \( -name "*${file_name}*trimmed.f*q.gz" -o -name "*${file_name}*_val_[12].f*q.gz" -o -name "*${file_name}*.f*q.gz" \) 2>/dev/null | wc -l)
    if [[ "$found_count" -eq 0 ]]; then
        missing_files+=("$file_name")
    fi
done < <(tail -n +2 "$SAMPLES_FILE" | tr -d '\r')

if [[ ${#missing_files[@]} -gt 0 ]]; then
    echo "ERROR: Required FASTQ files could not be found in '$FASTQ_DIR'." >&2
    for fname in "${missing_files[@]}"; do
        echo "  - Files for '$fname' are missing." >&2
    done
    exit 1
else
    echo "INFO:   - All required FASTQ files (or existing Salmon results) are present."
fi

# --- 6. GENCODE File Preparation ---
echo "----------------------------------------------------"
echo "INFO: Step 3/7 - Preparing GENCODE files"
mkdir -p "$GENCODE_DIR"

if [[ "$SPECIES" == "human" ]]; then
    FTP_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
    if [[ "$ANNOTATION_TYPE" == "basic" ]]; then
        ANNOTATION_GTF="gencode.v${GENCODE_RELEASE}.basic.annotation.gtf.gz"
        SOURCE_FASTA="gencode.v${GENCODE_RELEASE}.pc_transcripts.fa.gz"
    else
        TRANSCRIPT_FASTA="gencode.v${GENCODE_RELEASE}.${TRANSCRIPT_TYPE}transcripts.fa.gz"
        ANNOTATION_GTF="gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"
    fi
    SALMON_INDEX_DIR="${GENCODE_DIR}/gencode_v${GENCODE_RELEASE}_human_${INDEX_SUFFIX}_index"
else
    FTP_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${GENCODE_RELEASE}"
    if [[ "$ANNOTATION_TYPE" == "basic" ]]; then
        ANNOTATION_GTF="gencode.vM${GENCODE_RELEASE}.basic.annotation.gtf.gz"
        SOURCE_FASTA="gencode.vM${GENCODE_RELEASE}.pc_transcripts.fa.gz"
    else
        TRANSCRIPT_FASTA="gencode.vM${GENCODE_RELEASE}.${TRANSCRIPT_TYPE}transcripts.fa.gz"
        ANNOTATION_GTF="gencode.vM${GENCODE_RELEASE}.annotation.gtf.gz"
    fi
    SALMON_INDEX_DIR="${GENCODE_DIR}/gencode_vM${GENCODE_RELEASE}_mouse_${INDEX_SUFFIX}_index"
fi

TX2GENE_MAP="${GENCODE_DIR}/tx2gene_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"

# Download Annotation GTF
if [[ -s "${GENCODE_DIR}/${ANNOTATION_GTF}" ]]; then
    echo "INFO:   - Annotation GTF found. Skipping download."
else
    echo "INFO:   - Downloading annotation GTF..."
    wget -q -c -P "$GENCODE_DIR" "${FTP_BASE}/${ANNOTATION_GTF}"
fi

# Prepare FASTA
if [[ "$ANNOTATION_TYPE" == "basic" ]]; then
    GENERATED_FASTA="${GENCODE_DIR}/gencode_basic_${TRANSCRIPT_TYPE}transcripts.fa"
    
    if [[ -s "${GENCODE_DIR}/${SOURCE_FASTA}" ]]; then
        echo "INFO:   - Comprehensive protein-coding FASTA found. Skipping download."
    else
        echo "INFO:   - Downloading comprehensive protein-coding FASTA as source..."
        wget -q -c -P "$GENCODE_DIR" "${FTP_BASE}/${SOURCE_FASTA}"
    fi
    
    if [[ -s "${GENERATED_FASTA}.gz" ]]; then
        echo "INFO:   - Generated transcript FASTA found. Skipping filtering."
    else
        echo "INFO:   - Extracting basic transcript IDs from GTF..."
        BASIC_TRANSCRIPT_IDS="${GENCODE_DIR}/basic_transcript_ids.txt"
        zcat "${GENCODE_DIR}/${ANNOTATION_GTF}" | \
            awk -F'\t' '$3 == "transcript" { match($9, /transcript_id "([^"]+)"/); print substr($9, RSTART+15, RLENGTH-16) }' | \
            sed 's/\..*//' > "$BASIC_TRANSCRIPT_IDS"
        
        echo "INFO:   - Filtering FASTA sequences (keeping only basic transcripts)..."
        zcat "${GENCODE_DIR}/${SOURCE_FASTA}" | \
            awk -v ids_file="$BASIC_TRANSCRIPT_IDS" '
            BEGIN {
                while ((getline line < ids_file) > 0) {
                    ids[line] = 1
                }
                close(ids_file)
                print_seq = 0
            }
            /^>/ {
                match($0, /^>([^ |]+)/, arr)
                tid = arr[1]
                gsub(/\.[0-9]+$/, "", tid)
                if (tid in ids) {
                    print_seq = 1
                    print $0
                } else {
                    print_seq = 0
                }
                next
            }
            print_seq { print }
            ' | gzip > "${GENERATED_FASTA}.gz"
        
        echo "INFO:   - Generated filtered FASTA: ${GENERATED_FASTA}.gz"
        rm -f "$BASIC_TRANSCRIPT_IDS"
    fi
    TRANSCRIPT_FASTA="$(basename "$GENERATED_FASTA").gz"
else
    if [[ -s "${GENCODE_DIR}/${TRANSCRIPT_FASTA}" ]]; then
        echo "INFO:   - Transcript FASTA found. Skipping download."
    else
        echo "INFO:   - Downloading transcript FASTA ($TRANSCRIPT_TYPE_MSG)..."
        wget -q -c -P "$GENCODE_DIR" "${FTP_BASE}/${TRANSCRIPT_FASTA}"
    fi
fi

# Create Mapping Files
if [[ -s "$TX2GENE_MAP" ]]; then
    echo "INFO:   - tx2gene map found. Skipping creation."
else
    echo "INFO:   - Creating tx2gene map..."
    zcat "${GENCODE_DIR}/${ANNOTATION_GTF}" | awk -F'\t' '$3 == "transcript" { match($9, /transcript_id "([^"]+)"/); tid = substr($9, RSTART+15, RLENGTH-16); match($9, /gene_id "([^"]+)"/); gid = substr($9, RSTART+9, RLENGTH-10); print tid "\t" gid }' > "$TX2GENE_MAP"
fi

if [[ -s "$GENE_SYMBOL_MAP" ]]; then
    echo "INFO:   - gene_id to gene_symbol map found. Skipping creation."
else
    echo "INFO:   - Creating gene_id to gene_symbol map..."
    zcat "${GENCODE_DIR}/${ANNOTATION_GTF}" | awk -F'\t' '$3 == "gene" { match($9, /gene_id "([^"]+)"/); gid = substr($9, RSTART+9, RLENGTH-10); match($9, /gene_name "([^"]+)"/); gsymbol = substr($9, RSTART+11, RLENGTH-12); print gid "\t" gsymbol }' > "$GENE_SYMBOL_MAP"
fi

# --- 7. Build Salmon Index ---
echo "----------------------------------------------------"
echo "INFO: Step 4/7 - Building Salmon Index"
if [[ -d "$SALMON_INDEX_DIR" && -f "${SALMON_INDEX_DIR}/pos.bin" ]]; then
    echo "INFO:   - Salmon index for transcript type '$TRANSCRIPT_TYPE_MSG' found at '$SALMON_INDEX_DIR'. Skipping."
else
    echo "INFO:   - Salmon index not found or incomplete. Building index..."
    if [[ "$EXEC_MODE" == "docker" ]]; then
        GENCODE_DIR_ABS=$(cd "$GENCODE_DIR" && pwd)
        docker run --rm -u "$(id -u):$(id -g)" -v "${GENCODE_DIR_ABS}":/data "$SALMON_CONTAINER" salmon index --threads "$THREADS" -t "/data/${TRANSCRIPT_FASTA}" -i "/data/$(basename "$SALMON_INDEX_DIR")" --gencode
    else
        salmon index --threads "$THREADS" -t "${GENCODE_DIR}/${TRANSCRIPT_FASTA}" -i "$SALMON_INDEX_DIR" --gencode
    fi
fi

# --- 8. Salmon Quantification ---
echo "----------------------------------------------------"
echo "INFO: Step 5/7 - Running Salmon Quantification"
while IFS=, read -r file_name sample_name condition rest; do
    [[ -z "$file_name" ]] && continue
    echo "INFO:   -> Processing sample: $sample_name ($file_name)"
    
    if [[ -f "${SALMON_OUT_DIR}/${file_name}/quant.sf" ]]; then
        echo "INFO:      - Output for $file_name found. Skipping."
        continue
    fi
    
    trimmed_files=()
    while IFS= read -r -d '' file; do
        trimmed_files+=("$file")
    done < <(find "$FASTQ_DIR" -maxdepth 1 \( -name "*${file_name}*_val_[12].f*q.gz" -o -name "*${file_name}*trimmed.f*q.gz" \) -print0 | sort -z)
    
    if [[ ${#trimmed_files[@]} -gt 0 ]]; then
        echo "INFO:      - Found trimmed FASTQ files."
        read_files=("${trimmed_files[@]}")
    else
        echo "INFO:      - No trimmed files found. Using raw FASTQ files."
        read_files=()
        while IFS= read -r -d '' file; do
            read_files+=("$file")
        done < <(find "$FASTQ_DIR" -maxdepth 1 -name "*${file_name}*.f*q.gz" -print0 | sort -z)
    fi
    num_files=${#read_files[@]}

    salmon_read_opts_host=""
    salmon_read_opts_docker=""

    if [[ "$num_files" -eq 1 ]]; then
        echo "INFO:      - Detected Single-End reads."
        salmon_read_opts_host="-r ${read_files[0]}"
        salmon_read_opts_docker="-r /fastq/$(basename "${read_files[0]}")"
    elif [[ "$num_files" -eq 2 ]]; then
        echo "INFO:      - Detected Paired-End reads."
        read1=$(printf '%s\n' "${read_files[@]}" | grep -E "(_R1|_1|\.1\.)" || echo "${read_files[0]}")
        read2=$(printf '%s\n' "${read_files[@]}" | grep -E "(_R2|_2|\.2\.)" || echo "${read_files[1]}")
        salmon_read_opts_host="-1 ${read1} -2 ${read2}"
        salmon_read_opts_docker="-1 /fastq/$(basename "$read1") -2 /fastq/$(basename "$read2")"
    else
        echo "ERROR: Found $num_files files for $file_name (expected 1 or 2). Exiting." >&2
        exit 1
    fi

    echo "INFO:      - Running Salmon..."
    if [[ "$EXEC_MODE" == "docker" ]]; then
        FASTQ_DIR_ABS=$(cd "$FASTQ_DIR" && pwd)
        GENCODE_DIR_ABS=$(cd "$GENCODE_DIR" && pwd)
        SALMON_OUT_DIR_ABS=$(cd "$SALMON_OUT_DIR" && pwd)
        
        docker run --rm -u "$(id -u):$(id -g)" \
            -v "${FASTQ_DIR_ABS}":/fastq:ro \
            -v "${GENCODE_DIR_ABS}":/gencode:ro \
            -v "${SALMON_OUT_DIR_ABS}":/output \
            -w /output \
            "$SALMON_CONTAINER" \
            salmon quant -i "/gencode/$(basename "$SALMON_INDEX_DIR")" -l A ${salmon_read_opts_docker} -p "$THREADS" --validateMappings --gcBias -o "/output/${file_name}"
    else
        salmon quant -i "$SALMON_INDEX_DIR" -l A ${salmon_read_opts_host} -p "$THREADS" --validateMappings --gcBias -o "${SALMON_OUT_DIR}/${file_name}"
    fi
    echo "INFO:      - Finished sample: $file_name"
done < <(tail -n +2 "$SAMPLES_FILE" | tr -d '\r')

# --- 9. Aggregation & Formatting (R/tximport) ---
echo "----------------------------------------------------"
echo "INFO: Step 6/7 - Aggregating Results and Formatting for edgeR"

if [[ -s "$TPM_FILE" && -s "$COUNTS_FILE" && -s "$SCALEDTPM_FILE" && -s "$EDGE_COUNTS_FILE" && -s "$TRANSCRIPT_COUNTS_FILE" ]]; then
    echo "INFO:   - All aggregated matrix files already exist. Skipping tximport step."
else
    R_SCRIPT_CONTENT=$(cat <<'EOF'
        suppressPackageStartupMessages({
            library(tximport)
            library(readr)
            library(dplyr)
            library(tibble)
        })

        args <- commandArgs(trailingOnly = TRUE)
        samples_file_path           <- args[1]
        salmon_out_dir              <- args[2]
        tx2gene_map_path            <- args[3]
        tpm_output_path             <- args[4]
        counts_output_path          <- args[5]
        scaledtpm_output_path       <- args[6]
        edge_counts_path            <- args[7]
        transcript_counts_output_path <- args[8]

        sample_info <- read.csv(samples_file_path, header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
        tx2gene <- read_tsv(tx2gene_map_path, col_names = c("TXNAME", "GENEID"), show_col_types = FALSE)

        files <- file.path(salmon_out_dir, sample_info$file_name, "quant.sf")
        names(files) <- sample_info$file_name
        if (!all(file.exists(files))) {
            stop("One or more Salmon 'quant.sf' files are missing.")
        }

        cat("INFO:      - Importing gene-level counts\n")
        txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
        write_tsv(as.data.frame(txi$abundance) %>% rownames_to_column("gene_id"), tpm_output_path)

        txi_scaled <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE)
        write_tsv(as.data.frame(txi_scaled$counts) %>% rownames_to_column("gene_id"), scaledtpm_output_path)

        counts_data <- as.data.frame(round(txi$counts)) %>% rownames_to_column("gene_id")
        write_tsv(counts_data, counts_output_path)

        if (!all(sample_info$file_name %in% colnames(counts_data))) {
            missing_cols <- setdiff(sample_info$file_name, colnames(counts_data))
            stop(paste("CRITICAL ERROR: The following 'file_name'(s) from your CSV are MISSING as columns in the count matrix:", paste(missing_cols, collapse=", ")))
        }

        formatted_counts <- counts_data[, c("gene_id", sample_info$file_name)]
        colnames(formatted_counts) <- c("Geneid", sample_info$sample_name)
        write_tsv(formatted_counts, edge_counts_path)

        cat("INFO:      - Importing transcript-level counts for DTU analysis\n")
        txi_tx <- tximport(files, type = "salmon", txOut = TRUE, ignoreAfterBar = TRUE)
        tx_counts_data <- as.data.frame(round(txi_tx$counts)) %>% rownames_to_column("transcript_id")

        if (!all(sample_info$file_name %in% colnames(tx_counts_data))) {
            missing_cols_tx <- setdiff(sample_info$file_name, colnames(tx_counts_data))
            stop(paste("CRITICAL ERROR (Transcript-level): The following 'file_name'(s) from your CSV are MISSING as columns in the count matrix:", paste(missing_cols_tx, collapse=", ")))
        }

        final_tx_counts <- tx_counts_data[, c("transcript_id", sample_info$file_name)]
        write_tsv(final_tx_counts, transcript_counts_output_path)

        cat("INFO:      - Successfully generated matrices.\n")
EOF
)

    if [[ "$EXEC_MODE" == "docker" ]]; then
        SAMPLES_DIR_ABS=$(cd "$(dirname "$SAMPLES_FILE")" && pwd)
        OUTPUT_BASE_DIR_ABS=$(cd "$OUTPUT_BASE_DIR" && pwd)
        
        R_OUTPUT=$(docker run --rm -i -u "$(id -u):$(id -g)" \
            -v "${SAMPLES_DIR_ABS}:/samples_dir:ro" \
            -v "${OUTPUT_BASE_DIR_ABS}:/workdir" \
            -w /workdir \
            "$TXIMPORT_CONTAINER" \
            Rscript - "/samples_dir/$(basename "$SAMPLES_FILE")" \
            "$(realpath --relative-to="$OUTPUT_BASE_DIR" "$SALMON_OUT_DIR")" \
            "$(realpath --relative-to="$OUTPUT_BASE_DIR" "$TX2GENE_MAP")" \
            "$(realpath --relative-to="$OUTPUT_BASE_DIR" "$TPM_FILE")" \
            "$(realpath --relative-to="$OUTPUT_BASE_DIR" "$COUNTS_FILE")" \
            "$(realpath --relative-to="$OUTPUT_BASE_DIR" "$SCALEDTPM_FILE")" \
            "$(realpath --relative-to="$OUTPUT_BASE_DIR" "$EDGE_COUNTS_FILE")" \
            "$(realpath --relative-to="$OUTPUT_BASE_DIR" "$TRANSCRIPT_COUNTS_FILE")" <<< "$R_SCRIPT_CONTENT" 2>&1) || {
            echo "ERROR: R script for aggregation/formatting failed." >&2
            echo "--- R Error Message ---" >&2
            echo "$R_OUTPUT" >&2
            exit 1
        }
    else
        R_OUTPUT=$(Rscript - "$SAMPLES_FILE" "$SALMON_OUT_DIR" "$TX2GENE_MAP" "$TPM_FILE" "$COUNTS_FILE" "$SCALEDTPM_FILE" "$EDGE_COUNTS_FILE" "$TRANSCRIPT_COUNTS_FILE" <<< "$R_SCRIPT_CONTENT" 2>&1) || {
            echo "ERROR: R script for aggregation/formatting failed." >&2
            echo "--- R Error Message ---" >&2
            echo "$R_OUTPUT" >&2
            exit 1
        }
    fi
    echo "$R_OUTPUT"
fi

# --- 10. MultiQC Report Generation ---
echo "----------------------------------------------------"
echo "INFO: Step 7/7 - Generating MultiQC Report"
mkdir -p "$MULTIQC_DIR"

if [[ -f "${MULTIQC_DIR}/multiqc_report.html" ]]; then
    echo "INFO:   - MultiQC report already exists at '${MULTIQC_DIR}/multiqc_report.html'. Skipping."
else
    MULTIQC_SEARCH_DIRS=()

    RAW_FASTQ_CANDIDATES=(
        "./fastq"                           
        "../fastq"                          
        "$(dirname "${FASTQ_DIR}")/fastq"  
        "."                                 
    )

    echo "INFO:   - Searching for raw FastQC results..."
    RAW_FASTQ_FOUND=""
    for candidate in "${RAW_FASTQ_CANDIDATES[@]}"; do
        if [[ -d "$candidate" ]]; then
            fastqc_count=$(find "$candidate" -maxdepth 1 \( -name "*_fastqc.html" -o -name "*_fastqc.zip" \) 2>/dev/null | wc -l)
            if [[ $fastqc_count -gt 0 ]]; then
                RAW_FASTQ_FOUND="$candidate"
                echo "INFO:      -> Found raw FastQC results in: $RAW_FASTQ_FOUND ($fastqc_count files)"
                break
            fi
        fi
    done

    if [[ -n "$RAW_FASTQ_FOUND" ]]; then
        MULTIQC_SEARCH_DIRS+=("$RAW_FASTQ_FOUND")
    else
        echo "INFO:      -> No raw FastQC results found"
    fi

    echo "INFO:   - Searching for trimmed FastQC results..."
    TRIMMED_FASTQC_CANDIDATES=(
        "${FASTQ_DIR}"                      
        "./trimed_fastq"                    
    )

    TRIMMED_FASTQC_FOUND=""
    for candidate in "${TRIMMED_FASTQC_CANDIDATES[@]}"; do
        if [[ -d "$candidate" ]]; then
            fastqc_count=$(find "$candidate" -maxdepth 1 \( -name "*_fastqc.html" -o -name "*_fastqc.zip" \) 2>/dev/null | wc -l)
            if [[ $fastqc_count -gt 0 ]]; then
                TRIMMED_FASTQC_FOUND="$candidate"
                echo "INFO:      -> Found trimmed FastQC results in: $TRIMMED_FASTQC_FOUND ($fastqc_count files)"
                break
            fi
        fi
    done

    if [[ -n "$TRIMMED_FASTQC_FOUND" ]]; then
        MULTIQC_SEARCH_DIRS+=("$TRIMMED_FASTQC_FOUND")
    else
        echo "INFO:      -> No trimmed FastQC results found"
    fi

    echo "INFO:   - Searching for Salmon results..."
    if [[ -d "${SALMON_OUT_DIR}" ]]; then
        echo "INFO:      -> Found Salmon results in: ${SALMON_OUT_DIR}"
        MULTIQC_SEARCH_DIRS+=("${SALMON_OUT_DIR}")
    fi

    if [[ ${#MULTIQC_SEARCH_DIRS[@]} -eq 0 ]]; then
        echo "ERROR: No QC results found for MultiQC aggregation." >&2
        exit 1
    fi

    echo "INFO:   - MultiQC will aggregate results from ${#MULTIQC_SEARCH_DIRS[@]} location(s):"
    for dir in "${MULTIQC_SEARCH_DIRS[@]}"; do
        echo "INFO:      * $dir"
    done

    if [[ "$EXEC_MODE" == "docker" ]]; then
        DOCKER_MOUNTS=()
        DOCKER_PATHS=()
        
        for i in "${!MULTIQC_SEARCH_DIRS[@]}"; do
            abs_path=$(cd "${MULTIQC_SEARCH_DIRS[$i]}" && pwd)
            mount_point="/data/source_${i}"
            DOCKER_MOUNTS+=("-v" "${abs_path}:${mount_point}:ro")
            DOCKER_PATHS+=("${mount_point}")
        done
        
        MULTIQC_DIR_ABS=$(cd "$MULTIQC_DIR" && pwd)
        
        docker run --rm -u "$(id -u):$(id -g)" \
            "${DOCKER_MOUNTS[@]}" \
            -v "${MULTIQC_DIR_ABS}:/output" \
            -w /output \
            "$MULTIQC_CONTAINER" \
            multiqc "${DOCKER_PATHS[@]}" -o /output --force \
            --fullnames \
            --dirs --dirs-depth 1 > /dev/null 2>&1
    else
        multiqc "${MULTIQC_SEARCH_DIRS[@]}" -o "$MULTIQC_DIR" --force \
            --fullnames \
            --dirs --dirs-depth 1 > /dev/null 2>&1
    fi

    echo "INFO:   - MultiQC report saved in '$MULTIQC_DIR' directory."
fi

# --- 11. Final Report ---
echo "========================================================================"
echo "INFO: Pipeline Finished Successfully"
echo "INFO: All outputs are located in '$OUTPUT_BASE_DIR' under versioned directories:"
echo "INFO:"
echo "INFO:   1. Gene/Transcript Matrices (in $SALMON_OUT_DIR):"
echo "INFO:      - Counts (Gene)    : $(basename "$COUNTS_FILE")"
echo "INFO:      - Counts (Transcript): $(basename "$TRANSCRIPT_COUNTS_FILE")"
echo "INFO:      - TPM (Gene)       : $(basename "$TPM_FILE")"
echo "INFO:      - ScaledTPM (Gene) : $(basename "$SCALEDTPM_FILE")"
echo "INFO:"
echo "INFO:   2. Annotation Files (in $GENCODE_DIR):"
echo "INFO:      - Transcript-to-Gene Map: $(basename "$TX2GENE_MAP")"
echo "INFO:      - GeneID-to-Symbol Map  : $(basename "$GENE_SYMBOL_MAP")"
echo "INFO:"
echo "INFO:   3. Formatted Data for edgeR (in $SALMON_OUT_DIR):"
echo "INFO:      - Counts:  $(basename "$EDGE_COUNTS_FILE")"
echo "INFO:"
echo "INFO:   4. Quality Control Report:"
echo "INFO:      - Located in: $MULTIQC_DIR/multiqc_report.html"
echo "INFO:"
echo "INFO:   5. Log File (in $SALMON_OUT_DIR):"
echo "INFO:      - File: $(basename "$LOG_FILE")"
echo "========================================================================"
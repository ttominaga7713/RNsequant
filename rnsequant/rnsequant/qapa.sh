#!/bin/bash
set -e
set -o pipefail

# ==============================================================================
# Script Name: qapa.sh
# Version: 1.0
# Description: Automates the complete workflow for Alternative Polyadenylation (APA) analysis.
#              1. Prepares QAPA reference files (Human/Mouse).
#              2. Quantifies expression using Salmon.
#              3. Aggregates results using QAPA.
#              4. Generates count matrices formatted for DRIMSeq (DTU analysis).
#              - Auto-accepts Paired CSV design if 'subject' column is present.
# Supports:    Docker (default) and Local execution modes.
# ==============================================================================

# --- Global Configuration (Container Images) ---
SALMON_CONTAINER="quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5"
QAPA_CONTAINER="ezojika7713/qapa:1.3.3"
DRIMSEQ_CONTAINER="ezojika7713/drimseq:1.22.0"
TXIMPORT_CONTAINER="ezojika7713/tximport:v1.22"

# --- Help Message ---
usage() {
  cat <<EOM
Usage: $(basename "$0") <samples.csv> <species> <gencode_release> [OPTIONS]

Required Arguments:
  <samples.csv>     Mandatory CSV file. Header must be exactly one of:
                      Unpaired: file_name,sample_name,condition
                      Paired:   file_name,sample_name,condition,subject
  <species>         Target species: "human" or "mouse".
  <gencode_release> GENCODE release version (e.g., 31 for human, 22 for mouse).

Optional Arguments:
  --input-dir <path>  Directory containing FASTQ files (default: current directory ".").
  --threads <int>     Number of threads for parallel processing (default: 8).
  --local             Execute tools locally (requires salmon, qapa, Rscript, wget in PATH).
                      Default is to use Docker containers.
  -h, --help          Show this help message.
EOM
  exit 1
}

# --- 1. Argument Parsing & Setup ---
if [[ "$1" == "-h" || "$1" == "--help" ]]; then usage; fi
if [[ $# -lt 3 ]]; then echo "ERROR: Missing required arguments." >&2; usage; fi

SAMPLES_FILE=$1
SPECIES=$2
GENCODE_RELEASE=$3
PREFIX=$(basename "$SAMPLES_FILE" .csv)
shift 3

EXEC_MODE="docker"
THREADS=8
FASTQ_DIR="."

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir)
      if [[ -z "$2" ]] || [[ "$2" == -* ]]; then echo "ERROR: --input-dir requires a path." >&2; exit 1; fi
      FASTQ_DIR="$2"; shift 2 ;;
    --local)
      EXEC_MODE="local"; shift ;;
    --threads)
      if [[ -z "$2" ]] || [[ "$2" == -* ]]; then echo "ERROR: --threads requires a number." >&2; exit 1; fi
      THREADS="$2"; shift 2 ;;
    *) echo "ERROR: Unknown argument '$1'" >&2; usage ;;
  esac
done

# --- 2. Versioning & Absolute Directory Setup ---
SPECIES_PREFIX=""
if [[ "$SPECIES" == "mouse" ]]; then SPECIES_PREFIX="M"; fi
VERSION_PREFIX="v${SPECIES_PREFIX}${GENCODE_RELEASE}"

WORK_DIR="$(pwd)"
FASTQ_DIR="$(cd "$FASTQ_DIR" && pwd)"
SAMPLES_FILE_ABS="$(cd "$(dirname "$SAMPLES_FILE")" && pwd)/$(basename "$SAMPLES_FILE")"

MAIN_OUT_DIR="${WORK_DIR}/${VERSION_PREFIX}_qapa"
QAPA_FASTA_DIR="${WORK_DIR}/${VERSION_PREFIX}_qapa_fasta"
SALMON_OUT_DIR="$MAIN_OUT_DIR"
SALMON_INDEX_DIR="${QAPA_FASTA_DIR}/qapa_gencode_${SPECIES}_salmon_index"
LOG_FILE="${SALMON_OUT_DIR}/${PREFIX}_qapa_pipeline.log"
DB_DIR="${QAPA_FASTA_DIR}/db"

# --- 3. Pre-flight Checks ---
if [ ! -d "$FASTQ_DIR" ]; then echo "ERROR: FASTQ directory '$FASTQ_DIR' not found." >&2; exit 1; fi
if [ "$EXEC_MODE" == "local" ]; then
    for cmd in salmon qapa wget Rscript; do
        if ! command -v $cmd &> /dev/null; then echo "ERROR: '$cmd' required for local mode." >&2; exit 1; fi
    done
fi

# --- 4. Initialize Logging ---
mkdir -p "$SALMON_OUT_DIR" "$QAPA_FASTA_DIR" "$DB_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================================================"
echo "INFO: QAPA & Salmon Pipeline Started (v1.0)"
echo "INFO: Date: $(date)"
echo "INFO: Target CSV: $SAMPLES_FILE"
echo "INFO: Species: $SPECIES, Release: $GENCODE_RELEASE ($VERSION_PREFIX)"
echo "INFO: Mode: $EXEC_MODE, Threads: $THREADS"
echo "INFO: FastQ Directory: $FASTQ_DIR"
echo "========================================================================"

# --- Strict Header Validation (Paired/Unpaired Auto-Detect) ---
EXPECTED_HEADER_UNPAIRED="file_name,sample_name,condition"
EXPECTED_HEADER_PAIRED="file_name,sample_name,condition,subject"
ACTUAL_HEADER=$(head -1 "${SAMPLES_FILE_ABS}" | tr -d '\r')

if [[ "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_UNPAIRED}" && "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_PAIRED}" ]]; then
  echo "ERROR: Invalid CSV header format in ${SAMPLES_FILE}" >&2
  echo "       Expected:" >&2
  echo "         Unpaired: ${EXPECTED_HEADER_UNPAIRED}" >&2
  echo "         Paired:   ${EXPECTED_HEADER_PAIRED}" >&2
  echo "       Actual:   ${ACTUAL_HEADER}" >&2
  exit 1
fi
echo "INFO: CSV header validated: ${ACTUAL_HEADER}"

# --- 5. Prepare QAPA Resources (Download & Format) ---
echo "----------------------------------------------------"
echo "INFO: Step 1/9 - Preparing Reference Files..."

if [[ "$SPECIES" == "human" ]]; then
    ASSEMBLY_NAME="hg38"
    GENCODE_VERSION_PREFIX_URL="v${GENCODE_RELEASE}"
    GENOME_FASTA_BASENAME="GRCh38.primary_assembly.genome_${GENCODE_VERSION_PREFIX_URL}.fa"
    BED_BASENAME="qapa_3utrs.gencode_V31.hg38.bed"
    FTP_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
    FTP_GENOME_FASTA="GRCh38.primary_assembly.genome.fa.gz"
    BED_URL="https://github.com/morrislab/qapa/releases/download/v1.3.0/qapa_3utrs.gencode_V31.hg38.bed.gz"
    BIOMART_DATASET="hsapiens_gene_ensembl"
elif [[ "$SPECIES" == "mouse" ]]; then
    ASSEMBLY_NAME="mm10"
    GENCODE_VERSION_PREFIX_URL="vM${GENCODE_RELEASE}"
    GENOME_FASTA_BASENAME="GRCm38.primary_assembly.genome_${GENCODE_VERSION_PREFIX_URL}.fa"
    BED_BASENAME="qapa_3utrs.gencode_VM22.mm10.bed"
    FTP_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${GENCODE_RELEASE}"
    FTP_GENOME_FASTA="GRCm38.primary_assembly.genome.fa.gz"
    BED_URL="https://github.com/morrislab/qapa/releases/download/v1.3.0/qapa_3utrs.gencode_VM22.mm10.bed.gz"
    BIOMART_DATASET="mmusculus_gene_ensembl"
else
    echo "ERROR: Species '$SPECIES' not supported." >&2; exit 1
fi

BIOMART_URL='https://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="'${BIOMART_DATASET}'" interface="default"><Attribute name="ensembl_transcript_id" /><Attribute name="ensembl_gene_id" /><Attribute name="external_gene_name" /><Attribute name="gene_biotype" /></Dataset></Query>'

BED_FILE="${QAPA_FASTA_DIR}/${BED_BASENAME}"
GENOME_FASTA_FIXED="${QAPA_FASTA_DIR}/${GENOME_FASTA_BASENAME}"
QAPA_TRANSCRIPTOME_FASTA="${QAPA_FASTA_DIR}/qapa.gencode_${SPECIES}_${ASSEMBLY_NAME}_${VERSION_PREFIX}.fa"
BED_GZ_FILE="${BED_FILE}.gz"
GENOME_GZ_FILE="${QAPA_FASTA_DIR}/${FTP_GENOME_FASTA}"
ANNOTATION_DB_FILE="${DB_DIR}/ensembl_identifiers.${ASSEMBLY_NAME}.txt"

if [ ! -f "$ANNOTATION_DB_FILE" ]; then 
    echo "INFO:   - Downloading Ensembl Annotation..."
    wget -qO "$ANNOTATION_DB_FILE" "$BIOMART_URL"
else 
    echo "INFO:   - Annotation DB already exists."
fi

if [ ! -f "$BED_FILE" ]; then
    if [ ! -f "$BED_GZ_FILE" ]; then 
        echo "INFO:   - Downloading 3'UTR BED..."
        wget -q -c -O "$BED_GZ_FILE" "$BED_URL"
    fi
    gunzip -c "$BED_GZ_FILE" > "$BED_FILE"
else 
    echo "INFO:   - 3'UTR BED already exists."
fi

if [ ! -f "$GENOME_FASTA_FIXED" ]; then
    if [ ! -f "$GENOME_GZ_FILE" ]; then 
        echo "INFO:   - Downloading Genome FASTA..."
        wget -q -c -O "$GENOME_GZ_FILE" "${FTP_BASE}/${FTP_GENOME_FASTA}"
    fi
    echo "INFO:   - Fixing FASTA headers (adding 'chr')..."
    gunzip -c "$GENOME_GZ_FILE" | awk '/^>/ {if (!/^>chr/) {print ">chr" substr($1,2)} else {print $0}} !/^>/ {print $0}' > "$GENOME_FASTA_FIXED"
else 
    echo "INFO:   - Fixed Genome FASTA already exists."
fi

# --- 6. Build QAPA Transcriptome Reference ---
echo "----------------------------------------------------"
echo "INFO: Step 2/9 - Generating QAPA Transcriptome FASTA..."
if [ ! -f "$QAPA_TRANSCRIPTOME_FASTA" ]; then
    echo "INFO:   - Running 'qapa fasta'..."
    if [ "$EXEC_MODE" == "docker" ]; then
        docker run --rm -u $(id -u):$(id -g) -v "$QAPA_FASTA_DIR":/data --entrypoint qapa "$QAPA_CONTAINER" fasta -f "/data/$(basename "$GENOME_FASTA_FIXED")" "/data/$(basename "$BED_FILE")" "/data/$(basename "$QAPA_TRANSCRIPTOME_FASTA")"
    else
        qapa fasta -f "$GENOME_FASTA_FIXED" "$BED_FILE" "$QAPA_TRANSCRIPTOME_FASTA"
    fi
else
    echo "INFO:   - QAPA Transcriptome FASTA already exists."
fi

# --- 7. Build Salmon Index ---
echo "----------------------------------------------------"
echo "INFO: Step 3/9 - Building Salmon Index..."
if [ ! -d "$SALMON_INDEX_DIR" ]; then
    echo "INFO:   - Running 'salmon index'..."
    if [ "$EXEC_MODE" == "docker" ]; then
        docker run --rm -u $(id -u):$(id -g) -v "$QAPA_FASTA_DIR":/data/refs "$SALMON_CONTAINER" salmon index --threads "$THREADS" -t "/data/refs/$(basename "$QAPA_TRANSCRIPTOME_FASTA")" -i "/data/refs/$(basename "$SALMON_INDEX_DIR")"
    else
        salmon index --threads "$THREADS" -t "$QAPA_TRANSCRIPTOME_FASTA" -i "$SALMON_INDEX_DIR"
    fi
else
    echo "INFO:   - Salmon index already exists."
fi

# --- 8. Salmon Quantification ---
echo "----------------------------------------------------"
echo "INFO: Step 4/9 - Running Salmon Quantification..."
# [FIXED]: Added 'rest' to catch the 4th column (subject) so it doesn't get mixed into 'condition'
tail -n +2 "$SAMPLES_FILE_ABS" | tr -d '\r' | while IFS=, read -r file_name sample_name condition rest; do
    echo "INFO:   -> Processing: $sample_name ($file_name)"
    if [ -f "${SALMON_OUT_DIR}/${file_name}/quant.sf" ]; then echo "INFO:      - Skipped (Already exists)."; continue; fi
    
    read_files=($(find "$FASTQ_DIR" -maxdepth 1 \( -name "*${file_name}*trimmed.f*q.gz" -o -name "*${file_name}*_val_[12].f*q.gz" \) | sort))
    if [ ${#read_files[@]} -eq 0 ]; then
        read_files=($(find "$FASTQ_DIR" -maxdepth 1 -name "*${file_name}*.f*q.gz" ! -name "*trimmed.f*q.gz" ! -name "*_val_[12].f*q.gz" | sort))
    fi
    
    num_files=${#read_files[@]}
    if [ "$num_files" -eq 1 ]; then
        read_opts_docker="-r /data/fastq/$(basename "${read_files[0]}")"; read_opts_local="-r ${read_files[0]}"
    elif [ "$num_files" -eq 2 ]; then
        read1=$(echo "${read_files[@]}" | tr ' ' '\n' | grep -E "(_R1|_1|\.1\.)" || echo "${read_files[0]}")
        read2=$(echo "${read_files[@]}" | tr ' ' '\n' | grep -E "(_R2|_2|\.2\.)" || echo "${read_files[1]}")
        read_opts_docker="-1 /data/fastq/$(basename "$read1") -2 /data/fastq/$(basename "$read2")"
        read_opts_local="-1 $read1 -2 $read2"
    else
        echo "ERROR: Found $num_files FASTQ files for $file_name. Expected 1 or 2." >&2; exit 1
    fi

    if [ "$EXEC_MODE" == "docker" ]; then
        docker run --rm -u $(id -u):$(id -g) -v "$FASTQ_DIR":/data/fastq -v "$SALMON_INDEX_DIR":/data/index -v "$SALMON_OUT_DIR":/data/output "$SALMON_CONTAINER" salmon quant -i /data/index -l A ${read_opts_docker} -p "$THREADS" --validateMappings --gcBias -o "/data/output/${file_name}" > /dev/null 2>&1
    else
        salmon quant -i "$SALMON_INDEX_DIR" -l A ${read_opts_local} -p "$THREADS" --validateMappings --gcBias -o "${SALMON_OUT_DIR}/${file_name}" > /dev/null 2>&1
    fi
done

# --- 9. QAPA Quantification (PAU Calculation) ---
echo "----------------------------------------------------"
echo "INFO: Step 5/9 - Aggregating PAU with 'qapa quant'..."
PAU_RESULTS_FILE="${SALMON_OUT_DIR}/${PREFIX}_pau_results.txt"

if [ ! -f "$PAU_RESULTS_FILE" ]; then
    QUANT_FILES=$(find "${SALMON_OUT_DIR}" -mindepth 2 -name "quant.sf" | sort)
    if [ -z "$QUANT_FILES" ]; then echo "ERROR: No quant.sf files found." >&2; exit 1; fi
    
    echo "INFO:   - Running 'qapa quant'..."
    if [ "$EXEC_MODE" == "docker" ]; then
        DB_MNT="/data/db"; QUANT_MNT="/data/quant"
        QUANT_FILES_DOCKER=$(echo "$QUANT_FILES" | sed "s|${SALMON_OUT_DIR}|$QUANT_MNT|")
        docker run --rm -u $(id -u):$(id -g) -v "$DB_DIR":${DB_MNT} -v "$SALMON_OUT_DIR":${QUANT_MNT} --entrypoint qapa "$QAPA_CONTAINER" quant --db "${DB_MNT}/$(basename "$ANNOTATION_DB_FILE")" ${QUANT_FILES_DOCKER} > "$PAU_RESULTS_FILE"
    else
        qapa quant --db "$ANNOTATION_DB_FILE" ${QUANT_FILES} > "$PAU_RESULTS_FILE"
    fi
else
    echo "INFO:   - PAU results file already exists."
fi

# --- 10. Generate Mapping Files ---
echo "----------------------------------------------------"
echo "INFO: Step 6/9 - Creating ID Mapping Files..."
TX2APA_ID_FILE="${QAPA_FASTA_DIR}/3utrtx2apaID.tsv"
TX2GENE_FILE="${QAPA_FASTA_DIR}/3utrtx2apaiso.tsv"
TX2SYMBOL_FILE="${QAPA_FASTA_DIR}/apaiso2genesymbol.tsv"

if [ ! -f "$TX2APA_ID_FILE" ] || [ ! -f "$TX2GENE_FILE" ] || [ ! -f "$TX2SYMBOL_FILE" ]; then
    awk -F'\t' \
        -v tx2apaID_file="$TX2APA_ID_FILE" \
        -v tx2gene_file="$TX2GENE_FILE" \
        -v tx2symbol_file="$TX2SYMBOL_FILE" \
        '
        BEGIN { OFS="\t"; print "TXNAME", "GENEID" > tx2apaID_file; print "TXNAME", "GENEID" > tx2gene_file; print "GENEID", "SYMBOL" > tx2symbol_file }
        FILENAME == ARGV[1] { if (FNR>1) { n=split($2, t, ","); for(i=1; i<=n; i++) k[t[i] "_" $9 "_" $10]=$1 }; next }
        FILENAME == ARGV[2] { if (FNR>1) { if($1!="") { g[$1]=$2; s[$1]=$3 } }; next }
        FILENAME == ARGV[3] {
            if (/^>/) {
                hdr = substr($0, 2)
                if (match(hdr, /^(ENS[A-Z]*T[0-9]+).*_utr_([0-9]+)_([0-9]+)\(.\)$/, m)) {
                    key = m[1] "_" m[2] "_" m[3]; apa = k[key]; gene = g[m[1]]; sym = s[m[1]]
                    if (apa != "") { print hdr, apa > tx2apaID_file; if(sym!="") a2s[apa]=sym }
                    if (gene != "") print hdr, gene > tx2gene_file
                }
            }
        }
        END { for(i in a2s) print i, a2s[i] > tx2symbol_file; close(tx2apaID_file); close(tx2gene_file); close(tx2symbol_file) }
        ' "$PAU_RESULTS_FILE" "$ANNOTATION_DB_FILE" "$QAPA_TRANSCRIPTOME_FASTA"
    echo "INFO:   - Mapping files generated."
else
    echo "INFO:   - Mapping files already exist."
fi

# --- 11. Generate Individual Isoform Count Matrix (tximport) ---
echo "----------------------------------------------------"
echo "INFO: Step 7/9 - Importing Counts (Individual Isoforms)..."
COUNT_MATRIX_INDIVIDUAL_ISOFORM_FILE="${SALMON_OUT_DIR}/qapa_individual_isoform_counts_matrix.tsv"

if [ ! -f "$COUNT_MATRIX_INDIVIDUAL_ISOFORM_FILE" ]; then
    R_SCRIPT_CONTENT=$(cat <<'EOF'
suppressPackageStartupMessages({ library(tximport); library(readr); library(dplyr); library(tibble) })
args <- commandArgs(trailingOnly = TRUE)
samples <- read.csv(args[1], header=TRUE, stringsAsFactors=FALSE)
files <- file.path(args[2], samples$file_name, "quant.sf")
names(files) <- samples$file_name
tx2gene <- read_tsv(args[3], col_types = "cc")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE)
counts_df <- as.data.frame(round(txi$counts)) %>% rownames_to_column(var = "transcript_id")
write_tsv(counts_df, args[4])
EOF
    )

    if [ "$EXEC_MODE" == "docker" ]; then
        docker run --rm -i -u $(id -u):$(id -g) -v "$SALMON_OUT_DIR":/data/quant -v "$QAPA_FASTA_DIR":/data/fasta -v "$SAMPLES_FILE_ABS":/data/s.csv "$TXIMPORT_CONTAINER" Rscript - "/data/s.csv" "/data/quant" "/data/fasta/$(basename "$TX2GENE_FILE")" "/data/quant/$(basename "$COUNT_MATRIX_INDIVIDUAL_ISOFORM_FILE")" <<< "$R_SCRIPT_CONTENT"
    else
        Rscript - "$SAMPLES_FILE_ABS" "$SALMON_OUT_DIR" "$TX2GENE_FILE" "$COUNT_MATRIX_INDIVIDUAL_ISOFORM_FILE" <<< "$R_SCRIPT_CONTENT"
    fi
    echo "INFO:   - Individual isoform matrix generated."
else
    echo "INFO:   - Individual isoform matrix already exists."
fi

# --- 12. P/D Aggregation for DRIMSeq ---
echo "----------------------------------------------------"
echo "INFO: Step 8/9 - Aggregating P/D Isoforms for DTU Analysis..."
COUNT_MATRIX_APA_AGGREGATED_FILE="${SALMON_OUT_DIR}/qapa_transcript_counts_matrix.tsv"
APAID_TO_GENEID_FILE="${QAPA_FASTA_DIR}/apaid2geneid.tsv"
GENEID_TO_SYMBOL_FILE="${QAPA_FASTA_DIR}/geneid2genesymbol.tsv"

if [ ! -f "$COUNT_MATRIX_APA_AGGREGATED_FILE" ]; then
    echo "INFO:   - Generating DRIMSeq-ready files..."

    awk -F'\t' 'BEGIN{OFS="\t"}
        NR==FNR { if(FNR>1) mapA[$1] = $2; next }
        FNR>1 { mapB[$1] = $2 }
        END {
            print "TXNAME", "GENEID"
            for (tx in mapA) {
                apa = mapA[tx]; geneID = mapB[tx]
                split(apa, p, "_"); type = p[length(p)]
                g = p[1]; for(i=2; i<length(p)-1; i++) g = g "_" p[i]
                if ((type=="P"||type=="D") && geneID!="") {
                    k = g "_" type; out[k] = geneID
                }
            }
            for (k in out) print k, out[k]
        }
    ' "$TX2APA_ID_FILE" "$TX2GENE_FILE" > "$APAID_TO_GENEID_FILE"

    awk -F'\t' 'BEGIN{OFS="\t"} FNR>1{ if($3!="" && $4!="") m[$3]=$4 } END{ print "GENEID","SYMBOL"; for(i in m) print i, m[i] }' "$PAU_RESULTS_FILE" > "$GENEID_TO_SYMBOL_FILE"

    awk -F'\t' 'BEGIN{OFS="\t"}
        FILENAME==ARGV[1] { if(FNR>1) mapA[$1]=$2; next }
        FILENAME==ARGV[2] { if(FNR>1) valid[$1]=1; next }
        FILENAME==ARGV[3] {
            if(FNR==1) { header=$0; NF_h=NF; next }
            apa = mapA[$1]
            if (apa == "") next
            
            split(apa, p, "_"); type = p[length(p)]
            g = p[1]; for(i=2; i<length(p)-1; i++) g = g "_" p[i]
            
            if (type == "P" || type == "D") {
                aggID = g "_" type
                if (aggID in valid) {
                    for(i=2; i<=NF; i++) counts[aggID, i] += $i
                }
            }
        }
        END {
            print header
            for (k in counts) { split(k, z, SUBSEP); ids[z[1]]=1 }
            n = asorti(ids, sorted)
            for (i=1; i<=n; i++) {
                id = sorted[i]; printf "%s", id
                for (j=2; j<=NF_h; j++) printf "%s%s", OFS, (counts[id, j]+0)
                printf "\n"
            }
        }
    ' "$TX2APA_ID_FILE" "$APAID_TO_GENEID_FILE" "$COUNT_MATRIX_INDIVIDUAL_ISOFORM_FILE" > "$COUNT_MATRIX_APA_AGGREGATED_FILE"

    echo "INFO:   - P/D Aggregation Complete."
else
    echo "INFO:   - Aggregated count matrix already exists."
fi

# --- 13. Cleanup ---
echo "----------------------------------------------------"
echo "INFO: Step 9/9 - Cleaning Intermediate Files..."
rm -f "$TX2GENE_FILE" "$TX2SYMBOL_FILE"
echo "INFO:   - Done."

# --- Final Report ---
echo "========================================================================"
echo "INFO: QAPA Pipeline Finished Successfully."
echo "INFO: Main Outputs:"
echo "INFO:   1. DRIMSeq Count Matrix (P/D Aggregated):"
echo "INFO:      -> ${COUNT_MATRIX_APA_AGGREGATED_FILE}"
echo "INFO:   2. Gene Mapping File:"
echo "INFO:      -> ${APAID_TO_GENEID_FILE}"
echo "INFO:   3. Symbol Mapping File:"
echo "INFO:      -> ${GENEID_TO_SYMBOL_FILE}"
echo "========================================================================"
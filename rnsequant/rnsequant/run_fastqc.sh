#!/bin/bash
set -e

# ==============================================================================
# Script Name: run_fastqc.sh
# Version: 1.0
# Description: Automated QC and Trimming pipeline for NGS data.
#              Performs FastQC (pre-trim), Trim Galore!, FastQC (post-trim),
#              and aggregates results with MultiQC.
# Supports:    Docker (default) and Local execution modes.
# ==============================================================================

# --- Configuration ---
# Docker images to be used
TRIMGALORE_IMG="quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_1"
FASTQC_IMG="quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
MULTIQC_IMG="quay.io/biocontainers/multiqc:1.10.1--py_0"

# --- Help Function ---
show_help() {
    echo "Usage: $0 [options] <sample_sheet.csv>"
    echo ""
    echo "Description:"
    echo "  This script performs QC (FastQC), trimming (Trim Galore!), and summary reporting (MultiQC)"
    echo "  for a list of FASTQ files specified in a CSV file."
    echo ""
    echo "Arguments:"
    echo "  <sample_sheet.csv>    Mandatory CSV file. Header must be exactly one of:"
    echo "                          Unpaired: file_name,sample_name,condition"
    echo "                          Paired:   file_name,sample_name,condition,subject"
    echo ""
    echo "Options:"
    echo "  -h, --help                Show this help message and exit."
    echo "  -v, --version             Show versions of both local and Docker tools and exit."
    echo "  -i, --input-fastq-dir <path>"
    echo "                            Path to the directory containing input FASTQ files. (Default: current directory)"
    echo "  -o, --output-dir <path>"
    echo "                            Path to the directory for trimmed FASTQ and QC results. (Default: ./trimed_fastq)"
    echo "  --local                   Use locally installed tools instead of Docker containers."
    echo ""
}

# --- Default Variables ---
RUN_MODE="Docker"
CSV_FILE=""
INPUT_DIR="."
OUTPUT_DIR="trimed_fastq"

# --- Argument Parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      show_help
      exit 0
      ;;
    -v|--version)
      echo "--- Docker Image Versions ---"
      docker run --rm "${FASTQC_IMG}" fastqc --version
      docker run --rm "${TRIMGALORE_IMG}" trim_galore --version
      docker run --rm "${MULTIQC_IMG}" multiqc --version
      echo -e "\n--- Local Tool Versions ---"
      command -v fastqc &> /dev/null && fastqc --version || echo "FastQC (Local): Not installed"
      command -v trim_galore &> /dev/null && trim_galore --version || echo "Trim Galore! (Local): Not installed"
      command -v cutadapt &> /dev/null && cutadapt --version || echo "Cutadapt (Local): Not installed"
      command -v multiqc &> /dev/null && multiqc --version || echo "MultiQC (Local): Not installed"
      exit 0
      ;;
    -i|--input-fastq-dir)
      INPUT_DIR="$2"
      shift 2
      ;;
    -o|--output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --local)
      RUN_MODE="Local"
      shift
      ;;
    *)
      if [ -z "${CSV_FILE}" ]; then
        CSV_FILE="$1"
      else
        echo "ERROR: Too many arguments. Please provide only one CSV file." >&2; show_help >&2; exit 1
      fi
      shift
      ;;
  esac
done

# --- Pre-flight Checks ---
if [ -z "${CSV_FILE}" ]; then echo "ERROR: Missing CSV file argument." >&2; show_help >&2; exit 1; fi
if [ ! -r "${CSV_FILE}" ]; then echo "ERROR: Cannot read CSV file at ${CSV_FILE}"; exit 1; fi
if [ ! -d "${INPUT_DIR}" ]; then echo "ERROR: Input FASTQ directory not found at ${INPUT_DIR}"; exit 1; fi

# --- Prepare Output Directory ---
echo "INFO: Output for trimmed files will be in '${OUTPUT_DIR}/'"
mkdir -p "${OUTPUT_DIR}"

echo "INFO: Running in ${RUN_MODE} mode."
echo "INFO: Looking for input FASTQ files in: ${INPUT_DIR}"

# --- Strict Header Validation (Paired/Unpaired Auto-Detect) ---
EXPECTED_HEADER_UNPAIRED="file_name,sample_name,condition"
EXPECTED_HEADER_PAIRED="file_name,sample_name,condition,subject"
# tr -d '\r' is used to handle potential Windows-style line endings (CRLF)
ACTUAL_HEADER=$(head -1 "${CSV_FILE}" | tr -d '\r')

if [[ "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_UNPAIRED}" && "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_PAIRED}" ]]; then
  echo "ERROR: Invalid CSV header format."
  echo "       Expected:"
  echo "         Unpaired: ${EXPECTED_HEADER_UNPAIRED}"
  echo "         Paired:   ${EXPECTED_HEADER_PAIRED}"
  echo "       Actual:   ${ACTUAL_HEADER}"
  exit 1
fi
echo "INFO: CSV header validated: ${ACTUAL_HEADER}"

# Extract IDs from the 1st column (skip header)
ID_LIST=$(tail -n +2 "${CSV_FILE}" | tr -d '\r' | cut -d, -f1)


# --- Step 0: Verify FASTQ Existence ---
echo "----------------------------------------------------"
echo "INFO: Step 0/5 - Pre-flight check: Verifying FASTQ files..."
for RUN_ID in ${ID_LIST}; do
  if [ -f "${INPUT_DIR}/${RUN_ID}_1.fastq.gz" ]; then
    if [ ! -f "${INPUT_DIR}/${RUN_ID}_2.fastq.gz" ]; then
      echo "ERROR: Read 2 file is missing for PE read ${RUN_ID} (${INPUT_DIR}/${RUN_ID}_2.fastq.gz)"
      exit 1
    fi
  elif [ ! -f "${INPUT_DIR}/${RUN_ID}.fastq.gz" ]; then
    echo "ERROR: FASTQ file not found for SE read ${RUN_ID} (${INPUT_DIR}/${RUN_ID}.fastq.gz or ${INPUT_DIR}/${RUN_ID}_1.fastq.gz)"
    exit 1
  fi
done
echo "INFO: Pre-flight check PASSED. All required FASTQ files are present."

# --- Main Processing Loop ---
for RUN_ID in ${ID_LIST}; do
  echo "===================================================="
  echo "INFO: Processing Run ID: ${RUN_ID}"
  
  IS_PE=false
  if [ -f "${INPUT_DIR}/${RUN_ID}_1.fastq.gz" ]; then IS_PE=true; fi

  # --- Step 1: Pre-trim FastQC ---
  echo "----------------------------------------------------"
  echo "INFO: Step 1/5 - Running FastQC (pre-trim) for ${RUN_ID}..."
  
  if ${IS_PE}; then 
    FILES_TO_QC=("${RUN_ID}_1.fastq.gz" "${RUN_ID}_2.fastq.gz")
  else 
    FILES_TO_QC=("${RUN_ID}.fastq.gz")
  fi

  for fq_file in "${FILES_TO_QC[@]}"; do
    base_name=$(basename "${fq_file}" .gz | sed -e 's/\.fastq$//' -e 's/\.fq$//')
    expected_output="${INPUT_DIR}/${base_name}_fastqc.html"
    
    if [ ! -f "${expected_output}" ]; then
      echo "  -> Running FastQC on ${fq_file}";
      if [ "${RUN_MODE}" == "Local" ]; then
        fastqc --outdir "${INPUT_DIR}" "${INPUT_DIR}/${fq_file}";
      else
        docker run --rm -u "$(id -u):$(id -g)" -v "${INPUT_DIR}":/data "${FASTQC_IMG}" fastqc --outdir /data "/data/${fq_file}";
      fi
    else 
      echo "  -> Skipping FastQC for ${fq_file} (cached)"
    fi
  done
  
  # --- Step 2: Trim Galore! ---
  echo "----------------------------------------------------"
  echo "INFO: Step 2/5 - Running Trim Galore! for ${RUN_ID}..."
  
  if ${IS_PE}; then
    expected_output="${OUTPUT_DIR}/${RUN_ID}_1_val_1.fq.gz"
    if [ ! -f "${expected_output}" ]; then
      echo "  -> Running Trim Galore! on pair: ${RUN_ID}"
      if [ "${RUN_MODE}" == "Local" ]; then
        trim_galore --paired --cores 4 -o "${OUTPUT_DIR}" "${INPUT_DIR}/${RUN_ID}_1.fastq.gz" "${INPUT_DIR}/${RUN_ID}_2.fastq.gz"
      else
        docker run --rm -u "$(id -u):$(id -g)" -v "${INPUT_DIR}":/input -v "$(pwd)/${OUTPUT_DIR}":/output "${TRIMGALORE_IMG}" \
          trim_galore --paired --cores 4 -o /output "/input/${RUN_ID}_1.fastq.gz" "/input/${RUN_ID}_2.fastq.gz"
      fi
    else 
      echo "  -> Skipping Trim Galore! for pair ${RUN_ID} (cached)"
    fi
  else # SE
    expected_output="${OUTPUT_DIR}/${RUN_ID}_trimmed.fq.gz"
    if [ ! -f "${expected_output}" ]; then
      echo "  -> Running Trim Galore! on single-end: ${RUN_ID}"
      if [ "${RUN_MODE}" == "Local" ]; then
        trim_galore --cores 4 -o "${OUTPUT_DIR}" "${INPUT_DIR}/${RUN_ID}.fastq.gz"
      else
        docker run --rm -u "$(id -u):$(id -g)" -v "${INPUT_DIR}":/input -v "$(pwd)/${OUTPUT_DIR}":/output "${TRIMGALORE_IMG}" \
          trim_galore --cores 4 -o /output "/input/${RUN_ID}.fastq.gz"
      fi
    else 
      echo "  -> Skipping Trim Galore! for single-end ${RUN_ID} (cached)"
    fi
  fi
  
  # --- Step 3: Post-trim FastQC ---
  echo "----------------------------------------------------"
  echo "INFO: Step 3/5 - Running FastQC (post-trim) for ${RUN_ID}..."
  
  if ${IS_PE}; then 
    FILES_TO_QC=("${RUN_ID}_1_val_1.fq.gz" "${RUN_ID}_2_val_2.fq.gz")
  else 
    FILES_TO_QC=("${RUN_ID}_trimmed.fq.gz")
  fi

  for fq_file in "${FILES_TO_QC[@]}"; do
    [ ! -e "${OUTPUT_DIR}/${fq_file}" ] && continue
    
    base_name=$(basename "${fq_file}" .gz | sed -e 's/\.fastq$//' -e 's/\.fq$//')
    expected_output="${OUTPUT_DIR}/${base_name}_fastqc.html"
    
    if [ ! -f "${expected_output}" ]; then
      echo "  -> Running FastQC on ${fq_file}"
      if [ "${RUN_MODE}" == "Local" ]; then
        fastqc --outdir "${OUTPUT_DIR}" "${OUTPUT_DIR}/${fq_file}";
      else
        docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)/${OUTPUT_DIR}":/data "${FASTQC_IMG}" fastqc --outdir /data "/data/${fq_file}";
      fi
    else 
      echo "  -> Skipping FastQC for ${fq_file} (cached)"
    fi
  done
done


# --- Step 4 & 5: Completion Check & MultiQC ---
echo "===================================================="
echo "INFO: Step 4/5 - Final check before running MultiQC..."

ALL_DONE=true
for RUN_ID in ${ID_LIST}; do
  if [ -f "${INPUT_DIR}/${RUN_ID}_1.fastq.gz" ]; then # PE
    if [ ! -f "${OUTPUT_DIR}/${RUN_ID}_1_val_1.fq.gz" ] || [ ! -f "${OUTPUT_DIR}/${RUN_ID}_1_val_1_fastqc.html" ]; then 
      echo "  - WARNING: Post-trim output missing for ${RUN_ID}"; ALL_DONE=false
    fi
  else # SE
    if [ ! -f "${OUTPUT_DIR}/${RUN_ID}_trimmed.fq.gz" ] || [ ! -f "${OUTPUT_DIR}/${RUN_ID}_trimmed_fastqc.html" ]; then 
      echo "  - WARNING: Post-trim output missing for ${RUN_ID}"; ALL_DONE=false
    fi
  fi
done

if ${ALL_DONE}; then
  echo "----------------------------------------------------"
  echo "INFO: Step 5/5 - Generating final report with MultiQC..."
  
  FINAL_REPORT="${OUTPUT_DIR}/multiqc_report.html"
  
  if [ -f "${FINAL_REPORT}" ]; then
    echo "  -> Skipping MultiQC. Report already exists at ${FINAL_REPORT} (cached)"
    echo "----------------------------------------------------"
    echo "INFO: All steps completed successfully."
  else
    if [ "${RUN_MODE}" == "Local" ]; then
      multiqc "${INPUT_DIR}" "${OUTPUT_DIR}" --force -o "${OUTPUT_DIR}"
    else
      docker run --rm -u "$(id -u):$(id -g)" \
        -v "${INPUT_DIR}":/input \
        -v "$(pwd)/${OUTPUT_DIR}":/output \
        "${MULTIQC_IMG}" multiqc /input /output --force --outdir /output
    fi
    
    echo "----------------------------------------------------"
    echo "INFO: All steps completed successfully."
    echo "Final report: ${FINAL_REPORT}"
  fi
else
  echo "----------------------------------------------------"
  echo "INFO: Step 5/5 - MultiQC skipped because some files were not fully processed."
  exit 1
fi
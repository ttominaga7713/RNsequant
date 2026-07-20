#!/bin/bash

# ==============================================================================
# Script Name: run_fastq_download.sh
# Version: 1.0
# Description: Download FASTQ files from NCBI SRA (default) or ENA.
#              - Auto-accepts Paired CSV design if 'subject' column is present.
#              - Includes advanced rollback, verification, and retry logic.
#              - Safely preserves 'subject' column during FASTQ merging.
# Supports:    Docker (default) and Local execution modes.
# ==============================================================================

# --- Configuration ---
SRA_TOOLKIT_IMAGE="ncbi/sra-tools:3.2.1"
PIGZ_IMAGE="genevera/docker-pigz"
OUTPUT_DIR="."
DEFAULT_THREADS=8
DEFAULT_SOURCE="ncbi"
MAX_RETRIES=100
TIMEOUT_SEC=3        # ENA Wget timeout
ROLLBACK_MB=5        # ENA Rollback size

# --- Signal Handling ---
function handle_interrupt {
    echo -e "\n\nERROR: [!] Script interrupted by Ctrl+C." >&2
    rm -f wget_error.log 
    if [ -n "$TEMP_CSV_FILE" ] && [ -f "$TEMP_CSV_FILE" ]; then
        rm -f "$TEMP_CSV_FILE"
    fi
    exit 1
}
trap handle_interrupt SIGINT

# --- Help Function ---
show_help() {
    echo "Usage: $(basename "$0") [options] <sample_sheet.csv>"
    echo ""
    echo "Description:"
    echo "  Downloads FASTQ files from ENA or NCBI SRA based on a CSV file."
    echo ""
    echo "Arguments:"
    echo "  <sample_sheet.csv>    Mandatory CSV file. Header must be exactly one of:"
    echo "                          Unpaired: file_name,sample_name,condition"
    echo "                          Paired:   file_name,sample_name,condition,subject"
    echo ""
    echo "Options:"
    echo "  -h, --help                Show this help message and exit."
    echo "  -v, --version             Show tool versions and exit."
    echo "  -o, --output-dir <path>   Directory to save downloaded FASTQ files. (Default: .)"
    echo "  -t, --threads <int>       Number of threads for fasterq-dump/pigz. (Default: 8)"
    echo "  -s, --source <ena|ncbi>   Data source to use. (Default: ncbi)"
    echo "  --local                   Use locally installed tools instead of Docker containers."
    echo ""
}

# --- Argument Parsing ---
USE_DOCKER=true
CSV_FILE=""
THREADS=$DEFAULT_THREADS
SOURCE=$DEFAULT_SOURCE

while (( "$#" )); do
  case "$1" in
    -h|--help) show_help; exit 0 ;;
    -v|--version)
      echo "--- Docker Image Versions ---"
      docker run --rm "${SRA_TOOLKIT_IMAGE}" fasterq-dump --version
      docker run --rm "${PIGZ_IMAGE}" pigz --version
      echo -e "\n--- Local Tool Versions ---"
      command -v fasterq-dump &> /dev/null && fasterq-dump --version || echo "fasterq-dump (Local): Not installed"
      command -v pigz &> /dev/null && pigz --version || echo "pigz (Local): Not installed"
      exit 0
      ;;
    --local) USE_DOCKER=false; shift ;;
    -o|--output-dir) OUTPUT_DIR=$2; shift 2 ;;
    -t|--threads) THREADS=$2; shift 2 ;;
    -s|--source)
      if [[ "$2" =~ ^(ncbi|ena)$ ]]; then SOURCE=$2; shift 2; else echo "ERROR: Invalid source." >&2; exit 1; fi ;;
    -*) echo "ERROR: Unknown option '$1'" >&2; show_help >&2; exit 1 ;;
    *) 
      if [ -z "${CSV_FILE}" ]; then CSV_FILE="$1"; else echo "ERROR: Too many arguments." >&2; show_help >&2; exit 1; fi
      shift ;;
  esac
done

# --- Pre-flight Checks ---
if [ -z "$CSV_FILE" ]; then echo "ERROR: Missing CSV file argument." >&2; show_help >&2; exit 1; fi
if [ ! -r "$CSV_FILE" ]; then echo "ERROR: Cannot read CSV file at ${CSV_FILE}" >&2; exit 1; fi

# --- Path Initialization ---
if [ ! -d "$OUTPUT_DIR" ]; then mkdir -p "$OUTPUT_DIR"; fi
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)

echo "========================================================================"
echo "INFO: FastQ Download Pipeline Started (v1.0)"
echo "INFO: Target CSV: ${CSV_FILE}"
echo "INFO: Source: ${SOURCE} | Threads: ${THREADS} | Docker: ${USE_DOCKER}"
echo "INFO: Output directory: ${OUTPUT_DIR}/"
echo "========================================================================"

# --- Strict Header Validation (Paired/Unpaired Auto-Detect) ---
EXPECTED_HEADER_UNPAIRED="file_name,sample_name,condition"
EXPECTED_HEADER_PAIRED="file_name,sample_name,condition,subject"
ACTUAL_HEADER=$(head -1 "${CSV_FILE}" | tr -d '\r')

if [[ "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_UNPAIRED}" && "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_PAIRED}" ]]; then
  echo "ERROR: Invalid CSV header format." >&2
  echo "       Expected:" >&2
  echo "         Unpaired: ${EXPECTED_HEADER_UNPAIRED}" >&2
  echo "         Paired:   ${EXPECTED_HEADER_PAIRED}" >&2
  echo "       Actual:   ${ACTUAL_HEADER}" >&2
  exit 1
fi

IS_PAIRED=false
if [ "${ACTUAL_HEADER}" == "${EXPECTED_HEADER_PAIRED}" ]; then
    IS_PAIRED=true
fi
echo "INFO: CSV header validated: ${ACTUAL_HEADER} (Paired: ${IS_PAIRED})"

# CSV Loading (Extract SRR IDs from first column)
SRR_LIST=($(tail -n +2 "$CSV_FILE" | cut -d ',' -f 1 | tr -d '\r' | sed '/^\s*$/d'))

if [ ${#SRR_LIST[@]} -eq 0 ]; then
    echo "ERROR: No valid SRR IDs found in $CSV_FILE" >&2
    exit 1
fi

# ==============================================================================
# Helper Functions
# ==============================================================================

function compress_fastq_files {
    local target_dir=$1
    local srr_id=$2
    if [ "$USE_DOCKER" = true ]; then
        local files=$(find "$target_dir" -maxdepth 1 -name "${srr_id}*.fastq")
        for file in $files; do
            echo "INFO:   -> Compressing $(basename "$file")..."
            docker run --rm -u "$(id -u):$(id -g)" -v "${target_dir}:/data" -w /data "$PIGZ_IMAGE" pigz -p "$THREADS" "$(basename "$file")"
        done
    else
        local files=("${target_dir}/${srr_id}"*.fastq)
        if [ -e "${files[0]}" ]; then
            for file in "${files[@]}"; do
                echo "INFO:   -> Compressing $(basename "$file")..."
                [ -f "$file" ] && (command -v pigz &> /dev/null && pigz -p "$THREADS" "$file" || gzip "$file")
            done
        fi
    fi
}

function rollback_file {
    local file=$1
    if [ -f "$file" ]; then
        local fsize=$(stat -c%s "$file")
        local cut_bytes=$((ROLLBACK_MB * 1024 * 1024))
        if [ "$fsize" -gt "$cut_bytes" ]; then
            echo "WARNING:  [ROLLBACK] Cutting last ${ROLLBACK_MB}MB..."
            truncate -s -${ROLLBACK_MB}M "$file"
        else
            echo "WARNING:  [ROLLBACK] File too small. Deleting..."
            rm -f "$file"
        fi
    fi
}

function diagnose_wget_error {
    local log_file=$1
    if grep -a -q "unable to resolve host address" "$log_file"; then
        echo "WARNING:  [DIAGNOSIS] DNS Error. Check /etc/resolv.conf"
    elif grep -a -q "404 Not Found" "$log_file"; then
        echo "WARNING:  [DIAGNOSIS] 404 Not Found on ENA."
    elif grep -a -q "No space left on device" "$log_file"; then
        echo "ERROR:    [FATAL] Disk Full!" >&2
        exit 1
    fi
}

# ---------------------------------------------------------
# Download Core Logic
# ---------------------------------------------------------
function perform_download {
    local srr_id=$1
    local filter_pattern=$2
    local outdir=$OUTPUT_DIR
    local log_file="wget_error.log"
    
    if [ "$SOURCE" == "ena" ]; then
        local api_url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr_id}&result=read_run&fields=fastq_ftp&format=tsv&download=true&limit=0"
        
        local ftp_links=$(curl -s "$api_url" | tail -n +2 | cut -f 2)
        if [ -z "$ftp_links" ]; then
            echo "ERROR:    No FTP links found for $srr_id via ENA API." >&2
            return 1 
        fi
        
        IFS=';' read -ra LINKS <<< "$ftp_links"
        for link in "${LINKS[@]}"; do
            if [[ "$link" != http* && "$link" != ftp* ]]; then link="ftp://${link}"; fi
            local filename=$(basename "$link")
            
            if [[ -n "$filter_pattern" && "$filename" != *"$filter_pattern"* ]]; then continue; fi
            [ -f "${outdir}/${filename}.verified" ] && continue

            local attempt=1
            local success=false
            while [ $attempt -le $MAX_RETRIES ]; do
                wget -c -q --show-progress --progress=bar:force:noscroll --read-timeout=$TIMEOUT_SEC --connect-timeout=$TIMEOUT_SEC -t 1 -P "$outdir" "$link" 2>&1 | tee "$log_file"
                if [ ${PIPESTATUS[0]} -eq 0 ]; then
                    success=true; rm -f "$log_file"; break
                else
                    echo "" 
                    diagnose_wget_error "$log_file"
                    rollback_file "${outdir}/${filename}"
                fi
                ((attempt++))
            done
            [ "$success" = false ] && return 1
        done

    elif [ "$SOURCE" == "ncbi" ]; then
        if ls "${outdir}/${srr_id}"*.fastq.gz* 1> /dev/null 2>&1; then
            echo "INFO:   [SKIP] NCBI mode: Files for $srr_id already exist."
            return 0
        fi

        # Disable SRA Toolkit cloud delivery to avoid AWS IMDS errors on non-EC2 hosts.
        # (fasterq-dump >=3.x tries to resolve via AWS by default, which fails locally.)
        local SRA_MKFG_KEYS=(
            '/libs/cloud/report_instance_identity = "false"'
            '/repository/user/main/public/cache-enabled = "false"'
        )
        # -- Docker mode: inject the config file at container startup --
        local MKFG_INJECT=""
        local MKFG_PATH="/root/.ncbi/user-settings.mkfg"
        for key in "${SRA_MKFG_KEYS[@]}"; do
            MKFG_INJECT+="mkdir -p \$(dirname ${MKFG_PATH}); grep -qF '${key%%=*}' ${MKFG_PATH} 2>/dev/null || echo '${key}' >> ${MKFG_PATH}; "
        done

        # -- Local mode: write config to the host's ~/.ncbi --
        function _ensure_local_sra_config {
            local cfg="${HOME}/.ncbi/user-settings.mkfg"
            mkdir -p "$(dirname "$cfg")"
            for key in "${SRA_MKFG_KEYS[@]}"; do
                grep -qF "${key%%=*}" "$cfg" 2>/dev/null || echo "$key" >> "$cfg"
            done
        }

        local attempt=1
        local success=false
        while [ $attempt -le 5 ]; do
            if [ "$USE_DOCKER" = true ]; then
                docker run --rm -u "$(id -u):$(id -g)" -v "${outdir}:/data" -w /data \
                    --entrypoint sh "$SRA_TOOLKIT_IMAGE" \
                    -c "${MKFG_INJECT} fasterq-dump --progress --split-files --threads ${THREADS} ${srr_id}"
            else
                _ensure_local_sra_config
                cd "$outdir" && fasterq-dump --progress --split-files --threads "$THREADS" "$srr_id" && cd - > /dev/null
            fi
            if [ $? -eq 0 ]; then
                compress_fastq_files "$outdir" "$srr_id"
                success=true; break
            fi
            ((attempt++))
        done
        [ "$success" = false ] && return 1
    fi
}

function check_and_mark_verified {
    local file_path=$1
    local verified_path="${file_path}.verified"
    [ -f "$verified_path" ] && return 0
    
    if [ -f "$file_path" ]; then
        if [ "$SOURCE" == "ncbi" ]; then
            echo "INFO:   -> Skipping verify (gzip -t) for NCBI: $(basename "$file_path") ... [OK]"
            mv "$file_path" "$verified_path"
            sync
            return 0
        fi

        echo -n "INFO:   -> Verifying (gzip -t): $(basename "$file_path") ... "
        if gzip -t "$file_path" 2>/dev/null; then
            echo "[OK]"
            mv "$file_path" "$verified_path"
            sync
            return 0
        else
            echo "[FAIL] Corrupted. Deleting file to force redownload."
            rm -f "$file_path"
            return 1
        fi
    fi
    return 1
}

# ---------------------------------------------------------
# Download and Verify Wrapper with Auto-Retry Logic
# ---------------------------------------------------------
function download_and_verify {
    local srr_id=$1
    local filter_pattern=$2
    local expected_file=$3
    local attempt=1
    local max_attempts=3
    
    if [ -f "${expected_file}.verified" ]; then
         echo "INFO:   [SKIP] $(basename "$expected_file") already verified."
         return 0
    fi

    while [ $attempt -le $max_attempts ]; do
        perform_download "$srr_id" "$filter_pattern"
        
        if [ ! -f "$expected_file" ] && [ ! -f "${expected_file}.verified" ]; then
            return 1 
        fi

        if check_and_mark_verified "$expected_file"; then
            return 0
        fi
        
        echo "WARNING: Retrying download for $(basename "$expected_file") (Attempt $attempt/$max_attempts)..."
        ((attempt++))
    done
    
    echo "ERROR: Corrupted file persists for $(basename "$expected_file") after $max_attempts verification attempts." >&2
    return 1
}

# ==============================================================================
# Phase 0 & 1: Fast Track & Layout Detection
# ==============================================================================
echo "----------------------------------------------------"
echo "INFO: Phase 0: Checking completion..."
ALL_COMPLETE=true
for srr_id in "${SRR_LIST[@]}"; do
    if [ -f "${OUTPUT_DIR}/${srr_id}.fastq.gz" ] || ([ -f "${OUTPUT_DIR}/${srr_id}_1.fastq.gz" ] && [ -f "${OUTPUT_DIR}/${srr_id}_2.fastq.gz" ]); then continue; fi
    ALL_COMPLETE=false; break
done

if [ "$ALL_COMPLETE" = true ]; then
    echo "INFO: Basic files are already downloaded."
fi

if [ "$ALL_COMPLETE" = false ]; then
    echo "----------------------------------------------------"
    echo "INFO: Phase 1: Detecting Layout..."
    FIRST_SRR="${SRR_LIST[0]}"

    if [ -f "${OUTPUT_DIR}/${FIRST_SRR}_1.fastq.gz" ] || [ -f "${OUTPUT_DIR}/${FIRST_SRR}_1.fastq.gz.verified" ]; then
        LAYOUT="PE"
    elif [ -f "${OUTPUT_DIR}/${FIRST_SRR}.fastq.gz" ] || [ -f "${OUTPUT_DIR}/${FIRST_SRR}.fastq.gz.verified" ]; then
        LAYOUT="SE"
    else
        echo "INFO:   -> Probing for PE layout (_1)..."
        if download_and_verify "$FIRST_SRR" "_1.fastq.gz" "${OUTPUT_DIR}/${FIRST_SRR}_1.fastq.gz"; then
            echo "INFO:   -> PE Detected. Downloading _2..."
            download_and_verify "$FIRST_SRR" "_2.fastq.gz" "${OUTPUT_DIR}/${FIRST_SRR}_2.fastq.gz" || exit 1
            LAYOUT="PE"
        else
            echo "INFO:   -> _1 not found. Probing for SE layout..."
            if download_and_verify "$FIRST_SRR" "" "${OUTPUT_DIR}/${FIRST_SRR}.fastq.gz"; then
                 LAYOUT="SE"
            else
                 echo "ERROR: Failed to download initial file. ENA API might be down or ID invalid." >&2
                 exit 1
            fi
        fi
    fi
    echo "INFO: Layout detected: $LAYOUT"

    # ==============================================================================
    # Phase 2: Main Loop
    # ==============================================================================
    echo "----------------------------------------------------"
    echo "INFO: Phase 2: Download & Verify Loop..."
    for srr_id in "${SRR_LIST[@]}"; do
        echo "INFO: Processing: $srr_id"
        if [ "$LAYOUT" == "PE" ]; then
            download_and_verify "$srr_id" "_1.fastq.gz" "${OUTPUT_DIR}/${srr_id}_1.fastq.gz" || exit 1
            download_and_verify "$srr_id" "_2.fastq.gz" "${OUTPUT_DIR}/${srr_id}_2.fastq.gz" || exit 1
        else
            download_and_verify "$srr_id" "" "${OUTPUT_DIR}/${srr_id}.fastq.gz" || exit 1
        fi
    done
else
    # Detect Layout from existing files
    FIRST_SRR="${SRR_LIST[0]}"
    if [ -f "${OUTPUT_DIR}/${FIRST_SRR}_1.fastq.gz" ]; then
        LAYOUT="PE"
    else
        LAYOUT="SE"
    fi
    echo "INFO: Layout detected from existing files: $LAYOUT"
fi

# ==============================================================================
# Phase 3: Finalize
# ==============================================================================
echo "----------------------------------------------------"
echo "INFO: Phase 3: Finalize..."
ALL_CLEARED=true
for srr_id in "${SRR_LIST[@]}"; do
    if [ "$LAYOUT" == "PE" ]; then
        ([ ! -f "${OUTPUT_DIR}/${srr_id}_1.fastq.gz.verified" ] && [ ! -f "${OUTPUT_DIR}/${srr_id}_1.fastq.gz" ]) && ALL_CLEARED=false
        ([ ! -f "${OUTPUT_DIR}/${srr_id}_2.fastq.gz.verified" ] && [ ! -f "${OUTPUT_DIR}/${srr_id}_2.fastq.gz" ]) && ALL_CLEARED=false
    else
        ([ ! -f "${OUTPUT_DIR}/${srr_id}.fastq.gz.verified" ] && [ ! -f "${OUTPUT_DIR}/${srr_id}.fastq.gz" ]) && ALL_CLEARED=false
    fi
done

if [ "$ALL_CLEARED" = true ]; then
    find "$OUTPUT_DIR" -name "*.fastq.gz.verified" | while read -r vfile; do mv "$vfile" "${vfile%.verified}"; done
    echo "INFO: Flushing disk cache (sync)..."
    sync

    # ==============================================================================
    # Phase 4: FastQ Merging & In-place CSV Update
    # ==============================================================================
    echo "----------------------------------------------------"
    echo "INFO: Phase 4: Checking if FASTQ merging is required..."
    
    MERGE_REQUIRED=false
    UNIQUE_SAMPLES=$(tail -n +2 "$CSV_FILE" | cut -d ',' -f 2 | tr -d '\r' | sed '/^\s*$/d' | sort -u)

    for sample in $UNIQUE_SAMPLES; do
        SRR_COUNT=$(awk -F',' -v s="$sample" 'NR>1 && $2==s {count++} END {print count}' "$CSV_FILE")
        if [ "$SRR_COUNT" -gt 1 ]; then
            MERGE_REQUIRED=true
            break
        fi
    done

    if [ "$MERGE_REQUIRED" = false ]; then
        echo "INFO: No duplicate sample_names found. Skipping merge phase."
        echo "========================================================================"
        echo "INFO: FastQ Download Pipeline Finished Successfully."
        echo "========================================================================"
        exit 0
    fi

    echo "INFO: Duplicate sample_names detected. Merging FASTQ files..."
    
    ORIGINAL_BACKUP="${CSV_FILE%.csv}_original.csv"
    if [ ! -f "$ORIGINAL_BACKUP" ]; then
        cp "$CSV_FILE" "$ORIGINAL_BACKUP"
    fi

    TEMP_CSV_FILE="${CSV_FILE}.tmp"
    echo "$ACTUAL_HEADER" > "$TEMP_CSV_FILE"

    echo "$UNIQUE_SAMPLES" | while IFS= read -r sample; do
        
        SRRS_FOR_SAMPLE=$(awk -F',' -v s="$sample" 'NR>1 && $2==s {gsub(/\r/,""); print $1}' "$CSV_FILE")
        SRR_COUNT=$(echo "$SRRS_FOR_SAMPLE" | wc -w)
        CONDITION=$(awk -F',' -v s="$sample" 'NR>1 && $2==s {gsub(/\r/,""); print $3; exit}' "$CSV_FILE")
        
        if [ "$IS_PAIRED" = true ]; then
            SUBJECT=$(awk -F',' -v s="$sample" 'NR>1 && $2==s {gsub(/\r/,""); print $4; exit}' "$CSV_FILE")
            ROW_SUFFIX="${sample},${CONDITION},${SUBJECT}"
        else
            ROW_SUFFIX="${sample},${CONDITION}"
        fi
        
        if [ "$SRR_COUNT" -eq 1 ]; then
            echo "INFO:   -> [SKIP] sample '$sample' has only 1 file ($SRRS_FOR_SAMPLE). No merge needed."
            echo "${SRRS_FOR_SAMPLE},${ROW_SUFFIX}" >> "$TEMP_CSV_FILE"
        else
            MERGED_PREFIX=$(echo $SRRS_FOR_SAMPLE | tr ' ' '_')
            
            ALREADY_MERGED=true
            if [ "$LAYOUT" == "PE" ]; then
                [ ! -f "${OUTPUT_DIR}/${MERGED_PREFIX}_1.fastq.gz" ] && ALREADY_MERGED=false
                [ ! -f "${OUTPUT_DIR}/${MERGED_PREFIX}_2.fastq.gz" ] && ALREADY_MERGED=false
            else
                [ ! -f "${OUTPUT_DIR}/${MERGED_PREFIX}.fastq.gz" ] && ALREADY_MERGED=false
            fi
            
            if [ "$ALREADY_MERGED" = true ]; then
                echo "INFO:   -> [SKIP] Merged files for sample '$sample' already exist."
                echo "${MERGED_PREFIX},${ROW_SUFFIX}" >> "$TEMP_CSV_FILE"
                continue
            fi

            echo "INFO:   -> [MERGE] Merging files for sample: $sample"
            echo "INFO:      -> Output prefix: $MERGED_PREFIX"
            
            if [ "$LAYOUT" == "PE" ]; then
                R1_FILES=""
                R2_FILES=""
                for srr in $SRRS_FOR_SAMPLE; do
                    R1_FILES="$R1_FILES ${OUTPUT_DIR}/${srr}_1.fastq.gz"
                    R2_FILES="$R2_FILES ${OUTPUT_DIR}/${srr}_2.fastq.gz"
                done
                
                echo "INFO:      -> Creating ${MERGED_PREFIX}_1.fastq.gz"
                cat $R1_FILES > "${OUTPUT_DIR}/${MERGED_PREFIX}_1.fastq.gz"
                echo "INFO:      -> Verifying integrity of ${MERGED_PREFIX}_1.fastq.gz..."
                if ! gzip -t "${OUTPUT_DIR}/${MERGED_PREFIX}_1.fastq.gz" 2>/dev/null; then
                    echo "ERROR: Merged file corrupted! Deleting ${MERGED_PREFIX}_1.fastq.gz" >&2
                    rm -f "${OUTPUT_DIR}/${MERGED_PREFIX}_1.fastq.gz"
                    rm -f "$TEMP_CSV_FILE"
                    exit 1
                fi
                echo "INFO:         [OK]"

                echo "INFO:      -> Creating ${MERGED_PREFIX}_2.fastq.gz"
                cat $R2_FILES > "${OUTPUT_DIR}/${MERGED_PREFIX}_2.fastq.gz"
                echo "INFO:      -> Verifying integrity of ${MERGED_PREFIX}_2.fastq.gz..."
                if ! gzip -t "${OUTPUT_DIR}/${MERGED_PREFIX}_2.fastq.gz" 2>/dev/null; then
                    echo "ERROR: Merged file corrupted! Deleting ${MERGED_PREFIX}_2.fastq.gz" >&2
                    rm -f "${OUTPUT_DIR}/${MERGED_PREFIX}_2.fastq.gz"
                    rm -f "${OUTPUT_DIR}/${MERGED_PREFIX}_1.fastq.gz" 
                    rm -f "$TEMP_CSV_FILE"
                    exit 1
                fi
                echo "INFO:         [OK]"

            else
                SE_FILES=""
                for srr in $SRRS_FOR_SAMPLE; do
                    SE_FILES="$SE_FILES ${OUTPUT_DIR}/${srr}.fastq.gz"
                done
                
                echo "INFO:      -> Creating ${MERGED_PREFIX}.fastq.gz"
                cat $SE_FILES > "${OUTPUT_DIR}/${MERGED_PREFIX}.fastq.gz"
                echo "INFO:      -> Verifying integrity of ${MERGED_PREFIX}.fastq.gz..."
                if ! gzip -t "${OUTPUT_DIR}/${MERGED_PREFIX}.fastq.gz" 2>/dev/null; then
                    echo "ERROR: Merged file corrupted! Deleting ${MERGED_PREFIX}.fastq.gz" >&2
                    rm -f "${OUTPUT_DIR}/${MERGED_PREFIX}.fastq.gz"
                    rm -f "$TEMP_CSV_FILE"
                    exit 1
                fi
                echo "INFO:         [OK]"
            fi
            
            echo "${MERGED_PREFIX},${ROW_SUFFIX}" >> "$TEMP_CSV_FILE"
        fi
    done

    mv "$TEMP_CSV_FILE" "$CSV_FILE"

    echo "========================================================================"
    echo "INFO: FastQ Download & Merge Pipeline Finished Successfully."
    echo "INFO: Original CSV backed up as: $(basename "$ORIGINAL_BACKUP")"
    echo "INFO: Target CSV ($(basename "$CSV_FILE")) successfully updated."
    echo "========================================================================"
else
    echo "ERROR: Download incomplete. Please run the script again." >&2
    exit 1
fi
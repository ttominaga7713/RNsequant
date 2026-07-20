#!/bin/bash
set -e
set -o pipefail
set -u

# ==============================================================================
# Script Name: run_edgeR_core.sh
# Version: 1.0
# Description: Performs Differential Expression Gene (DEG) analysis using edgeR.
#              - Generates ALL pairwise permutations.
#              - Auto-detects Paired design if 'subject' column is present.
#              - Requires exact CSV header: 
#                Unpaired: 'file_name,sample_name,condition'
#                Paired:   'file_name,sample_name,condition,subject'
# Supports:    Docker (default) and Local execution modes.
# ==============================================================================

# --- Configuration ---
EDGER_CONTAINER="ezojika7713/edger:3.36.0"

# --- Help Message ---
usage() {
  cat <<EOM
Usage: $(basename "$0") [OPTIONS] <samples.csv> <counts.tsv>

Description:
  Performs Differential Expression Gene (DEG) analysis using the edgeR package.
  It automatically creates ALL pairwise permutations from the 'condition' column.
  If a 'subject' column is detected, it automatically runs a Paired analysis.

Required Arguments:
  <samples.csv>    Sample metadata file in CSV format.
                   Header must be exactly one of the following:
                     Unpaired: file_name,sample_name,condition
                     Paired:   file_name,sample_name,condition,subject
                   
  <counts.tsv>     Raw gene count matrix in TSV (tab-separated) format.
                     - First column: Gene IDs
                     - Subsequent columns: Raw read counts for each sample
                     - Column names must match 'sample_name' in samples.csv

Optional Flags:
  --prefix <string>            Output file prefix.
  --gene-symbol-map <map.tsv>  Optional 2-column TSV file to ANNOTATE results.
  --local                      Use locally installed R instead of Docker.
  -h, --help                   Display this help message and exit.

Output:
  1. Pairwise comparison files: {prefix}_{Control}_vs_{Target}_edgeR_results.tsv
  2. ANOVA-like test (if 3+ groups): {prefix}_ANOVA_like_edgeR_results.tsv
EOM
  exit 1
}

# --- 1. Argument Parsing ---
PREFIX=""
GENE_SYMBOL_MAP_FILE=""
EXEC_MODE="docker"
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      PREFIX="$2"; shift 2 ;;
    --gene-symbol-map)
      if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
        echo "ERROR: --gene-symbol-map option requires a file path." >&2; exit 1
      fi
      GENE_SYMBOL_MAP_FILE="$2"; shift 2 ;;
    --local)
      EXEC_MODE="local"; shift ;;
    -h|--help)
      usage ;;
    -*)
      echo "ERROR: Unknown option '$1'" >&2; usage ;;
    *)
      POSITIONAL_ARGS+=("$1"); shift ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}"

if [[ $# -ne 2 ]]; then
  echo "ERROR: Missing required arguments." >&2
  usage
fi

SAMPLES_FILE=$1
COUNTS_FILE=$2

if [[ -z "$PREFIX" ]]; then
  PREFIX=$(basename "$COUNTS_FILE" | sed 's/_counts\.tsv$//; s/_counts\.txt$//')
fi

# Validation
for f in "$SAMPLES_FILE" "$COUNTS_FILE"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Required file not found: $f" >&2; exit 1
  fi
done

if [[ -n "$GENE_SYMBOL_MAP_FILE" && ! -f "$GENE_SYMBOL_MAP_FILE" ]]; then
  echo "ERROR: Gene symbol map file not found: $GENE_SYMBOL_MAP_FILE" >&2; exit 1
fi

if [[ "$EXEC_MODE" == "local" ]] && ! command -v Rscript &> /dev/null; then
  echo "ERROR: In --local mode, 'Rscript' must be installed and in your PATH." >&2
  exit 1
fi

echo "========================================================================"
echo "INFO: edgeR Differential Expression Pipeline Started (v1.0)"
echo "INFO: Mode:              $EXEC_MODE"
echo "INFO: Samples file:      $SAMPLES_FILE"
echo "INFO: Counts file:       $COUNTS_FILE"
echo "INFO: Output prefix:     $PREFIX"
echo "INFO: Naming Convention: {Control}_vs_{Target}"
echo "========================================================================"

# --- Strict Header Validation ---
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

# --- 2. Execute R Pipeline ---
echo "----------------------------------------------------"
echo "INFO: Executing R Analysis Script..."

# Define R script as a heredoc
read -r -d '' R_SCRIPT <<'RSCRIPT_EOF' || true
suppressPackageStartupMessages({
    library(edgeR)
    library(readr)
    library(dplyr)
    library(tibble)
})

# === 1. Parse Arguments ===
args <- commandArgs(trailingOnly = TRUE)
samples_file_path       <- args[1]
counts_file_path        <- args[2]
prefix                  <- args[3]
gene_symbol_map_path    <- args[4]

# === 2. Load Data ===
cat("INFO: --- Reading count and sample data ---\n")
sample_info <- read_csv(samples_file_path, show_col_types = FALSE)
DGE_counts_raw <- read_tsv(counts_file_path, show_col_types = FALSE)
unique_conditions <- unique(sample_info$condition)

if (length(unique_conditions) < 2) {
    stop("ERROR: Cannot generate contrasts. At least two unique conditions are required.")
}

# === 3. Generate Permutations (Control vs Target) ===
pairs <- expand.grid(Control=unique_conditions, Target=unique_conditions, stringsAsFactors=FALSE)
pairs <- pairs[pairs$Control != pairs$Target, ]

contrasts_list <- character()
valid_contrasts_list <- character()
file_suffixes <- character()

for (i in 1:nrow(pairs)) {
    control <- pairs$Control[i]
    target <- pairs$Target[i]
    
    contrast_name <- paste(target, control, sep = "-")
    valid_name <- paste(make.names(target), make.names(control), sep = "-")
    suffix <- paste(control, "vs", target, sep = "_")
    
    contrasts_list <- c(contrasts_list, contrast_name)
    valid_contrasts_list <- c(valid_contrasts_list, valid_name)
    file_suffixes <- c(file_suffixes, suffix)
}

cat(paste("INFO: Generated", length(contrasts_list), "comparisons (Control_vs_Target):\n"))
for(i in seq_along(contrasts_list)) {
    cat(paste0("INFO:   - File: ", file_suffixes[i], "  [Calc: ", contrasts_list[i], "]\n"))
}

# === 4. Load Symbol Map ===
symbol_map <- NULL
if (gene_symbol_map_path != "NULL") {
    cat(paste("INFO: --- Reading gene symbol map for annotation ---\n"))
    map_col_name <- colnames(DGE_counts_raw)[1]
    symbol_map <- read_tsv(gene_symbol_map_path, col_names = c(map_col_name, "GeneSymbol"), show_col_types = FALSE) %>%
      filter(!is.na(!!sym(map_col_name)), !is.na(GeneSymbol)) %>%
      distinct(!!sym(map_col_name), .keep_all = TRUE)
}

# === 5. Create DGEList ===
cat("INFO: --- Creating DGEList object ---\n")
id_col_name <- colnames(DGE_counts_raw)[1]
DGE_counts <- DGE_counts_raw %>%
    as.data.frame() %>%
    column_to_rownames(var = id_col_name)
DGE_counts <- DGE_counts[, sample_info$sample_name]

group <- factor(sample_info$condition)
levels(group) <- make.names(levels(group))
dge <- DGEList(counts = DGE_counts, group = group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# === 6. Calculate CPM ===
cpm_df <- as.data.frame(cpm(dge, log = FALSE)) %>%
    rownames_to_column(id_col_name)
colnames(cpm_df)[-1] <- paste0("CPM_", colnames(cpm_df)[-1])

# === 7. Estimate Dispersion (Paired / Unpaired 自動判別) ===
if ("subject" %in% colnames(sample_info)) {
    cat("INFO: --- 'subject' column detected: Running PAIRED analysis ---\n")
    subject <- factor(sample_info$subject)
    design <- model.matrix(~0 + group + subject)
    # Contrast用に前半の列名を揃える
    colnames(design)[1:length(levels(group))] <- levels(group)
} else {
    cat("INFO: --- Running UNPAIRED analysis ---\n")
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
}

dge <- estimateDisp(dge, design)

# === 8. ANOVA-like F-test (3+ groups) ===
if (length(unique_conditions) >= 3) {
    cat("INFO: --- Performing ANOVA-like F-test ---\n")
    fit <- glmQLFit(dge, design)
    
    # 基準となる1群に対する残りの全群の差を検定（全てに差がないという帰無仮説）
    base_group <- levels(group)[1]
    anova_contrasts <- list()
    for (g in levels(group)[-1]) {
        anova_contrasts[[paste0(g, "_vs_", base_group)]] <- paste(make.names(g), "-", make.names(base_group))
    }
    
    contrast_matrix_anova <- makeContrasts(contrasts = unlist(anova_contrasts), levels = design)
    qlf_anova <- glmQLFTest(fit, contrast = contrast_matrix_anova)
    
    res_anova <- topTags(qlf_anova, n=Inf)$table %>%
        rownames_to_column(id_col_name) %>%
        left_join(cpm_df, by = id_col_name)
    
    if (!is.null(symbol_map)) {
        res_anova <- res_anova %>% left_join(symbol_map, by = id_col_name)
    }
    
    base_cols <- c(id_col_name)
    if ("GeneSymbol" %in% colnames(res_anova)) {
        base_cols <- c(base_cols, "GeneSymbol")
    }
    stat_cols <- c("logCPM", "F", "PValue", "FDR")
    cpm_cols <- grep("^CPM_", colnames(res_anova), value = TRUE)
    
    res_anova <- res_anova %>% select(any_of(base_cols), any_of(stat_cols), all_of(cpm_cols))
    
    output_anova_filename <- paste0(prefix, "_ANOVA_like_edgeR_results.tsv")
    write_tsv(res_anova, output_anova_filename)
    cat(paste("INFO:     -> Saved ANOVA-like results to:", output_anova_filename, "\n"))
}

# === 9. Pairwise Tests ===
cat("INFO: --- Performing pairwise DEG tests ---\n")
fit <- glmQLFit(dge, design)

for (i in seq_along(contrasts_list)) {
    file_suffix <- gsub("\\s", "", file_suffixes[i])
    val_contrast <- valid_contrasts_list[i]

    cat(paste("INFO:   - Processing:", file_suffix, "\n"))
    contrast.matrix <- makeContrasts(contrasts=val_contrast, levels=design)
    qlf <- glmQLFTest(fit, contrast=contrast.matrix)
    
    res <- topTags(qlf, n=Inf)$table %>%
        rownames_to_column(id_col_name) %>%
        left_join(cpm_df, by = id_col_name)

    if (!is.null(symbol_map)) {
        res <- res %>% left_join(symbol_map, by = id_col_name)
    }

    base_cols <- c(id_col_name)
    if ("GeneSymbol" %in% colnames(res)) {
        base_cols <- c(base_cols, "GeneSymbol")
    }
    stat_cols <- c("logFC", "logCPM", "F", "PValue", "FDR")
    cpm_cols <- grep("^CPM_", colnames(res), value = TRUE)
    
    res <- res %>% select(any_of(base_cols), any_of(stat_cols), all_of(cpm_cols))
    
    output_filename <- paste0(prefix, "_", file_suffix, "_edgeR_results.tsv")
    write_tsv(res, output_filename)
}

cat("INFO: --- edgeR R script finished successfully! ---\n")
RSCRIPT_EOF

R_GENE_SYMBOL_MAP_FILE=${GENE_SYMBOL_MAP_FILE:-"NULL"}

if [[ "$EXEC_MODE" == "docker" ]]; then
  echo "INFO: Running R script inside Docker container..."
  docker run --rm -i -u "$(id -u):$(id -g)" -v "$(pwd)":/workdir -w /workdir "$EDGER_CONTAINER" \
    Rscript - "$SAMPLES_FILE" "$COUNTS_FILE" "$PREFIX" "$R_GENE_SYMBOL_MAP_FILE" <<< "$R_SCRIPT"
else
  echo "INFO: Running R script locally..."
  Rscript - "$SAMPLES_FILE" "$COUNTS_FILE" "$PREFIX" "$R_GENE_SYMBOL_MAP_FILE" <<< "$R_SCRIPT"
fi

echo "----------------------------------------------------"
echo "INFO: edgeR Pipeline Finished Successfully."
echo "========================================================================"
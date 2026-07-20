#!/bin/bash
# ==============================================================================
# DRIMSeq Alternative Splicing Analysis Pipeline (v1.0 - Integrated)
# ==============================================================================
# Description:
#   - Fixed path resolution issues for Docker/WSL.
#   - Includes full pairwise testing (Outputs both A vs B and B vs A).
#   - Supports 3+ groups (ANOVA + all pairwise combinations automatically).
# ==============================================================================

set -e
set -o pipefail

# --- Global Variables (Container Image) ---
DRIMSEQ_CONTAINER="ezojika7713/drimseq:1.22.0"

# --- Utility Functions ---
get_gtf_info() {
    if [[ "$SPECIES" == "human" ]]; then
        FTP_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
        if [ "$USE_BASIC" = true ]; then
            ANNOTATION_GTF_BASE="gencode.v${GENCODE_RELEASE}.basic.annotation.gtf.gz"
        else
            ANNOTATION_GTF_BASE="gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"
        fi
    else
        FTP_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${GENCODE_RELEASE}"
        if [ "$USE_BASIC" = true ]; then
            ANNOTATION_GTF_BASE="gencode.vM${GENCODE_RELEASE}.basic.annotation.gtf.gz"
        else
            ANNOTATION_GTF_BASE="gencode.vM${GENCODE_RELEASE}.annotation.gtf.gz"
        fi
    fi
    ANNOTATION_GTF="${DRIMSEQ_GENCODE_DIR}/${ANNOTATION_GTF_BASE}"
}

download_gtf_if_needed() {
    get_gtf_info
    if [ ! -f "$ANNOTATION_GTF" ]; then
        echo "INFO: Downloading GENCODE GTF file..."
        wget -c -P "$DRIMSEQ_GENCODE_DIR" "${FTP_BASE}/${ANNOTATION_GTF_BASE}"
    else
        echo "INFO: GTF file already exists: $ANNOTATION_GTF"
    fi
}

validate_file() {
    local file=$1; local description=$2; local help_msg=$3
    if [ ! -f "$file" ]; then
        echo "ERROR: $description file not found at: $file" >&2
        if [ -n "$help_msg" ]; then echo "$help_msg" >&2; fi
        return 1
    fi
    return 0
}

validate_positive_integer() {
    if ! [[ "$1" =~ ^[0-9]+$ ]] || [ "$1" -lt 1 ]; then
        echo "ERROR: $2 must be a positive integer (got: $1)" >&2
        return 1
    fi
    return 0
}

# --- Help Message ---
usage() {
  cat <<EOM
Usage: $(basename "$0") <samples.csv> <species> <gencode_release> [count_matrix.tsv] [OPTIONS]

Required Arguments:
  <samples.csv>         Sample information CSV file.
  <species>             "human" or "mouse".
  <gencode_release>     GENCODE release version.

Optional Positional Argument:
  [count_matrix.tsv]    Path to count matrix TSV file.

Options:
  --output-dir <path>         Output directory.
  --tx2gene <path>            Path to tx2gene mapping file.
  --add-gene-symbol [path]    Add gene symbol column.
  --threads <n>               Number of threads (default: 8).
  --pvalue <n>                P-value threshold (default: 0.05).
  --ofdr <n>                  Overall FDR (default: 0.05).
  --make-tx2gene              Create tx2gene mapping file from GTF.
  --make-geneid2symbol        Create geneid2symbol mapping file from GTF.
  --all-transcripts           Use comprehensive annotation.
  --basic-protein-coding      Use GENCODE basic protein-coding.
  --local                     Use locally installed R.
  -h, --help                  Display this help message.
EOM
  exit 0
}

# --- 1. Argument Parsing ---
if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then usage; fi
if [[ "$1" == "-v" || "$1" == "--version" ]]; then echo "Container: ${DRIMSEQ_CONTAINER}"; exit 0; fi

if [[ $# -lt 3 ]]; then echo "ERROR: Missing required arguments."; exit 1; fi

SAMPLES_FILE=$1
SPECIES=$2
GENCODE_RELEASE=$3
shift 3

if [ -f "$SAMPLES_FILE" ]; then SAMPLES_FILE=$(realpath "$SAMPLES_FILE"); fi
if [ ! -f "$SAMPLES_FILE" ]; then echo "ERROR: Sample file not found: $SAMPLES_FILE"; exit 1; fi
if ! [[ "$SPECIES" =~ ^(human|mouse)$ ]]; then echo "ERROR: Species must be 'human' or 'mouse'."; exit 1; fi

COUNT_MATRIX_PATH=""
if [[ $# -gt 0 ]] && [[ "$1" != -* ]]; then
    COUNT_MATRIX_PATH="$1"
    shift
fi
if [[ -n "$COUNT_MATRIX_PATH" && -f "$COUNT_MATRIX_PATH" ]]; then
    COUNT_MATRIX_PATH=$(realpath "$COUNT_MATRIX_PATH")
fi

EXEC_MODE="docker"
INDEX_SUFFIX="pc"
USE_BASIC=false
OUTPUT_DIR=""
TX2GENE_PATH=""
GENE_SYMBOL_PATH=""
ADD_GENE_SYMBOL=false
MAKE_TX2GENE=false
MAKE_GENEID2SYMBOL=false
THREADS=8
PVALUE_THRESHOLD=0.05
OFDR=0.05

MIN_SAMPS_GENE_EXPR="all"
MIN_GENE_EXPR=10
MIN_SAMPS_FEATURE_EXPR="auto"
MIN_FEATURE_EXPR=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local) EXEC_MODE="local"; shift ;;
    --all-transcripts) INDEX_SUFFIX="all"; shift ;;
    --basic-protein-coding) INDEX_SUFFIX="pc"; USE_BASIC=true; shift ;;
    --basic-all-transcripts) INDEX_SUFFIX="all"; USE_BASIC=true; shift ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --tx2gene) TX2GENE_PATH="$2"; shift 2 ;;
    --add-gene-symbol)
        ADD_GENE_SYMBOL=true
        if [[ -n "$2" ]] && [[ "$2" != -* ]]; then GENE_SYMBOL_PATH="$2"; shift 2; else shift; fi
        ;;
    --make-tx2gene) MAKE_TX2GENE=true; shift ;;
    --make-geneid2symbol) MAKE_GENEID2SYMBOL=true; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --pvalue) PVALUE_THRESHOLD="$2"; shift 2 ;;
    --ofdr) OFDR="$2"; shift 2 ;;
    --min-samps-gene-expr) MIN_SAMPS_GENE_EXPR="$2"; shift 2 ;;
    --min-gene-expr) MIN_GENE_EXPR="$2"; shift 2 ;;
    --min-samps-feature-expr) MIN_SAMPS_FEATURE_EXPR="$2"; shift 2 ;;
    --min-feature-expr) MIN_FEATURE_EXPR="$2"; shift 2 ;;
    *) echo "ERROR: Unknown option: $1"; exit 1 ;;
  esac
done

PREFIX=$(basename "$SAMPLES_FILE" .csv)

# --- 2. Construct Version Prefix & Paths ---
SPECIES_PREFIX=""; if [[ "$SPECIES" == "mouse" ]]; then SPECIES_PREFIX="M"; fi
if [ "$USE_BASIC" = true ]; then
    VERSION_PREFIX="v${SPECIES_PREFIX}${GENCODE_RELEASE}.${INDEX_SUFFIX}.basic"
else
    VERSION_PREFIX="v${SPECIES_PREFIX}${GENCODE_RELEASE}.${INDEX_SUFFIX}"
fi

BASE_OUT_DIR="."
BASE_OUT_DIR=$(cd "$BASE_OUT_DIR" && pwd)

DRIMSEQ_GENCODE_DIR="${BASE_OUT_DIR}/${VERSION_PREFIX}_gencode_files"

if [ -z "$OUTPUT_DIR" ]; then
    DRIMSEQ_OUT_DIR="${BASE_OUT_DIR}/${VERSION_PREFIX}_drimseq"
else
    mkdir -p "$OUTPUT_DIR"
    DRIMSEQ_OUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
fi
mkdir -p "$DRIMSEQ_GENCODE_DIR"

if [[ -n "$TX2GENE_PATH" && -f "$TX2GENE_PATH" ]]; then TX2GENE_PATH=$(realpath "$TX2GENE_PATH"); fi
if [[ -n "$GENE_SYMBOL_PATH" && -f "$GENE_SYMBOL_PATH" ]]; then GENE_SYMBOL_PATH=$(realpath "$GENE_SYMBOL_PATH"); fi


# --- 3. GTF Preparation (Helper Mode) ---
if [ "$MAKE_TX2GENE" = true ] || [ "$MAKE_GENEID2SYMBOL" = true ]; then
    echo "INFO: === GENCODE Annotation File Preparation ==="
    download_gtf_if_needed
    
    if [ "$MAKE_TX2GENE" = true ]; then
        TX2GENE_OUTPUT="${DRIMSEQ_GENCODE_DIR}/tx2gene_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"
        echo "INFO: Creating tx2gene: $TX2GENE_OUTPUT"
        zcat "$ANNOTATION_GTF" | awk -F'\t' '$3 == "transcript" {match($9, /transcript_id "([^"]+)"/); t=substr($9,RSTART+15,RLENGTH-16); match($9, /gene_id "([^"]+)"/); g=substr($9,RSTART+9,RLENGTH-10); print t "\t" g}' > "$TX2GENE_OUTPUT"
    fi
    if [ "$MAKE_GENEID2SYMBOL" = true ]; then
        GENEID2SYMBOL_OUTPUT="${DRIMSEQ_GENCODE_DIR}/geneid2symbol_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"
        echo "INFO: Creating geneid2symbol: $GENEID2SYMBOL_OUTPUT"
        zcat "$ANNOTATION_GTF" | awk -F'\t' '$3 == "gene" {match($9, /gene_id "([^"]+)"/); g=substr($9,RSTART+9,RLENGTH-10); match($9, /gene_name "([^"]+)"/); n=substr($9,RSTART+11,RLENGTH-12); if(n!="") print g "\t" n}' > "$GENEID2SYMBOL_OUTPUT"
    fi
    exit 0
fi

# --- 4. Analysis Setup ---
mkdir -p "$DRIMSEQ_OUT_DIR"
LOG_FILE="${DRIMSEQ_OUT_DIR}/run_drimseq_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================================================"
echo "INFO: DRIMSeq Analysis Pipeline Started (v1.0 - Integrated)"
echo "INFO: Date: $(date)"
echo "INFO: Target: $PREFIX"
echo "========================================================================"

if [ -z "$COUNT_MATRIX_PATH" ]; then
    POSSIBLE_PATH="${BASE_OUT_DIR}/${VERSION_PREFIX}_salmon_quant/${PREFIX}_transcript_counts_matrix.tsv"
    if [ -f "$POSSIBLE_PATH" ]; then
        COUNT_MATRIX_PATH="$POSSIBLE_PATH"
        echo "INFO: Auto-detected count matrix: $COUNT_MATRIX_PATH"
    fi
fi
validate_file "$COUNT_MATRIX_PATH" "Count matrix" || exit 1
COUNT_MATRIX_PATH=$(realpath "$COUNT_MATRIX_PATH")

if [ -z "$TX2GENE_PATH" ]; then
    TX2GENE_PATH="${DRIMSEQ_GENCODE_DIR}/tx2gene_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"
    if [ ! -f "$TX2GENE_PATH" ]; then
        echo "INFO: Creating auto-detected tx2gene map..."
        download_gtf_if_needed
        zcat "$ANNOTATION_GTF" | awk -F'\t' '$3 == "transcript" {match($9, /transcript_id "([^"]+)"/); t=substr($9,RSTART+15,RLENGTH-16); match($9, /gene_id "([^"]+)"/); g=substr($9,RSTART+9,RLENGTH-10); print t "\t" g}' > "$TX2GENE_PATH"
    fi
fi
validate_file "$TX2GENE_PATH" "tx2gene" || exit 1
TX2GENE_PATH=$(realpath "$TX2GENE_PATH")

GENE_SYMBOL_PATH_FOR_R="NA"
if [ "$ADD_GENE_SYMBOL" = true ]; then
    if [ -z "$GENE_SYMBOL_PATH" ]; then
        GENE_SYMBOL_PATH="${DRIMSEQ_GENCODE_DIR}/geneid2symbol_gencode_${SPECIES}_v${GENCODE_RELEASE}.tsv"
        if [ ! -f "$GENE_SYMBOL_PATH" ]; then
            echo "INFO: Creating auto-detected gene symbol map..."
            download_gtf_if_needed
            zcat "$ANNOTATION_GTF" | awk -F'\t' '$3 == "gene" {match($9, /gene_id "([^"]+)"/); g=substr($9,RSTART+9,RLENGTH-10); match($9, /gene_name "([^"]+)"/); n=substr($9,RSTART+11,RLENGTH-12); if(n!="") print g "\t" n}' > "$GENE_SYMBOL_PATH"
        fi
    fi
    validate_file "$GENE_SYMBOL_PATH" "Gene symbol" || exit 1
    GENE_SYMBOL_PATH=$(realpath "$GENE_SYMBOL_PATH")
    GENE_SYMBOL_PATH_FOR_R="$GENE_SYMBOL_PATH"
fi

# --- 5. R Script Content ---
FINAL_RDS_FILE="${DRIMSEQ_OUT_DIR}/drimseq_object_fitted.Rds"

R_SCRIPT_CONTENT=$(cat <<'EOF'
# --- DRIMSeq Analysis R Script ---
options(error = function() {
    cat("\nERROR: === R ERROR TRACEBACK ===\n")
    traceback(2)
    cat("================================\n")
    quit(status = 1)
})

suppressPackageStartupMessages({
    library(DRIMSeq)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(BiocParallel)
    library(stageR)
    library(rlang)
})

args <- commandArgs(trailingOnly = TRUE)
samples_file_path <- args[1]; count_matrix_path <- args[2]; tx2gene_map_path <- args[3]
drimseq_out_dir <- args[4]; gene_symbol_path <- args[5]; threads <- as.integer(args[6])
min_samps_gene_expr <- args[7]; min_gene_expr <- as.integer(args[8])
min_samps_feature_expr <- args[9]; min_feature_expr <- as.integer(args[10])
pvalue_threshold <- as.numeric(args[11]); ofdr <- as.numeric(args[12])
prefix <- args[13]

available_cores <- parallel::detectCores()
optimal_threads <- min(threads, max(1, available_cores - 1))
if (optimal_threads < threads) threads <- optimal_threads

cat("INFO: --- Loading inputs ---\n")
samples_full <- read.csv(samples_file_path, header = TRUE, strip.white = TRUE)
count_matrix <- read_tsv(count_matrix_path, show_col_types = FALSE)
tx2gene <- read_tsv(tx2gene_map_path, col_names = c("TXNAME", "GENEID"), show_col_types = FALSE)

add_gene_symbol <- gene_symbol_path != "NA"
gene2symbol <- NULL
if (add_gene_symbol) {
    gene2symbol <- read_tsv(gene_symbol_path, col_names = c("GENEID", "GENESYMBOL"), 
                            show_col_types = FALSE, col_types = "cc")
    cat(paste0("INFO:   Loaded ", nrow(gene2symbol), " gene symbols\n"))
}

transcript_ids <- count_matrix[[1]]
cts <- as.matrix(count_matrix[, -1])
rownames(cts) <- transcript_ids
cts <- cts[, samples_full$file_name, drop = FALSE]
colnames(cts) <- samples_full$sample_name
cts <- cts[rowSums(cts) > 0, ]

counts <- data.frame(
    gene_id = tx2gene$GENEID[match(rownames(cts), tx2gene$TXNAME)],
    feature_id = rownames(cts),
    cts, check.names = FALSE
)
counts <- counts[!is.na(counts$gene_id), ]

cat("INFO: --- Creating DRIMSeq object ---\n")
# Always use Unpaired design (sample_name and condition) for mathematical stability
samples_drim <- samples_full[, c("sample_name", "condition")]
colnames(samples_drim) <- c("sample_id", "group")

original_group_names <- levels(as.factor(samples_drim$group))
sanitized_group_names <- make.names(original_group_names)
names(original_group_names) <- sanitized_group_names
samples_drim$group <- factor(make.names(samples_drim$group), levels = sanitized_group_names)

d <- dmDSdata(counts = counts, samples = samples_drim)

cat("INFO: --- Filtering & Modeling ---\n")
n_min <- min(table(samples_drim$group))
msge <- if (min_samps_gene_expr == "all") nrow(samples_drim) else if (min_samps_gene_expr == "auto") n_min else as.integer(min_samps_gene_expr)
msfe <- if (min_samps_feature_expr == "auto") n_min else as.integer(min_samps_feature_expr)

d_filtered <- dmFilter(d, min_samps_gene_expr = msge, min_gene_expr = min_gene_expr, 
                       min_samps_feature_expr = msfe, min_feature_expr = min_feature_expr, 
                       run_gene_twice = TRUE)

set.seed(123)

cat("INFO: --- Fitting ---\n")
design_full <- model.matrix(~ group, data = DRIMSeq::samples(d_filtered))

d_precision <- dmPrecision(d_filtered, design = design_full, BPPARAM = MulticoreParam(workers = threads))
pdf(file.path(drimseq_out_dir, "precision_plot.pdf")); plotPrecision(d_precision); dev.off()
d_fitted <- dmFit(d_precision, design = design_full, BPPARAM = MulticoreParam(workers = threads))

# --- Proportions ---
props <- as.data.frame(DRIMSeq::proportions(d_fitted))
colnames(props) <- c("gene_id", "feature_id", make.names(colnames(props)[-c(1,2)]))
props_long <- props %>% pivot_longer(-c(gene_id, feature_id), names_to = "sample_id", values_to = "proportion")
samps_info <- DRIMSeq::samples(d_fitted)
samps_info$sample_id <- make.names(samps_info$sample_id)
props_long <- props_long %>% left_join(samps_info, by = "sample_id")
mean_props <- props_long %>% group_by(gene_id, feature_id, group) %>% summarize(mean_proportion = mean(proportion), .groups = "drop")
mean_props_wide <- mean_props %>% pivot_wider(names_from = group, values_from = mean_proportion, names_prefix = "mean_prop_")

for (i in seq_along(sanitized_group_names)) {
    old_name <- paste0("mean_prop_", sanitized_group_names[i])
    new_name <- paste0("mean_prop_", original_group_names[sanitized_group_names[i]])
    if (old_name %in% colnames(mean_props_wide)) colnames(mean_props_wide)[colnames(mean_props_wide) == old_name] <- new_name
}

# === Helper: Robust Symbol & Reorder ===
add_symbol_if_requested <- function(df) {
    if (add_gene_symbol && !is.null(gene2symbol)) {
        if (!"gene_symbol" %in% colnames(df)) {
            df$gene_id <- as.character(df$gene_id)
            gene2symbol$GENEID <- as.character(gene2symbol$GENEID)
            df_joined <- df %>% left_join(dplyr::distinct(gene2symbol), by = c("gene_id" = "GENEID"))
            if ("GENESYMBOL" %in% colnames(df_joined)) {
                df <- df_joined %>% dplyr::rename(gene_symbol = GENESYMBOL)
            } else {
                df <- df_joined
            }
        }
    }
    
    base_cols <- c("gene_id")
    if (add_gene_symbol && "gene_symbol" %in% colnames(df)) base_cols <- c("gene_id", "gene_symbol")

    if ("feature_id" %in% colnames(df)) {
        stats_cols <- c("feature_id", "lr", "df", "pvalue", "adj_pvalue")
        if ("stageR_adjP" %in% colnames(df)) stats_cols <- c(stats_cols, "stageR_adjP")
        delta_prop_cols <- grep("^delta_prop_", colnames(df), value = TRUE)
        if (length(delta_prop_cols) > 0) stats_cols <- c(stats_cols, delta_prop_cols)
        present_cols <- stats_cols[stats_cols %in% colnames(df)]
        df <- df %>% relocate(all_of(base_cols)) %>% relocate(all_of(present_cols), .after = last(base_cols))
    } else { 
        stats_cols <- c("lr", "df", "pvalue", "adj_pvalue")
        if ("stageR_gene_adjP" %in% colnames(df)) stats_cols <- c(stats_cols, "stageR_gene_adjP")
        present_cols <- stats_cols[stats_cols %in% colnames(df)]
        df <- df %>% relocate(all_of(base_cols)) %>% relocate(all_of(present_cols), .after = last(base_cols))
    }
    return(df)
}

cat("INFO: --- Testing ---\n")
n_groups <- nlevels(DRIMSeq::samples(d_fitted)$group)

if (n_groups > 2) {
    cat("INFO: Detected >2 groups. Running ANOVA logic...\n")
    # --- 1. Global Test (ANOVA) ---
    d_tested_main <- dmTest(d_fitted, coef = 2:n_groups, BPPARAM = MulticoreParam(workers = threads))
    res_main_gene <- results(d_tested_main, level = "gene")
    res_main_txp <- results(d_tested_main, level = "feature") %>% 
        left_join(mean_props_wide, by = c("gene_id", "feature_id")) %>% 
        left_join(props, by = c("gene_id", "feature_id"))
    
    cat(paste0("INFO: Running stageR for ANOVA (OFDR = ", ofdr, ")\n"))
    pScreen <- res_main_gene$pvalue; names(pScreen) <- res_main_gene$gene_id
    pConfirmation <- matrix(res_main_txp$pvalue, ncol = 1)
    rownames(pConfirmation) <- res_main_txp$feature_id; colnames(pConfirmation) <- "transcript"
    tx2gene_stager <- res_main_txp[, c("feature_id", "gene_id")]; colnames(tx2gene_stager) <- c("transcript", "gene")
    
    stageRObj <- tryCatch({ stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=FALSE, tx2gene=tx2gene_stager) }, error=function(e) NULL)
    
    if(!is.null(stageRObj)) {
        stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=ofdr, allowNA=TRUE)
        padj_stager <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
        gene_padj <- padj_stager %>% group_by(geneID) %>% summarize(stageR_gene_adjP = min(gene, na.rm=TRUE), .groups="drop")
        res_main_gene <- res_main_gene %>% left_join(gene_padj, by=c("gene_id"="geneID"))
        tx_cols <- setdiff(colnames(padj_stager), c("txID", "geneID", "gene"))
        if(length(tx_cols) > 0) {
             res_main_txp <- res_main_txp %>% left_join(padj_stager %>% select(txID, stageR_adjP = all_of(tx_cols[1])), by=c("feature_id"="txID"))
        }
    } 

    write.table(res_main_gene %>% add_symbol_if_requested(), file.path(drimseq_out_dir, paste0(prefix, "_drimseq_results_gene_level_ANOVA.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(res_main_txp %>% add_symbol_if_requested(), file.path(drimseq_out_dir, paste0(prefix, "_drimseq_results_transcript_level_ANOVA.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
    
    # --- 2. Pairwise Combinations Loop ---
    cat("INFO: Running all pairwise comparisons using global fit...\n")
    pair_combs <- combn(sanitized_group_names, 2, simplify = FALSE)
    
    for (pair in pair_combs) {
        g1_name <- pair[1]; g2_name <- pair[2]
        orig_g1_name <- original_group_names[g1_name]; orig_g2_name <- original_group_names[g2_name]
        cat(sprintf("INFO:   -> Testing %s vs %s...\n", orig_g1_name, orig_g2_name))
        
        idx1 <- which(sanitized_group_names == g1_name)
        idx2 <- which(sanitized_group_names == g2_name)
        contrast_vec <- rep(0, n_groups)
        if (idx1 > 1) contrast_vec[idx1] <- -1
        if (idx2 > 1) contrast_vec[idx2] <- 1
        
        d_tested_pair <- dmTest(d_fitted, contrast = contrast_vec, BPPARAM = MulticoreParam(workers = threads))
        res_pair_gene <- results(d_tested_pair, level = "gene")
        res_pair_txp <- results(d_tested_pair, level = "feature")
        
        mean_props_pair <- mean_props_wide %>% select(gene_id, feature_id, !!paste0("mean_prop_", orig_g1_name), !!paste0("mean_prop_", orig_g2_name))
        res_pair_txp <- res_pair_txp %>% left_join(mean_props_pair, by = c("gene_id", "feature_id"))
        res_pair_txp[[paste0("delta_prop_", orig_g1_name, "-", orig_g2_name)]] <- res_pair_txp[[paste0("mean_prop_", orig_g1_name)]] - res_pair_txp[[paste0("mean_prop_", orig_g2_name)]]
        
        pScreen <- res_pair_gene$pvalue; names(pScreen) <- res_pair_gene$gene_id
        pConfirmation <- matrix(res_pair_txp$pvalue, ncol = 1)
        rownames(pConfirmation) <- res_pair_txp$feature_id; colnames(pConfirmation) <- "transcript"
        
        stageRObj <- tryCatch({ stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=FALSE, tx2gene=tx2gene_stager) }, error=function(e) NULL)
        
        if(!is.null(stageRObj)) {
            stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=ofdr, allowNA=TRUE)
            padj_stager <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
            gene_padj <- padj_stager %>% group_by(geneID) %>% summarize(stageR_gene_adjP = min(gene, na.rm=TRUE), .groups="drop")
            res_pair_gene <- res_pair_gene %>% left_join(gene_padj, by=c("gene_id"="geneID"))
            tx_cols <- setdiff(colnames(padj_stager), c("txID", "geneID", "gene"))
            if(length(tx_cols) > 0) {
                 res_pair_txp <- res_pair_txp %>% left_join(padj_stager %>% select(txID, stageR_adjP = all_of(tx_cols[1])), by=c("feature_id"="txID"))
            }
        }
        
        target_samples <- samps_info %>% filter(group %in% c(g1_name, g2_name)) %>% pull(sample_id)
        res_pair_txp <- res_pair_txp %>% left_join(props %>% select(gene_id, feature_id, any_of(target_samples)), by = c("gene_id", "feature_id"))
        
        # Output A vs B
        fname_gene <- paste0(prefix, "_", orig_g1_name, "_vs_", orig_g2_name, "_drimseq_results_gene_level_pairwise.tsv")
        fname_txp  <- paste0(prefix, "_", orig_g1_name, "_vs_", orig_g2_name, "_drimseq_results_transcript_level_with_proportions.tsv")
        write.table(res_pair_gene %>% add_symbol_if_requested(), file.path(drimseq_out_dir, fname_gene), sep="\t", quote=FALSE, row.names=FALSE)
        write.table(res_pair_txp %>% add_symbol_if_requested(), file.path(drimseq_out_dir, fname_txp), sep="\t", quote=FALSE, row.names=FALSE)
        
        # Output B vs A
        res_txp_rev <- res_pair_txp
        delta_col_orig <- paste0("delta_prop_", orig_g1_name, "-", orig_g2_name)
        delta_col_rev <- paste0("delta_prop_", orig_g2_name, "-", orig_g1_name)
        if (delta_col_orig %in% colnames(res_txp_rev)) {
            res_txp_rev[[delta_col_rev]] <- -res_txp_rev[[delta_col_orig]]
            res_txp_rev[[delta_col_orig]] <- NULL
        }
        fname_txp_rev <- paste0(prefix, "_", orig_g2_name, "_vs_", orig_g1_name, "_drimseq_results_transcript_level_with_proportions.tsv")
        write.table(res_txp_rev %>% add_symbol_if_requested(), file.path(drimseq_out_dir, fname_txp_rev), sep="\t", quote=FALSE, row.names=FALSE)
    }
    
} else {
    cat("INFO: Detected 2 groups. Running pairwise logic.\n")
    g1_name <- sanitized_group_names[1]; g2_name <- sanitized_group_names[2]
    orig_g1_name <- original_group_names[g1_name]; orig_g2_name <- original_group_names[g2_name]
    
    d_tested_main <- dmTest(d_fitted, coef = 2, BPPARAM = MulticoreParam(workers = threads))
    res_main_gene <- results(d_tested_main, level = "gene")
    res_main_txp <- results(d_tested_main, level = "feature")
    
    mean_props_pair <- mean_props_wide %>% select(gene_id, feature_id, !!paste0("mean_prop_", orig_g1_name), !!paste0("mean_prop_", orig_g2_name))
    res_main_txp <- res_main_txp %>% left_join(mean_props_pair, by = c("gene_id", "feature_id"))
    res_main_txp[[paste0("delta_prop_", orig_g1_name, "-", orig_g2_name)]] <- res_main_txp[[paste0("mean_prop_", orig_g1_name)]] - res_main_txp[[paste0("mean_prop_", orig_g2_name)]]
    
    cat(paste0("INFO: Running stageR (OFDR = ", ofdr, ")\n"))
    pScreen <- res_main_gene$pvalue; names(pScreen) <- res_main_gene$gene_id
    pConfirmation <- matrix(res_main_txp$pvalue, ncol = 1)
    rownames(pConfirmation) <- res_main_txp$feature_id; colnames(pConfirmation) <- "transcript"
    tx2gene_stager <- res_main_txp[, c("feature_id", "gene_id")]; colnames(tx2gene_stager) <- c("transcript", "gene")
    
    stageRObj <- tryCatch({ stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=FALSE, tx2gene=tx2gene_stager) }, error=function(e) NULL)
    
    if(!is.null(stageRObj)) {
        stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=ofdr, allowNA=TRUE)
        padj_stager <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE)
        gene_padj <- padj_stager %>% group_by(geneID) %>% summarize(stageR_gene_adjP = min(gene, na.rm=TRUE), .groups="drop")
        res_main_gene <- res_main_gene %>% left_join(gene_padj, by=c("gene_id"="geneID"))
        tx_cols <- setdiff(colnames(padj_stager), c("txID", "geneID", "gene"))
        if(length(tx_cols) > 0) {
             res_main_txp <- res_main_txp %>% left_join(padj_stager %>% select(txID, stageR_adjP = all_of(tx_cols[1])), by=c("feature_id"="txID"))
        }
    } 
    
    cols_to_select <- c("gene_id", "feature_id", samples_full$sample_name)
    props_pair <- props %>% select(any_of(cols_to_select))
    res_main_txp <- res_main_txp %>% left_join(props_pair, by = c("gene_id", "feature_id"))

    # Output A vs B
    fname_gene <- paste0(prefix, "_", orig_g1_name, "_vs_", orig_g2_name, "_drimseq_results_gene_level_pairwise.tsv")
    fname_txp  <- paste0(prefix, "_", orig_g1_name, "_vs_", orig_g2_name, "_drimseq_results_transcript_level_with_proportions.tsv")
    write.table(res_main_gene %>% add_symbol_if_requested(), file.path(drimseq_out_dir, fname_gene), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(res_main_txp %>% add_symbol_if_requested(), file.path(drimseq_out_dir, fname_txp), sep="\t", quote=FALSE, row.names=FALSE)
    
    # Output B vs A
    cat(paste0("INFO: Creating reversed comparison: ", orig_g2_name, " vs ", orig_g1_name, "\n"))
    res_txp_rev <- res_main_txp
    delta_col_orig <- paste0("delta_prop_", orig_g1_name, "-", orig_g2_name)
    delta_col_rev <- paste0("delta_prop_", orig_g2_name, "-", orig_g1_name)
    
    if (delta_col_orig %in% colnames(res_txp_rev)) {
        res_txp_rev[[delta_col_rev]] <- -res_txp_rev[[delta_col_orig]]
        res_txp_rev[[delta_col_orig]] <- NULL
    }
    fname_txp_rev <- paste0(prefix, "_", orig_g2_name, "_vs_", orig_g1_name, "_drimseq_results_transcript_level_with_proportions.tsv")
    write.table(res_txp_rev %>% add_symbol_if_requested(), file.path(drimseq_out_dir, fname_txp_rev), sep="\t", quote=FALSE, row.names=FALSE)
}

saveRDS(d_fitted, file.path(drimseq_out_dir, "drimseq_object_fitted.Rds"))
cat("INFO: --- DRIMSeq analysis complete! ---\n")
EOF
)

# --- 6. Execution Block ---
if [ "$EXEC_MODE" == "docker" ]; then
    echo "INFO: Running DRIMSeq in Docker..."
    
    H_SAMPLES_DIR=$(dirname "$SAMPLES_FILE"); H_SAMPLES_NAME=$(basename "$SAMPLES_FILE")
    H_COUNTS_DIR=$(dirname "$COUNT_MATRIX_PATH"); H_COUNTS_NAME=$(basename "$COUNT_MATRIX_PATH")
    H_TX2GENE_DIR=$(dirname "$TX2GENE_PATH"); H_TX2GENE_NAME=$(basename "$TX2GENE_PATH")
    H_OUT_DIR="$DRIMSEQ_OUT_DIR"
    
    DOCKER_SYMBOL_MOUNT=""
    R_SYMBOL_PATH="NA"
    if [ "$GENE_SYMBOL_PATH_FOR_R" != "NA" ]; then
        H_SYMBOL_DIR=$(dirname "$GENE_SYMBOL_PATH_FOR_R"); H_SYMBOL_NAME=$(basename "$GENE_SYMBOL_PATH_FOR_R")
        DOCKER_SYMBOL_MOUNT="-v ${H_SYMBOL_DIR}:/data/symbols"
        R_SYMBOL_PATH="/data/symbols/${H_SYMBOL_NAME}"
    fi

    docker run --rm -i -u $(id -u):$(id -g) \
        -v "${H_SAMPLES_DIR}:/data/samples" \
        -v "${H_COUNTS_DIR}:/data/counts" \
        -v "${H_TX2GENE_DIR}:/data/tx2gene" \
        -v "${H_OUT_DIR}:/data/output" \
        $DOCKER_SYMBOL_MOUNT \
        -w /data/output \
        "$DRIMSEQ_CONTAINER" Rscript - \
        "/data/samples/${H_SAMPLES_NAME}" \
        "/data/counts/${H_COUNTS_NAME}" \
        "/data/tx2gene/${H_TX2GENE_NAME}" \
        "/data/output" \
        "$R_SYMBOL_PATH" \
        "$THREADS" "$MIN_SAMPS_GENE_EXPR" "$MIN_GENE_EXPR" "$MIN_SAMPS_FEATURE_EXPR" "$MIN_FEATURE_EXPR" \
        "$PVALUE_THRESHOLD" "$OFDR" "$PREFIX" <<< "$R_SCRIPT_CONTENT"
else
    echo "INFO: Running DRIMSeq Locally..."
    Rscript - \
        "$SAMPLES_FILE" "$COUNT_MATRIX_PATH" "$TX2GENE_PATH" "$DRIMSEQ_OUT_DIR" \
        "$GENE_SYMBOL_PATH_FOR_R" "$THREADS" \
        "$MIN_SAMPS_GENE_EXPR" "$MIN_GENE_EXPR" "$MIN_SAMPS_FEATURE_EXPR" "$MIN_FEATURE_EXPR" \
        "$PVALUE_THRESHOLD" "$OFDR" "$PREFIX" <<< "$R_SCRIPT_CONTENT"
fi

echo "========================================================================"
echo "INFO: Pipeline Finished Successfully."
echo "INFO: Results in: $DRIMSEQ_OUT_DIR"
echo "========================================================================"
# RNsequant
Automated RNA-seq analysis pipeline for human and mouse data (v0.9)

## Pipeline Architecture and File Placement
The main executable `rnsequant.sh` serves as a wrapper script. Rather than performing the data processing directly, it controls the execution flow and automatically calls the specialized submodule scripts (such as `run_fastq_download.sh`, `run_salmon.sh`, `run_edgeR_core.sh`, etc.) located inside the `rnsequant` directory.

**File Placement:**
To run this pipeline, you must place the `rnsequant.sh` file and the `rnsequant` directory (containing the sub-scripts) in the **same hierarchical level** within your local working directory.

Please set up your directory structure as follows before executing `rnsequant.sh`:

```text
Your_Working_Directory/
├── rnsequant.sh         # The wrapper script to execute
├── rnsequant/           # Directory containing sub-scripts (Required)
│   ├── run_fastq_download.sh
│   ├── run_salmon.sh
│   └── ...
└── sample_meta.csv      # The sample metadata you provide (CSV)
```

You only need to interact with and run the outermost `rnsequant.sh` wrapper. The pipeline will automatically pass the correct arguments and environment variables to each script inside the `rnsequant` directory, so there is no need to execute the internal scripts individually.

## Environment (Docker Supported)
This pipeline operates using **Docker containers by default**. You do not need to install complex bioinformatics tools (e.g., Salmon, edgeR, FastQC) on your local machine. As long as you have a **Docker-enabled environment** and a **metadata CSV file**, you can start your analysis immediately.

*Note: You can also use locally installed tools by appending the --local flag.*

## Required Files (Sample CSV)
To run the pipeline, you must provide a sample metadata CSV file (or a directory containing CSV files) as an argument. Depending on your experimental design, **you must use one of the following exact headers (first row)**:

**1. Unpaired Design**
```csv
file_name,sample_name,condition
SRR1234567,Sample_A1,Control
SRR1234568,Sample_A2,Control
SRR1234569,Sample_B1,Treatment
```

**2. Paired Design**
By adding a subject column, the pipeline will automatically process the data as a paired design.
```csv
file_name,sample_name,condition,subject
SRR1234567,Patient1_Pre,Control,Patient1
SRR1234568,Patient1_Post,Treatment,Patient1
SRR1234569,Patient2_Pre,Control,Patient2
```

## Default Workflow
If executed without special options, the wrapper script runs the following **default flow**:

1. **FastQ Preparation:** Downloads FastQ data from ENA/NCBI (local files can also be used).
2. **FastQC & Trimming:** Performs read quality control and trimming.
3. **Salmon:** Quantifies expression at both the gene and transcript levels.
4. **edgeR:** Performs Differential Expression Gene (DEG) analysis, after which the pipeline stops.

*Note: If you want to run DRIMSeq (DTU analysis) or QAPA (APA analysis) downstream, you need to add the respective option flags.*

## Usage and Arguments

### Basic Syntax
```bash
bash rnsequant.sh <csv_file|directory> <species> <gencode_release> [OPTIONS]
```

### Required Arguments
1. `<csv_file|directory>`: Path to a single .csv sample metadata file, or a directory containing .csv files.
2. `<species>`: Target species (human or mouse).
3. `<gencode_release>`: GENCODE release version to use (e.g., 46).

### Key Options
Flags to control the analysis flow and inputs:
* `-o, --output-dir <dir>`: Specifies the base directory for output results.
* `--start-at <step>`: Starts or resumes the pipeline from a specific step. *(Valid options: download, fastqc, salmon, edger, drimseq, qapa)*
* `--local-fastq <dir>`: Skips the download step and uses local FastQ files located in `<dir>`.
* `--run-drimseq`: Executes Differential Transcript Usage (DTU) analysis using DRIMSeq after the default edgeR analysis.
* `--run-qapa`: Executes Alternative Polyadenylation (APA) analysis using QAPA after the edgeR analysis (includes a QAPA-specific DRIMSeq step).
  **Note: QAPA analysis has strict GENCODE release version requirements. It is only supported for human release 31 and mouse release 22.**
* `--retry [num]`: Sets the number of automatic retries on failure (defaults to 2 if omitted).
* `--local`: Executes using locally installed tools instead of Docker containers.

### Transcript Type Options
Specifies the type of reference to use for quantification (Defaults to GENCODE basic protein-coding transcripts).
* `--all-transcripts`: Use all transcripts.
* `--basic-protein-coding`: Use GENCODE basic protein-coding transcripts.
* `--basic-all-transcripts`: Use GENCODE basic all transcripts.

## Tool Versions
The following tool versions are utilized when running the pipeline in Docker mode:
* **SRA-Toolkit:** 3.2.1
* **FastQC:** 0.12.1
* **Trim Galore!:** 0.6.10
* **Salmon:** 1.10.3
* **edgeR:** 3.36.0
* **DRIMSeq:** 1.22.0
* **QAPA:** 1.3.3
* **tximport:** 1.22
* **MultiQC:** 1.10.1

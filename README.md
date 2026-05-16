# merip
# RNA-seq & MeRIP-seq Analysis Pipeline

This repository contains R scripts for identifying differential m6A methylation peaks, performing differential expression (DE) analysis, and visualizing GO enrichment results.

## 1. Prerequisites
Ensure the following R libraries are installed:
* **Analysis**: `exomePeak`, `DESeq2`, `data.table`
* **Visualization**: `ggplot2`, `gridExtra`

## 2. Workflow

### Step 1: Differential m6A Peak Calling
Using the `exomePeak` package to identify differential methylation peaks between treated and control conditions using BAM files and a GTF annotation.
* **Input**: IP and Input BAM files (2 replicates each per condition).
* **Output**: `diff_peak.txt` (contains peak coordinates and fold changes).

### Step 2: Differential Expression (DE) Analysis
Using `DESeq2` to identify differentially expressed genes (DEGs) from raw count matrices.
* **Filter**: Removes genes with 0 reads across all samples.
* **Normalization**: Median-of-ratios method via `DESeq()`.
* **Output**: `C4_diff.txt` (standardized results including $log_2$ FoldChange and $p$-adjusted values).

### Step 3: Data Extraction
Extracted gene names associated with significantly "Up" regulated m6A peaks for downstream functional analysis.
* **Output**: `name_up.txt`.

### Step 4: GO Enrichment & Visualization
Processes GO term enrichment data (Biological Process, Molecular Function, Cellular Component).

* **Data Processing**: 
    * Filters the top 10 terms per category based on FDR.
    * Reorders terms by gene count for better visualization.
* **Visualizations**:
    * **Bar Charts**: Faceted by category (`BP`, `MF`, `CC`) showing gene counts and significance (P-value/FDR).
    * **Bubble Plots**: Displaying Gene Ratio vs. Description, colored by $p$-adj.
    * **Combined Plots**: Uses `gridExtra` to merge multiple plots into a single output (e.g., `GO_fdr_combined.png`).

## 3. Key Files & Directories
* `/BAM/`: Input alignment files.
* `exomePeak_output/`: Differential peak results.
* `C4_diff.txt`: DESeq2 differential expression results.
* `GO_FDR10_up.txt`: Top 10 enriched GO terms for upregulated genes.
* `GO.png` / `GO_fdr_combined.png`: Final visualization plots.

## 4. Usage Note
Paths in the scripts are currently set to absolute local directories (e.g., `/data1/zhouyiming/...`). Ensure these are updated to your local environment before running.

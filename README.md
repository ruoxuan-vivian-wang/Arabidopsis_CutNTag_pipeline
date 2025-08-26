# Arabidopsis_CutNTag_pipeline

**Author:** Ruoxuan  
**Date:** 2025-07-03  
**Description:** A general pipeline for processing Cut&Tag data.  
It performs quality control, genome alignment, peak calling, bigwig generation, peak merging, and raw counts summarization.

---

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Pipeline Overview](#pipeline-overview)
- [Usage](#usage)
- [Outputs](#outputs)
- [Notes](#notes)

---

## Requirements

- **Operating System:** Linux / macOS  
- **Tools:**  
  - `micromamba` or `conda` for environment management  
  - `fastp` for QC and adapter trimming  
  - `bowtie2` for genome alignment  
  - `samtools` for BAM file processing  
  - `MACS2` for peak calling  
  - `deepTools` for bigwig generation and heatmaps  
  - `bedtools` for peak merging  

---

## Installation

```bash
# Using micromamba (recommended)
eval "$(~/bin/micromamba shell hook --shell bash)"
micromamba create -n cnt_env python=3.10 -y
micromamba activate cnt_env
micromamba install -y -c bioconda fastp bowtie2 samtools macs2 bedtools
micromamba install -y -c bioconda -c conda-forge deeptools=3.5.4
```

---

## Pipeline Overview

The pipeline performs the following steps:

1. **QC & Adapter Trimming** (`fastp`)

   * Removes adapters, trims low-quality bases, and generates HTML/JSON QC reports.

2. **Genome Alignment** (`bowtie2`)

   * Aligns reads to the reference genome, outputs sorted BAM files.

3. **Peak Calling** (`MACS2`)

   * Calls narrow peaks using treatment and control BAM files.

4. **BigWig Generation** (`bamCoverage` from `deepTools`)

   * Creates normalized signal tracks for visualization.

5. **Peak Merging** (`bedtools merge`)

   * Combines peaks from multiple samples into a single non-redundant peak set.

6. **Raw Counts Summarization** (`multiBamSummary`)

   * Generates a table of read counts per peak for downstream analysis.

---

## Usage

1. **Clone the repository:**

```bash
git clone https://github.com/ruoxuan-vivian-wang/Arabidopsis_CutNTag_pipeline.git
cd Arabidopsis_CutNTag_pipeline
```

2. **Edit the pipeline variables in `pipeline.sh`:**

```bash
# Example
R1_IN=/path/to/control_R1.fq.gz
R2_IN=/path/to/control_R2.fq.gz
R1_AC=/path/to/treatment_R1.fq.gz
R2_AC=/path/to/treatment_R2.fq.gz

SAMPLE_IN=control
SAMPLE_AC=treatment

GENOME_INDEX=/path/to/genome_index
GENOME_SIZE=1.35e8
THREADS=16
OUTDIR=/path/to/output
```

3. **Run the pipeline:**

```bash
bash cutntag_pipeline.sh
```

---

## Outputs

| File / Folder                     | Description                            |
| --------------------------------- | -------------------------------------- |
| `*_fastp_*.html / *_fastp_*.json` | QC reports after trimming              |
| `*.sorted.bam`                    | Sorted BAM files                       |
| `*.sorted.bam.bai`                | BAM index files                        |
| `macs2_output/*_peaks.narrowPeak` | Called narrow peaks                    |
| `*.bw`                            | Normalized bigwig tracks               |
| `matrix_*.gz`                     | ComputeMatrix output for heatmaps      |
| `heatmap_*.png`                   | Heatmaps of signal around peaks        |
| `all_peaks_sorted.bed`            | Sorted combined peaks from all samples |
| `combined_peaks_final.bed`        | Merged non-redundant peak set          |
| `counts_all.tab`                  | Raw read counts per peak               |

---

## Notes

* Adjust `GENOME_SIZE` according to the organism.
* NarrowPeak is used by default; remove `--broad` from MACS2 if previously included.
* For multiple samples, loop through each pair of control/treatment BAMs.
* Make sure sufficient CPU threads (`THREADS`) are allocated for faster processing.

---




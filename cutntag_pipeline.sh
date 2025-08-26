#!/bin/bash
# ==========================================
#   Cut&Tag Pipeline
#   Author: Ruoxuan
#   Date: 2025-07-03
#   Description: General pipeline for Cut&Tag data
# ==========================================

# -------- User-defined variables ----------
# Environment
ENV_NAME=cnt_env
PYTHON_VERSION=3.10

# Input files (paired-end)
# You can add multiple samples later in a loop
R1_IN=/path/to/sample_input_R1.fq.gz
R2_IN=/path/to/sample_input_R2.fq.gz
R1_AC=/path/to/sample_treatment_R1.fq.gz
R2_AC=/path/to/sample_treatment_R2.fq.gz

# Sample names (prefixes for output)
SAMPLE_IN=control_sample
SAMPLE_AC=treatment_sample

# Reference genome
GENOME_INDEX=/path/to/genome_index
GENOME_SIZE=1.35e8 # Example: Arabidopsis
THREADS=16

# Output directory
OUTDIR=/path/to/output

# Peaks BED file for multiBamSummary
COMBINED_PEAKS_BED=${OUTDIR}/combined_peaks.bed
# ------------------------------------------

# -------------------------
# Step 0: Setup environment
# -------------------------
eval "$(~/bin/micromamba shell hook --shell bash)"
micromamba create -n $ENV_NAME -c bioconda python=$PYTHON_VERSION -y
micromamba activate $ENV_NAME
micromamba install -y -c bioconda fastp bowtie2 samtools macs2
micromamba install -y -n $ENV_NAME -c bioconda -c conda-forge deeptools=3.5.4 bedtools

mkdir -p $OUTDIR

# -------------------------
# Step 1: QC & Adapter trimming
# -------------------------
fastp \
  -i $R1_IN -I $R2_IN \
  -o ${OUTDIR}/cleaned_${SAMPLE_IN}_R1.fq.gz \
  -O ${OUTDIR}/cleaned_${SAMPLE_IN}_R2.fq.gz \
  --detect_adapter_for_pe \
  -w $THREADS \
  -h ${OUTDIR}/fastp_${SAMPLE_IN}.html \
  -j ${OUTDIR}/fastp_${SAMPLE_IN}.json

fastp \
  -i $R1_AC -I $R2_AC \
  -o ${OUTDIR}/cleaned_${SAMPLE_AC}_R1.fq.gz \
  -O ${OUTDIR}/cleaned_${SAMPLE_AC}_R2.fq.gz \
  --detect_adapter_for_pe \
  -w $THREADS \
  -h ${OUTDIR}/fastp_${SAMPLE_AC}.html \
  -j ${OUTDIR}/fastp_${SAMPLE_AC}.json

# -------------------------
# Step 2: Genome Alignment
# -------------------------
for SAMPLE in $SAMPLE_IN $SAMPLE_AC
do
  bowtie2 -p $THREADS -x $GENOME_INDEX \
    -1 ${OUTDIR}/cleaned_${SAMPLE}_R1.fq.gz \
    -2 ${OUTDIR}/cleaned_${SAMPLE}_R2.fq.gz | \
  samtools view -bS - > ${OUTDIR}/${SAMPLE}.bam

  samtools sort -o ${OUTDIR}/${SAMPLE}.sorted.bam ${OUTDIR}/${SAMPLE}.bam -@ $THREADS
  samtools index ${OUTDIR}/${SAMPLE}.sorted.bam
done

# -------------------------
# Step 3: Peak Calling with MACS2
# -------------------------
macs2 callpeak \
  -t ${OUTDIR}/${SAMPLE_AC}.sorted.bam \
  -c ${OUTDIR}/${SAMPLE_IN}.sorted.bam \
  -f BAMPE -g $GENOME_SIZE -q 0.01 \
  --keep-dup all \
  -n ${SAMPLE_AC} \
  --outdir ${OUTDIR}/macs2_output/

# -------------------------
# Step 4: Generate BigWig
# -------------------------
bamCoverage -b ${OUTDIR}/${SAMPLE_AC}.sorted.bam \
  -o ${OUTDIR}/${SAMPLE_AC}.bw \
  --normalizeUsing CPM \
  --binSize 10

# -------------------------
# Step 5: Combine Peaks
# -------------------------
cat ${OUTDIR}/macs2_output/*_peaks.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n > ${OUTDIR}/all_peaks_raw.bed
bedtools merge -i ${OUTDIR}/all_peaks_raw.bed > ${OUTDIR}/combined_peaks.bed

# -------------------------
# Step 6: Generate raw counts for all samples
# -------------------------
# Index all BAMs before counting
for BAM in ${OUTDIR}/*.sorted.bam
do
  samtools index $BAM
done

multiBamSummary BED-file \
  --BED $COMBINED_PEAKS_BED \
  --bamfiles ${OUTDIR}/${SAMPLE_IN}.sorted.bam ${OUTDIR}/${SAMPLE_AC}.sorted.bam \
  -out ${OUTDIR}/counts_all.npz \
  --outRawCounts ${OUTDIR}/counts_all.tab \
  --smartLabels \
  --numberOfProcessors $THREADS


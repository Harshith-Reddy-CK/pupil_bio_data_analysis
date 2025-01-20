#!/bin/bash

# Custom Somatic Mutation Calling Pipeline
# Author: Harshith_Reddy
# Date: 17-01-2025

# Define paths and filenames
RAW_DATA_DIR="/path/to/raw_data"
OUTPUT_DIR="/path/to/output"
REFERENCE_GENOME="/path/to/reference_genome.fa"
FASTQC_OUTPUT="$OUTPUT_DIR/fastqc_results"
TRIMMING_OUTPUT="$OUTPUT_DIR/trimming_results"
ALIGNMENT_OUTPUT="$OUTPUT_DIR/alignment_results"
MUTATION_OUTPUT="$OUTPUT_DIR/mutation_results"
CANCER_FASTQ="$RAW_DATA_DIR/cancer_sample.fastq"
NORMAL_FASTQ="$RAW_DATA_DIR/normal_sample.fastq"
CANCER_BAM="$ALIGNMENT_OUTPUT/cancer_sample.bam"
NORMAL_BAM="$ALIGNMENT_OUTPUT/normal_sample.bam"

# Create necessary directories
mkdir -p $FASTQC_OUTPUT $TRIMMING_OUTPUT $ALIGNMENT_OUTPUT $MUTATION_OUTPUT

# 1. Quality Control
echo "Performing Quality Control using FastQC..."
fastqc -o $FASTQC_OUTPUT $CANCER_FASTQ $NORMAL_FASTQ
echo "Quality Control completed. Results saved in $FASTQC_OUTPUT."

# 2. Trimming Reads
echo "Trimming low-quality reads using fastp..."
fastp -i $CANCER_FASTQ -o $TRIMMING_OUTPUT/cancer_sample_trimmed.fastq -h $TRIMMING_OUTPUT/cancer_fastp_report.html
fastp -i $NORMAL_FASTQ -o $TRIMMING_OUTPUT/normal_sample_trimmed.fastq -h $TRIMMING_OUTPUT/normal_fastp_report.html
echo "Trimming completed. Trimmed reads saved in $TRIMMING_OUTPUT."

# 3. Alignment using BWA
echo "Aligning samples to the reference genome using BWA..."
bwa mem $REFERENCE_GENOME $TRIMMING_OUTPUT/cancer_sample_trimmed.fastq > $ALIGNMENT_OUTPUT/cancer_sample.sam
bwa mem $REFERENCE_GENOME $TRIMMING_OUTPUT/normal_sample_trimmed.fastq > $ALIGNMENT_OUTPUT/normal_sample.sam

# Convert SAM to BAM, sort, and index
samtools view -Sb $ALIGNMENT_OUTPUT/cancer_sample.sam | samtools sort -o $CANCER_BAM
samtools view -Sb $ALIGNMENT_OUTPUT/normal_sample.sam | samtools sort -o $NORMAL_BAM
samtools index $CANCER_BAM
samtools index $NORMAL_BAM
echo "Alignment completed. BAM files are ready and indexed."

# 4. Mutation Calling using Samtools and Bcftools
echo "Generating mpileup for mutation calling..."
samtools mpileup -f $REFERENCE_GENOME $CANCER_BAM $NORMAL_BAM -o $MUTATION_OUTPUT/samples.mpileup

# Call variants with bcftools
echo "Calling variants with bcftools..."
bcftools call -mv -Ov -o $MUTATION_OUTPUT/raw_variants.vcf $MUTATION_OUTPUT/samples.mpileup

# 5. Filtering variants
echo "Filtering variants to identify somatic mutations..."
bcftools filter -e 'QUAL<30 || DP<10' $MUTATION_OUTPUT/raw_variants.vcf -Ov -o $MUTATION_OUTPUT/filtered_variants.vcf

# Separate somatic mutations (present in cancer, absent in normal)
vcftools --vcf $MUTATION_OUTPUT/filtered_variants.vcf --keep-only-indels --out $MUTATION_OUTPUT/somatic_indels
vcftools --vcf $MUTATION_OUTPUT/filtered_variants.vcf --keep-only-snps --out $MUTATION_OUTPUT/somatic_snps

echo "Mutation calling and filtering completed. Results saved in $MUTATION_OUTPUT."

# 6. (Optional) Annotation or further analysis can be added here

echo "Pipeline execution completed. Check $OUTPUT_DIR for all results."


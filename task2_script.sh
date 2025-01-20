#!/bin/bash

# Script for Quality Control, Alignment, and Mutation Calling
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
CANCER_BAM="$ALIGNMENT_OUTPUT/cancer_sample.bam"
NORMAL_BAM="$ALIGNMENT_OUTPUT/normal_sample.bam"


# Flag to use trimming tools (set to 1 to enable)
USE_TRIMMING=1

# Choose the mutation calling tool: "Strelka2", "Mutect2", or "Custom"
TOOL_CHOICE="Mutect2"

# Create necessary directories
mkdir -p $FASTQC_OUTPUT $TRIMMING_OUTPUT $ALIGNMENT_OUTPUT $MUTATION_OUTPUT

# 1. Quality Control 
echo "Performing Quality Control using FastQC..."
fastqc -o $FASTQC_OUTPUT $RAW_DATA_DIR/*.fastq
echo "Quality Control completed. Results saved in $FASTQC_OUTPUT."

# Conditional Trimming Step
if [ $USE_TRIMMING -eq 1 ]; then
    echo "Trimming low-quality reads using fastp and Trimmomatic..."
    
    # Using fastp
    fastp -i $RAW_DATA_DIR/cancer_sample.fastq -o $TRIMMING_OUTPUT/cancer_sample_trimmed.fastq -h $TRIMMING_OUTPUT/cancer_fastp_report.html
    fastp -i $RAW_DATA_DIR/normal_sample.fastq -o $TRIMMING_OUTPUT/normal_sample_trimmed.fastq -h $TRIMMING_OUTPUT/normal_fastp_report.html

    # Using Trimmomatic
    java -jar trimmomatic.jar SE -phred33 $RAW_DATA_DIR/cancer_sample.fastq $TRIMMING_OUTPUT/cancer_sample_trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    java -jar trimmomatic.jar SE -phred33 $RAW_DATA_DIR/normal_sample.fastq $TRIMMING_OUTPUT/normal_sample_trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "Trimming completed. Trimmed reads saved in $TRIMMING_OUTPUT."
    RAW_DATA_DIR=$TRIMMING_OUTPUT
else
    echo "Skipping trimming step as USE_TRIMMING is set to 0."
fi

# 2. Alignment to the human genome (10 points)
echo "Aligning samples using BWA..."
bwa mem $REFERENCE_GENOME $RAW_DATA_DIR/cancer_sample_trimmed.fastq > $ALIGNMENT_OUTPUT/cancer_sample.sam
bwa mem $REFERENCE_GENOME $RAW_DATA_DIR/normal_sample_trimmed.fastq > $ALIGNMENT_OUTPUT/normal_sample.sam

# Convert SAM to BAM, sort, and index
samtools view -Sb $ALIGNMENT_OUTPUT/cancer_sample.sam | samtools sort -o $CANCER_BAM
samtools view -Sb $ALIGNMENT_OUTPUT/normal_sample.sam | samtools sort -o $NORMAL_BAM
samtools index $CANCER_BAM
samtools index $NORMAL_BAM
echo "Alignment completed. BAM files generated and indexed."

# 3. Mutation Calling (40 points)
echo "Starting mutation calling using $TOOL_CHOICE..."

case $TOOL_CHOICE in
    "Strelka2")
        echo "Calling somatic mutations using Strelka2..."
        strelka2 --normalBam $NORMAL_BAM --tumorBam $CANCER_BAM --referenceFasta $REFERENCE_GENOME --runDir $MUTATION_OUTPUT/strelka
        $MUTATION_OUTPUT/strelka/runWorkflow.py -m local -j 4
        ;;
    "Mutect2")
        echo "Calling somatic mutations using GATK/Mutect2..."
        gatk Mutect2 -R $REFERENCE_GENOME -I $CANCER_BAM -tumor CancerSample -I $NORMAL_BAM -normal NormalSample -O $MUTATION_OUTPUT/mutect2_somatic.vcf
        gatk FilterMutectCalls -V $MUTATION_OUTPUT/mutect2_somatic.vcf -O $MUTATION_OUTPUT/mutect2_filtered.vcf
        ;;
    "Custom")
        echo "Identifying somatic mutations with custom scripts..."
        samtools mpileup -f $REFERENCE_GENOME $CANCER_BAM $NORMAL_BAM | bcftools call -mv -Oz -o $MUTATION_OUTPUT/somatic_variants.vcf.gz
        bcftools index $MUTATION_OUTPUT/somatic_variants.vcf.gz
        ;;
    *)
        echo "Invalid TOOL_CHOICE selected. Please choose 'Strelka2', 'Mutect2', or 'Custom'."
        exit 1
        ;;
esac

echo "Mutation calling completed using $TOOL_CHOICE. Results saved in $MUTATION_OUTPUT."

# Additional analysis with Python/R can be added here

echo "All tasks completed successfully."
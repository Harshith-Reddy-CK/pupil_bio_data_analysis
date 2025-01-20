#🧬 Pupil Bio Data Analysis
Welcome to the Pupil Bio Data Analysis repository! This repository contains code and scripts used for the bioinformatics challenge related to identifying somatic mutations and performing quality control for cancer sample analysis.

##⚡ Overview of Tasks
<details> <summary>🧑‍🔬 Task 1: Coverage Analysis and Biomarker Identification</summary>
Coverage Analysis

📊 Calculate the median and coefficient of variation (CV) for single CpG coverage in each tissue.
📈 Generate plots summarizing the coverage statistics.
🔬 Biomarker Identification (20 points)

🧬 Identify PMPs (Potential Mutant Pairs) with high specificity for tissue differentiation. The goal is to minimize false positives for Tissue #1 while allowing some false negatives. Statistical and machine learning approaches are used to assign confidence (e.g., p-values) to each PMP.
🧮 Calculate the mean variant read fraction (VRF) for each PMP in both tissues.

</details> <details> <summary>🧪 Task 2: Quality Control, Alignment, and Mutation Calling</summary>
⚙️ Quality Control 

Perform quality checks using tools like FastQC, and summarize quality metrics (e.g., sequence counts, per-base quality, and read duplication levels).
💥 Alignment and Mutation Calling 

🧬 Align the samples to the human genome using tools like Bowtie2 or BWA.
🔍 Identify somatic mutations present in the cancer sample but absent in the normal tissue using custom scripts and open-source tools like Samtools, bcftools, and vcftools.
</details>
🖥️ Scripts Folder
The scripts folder contains Python and R scripts that automate the following tasks:

🧪 Quality control analysis using FastQC and other tools.
🧬 Alignment of sequencing data using tools such as Bowtie2 and BWA.
🔬 Mutation calling and somatic mutation identification using custom Python scripts leveraging Samtools and bcftools.
📊 Plots and visualizations for quality control and biomarker identification.
The scripts are organized into separate folders for each task, allowing users to easily navigate through the pipeline.

🚀 Usage
To use the scripts, follow these steps:

Clone the repository:

'''bash
git clone https://github.com/Harshith-Reddy-CK/pupil_bio_data_analysis.git
cd pupil_bio_data_analysis


Install dependencies listed in the requirements (e.g., Python, R packages, etc.).

Navigate to the appropriate folder (e.g., scripts/task1 or scripts/task2) and run the scripts as needed for each task.

💾 Supplementary and Output Files
The repository currently contains scripts for performing the tasks outlined above. Due to the size limitations of GitHub, large files, including supplementary data and output files, are not hosted directly on this repository.

To request the large files (e.g., BAM files, VCF files, output data), you can contact the repository owner or download them from the provided Dropbox or Google Drive links below:

<details> <summary>📥 Download Links</summary>
Dropbox Link to Output Files
Google Drive Link to Supplementary Files
</details>
Please reach out to me if you would like access to these files.

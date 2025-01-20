
#My Analysis codes Date: 15-01-2025

#QC
 fastqc *

 #trimming 
 fastp -i PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz -I PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz -o trimmed_PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz -O trimmed_PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz -h Norm_fastp.html

#Alignment 
 bwa mem ./human_ref/hg38.fna trimmed_PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz trimmed_PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz > aligned_tumor.sam

  bwa mem ./human_ref/hg38.fna trimmed_PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz trimmed_PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz > aligned_normal.sam

#samtools 
   samtools view -bS aligned_tumor.sam > tumor.bam
   samtools view -bS aligned_normal.sam > normal.bam
   
   samtools view -H normal.bam
    samtools sort -o normal_sorted.bam normal.bam
    samtools sort -o tumor_sorted.bam tumor.bam

    samtools index normal_sorted.bam
    samtools index tumor_sorted.bam

    samtools faidx hg38.fna

#Strelka
   git clone https://github.com/Illumina/strelka.git
   strelka
    conda install -c bioconda strelka
     cd strelka/
     wget https://github.com/Illumina/strelka/releases/download/v${STRELKA_VERSION}/strelka-${STRELKA_VERSION}.release_src.tar.bz2
      tar -xjf strelka-${STRELKA_VERSION}.release_src.tar.bz2

      ./bin/configureStrelkaSomaticWorkflow.py   --normalBam normal_sorted.bam   --tumorBam tumor_sorted.bam   --referenceFasta ./human_ref/hg38.fna   --runDir strelka_output
    
#Mutation calling 

    tabix -l somatic.snvs.vcf.gz  
    tabix -l somatic.indels.vcf.gz

    bcftools view somatic.snvs.vcf.gz | head

    bcftools mpileup -f ./human_ref/hg38.fna normal_sorted.bam
    bcftools mpileup -f ./human_ref/hg38.fna normal_sorted.bam | bcftools call -mv -Oz -o normal.vcf.gz


    bcftools isec -p output_dir -n=2 ./strelka_output/results/variants/somatic.snvs.vcf.gz normal.vcf.gz

    bg_mutation_rate=$(bcftools stats normal_pass.vcf | grep "number of SNPs" | awk '{print $NF}')

    bcftools query -f '%QUAL\t%DP\n' /data/normal_variants.vcf.gz > /data/normal_quality_depth.tsv


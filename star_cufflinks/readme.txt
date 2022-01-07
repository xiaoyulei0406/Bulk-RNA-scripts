###This package is a set of scripts that perform analysis of RNA-Seq data for STEP 1: alignment

###########################################
###input files:
###########################################

1 Sample infomation files: sample.txt and sample_pairs.txt
Example of sample.txt :
001_13S_S7
001_13S_S1
003_13S_S8
003_13S_S3
004_13S2_S9
004_13S_S4
005_13S_S5
005_S13_S10
006_13S_S11
006_13S8_S6


2 fastq files

3 references files and tools:

ref = '/home/ec2-user/tools/gatk_bundle/b37/fasta/human_g1k_v37.fasta'
STAR='/home/cyu/tools/STAR-2.7.5a/bin/Linux_x86_64/STAR '
CUFFLINKS='/home/ec2-user/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks'
GATK = 'java -jar /home/ec2-user/tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
SAMTOOLS ='/home/cyu/tools/samtools-1.9/samtools'
BEDTOOLS ='/home/ec2-user/anaconda3/bin/bedtools'
PICARD ='java -jar /home/ec2-user/tools/picard.jar'
interval_bed= '/home/ec2-user/tools/gatk_bundle/b37/b37.interval_list_removevirus.bed'
HTseq='/home/ec2-user/anaconda3/bin/htseq-count'

###########################################
###scripts:
###########################################
#pre-quality control

1 run_fastqc.py This script get quality control reports from fastqc
Example command:
python /data/cyu/scripts/rna/run_fastqc.py \
-i /data/cyu/RNA_Bea/fastq/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 1 -e 10 -o /data/cyu/RNA_Bea/

2 run_trim_galore.py (optional) This script uses trim_galore to trim adaptor

Example command:
python /data/cyu/scripts/rna/run_trim_galore.py \
-i /data/cyu/RNA_Bea/fastq/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 10 -e 10 -o /data/cyu/RNA_Bea/


3 run_STAR.py This script run STAR to alignment. We can use raw fastq files if the fastqc is good or we can use trimmed fastq files if we do run_trim_galore.py step.

Example command:
python /data/cyu/scripts/rna/run_STAR.py \
-i /data/cyu/RNA_Bea/fastq/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 1 -e 1 -f raw_fq -o /data/cyu/RNA_Bea/

python /data/cyu/scripts/rna/run_STAR.py \
-i /data/cyu/RNA_Bea/fastq/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 1 -e 1 -f trim_fq -o /data/cyu/RNA_Bea/


##post-quality control
run_rseqc.py: Please try run all the samples in one, it will return a summary pdf file.
Example command:
python /data/cyu/scripts/rna/run_rseqc.py \
-i /data/cyu/RNA_Bea/data/STAR_v1/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 1 -e 10 -o /data/cyu/RNA_Bea/

##Quatification
4 run_htseqcount.py This script uses htseqcount to get gene counts
Example command:
python /data/cyu/scripts/rna/run_htseqcount.py \
-i /data/gbm/rna/data/bam/ -p /data/cyu/gbm_rna/samples.txt -s 1 -e 6 -o /data/cyu/gbm_rna/


5 run_cufflinks.py This script uses cufflinks to calculate FPKM
Example command:
python /data/cyu/scripts/rna/run_cufflinks.py \
-i /data/cyu/gbm_rna/data/bam/ -p /data/cyu/gbm_rna/samples.txt -s 1 -e 2 -o /data/cyu/gbm_rna/





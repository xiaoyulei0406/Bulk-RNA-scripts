#!/bin/bash

###################
# The input parameters of this pipeline:
#(1) R1.fastq.gz, R2.fastq.gz and output folder (Please provide the absolute path of the file)
#(2) Now the pipeline supports species for: human; mouse; turkey; cow
#(3) RNAseq library strandedness (as input for HTseq-count):reverse; yes or no

###################

fastqc=/home/ec2-user/bin/FastQC/fastqc
picard=/home/ec2-user/tools/picard.jar
STAR=STAR
htseq=htseq-count
File1=$1   #need absolute path
File2=$2
OutFolder=$3    #need absolute path
org=$4
strandedness=$5


case $org in
  human)
    STARGenomeFolder=/home/ec2-user/anaconda3/STAR
    GTF=/home/ec2-user/refs/human/gencode.v19.annotation.gtf
    RefFlat=/home/ec2-user/refs/human/refFlat.txt
    rRNA_interval=/home/ec2-user/refs/human/hg19_rRNA18.interval
  ;;
  mouse)
    STARGenomeFolder=/home/ec2-user/refs/mm10
    GTF=/home/ec2-user/refs/mm10/mm10.gtf
    Genome=/home/ec2-user/refs/mm10/mm10.fa
    RefFlat=/home/ec2-user/refs/mm10/mm10_refFlat.txt
    rRNA_interval=/home/ec2-user/refs/mm10/mm10_rRNA.interval
    ;;
  turkey)
      STARGenomeFolder=/home/ec2-user/refs//turkey/
      GTF=/home/ec2-user/refs/turkey/Meleagris_gallopavo.UMD2.94.gtf
      Genome=/home/ec2-user/refs/turkey/Meleagris_gallopavo.UMD2.dna_sm.chromosome.all.fa
      rRNA_interval=/home/ec2-user/refs/turkey/melGal1-rRNA.interval
  ;;
  cow)
    STARGenomeFolder=/home/ec2-user/refs/cow
    GTF=/home/ec2-user/refs/cow/Bos_taurus.UMD3.1.90.gtf
    RefFlat=/home/ec2-user/refs/cow/STAR/refFlat.txt
    rRNA_interval=/home/ec2-user/refs/cow/STAR/bosTau8_RNArepeat0.interval
  ;;
  other)
    STARGenomeFolder=/home/ec2-user/refs/Pseudomonas_aeruginosa
    GTF=/home/ec2-user/refs/Pseudomonas_aeruginosa/Pseudomonas_aeruginosa_PAO1_107.gtf
  ''
esac

echo ============== STARTING ${name} =================

echo Output file be stored at:  ${OutFolder}

i=${File1##*/}   #unabslute path
j=${File2##*/}
name=$( echo $i|cut -d '_' -f 1 )   #file name before _S

if [ ! -d "$OutFolder" ]
then
  mkdir -p ${OutFolder}
fi
cd ${OutFolder}


 #Run fastQC
echo "  Runing fastQC on ${i}"
$fastqc -o ${OutFolder} ${File1}
echo "  Runing fastQC on ${j}"
$fastqc -o ${OutFolder} ${File2}

# Run STAR
#date 
echo "  Runing STAR on ${name}"
$STAR --limitSjdbInsertNsj 1800000 --genomeDir ${STARGenomeFolder} --runThreadN 5 --readFilesCommand zcat --readFilesIn ${File1} ${File2}  --sjdbGTFfile ${GTF} --outFileNamePrefix ${name}_ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotifi --outReadsUnmapped Fastx 

if [[ -a $RefFlat ]] && [ -a "$rRNA_interval" ]
then
#  date 
  echo "  Runing QC on Bam"    #This only calculate statistics from mapped reads. Check mapping precentage first from STAR statistics
  $java -jar ${picard} CollectRnaSeqMetrics \
  I=${name}_Aligned.sortedByCoord.out.bam \
  O=${name}_RNA_Metrics \
  REF_FLAT=${RefFlat} \
  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=${rRNA_interval} \
  CHART=${name}_Aligned.sortedByCoord.out.geneBodyCov.pdf
else
  echo "  REF_FLAT and RIBOSOMAL_INTERVALS are not available. Skipping CollectRnaSeqMetrics"
fi

# Run HT-Seq
echo "  Runing HT-Seq on ${name}"
$htseq -s ${strandedness} -f bam  ${name}_Aligned.sortedByCoord.out.bam ${GTF} > ${name}_HTseq.txt
#date 

echo ==================== ENDING ${name} =======================
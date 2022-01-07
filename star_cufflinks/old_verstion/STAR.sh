fq1=$1
fq2=$2
sampleName=$3
outdir=$4


echo ======================Alignment ===================


/home/cyu/tools/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir /home/cyu/RNASeq/index/ICON/STAR/GRCh37_ensembl/ \
--readFilesIn ${fq1} ${fq2} \
--readFilesCommand zcat \
--outSAMstrandField intronMotif \
--outFileNamePrefix ${outdir}/${sampleName} \
--outReadsUnmapped Fastx \
--twopassMode Basic \
--limitSjdbInsertNsj 1800000 \
--outSAMtype BAM Unsorted \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMattrRGline ID:${sampleName} PL:illumina PU:CCD LIB:${sampleName} SM:${sampleName}

echo ======================bam sort ===================
echo "Samtools start"
/home/cyu/tools/samtools-1.9/samtools sort -@ 8 -n ${outdir}/${sampleName}Aligned.out.bam -o ${outdir}/${sampleName}_name_sorted.bam #for HTSeq
/home/cyu/tools/samtools-1.9/samtools sort -@ 8 ${outdir}/${sampleName}Aligned.out.bam -o ${outdir}/${sampleName}_sorted.bam # for Cufflinks
/home/cyu/tools/samtools-1.9/samtools index -@ 8 ${outdir}/${sampleName}_sorted.bam # for Cufflinks


import os,sys, glob, subprocess
import argparse

STAR='/home/cyu/tools/STAR-2.7.5a/bin/Linux_x86_64/STAR '
genomeindex='/home/cyu/RNASeq/index/ICON/STAR/GRCh37_ensembl/'
SAMTOOLS = '/home/cyu/tools/samtools-1.9/samtools'

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
		description='Run STAR for RNA-Seq in FastQC STAR-2.7.5a, \n\n ex: python run_STAR.py -i <input_dir> -o <out_dir> <patient file> <start> <end>')
	parser.add_argument ('-i', dest='input_dir', help='Input directionary of RNA-Seq fastq files', type=str)
	parser.add_argument ('-o', dest='out_dir', help='Output directionary of STAR results', type=str)
	parser.add_argument ('-p', dest='pt_file', help='Input patient information file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	args = parser.parse_args ()
	return args

def main():
	args = parse_arguments ()
	input_dir = args.input_dir
	outdir = args.out_dir
	pt_file = args.pt_file
	start = args.start
	end = args.end

	sh_dir = outdir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = outdir + "/log/"  # The script log/alignment_start_end.sh will be executed
	subprocess.call ('mkdir -p ' + outdir + "/{sh,log}/", shell=True)

	fout = sh_dir + 'STAR_' + start + '_' + end + '.sh'
	sh = open (fout, 'w')
	pts = open (pt_file)  # note there is a header
	log = log_dir + 'STAR_' + start + '_' + end + '.log'
	lns = pts.readlines ()
	cmd = 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n'
	cmd += 'echo Run STAR for RNA-Seq for RNA-Seq >>' + log + ' 2>&1 &'
	cmd += 'echo STAR Version: STAR-2.7.5a >>' + log + ' 2>&1 &'
	cmd +='\n'
	cmd += 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n\n\n'

	for i in range ((int) (start), (int) (end) + 1):
		print (lns[i])
		tmp = lns[i].strip ('\n').split (',')
		pID = tmp[0]  #patientID
		subprocess.call ('mkdir -p ' + outdir + '/data/STAR/', shell=True)

		fq1 = input_dir + "/" + pID + '_R1_*.gz'
		fq2 = input_dir + "/" + pID + '_R2_*.gz'
		out_dir1 = outdir + '/data/STAR/'

		cmd += STAR + ' --runThreadN 4 --genomeDir ' + genomeindex + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' '
		cmd += '--readFilesCommand zcat '
		cmd += '--outSAMstrandField intronMotif '
		cmd += '--outFileNamePrefix ' + out_dir1 + '/' + pID + '.'
		cmd +=' --outReadsUnmapped Fastx --twopassMode Basic --limitSjdbInsertNsj 1800000 '
		cmd +='--outSAMtype BAM Unsorted --outFilterIntronMotifs RemoveNoncanonical '
		cmd +='--outSAMattrRGline ID:' + pID + ' PL:illumina PU:CCD LIB:' + pID +' SM:' + pID + ' \n'
		cmd += 'echo $(date) STAR is Done. \n\n'

		cmd +='Start bam sort for HTSeq-count \n\n'
		cmd += SAMTOOLS + ' sort -@ 8 -n ' + out_dir1 + '/' + pID + '.Aligned.out.bam -o ' + out_dir1 + '/' + pID+ '_name_sorted.bam \n'
		cmd +='Start bam sort for cufflinks \n\n'
		cmd += SAMTOOLS + ' sort -@ 8 ' + out_dir1 + '/' + pID + '.Aligned.out.bam -o ' + out_dir1 + '/' + pID+ '_sorted.bam \n'
		cmd += 'echo Finished at  $(date)\n\n'
		sh.write (cmd)
	sh.close ()
	pts.close ()
	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)



if __name__ == '__main__':
	main()

'''
python /data/cyu/scripts/rna/run_STAR.py \
-i /data/cyu/RNA_Bea/data/trim_fq/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 1 -e 1 -o /data/cyu/RNA_Bea/

python /data/cyu/scripts/rna/run_STAR.py \
-i /data/cyu/RNA_Bea/data/trim_fq/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 2 -e 10 -o /data/cyu/RNA_Bea/

'''

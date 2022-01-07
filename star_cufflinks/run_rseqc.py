import os,sys, glob, subprocess
import argparse

FASTQC = '/home/cyu/tools/FastQC/fastqc'

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
		description='Run Rseqcfor RNA-Seq in post alignmet, ex: python run_fastqc.py -i /data/cyu/RNA_Bea/fastq/ -p /data/cyu/RNA_Bea/sample.txt -s 1 -e 10 -o /data/cyu/RNA_Bea/')
	parser.add_argument ('-i', dest='input_dir', help='Input directionary of RNA-Seq fastq files', type=str)
	parser.add_argument ('-o', dest='out_dir', help='Output directionary of FASTQC results', type=str)
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

	fout = sh_dir + 'rseqc_' + start + '_' + end + '.sh'
	sh = open (fout, 'w')
	pts = open (pt_file)  # note there is a header
	log = log_dir + 'rseqc_' + start + '_' + end + '.log'
	lns = pts.readlines ()


	cmd = 'echo \'###########################################\' >>' + log + ' 2>&1 & wait; \n'
	cmd += '\n'
	cmd += 'echo Run rseqc for RNA-Seq >>' + log + ' 2>&1 & wait; \n'
	cmd += 'echo rseqc version 3.3.0 >>' + log + ' 2>&1 & wait; \n'
	cmd +='\n'
	cmd += 'echo \'###########################################\' >>' + log + ' 2>&1 & wait; \n'
	cmd += '\n\n\n'

	for i in range ((int) (start)-1, (int) (end) ):
		print (lns[i])
		tmp = lns[i].strip ('\n').split (',')
		pID = tmp[0]  #patientID
		subprocess.call ('mkdir -p ' + outdir + '/data/rseqc/', shell=True)

		srtbam = input_dir + "/" + pID + 'sorted.bam'
		alnbam= input_dir + "/" + pID + '.Aligned.out.bam'
		outdir1 = outdir + '/data/rseqc/'
		cmd += 'echo Run rseqc geneBody_coverage.py \n'
		cmd += '/home/ec2-user/anaconda3/bin/python /home/ec2-user/anaconda3/bin/geneBody_coverage.py '
		cmd += ' -r /home/cyu/RNASeq/refs/Annotations/hg19.HouseKeepingGenes.rmchr.bed '
		cmd += ' -i ' + srtbam
		cmd += ' -o ' + outdir1 + '/' + pID + ' & wait; \n'

		cmd += 'echo Run rseqc read_distribution.py \n'
		cmd += '/home/ec2-user/anaconda3/bin/python /home/ec2-user/anaconda3/bin/read_distribution.py '
		cmd += ' -i ' + alnbam
		cmd += ' -r /home/cyu/RNASeq/refs/Annotations/hg19_Ensembl_gene_removechr.bed > '
		cmd += outdir1 + '/' + pID + '_read_distribution.txt & wait; \n'

		cmd += 'echo $(date) rseqc is Done. \n\n'
		sh.write (cmd)
	sh.close ()
	pts.close ()
	cmd = 'nohup bash ' + fout + '>>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)



if __name__ == '__main__':
	main()

'''
python /data/cyu/scripts/rna/run_rseqc.py \
-i /data/cyu/RNA_Bea/data/STAR_v1/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 1 -e 1 -o /data/cyu/RNA_Bea/
'''

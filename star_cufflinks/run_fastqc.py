import os,sys, glob, subprocess
import argparse

FASTQC = '/home/cyu/tools/FastQC/fastqc'

def parse_arguments():
	'''
	Gives options
	'''
	parser = argparse.ArgumentParser (
		description='Run FASTQC for RNA-Seq in FastQC v0.11.8, ex: python run_fastqc.py -i /data/cyu/RNA_Bea/fastq/ -p /data/cyu/RNA_Bea/sample.txt -s 1 -e 10 -o /data/cyu/RNA_Bea/')
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

	fout = sh_dir + 'fastqc_' + start + '_' + end + '.sh'
	sh = open (fout, 'w')
	pts = open (pt_file)  # note there is a header
	log = log_dir + 'fastqc_' + start + '_' + end + '.log'
	lns = pts.readlines ()


	cmd = 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n'
	cmd += 'echo Run FASTQC for RNA-Seq >>' + log + ' 2>&1 &'
	cmd += 'echo FASTQC version v0.11.8 >>' + log + ' 2>&1 &'
	cmd +='\n'
	cmd += 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n\n\n'

	for i in range ((int) (start)-1, (int) (end) ):
		print (lns[i])
		tmp = lns[i].strip ('\n').split (',')
		pID = tmp[0]  #patientID
		subprocess.call ('mkdir -p ' + outdir + '/data/fastqc/', shell=True)

		fq1 = input_dir + "/" + pID + '_R1*.gz'
		fq2 = input_dir + "/" + pID + '_R2*.gz'
		outdir1 = outdir + '/data/fastqc/'

		cmd += FASTQC + ' -t 8 '
		cmd += ' -o ' + outdir1 + ' ' + fq1 + '& wait; \n'
		cmd += FASTQC + ' -t 8 '
		cmd += ' -o ' + outdir1 + ' ' + fq2 + '& wait; \n'

		cmd += 'echo $(date) fastqc is Done. \n\n'
		cmd += 'echo Finished at  $(date)\n\n'
		sh.write (cmd)
	sh.close ()
	pts.close ()
	cmd = 'nohup bash ' + fout + '>>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)



if __name__ == '__main__':
	main()

'''
python /data/cyu/scripts/rna/run_fastqc.py \
-i /data/cyu/RNA_Bea/fastq/ \
-p /data/cyu/RNA_Bea/sample.txt \
-s 1 -e 10 -o /data/cyu/RNA_Bea/
'''

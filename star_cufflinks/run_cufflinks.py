import argparse, sys, re, subprocess
CUFFLINKS='/home/ec2-user/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks'
import argparse
def parse_arguments():
	'''
	Gives options
	'''
	# Parse the user supplied arguments
	parser = argparse.ArgumentParser (description='Running cufflinks')
	parser.add_argument ('-i', dest='input_dir', help='The input path')
	parser.add_argument ('-p', dest='pt_file', help='Tnput patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-l', dest='library_types', nargs='?', default='fr-unstranded',
						 help='Supported library types:ff-firststrand,ff-secondstrand,ff-unstranded,fr-firststrand,fr-secondstrand,fr-unstranded (default),transfrags')

	args = parser.parse_args ()
	return args

def main():

	args = parse_arguments ()
	input_dir = args.input_dir
	pt_file = args.pt_file
	start = args.start
	end = args.end
	out_dir = args.out_dir
	library_types =args.library_types

	sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed
	subprocess.call ('mkdir -p ' + out_dir + "/{sh,log}/", shell=True)
	out_dir = out_dir + "/data/cufflinks/"
	subprocess.call ('mkdir -p ' + out_dir, shell=True)

	fout = sh_dir + 'cufflinks_' + start + '_' + end + '.sh'
	sh = open (fout, 'w')
	pts = open (pt_file)  # note there is a header
	log = log_dir + 'cufflinks_' + start + '_' + end + '.log'
	lns = pts.readlines ()

	cmd = 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n'
	cmd += 'echo Run Cufflinks for RNA-Seq for RNA-Seq >>' + log + ' 2>&1 &'
	cmd += 'echo cufflinks Version: cufflinks-2.2.1.Linux_x86_64 >>' + log + ' 2>&1 &'
	cmd += '\n'
	cmd += 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n\n\n'

	for i in range ((int) (start)-1, (int) (end)):
		tmp = lns[i].strip ("\n").split ('\t')
		inputdir = input_dir
		sample_id = tmp[0]
		out_dir1= out_dir + sample_id + "/"
		subprocess.call ('mkdir -p ' + out_dir1, shell=True)
		bam = inputdir + sample_id + '_sorted.bam'

		cmd += 'echo 1.Running Cufflinks commandline on patient ' + sample_id + '...\n\n'

		cmd += CUFFLINKS + ' -p 1 -o ' + out_dir1 + ' --multi-read-correct '
		cmd += ' --library-type '+ library_types +' -G /home/cyu/RNASeq/refs/Annotations/Homo_sapiens.GRCh37.75.gtf '
		cmd += bam + ' \n'

		cmd += 'echo $(date) cufflinks is done. \n\n'
		sh.write (cmd)
	sh.close ()
	pts.close ()
	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	subprocess.call (cmd, shell=True)

if __name__ == '__main__':
	main ()

'''
python /data/cyu/scripts/rna/run_cufflinks.py \
-i /data/cyu/gbm_rna/data/bam/ -p /data/cyu/gbm_rna/samples.txt -s 1 -e 2 -o /data/cyu/gbm_rna/

python /data/cyu/scripts/rna/run_cufflinks.py \
-i /data/cyu/RNA_Bea/data/STAR_v1/ -p /data/cyu/RNA_Bea/sample.txt -s 1 -e 1 -o //data/cyu/RNA_Bea/


'''

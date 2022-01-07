import argparse, sys, re, subprocess

ref = '/home/ec2-user/tools/gatk_bundle/b37/fasta/human_g1k_v37.fasta'
GATK = 'java -jar /home/ec2-user/tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
GATKBUNDLE = '/home/ec2-user/tools/gatk_bundle/'
BWA = '/home/ec2-user/anaconda3/bin/bwa'
SAMTOOLS = '/home/cyu/tools/samtools-1.9/samtools'
BEDTOOLS = '/home/ec2-user/anaconda3/bin/bedtools'
PICARD = 'java -jar /home/ec2-user/tools/picard.jar'
g1000_indel = '/home/ec2-user/tools/gatk_bundle/b37/1000G_phase1.indels.b37.vcf'
mills = '/home/ec2-user/tools/gatk_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf'
interval_bed = '/home/ec2-user/tools/gatk_bundle/b37/b37.interval_list_removevirus.bed'
STAR='/home/cyu/tools/STAR-2.7.5a/bin/Linux_x86_64/STAR'
genomeindex='/home/cyu/RNASeq/index/ICON/STAR/GRCh37_ensembl/'
HTseq='/home/ec2-user/anaconda3/bin/htseq-count'
RSEM='/home/ec2-user/anaconda3/bin/rsem-calculate-expression'



def main(args):
	# Parse the user supplied arguments
	parser = argparse.ArgumentParser (description='Ziopharm alignment process')
	parser.add_argument ('-i', dest='input_dir', help='The input path')
	parser.add_argument ('-p', dest='pt_file', help='Tnput patient file')
	parser.add_argument ('-s', dest='start', help='The start job id')
	parser.add_argument ('-e', dest='end', help='The end job id')
	parser.add_argument ('-o', dest='out_dir', help='The output path')
	parser.add_argument ('-l', dest='stranded', nargs='?', default='no',
						 help='Supported library types:yes,no,reverse')

	args = parser.parse_args ()
	input_dir = args.input_dir
	pt_file = args.pt_file
	start = args.start
	end = args.end
	out_dir = args.out_dir
	stranded=args.stranded

	sh_dir = out_dir + "/sh/"  # The script sh/htseq_start_end.sh will be executed
	log_dir = out_dir + "/log/"  # The script log/htseq_start_end.sh will be executed
	subprocess.call ('mkdir -p ' + out_dir + "/{sh,log}/", shell=True)
	out_dir = out_dir + "/data/htseq/"

	subprocess.call ('mkdir -p ' + out_dir, shell=True)

	fout = sh_dir + 'htseq_' + start + '_' + end + '.sh'
	sh = open (fout, 'w')
	pts = open (pt_file)  # note there is a header
	log = log_dir + 'htseq_' + start + '_' + end + '.log'
	lns = pts.readlines ()
	cmd = 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n'
	cmd += 'echo Run htseq-count for RNA-Seq >>' + log + ' 2>&1 &'
	cmd += 'echo htseq-count Version: 0.12.4 >>' + log + ' 2>&1 &'
	cmd += '\n'
	cmd += 'echo \'###########################################\' >>' + log + ' 2>&1 &'
	cmd += '\n\n\n'
	for i in range ((int) (start)-1, (int) (end)):
		tmp = lns[i].strip ("\n").split ('\t')
		inputdir = input_dir
		sample_id = tmp[0]

		bam = inputdir + sample_id + '_name_sorted.bam'

		cmd = 'echo starting expression at time `date` \n\n'
		cmd += 'echo 1.Running HTSeq-count commandline on patient ' + sample_id + '...\n\n'

		cmd += HTseq + ' -s '+ stranded +' -r name --additional-attr gene_name --additional-attr transcript_id '
		cmd += '-f bam ' + bam
		cmd += ' /home/cyu/RNASeq/refs/Annotations/Homo_sapiens.GRCh37.75.gtf > '
		cmd += out_dir + '/' + sample_id +'_HTSeq.txt \n'

		cmd += 'echo Finished HTSeq-count \n\n'
		sh.write (cmd)
	sh.close ()
	pts.close ()

	cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
	print ('Submitting Alignment jobs to process ' + sample_id + '\n\n')
	subprocess.call (cmd, shell=True)


if __name__ == '__main__':
	main (sys.argv)

'''
python /data/cyu/scripts/rna/run_htseqcount.py \
-i /data/gbm/rna/data/bam/ -p /data/cyu/gbm_rna/samples.txt -s 1 -e 6 -o /data/cyu/gbm_rna/

python /data/cyu/scripts/rna/run_htseqcount.py \
-i /data/cyu/RNA_Bea/data/STAR_v1/ -p /data/cyu/RNA_Bea/sample.txt -s 1 -e 1 -o /data/cyu/RNA_Bea/

'''

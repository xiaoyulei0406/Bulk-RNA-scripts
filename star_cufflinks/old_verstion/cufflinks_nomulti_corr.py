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
CUFFLINKS='/home/ec2-user/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks'


def main(args):
    # Parse the user supplied arguments
    parser = argparse.ArgumentParser (description='Ziopharm alignment process')
    parser.add_argument ('-i', dest='input_dir', help='The input path')
    parser.add_argument ('-p', dest='pt_file', help='Tnput patient file')
    parser.add_argument ('-s', dest='start', help='The start job id')
    parser.add_argument ('-e', dest='end', help='The end job id')
    parser.add_argument ('-o', dest='out_dir', help='The output path')

    args = parser.parse_args ()
    input_dir = args.input_dir
    pt_file = args.pt_file
    start = args.start
    end = args.end
    out_dir = args.out_dir

    sh_dir = out_dir + "/sh/"  # The script sh/alignment_start_end.sh will be executed
    log_dir = out_dir + "/log/"  # The script log/alignment_start_end.sh will be executed
    subprocess.call ('mkdir -p ' + out_dir + "/{sh,log}/", shell=True)
    out_dir = out_dir + "/data/cufflinks/no_multi_corr/"
    subprocess.call ('mkdir -p ' + out_dir, shell=True)

    fout = sh_dir + 'cufflinks_no_multi_corr_' + start + '_' + end + '.sh'
    sh = open (fout, 'w')
    pts = open (pt_file)  # note there is a header
    log = log_dir + 'cufflinks_no_multi_corr_' + start + '_' + end + '.log'
    lns = pts.readlines ()

    for i in range ((int) (start), (int) (end) + 1):
        tmp = lns[i].strip ("\n").split ('\t')
        inputdir = input_dir
        sample_id = tmp[0]
        out_dir1= out_dir + sample_id + "/"
        subprocess.call ('mkdir -p ' + out_dir1, shell=True)
        bam = inputdir + sample_id + '_sorted.bam'

        cmd = 'echo starting alignment at time `date` \n\n'
        cmd += 'echo 1.Running STAR commandline on patient ' + sample_id + '...\n\n'

        cmd += CUFFLINKS + ' -p 8 -o ' + out_dir1 + ' --max-bundle-length 10000000 '
        cmd += ' --library-type fr-firststrand -G /home/cyu/RNASeq/refs/Annotations/Homo_sapiens.GRCh37.75.gtf '
        cmd += bam + ' \n'

        cmd += 'echo $(date) cufflinks is done. \n\n'
        sh.write (cmd)
    sh.close ()
    pts.close ()

    cmd = 'nohup bash ' + fout + '>' + log + ' 2>&1 &'
    print ('Submitting Alignment jobs to process ' + sample_id + '\n\n')
    subprocess.call (cmd, shell=True)


if __name__ == '__main__':
    main (sys.argv)

'''
python /data/cyu/scripts/rna/cufflinks_nomulti_corr.py \
-i /data/cyu/gbm_rna/data/bam/ -p /data/cyu/gbm_rna/samples.txt -s 3 -e 4 -o /data/cyu/gbm_rna/

python /data/cyu/scripts/rna/cufflinks_nomulti_corr.py \
-i /data/cyu/gbm_rna/data/bam/ -p /data/cyu/gbm_rna/samples.txt -s 5 -e 6 -o /data/cyu/gbm_rna/

python /data/cyu/scripts/rna/cufflinks_nomulti_corr.py \
-i /data/cyu/gbm_rna/data/bam/ -p /data/cyu/gbm_rna/samples.txt -s 1 -e 2 -o /data/cyu/gbm_rna/


'''


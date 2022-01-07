#source("http://bioconductor.org/biocLite.R")
#biocLite("ensembldb")
#biocLite("EnsDb.Hsapiens.v79")
#Rscript  ./sym2id.r -i /Users/chunleiyu/Work/data/rnaseq_p03012018/DE/test-sigEnsemble.txt -o /Users/chunleiyu/Work/data/rnaseq_p03012018/GO/sym2id.txt
suppressMessages(library("optparse"))

optionList = list(make_option(c("-i", "--input"), type="character", default=NULL, help="count input: file or working folder.\n(1)input file needs the first header line(delimited by \",\");\nor (2) the file in the working folder does not need first head line and count number is in the second column]"),
                  make_option(c("-o","--output"),type="character",default=NULL,help="output folder"));

opt_parser = OptionParser(option_list=optionList)
opt = parse_args(opt_parser)
suppressMessages(library(ensembldb))
suppressMessages(library(EnsDb.Hsapiens.v79))
genes=read.table(opt$input,sep="\n",header=F,stringsAsFactors = F)
genes<-genes[,1]
#genes(EnsDb.Hsapiens.v79, filter=list(GenenameFilter(genes),GeneIdFilter("ENSG", "startsWith")), return.type="data.frame", columns=c("gene_id"))
f<-genes(EnsDb.Hsapiens.v79, filter=list(GenenameFilter(genes),GeneIdFilter("ENSG", "startsWith")), return.type="data.frame", columns=c("gene_id"))[,1]
write.table(f,opt$output,quote = F, row.names = F,col.names = F)


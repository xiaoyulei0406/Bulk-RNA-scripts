#!/usr/bin/env Rscript
##Rscirpt
suppressMessages(library("optparse"))

optionList = list(make_option(c("-i", "--input"), type="character", default=NULL, help="count input: file or working folder.\n(1)input file needs the first header line(delimited by \",\");\nor (2) the file in the working folder does not need first head line and count number is in the second column]"),
                  make_option(c("-p", "--prefix"), type="character", default="output", help="output prefix [default= %default]"),
                  make_option(c("-m","--meta"),type="character",default=NULL,help="metadata file [the format of the file is shown in \"metadata.txt\"]"),
                  make_option(c("-d","--dir"),type="character",default=NULL,help="output folder"),
                  make_option(c("-r","--ref"),type="character",default=NULL,help="reference condition [e.g. \"control\",which must be shwon in your condition level in the metadata file]"),
                  make_option(c("-t", "--thread"), type="integer", default="1", help="number of threads.1: single thread [default= %default]"));

opt_parser = OptionParser(option_list=optionList)
opt = parse_args(opt_parser)

if (is.null(opt$input)||is.null(opt$dir)||is.null(opt$meta)||is.null(opt$ref)){
  print_help(opt_parser)
  stop("Parameter errors!", call.=FALSE)
}

#.libPaths("/home/ec2-user/R/x86_64-pc-linux-gnu-library/3.3/")
options(java.parameters = "-Xmx8000m")
#.libPaths() #show the lib path for debug purpose 

suppressMessages(library(BiocParallel))
suppressMessages(library("reshape2"))
suppressMessages(library("gplots"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))
suppressMessages(library("DESeq2") )
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library("pheatmap"))
suppressMessages(library('xlsx'))
suppressMessages(library('png'))

readCount<- function(line)
{
  temp = read.csv(line[2],row.names =1,sep="\t")
  return(temp[,1])
}

InputFileCount=opt$input
InputFileMeta=opt$meta
outputPrefix=paste0(opt$dir,"/",opt$prefix)
Metadata=read.csv(file=InputFileMeta,sep="\t",header=T,check.names=FALSE,row.names=1)

if (file_test("-f",InputFileCount)){
  inputCount=read.csv(file=InputFileCount,header=T,check.names=FALSE)
  inputCount = read.csv(file=InputFileCount,header=T,check.names=FALSE,row.names=colnames(inputCount)[1])
  inputCount<-inputCount[, rownames(Metadata)]
  
  ddspre <- DESeqDataSetFromMatrix(countData = inputCount,colData=Metadata,design=~condition)
}else{
  sampleTable <- data.frame(sampleName=Metadata[,1],fileName=Metadata[,1],condition=Metadata[,2])
  
  ddspre <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory=opt$input,design=~condition)
  meta=read.table(file=InputFileMeta,header = T,sep="\t",check.names=FALSE)
  print (meta)
  #newFiles<-outer(Metadata[,1],InputFileCount,FUN=function(r,c) paste(trimws(c),r,sep="/"))
  originalDir<-getwd()
  setwd(InputFileCount)
  inputCount = apply(meta,1,readCount)
  inputCount = as.data.frame(inputCount)
  colnames(inputCount) = meta[,1]
  setwd(originalDir)
}


ddspre$condition <-relevel(ddspre$condition,ref=opt$ref)
if (opt$thread==1){
  dds<-DESeq(ddspre)
  resultsDDS<-results(dds)
}else{
  dds<-DESeq(ddspre,parallel=TRUE, BPPARAM=MulticoreParam(opt$thread))
  resultsDDS<-results(dds,parallel=TRUE)
}

#resultsDDS<-results(dds)
orderResultsDDS<-resultsDDS[order(resultsDDS$padj),]
write.xlsx(orderResultsDDS, file=paste0(outputPrefix,"-AnalysisResult.xlsx"), sheetName="DE_All")

Sig<-subset(resultsDDS, resultsDDS$pvalue < 0.05)
if (nrow(Sig)>0){
  OrderSig = Sig[order(Sig$padj),]
  write.xlsx(OrderSig,file=paste0(outputPrefix,"-AnalysisResult.xlsx"),sheetName="DE_Significant", append =TRUE)
  write(rownames(OrderSig),file=paste0(outputPrefix,"-sigEnsemble.txt"))
}





rld=rlogTransformation(ddspre)
vsd <- varianceStabilizingTransformation(ddspre)

## plot count of each sample ###
pseudoCount = log2(inputCount + 1)
pseudoCount = as.data.frame(pseudoCount)
df = melt(pseudoCount) # reshape the matrix
png(file=paste0(outputPrefix,"-BetweenSampleDis.png"),type="cairo")
ggplot(df,aes(x=variable,y=value))+ geom_boxplot(fill="pink") +  ylab(expression(log[2](count+ 1)))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

##plot heat map
png(file=paste0(outputPrefix,"-heatmap.png"),type="cairo")
distsRL <- dist(t(assay(vsd))) # Calculate distances using transformed (and normalized) counts

mat <- as.matrix(distsRL) # convert to matrix
#mat
rownames(mat) <- colnames(mat) <- with(colData(ddspre), paste(condition,Metadata[,1],sep=":")) # set rownames in the matrix
colnames(mat) = NULL # remove column names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # set colors
pheatmap(mat,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors)

dev.off()



#p-value plot
png(file=paste0(outputPrefix,"-pval.png"),type="cairo")
pvalueDf<-data.frame(pvalue=resultsDDS$pvalue)
ggplot(pvalueDf, aes(x=pvalue)) + geom_histogram(color="black", fill="blue",binwidth=.05)+xlab("p-values")
dev.off()

#MA-plot
png(file=paste0(outputPrefix,"-DE_MAplot.png"),type="cairo")
plotMA(resultsDDS, ylim=c(-7,7))
dev.off()

#PCA
png(file=paste0(outputPrefix,"-DE_PCA.png"),type="cairo")
plotPCA(rld, intgroup=c("condition"))
dev.off()

#volcano plot
png(file=paste0(outputPrefix,"-DE_VolcanoPlot.png"), type='cairo')
pvalue.cutoff<-.01 # threshold on the adjust p-value
log2FolderChange.cutoff<-2 #threshold on the log2 folder change

sigDDS<-subset(resultsDDS, resultsDDS$pvalue<pvalue.cutoff & abs(resultsDDS$log2FoldChange)>log2FolderChange.cutoff)
with(resultsDDS, plot(resultsDDS$log2FoldChange, -log10(resultsDDS$pvalue), pch=20, main="DE_VolcanoPlot",col="black"))
with(sigDDS, points(sigDDS$log2FoldChange, -log10(sigDDS$pvalue), pch=20, col="red"))
dev.off()

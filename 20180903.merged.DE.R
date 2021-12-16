#args = commandArgs(trailingOnly=T)
setwd('/vol1/Working/YanLab/Meiling/Project_Mz424/20180903.Bowtie2Deseq2')

# annotate the samples manually with meiling
sampleInfo = read.csv('sampleInfo.csv', comment='#', row.names=1)
sele = sampleInfo$pick=='Y'
sampleInfo=sampleInfo[sele,]
# add metaTissue to sampleInfo
#metaTissue=as.character(sampleInfo$tissue)
#metaTissue[seq(1,nrow(sampleInfo),by=2)] = metaTissue[seq(2,nrow(sampleInfo),by=2)]
#write.csv(sampleInfo, 'sampleInfo.csv')

# read merged gene count from the four batches, rename and align the
# samplenames with sampleInfo
raw = read.csv('mergedGeneCount.csv', row.names=1, check.names=F)
new_names = sub('M10_28_[0-9]+_', 'M10_',colnames(raw))
rownames(sampleInfo) %in% new_names
colnames(raw) = new_names
geneCount = raw[,rownames(sampleInfo)]
dim(geneCount)
dim(sampleInfo)

# minimal gene filter
sele = rowSums(geneCount)>0
geneCount = geneCount[sele,]
write.csv(geneCount, 'geneCount.csv')

# deseq2 normalization and plot
require(DESeq2)
colData = sampleInfo
countData = geneCount
all(colnames(countData) == rownames(colData))
rownames(colData)=paste0(colData$patientIndex, '.', colData$metaStatus)
colnames(countData) = rownames(colData)
#sele = which(colData$metaTissue!='liver')
# setup obj from count data
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~metaStatus)
vst = varianceStabilizingTransformation(dds, blind=T)
expr = assay(vst)
write.csv(expr, 'expr.vst.csv')

library("RColorBrewer")
require('ggplot2')
require('pheatmap')

# color by metaTissue
data <- plotPCA(vst, intgroup=c('metaTissue', 'metaStatus'), ntop=5000, returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=metaTissue, shape=metaStatus)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance"))# +
dev.copy(pdf, 'plotPCA_top5000.pdf'); dev.off()
###################################################################

## read sampleInfo
sampleInfo = read.csv('sampleInfo.csv', row.names=1)
sampleInfo$tissue = relevel(sampleInfo$tissue, 'breast')
sampleInfo$metaTissue = as.factor(sampleInfo$metaTissue)
sampleInfo$metaStatus = relevel(sampleInfo$metaStatus, 'P')
sampleInfo$patientIndex = as.factor(sampleInfo$patientIndex)
sampleInfo$Batch=as.factor(sampleInfo$Batch)
colData = sampleInfo


# DEGs
# unpaired
dds = DESeq(dds)
group = colData$metaStatus
y = as.data.frame(expr)
x = results(dds)
x = x[order(x$stat),]
xx = cbind(x, y[rownames(x), c(which(group=='P'), which(group=='M'))])
write.csv(xx, 'DE_unpaired.M_vs_P.csv')
sele = which(xx$padj<.05)
res = xx[sele,]
write.csv(res, 'DE_unpaired.M_vs_P.significant.csv')

## paired
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData,
        design=~ patientIndex + metaStatus)
dds = DESeq(dds)
y = as.data.frame(expr)
#x = results(dds)
x = results(dds, contrast=c('metaStatus', 'M', 'P'))
x = x[order(x$stat),]
#x = tmp
xx = cbind(x, y[rownames(x),])
write.csv(xx, 'DE_paired.M_vs_P.csv')
sele = which(xx$padj<.05)
res = xx[sele,]
write.csv(res, 'DE_paired.M_vs_P.significant.csv')

DESeq2::plotMA(x, alpha=0.05, ylim=c(-4,3))
dev.copy(pdf, 'DE_paired.M_vs_P.plotMA.pdf'); dev.off()


#==============================================================
expr = as.data.frame(expr)
## tissue specific
ddsTissue = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~ patientIndex + tissue)
ddsTissue = DESeq(ddsTissue)
tissueContrasts = list(
    brain=c('brain', 'breast'),
    liver=c('liver', 'breast'),
    lung=c('lung', 'breast'),
    ovary=c('ovary', 'breast'))

contrasts = tissueContrasts
groupCol = 'tissue'; group = colData[,groupCol]
dds = ddsTissue
y = expr

de = list()
for (i in names(contrasts)) {
    cc = contrasts[[i]]
    res = results(ddsTissue, contrast=c(groupCol, cc))
    DESeq2::plotMA(res, alpha=.05, ylim=c(-3,3), main=i)
    #identify(res$baseMean, res$log2FoldChange, labels=rownames(res))
    dev.copy(pdf, paste0(i, '.plotMA.pdf')); dev.off()
    x = as.data.frame(res)
    x = x[order(x$stat),]
	seleCols = which(colData$metaTissue==cc[1])
    xx = cbind(x, y[rownames(x), seleCols])
    de[[i]] = xx
}

for (i in names(de)) {
    curr = de[[i]]
    write.csv(curr, paste0(i, '.StatWithVst.csv'))
    curr = as.data.frame(curr)
    sele = which(curr$padj<.05)
    write.csv(curr[sele,], paste0(i, '.StatWithVst.padj_05.csv'))
}

## compare of contrasts
## plot venndiagrams: p=.05, fc=1.5
source ('~/code/working/venndiagram.R')
runVenns(de)
de[['paired']] = read.csv('DE_paired.M_vs_P.csv', row.names=1)
runVenns(de[c('paired', 'liver')], 'paired_liver')
runVenns(de[c('paired', 'brain')], 'paired_brain')
runVenns(de[c('paired', 'ovary')], 'paired_ovary')
runVenns(de[c('paired', 'lung')], 'paired_lung')



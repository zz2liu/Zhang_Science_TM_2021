#args = commandArgs(trailingOnly=T)
sampleInfo = read.csv('paired_meta.csv', row.names=2)
expr = read.csv('/vol/Working/NguyenLab/TCGA_BRCA/20150924.RnaseqV2.Level3/gene_vst.csv',
	 row.names=1, check.names=F)
exprPaired = expr[,rownames(sampleInfo)]

meilineSig = read.csv('/vol/Working/YanLab/Meiling/Project_Mz424/Bowtie2Deseq2/DE_paired.M_vs_P.significant.csv', row.names=1)
sigGenes = rownames(meilineSig)

tcgaGenes = sapply(rownames(expr), function(x) {strsplit(x, split='|', fixed=T)[[1]][1]})
sele = match(sigGenes, tcgaGenes)
tmp = exprPaired[sele,]
priCols = seq(1,ncol(exprPaired), by=2)
metCols = seq(2, ncol(exprPaired), by=2)
diffs = tmp[,metCols] - tmp[,priCols]
diffs$log2FC = rowMeans(diffs)
colnames(diffs) = substr(colnames(diffs), 1, 12)
tcgaExpr = tmp
group = rep(c(1,2), ncol(tcgaExpr)-1)
i = 1
x = as.numeric(tcgaExpr[i,])
.testOne = function(x) {
	tmp=t.test(x[group==2], x[group==1], paired=T) 
	c(tmp$statistic, pval=tmp$p.value,
	logFC=mean(x[group==2]-x[group==1]),
	baseMean=mean(x))
}

tcgaExpr = as.matrix(tcgaExpr)
res = t(apply(tcgaExpr, 1, .testOne)) 


#########
raw = read.csv('/vol/Working/NguyenLab/TCGA_BRCA/20150924.RnaseqV2.Level3/gene_count.csv',
	 row.names=1, check.names=F)
tmp = raw[,rownames(sampleInfo)]
geneCount = tmp[rowSums(tmp)>0,]
sampleInfo$patient = substr(rownames(sampleInfo), 1, 12)
write.csv(sampleInfo, 'sampleInfo.csv')
write.csv(geneCount, 'geneCount.csv')

require(DESeq2)
tmp = as.matrix(geneCount)
countData = matrix(as.integer(tmp), nrow=nrow(tmp), dimnames=dimnames(tmp))
colData = sampleInfo
colData$MetaStatus = relevel(colData$MetaStatus, 'Primary Tumor')
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~patient+MetaStatus)
vst = varianceStabilizingTransformation(dds, blind=T)
expr = assay(vst)
write.csv(expr, 'expr.vst.csv')

library("RColorBrewer")
require('ggplot2')
require('pheatmap')

#data <- plotPCA(vst, intgroup=c('metaTissue', 'metaStatus'), ntop=500, returnData=TRUE)
data <- plotPCA(vst, intgroup=c('patient', 'MetaStatus'), ntop=5000, returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=patient, shape=MetaStatus)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance"))# +
dev.copy(pdf, 'plotPCA_top5000.pdf'); dev.off()

sampleDists = dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(ncol(vst), 'Blues')))(255)
require('gplots')
colors = greenred(99)
labels = paste0(substr(sampleInfo$MetaStatus, 1,4), '.', substring(sampleInfo$patient, 6))
dimnames(sampleDistMatrix) = list(labels, labels)
pheatmap(sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors)
dev.copy(pdf, 'sampleClustering.pdf'); dev.off()

## DEGs
dds = DESeq(dds)
y = as.data.frame(expr)
x = results(dds)
x = x[order(x$stat),]
xx = cbind(x, y[rownames(x),])
xx$geneSymbol = sapply(rownames(xx), function(x) { strsplit(x, split='|', fixed=T)[[1]][1] })
write.csv(xx, 'DE_paired.M_vs_P.csv')
sele = which(xx$padj<.05)

res = xx[sele,]
write.csv(res, 'DE_paired.M_vs_P.significant.csv')

## compare DEGs
#read and filter indepently
mRes = read.csv('/vol/Working/YanLab/Meiling/Project_Mz424/Bowtie2Deseq2/DE_paired.M_vs_P.csv',
    row.names=1)
mRes = mRes[!is.na(mRes$padj),]

tRes = read.csv('DE_paired.M_vs_P.csv', row.names=1)
tRes = tRes[!is.na(tRes$padj),]

#align the two DEGs
commGenes = intersect(rownames(mRes), tRes$geneSymbol)
#pick most significant geneId for each symbol
tOrd = order(-abs(tRes$stat))
tmp = tRes[tOrd,]
sele = match(commGenes, tRes$geneSymbol)
tRes_ = tRes[sele,]
mRes_ = mRes[commGenes,]

write.csv(tRes_, 'tcgaDe.csv')
write.csv(mRes_, 'myDe.csv')

## compare stats by correlation
rownames(tRes_) = rownames(mRes_)
cor.test(tRes_$stat, mRes_$stat)
cor.test(tRes_$stat, mRes_$stat, method='spearman')
require(LSD)
heatscatter(tRes_$stat, mRes_$stat)
abline(lm(mRes_$stat~tRes_$stat))
title(main='cor.test p-value < 2.2e-16')
dev.copy(pdf, 'Stat_comparision.pdf'); dev.off()

## compare sig genes by vennDiagram
source('~/code/working/venndiagram.R')
dfs = list(m=mRes_[,1:6], t=tRes_[,1:6])
runVenns(dfs)
overlapping.hyper(39, 668+39, 39+215, 16000)[[1]]
overlapping.hyper(5, 236, 65, 16000)[[1]]
overlapping.hyper(34, 34+437, 34+155, 16000)[[1]]



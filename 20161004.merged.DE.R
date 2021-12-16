#args = commandArgs(trailingOnly=T)
setwd('/vol1/Working/YanLab/Meiling/Project_Mz424/Bowtie2Deseq2')

# annotate the samples manually with meiling
sampleInfo = read.csv('../sampleInfo.csv', row.names=1)
sele = sampleInfo$pick=='Y'
sampleInfo=sampleInfo[sele,]
# add metaTissue to sampleInfo
metaTissue=as.character(sampleInfo$tissue)
metaTissue[seq(1,nrow(sampleInfo),by=2)] = metaTissue[seq(2,nrow(sampleInfo),by=2)]
write.csv(sampleInfo, 'sampleInfo.csv')

## read sampleInfo
sampleInfo = read.csv('sampleInfo.csv', row.names=1)
sampleInfo$tissue = relevel(sampleInfo$tissue, 'breast')
sampleInfo$metaTissue = as.factor(sampleInfo$metaTissue)
sampleInfo$metaStatus = relevel(sampleInfo$metaStatus, 'P')
sampleInfo$patientIndex = as.factor(sampleInfo$patientIndex)
sampleInfo$Batch=as.factor(sampleInfo$Batch)

## read the count data and aligned to sampleInfo
tmp2 = read.delim('../Batch2/bowtie2hg19/Sample.featureCounts.count_', row.names=1)
tmp3 = read.csv('../Batch3/bowtie2hg19/featureCounts.csv', row.names=1)
tmp = cbind(tmp2, tmp3)
geneCount = tmp[,rownames(sampleInfo)]
nz = rowSums(geneCount)>0
sum(nz) #20K
geneCountNz = geneCount[nz,]
head(geneCountNz)
write.csv(geneCountNz, 'geneCountNz.csv')
geneCountNz = read.csv('geneCountNz.csv', row.names=1)


## Deseq2
require(DESeq2)
colData = sampleInfo
countData = geneCountNz
all(colnames(countData) == rownames(colData))
rownames(colData)=paste0(colData$patientIndex, '.', colData$metaStatus)
colnames(countData) = rownames(colData)
sele = which(colData$metaTissue!='liver')
# setup obj from count data
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~metaStatus)
vst = varianceStabilizingTransformation(dds, blind=T)
expr = assay(vst)
write.csv(expr, 'expr.vst.csv')

library("RColorBrewer")
require('ggplot2')
require('pheatmap')

#PCA
#plotPCA(vst, intgroup=c('patientIndex','metaStatus' ))
# color by metaTissue
#data <- plotPCA(vst, intgroup=c('metaTissue', 'metaStatus'), ntop=500, returnData=TRUE)
data <- plotPCA(vst, intgroup=c('metaTissue', 'metaStatus'), ntop=5000, returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=metaTissue, shape=metaStatus)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance"))# +
dev.copy(pdf, 'plotPCA_top5000.withoutLiver.pdf'); dev.off()

sampleDists = dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(ncol(vst), 'Blues')))(255)
require('gplots')
colors = greenred(99)
pheatmap(sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors)
dev.copy(pdf, 'sampleClustering.withoutLiver.pdf'); dev.off()
#coord_fixed()
## patient 8 is an outlier, P->M in an oposite direction comparing with the other six patients
## > exclude it and redo the analysis
#PCA with all genes and also plot the third axis
plotPcaFromVst = function(object, intgroup = "condition", ntop = 500,
        returnData = FALSE)
{
	rv <- rowVars(assay(object))
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
		length(rv)))]
	pca <- prcomp(t(assay(object)[select, ]))
	percentVar <- pca$sdev^2/sum(pca$sdev^2)
	if (!all(intgroup %in% names(colData(object)))) {
		stop("the argument 'intgroup' should specify columns of colData(dds)")
	}
	intgroup.df <- as.data.frame(colData(object)[, intgroup,
		drop = FALSE])
	group <- if (length(intgroup) > 1) {
		factor(apply(intgroup.df, 1, paste, collapse = " : "))
	}
	else {
		colData(object)[[intgroup]]
	}
	d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
		intgroup.df, name = colnames(object))
	if (returnData) {
		attr(d, "percentVar") <- percentVar[1:2]
		return(d)
	}
	ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
		geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
		100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *
		100), "% variance")) + coord_fixed()
}

# cal rpkm
#sampleM = colSums(countNz) / 1000000
#geneK = geneInfo[rownames(expr), 'medianTxLength'] / 1000

## Deseq2 DE #using expr
#dds = DESeqDataSetFromMatrix(countData=countNz, colData=sampleInfo, design=~treatment)
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
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~ patientIndex + metaStatus)
dds = DESeq(dds)
y = as.data.frame(expr)
x = results(dds)
x = x[order(x$stat),]
xx = cbind(x, y[rownames(x),])
write.csv(xx, 'DE_paired.M_vs_P.csv')
sele = which(xx$padj<.05)
res = xx[sele,]
write.csv(res, 'DE_paired.M_vs_P.significant.csv')

DESeq2::plotMA(x, alpha=0.05, ylim=c(-4,3))
dev.copy(pdf, 'DE_paired.M_vs_P.plotMA.pdf'); dev.off()
plotGene = function(gene, patients, col=1:12) {
	row = expr[gene,]	
	mat = matrix(row, nrow=2)
	colnames(mat) = sampleInfo$patientIndex[seq(1, nrow(sampleInfo), by=2)]
	rownames(mat) = c('P', 'M')
	matplot(mat, type='l', lty=1, col=col,xaxt='n', ylab='vst', main=gene)
}	



par(mfrow=c(2,5))
dnGenes = head(rownames(res), 5)
upGenes = tail(rownames(res), 5)
patients = as.character(sampleInfo$patientIndex[seq(1, nrow(sampleInfo), by=2)])
col = 1:12
for (i in c(dnGenes, upGenes)) {
	plotGene(i, patients, col)
}
source('~/code/r_util/legend.R')
outerLegend('topright', legend=patients, fill=col)
dev.copy(pdf, 'DE_paired.M_vs_P.top5.plots.pdf')
dev.off()

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


## END
x = results(dds)
x = x[order(x$stat),]
xx = cbind(x, expr[rownames(x),])
write.csv(xx, 'DE_paired.M_vs_P.csv')
sele = which(xx$padj<.05)
res = xx[sele,]
write.csv(res, 'DE_paired.M_vs_P.significant.csv')

plotMA(x, alpha=0.05, ylim=c(-4,3))
dev.copy(pdf, 'DE_paired.M_vs_P.plotMA.pdf'); dev.off()
plotGene = function(gene, patients, col=1:12) {
	row = expr[gene,]	
	mat = matrix(row, nrow=2)
	colnames(mat) = sampleInfo$patientIndex[seq(1, nrow(sampleInfo), by=2)]
	rownames(mat) = c('P', 'M')
	matplot(mat, type='l', lty=1, col=col,xaxt='n', ylab='vst', main=gene)
}	



par(mfrow=c(2,5))
dnGenes = head(rownames(res), 5)
upGenes = tail(rownames(res), 5)
patients = as.character(sampleInfo$patientIndex[seq(1, nrow(sampleInfo), by=2)])
col = 1:12
for (i in c(dnGenes, upGenes)) {
	plotGene(i, patients, col)
}
source('~/code/r_util/legend.R')
outerLegend('topright', legend=patients, fill=col)
dev.copy(pdf, 'DE_paired.M_vs_P.top5.plots.pdf')
dev.off()

## batch effects?

#args = commandArgs(trailingOnly=T)
setwd('/vol1/Working/YanLab/Meiling/Project_Mz424/bowtie2hg19/deseq2')

# annotate the samples manually with meiling
sampleInfo = read.csv('../../sampleInfo.csv', row.names=1)
sele = sampleInfo$pick=='Y'
sampleInfo=sampleInfo[sele,]
sampleInfo$metaStatus = relevel(sampleInfo$metaStatus, 'P')
sampleInfo$patientIndex = as.factor(sampleInfo$patientIndex)
summary(sampleInfo)

tmp = read.delim('../Sample.featureCounts.count_', row.names=1)
geneCount = tmp[,rownames(sampleInfo)]
nz = rowSums(geneCount)>0
sum(nz) #20K
geneCountNz = geneCount[nz,]
head(geneCountNz)
write.csv(geneCountNz, 'geneCountNz.csv')

## Deseq2
require(DESeq2)
colData = sampleInfo
countData = geneCountNz
all(colnames(countData) == rownames(colData))
rownames(colData)=paste0(colData$patientIndex, '.', colData$metaStatus)
colnames(countData) = rownames(colData)
# setup obj from count data
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~metaStatus)
vst = varianceStabilizingTransformation(dds, blind=T)
expr = assay(vst)
write.csv(expr, 'expr.vst.csv')

library("RColorBrewer")
require('ggplot2')

#PCA
#plotPCA(vst, intgroup=c('patientIndex','metaStatus' ))
plotPCA(vst, intgroup=c('patientIndex', 'metaStatus'))
data <- plotPCA(vst, intgroup=c('patientIndex', 'metaStatus'), ntop=500, returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=patientIndex, shape=metaStatus)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance"))# +
dev.copy(pdf, 'plotPCA_top500.pdf'); dev.off()
#coord_fixed()
## patient 8 is an outlier, P->M in an oposite direction comparing with the other six patients
## > exclude it and redo the analysis
#PCA with all genes and also plot the third axis
rowVars = function(x){ apply(x, 1, var) }
ord = order(rowVars(expr), decreasing=T)
topExpr = expr[head(ord, n=500),]
pcaRes = prcomp(topExpr)

par(mfrow=c(2,2))
biplot(pcaRes)
plot(pcaRes$rotation[,1:2])
text(pcaRes$rotation[,1:2], colnames(expr))
plot(pcaRes)
dev.copy(pdf, 'myplot.pdf'); dev.off()

percentVar = summary(pcaRes)$importance[2,] * 100
data = cbind(sampleInfo, as.data.frame(pcaRes$rotation)) 
ggplot(data, aes(PC1, PC2, color=patientIndex, shape=metaStatus)) +
	geom_point(size=3) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))

with(as.data.frame(pcaRes$rotation), plot(PC1, PC2, col=sampleInfo$patientIndex, pch=c(16,17,18,15)[sampleInfo$metaStatus]))
legend('center', 

pca = prcomp(t(topExpr))

par(mfrow=c(2,2))
biplot(pca)
plot(pca$rotation[,1:2])
text(pca$rotation[,1:2], colnames(expr))
plot(pca)
dev.copy(pdf, 'myplotPca.pdf'); dev.off()

percentVar = summary(pca)$importance[2,] * 100
data = cbind(sampleInfo, as.data.frame(pca$x)) 
ggplot(data, aes(PC1, PC3, color=patientIndex, shape=metaStatus)) +
	geom_point(size=3) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC3: ",percentVar[3],"% variance"))

with(as.data.frame(pca$rotation), plot(PC1, PC2, col=sampleInfo$patientIndex, pch=c(16,17,18,15)[sampleInfo$metaStatus]))
legend('center', 


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

## Deseq2 DE
#dds = DESeqDataSetFromMatrix(countData=countNz, colData=sampleInfo, design=~treatment)
# unpaired
dds = DESeq(dds)
group = colData$metaStatus
y = as.data.frame(expr)
x = results(dds)
x = x[order(x$stat),]
xx = cbind(x, y[rownames(x), c(which(group=='P'), which(group=='M'))])
write.csv(xx, 'DE_unpaired.M_vs_P.csv')
sele = xx$padj<.05
sele[is.na(sele)] = F
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
sele = xx$padj<.05
sele[is.na(sele)] = F
res = xx[sele,]
write.csv(res, 'DE_paired.M_vs_P.significant.csv')

plotMA(x, alpha=0.05, ylim=c(-4,3))
dev.copy(pdf, 'DE_paired.M_vs_P.plotMA.pdf'); dev.off()
plotGene = function(gene) {
	row = expr[gene,]	
	mat = matrix(row, nrow=2)
	colnames(mat) = c(1,2,5,6,7,9)
	rownames(mat) = c('P', 'M')
	matplot(mat, type='l', lty=1, col=1:6)
	legend('topright', legend=colnames(mat), fill=1:6)
}	

par(mfrow=c(2,5))
dnGenes = head(rownames(res), 5)
upGenes = tail(rownames(res), 5)
for (i in c(dnGenes, upGenes)) {
	plotGene(i)
}
#contrasts = list(
#        hmmr2D=c('B', 'A'),
#        hmmr3D=c('D', 'C'),
#        hmmrCo=c('G', 'F'))
#groupCol = 'treatment'; group = sampleInfo[,groupCol]
#de = list()
#for (i in names(contrasts)) {
#    cc = contrasts[[i]]
#    x = results(dds, contrast=c(groupCol, cc))
#    x = x[order(x$stat),]
#    xx = cbind(x, y[rownames(x), c(which(group==cc[2]), which(group==cc[1]))])
#    de[[i]] = xx
#}




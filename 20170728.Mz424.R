#args = commandArgs(trailingOnly=T)
expr = read.csv('DE_paired.M_vs_P. expr.csv', row.names=1)
genesPick = read.csv('DE_paired.M_vs_P. pickGenes.csv', as.is=T)[,1]
exprPicked = expr[genesPick,]
write.csv(exprPicked, 'exprPicked.csv')

source ('~/code/working/heatmap.R')
expr.z = t(scale(t(exprPicked)))
heatmap2(as.matrix(expr.z))
heatmap2(as.matrix(expr.z), Colv=NULL, Rowv=NULL)
heatmap2(as.matrix(expr.z), Colv=NULL)

# reorder the expr.z with meta tissue type
sampleInfo = read.csv('/vol/Working/YanLab/Meiling/Project_Mz424/sampleInfo.csv', as.is=T)
row.names = paste0('X', sampleInfo$patientIndex, '.', sampleInfo$metaStatus)
rownames(sampleInfo) = row.names
sampleInfo = sampleInfo[colnames(expr),]

tmp = sampleInfo$tissue
tmp[seq(1,length(tmp), 2)] = tmp[seq(2,length(tmp), 2)]
sampleInfo$metaTissue = tmp

ord = order(sampleInfo$metaTissue, sampleInfo$sampleIndex, sampleInfo$metaStatus)
sampleInfo[ord, c(2:5, ncol(sampleInfo))]
sampleInfo = sampleInfo[ord, ]
expr.z = expr.z[,ord]

# heatmap
heatmap2(as.matrix(expr.z), Colv=NULL)

#rename the samle names in expr.z
col.names = paste0(colnames(expr.z), '.', substr(sampleInfo$metaTissue, 1, 2))
colnames(expr.z) = col.names
heatmap2(as.matrix(expr.z), Colv=NULL)

dev.copy(pdf, 'pickedGenes.heatmap.pdf'); dev.off()



# adjust expr for patient effects (debatch)
require(sva)
expr = expr[,ord]
sampleInfo$patientIndex = as.factors(patientIndex)
expr.combat = ComBat(expr, batch=sampleInfo$patientIndex, mod=model.matrix(~sampleInfo$metaStatus))

# repick genes
genesPick = read.csv('pickGenes_new.csv', as.is=T)[,1]

exprPicked = expr.combat[genesPick,rownames(sampleInfo)]
expr.z = t(scale(t(exprPicked)))
col.names = paste0(colnames(expr.z), '.', substr(sampleInfo$metaTissue, 1, 2))
colnames(expr.z) = col.names
heatmap2(as.matrix(expr.z))
dev.copy(pdf, 'debatchedExpr.pickedGenes_new.heatmap.pdf'); dev.off()

heatmap2(as.matrix(expr.z), Colv=NULL)
dev.copy(pdf, 'debatchedExpr.pickedGenes_new.paired.heatmap.pdf'); dev.off()

# repick genes without debatch
exprPicked = expr[genesPick,rownames(sampleInfo)]
expr.z = t(scale(t(exprPicked)))
col.names = paste0(colnames(expr.z), '.', substr(sampleInfo$metaTissue, 1, 2))
colnames(expr.z) = col.names

heatmap2(as.matrix(expr.z))
dev.copy(pdf, 'originalExpr.pickedGenes_new.heatmap.pdf'); dev.off()

heatmap2(as.matrix(expr.z), Colv=NULL)
dev.copy(pdf, 'originalExpr.pickedGenes_new.paired.heatmap.pdf'); dev.off()

# reorder the expr.z with meta tissue type
sampleInfo = read.csv('/vol/Working/YanLab/Meiling/Project_Mz424/sampleInfo.csv', as.is=T)
row.names = paste0('X', sampleInfo$patientIndex, '.', sampleInfo$metaStatus)
rownames(sampleInfo) = row.names
sampleInfo = sampleInfo[colnames(expr),]

tmp = sampleInfo$tissue
tmp[seq(1,length(tmp), 2)] = tmp[seq(2,length(tmp), 2)]
sampleInfo$metaTissue = tmp

ord = order(sampleInfo$metaTissue, sampleInfo$sampleIndex, sampleInfo$metaStatus)
sampleInfo[ord, c(2:5, ncol(sampleInfo))]
sampleInfo = sampleInfo[ord, ]
expr.z = expr.z[,ord]

# heatmap
heatmap2(as.matrix(expr.z), Colv=NULL)

#rename the samle names in expr.z
col.names = paste0(colnames(expr.z), '.', substr(sampleInfo$metaTissue, 1, 2))
colnames(expr.z) = col.names
heatmap2(as.matrix(expr.z), Colv=NULL)

dev.copy(pdf, 'pickedGenes.heatmap.pdf'); dev.off()

heatmap2(as.matrix(expr.z))

readCsv = function(x) read.csv(x, row.names=1, check.names=F, stringsAsFactors=F)
source ('~/code/working/heatmap.R')
setwd('/vol/Working/YanLab/Meiling/Project_Mz424/20180903.Bowtie2Deseq2/20181101.draft_heatmap')
# read Expr
expr = read.csv('../expr.vst.csv', row.names=1, check.name=F)

# read sampleInfo and align with expr
sampleInfo = read.csv('../sampleInfo.csv', comment='#')
sampleInfo = sampleInfo[sampleInfo$pick=='Y',]
#tmp = paste(sampleInfo$patientIndex, sampleInfo$metaStatus, sep='.')
#tmp == colnames(expr)
rownames(sampleInfo) = colnames(expr)

# adjust expr for patient effects (debatch)
require(sva)
sampleInfo$patientIndex = factor(sampleInfo$patientIndex)
expr.combat = ComBat(expr, batch=sampleInfo$patientIndex, mod=model.matrix(~sampleInfo$metaStatus))
write.csv(expr.combat, 'expr.combat.csv')
write.csv(sampleInfo, 'sampleInfo.csv')

# pick epigenetic genes with fc>1.2
epigene = read.csv('/vol/Working/YanLab/Meiling/epigenetic_genes.csv', row.names=1)
siggene = read.csv('../DE_paired.M_vs_P.significant.csv', row.names=1)
common = intersect(rownames(epigene), rownames(siggene))
tmp = siggene[common,]
sele = abs(tmp$log2FoldChange)>log2(1.2)
genesPick = rownames(tmp)[sele]

# plot the combated expr
#colnames(expr.combat) == rownames(sampleInfo)
exprPicked = expr.combat[genesPick,]
expr.z = t(scale(t(exprPicked)))
#col.names = paste0(colnames(expr.z), '.', substr(sampleInfo$metaTissue, 1, 2))
#colnames(expr.z) = col.names
heatmap2(as.matrix(expr.z))
dev.copy(pdf, 'debatchedExpr.epigenes.padj_.05_FC_1.2_heatmap.pdf'); dev.off()

# reorder the samples by patient (the original order)
heatmap2(as.matrix(expr.z), Colv=F)
dev.copy(pdf, 'debatchedExpr.epigenes.padj_.05_FC_1.2_heatmap.orig_order.pdf'); dev.off()



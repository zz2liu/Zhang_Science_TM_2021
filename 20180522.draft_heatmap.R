readCsv = function(x) read.csv(x, row.names=1, check.names=F, stringsAsFactors=F)
source ('~/code/working/heatmap.R')
setwd('/vol/Working/YanLab/Meiling/Project_Mz424/20170728.heatmap')

# read Expr
expr = read.csv('DE_paired.M_vs_P. expr.csv', row.names=1)

# read and manipulate sampleInfo 
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

# adjust expr for patient effects (debatch)
require(sva)
expr = expr[,ord]
sampleInfo$patientIndex = as.factors(patientIndex)
expr.combat = ComBat(expr, batch=sampleInfo$patientIndex, mod=model.matrix(~sampleInfo$metaStatus))
write.csv(expr.combat, 'expr.combat.csv')
write.csv(sampleInfo, 'sampleInfo.csv')

# repick genes
#TET3 removed because of small FC
genesPick = read.csv('pickGenes_new.csv', as.is=T)[,1]
exprPicked = expr.combat[genesPick,rownames(sampleInfo)]
expr.z = t(scale(t(exprPicked)))
col.names = paste0(colnames(expr.z), '.', substr(sampleInfo$metaTissue, 1, 2))
colnames(expr.z) = col.names
heatmap2(as.matrix(expr.z))
dev.copy(pdf, 'debatchedExpr.pickedGenes_new.heatmap.pdf'); dev.off()

# 20180522 reorder the genes with Stat
de = readCsv('../Bowtie2Deseq2/DE_paired.M_vs_P.csv')
dePicked= de[genesPick,]
ord = order(dePicked$stat)
#heatmap2(as.matrix(expr.z)) #, Rowv=F)
heatmap2(as.matrix(expr.z[ord,]), Rowv=F)
dev.copy(pdf, 'debatchedExpr.pickedGenes_new.heatmap.orederByStat.pdf'); dev.off()




#heatmap2(as.matrix(expr.z), Colv=NULL)
#dev.copy(pdf, 'debatchedExpr.pickedGenes_new.paired.heatmap.pdf'); dev.off()

## reorder the expr.z with meta tissue type
#sampleInfo = read.csv('/vol/Working/YanLab/Meiling/Project_Mz424/sampleInfo.csv', as.is=T)
#row.names = paste0('X', sampleInfo$patientIndex, '.', sampleInfo$metaStatus)
#rownames(sampleInfo) = row.names
#sampleInfo = sampleInfo[colnames(expr),]
#
#tmp = sampleInfo$tissue
#tmp[seq(1,length(tmp), 2)] = tmp[seq(2,length(tmp), 2)]
#sampleInfo$metaTissue = tmp
#
#ord = order(sampleInfo$metaTissue, sampleInfo$sampleIndex, sampleInfo$metaStatus)
#sampleInfo[ord, c(2:5, ncol(sampleInfo))]
#sampleInfo = sampleInfo[ord, ]
#expr.z = expr.z[,ord]
#
## heatmap
#heatmap2(as.matrix(expr.z), Colv=NULL)
#
##rename the samle names in expr.z
#col.names = paste0(colnames(expr.z), '.', substr(sampleInfo$metaTissue, 1, 2))
#colnames(expr.z) = col.names
#heatmap2(as.matrix(expr.z), Colv=NULL)
#
#dev.copy(pdf, 'pickedGenes.heatmap.pdf'); dev.off()
#
#heatmap2(as.matrix(expr.z))

# detect the genes with a bigger/smaller change in Meta-matchedNormal vs Prim-matchedNormal
setwd("/tank/home/zongzhi/sshfs/ruddleProject/tmpMeiling/hg38.deseq2.paired_NM")
vst = read.csv('gene.vst_blind.csv', row.names=1, check.names=F)
vst = as.matrix(vst)
sampleInfo = read.csv('sampleInfo.filtered.csv', row.names=1, check.names=F)
#play with data.table
require(data.table)
require(dplyr)
dt = data.frame(sampleInfo, keep.rownames='rownames')

# get a matrix of [gene x tissue] difference (tumor-normal)
normCol = seq(1, nrow(sampleInfo), by=2)
tumorCol = normCol+1
tissue = sampleInfo[normCol,'batch']
delta = vst[,tumorCol] - vst[, normCol]
colnames(delta) = tissue

# perform a t test for each gene: breast vs others
require(genefilter)
fac = factor(c('primary', rep('metasite', length(tissue)-1)), levels=c('metasite', 'primary'))
res = rowttests(elta), fac)
res$padj = p.adjust(res$p.value)
result = cbind(res, as.matrix(delta))
write.csv(result, 'differtial_delta.metasite-primary.ttests.csv')

# output significant genes
sele = with(result, padj < .05)
tmp = result[sele,]
write.csv(sig, 'differtial_delta.metasite-primary.ttests.sig.csv')
    # not quite right, all the sigones are almost all zeros.

#> use log2(M_i/N - P/N) to test instead
exprRaw = vst^2
tumorFC = exprRaw[,tumorCol] / exprRaw[,normCol]
colnames(tumorFC) = tissue
primCol=1
deltaFC = tumorFC[,-primCol] - tumorFC[,primCol]
logDeltaFC = log2(deltaFC)
logDelta
    # not works

#> use linear model with a interaction item tumor:meta, and tissue as random variable in additive model
info = data.frame(
    tumor = factor(rep(c(0,1), times=nrow(sampleInfo)/2)),
    meta = factor(c(0, 0, rep(1, nrow(sampleInfo)-2))),
    tissue = relevel(sampleInfo$batch, 'breast'))

# test a gene
i=1
y = vst[i,]
fit = lm(y ~ tumor + meta + tumor:meta + tissue, data=info)
tmp = summary(fit)$coefficients
res = tmp[nrow(tmp),]
gene = rownames(vst)[i]
pdf(paste0(gene, '.pdf'))
g = ggplot(info, aes(tumor, y, color=meta)) + geom_line(aes(group=tissue)) + 
    ggtitle(paste(gene, ' pvalue=', res[length(res)]))
plot(g)
par(mfrow=c(2,2))
plot(fit)
dev.off()

# gather results from all genes; using info
fitOneRow = function(y) {
    fit = lm(y ~ tumor + meta + tumor:meta + tissue, data=info)
    tmp = summary(fit)$coefficients
    res = tmp[nrow(tmp),]
}
#tmp = apply(head(vst), 1, fitOneRow)
res = t(apply(vst, 1, fitOneRow))
padj = p.adjust(res[,4])
result = cbind(res, padj)
write.csv(result, 'tumor_meta_interaction.csv')

result = as.data.frame(result)
sele = with(result, padj < .05)
tmp = result[sele,]
ord = order(tmp[,3])
sig = tmp[ord,]
sigif = cbind(sig, vst[rownames(sig),])
write.csv(sig, 'tumor_meta_interaction.sig.csv')

# plot one gene: using vst, info
plotOneGene = function(gene) {
    y = vst[gene,]
    write.csv(cbind(info, y), paste0(gene, '.input.csv'))

    fit = lm(y ~ tumor + meta + tumor:meta + tissue, data=info)
    tmp = summary(fit)$coefficients
    res = tmp[nrow(tmp),]
    write.csv(tmp, paste0(gene, '.summary.csv'))

    pdf(paste0(gene, '.pdf'))
    g = ggplot(info, aes(tumor, y, color=meta)) + geom_line(aes(group=tissue)) + 
        ggtitle(paste(gene, ' pvalue=', res[length(res)]))
    plot(g)
    par(mfrow=c(1,1))
    textplot(tmp)
    par(mfrow=c(2,2))
    plot(fit)
    dev.off()
}



#  For each meta tissue, perform a z test for each gene: one meta tissue vs other meta tissues


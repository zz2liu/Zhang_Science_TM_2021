require (ggplot2)
require (reshape2)
gplotGenes = function(genes, col=1:12,
        title='Up in Metastatic tissue', height=5, width=9) {
    mat = expr[genes,]
    dat = melt(mat)
    colnames(dat) = c('gene', 'patient_status', 'expr')
    dat$patient = sub('.[PM]','', dat$patient_status)
    dat$metaStatus = sub('[0-9]+.','', dat$patient_status)
    dat$metaStatus = factor(dat$metaStatus, levels=c('P', 'M'))
    g = ggplot(dat, aes(x=metaStatus, y=expr, group=patient, color=patient))
    g + geom_point() + geom_line() + facet_grid(~gene) + labs(title=title)
    ggsave(paste0(title, '.pdf'), height=height)
}

expr = read.csv('expr.vst.csv', row.names=1, check.names=F)
head(res)
dnGenes = head(rownames(res), 10)
gplotGenes(c('CECR2'), title='tmp.Down_in_Meta', height=4)
gplotGenes(dnGenes, title='Down_in_Meta', height=4)
upGenes = tail(rownames(res), 10)
gplotGenes(upGenes, title='Up_in_Meta', height=4)




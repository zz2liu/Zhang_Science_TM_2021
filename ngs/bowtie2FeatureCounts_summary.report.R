#usage: PROG bowtie2FeatureCounts.summary.report <bowtie2FeatureCounts.summary.csv 
# outPrefix: used to output .csv and .pdf
# stdin: the csv file generated from bowtie2Summary.py
args = commandArgs(trailingOnly=T)
if (length(args) != 1) {
    print ('USAGE: PROG <outPrefix>')
}
require(reshape2)
inFile = file('stdin')
outPrefix = args[1]
outCsv = paste0(outPrefix, '.csv')
outPdf = paste0(outPrefix, '.pdf')

tmp = read.csv(inFile, row.names=1, check.names=F) #'bowtie2Summary.csv')
tmp$poorMapq = tmp[,1] - rowSums(tmp[,-1]) # first column is total
res = tmp[,-1]
ord = c('oneGene', 'manyGene', 'noGene', 'poorMapq', 'unaligned')
res = res[,ord]
write.csv(res, outCsv)

source('~/code/r_util/legend.R')
pdf(outPdf, width=12)
col = rev(1:ncol(res))
par(oma=c(2,2,2,8))
tmp=barplot(t(res), col=col, legend=F, las=2) #, ylab='# of reads')
mtext('# of reads', side=2, line=4)
#plot(0,0, oma=0, mar=rep(0,4),fig=c(0,1,0,1), new=T)
browser()
outerLegend('topright', rev(colnames(res)), fill=rev(col))

resFrac = res/rowSums(res)
barplot(t(resFrac), col=col, legend=F, las=2, ylab='Fraction')
outerLegend('topright', rev(colnames(res)), fill=rev(col))
dev.off()


## alternative way using ggplot2
require(ggplot2)
res$sampleId = rownames(res)
data = data.table::melt(res, 'sampleId')
g = ggplot(data, aes(x=sampleId, y=value, fill=variable))
g + geom_bar(stat='identity')
ggsave(paste0(outPrefix, '.svg'))

require(sqldf)
require(scales)
tmp = sqldf('select sampleId, sum(value) total from data group by sampleId') 
dataFrac = sqldf('select sampleId, variable, 1.*value/total fraction 
    from data join tmp using (sampleId)')
g = ggplot(dataFrac, aes(x=sampleId, y=fraction, fill=variable))
g + geom_bar(stat='identity') + scale_y_continuous(labels=percent)
ggsave(paste0(outPrefix, '.percent.svg'))

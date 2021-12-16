args = commandArgs(trailingOnly=T)

##R
counts = read.delim('Sample.featureCounts.count_', row.names=1)
#counts = as.matrix(df[,6:9])
#colnames(counts) = substr(colnames(counts), 8, 12)
rpm = t(t(counts)/colSums(counts))*1000000
goodRpm = rpm[rowSums(rpm)>1,]
pdf('log10 goodRpm.pdf'); hist(log10(goodRpm)); dev.off()
expr = log2(goodRpm+.0001)
pdf('log2 goodRpm.pdf'); hist(expr); dev.off()
source('~/code/working/heatmap.R')
runHeatmap(expr, 'log2rpm.heatmap.pdf')

## numGenes covered by numReads or RPM

## insert size of the pairs

## count summary
countSummary = read.delim('Sample.featureCounts.count.summary_', row.names=1)
x = as.matrix(countSummary)
x = x[rowSums(x)>0,]
countSummary = x
countSummaryPercentage= t(t(x)/colSums(x))
write.csv(countSummaryPercentage, 'Sample.featureCounts.countSummaryPercentage.csv')
write.csv(countSummary, 'Sample.featureCounts.countSummary.csv')
##
pdf('Sample.featureCounts.countSummaryPercentage.pdf'); barplot(countSummaryPercentage, las=2, legend=T, args.legend=list(x='bottomleft')); dev.off()
pdf('Sample.featureCounts.countSummary.pdf'); barplot(countSummary, las=2, legend=T, args.legend=list(x='bottomleft')); dev.off()



### bowtie2 Summary
df = read.csv('bowtie2Summary_.csv')
require(reshape2)
x = acast(df,alignCatagory~fileName)
write.csv(x, 'bowtie2SummaryMatrix.csv')
bowtie2Summary = x
ord = c(2,1,3)
bowtie2Summary=bowtie2Summary[ord,]

pdf('bowtie2Summary.pdf'); 
barplot(bowtie2Summary, las=2, legend=T, args.legend=list(x='bottomleft'),col=2:7);
barplot(t(t(bowtie2Summary)/colSums(bowtie2Summary)), las=2, legend=T, args.legend=list(x='bottomleft'),col=2:7);
dev.off()



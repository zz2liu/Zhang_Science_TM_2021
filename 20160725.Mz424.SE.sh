ppn=8
mem=32G
qsub -v ppn=$ppn,mem=$mem -I -q interactive -l nodes=1:ppn=$ppn -l mem=$mem

genomeDir=~/local/data/hg19
export genomePrefix=$genomeDir/Sequence/Bowtie2Index/hg19
export geneGtf=$genomeDir/Annotation/Genes/genes.gtf 

workDir=~/scratch60/tmpMeiling/Project_Mz424
cd $workDir
mkdir bowtie2hg19
cd bowtie2hg19

for i in $(ls -d ../Sample*); do
    echo $i; date
    sampleDir=$i bash -x ~/code/seqyale/bowtie2.sh
    #break
done    
# count with featureCounts# much faster
~/local/install/subread-1.5.0-p1-Linux-x86_64/bin/featureCounts -T $ppn -s 2 -t exon -g gene_id -a $geneGtf -o Sample.featureCounts.count *.bam 
    #s2 maps a bit less that s0 

#summary bowtie
python ~/code/ngs/bowtie2Summary.py *.bowtie2.log > bowtie2Summary.csv
cat bowtie2Summary.csv | sed 's/Sample_\|MZ424_\|.bowtie2.log//g' | sed '/totalReads/d' > bowtie2Summary_.csv
cat Sample.featureCounts.count | sed '1d' | sed '1s/Sample_\|MZ424_\|.q.10.bam//g' | cut -f1,7- > Sample.featureCounts.count_
cat Sample.featureCounts.count.summary | sed '1s/Sample_\|MZ424_\|.q.10.bam//g' > Sample.featureCounts.count.summary_
#summary and plot count
Rscript $scriptDir/20160725.Mz424.R

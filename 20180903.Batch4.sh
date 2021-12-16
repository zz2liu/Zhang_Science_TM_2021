cd ~/scratch60/tmpMeiling/Batch4
mkdir bowtie2hg19
cd bowtie2hg19

genomeDir=~/local/data/hg19
export genomePrefix=$genomeDir/Sequence/Bowtie2Index/hg19
export gtfFile=$genomeDir/Annotation/Genes/genes.gtf 
export chromSizeFile=$genomeDir/chrom.sizes

workDir=$PWD
ppn=8
mem=32G

##test interactively
srunI
i=../input/Sample_M7_28_2
outDir=$(echo $i | sed 's/.*Sample_//')
mkdir $outDir
echo $i; date
export sampleDir=$(realpath $i)
export outDir=$(realpath $outDir)
##test interactively
bash -x \
    ~/code/ngs/bowtie2localSe.rpm.featureCount.pipeline.sh

##test batch
sbatch -c8 -N1 --mem-per-cpu=4000 \
    ~/code/ngs/bowtie2localSe.rpm.featureCount.pipeline.sh

## run batch
for i in $(ls -d ../input/Sample*); do
    outDir=$(echo $i | sed 's/.*Sample_//')
    mkdir $outDir
    echo $i; date
    export sampleDir=$(realpath $i)
    export outDir=$(realpath $outDir)
    sbatch -c8 -N1 --mem-per-cpu=4000 \
        ~/code/ngs/bowtie2localSe.rpm.featureCount.pipeline.sh
done
# squeue -u zl99 #monitor

##wait until all jobs finished.
#summary bowtie and featureCounts
find -name "bowtie2.log" | xargs python ~/code/ngs/bowtie2Summary.py > tmp.csv
Rscript ~/code/ngs/bowtie2Summary.report.R bowtie2Summary < tmp.csv

bash ~/code/ngs/bowtie2FeatureCounts_summary.sh */ | sed 's%/%%' > tmp.csv
Rscript ~/code/ngs/bowtie2FeatureCounts_summary.report.R bowtie2FeatureCounts_summary.report <tmp.csv

#join the count files
files=$(echo */ | sed 's%/%/featureCounts.txt%g')
python3 ~/code/ngs/join_count_files.py $files > tmp.csv
samples=$(echo */ | sed 's%/%%g;s/ /,/g')
cat tmp.csv | sed "1d;2s/,.*/,$samples/" > featureCounts.csv

## DE analysis

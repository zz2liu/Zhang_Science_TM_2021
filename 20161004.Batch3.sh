cd ~/scratch60/tmpMeiling/Batch3
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
qsub -I -q interactive \
    -V -d . -l nodes=1:ppn=$ppn -l mem=$mem

    i=../rawData/Sample_MZ424_M4_007
    outDir=$(echo $i | sed 's/.*Sample_MZ424_//')
    mkdir $outDir
    echo $i; date
    export sampleDir=$(realpath $i)
    export outDir=$(realpath $outDir)
    bash -x ~/code/ngs/bowtie2localSe.rpm.featureCount.pipeline.sh

## run batch
for i in $(ls -d ../rawData/Sample*); do
    id=$(echo $i | sed 's/.*Sample_MZ424_//')
    mkdir $id
    echo $id; date
    export sampleDir=$(realpath $i)
    export outDir=$(realpath $id)
    qsub -N $id \
        -V -d . -l nodes=1:ppn=$ppn -l mem=$mem \
        ~/code/ngs/bowtie2localSe.rpm.featureCount.pipeline.sh
done    

#summary bowtie and featureCounts
find -name "bowtie2.log" | xargs python ~/code/ngs/bowtie2Summary.py | sed 's%./%%;s%/bowtie2.log%%;/totalReads/d' > tmp.csv
Rscript ~/code/ngs/bowtie2Summary.report.R bowtie2Summary < tmp.csv

bash ~/code/ngs/bowtie2FeatureCounts_summary.sh */ | sed 's%/%%' > tmp.csv
Rscript ~/code/ngs/bowtie2FeatureCounts_summary.report.R bowtie2FeatureCounts_summary.report <tmp.csv

#join the count files
files=$(echo */ | sed 's%/%/featureCounts.txt%g')
python3 ~/code/ngs/join_count_files.py $files > tmp.csv
samples=$(echo */ | sed 's%/%%g;s/ /,/g')
cat tmp.csv | sed "1d;2s/,.*/,$samples/" > featureCounts.csv

## DE analysis

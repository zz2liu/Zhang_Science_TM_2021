#!/bin/bash
### build the SE pipeline for each sample
#create test case with 10k reads
#qsub -I -q interactive -d . -lnodes=1:ppn=8 -lmem=32g
## require
#genomePrefix=~/local/data/hg38/genome
#gtfFile=~/local/data/hg38/gencode.v24.annotation.level12.gtf 
#chromSizeFile=~/local/data/hg38.chrom.sizes
#sampleDir, outDir
memEach=4 #to calc from $mem
nThread=8 #=$ppn
minMapq=10 #=$minMapq

PATH=~/local/install/subread-1.5.0-p1-Linux-x86_64/bin:$PATH
#module load Apps/ucsc/2016-05-27

## run one sample
function main { #using sampleDir,outDir,minMapq, memEach, nThread
    echo $sampleDir; date
    tmpDir=$(mktemp -d)
    fastqR1=$(ls  $sampleDir/*_R1_???.fastq.gz | tr '\n' , | rev | cut -c 2- |rev)
    # map to hgmm genome
    cd $tmpDir
    bowtie2 -p $nThread --local --non-deterministic \
            -x $genomePrefix -U $fastqR1 \
            2> bowtie2.log \
            | samtools view -Sbh -q $minMapq - > tmp.bam
    #sort
    samtools sort -m ${memEach}G -@ $nThread tmp.bam > sorted.bam \
        && rm tmp.bam
    samtools index sorted.bam

    #rpm.bigwig
    bamFile=sorted.bam
    totalReads=$(samtools view -c $bamFile)
    echo $totalReads > $bamFile.totalReads
    scale=$(bc <<< "scale=6; 1000000/$totalReads")
    bedtools genomecov -ibam $bamFile -scale $scale -g $chromSizeFile -bg > tmp.bedgraph
    bedSort tmp.bedgraph tmp.bedgraph
    bedGraphToBigWig tmp.bedgraph $chromSizeFile $bamFile.rpm.bw

    #gene count #s0 unstranded
    featureCounts -T $nThread -a $gtfFile -o featureCounts.txt -t exon -g gene_id -s 0 sorted.bam 2> featureCounts.log

    #clean up
    #rm  *.sam *.bedgraph #*.bam *.bam.bai
    mkdir $outDir
    rsync -azu --exclude='tmp.*' --remove-source-files $tmpDir/* $outDir
    rm -r $tmpDir
}

main

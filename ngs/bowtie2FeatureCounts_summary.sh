# Example PROG */ 
# under each sampleDir should have a bowtie2.log and a featureCounts.txt.summary

tmpDir=$(mktemp -d)

#for i in ${sampleDirs[@]}; do 
for i in $@; do 
    echo $i >> $tmpDir/sampleId
    cat $i/bowtie2.log | head -n1 | cut -f1 -d' '>>$tmpDir/total
    cat $i/bowtie2.log | grep ' aligned 0 times' | sed 's/^ *//' | cut -f1 -d' ' >> $tmpDir/unaligned
    cat $i/featureCounts.txt.summary | grep '^Assigned' | cut -f2 >> $tmpDir/oneGene
    cat $i/featureCounts.txt.summary | grep '^Unassigned_Ambiguity' | cut -f2 >> $tmpDir/manyGene
    cat $i/featureCounts.txt.summary | grep '^Unassigned_NoFeatures' | cut -f2 >> $tmpDir/noGene
done
paste -d',' $tmpDir/sampleId $tmpDir/total $tmpDir/unaligned $tmpDir/oneGene $tmpDir/manyGene $tmpDir/noGene \
| sed '1i\
sampleId,total,unaligned,oneGene,manyGene,noGene'

rm $tmpDir/*; rmdir $tmpDir

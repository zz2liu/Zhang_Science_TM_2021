# driving_breast_cancer_metastasis
Computing code for the draft: CECR2 Drives Breast Cancer Metastasis by Suppressing Macrophage Inflammatory Responses
## Alignment to genome, count to each annotated gene
packages: bowtie2, FeatureCounts, samtools, bedtools.
scripts: Mz424.SE.sh, Batch3.sh, Batch4.sh.

## Differential Gene Expression
packages: Deseq2
scripts: merged.DE.R

To calculate differential gene expression, the uniquely mapped reads (MAPQ > = 10) were counted to gene annotations using featureCounts (Liao et al., 2014), thus generating a count matrix of gene by sample. The matrix was then used as input to DESeq2 R package (Love et al., 2014). 

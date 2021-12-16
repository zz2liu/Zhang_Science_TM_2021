# Zhang Science TM 2021
Computing code for the draft: CECR2 Drives Breast Cancer Metastasis by Promoting NF-κB Signaling and Macrophage-mediated Immune Suppression

## Alignment to genome, count to each annotated gene
packages: bowtie2, FeatureCounts, samtools, bedtools.

scripts: Mz424.SE.sh, Batch3.sh, Batch4.sh.

## Differential Gene Expression
R packages: Deseq2

scripts: merged.DE.R

To calculate differential gene expression, the uniquely mapped reads (MAPQ > = 10) were counted to gene annotations using featureCounts (Liao et al., 2014), thus generating a count matrix of gene by sample. The matrix was then used as input to DESeq2 R package (Love et al., 2014). 

To mitigate the influence of individual participants, the patient id is integrated into the DESEQ regression model, so that each primary tumor is kind of paired with the counterpart of metastatic tumor. 

## Heatmap
R packages: heatmaply

scripts: meiling_heatmap_again.R

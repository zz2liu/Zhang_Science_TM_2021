# redo the heatmap in figure 1, using data from table S2
# expr = read.csv('/Users/zongzhiliu/OneDrive_hotmail/OneDrive/yanlab/meiling_heatmap.csv', row.names=1)
# source ('~/git/code/working/heatmap.R')
# heatmap2(as.matrix(expr))
# cor_mat = cor(expr)
# heatmap2(cor_mat)
# heatmaply_cor(x=cor_mat,
#               k_col=2,
#               k_row=2,
#               limits=c(min(cor_mat),1))
# vst = read.csv('/Users/zongzhiliu/OneDrive_hotmail/OneDrive/yanlab/regeneratefigures1/DE_paired.M_vs_P.csv',
#                row.names=1)
vst = read.csv('/Users/zongzhiliu/OneDrive_hotmail/OneDrive/yanlab/regeneratefigures1/GSE147995_vst.csv',
                row.names=1)
vst.z = t(scale(t(vst)))
cor_vst = cor(vst.z)
heatmaply_cor(x=cor_vst,
              k_col=2,
              k_row=2,
              limits=c(min(cor_vst),1))

# map name to new set
info = read.csv('/Users/zongzhiliu/OneDrive_hotmail/OneDrive/yanlab/patient_info.csv',
                )
ori = t(data.frame(strsplit(colnames(cor_vst), '.', fixed=T)))
ori = data.frame(ori)
colnames(ori) = c('original', 'status')
require(sqldf)
labels = sqldf("select cast(patient as varchar) || status
               from ori join info using (original)")[[1]]
er = sqldf("select 'ER'||ER
           from ori join info using (original)")[[1]]
x = cor_vst
colnames(x) = labels
rownames(x) = labels
?cor
#plot
heatmaply_cor(x, symm=T,
              #limits=c(min(x), 1),
              col_side_colors=data.frame(ER=er),
              col_side_palette=c('ER+'='blue', 'ER-'='red'),
              row_side_colors=data.frame(ER=er),
              row_side_palette=c('ER+'='blue', 'ER-'='red'),
              show_dendrogram=c(T, F),
              #hclust_method='average',
              #yaxis = list(side = "right")
              )
write.csv(x, 'expr_vst_z_pearson.csv')

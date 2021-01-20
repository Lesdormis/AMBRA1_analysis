# library(devtools)
# install_github("miccec/yaGST")
library(yaGST)
library(doParallel)
library(pheatmap)

load("TCGA_DLBC_fpkm_counts.Robj")
exprData <- t(fpkm_count)

ambra_del <- "TCGA-FF-A7CX"; ambra_mut <- c("TCGA-FA-A86F","TCGA-FF-8041")
ccnd3_mut <- "TCGA-GS-A9TU"; ccnd3_amp <- c("TCGA-GS-A9TQ", "TCGA-GS-A9TY")
cl <- makePSOCKcluster(7)
clusterApply(cl, floor(runif(length(cl), max = 10000000)), set.seed)
registerDoParallel(cl)
ans_eeMWW <- eeMWW(exprData, c(ambra_del,ambra_mut,ccnd3_mut,ccnd3_amp))
#save(ans_eeMWW, file = "TCGA_AMBRA1_CCND3_mut3_WMM.Robj")
int_gene <- c(head(names(ans_eeMWW),n = 150), tail(names(ans_eeMWW),n=150))
palette <- 200
cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(palette)
hp_df <- fpkm_count[int_gene,]
all_sampe <- colnames(fpkm_count); rest_sample <- all_sampe[!all_sampe %in% c(ambra_del,ambra_mut,ccnd3_mut,ccnd3_amp,rb1_del,cdk4_amp,cdk6_amp)]
annotate_col <- data.frame(sample = rep(c("AMBRA1_Dele","AMBRA1_Mut","CCND3_Mut","CCND3_Amp","rest_samples"),c(1,2,1,2,length(rest_sample))))
rownames(annotate_col) <- c(ambra_del,ambra_mut,ccnd3_mut,ccnd3_amp,rest_sample)
color_annotation = list(sample = c(AMBRA1_Dele = "red",AMBRA1_Mut = "purple",CCND3_Mut = "blue",CCND3_Amp = "green",rest_samples = "white"))

hp_df_log <- log2(hp_df + 0.01) #add pesudocount
hp_df_log_z <- scale(t(hp_df_log))
hp_df_log_z <- t(hp_df_log_z)
hp_df_log_z <- as.data.frame(hp_df_log_z)

hp_df_log_z_fix <- FixScale(mydata = hp_df_log_z, min = -2.5, max = 2.5)

ph <- pheatmap(hp_df_log_z_fix, color = cell_colors, cluster_rows = F, cluster_cols = T, 
               border_color = NA, fontsize_row = 1.2, fontsize_col = 1, 
               fontsize = 8, show_rownames =T,show_colnames = F, scale = "none",
               annotation_col = annotate_col,annotation_colors = color_annotation[1],filename = "Fig4e.pdf")
# get the heatmap order
ord <- ph$tree_col$order
hp_df <- hp_df[,ord]

write.csv(hp_df, file = "heatmap_gene_signatures_FPKM_matrix.csv", row.names = T)
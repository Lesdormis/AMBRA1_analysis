options(stringsAsFactors = F)
library(doParallel)
library(pheatmap)
library(openxlsx)
load("TCGA_DLBC_fpkm_counts.Robj")
ambra_del <- "TCGA-FF-A7CX"; ambra_mut <- c("TCGA-FA-A86F","TCGA-FF-8041")
ccnd3_mut <- "TCGA-GS-A9TU"; ccnd3_amp <- c("TCGA-GS-A9TQ", "TCGA-GS-A9TY")
palette <- 200
cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(palette)

all_sampe <- colnames(fpkm_count); rest_sample <- all_sampe[!all_sampe %in% c(ambra_del,ambra_mut,ccnd3_mut,ccnd3_amp)]
annotate_col <- data.frame(sample = rep(c("AMBRA1_Dele","AMBRA1_Mut","CCND3_Mut","CCND3_Amp","rest_samples"),c(1,2,1,2,length(rest_sample))))
rownames(annotate_col) <- c(ambra_del,ambra_mut,ccnd3_mut,ccnd3_amp,rest_sample)
color_annotation = list(sample = c(AMBRA1_Dele = "red",AMBRA1_Mut = "purple",CCND3_Mut = "blue",CCND3_Amp = "green",rest_samples = "white"))

hp_df_log_z_fix = read.xlsx("../Extended Data Figure 10.xlsx",sheet = 3,rowNames = T)
ph <- pheatmap(hp_df_log_z_fix, color = cell_colors, cluster_rows = F, cluster_cols = T, 
               border_color = NA, fontsize_row = 1.2, fontsize_col = 1, 
               fontsize = 8, show_rownames =T,show_colnames = F, scale = "none",
               annotation_col = annotate_col,annotation_colors = color_annotation[1],filename = "Ext_Fig10f.pdf")

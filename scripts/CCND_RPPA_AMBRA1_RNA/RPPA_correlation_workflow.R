library(ggplot2)
library(ggpubr)
library(openxlsx)
options(stringsAsFactors = F)
# LGG
temp_rppa = read.xlsx("../Extended Data Figure 7.xlsx",sheet = 1)
pdf(file = "Ext_Fig7b_LGG.pdf",width = 7,height = 7)
gp <- ggscatter(temp_rppa, x = "AMBRA1.(FPKM-UQ)", y = "cyclin.D1.(RPPA)",
                add = "reg.line",
                conf.int = T,
                add.params = list(color = "blue",
                                  fill = "lightgray"), title = paste("TCGA_LGG")) +
  stat_cor(method = "pearson")
print(gp)
dev.off()
# OV
temp_rppa = read.xlsx("../Extended Data Figure 7.xlsx",sheet = 2)
pdf(file = "Ext_Fig7b_OV.pdf",width = 7,height = 7)
gp <- ggscatter(temp_rppa, x = "AMBRA1.(FPKM-UQ)", y = "cyclin.D1.(RPPA)",
                add = "reg.line",
                conf.int = T,
                add.params = list(color = "blue",
                                  fill = "lightgray"), title = paste("TCGA_OV")) +
  stat_cor(method = "pearson")
print(gp)
dev.off()
# PRAD
temp_rppa = read.xlsx("../Extended Data Figure 7.xlsx",sheet = 3)
pdf(file = "Ext_Fig7b_PRAD.pdf",width = 7,height = 7)
gp <- ggscatter(temp_rppa, x = "AMBRA1.(FPKM-UQ)", y = "cyclin.D1.(RPPA)",
                add = "reg.line",
                conf.int = T,
                add.params = list(color = "blue",
                                  fill = "lightgray"), title = paste("TCGA_PRAD")) +
  stat_cor(method = "pearson")
print(gp)
dev.off()
# TGCT
temp_rppa = read.xlsx("../Extended Data Figure 7.xlsx",sheet = 4)
pdf(file = "Ext_Fig7b_TGCT.pdf",width = 7,height = 7)
gp <- ggscatter(temp_rppa, x = "AMBRA1.(FPKM-UQ)", y = "cyclin.D1.(RPPA)",
                add = "reg.line",
                conf.int = T,
                add.params = list(color = "blue",
                                  fill = "lightgray"), title = paste("TCGA_TGCT")) +
  stat_cor(method = "pearson")
print(gp)
dev.off()

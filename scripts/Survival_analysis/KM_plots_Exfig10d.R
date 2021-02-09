library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(openxlsx)
options(stringsAsFactors = F)
lenz <- read.xlsx("../Extended Data Figure 10.xlsx",sheet = 1)
stank <- read.xlsx("../Extended Data Figure 10.xlsx",sheet = 2)
# scale each expression on each cohort
lenz$AMBRA1 <- scale(lenz$AMBRA1)

stank$AMBRA1 <- scale(stank$AMBRA1)
cn <- c("GEO_ID","Surv_status","surv_years","AMBRA1")

colnames(stank) <- cn; colnames(lenz) <- cn

pdf("Ext_Fig10d.pdf")
  q = 0.1
  # lenz
  qua = quantile(lenz[,"AMBRA1"],c(q,1-q)) 
  df <- lenz
  df$rna <- ifelse(df[4] >= qua[[2]], "high",ifelse(df[4] <= qua[[1]],"low","med"))
  df = subset(df, rna != "med")
  df$rna <- as.character(df$rna)
  print(table(df$rna))
  lenz_q <- df
  
  # stank
  qua = quantile(stank[,"AMBRA1"],c(q,1-q)) 
  df <- stank
  df$rna <- ifelse(df[4] >= qua[[2]], "high",ifelse(df[4] <= qua[[1]],"low","med"))
  df = subset(df, rna != "med") 
  df$rna <- as.character(df$rna)
  print(table(df$rna))
  stank_q <- df
  
  # combine
  df <- rbind(lenz_q, stank_q)
  print(table(df$rna))
  surv_t <- Surv(df$surv_years, df$Surv_status)
  fit_t <- survfit(surv_t ~ rna , data = df)
  pg.1 = ggsurvplot(fit_t, pval = T,title=colnames(df)[2])
  pp.1 = pg.1$plot
  pp.1 = pp.1 + labs(subtitle = paste("quantile",q,1-q))
  print(pp.1)
dev.off()
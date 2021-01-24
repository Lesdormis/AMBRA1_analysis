options(stringsAsFactors = F)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(openxlsx)
survplot_table = read.csv("survplot_table_TCGAPAN.csv")

fisherIntegration  <- function (vector){
  my_length=length(vector)
  deg_free=my_length*2
  y=-2*sum(log(vector))
  p.val <- 1-pchisq(y, df = deg_free);
  p.val=as.numeric(p.val);
  o = c(y,p.val)
  return(o)
}

# make a table
cancer_types = sort(unique(survplot_table$CANCER_TYPE_ACRONYM))
cancer_tab = as.data.frame(table(survplot_table$CANCER_TYPE_ACRONYM),stringsAsFactors = F)
qua_vec <- seq(from = 0.05, to = 0.5, by = 0.05)
low_qua_name = paste("low",qua_vec,sep  = "_")
high_qua_name = paste("high",1-qua_vec,sep  = "_")
qua_name = paste(high_qua_name,low_qua_name,sep="_vs_")
sur_p_table = as.data.frame(matrix(data = NA, nrow = length(cancer_types),ncol = length(qua_vec),dimnames = list(cancer_types,qua_name)))
p_table = as.data.frame(matrix(data = NA, nrow = length(cancer_types),ncol = length(qua_vec),dimnames = list(cancer_types,qua_name)))

n.cutoff = 80
in_cancers = cancer_tab$Var1[cancer_tab$Freq>=n.cutoff]
# filter with in_cancers
filt_survplot_table = survplot_table %>% filter(CANCER_TYPE_ACRONYM %in% in_cancers)
qua_vec <- seq(from = 0.05, to = 0.5, by = 0.05)
low_qua_name = paste("low",qua_vec,sep  = "_")
high_qua_name = paste("high",1-qua_vec,sep  = "_")
qua_name = paste(high_qua_name,low_qua_name,sep="_vs_")
sur_p_table = as.data.frame(matrix(data = NA, nrow = length(in_cancers),ncol = length(qua_vec),dimnames = list(in_cancers,qua_name)))
p_table = as.data.frame(matrix(data = NA, nrow = length(in_cancers),ncol = length(qua_vec),dimnames = list(in_cancers,qua_name)))
# calculate each p_value
for (i in 1:length(qua_vec)) {
  q=qua_vec[i]
  
  for (ca in in_cancers) {
    temp_surv_table = filt_survplot_table %>% filter(CANCER_TYPE_ACRONYM == ca)
    qua = quantile(temp_surv_table[,"AMBRA1"],c(q,1-q)) 
    df <- temp_surv_table
    df$rna <- ifelse(df[2] >= qua[[2]], "high",ifelse(df[2] <= qua[[1]],"low","mid"))
    df = subset(df, rna != "mid")
    surv_t <- Surv(df$Overall.Survival.years, df$Surv_status)
    fit_t <- survfit(surv_t ~ rna , data = df)
    
    high_m = surv_median(fit_t)[1,"median"]
    low_m = surv_median(fit_t)[2,"median"]
    if (is.na(high_m)&is.na(low_m)) {
      op = "NA"
    } else if(is.na(high_m)&!is.na(low_m)){
      op = "high>low"
    } else if(!is.na(high_m)&is.na(low_m)){
      op = "high<low"
    } else if(high_m > low_m){
      op = "high>low"
    } else if(high_m < low_m){
      op = "high<low"
    }
    # get p value
    p = surv_pvalue(fit_t)[1,2]
    p_table[ca,i] = p
    p_text = paste("p=",round(p,3))
    tet = paste(op,p_text,sep = ";")
    
    sur_p_table[ca,i] = tet
  } 
  #dev.off()
}
# calculate fisher method
fisher.p_05 = fisherIntegration(p_table$high_0.5_vs_low_0.5)[2]
fisher.p_02 = fisherIntegration(p_table$high_0.8_vs_low_0.2)[2]

# overall survival
# 0.5 quantile
pdf(file = "Ext_fig7f.pdf",width = 7,height = 7)
q = 0.5
qua = quantile(filt_survplot_table[,"AMBRA1"],c(0.5,0.5)) 
df <- filt_survplot_table
df$rna <- ifelse(df[2] >= qua[[2]], "high",ifelse(df[2] <= qua[[1]],"low","mid"))
df = subset(df, rna != "mid")
surv_t <- Surv(df$Overall.Survival.years, df$Surv_status)
fit_t <- survfit(surv_t ~ rna , data = df)
pg.1 = ggsurvplot(fit_t, pval = T,title=colnames(df)[2])
os_fit = coxph(surv_t ~ rna, data = df)
p = summary(os_fit)$coefficients[1,5]
pp.1 = pg.1$plot
pp.1 = pp.1 + labs(subtitle = paste("quantile",q,1-q),caption = paste("p-value",p,"fisher.p.value",fisher.p_05))
print(pp.1)
dev.off()

# make table
# filt_survplot_table$rna = ifelse(filt_survplot_table[2] >= qua[[2]], "high",ifelse(filt_survplot_table[2] <= qua[[1]],"low","med"))
# print(table(filt_survplot_table$rna))
# write.xlsx(filt_survplot_table, file = paste0("number_cutoff",n.cutoff,"_OS_survtable_5050.xlsx"))

# 0.2 0.8 quantile
pdf(file = "Fig4b.pdf",width = 7,height = 7)
q = 0.2
qua = quantile(filt_survplot_table[,"AMBRA1"],c(0.2,0.8)) 
df <- filt_survplot_table
df$rna <- ifelse(df[2] >= qua[[2]], "high",ifelse(df[2] <= qua[[1]],"low","mid"))
df = subset(df, rna != "mid")
surv_t <- Surv(df$Overall.Survival.years, df$Surv_status)
fit_t <- survfit(surv_t ~ rna , data = df)
pg.1 = ggsurvplot(fit_t, pval = T,title=colnames(df)[2])
os_fit = coxph(surv_t ~ rna, data = df)
p = summary(os_fit)$coefficients[1,5]
pp.1 = pg.1$plot
pp.1 = pp.1 + labs(subtitle = paste("quantile",q,1-q),caption = paste("p-value",p,"fisher.p.value",fisher.p_02))
print(pp.1)
dev.off()
# make table
# filt_survplot_table$rna = ifelse(filt_survplot_table[2] >= qua[[2]], "high",ifelse(filt_survplot_table[2] <= qua[[1]],"low","mid"))
# print(table(filt_survplot_table$rna))
# write.xlsx(filt_survplot_table, file = paste0("number_cutoff",n.cutoff,"_OS_survtable_2080.xlsx"))




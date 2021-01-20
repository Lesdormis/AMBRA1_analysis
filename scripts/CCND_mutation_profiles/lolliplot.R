library(ggplot2)
options(stringsAsFactors = F)
all_mut <- read.table("CCND_ALL_MUT_Cbio_COSMIC_BL.tsv",header = T, sep = "\t")
ccnd1 <- subset(all_mut, Gene == "CCND1" & Mutation.Type %in% c("Missense_Mutation","Substitution - Missense"))
ccnd2 <- subset(all_mut, Gene == "CCND2" & Mutation.Type %in% c("Missense_Mutation","Substitution - Missense"))
ccnd3 <- subset(all_mut, Gene == "CCND3" & Mutation.Type %in% c("Missense_Mutation","Substitution - Missense"))
# CCND1
freq <- as.data.frame(table(ccnd1$Protein.Position))
data <- data.frame(x = seq(1,295), y = rep(0,295))
data <- merge(data,freq,by.x = "x",by.y = "Var1", all.x = T)
plot <- ggplot(data, aes(x = x, y = Freq)) + 
  geom_segment(aes(x = x, xend = x, y =0, yend = Freq), color = "grey") + 
  geom_point(color = "dark green", size= 3) + scale_x_continuous(trans = "reverse") + 
  coord_flip() + theme_classic()
#plot
ggsave(plot, filename = "Ext_Fig8a.pdf",height = 18,width = 5)

# CCND2
freq <- as.data.frame(table(ccnd2$Protein.Position))
data <- data.frame(x = seq(1,295), y = rep(0,295))
data <- merge(data,freq,by.x = "x",by.y = "Var1", all.x = T)
plot <- ggplot(data, aes(x = x, y = Freq)) + 
  geom_segment(aes(x = x, xend = x, y =0, yend = Freq), color = "grey") + 
  geom_point(color = "dark green", size= 3) + scale_x_continuous(trans = "reverse") + 
  coord_flip() + theme_classic()
plot
ggsave(plot, filename = "Ext_Fig8b.pdf",height = 18,width = 5)

# CCND3
freq <- as.data.frame(table(ccnd3$Protein.Position))
data <- data.frame(x = seq(1,295), y = rep(0,295))
data <- merge(data,freq,by.x = "x",by.y = "Var1", all.x = T)
plot <- ggplot(data, aes(x = x, y = Freq)) + 
  geom_segment(aes(x = x, xend = x, y =0, yend = Freq), color = "grey") + 
  geom_point(color = "dark green", size= 3) + scale_x_continuous(trans = "reverse")  +
  coord_flip() + theme_classic() 
plot
ggsave(plot, filename = "Ext_Fig8c.pdf",height = 18,width = 5)

library(ggplot2)
library(dplyr)
library(readr)

#import files for PCAone plotting
#meta <- read.delim("~/Downloads/model_gf_0.3_gen_adm_500_meta.tsv")
#eigenvec <- read.delim("~/Downloads/model_gf_0.3_gen_adm_500_filt.eigvecs", header=FALSE)
#eigenval <- read.table("~/Downloads/model_gf_0.3_gen_adm_500_filt.eigvals", quote="\"", comment.char="")

meta <- read.delim(snakemake@input[[1]])
eigenvec <- read.delim(snakemake@input[[2]], header=FALSE)
eigenval <- read.table(snakemake@input[[3]], quote="\"", comment.char="")


#merge metadata with pcaone output:
all_info <- cbind(meta, eigenvec)

#make list with % of each PC:
mylist<-c()
for (s in eigenval$V1) {
  print(s / sum(eigenval$V1))
  a<-(s / sum(eigenval$V1))
  mylist <- c(mylist, a)
}

#round numbers:
PC1 <- round(mylist[1]*100, digits=2)
PC2 <- round(mylist[2]*100, digits=2)
PC3 <- round(mylist[3]*100, digits=2)
PC4 <- round(mylist[4]*100, digits=2)


#plot
png(snakemake@output[[1]], width=8, height=6, units='in', res=200, pointsize=4)
ggplot(all_info, aes(x=V1, y=V2, col=time, shape=pop), size = 1, alpha=0.3) +
  geom_point()+
  scale_colour_viridis_c(name = "Generations bp")+
  theme_classic()+                                           
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  labs(x = paste("PC1 (",PC1,"%)", sep = ""), y = paste("PC2 (",PC2,"%)", sep = ""))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"))
dev.off()

png(snakemake@output[[2]], width=8, height=6, units='in', res=200, pointsize=4)
ggplot(all_info, aes(x=V2, y=V3, col=time, shape=pop), size = 1, alpha=0.3) +
  geom_point()+
  scale_colour_viridis_c(name = "Generations bp")+
  theme_classic()+                                           
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  labs(x = paste("PC2 (",PC2,"%)", sep = ""), y = paste("PC3 (",PC3,"%)", sep = ""))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"))
dev.off()

png(snakemake@output[[3]], width=8, height=6, units='in', res=200, pointsize=4)
ggplot(all_info, aes(x=V1, y=V4, col=time, shape=pop), size = 1, alpha=0.3) +
  geom_point()+
  scale_colour_viridis_c(name = "Generations bp")+
  theme_classic()+                                           
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  labs(x = paste("PC1 (",PC1,"%)", sep = ""), y = paste("PC4 (",PC4,"%)", sep = ""))+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(0.5, "cm"))
dev.off()


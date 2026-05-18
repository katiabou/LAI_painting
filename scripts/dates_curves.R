library(ggplot2)
library(dplyr)
library(scales)

target <- read.table(snakemake@input[[1]], quote="\"")
#target <- read.table("~/Downloads/pop_mix:expfit.out", quote="\"") #this is simulations

mean_date <- read.table(snakemake@input[[2]], quote="\"", comment.char="")
#mean_date <- read.table("~/Downloads/model_gf_0.3_source_50_gen_adm_50_source_time_1-mean_dates.txt", quote="\"", comment.char="")

nrmsd <- read.table(snakemake@input[[3]], quote="\"", comment.char="")

#get mean date estimate for that run:
Mean_date = paste('mean date=', mean_date$V7, sep = "")
NRMSD = paste('NRMSD=', nrmsd$V6, sep = "")

#input run details
gf_rate <- paste('gf_rate=', snakemake@params[1], sep = "")
gen_adm <- paste('gen_adm=', snakemake@params[2], sep = "")
source <- paste('source_size=', snakemake@params[3], sep = "")
source_time <- paste('source_time=', snakemake@params[4], sep = "")
seed <- paste('seed=', snakemake@params[5], sep = "")

#title
title_plot <- paste(gf_rate,gen_adm,source,source_time,seed,Mean_date,NRMSD,sep = " | ")

#export plot
png(snakemake@output[[1]], width=8, height=6, units='in', res=200, pointsize=4)
#png('~/Downloads/test.png', width=8, height=6, units='in', res=200, pointsize=4)
ggplot(target) +
  geom_point(aes(V1, V2), shape = 3, colour='purple') +
  geom_line(aes(V1, V3), colour='green4') +
  theme_classic()+
  labs(x='Genetic Distance (cM)', y='Weighted Covariance', title=title_plot)+
  theme(plot.title = element_text(size=8, hjust = 0.5))+
  xlim(0,20)
  #scale_x_continuous(breaks=seq(0, 20, 2), limits = c(0.5,20))+
  #scale_y_continuous(breaks=seq(-0.05, 0.25, 0.05), limits = c(-0.05, 0.25))
dev.off()



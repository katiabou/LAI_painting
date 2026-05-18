#! /usr/bin/env Rscript

#import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(data.table)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

list_of_tracts <- args[1]

a <- str_split(list_of_tracts, pattern=',')


tracts <- c()

for (i in 1:length(a[[1]])){
  tmp = read.delim(a[[1]][i])
  colnames(tmp)[9:11] <- c('gen_adm', 'gf','seed') 
  tracts <- rbind(tracts,tmp)
}

tracts$length <- as.numeric(tracts$length)
tracts$gen_adm <- as.factor(tracts$gen_adm)

options(scipen=999)

#tracts <- read.delim("~/Downloads/model_gf_0.3_gen_adm_50_tracts.tsv")

#plot tracks 
p1 <- tracts %>% filter(source_pop=='pop_c') %>%
  ggplot(aes(x=length/1000000, colour=gen_adm, fill=gen_adm)) + 
  geom_density(alpha = 0.15)+ 
  xlab('Mb')+
  theme_bw()+
  scale_color_brewer(name=expression(t[admix]), palette = "Dark2")+
  scale_fill_brewer(name=expression(t[admix]), palette = "Dark2")+
  ylab('Density')



#get average track length across seeds and individuals
average_per_sample <- tracts %>%
  filter(source_pop=='pop_c') %>%
  group_by(gen_adm) %>% 
  summarise(avg_length_bp=mean(length))%>%
  mutate(avg_length_mb=avg_length_bp/1000000)


#export png
#png('~/Downloads/hist.png', width=7, height=7, units='in', res=200, pointsize=4)
png(args[2], width=8, height=5, units='in', res=200, pointsize=4)
p1
dev.off()

#export tsv
write.table(average_per_sample, file=args[3], quote=FALSE, sep='\t', row.names = FALSE)

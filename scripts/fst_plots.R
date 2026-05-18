#! /usr/bin/env Rscript

#import libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(data.table)
library(tidyr)
library(MetBrewer)

args <- commandArgs(trailingOnly = TRUE)

list_of_fst <- args[1]

a <- str_split(list_of_fst, pattern=',')

fst_all <- c()

for (i in 1:length(a[[1]])){
  tmp = read.delim(a[[1]][i])
  fst_all <- rbind(fst_all,tmp)
}

# Sources, across all gen_adm:
sources <- fst_all %>% filter(comparison=='source1_source2')

s <- ggplot(sources, aes(x=as.factor(source_time), y=fst)) + 
  geom_boxplot(aes(colour = as.factor(gen_adm)))+
  theme_bw()+
  scale_colour_manual(values = met.brewer("Gauguin")[2:7])+
  ggtitle('Source1 vs Source2')+
  theme(plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=8))+
  ylim(0.01,0.04)+
  labs(y= 'Fst populations B and C',  x='Source sampling time (generations bp)', colour=bquote(t[admix]))+
  facet_wrap(vars(gen_adm), scales = "free_x", nrow = 2)



# Source1 and target, across all gen_adm:
sources1 <- fst_all %>% filter(comparison=='source1_target')

s1 <- ggplot(sources1, aes(x=as.factor(source_time), y=fst)) + 
  geom_boxplot(aes(colour = as.factor(gen_adm)))+
  theme_bw()+
  scale_colour_manual(values = met.brewer("Gauguin")[2:7])+
  ggtitle('Source1 vs Target')+
  theme(plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=8))+
  ylim(0.001,0.04)+
  labs(y= 'Fst target and population B',  x='Source sampling time (generations bp)', colour=bquote(t[admix]))+
  facet_wrap(vars(gen_adm), scales = "free_x", nrow = 2)


# Source2 and target, across all gen_adm:
sources2 <- fst_all %>% filter(comparison=='source2_target')

s2 <- ggplot(sources2, aes(x=as.factor(source_time), y=fst)) + 
  geom_boxplot(aes(colour = as.factor(gen_adm)))+
  theme_bw()+
  scale_colour_manual(values = met.brewer("Gauguin")[2:7])+
  ggtitle('Source2 vs Target')+
  theme(plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=8))+
  ylim(0.001,0.04)+
  labs(y= 'Fst target and population C',  x='Source sampling time (generations bp)', colour=bquote(t[admix]))+
  facet_wrap(vars(gen_adm), scales = "free_x", nrow = 2)
  

library("cowplot")
ggdraw() +
  draw_plot(s, x = 0, y = 0.666, width = 1, height = 0.333) +
  draw_plot(s1, x = 0, y = 0.333, width = 1, height = 0.333) +
  draw_plot(s2, x = 0, y = 0, width = 1, height = 0.333) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0), y = c(1, 0.666, 0.333))

#ggsave('~/Desktop/test/example_fig.png', width=10, height=11, units='in', dpi=300)
ggsave(args[2], width=10, height=11, units='in', dpi=300)



# Plot with just sources
s
ggsave(args[3], width=9, height=5, units='in', dpi=300)

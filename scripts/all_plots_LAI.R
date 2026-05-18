#import libraries
library(dplyr)
library(viridis)
library(ggplot2)
library(MetBrewer)
library(stringr)
library(tidyr)
library(data.table)
library(ggpubr)
library(readr)
library("cowplot")

options(scipen=999)

#import prob cutoff
prob <- snakemake@params[[1]]
gf <- snakemake@params[[2]]


#Martin's magic
seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)


#These are the overlap stats (nMCC etc)
stat_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/overlap/"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_.*_.*_tracts_overlap_stats_prob_", prob, ".tsv"), full.names = TRUE)
}))


c <- grep("flare_tracts_overlap_stats", stat_outputs, value = TRUE)
d <- grep("mosaic_tracts_overlap_stats", stat_outputs, value = TRUE)
e <- grep("rfmix_tracts_overlap_stats", stat_outputs, value = TRUE)
f <- grep("simplai_tracts_overlap_stats", stat_outputs, value = TRUE)

##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

resultb <- future_lapply(c(c, d, e, f), function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)


########################### 
#PLOTTING STARTS FROM HERE
########################### 

#resultb are the overlap files (a value for all individuals, will average below)

#fix the method names to capital
resultb$method <- toupper(resultb$method)


#######################
#
# nMCC heatmaps
#
#######################

#group by seed, method, pop, gen_adx, source_size, gf_rate, source_time and average across the 10 individuals, then group by everything apart from seed to get average across seeds
resultb_group <- resultb %>%
  group_by(seed, method, pop, gf_rate, gen_adm, source_pop_size, source_time) %>%
  summarise(F1_mean=mean(F1, na.rm=TRUE), nMCC_mean=mean(nMCC, na.rm=TRUE), recall_mean=mean(recall, na.rm=TRUE), precision_mean=mean(precision, na.rm=TRUE), specificity_mean=mean(specificity, na.rm=TRUE), accuracy_mean=mean(accuracy, na.rm=TRUE)) %>%
  group_by(method, pop, gf_rate, gen_adm, source_pop_size, source_time) %>%
  summarise(F1_mean2=mean(F1_mean, na.rm=TRUE), nMCC_mean2=mean(nMCC_mean, na.rm=TRUE), recall_mean2=mean(recall_mean, na.rm=TRUE), precision_mean2=mean(precision_mean, na.rm=TRUE), specificity_mean2=mean(specificity_mean, na.rm=TRUE), accuracy_mean2=mean(accuracy_mean, na.rm=TRUE))

#make variables as factor
resultb_group$gen_adm <- as.factor(resultb_group$gen_adm)
resultb_group$source_pop_size <- as.factor(resultb_group$source_pop_size)
resultb_group$source_time <- as.factor(resultb_group$source_time)

pop_names <- c(
  'pop_b' = "Population B tracts",
  'pop_c' = "Population C tracts")

#---------------------------------------------------------------------------------------------------------------------------

#all source times against gen_adm for nMCC

#fix the y axis labels to have a subscript
resultb_group$source_time2 <- resultb_group$source_time
resultb_group$source_time2 <- factor(resultb_group$source_time2,levels=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))
levels(resultb_group$source_time2) <- c("-t[5]","-t[4]", "-t[3]", "-t[2]",
                                       "-t[1]","t[admix]","t[1]","t[2]",
                                       "t[3]","t[4]","t[5]")

resultb_group$gen_adm2 <- resultb_group$gen_adm
levels(resultb_group$gen_adm2) <- c("t[admix]==50","t[admix]==100", 
                                   "t[admix]==200","t[admix]==300",
                                   "t[admix]==400","t[admix]==500")

#png('~/Downloads/nMCC_all.png', width=20, height=12, units='in', res=200, pointsize=4)
png(snakemake@output[[1]], width=20, height=12, units='in', res=200, pointsize=4)
resultb_group %>% filter(pop=='pop_c') %>%
  ggplot(aes(source_pop_size, source_time, fill= nMCC_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "Spectral", na.value='grey', breaks=seq(0.5,1,0.1), labels=seq(0.5,1,0.1), limits=c(0.49,1), direction=-1) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='nMCC')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)
dev.off()


#plot other stuff:

#plot Sensitivity (otherwise known as recall) for pop_c
#png('sensitivity_gf_0.3.png', width=20, height=12, units='in', res=200, pointsize=4)
#png(snakemake@output[[2]], width=20, height=12, units='in', res=200, pointsize=4)
p1 <- resultb_group %>% filter(pop=='pop_b') %>%
  ggplot(aes(source_pop_size, source_time, fill= recall_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), direction=-1) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Sensitivity', title='Population B tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)
#dev.off()

p2 <- resultb_group %>% filter(pop=='pop_c') %>%
  ggplot(aes(source_pop_size, source_time, fill= recall_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), direction=-1) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Sensitivity', title='Population C tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)

ggdraw() +
  draw_plot(p1, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p2, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5))

#ggsave('~/Desktop/test/example_fig.png', width=10, height=11, units='in', dpi=300)
ggsave(snakemake@output[[2]], width=16, height=20, units='in', dpi=300)


#plot Specificity for pop_c
#png('Specificity_gf_0.3.png', width=20, height=12, units='in', res=200, pointsize=4)
#png(snakemake@output[[3]], width=20, height=12, units='in', res=200, pointsize=4)
p3 <- resultb_group %>% filter(pop=='pop_b') %>%
  ggplot(aes(source_pop_size, source_time, fill= specificity_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), direction=-1) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Specificity', title='Population B tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)

p4 <- resultb_group %>% filter(pop=='pop_c') %>%
  ggplot(aes(source_pop_size, source_time, fill= specificity_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), direction=-1) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Specificity', title='Population C tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)
#dev.off()

ggdraw() +
  draw_plot(p3, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p4, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5))

#ggsave('~/Desktop/test/example_fig.png', width=10, height=11, units='in', dpi=300)
ggsave(snakemake@output[[3]], width=16, height=20, units='in', dpi=300)



#plot precision
#png('precision_gf_0.3.png', width=20, height=12, units='in', res=200, pointsize=4)
#png(snakemake@output[[4]], width=20, height=12, units='in', res=200, pointsize=4)
p5 <- resultb_group %>% filter(pop=='pop_b') %>%
  ggplot(aes(source_pop_size, source_time, fill= precision_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1)) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Precision', title='Population B tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)

p6 <- resultb_group %>% filter(pop=='pop_c') %>%
  ggplot(aes(source_pop_size, source_time, fill= precision_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1)) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Precision', title='Population C tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)
#dev.off()

ggdraw() +
  draw_plot(p5, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p6, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5))

#ggsave('~/Desktop/test/example_fig.png', width=10, height=11, units='in', dpi=300)
ggsave(snakemake@output[[4]], width=16, height=20, units='in', dpi=300)


#plot Accuracy
#png('accuracy_gf_0.3.png', width=20, height=12, units='in', res=200, pointsize=4)
#png(snakemake@output[[5]], width=20, height=12, units='in', res=200, pointsize=4)
p7 <- resultb_group %>% filter(pop=='pop_b') %>%
  ggplot(aes(source_pop_size, source_time, fill= accuracy_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), direction=-1) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Accuracy', title='Population B tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)
#dev.off()

p8 <- resultb_group %>% filter(pop=='pop_c') %>%
  ggplot(aes(source_pop_size, source_time, fill= accuracy_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(resultb_group$source_time2)))+
  scale_fill_distiller(palette = "YlGnBu", na.value='grey', breaks=seq(0,1,0.2), labels=seq(0,1,0.2), limits=c(0,1), direction=-1) +
  labs(y=bquote(t[sampling]), x='Source sample size', fill='Accuracy', title='Population C tracts')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = 90, size = 13, vjust = 0.5),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=16), 
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        strip.text = element_text(size=17),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(0.7, 'cm'),
        plot.title = element_text(size=16))+
  facet_grid(method~gen_adm2, labeller = label_parsed)

ggdraw() +
  draw_plot(p7, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p8, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5))

#ggsave('~/Desktop/test/example_fig.png', width=10, height=11, units='in', dpi=300)
ggsave(snakemake@output[[5]], width=16, height=20, units='in', dpi=300)


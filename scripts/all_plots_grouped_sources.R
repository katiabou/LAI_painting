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

#args <- commandArgs(trailingOnly = TRUE)

#import prob cutoff
prob1 <- as.numeric(snakemake@params[[1]])
prob <- sprintf("%.1f", prob1) #this is used to input the prob as 0.0 and not 0, to perfectly match the file names

gf <- snakemake@params[[2]]

#Martin's magic
seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)


### Input files for grouped sources

#These are the admixture time and admixture proportions
# time_outputs_s3 <- unlist(lapply(seed_dirs, function(seed_dir) {
#     list.files(path = file.path(seed_dir, "compare_software/files_scenario3/admixture_time_proportion/"), pattern = paste0("model_gf_.*_source_.*_gen_adm_.*_scenario3_.*_adm_time_proportions", ".txt"), full.names = TRUE)
# }))


#These are the overlap stats (nMCC etc) for the grouped run
stat_outputs_grouped_sources <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files_grouped_sources/overlap"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_grouped_sources.*_tracts_overlap_stats_prob_", prob, ".tsv"), full.names = TRUE)
}))


### Input files for initial scenario

#These are the admixture time and admixture proportions
# time_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
#     list.files(path = file.path(seed_dir, "compare_software/files/admixture_time_proportion/"), pattern = paste0("model_gf_.*_source_.*_gen_adm_.*_source_time_.*_.*_adm_time_proportions", ".txt"), full.names = TRUE)
# }))

#These are the overlap stats (nMCC etc)
# Selected time points and source sizes from individual runs
source_time_sub <- c(0.2,0.6,1.4,1.8)
source_sub <- c(4,8,12,16)

pattern <- paste0(
  "model_gf_",paste(gf, collapse = "|"),"_source_(",
  paste(source_sub, collapse = "|"),
  ")_gen_adm_.*_source_time_(",
  paste(source_time_sub, collapse = "|"),
  ")_.*_tracts_overlap_stats_prob_", paste(prob),".tsv"
)

stat_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    #list.files(path = file.path(seed_dir, "compare_software/files/overlap/"), pattern = paste0("model_gf_.*_source_.*_gen_adm_.*_source_time_.*_.*_tracts_overlap_stats_prob_", prob, ".tsv"), full.names = TRUE)
    list.files(path = file.path(seed_dir, "compare_software/files/overlap"), pattern = pattern, full.names = TRUE)
}))

# #These are the percentage of sequence left after filtering for posterior probability
# filt_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
#     list.files(path = file.path(seed_dir, "compare_software/files/tracts/"), pattern = paste0("model_gf_.*_source_.*_gen_adm_.*_source_time_.*_.*_tracts_prob_", prob, "_filt_prop.tsv"), full.names = TRUE)
# }))

a <- grep("grouped_sources_flare_tracts_overlap_stats", stat_outputs_grouped_sources, value = TRUE)
b <- grep("grouped_sources_mosaic_tracts_overlap_stats", stat_outputs_grouped_sources, value = TRUE)
c <- grep("grouped_sources_rfmix_tracts_overlap_stats", stat_outputs_grouped_sources, value = TRUE)
d <- grep("grouped_sources_simplai_tracts_overlap_stats", stat_outputs_grouped_sources, value = TRUE)

e <- grep("flare_tracts_overlap_stats", stat_outputs, value = TRUE)
f <- grep("mosaic_tracts_overlap_stats", stat_outputs, value = TRUE)
g <- grep("rfmix_tracts_overlap_stats", stat_outputs, value = TRUE)
h <- grep("simplai_tracts_overlap_stats", stat_outputs, value = TRUE)


##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

resultb_grouped <- future_lapply(c(a,b,c,d), function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
resultb <- future_lapply(c(e,f,g,h), function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)

#fix the method names to capital
resultb_grouped$method <- toupper(resultb_grouped$method)
resultb$method <- toupper(resultb$method)



#### #### #### #### #### #### 
#### HEAT MAP FIGURES ####
#### #### #### #### #### #### 

#group by seed, method, pop, gen_adx, source_size, gf_rate, source_time and average across the 10 individuals, then group by everyting apart from seed to get average across seeds
result_group <- resultb %>%
  group_by(seed, method, pop, gf_rate, gen_adm, source_pop_size, source_time) %>%
  summarise(F1_mean=mean(F1, na.rm=TRUE), nMCC_mean=mean(nMCC, na.rm=TRUE), recall_mean=mean(recall, na.rm=TRUE), precision_mean=mean(precision, na.rm=TRUE), specificity_mean=mean(specificity, na.rm=TRUE), accuracy_mean=mean(accuracy, na.rm=TRUE)) %>%
  group_by(method, pop, gf_rate, gen_adm, source_pop_size, source_time) %>%
  summarise(F1_mean2=mean(F1_mean, na.rm=TRUE), nMCC_mean2=mean(nMCC_mean, na.rm=TRUE), recall_mean2=mean(recall_mean, na.rm=TRUE), precision_mean2=mean(precision_mean, na.rm=TRUE), specificity_mean2=mean(specificity_mean, na.rm=TRUE), accuracy_mean2=mean(accuracy_mean, na.rm=TRUE))

#result_group$source_time <- as.factor(result_group$source_time)

result_group_s3 <- resultb_grouped %>%
  group_by(seed, method, pop, gf_rate, gen_adm, source_pop_size, source_time) %>%
  summarise(F1_mean=mean(F1, na.rm=TRUE), nMCC_mean=mean(nMCC, na.rm=TRUE), recall_mean=mean(recall, na.rm=TRUE), precision_mean=mean(precision, na.rm=TRUE), specificity_mean=mean(specificity, na.rm=TRUE), accuracy_mean=mean(accuracy, na.rm=TRUE)) %>%
  group_by(method, pop, gf_rate, gen_adm, source_pop_size, source_time) %>%
  summarise(F1_mean2=mean(F1_mean, na.rm=TRUE), nMCC_mean2=mean(nMCC_mean, na.rm=TRUE), recall_mean2=mean(recall_mean, na.rm=TRUE), precision_mean2=mean(precision_mean, na.rm=TRUE), specificity_mean2=mean(specificity_mean, na.rm=TRUE), accuracy_mean2=mean(accuracy_mean, na.rm=TRUE))

result_group_s3$source_time <- 'grouped_sources'
result_group$source_time <- as.factor(result_group$source_time)
#result_group_s3$source_time <- ifelse(result_group_s3$source_time=='grouped_sources', 1245, result_group_s3$source_time)
#result_group_s3$source_time <- as.factor(result_group_s3$source_time)

#merge the two dfs
result_group_merge <- rbind(result_group, result_group_s3)

# make variables as factor
result_group_merge$gen_adm <- as.factor(result_group_merge$gen_adm)
result_group_merge$source_pop_size <- as.factor(result_group_merge$source_pop_size)

#result and resultb need to fix the method names to capital
result_group_merge$method <- toupper(result_group_merge$method)

#subset for what we want:
result_group_merge_sub <- result_group_merge %>%
  filter(source_pop_size %in% c(4,8,12,16) & source_time %in% c('0.2','0.6','1.4','1.8','grouped_sources'))


pop_names <- c(
  'pop_b' = "Population B tracts",
  'pop_c' = "Population C tracts")


#all four source times against gen_adm for nMCC

#fix the y axis labels to have a subscript source_times_1 <- c(0.2,0.6,1.4,1.8) #sample from these specific time periods
result_group_merge_sub$source_time2 <- as.factor(result_group_merge_sub$source_time)
result_group_merge_sub$source_time2 <- factor(result_group_merge_sub$source_time2,levels=c('0.2', '0.6', '1.4', '1.8','grouped_sources'))
levels(result_group_merge_sub$source_time2) <- c("-t[4]","-t[2]","t[2]","t[4]","t[grouped]")

result_group_merge_sub$gen_adm2 <- result_group_merge_sub$gen_adm
levels(result_group_merge_sub$gen_adm2) <- c("t[admix]==50","t[admix]==100", 
                                   "t[admix]==200","t[admix]==300",
                                   "t[admix]==400","t[admix]==500")


#png('nMCC_grouped.png', width=20, height=11, units='in', res=200, pointsize=4)
png(snakemake@output[[1]], width=20, height=11, units='in', res=200, pointsize=4)
result_group_merge_sub %>% filter(pop=='pop_c') %>%
  ggplot(aes(source_pop_size, source_time, fill= nMCC_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(result_group_merge_sub$source_time2)))+
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


#png(snakemake@output[[2]], width=20, height=11, units='in', res=200, pointsize=4)
p1 <- result_group_merge_sub %>% filter(pop=='pop_b') %>%
  ggplot(aes(source_pop_size, source_time, fill= recall_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(result_group_merge_sub$source_time2)))+
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

p2 <- result_group_merge_sub %>% filter(pop=='pop_c') %>%
  ggplot(aes(source_pop_size, source_time, fill= recall_mean2)) + 
  geom_tile() +
  scale_y_discrete(labels = parse(text = levels(result_group_merge_sub$source_time2)))+
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

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

options(scipen=999)

# Martin's magic
seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)

gf <- as.numeric(snakemake@params[[1]])

# These are the admixture times for the grouped
time_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files_grouped_sources/admixture_times"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_grouped_sources_.*_adm_time", ".txt"), full.names = TRUE)
}))

# Selected time points and source sizes from individual runs
source_time_sub <- c(0.2,0.6,1.4,1.8)
source_sub <- c(4,8,12,16)

pattern <- paste0(
  "model_gf_",paste(gf, collapse = "|"),"_source_(",
  paste(source_sub, collapse = "|"),
  ")_gen_adm_.*_source_time_(",
  paste(source_time_sub, collapse = "|"),
  ")_.*_adm_time.txt"
)

# These are the admixture times for the individual time slots [0.2,0.6,1.4,1.8] and sizes [4,8,12,16]
time_outputs2 <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/admixture_times"), pattern = pattern, full.names = TRUE)
}))



# #These are the percentage of sequence left after filtering for posterior probability
# filt_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
#     list.files(path = file.path(seed_dir, "compare_software/files/tracts/"), pattern = paste0("model_gf_.*_source_.*_gen_adm_.*_source_time_.*_.*_tracts_prob_", prob, "_filt_prop.tsv"), full.names = TRUE)
# }))


a <- grep("flare_adm_time", time_outputs, value = TRUE)
b <- grep("mosaic_adm_time", time_outputs, value = TRUE)
c <- grep("dates_adm_time", time_outputs, value = TRUE)
d <- grep("flare_adm_time", time_outputs2, value = TRUE)
e <- grep("mosaic_adm_time", time_outputs2, value = TRUE)
f <- grep("dates_adm_time", time_outputs2, value = TRUE)

##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

result1 <- future_lapply(c(a, b, c, d, e, f), function(path) read.table(path, quote="\"", comment.char="")) %>% do.call(rbind, .)

# add column names to result df
colnames(result1) <- c('method','gf_rate','source_pop_size','gen_adm','source_time','seed', 'inferred_gen_adm')

# Filter for gene flow rate
result <- result1 %>% filter(gf_rate==gf)


################################
#                              #
#  PLOTS FOR INFERRED TIMES #
#                              #
################################

# result and resultb need to fix the method names to capital
result$method <- toupper(result$method)

# make source pop size as factor
result$source_pop_size <- as.factor(result$source_pop_size)

#fix the y axis labels to have a subscript
result$source_time2 <- result$source_time
result$source_time2 <- factor(result$source_time2,levels=c('grouped_sources',0.2,0.6,1.4,1.8))
levels(result$source_time2) <- c("t[grouped]","-t[4]","-t[2]", "t[2]","t[4]")

result$gen_adm2 <- result$gen_adm
result$gen_adm2 <- factor(result$gen_adm2,levels=c(50, 100, 200, 300, 400, 500))
levels(result$gen_adm2) <- c("t[admix]==50","t[admix]==100", 
                                   "t[admix]==200","t[admix]==300",
                                   "t[admix]==400","t[admix]==500")

#png('adm_time_grouped.png', width=15, height=12, units='in', res=200, pointsize=4)
png(snakemake@output[[1]], width=15, height=12, units='in', res=200, pointsize=4)
result %>%
  group_by(gen_adm) %>%
  mutate(inferred_gen_adm = ifelse(inferred_gen_adm >= 2*gen_adm, NA, inferred_gen_adm)) %>%
  ungroup() %>%
  ggplot(aes(x=source_pop_size, y=inferred_gen_adm))+
  geom_point(aes(colour= method), alpha=0.1, size=3)+
  geom_hline(aes(yintercept=gen_adm), linetype="dashed", color = "black", linewidth=0.8)+
  scale_color_manual(values = c("DATES" = "deepskyblue4",
                                "FLARE" = "firebrick3", 
                                "MOSAIC" = "darkolivegreen4"))+  labs(x = 'Source sample size', y = 'Inferred admixture time', colour='Method') +
  theme_bw()+
  theme(axis.text.x=element_text(angle = 90, size = 10, vjust = 0.5),  
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        plot.title = element_text(hjust = 0.5, size=14),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=12),
        legend.key.width= unit(1.5, 'cm'),
        legend.position="right")+
  facet_grid(gen_adm2~source_time2, scales="free", 
             labeller=label_parsed)+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
dev.off()


# export adm date file
write.table(result, file=snakemake@output[[2]], quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)




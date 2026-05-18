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


#import data
# d1 <- read.delim("~/Downloads/seed101/model_gf_0.3_source_8_gen_adm_50_source_time_0-all_simulated_post_prob.txt")
# d2 <- read.delim("~/Downloads/seed101/model_gf_0.3_source_8_gen_adm_200_source_time_0-all_simulated_post_prob.txt")
# d3 <- read.delim("~/Downloads/seed102/model_gf_0.3_source_8_gen_adm_50_source_time_0-all_simulated_post_prob.txt")
# d4 <- read.delim("~/Downloads/seed102/model_gf_0.3_source_8_gen_adm_200_source_time_0-all_simulated_post_prob.txt")

# all <- rbind(d1,d2,d3,d4)

gf <- snakemake@params[[1]]

#Martin's magic
seed_dirs <- dir(path = "output", pattern = "seed_101", full.names = TRUE)


#These are the overlap stats (nMCC etc)
stat_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/posterior_probabilities"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_.*-all_simulated_post_prob.txt"), full.names = TRUE)
}))


##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

resultb <- future_lapply(stat_outputs, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)


#import list of files
#pp <- args[1]

#split into single files
#a <- str_split(pp, pattern=',')

#merge all flare dfs
# all <- c()

# for (i in 1:length(a[[1]])){
#   tmp = read.delim(a[[1]][i])
#   all <- rbind(all,tmp)
# }


# Fix the method names to capital
resultb$method <- toupper(resultb$method)

#plot
png(snakemake@output[[1]], width=15, height=8, units='in', res=200, pointsize=4)
resultb %>%
  ggplot(aes(x=prob_anc1))+
  geom_density(aes(fill=method),alpha=0.4, color=NA)+
  scale_fill_manual(values = c("RFMIX" = "darkgoldenrod2",
                                "FLARE" = "firebrick3", 
                                "MOSAIC" = "darkolivegreen4"))+  
  labs(x = 'Posterior probability population B', y = 'Density', fill='Method') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 90, size = 10, vjust = 0.5),  
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        plot.title = element_text(hjust = 0.5, size=14),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white'),
        legend.key.width= unit(1.5, 'cm'),
        legend.position="right")+
  facet_grid(source~gen_adm, scales="free_y", 
            labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source)))+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  coord_cartesian(ylim = c(0,25))
dev.off()


png(snakemake@output[[2]], width=15, height=8, units='in', res=200, pointsize=4)
resultb %>%
  ggplot(aes(x=prob_anc2))+
  geom_density(aes(fill=method),alpha=0.4, color=NA)+
  scale_fill_manual(values = c("RFMIX" = "darkgoldenrod2",
                               "FLARE" = "firebrick3", 
                               "MOSAIC" = "darkolivegreen4"))+  
  labs(x = 'Posterior probability population C', y = 'Density', fill='Method') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 90, size = 10, vjust = 0.5),  
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14), 
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        plot.title = element_text(hjust = 0.5, size=14),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white'),
        legend.key.width= unit(1.5, 'cm'),
        legend.position="right")+
  facet_grid(source~gen_adm, scales="free_y", 
            labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source)))+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  coord_cartesian(ylim = c(0,25))
dev.off()


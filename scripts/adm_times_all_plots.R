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

# These are the admixture times
time_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/admixture_times"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_.*_.*_adm_time", ".txt"), full.names = TRUE)
}))

a <- grep("flare_adm_time", time_outputs, value = TRUE)
b <- grep("mosaic_adm_time", time_outputs, value = TRUE)
c <- grep("dates_adm_time", time_outputs, value = TRUE)

##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

result1 <- future_lapply(c(a, b, c), function(path) read.table(path, quote="\"", comment.char="")) %>% do.call(rbind, .)

# add column names to result df
colnames(result1) <- c('method','gf_rate','source_pop_size','gen_adm','source_time','seed', 'inferred_gen_adm')

# Filter for gene flow rate
result <- result1 %>% filter(gf_rate==gf)


################################
#                              #
#   PLOTS FOR INFERRED TIMES   #
#                              #
################################

# result and resultb need to fix the method names to capital
result$method <- toupper(result$method)

# make source pop size as factor
result$source_pop_size <- as.factor(result$source_pop_size)

#fix the y axis labels to have a subscript
result$source_time2 <- result$source_time
result$source_time2 <- factor(result$source_time2,levels=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))
levels(result$source_time2) <- c("-t[5]","-t[4]", "-t[3]", "-t[2]",
                                       "-t[1]","t[admix]","t[1]","t[2]",
                                       "t[3]","t[4]","t[5]")

result$gen_adm2 <- result$gen_adm
result$gen_adm2 <- factor(result$gen_adm2,levels=c(50, 100, 200, 300, 400, 500))
levels(result$gen_adm2) <- c("t[admix]==50","t[admix]==100", 
                                   "t[admix]==200","t[admix]==300",
                                   "t[admix]==400","t[admix]==500")


##### PLOT ALL INFERRED ADM TIME AGAINST SOURCE SAMPLE SIZE
#png('adm_time_all.png', width=30, height=15, units='in', res=200, pointsize=4)
png(snakemake@output[[1]], width=30, height=15, units='in', res=200, pointsize=4)
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
        strip.text = element_text(colour = 'white', size=14),
        legend.key.width= unit(1.5, 'cm'),
        legend.position="bottom")+
  facet_grid(gen_adm2~source_time2, scales="free_y", 
             labeller=label_parsed)+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
dev.off()



#png('adm_time_all_sub2.png', width=12, height=7, units='in', res=200, pointsize=4)
png(snakemake@output[[2]], width=12, height=7, units='in', res=200, pointsize=4)
result %>%
  group_by(gen_adm) %>%
  mutate(inferred_gen_adm = ifelse(inferred_gen_adm >= 2*gen_adm, NA, inferred_gen_adm)) %>%
  filter(source_time==0) %>%
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
        legend.position="bottom")+
  facet_wrap(~gen_adm2,
             labeller=label_parsed, ncol = 3, scales = "free")+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
dev.off()


# export adm date file
write.table(result, file=snakemake@output[[3]], quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)


### Plot geom_ribbon

# all scenarios
abcd <- result %>%
  group_by(method, gen_adm, source_pop_size, source_time) %>%
  summarise(
    avg_time = mean(inferred_gen_adm),
    sd  = sd(inferred_gen_adm),
    median_time = median(inferred_gen_adm),
    min = quantile(inferred_gen_adm, probs=0.25),
    max = quantile(inferred_gen_adm, probs=0.75),
    min2 = quantile(inferred_gen_adm, probs=0.1),
    max2 = quantile(inferred_gen_adm, probs=0.9),
    .groups = "drop"
  ) 

# for extreme estimates, cap at 2*adm_time, to make plot easier to read
abcd$max[abcd$max > (abcd$gen_adm * 2)] <- abcd$gen_adm[abcd$max > (abcd$gen_adm * 2)] *2
abcd$median_time[abcd$median_time > (abcd$gen_adm * 2)] <- abcd$gen_adm[abcd$median_time > (abcd$gen_adm * 2)] *2
abcd$max2[abcd$max2 > (abcd$gen_adm * 2)] <- abcd$gen_adm[abcd$max2 > (abcd$gen_adm * 2)] *2
abcd$min[abcd$min > (abcd$gen_adm * 2)] <- abcd$gen_adm[abcd$min > (abcd$gen_adm * 2)] *2
abcd$min2[abcd$min2 > (abcd$gen_adm * 2)] <- abcd$gen_adm[abcd$min2 > (abcd$gen_adm * 2)] *2

#fix the y axis labels to have a subscript
abcd$source_time2 <- abcd$source_time
abcd$source_time2 <- factor(abcd$source_time2,levels=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))
levels(abcd$source_time2) <- c("-t[5]","-t[4]", "-t[3]", "-t[2]",
                                 "-t[1]","t[admix]","t[1]","t[2]",
                                 "t[3]","t[4]","t[5]")

abcd$gen_adm2 <- abcd$gen_adm
abcd$gen_adm2 <- factor(abcd$gen_adm2,levels=c(50, 100, 200, 300, 400, 500))
levels(abcd$gen_adm2) <- c("t[admix]==50","t[admix]==100", 
                             "t[admix]==200","t[admix]==300",
                             "t[admix]==400","t[admix]==500")


##### PLOT ALL INFERRED ADM TIME AGAINST SOURCE SAMPLE SIZE
#png('adm_time_all.png', width=30, height=15, units='in', res=200, pointsize=4)
png(snakemake@output[[4]], width=30, height=15, units='in', res=200, pointsize=4)
abcd %>%
  ggplot()+
  geom_hline(aes(yintercept=gen_adm), linetype="dashed", color = "black", linewidth=0.8)+
  geom_point(aes(x = source_pop_size, y = median_time, colour=method), size = 3, alpha=0) +
  geom_line(aes(x = as.numeric(source_pop_size), y = median_time, colour=method)) +
  geom_ribbon(aes(x = as.numeric(source_pop_size), ymax = max, ymin = min, fill = method), alpha = 0.3) +
  geom_ribbon(aes(x = as.numeric(source_pop_size), ymax = max2, ymin = min2, fill = method), alpha = 0.3) +
  scale_fill_manual(values = c("DATES" = "deepskyblue4",
                               "FLARE" = "firebrick3", 
                               "MOSAIC" = "darkolivegreen4"))+
  scale_color_manual(values = c("DATES" = "deepskyblue4",
                                "FLARE" = "firebrick3", 
                                "MOSAIC" = "darkolivegreen4"))+  labs(x = 'Source sample size', y = 'Inferred admixture time', fill='Method') +
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
        legend.position="bottom")+
  facet_grid(gen_adm2~source_time2, scales="free_y", 
             labeller=label_parsed)+
  guides(fill = guide_legend(override.aes = list(alpha = 1)), colour = 'none')
dev.off()



#png('adm_time_all_sub2.png', width=12, height=7, units='in', res=200, pointsize=4)
png(snakemake@output[[5]], width=12, height=7, units='in', res=200, pointsize=4)
abcd %>%
  filter(source_time==0) %>%
  ggplot()+
  geom_hline(aes(yintercept=gen_adm), linetype="dashed", color = "black", linewidth=0.8)+
  geom_point(aes(x = source_pop_size, y = median_time, colour=method), size = 3, alpha=0) +
  geom_line(aes(x = as.numeric(source_pop_size), y = median_time, colour=method)) +
  geom_ribbon(aes(x = as.numeric(source_pop_size), ymax = max, ymin = min, fill = method), alpha = 0.3) +
  geom_ribbon(aes(x = as.numeric(source_pop_size), ymax = max2, ymin = min2, fill = method), alpha = 0.3) +
  scale_fill_manual(values = c("DATES" = "deepskyblue4",
                               "FLARE" = "firebrick3", 
                               "MOSAIC" = "darkolivegreen4"))+
  scale_color_manual(values = c("DATES" = "deepskyblue4",
                                "FLARE" = "firebrick3", 
                                "MOSAIC" = "darkolivegreen4"))+  labs(x = 'Source sample size', y = 'Inferred admixture time', fill='Method') +
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
        legend.position="bottom")+
  facet_wrap(~gen_adm2, scales="free_y", 
             labeller=label_parsed)+
  guides(fill = guide_legend(override.aes = list(alpha = 1)), colour = 'none')
dev.off()


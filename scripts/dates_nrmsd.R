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

#Martin's magic
seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)

gf <- as.numeric(snakemake@params[[1]])

#These are the mean date estimates 
stat_outputs_nrmsd <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "dates/output/"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_.*-nrmsd_dates.txt"), full.names = TRUE)
}))

stat_outputs_date <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/admixture_times/"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_.*_dates_adm_time.txt"), full.names = TRUE)
}))


a <- grep("-nrmsd_dates.txt", stat_outputs_nrmsd, value = TRUE)
b <- grep("_dates_adm_time.txt", stat_outputs_date, value = TRUE)


##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

resulta <- future_lapply(a, function(path) read.table(path, quote="\"", comment.char="")) %>% do.call(rbind, .)
resultb <- future_lapply(b, function(path) read.table(path, quote="\"", comment.char="")) %>% do.call(rbind, .)

#fix column names 
colnames(resulta) <- c('seed','gf_rate','source_sample_size','gen_adm','source_time','nrmsd')
colnames(resultb) <- c('method','gf_rate','source_sample_size','gen_adm','source_time', 'seed', 'dates_estimate')


# Merge both datasets
all_dates <- merge(resulta, resultb, by=c('seed','gf_rate','source_sample_size','gen_adm', 'source_time'))
all_dates$source_sample_size <- as.factor(all_dates$source_sample_size)
all_dates$nrmsd <- as.numeric(all_dates$nrmsd)


#fix the y axis labels to have a subscript
all_dates$source_time2 <- all_dates$source_time
all_dates$source_time2 <- factor(all_dates$source_time2,levels=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))
levels(all_dates$source_time2) <- c("-t[5]","-t[4]", "-t[3]", "-t[2]",
                                       "-t[1]","t[admix]","t[1]","t[2]",
                                       "t[3]","t[4]","t[5]")

all_dates$gen_adm2 <- all_dates$gen_adm
all_dates$gen_adm2 <- factor(all_dates$gen_adm2,levels=c(50, 100, 200, 300, 400, 500))
levels(all_dates$gen_adm2) <- c("t[admix]==50","t[admix]==100", 
                                   "t[admix]==200","t[admix]==300",
                                   "t[admix]==400","t[admix]==500")


# Plot 
#png('adm_time_nrmsd_dates.png', width=30, height=15, units='in', res=200, pointsize=4)
png(snakemake@output[[1]], width=30, height=15, units='in', res=200, pointsize=4)
all_dates %>%
  group_by(gen_adm2) %>%
  mutate(dates_estimate = ifelse(dates_estimate >= 2*gen_adm, NA, dates_estimate)) %>%
  ungroup() %>%
  mutate(nrmsd = ifelse(nrmsd >= 0.7, NA, nrmsd)) %>%
  ggplot(aes(x=source_sample_size, y=dates_estimate))+
  geom_point(aes(colour= nrmsd), alpha=0.6, size=3)+
  scale_colour_distiller(palette = "Spectral", na.value='grey', limits=c(0,0.7), breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
  geom_hline(aes(yintercept=gen_adm), linetype="dashed", color = "black", linewidth=0.8)+
  labs(y=expression("Inferred " *t[admix]), x='Source sample size', colour='NRMSD')+
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
        legend.key.width= unit(1.2, 'cm'),
        legend.key.height= unit(1.5, 'cm'),
        legend.position="right")+
  facet_grid(gen_adm2~source_time2, scales="free_y", labeller=label_parsed) 
  #guides(color = guide_legend(override.aes = list(alpha = 1)))
dev.off()








#### Old code


#Martin's magic
# seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)


# #These are the mean date estimates 
# stat_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
#     list.files(path = file.path(seed_dir, "dates/output/"), pattern = paste0("model_gf_.*_source_.*_gen_adm_.*_source_time_.*-nrmsd_dates.txt"), full.names = TRUE)
# }))


# a <- grep("-nrmsd_dates.txt", stat_outputs, value = TRUE)


# ##### MAKE IMPORTING FILES FASTER
# library(future.apply)

# plan(multicore, workers = availableCores())

# resultb <- future_lapply(a, function(path) read.table(path, quote="\"", comment.char="")) %>% do.call(rbind, .)

# #fix column names 
# colnames(resultb) <- c('seed','gf_rate','source_sample_size','gen_adm','source_time','dates_estimate')

# resultb$source_time <- as.numeric(resultb$source_time)

# #start plotting

# #subset for only source time=1:
# # g <- resultb %>%
# #   group_by(gen_adm, source_time, source_sample_size) %>%
# #   dplyr::summarise(mean.dates = mean(dates_estimate, na.rm = TRUE),
# #             sd.dates = sd(dates_estimate, na.rm = TRUE),
# #             n.dates = n()) %>%
# #   mutate(se.dates = sd.dates / sqrt(n.dates),
# #          lower.ci.dates = mean.dates - qt(1 - (0.05 / 2), n.dates - 1) * se.dates,
# #          upper.ci.dates = mean.dates + qt(1 - (0.05 / 2), n.dates - 1) * se.dates)


# # png(snakemake@output[[1]], width=8, height=15, units='in', res=200, pointsize=4)
# # g %>% filter(source_time==1)%>% 
# #   ggplot(aes(x = as.factor(source_sample_size), y = (mean.dates))) +
# #   geom_point() +
# #   geom_hline(aes(yintercept=gen_adm), linetype="dashed", color = "red", size=1)+
# #   labs(y='DATES estimated admixture time (generations)', x='Source sample size')+
# #   geom_errorbar(aes(ymin =  lower.ci.dates, ymax = upper.ci.dates), width=.1)+
# #   theme_bw()+
# #   facet_grid(rows = vars(gen_adm), scales="free")+
# #   ggtitle(paste("Source time=1"))+
# #   theme(plot.title = element_text(size=10, hjust = 0.5))
# #   #facet_wrap(.~gen_adm, scales="free")
# # dev.off()



# #raw data no summary metrics
# #subset for only source time=1:

# dates_means <- resultb %>% 
#   group_by(gen_adm, source_time, source_sample_size) %>%
#   dplyr::summarise(mean.dates = mean(dates_estimate, na.rm = TRUE))

# # png(snakemake@output[[1]], width=8, height=15, units='in', res=200, pointsize=4)
# # resultb %>% filter(source_time==1) %>% 
# #   ggplot(aes(x = as.factor(source_sample_size), y = dates_estimate)) +
# #   geom_point(alpha=0.7) +
# #   geom_point(data = dates_means %>% filter(source_time==1), mapping=aes(x = as.factor(source_sample_size), y = mean.dates), col="mediumturquoise", shape=19)+
# #   geom_hline(aes(yintercept=gen_adm), linetype="dashed", color = "red", size=0.5)+
# #   labs(y='DATES estimated admixture time (generations)', x='Source sample size')+
# #   theme_bw()+
# #   facet_grid(rows = vars(gen_adm), scales="free")+
# #   ggtitle(paste("Source time=1"))+
# #   theme(plot.title = element_text(size=10, hjust = 0.5))+
# #   scale_y_continuous(trans='log10')
# # dev.off()

# png(snakemake@output[[1]], width=8, height=15, units='in', res=200, pointsize=4)
# resultb %>% filter(source_time==0)%>% 
#   ggplot(aes(x = as.factor(source_sample_size), y = dates_estimate, fill=as.factor(gen_adm))) +
#   geom_violin(trim=FALSE)+
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="black")+
#   scale_fill_manual(values=met.brewer("Peru1", 6), name = "Admixture time (gen)")+
#   geom_hline(aes(yintercept=0.7), linetype="dashed", color = "red", size=0.5)+
#   labs(y='DATES NRMSD', x='Source sample size')+
#   theme_bw()+
#   facet_grid(rows = vars(gen_adm), scales="free")+
#   ggtitle(paste("Source time=0"))+
#   theme(plot.title = element_text(size=10, hjust = 0.5))+
#   scale_y_continuous(trans='log10')
# dev.off()


# png(snakemake@output[[2]], width=8, height=15, units='in', res=200, pointsize=4)
# resultb %>% filter(source_time==1)%>% 
#   ggplot(aes(x = as.factor(source_sample_size), y = dates_estimate, fill=as.factor(gen_adm))) +
#   geom_violin(trim=FALSE)+
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="black")+
#   scale_fill_manual(values=met.brewer("Peru1", 6), name = "Admixture time (gen)")+
#   geom_hline(aes(yintercept=0.7), linetype="dashed", color = "red", size=0.5)+
#   labs(y='DATES NRMSD', x='Source sample size')+
#   theme_bw()+
#   facet_grid(rows = vars(gen_adm), scales="free")+
#   ggtitle(paste("Source time=1"))+
#   theme(plot.title = element_text(size=10, hjust = 0.5))+
#   scale_y_continuous(trans='log10')
# dev.off()

# png(snakemake@output[[3]], width=8, height=15, units='in', res=200, pointsize=4)
# resultb %>% filter(source_time==2)%>% 
#   ggplot(aes(x = as.factor(source_sample_size), y = dates_estimate, fill=as.factor(gen_adm))) +
#   geom_violin(trim=FALSE)+
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="black")+
#   scale_fill_manual(values=met.brewer("Peru1", 6), name = "Admixture time (gen)")+
#   geom_hline(aes(yintercept=0.7), linetype="dashed", color = "red", size=0.5)+
#   labs(y='DATES NRMSD', x='Source sample size')+
#   theme_bw()+
#   facet_grid(rows = vars(gen_adm), scales="free")+
#   ggtitle(paste("Source time=2"))+
#   theme(plot.title = element_text(size=10, hjust = 0.5))+
#   scale_y_continuous(trans='log10')
# dev.off()

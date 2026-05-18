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

#Import parameters
sequence_length <- as.numeric(snakemake@params[[1]])
#sequence_length <- as.numeric(300000000)

gf <- as.numeric(snakemake@params[[2]])
#gf <- 0.3

#Martin's magic
seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)

#These are the tracts with pp cutoffs 
pp_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files_pp_cutoff/tracts/"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_0_.*_tracts_probcutoff_.*.tsv"), full.names = TRUE)
}))

c <- grep("probcutoff_0.6.tsv", pp_outputs, value = TRUE)
d <- grep("probcutoff_0.7.tsv", pp_outputs, value = TRUE)
e <- grep("probcutoff_0.8.tsv", pp_outputs, value = TRUE)
f <- grep("probcutoff_0.9.tsv", pp_outputs, value = TRUE)
g <- grep("probcutoff_0.999.tsv", pp_outputs, value = TRUE)


#These are the tracts with no pp cutoffs 
source_sub <- c(4,10,20,100)

pattern <- paste0(
  "model_gf_",paste(gf, collapse = "|"),"_source_(",
  paste(source_sub, collapse = "|"),
  ")_gen_adm_.*_source_time_0_.*_tracts_prob_0.0.tsv"
)

outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/tracts/"), pattern = pattern, full.names = TRUE)
}))

c2 <- grep("flare_tracts_prob_0.0.tsv", outputs, value = TRUE)
d2 <- grep("mosaic_tracts_prob_0.0.tsv", outputs, value = TRUE)
e2 <- grep("rfmix_tracts_prob_0.0.tsv", outputs, value = TRUE)
f2 <- grep("simplai_tracts_prob_0.0.tsv", outputs, value = TRUE)


# nMCC estimates for different pp cutoffs
nmcc_outputs_pp <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files_pp_cutoff/overlap/"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_0_.*_tracts_overlap_stats_probcutoff_.*.tsv"), full.names = TRUE)
}))

c3 <- grep("probcutoff_0.6.tsv", nmcc_outputs_pp, value = TRUE)
d3 <- grep("probcutoff_0.7.tsv", nmcc_outputs_pp, value = TRUE)
e3 <- grep("probcutoff_0.8.tsv", nmcc_outputs_pp, value = TRUE)
f3 <- grep("probcutoff_0.9.tsv", nmcc_outputs_pp, value = TRUE)
g3 <- grep("probcutoff_0.999.tsv", nmcc_outputs_pp, value = TRUE)


# nMCC estimates for no pp cutoffs
pattern <- paste0(
  "model_gf_",paste(gf, collapse = "|"),"_source_(",
  paste(source_sub, collapse = "|"),
  ")_gen_adm_.*_source_time_0_.*_tracts_overlap_stats_prob_0.0.tsv"
)

nmcc_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/overlap/"), pattern = pattern, full.names = TRUE)
}))

c4 <- grep("flare_tracts_overlap_stats_prob_0.0.tsv", nmcc_outputs, value = TRUE)
d4 <- grep("mosaic_tracts_overlap_stats_prob_0.0.tsv", nmcc_outputs, value = TRUE)
e4 <- grep("rfmix_tracts_overlap_stats_prob_0.0.tsv", nmcc_outputs, value = TRUE)
f4 <- grep("simplai_tracts_overlap_stats_prob_0.0.tsv", nmcc_outputs, value = TRUE)



##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

# pp files
c_df <- future_lapply(c, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
d_df <- future_lapply(d, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
e_df <- future_lapply(e, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
f_df <- future_lapply(f, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
g_df <- future_lapply(g, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)

# add pp info to dfs
c_df$pp <- 0.6
d_df$pp <- 0.7
e_df$pp <- 0.8
f_df$pp <- 0.9
g_df$pp <- 0.999

# merge dfs
all_pp <- rbind(c_df, d_df, e_df, f_df, g_df)


# no pp files
c2_df <- future_lapply(c2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
d2_df <- future_lapply(d2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
e2_df <- future_lapply(e2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
f2_df <- future_lapply(f2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)

# merge dfs
all <- rbind(c2_df, d2_df, e2_df, f2_df)
all$pp <- 0.0 


# nMCC pp files
c3_df <- future_lapply(c3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
d3_df <- future_lapply(d3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
e3_df <- future_lapply(e3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
f3_df <- future_lapply(f3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
g3_df <- future_lapply(g3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)

c3_df$pp <- 0.6
d3_df$pp <- 0.7
e3_df$pp <- 0.8
f3_df$pp <- 0.9
g3_df$pp <- 0.999

# merge dfs
nmcc_pp <- rbind(c3_df, d3_df, e3_df, f3_df, g3_df)

# nMCC pp files
c3_df <- future_lapply(c3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
d3_df <- future_lapply(d3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
e3_df <- future_lapply(e3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
f3_df <- future_lapply(f3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
g3_df <- future_lapply(g3, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)

# merge dfs
nmcc <- rbind(c3_df, d3_df, e3_df, f3_df, g3_df)
nmcc$pp <- 0.0 


# Merge the tracts for all pp cutoffs
all_tract <- rbind(all_pp, all)


# Group by and then estimate adm proportions
df2 <- all_tract %>% 
       filter(source_pop!='pop_no') %>%
       filter(!(type=='simplai' & pp!=0.0)) %>%
       separate(range, c("right","left"), sep="-", convert = TRUE) %>%
       mutate(length=(left-right)) %>%
       group_by(seed, name, source_pop, source_pop_size, pp, type, gen_adm, gf_rate, source_time) %>%
       summarise(proportions = (sum(length))/(2*sequence_length)) 

# df3 <- df2 %>%
#        group_by(seed, name, source_pop, pp, type, source_pop_size, gen_adm) %>%
#        summarise(prop_all = mean(proportions)) 


### Test plot with no nMCC

# Make pp a factor
df2$pp <- as.factor(df2$pp)

# Fix the method names to capital
df2$type <- toupper(df2$type)

#This is without nmcc info
# png('test_pp2.png', width=16, height=9, units='in', res=200, pointsize=4)
# df2 %>% filter(source_pop =='pop_c') %>%
#   ggplot(aes(pp, proportions)) + 
#   geom_hline(aes(yintercept = (gf)), col = "black", linetype="dashed", size=0.3) +
#   geom_point(aes(colour= type), alpha=0.5, size=3)+
#   #scale_shape_manual(values=c(4, 16))+
#   scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
#                                 "FLARE" = "firebrick3", 
#                                 "MOSAIC" = "darkolivegreen4",
#                                 "SIMPLAI" = "#74A089"))+  
#   labs(x = 'Posterior probability cutoff', y = 'Admixture proportions pop_c', colour='Method') +
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text.x=element_text(size = 10),  
#         axis.text.y=element_text(size=10),
#         axis.title.x=element_text(size=14), 
#         axis.title.y=element_text(size=14), 
#         legend.text = element_text(size=14),
#         legend.title = element_text(size=14),
#         plot.title = element_text(hjust = 0.5, size=14),
#         panel.grid.minor = element_blank(),
#         strip.background =element_rect(fill="gray28"),
#         strip.text = element_text(colour = 'white'),
#         legend.key.width= unit(1.5, 'cm'),
#         legend.position="right")+
#   facet_grid(source_pop_size~gen_adm,
#             labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source_pop_size)))+
#   guides(color = guide_legend(override.aes = list(alpha = 1)))+
#   ylim(0,1)

# dev.off()



# Merge the nMCC files with pp and no pp
all_nmcc <- rbind(nmcc_pp, nmcc)

filt_nmcc <- all_nmcc %>% 
       filter(pop!='pop_no') %>% 
       filter(!(method=='simplai' & pp!=0.0)) %>%
       select(sample, method, pop, gf_rate, source_pop_size, gen_adm, source_time, seed, pp, nMCC)

# Make pp a factor
filt_nmcc$pp <- as.factor(filt_nmcc$pp)

# Fix the method names to capital
filt_nmcc$method <- toupper(filt_nmcc$method)


# rename columns to match %sequence df
names(filt_nmcc)[names(filt_nmcc) == 'method'] <- 'type'
names(filt_nmcc)[names(filt_nmcc) == 'sample'] <- 'name'
names(filt_nmcc)[names(filt_nmcc) == 'pop'] <- 'source_pop'

# Merge two df based on common columns
merged_df <- merge(filt_nmcc, df2, by=c("seed","name","type","source_pop_size", "gen_adm", "pp", "source_pop", "gf_rate", "source_time"), all.x=TRUE)

# replace NA in prop_length column with 0
merged_df$proportions <- ifelse(is.na(merged_df$proportions), 0, merged_df$proportions)


merged_df$nMCC_cutoff <- ifelse(merged_df$nMCC>0.8, '>0.8', '<0.8')
merged_df <- as.data.table(merged_df)

# create order
my_order <- c(">0.8", "<0.8")

# reversed to get correct ordering
merged_df[, order := match(nMCC_cutoff, rev(my_order))]

# sort the order
setorder(merged_df, proportions, order)


# plot using nMCC cutoff:
png(snakemake@output[[1]], width=16, height=8, units='in', res=200, pointsize=4)
#png('test_pp3.png', width=16, height=8, units='in', res=200, pointsize=4)
merged_df %>% filter(source_pop =='pop_b') %>%
  ggplot(aes(pp, proportions)) + 
  geom_hline(aes(yintercept = (1-gf)), col = "black", linetype="dashed", size=0.3) +
  geom_point(aes(colour= type, shape=nMCC_cutoff), alpha=0.2, size=3)+
  scale_shape_manual(values=c(4, 16))+
  scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
                                "FLARE" = "firebrick3", 
                                "MOSAIC" = "darkolivegreen4",
                                "SIMPLAI" = "#74A089"))+  
  labs(x = 'Posterior probability cutoff', y = 'Admixture proportions population B', colour='Method', shape='nMCC') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 10),  
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
  facet_grid(source_pop_size~gen_adm,
            labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source_pop_size)))+
  guides(color = guide_legend(override.aes = list(alpha = 1)), shape = guide_legend(override.aes = list(alpha = 1)))+
  ylim(0,1)
dev.off()


png(snakemake@output[[2]], width=16, height=8, units='in', res=200, pointsize=4)
#png('test_pp3.png', width=16, height=8, units='in', res=200, pointsize=4)
merged_df %>% filter(source_pop =='pop_c') %>%
  ggplot(aes(pp, proportions)) + 
  geom_hline(aes(yintercept = (gf)), col = "black", linetype="dashed", size=0.3) +
  geom_point(aes(colour= type, shape=nMCC_cutoff), alpha=0.2, size=3)+
  scale_shape_manual(values=c(4, 16))+
  scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
                                "FLARE" = "firebrick3", 
                                "MOSAIC" = "darkolivegreen4",
                                "SIMPLAI" = "#74A089"))+  
  labs(x = 'Posterior probability cutoff', y = 'Admixture proportions population C', colour='Method', shape='nMCC') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 10),  
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
  facet_grid(source_pop_size~gen_adm,
            labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source_pop_size)))+
  guides(color = guide_legend(override.aes = list(alpha = 1)), shape = guide_legend(override.aes = list(alpha = 1)))+
  ylim(0,0.5)
dev.off()




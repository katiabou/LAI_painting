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

gf <- snakemake@params[[2]]
#gf <- 0.3

#Martin's magic
seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)

#These are the tracts with pp cutoffs 
pp_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files_pp_cutoff/tracts/"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_.*_.*_tracts_probcutoff_.*.tsv"), full.names = TRUE)
}))

c <- grep("probcutoff_0.6.tsv", pp_outputs, value = TRUE)
d <- grep("probcutoff_0.7.tsv", pp_outputs, value = TRUE)
e <- grep("probcutoff_0.8.tsv", pp_outputs, value = TRUE)
f <- grep("probcutoff_0.9.tsv", pp_outputs, value = TRUE)
g <- grep("probcutoff_0.999.tsv", pp_outputs, value = TRUE)



##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

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

#take out regions with low pp (pop_no regions) and then sum up the remaining bp
tracts_pass <- all_pp %>%
               filter(source_pop!='pop_no') %>%
               separate(range, c("start","end")) %>%
               mutate(length_bp=as.numeric(end)-as.numeric(start))

tracts_pass_prop <- tracts_pass %>%
                    group_by(seed, type, name, source_pop_size, gen_adm, pp) %>%
                    summarise(length_sum = sum(length_bp)) %>%
                    mutate(prop_length = (length_sum/(2*sequence_length))*100) #%>%
                    #summarise(avg_prop = mean(prop_length)) %>%
                    #mutate(avg_prop = ifelse(is.nan(avg_prop), 0, avg_prop))
                    #mutate(gf_rate=GF_RATE, source_pop_size=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIME, seed=SEED, method='FLARE')


# Make pp a factor
tracts_pass_prop$pp <- as.factor(tracts_pass_prop$pp)

# Fix the method names to capital
tracts_pass_prop$type <- toupper(tracts_pass_prop$type)


#plot per gen_adm and source_pop_size, coloured by type, but remove SIMPLAI since no pp are reported!!!
# png(snakemake@output[[1]], width=15, height=8, units='in', res=200, pointsize=4)
# #png('test_pp.png', width=15, height=8, units='in', res=200, pointsize=4)
# tracts_pass_prop %>% filter(type!='SIMPLAI') %>%
#   ggplot(aes(pp, prop_length)) + 
#   geom_point(aes(colour= type), alpha=0.1, size=3)+
#   scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
#                                 "FLARE" = "firebrick3", 
#                                 "MOSAIC" = "darkolivegreen4"))+  
#   labs(x = 'Posterior probability cutoff', y = '% of sequence retained', colour='Method') +
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
#   guides(color = guide_legend(override.aes = list(alpha = 1)))
# dev.off()



# Add nMCC estimates for different pp cutoffs
nmcc_outputs <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files_pp_cutoff/overlap/"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_.*_.*_tracts_overlap_stats_probcutoff_.*.tsv"), full.names = TRUE)
}))

c2 <- grep("probcutoff_0.6.tsv", nmcc_outputs, value = TRUE)
d2 <- grep("probcutoff_0.7.tsv", nmcc_outputs, value = TRUE)
e2 <- grep("probcutoff_0.8.tsv", nmcc_outputs, value = TRUE)
f2 <- grep("probcutoff_0.9.tsv", nmcc_outputs, value = TRUE)
g2 <- grep("probcutoff_0.999.tsv", nmcc_outputs, value = TRUE)

c2_df <- future_lapply(c2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
d2_df <- future_lapply(d2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
e2_df <- future_lapply(e2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
f2_df <- future_lapply(f2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
g2_df <- future_lapply(g2, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)

# add pp info to dfs
c2_df$pp <- 0.6
d2_df$pp <- 0.7
e2_df$pp <- 0.8
f2_df$pp <- 0.9
g2_df$pp <- 0.999

# merge dfs
all_pp2 <- rbind(c2_df, d2_df, e2_df, f2_df, g2_df)

# Make pp a factor
all_pp2$pp <- as.factor(all_pp2$pp)

# Fix the method names to capital
all_pp2$method <- toupper(all_pp2$method)


# Select one of the two ancestries (nMCC estimates are the same for both)
all_pp2_popc <- all_pp2 %>% filter(pop=='pop_c')

# rename columns to match %sequence df
names(all_pp2_popc)[names(all_pp2_popc) == 'method'] <- 'type'
names(all_pp2_popc)[names(all_pp2_popc) == 'sample'] <- 'name'

# Merge two df based on common columns
merged_df <- merge(all_pp2_popc, tracts_pass_prop, by=c("seed","name","type","source_pop_size", "gen_adm", "pp"), all.x=TRUE)

# replace NA in prop_length column with 0
merged_df$prop_length <- ifelse(is.na(merged_df$prop_length), 0, merged_df$prop_length)

# plot using nMCC cutoff:
# png('test_pp.png', width=15, height=8, units='in', res=200, pointsize=4)
# merged_df %>% filter(type!='SIMPLAI') %>%
#   ggplot(aes(pp, prop_length)) + 
#   geom_point(aes(colour= type), alpha=0.1, size=3)+
#   scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
#                                 "FLARE" = "firebrick3", 
#                                 "MOSAIC" = "darkolivegreen4"))+  
#   labs(x = 'Posterior probability cutoff', y = '% of sequence retained', colour='Method') +
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
#   guides(color = guide_legend(override.aes = list(alpha = 1)))
# dev.off()


merged_df$nMCC_cutoff <- ifelse(merged_df$nMCC>0.8, '>0.8', '<0.8')
merged_df <- as.data.table(merged_df)

# create order
my_order <- c(">0.8", "<0.8")

# reversed to get correct ordering
merged_df[, order := match(nMCC_cutoff, rev(my_order))]

# sort the order
setorder(merged_df, prop_length, order)


# plot using nMCC cutoff:
png(snakemake@output[[1]], width=15, height=8, units='in', res=200, pointsize=4)
#png('test_pp3.png', width=15, height=8, units='in', res=200, pointsize=4)
merged_df %>% filter(type!='SIMPLAI') %>%
  ggplot(aes(pp, prop_length)) + 
  geom_point(aes(colour= type, shape=nMCC_cutoff), alpha=0.2, size=3)+
  scale_shape_manual(values=c(4, 16))+
  scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
                                "FLARE" = "firebrick3", 
                                "MOSAIC" = "darkolivegreen4"))+  
  labs(x = 'Posterior probability cutoff', y = '% of sequence retained', colour='Method', shape='nMCC') +
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
  guides(color = guide_legend(override.aes = list(alpha = 1)), shape = guide_legend(override.aes = list(alpha = 1)))
dev.off()




# The code below is for adding standard error bars to the plot, but they end barely being there

# tracts_pass_prop <- tracts_pass %>%
#                     group_by(seed, type, name, source_pop_size, gen_adm, pp) %>%
#                     summarise(length_sum = sum(length_bp)) %>%
#                     mutate(prop_length = (length_sum/(2*sequence_length))*100) 

# tracts_pass_prop_avg <- tracts_pass_prop %>% 
#                         group_by(type, source_pop_size, gen_adm, pp) %>%
#                         summarise(avg_freq=mean(prop_length), sd_delta=sd(prop_length), n_delta=n()) %>%
#                         mutate(se_delta=sd_delta / sqrt(n_delta),
#                         lower_ci=avg_freq - qt(1 - (0.05 / 2), n_delta - 1) * se_delta,
#                         upper_ci=avg_freq + qt(1 - (0.05 / 2), n_delta - 1) * se_delta)

# # Make pp a factor
# tracts_pass_prop_avg$pp <- as.factor(tracts_pass_prop_avg$pp)

# # Fix the method names to capital
# tracts_pass_prop_avg$type <- toupper(tracts_pass_prop_avg$type)


# #plot per gen_adm and source_pop_size, coloured by type, but remove SIMPLAI since no pp are reported!!!
# png(snakemake@output[[1]], width=15, height=8, units='in', res=200, pointsize=4)
# #png('test_pp.png', width=15, height=8, units='in', res=200, pointsize=4)
# tracts_pass_prop_avg %>% filter(type!='SIMPLAI') %>%
#   ggplot(aes(pp, avg_freq)) + 
#   geom_point(aes(colour= type), alpha=0.7, size=3)+
#   geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, colour= type), width=.3,
#                 position=position_dodge(0.05)) +
#   scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
#                                 "FLARE" = "firebrick3", 
#                                 "MOSAIC" = "darkolivegreen4"))+  
#   labs(x = 'Posterior probability cutoff', y = '% of sequence retained', colour='Method') +
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
#             labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source_pop_size)))
# dev.off()


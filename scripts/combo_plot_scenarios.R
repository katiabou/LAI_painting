#import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(viridis)
library(readr)
library(ggpubr)

sequence_length <- as.numeric(snakemake@params[["seq_length"]]) # megabases
#sequence_length <- 300000000
source_time <- as.numeric(snakemake@params[["source_time"]])
source_pop_size <- as.numeric(snakemake@params[["source_pop_size"]])
gf <- as.numeric(snakemake@params[["gf_rate"]])
gen_adm <- as.numeric(snakemake@params[["gen_adm"]])
seed <- as.numeric(snakemake@params[["seed"]])

########## TRACT PLOT ############
# true_tracts <- read.delim('~/Desktop/test/model_gf_0.3_gen_adm_50_tracts.tsv')
# flare <- read.delim("~/Desktop/test/model_gf_0.3_source_20_gen_adm_50_source_time_0_flare_tracts_prob_0.0.tsv") 
# rfmix <- read.delim("~/Desktop/test/model_gf_0.3_source_20_gen_adm_50_source_time_0_rfmix_tracts_prob_0.0.tsv") 
# mosaic <- read.delim("~/Desktop/test/model_gf_0.3_source_20_gen_adm_50_source_time_0_mosaic_tracts_prob_0.0.tsv") 
# simplai <- read.delim("~/Desktop/test/model_gf_0.3_source_20_gen_adm_50_source_time_0_simplai_tracts_prob_0.0.tsv") 

# Tracts input files
true_tracts <- read.delim(snakemake@input[["tracts1"]])
flare <- read.delim(snakemake@input[["flare_tracts1"]]) 
rfmix <- read.delim(snakemake@input[["rfmix_tracts1"]]) 
mosaic <- read.delim(snakemake@input[["mosaic_tracts1"]]) 
simplai <- read.delim(snakemake@input[["simplai_tracts1"]]) 

# Posterior probability file (one seed)
pp <- read.delim(snakemake@input[["all_adm_prop"]])

# Posterior probability cutoffs (one seed)
pp_outputs <- unlist(
    list.files(path = file.path("output/seed_101/compare_software/files_pp_cutoff/tracts/"), pattern = paste0("model_gf_", gf, "_source_", source_pop_size, "_gen_adm_", gen_adm, "_source_time_", source_time, "_.*_tracts_probcutoff_.*.tsv"), full.names = TRUE)
  )

# nMCC with posterior probability cutoffs (one seed)
nmcc_outputs <- unlist(
  list.files(path = file.path("output/seed_101/compare_software/files_pp_cutoff/overlap/"), pattern = paste0("model_gf_", gf, "_source_", source_pop_size, "_gen_adm_", gen_adm, "_source_time_", source_time, "_.*_tracts_overlap_stats_probcutoff_.*.tsv"), full.names = TRUE)
)

# Dates input files
target <- read.table(snakemake@input[["dates_mix_expfit"]], quote="\"")
mean_date <- read.table(snakemake@input[["dates_mean"]], quote="\"", comment.char="")
nrmsd <- read.table(snakemake@input[["dates_nrmsd"]], quote="\"", comment.char="")



#### Start with tracts ####

# Take only one sample from true tracts
s1 <- true_tracts %>% filter(name == 'pop_mix_1')
      
# get first haplotype
s1_hap1 <- s1 %>%
  group_by(name) %>%
  filter(node_id==min(node_id))  %>%
  select('name','source_pop','left','right') %>%
  mutate(type='truth')

# get second haplotype
s1_hap2 <- s1 %>%
  group_by(name) %>%
  filter(node_id==max(node_id)) %>% 
  mutate(left_2=left+sequence_length+1) %>%
  mutate(right_2=right+sequence_length+1) %>%
  select('name','source_pop','left_2','right_2') %>%
  rename('left' = "left_2", 'right' = "right_2") %>%
  mutate(type='truth')

true_all <- rbind(s1_hap1, s1_hap2)


all_methods <- rbind(flare, rfmix, mosaic, simplai) %>%
               separate(range, into = c("left", "right"), sep = "-") %>%
               mutate(across(all_of(c('left','right')), as.numeric)) %>%
               select(name, source_pop, left, right, type)


all_methods_true <- rbind(all_methods, true_all)

# make method uppercase
all_methods_true$type <- toupper(all_methods_true$type)

# values = c("pop_b" = "chartreuse4", 
#            "pop_c" = "cyan3",
#            "pop_no" = "grey70")
values = c("pop_b" = "pink4", 
           "pop_c" = "darkslategray3",
           "pop_no" = "grey70")
           
p1 <- ggplot()+
  geom_segment(data=all_methods_true %>% filter(name=='pop_mix_1') %>% mutate(type=factor(type, levels=c("SIMPLAI","RFMIX","MOSAIC","FLARE", "TRUTH"))), aes(y = type, yend = type, x = left/1000000, xend = right/1000000, colour=source_pop), linewidth = 9) +
  theme_bw()+
  theme(panel.spacing.y=unit(0.1, "lines"),
        axis.title.x = element_text(margin=margin(t=5)),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(size = 10, angle = 0),
        legend.position = 'none',
        # strip.background =element_rect(fill="gray28"),
        # strip.text = element_text(colour = 'white')
  ) +
  scale_color_manual(name = "Ancestry", values = values) +
  labs(x="Genome position (Mb)", y="Sample")


########## TRACT LENGTH PLOT ############

inferred_tracts <- rbind(flare, rfmix, mosaic, simplai)

# fixing ranges for true tracts for both haplotypes and all samples
true_tracts_hap1_range <- true_tracts %>%
  group_by(gen_adm, name) %>%
  #group_by(name) %>%
  filter(node_id==min(node_id)) %>%
  mutate(range=paste(left,'-',right, sep=""), type="TRUTH", gf_rate=gf, source_time=source_time) %>%
  filter(left != 0 | right != 0) %>%
  select(name, source_pop, range, type, gf_rate, gen_adm, source_time, seed)

true_tracts_hap2_range <- true_tracts %>%
  group_by(gen_adm, name) %>%
  #group_by(name) %>%
  filter(node_id==max(node_id)) %>% 
  mutate(left_2=left+sequence_length+1) %>%
  mutate(right_2=right+sequence_length+1) %>%
  mutate(range=paste(left_2,'-',right_2, sep=""), type="TRUTH", gf_rate=gf, source_time=source_time) %>%  #using left_2 and right_2 since it's probably easier to consider one continuous chromosome
  filter(left != 0 | right != 0) %>%
  select(name, source_pop, range, type, gf_rate, gen_adm, source_time, seed)


# merge both haplotypes
true_tracts_range <- rbind(true_tracts_hap1_range, true_tracts_hap2_range)


# Function to filter for source pop sizes I want from inferred and add equivalent column to true, to be able to merge together:
get_source <- function(df_inferred, df_true, x){
  sub_inferred <- df_inferred %>%
    filter(source_pop_size==x)
  
  sub_true <- df_true %>%
    mutate(source_pop_size=x)
  
  sub_all <- rbind(sub_inferred, sub_true)
}

source_sub <- get_source(inferred_tracts, true_tracts_range, source_pop_size)


all <- source_sub %>% 
  separate(range,into=c("start", "end"), sep="-", convert = TRUE) %>%
  mutate(length=end-start) %>% 
  filter(source_pop!="pop_no")


# make method uppercase
all$type <- toupper(all$type)

options(scipen = 999)

values <- c(
  FLARE = "firebrick3", 
  RFMIX = "darkgoldenrod2",
  SIMPLAI = "#74A089",
  MOSAIC = "darkolivegreen4", 
  "TRUTH" = "black")

line_type <- c(
  FLARE = "longdash",
  RFMIX = "longdash",
  SIMPLAI = "longdash",
  MOSAIC = "longdash",
  "TRUTH" = "solid")

# Some NaNs are created using the log10 scale, these are tracts with 0 or -1 length (created when added the ends of the chromosomes to fill up the ancestries across the chromosome)
# They are not an issue for the plot

all$type <- as.factor(all$type)

p2 <- ggplot(all %>% filter(source_pop=='pop_c'), aes(x=(length/1000000), colour=type, linetype=type))+
  scale_linetype_manual(values=line_type, name="Method", drop=FALSE) +
  scale_color_manual(values=values, name="Method", drop=FALSE)+
  geom_density(linewidth = 0.8, show.legend = TRUE) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 10),  
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=11), 
        axis.title.y=element_text(size=11), 
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        plot.title = element_text(hjust = 0.5, size=14),
        panel.grid.minor = element_blank(),
        legend.key.width= unit(0.7, 'cm'),
        legend.position = c(0.15, 0.78)) +
  labs(x="Tract length population C (Mb)", y="Tract density")+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(limits = c(0,3)) 


########## DENSITY PLOT ############
#pp <- read.delim("~/Desktop/test/model_gf_0.3_source_20_gen_adm_50_source_time_0-all_simulated_post_prob.txt")

pp$method <- toupper(pp$method)

p3 <- ggplot(pp, aes(x=prob_anc2))+
  geom_density(aes(fill=method),alpha=0.4, color=NA)+
  scale_fill_manual(values = c("RFMIX" = "darkgoldenrod2",
                               "FLARE" = "firebrick3", 
                               "MOSAIC" = "darkolivegreen4"))+  
  labs(x = 'Posterior probability population C', y = 'Density', fill='Method') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 10),  
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=11), 
        axis.title.y=element_text(size=11), 
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        plot.title = element_text(hjust = 0.5, size=14),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white'),
        legend.key.width= unit(0.7, 'cm'),
        legend.position = c(0.25, 0.84))+
  guides(fill = guide_legend(override.aes = list(alpha = 1)))


  



########## PERCENTAGE SEQUENCE PP PLOT ############
  
#These are the tracts with pp cutoffs 
# pp_outputs <- unlist(
#     list.files(path = file.path("~/Desktop/test/"), pattern = paste0("model_gf_", gf, "_source_20_gen_adm_50_source_time_0_.*_tracts_probcutoff_.*.tsv"), full.names = TRUE)
#   )

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
  

# Add nMCC estimates for different pp cutoffs
# nmcc_outputs <- unlist(
#   list.files(path = file.path("~/Desktop/test/"), pattern = paste0("model_gf_", gf, "_source_20_gen_adm_50_source_time_0_.*_tracts_overlap_stats_probcutoff_.*.tsv"), full.names = TRUE)
# )

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
common_col_names <- intersect(names(all_pp2_popc), names(tracts_pass_prop))
merged_df <- merge(all_pp2_popc, tracts_pass_prop, by=common_col_names, all.x=TRUE)

# replace NA in prop_length column with 0
merged_df$prop_length <- ifelse(is.na(merged_df$prop_length), 0, merged_df$prop_length)

# AddnMCC label
merged_df$nMCC_cutoff <- ifelse(merged_df$nMCC>0.8, '>0.8', '<0.8')
merged_df <- as.data.table(merged_df)

# create order
my_order <- c(">0.8", "<0.8")

# reversed to get correct ordering
merged_df[, order := match(nMCC_cutoff, rev(my_order))]

# sort the order
setorder(merged_df, prop_length, order)

# plot using nMCC cutoff:
p4 <- merged_df %>% filter(type!='SIMPLAI') %>%
  ggplot(aes(pp, prop_length)) + 
  geom_point(aes(colour= type, shape=nMCC_cutoff), alpha=0.4, size=3)+
  scale_shape_manual(values=c(4, 16))+
  scale_colour_manual(values = c("RFMIX" = "darkgoldenrod2",
                                 "FLARE" = "firebrick3", 
                                 "MOSAIC" = "darkolivegreen4"))+  
  labs(x = 'Posterior probability cutoff', y = '% of sequence retained', colour='Method', shape='nMCC') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 10),  
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=11), 
        axis.title.y=element_text(size=11), 
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        plot.title = element_text(hjust = 0.5, size=14),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white'),
        legend.key.width= unit(0.7, 'cm'),
        legend.position = c(0.15, 0.3)) +
  guides(color = guide_legend(override.aes = list(alpha = 1)), shape = guide_legend(override.aes = list(alpha = 1)))



########## DATES PLOT ############
# target <- read.table("~/Desktop/test/pop_mix:expfit.out", quote="\"") #this is simulations
# mean_date <- read.table("~/Desktop/test/model_gf_0.3_source_20_gen_adm_50_source_time_0_dates_adm_time.txt", quote="\"", comment.char="")
# nrmsd <- read.table("~/Desktop/test/model_gf_0.3_source_20_gen_adm_50_source_time_0-nrmsd_dates.txt", quote="\"", comment.char="")

#get mean date estimate for that run:
Mean_date = paste('Inferred admixture time=', mean_date$V7, sep = "")
NRMSD = paste('NRMSD=', nrmsd$V6, sep = "")

#input run details
#gf_rate <- paste('gf_rate=', snakemake@params[1], sep = "")
#gen_adm <- paste('gen_adm=', snakemake@params[2], sep = "")
#source <- paste('source_size=', snakemake@params[3], sep = "")
#source_time <- paste('source_time=', snakemake@params[4], sep = "")
#seed <- paste('seed=', snakemake@params[5], sep = "")

#title
#title_plot <- paste(gf_rate,gen_adm,source,source_time,seed,Mean_date,NRMSD,sep = " | ")
title_plot <- paste(Mean_date,NRMSD,sep='\n')

#export plot
#png(snakemake@output[[1]], width=8, height=6, units='in', res=200, pointsize=4)
#png('~/Downloads/test.png', width=8, height=6, units='in', res=200, pointsize=4)
p5 <- ggplot(target) +
  geom_point(aes(V1, V2), shape = 3, colour='purple') +
  geom_line(aes(V1, V3), colour='green4') +
  annotate(
    geom = "text", x = 4.8, y = 0.017, 
    label = title_plot, hjust = 0, vjust = 1, size = 4
  ) +
  #labs(x='Genetic Distance (cM)', y='Weighted Covariance', title=title_plot)+
  labs(x='Genetic distance (cM)', y='Weighted covariance')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 10),  
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=11), 
        axis.title.y=element_text(size=11), 
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        plot.title = element_text(hjust = 0.5, size=14),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white'),
        legend.key.width= unit(0.7, 'cm'),
        legend.position = c(0.05, 0.9)) +
  xlim(0,20)
#scale_x_continuous(breaks=seq(0, 20, 2), limits = c(0.5,20))+
#scale_y_continuous(breaks=seq(-0.05, 0.25, 0.05), limits = c(-0.05, 0.25))
#dev.off()


# Put them all together
# ggarrange(p1,
#           ggarrange(p2,p3,p4,p5, labels = c("B","C","D","E")),
#           nrow=2,
#           labels="A")

#ggsave('~/Desktop/test/example_fig.png', width=10, height=15, units='in', dpi=300)


library("cowplot")
ggdraw() +
  draw_plot(p1, x = 0, y = 0.76, width = 1, height = 0.2) +
  draw_plot(p2, x = 0, y = 0.38, width = 0.5, height = 0.38) +
  draw_plot(p3, x = 0.5, y = 0.38, width = 0.5, height = 0.38) +
  draw_plot(p4, x = 0, y = 0, width = 0.5, height = 0.38) +
  draw_plot(p5, x = 0.5, y = 0, width = 0.5, height = 0.38) +
  draw_plot_label(label = c("i", "ii", "iii", "iv", "v"), size = 15,
                  x = c(0.02, 0, 0.50, 0, 0.50), y = c(1, 0.77, 0.77, 0.39, 0.39))

#ggsave('~/Desktop/test/example_fig.png', width=10, height=11, units='in', dpi=300)
ggsave(snakemake@output[[1]], width=10, height=11, units='in', dpi=300)


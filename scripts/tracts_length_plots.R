library(GenomicRanges)
library(magrittr)
library(dplyr)
library(tidyr)
library(mltools)
library(ggplot2)
library(readr)
library("cowplot")

# import parameters
sequence_length <- as.numeric(snakemake@params[["seq_length"]])
#sequence_length <- as.numeric(300000000)
gf=as.numeric(snakemake@params[["gf_rate"]])
#gf=as.numeric(0.3)


# THESE ARE THE SOURCE SIZES I WANT TO PLOT. I'VE ONLY CHOSEN SOURCE TIME AT 0. I select them in the unlist function below.
#source <- c(4,10,20,100) #Do not need to show here
#prob1 <- 0.0
prob1 <- as.numeric(snakemake@params[["prob_cutoff"]])
prob <- sprintf("%.1f", prob1) #this is used to input the prob as 0.0 and not 0, to perfectly match the file names
source_time <- as.numeric(snakemake@params[["source_time"]])
#source_time <- 0

# Set path pattern
seed_dirs <- dir(path = "output", pattern = "seed_.*", full.names = TRUE)


# These are the inferred tracts for all methods, but only taking source_time equal to 0 (I'll subset for the source numbers I want later)
tract_inferred <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "compare_software/files/tracts"), pattern = paste0("model_gf_", gf, "_source_.*_gen_adm_.*_source_time_", source_time, "_.*_tracts_prob_", prob, ".tsv"), full.names = TRUE)
}))


# Extract pattern per method
flare <- grep("flare_tracts", tract_inferred, value = TRUE)
mosaic <- grep("mosaic_tracts", tract_inferred, value = TRUE)
rfmix <- grep("rfmix_tracts", tract_inferred, value = TRUE)
simplai <- grep("simplai_tracts", tract_inferred, value = TRUE)

# These are the true tracts
tract_true <- unlist(lapply(seed_dirs, function(seed_dir) {
    list.files(path = file.path(seed_dir, "slendr/sim_data"), pattern = paste0("model_gf_", gf, "_gen_adm_.*_tracts.tsv"), full.names = TRUE)
}))



##### MAKE IMPORTING FILES FASTER
library(future.apply)

plan(multicore, workers = availableCores())

inferred_tracts <- future_lapply(c(flare, mosaic, rfmix, simplai), function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)
true_tracts <- future_lapply(tract_true, function(path) read_tsv(path, col_types = cols())) %>% do.call(rbind, .)

###### Looking only at lengths of inferred segments and lengths of true segments for pop_c and pop_b separately, for all samples

# fixing ranges for true tracts for both haplotypes and all samples
true_tracts_hap1_range <- true_tracts %>%
  group_by(gen_adm, name) %>%
  #group_by(name) %>%
  filter(node_id==min(node_id)) %>%
  mutate(range=paste(left,'-',right, sep=""), type="true", gf_rate=gf, source_time=source_time) %>%
  filter(left != 0 | right != 0) %>%
  select(name, source_pop, range, type, gf_rate, gen_adm, source_time, seed)

true_tracts_hap2_range <- true_tracts %>%
  group_by(gen_adm, name) %>%
  #group_by(name) %>%
  filter(node_id==max(node_id)) %>% 
  mutate(left_2=left+sequence_length+1) %>%
  mutate(right_2=right+sequence_length+1) %>%
  mutate(range=paste(left_2,'-',right_2, sep=""), type="true", gf_rate=gf, source_time=source_time) %>%  #using left_2 and right_2 since it's probably easier to consider one continuous chromosome
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

source_4 <- get_source(inferred_tracts, true_tracts_range, 4)
source_10 <- get_source(inferred_tracts, true_tracts_range, 10)
source_20 <- get_source(inferred_tracts, true_tracts_range, 20)
source_100 <- get_source(inferred_tracts, true_tracts_range, 100)

# Merge all into 1 df:
all <- rbind(source_4, source_10, source_20, source_100)

all <- all %>% 
        separate(range,into=c("start", "end"), sep="-", convert = TRUE) %>%
        mutate(length=end-start) %>% 
        filter(source_pop!="pop_no")

# Make method uppercase
all$type <- toupper(all$type)


options(scipen = 999)

values <- c(
  FLARE = "firebrick3", 
  RFMIX = "darkgoldenrod2",
  SIMPLAI = "#74A089",
  MOSAIC = "darkolivegreen4", 
  "TRUE" = "black")

line_type <- c(
  FLARE = "longdash",
  RFMIX = "longdash",
  SIMPLAI = "longdash",
  MOSAIC = "longdash",
  "TRUE" = "solid")

# Some NaNs are created using the log10 scale, these are tracts with 0 or -1 length (created when added the ends of the chromosomes to fill up the ancestries across the chromosome)
# They are not an issue for the plot

p1 <- ggplot(all %>% filter(source_pop=='pop_b'), aes(x=(length/1000000), colour=type, linetype=type))+
  scale_linetype_manual(values = line_type, name="Method") +
  scale_color_manual(values=values, name="Method")+  geom_density() +
  geom_density(linewidth = 0.6) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 7),  
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
        legend.position="right")+
  labs(x="Tract length population B (Mb)", y="Tract density")+
  facet_grid(source_pop_size~gen_adm,
            labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source_pop_size)))+
  scale_x_continuous(trans='log10', breaks=c(0.0001, 0.01, 1, 100), labels=c("0.0001", "0.01", "1", "100"))

#ggsave("plot1.png", device='png', width=15, height = 8)
#ggsave(snakemake@output[[1]], device='png', width=15, height = 8)



p2 <- ggplot(all %>% filter(source_pop=='pop_c'), aes(x=(length/1000000), colour=type, linetype=type))+
  scale_linetype_manual(values = line_type, name="Method") +
  scale_color_manual(values=values, name="Method")+  geom_density() +
  geom_density(linewidth = 0.6) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size = 7),  
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
        legend.position="right")+
  labs(x="Tract length population C (Mb)", y="Tract density")+
  facet_grid(source_pop_size~gen_adm,
            labeller=label_bquote(cols=t[admix]==.(gen_adm), rows=n==.(source_pop_size)))+
  scale_x_continuous(trans='log10', breaks=c(0.00001, 0.001, 0.1, 10), labels=c("0.00001", "0.001", "0.1", "10")) +
  ylim(0,3) #to limit the massive peak introduced by MOSAIC pushing everything else down

#ggsave("plot2.png", device='png', width=15, height = 8)
#ggsave(snakemake@output[[2]], device='png', width=15, height = 8)

ggdraw() +
  draw_plot(p1, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(p2, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5))

ggsave(snakemake@output[[1]], device='png', width=16, height = 20)

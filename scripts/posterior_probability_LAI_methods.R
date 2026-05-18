#import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

#import run info
gf <- as.numeric(snakemake@params[[1]])
source <- as.numeric(snakemake@params[[2]])
gen_adm <- as.numeric(snakemake@params[[3]])
source_time <- as.numeric(snakemake@params[[4]])
seed <- as.numeric(snakemake@params[[5]])

#import flare source and target names (first column sample name, second is pop name, does not seem like it's used outside of the uniq_ancestry part)
#sources <- read.table("~/Downloads/source_name_pop_model_gf_0.3_source_50_gen_adm_50_source_time_0.txt", quote="\"") 
sources <- read.table(snakemake@input[[1]], quote="\"") 

#one column with sample name of target
#targets <- read.table("~/Downloads/target_name_model_gf_0.3_source_50_gen_adm_50_source_time_0.txt", quote="\"")
targets <- read.table(snakemake@input[[2]], quote="\"")

#extract unique ancestries
uniq_ancestry <- unique(sources$V2)
Anc0 <- uniq_ancestry[1]
Anc1 <- uniq_ancestry[2]


#################
#               #
#     FLARE     #
#               #
#################

options(scipen=999)

#import flare output (contains all samples)
#flare <- read.delim("~/Downloads/model_gf_0.3_source_50_gen_adm_50_source_time_0.txt", header=FALSE)
flare <- read.delim(snakemake@input[[3]], header=FALSE)

#replace chr in chrom column
flare$V1 <- str_replace(flare$V1, "chr", "")
flare$V1 <- as.numeric(flare$V1)

#remove unwanted column
flare_1 <- flare %>% select(-3)

#create dummy list
datalist = list()

#for each sample, split the column with the ancestry and ancestry probabilities
for (i in 3:ncol(flare_1)){
  flare_n <- flare_1 %>% select(1,2,i)
  
  #seperate and create probability columns for each haplotype and each ancestry
  flare_n2 <- separate(data=flare_n, col = 3, into = c("GT", "AN1", "AN2", "ANP1", "ANP2"), sep = ":") %>%
    separate(ANP1, c("PROB_HAPL1_ANC1", "PROB_HAPL1_ANC2"), sep=',') %>%
    separate(ANP2, c("PROB_HAPL2_ANC1", "PROB_HAPL2_ANC2"), sep=',') 
  
  sampleno <- i-2
  flare_n2$sample <- paste("pop_mix_", sampleno, sep = "")
  flare_n2$AN1[flare_n2$AN1=='0'] <- Anc0
  flare_n2$AN1[flare_n2$AN1=='1'] <- Anc1
  flare_n2$AN2[flare_n2$AN2=='0'] <- Anc0
  flare_n2$AN2[flare_n2$AN2=='1'] <- Anc1
  datalist[[i]] <- flare_n2
} 

big_data = do.call(rbind, datalist) #combine all together


#get PP from each ancestry and haplotype
big_data$PROB_HAPL1_ANC1 <- as.numeric(big_data$PROB_HAPL1_ANC1)
big_data$PROB_HAPL1_ANC2 <- as.numeric(big_data$PROB_HAPL1_ANC2)
big_data$PROB_HAPL2_ANC1 <- as.numeric(big_data$PROB_HAPL2_ANC1)
big_data$PROB_HAPL2_ANC2 <- as.numeric(big_data$PROB_HAPL2_ANC2)

#combine the pp of each ancestry from both haplotypes:
pp_popb_flare <- data.frame(prob_anc1 = c(big_data[,"PROB_HAPL1_ANC1"], big_data[,"PROB_HAPL2_ANC1"]))
pp_popc_flare <- data.frame(prob_anc2 = c(big_data[,"PROB_HAPL1_ANC2"], big_data[,"PROB_HAPL2_ANC2"]))

pp_flare <- cbind(pp_popb_flare, pp_popc_flare)
pp_flare$seed <- seed
pp_flare$gf <- gf
pp_flare$source <- source
pp_flare$gen_adm <- gen_adm
pp_flare$source_time <- source_time
pp_flare$method <- 'FLARE'

# plot_pp_popb_flare <- ggplot(pp_popb_flare) + 
#   geom_histogram(aes(x = prob_anc1), fill = "navy", binwidth = 0.1) +
#   theme_bw()+
#   theme(panel.spacing.y=unit(0.1, "lines"),
#         axis.title.x = element_text(margin=margin(t=5)),
#         axis.title.y = element_text(margin=margin(r=5)),
#         axis.text.x = element_text(color = "black"),
#         axis.ticks.x = element_line(linewidth = 0.3),
#         plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.text.y.left = element_text(size = 10, angle = 0),
#         legend.position = 'none',
#         strip.background =element_rect(fill="gray28"),
#         strip.text = element_text(colour = 'white')) +
#   xlab("Posterior probability for population B ancestry")



#################
#               #
#     MOSAIC    #
#               #
#################

library(MOSAIC)

la_results = snakemake@input[[4]]
#la_results = '~/Downloads/localanc_pop_mix.RData'

model_results = snakemake@input[[5]]
#model_results = '~/Downloads/pop_mix.RData'

mosaic_input_dir = snakemake@params[[6]]
#mosaic_input_dir = '~/Downloads/'

snp=read.table(snakemake@input[[6]], quote="\"", comment.char="")
#snp=read.table('~/Downloads/snpfile.1', quote="\"", comment.char="")

#load files
load(model_results)
load(la_results)

# localanc gives the local ancestry at each grid point
# get local ancestry probabilities at each SNP
local_pos=grid_to_pos(localanc, mosaic_input_dir, g.loc, chrnos)  #look here for code https://github.com/rwaples/lai-sim/blob/main/workflow/scripts/export_mosaic_results.R

local_pos_df <- reshape2::melt(local_pos[[1]])
colnames(local_pos_df) <- c("ancestry", "haplotype", "number", "prob")

#create dummy list
datalist_mosaic_hapl1 = list()
datalist_mosaic_hapl2 = list()

for(i in seq(1, max(local_pos_df$haplotype), 2)){
  
  mosaic_sample_hap1 <- local_pos_df %>%
    filter(haplotype %in% c(i)) %>% 
    pivot_wider(names_from = ancestry, values_from = prob, names_prefix="prob_anc") %>%
    mutate(AN1 = ifelse(prob_anc1>prob_anc2, Anc0, Anc1)) %>%
    #mutate(AN1 = ifelse(prob_anc1<threshold & prob_anc2<threshold, 'pop_no', AN1)) %>%
    mutate(sample = paste("pop_mix_", (i+1)/2, sep = ""))
  mosaic_sample_hap1_snp <- cbind(mosaic_sample_hap1, snp[,4, drop=FALSE])
  
  datalist_mosaic_hapl1[[i]] <- mosaic_sample_hap1_snp
  
  mosaic_sample_hap2 <- local_pos_df %>%
    filter(haplotype %in% c(i+1)) %>% 
    pivot_wider(names_from = ancestry, values_from = prob, names_prefix="prob_anc") %>%
    mutate(AN2 = ifelse(prob_anc1>prob_anc2, Anc0, Anc1)) %>%
    #mutate(AN2 = ifelse(prob_anc1<threshold & prob_anc2<threshold, 'pop_no', AN2)) %>%
    mutate(sample = paste("pop_mix_", (i+1)/2, sep = ""))
  mosaic_sample_hap2_snp <- cbind(mosaic_sample_hap2, snp[,4, drop=FALSE])
  
  datalist_mosaic_hapl2[[i]] <- mosaic_sample_hap2_snp
}

big_data_mosaic_hap1 = do.call(rbind, datalist_mosaic_hapl1) #combine all together hapl1
big_data_mosaic_hap2 = do.call(rbind, datalist_mosaic_hapl2) #combine all together hapl2


#combine the pp of each ancestry from both haplotypes:
pp_popb_mosaic <- data.frame(prob_anc1 = c(big_data_mosaic_hap1[,"prob_anc1"], big_data_mosaic_hap2[,"prob_anc1"]))
pp_popc_mosaic <- data.frame(prob_anc2 = c(big_data_mosaic_hap1[,"prob_anc2"], big_data_mosaic_hap2[,"prob_anc2"]))

pp_mosaic <- cbind(pp_popb_mosaic, pp_popc_mosaic)
pp_mosaic$seed <- seed
pp_mosaic$gf <- gf
pp_mosaic$source <- source
pp_mosaic$gen_adm <- gen_adm
pp_mosaic$source_time <- source_time
pp_mosaic$method <- 'MOSAIC'



#################
#               #
#     RFMIX     #
#               #
#################

#import rfmix output (contains all samples)
fb <- read.delim(snakemake@input[[7]], comment.char="#")
#fb <- read.delim("~/Downloads/model_gf_0.3_source_50_gen_adm_50_source_time_0.fb.tsv", comment.char="#")

#replace chr in chrom column
fb$chromosome <- str_replace(fb$chromosome, "chr", "")
fb$chromosome <- as.numeric(fb$chromosome)

#remove unwanted column
fb_1 <- fb %>% select(c(-1:-4))

#take 1 df per sample
a <- split.default(fb_1, rep(1:nrow(targets), each = 4))

#create dummy list
datalist = list()

for (i in 1:length(a)) {
  fb_n <- a[[i]]
  sampleno <- i
  colnames(fb_n) <- c("PROB_HAPL1_ANC1", "PROB_HAPL1_ANC2", "PROB_HAPL2_ANC1", "PROB_HAPL2_ANC2")
  fb_n$sample <- paste("pop_mix_", sampleno, sep = "")
  fb_n$chr <- fb$chromosome
  fb_n$pos <- fb$physical_position
  datalist[[i]] <- fb_n
}

big_data_n = do.call(rbind, datalist) #combine all together


#combine the pp of each ancestry from both haplotypes:
pp_popb_rfmix <- data.frame(prob_anc1 = c(big_data_n[,"PROB_HAPL1_ANC1"], big_data_n[,"PROB_HAPL2_ANC1"]))
pp_popc_rfmix <- data.frame(prob_anc2 = c(big_data_n[,"PROB_HAPL1_ANC2"], big_data_n[,"PROB_HAPL2_ANC2"]))

pp_rfmix <- cbind(pp_popb_rfmix, pp_popc_rfmix)
pp_rfmix$seed <- seed
pp_rfmix$gf <- gf
pp_rfmix$source <- source
pp_rfmix$gen_adm <- gen_adm
pp_rfmix$source_time <- source_time
pp_rfmix$method <- 'rfmix'


#merge pp files for this run:
all_pp <- rbind(pp_flare, pp_mosaic, pp_rfmix)

#export file
#write.table(all_pp, file='~/Downloads/all_pp.tsv', quote=FALSE, sep='\t', row.names=FALSE)
write.table(all_pp, file=snakemake@output[[1]], quote=FALSE, sep='\t', row.names=FALSE)




#import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(viridis)

# import parameters
sequence_length <- as.numeric(snakemake@params[["seq_length"]]) # megabases
threshold <- as.numeric(snakemake@params[["prob_cutoff"]])
GF_RATE=as.numeric(snakemake@params[["gf_rate"]])
#GF_RATE=as.numeric(0.3)
SOURCE_POP_SIZE=as.numeric(snakemake@params[["source_pop_size"]])
#SOURCE_POP_SIZE=as.numeric(50)
GEN_ADM=as.numeric(snakemake@params[["gen_adm"]])
#GEN_ADM=as.numeric(50)
SOURCE_TIME=snakemake@params[["source_time"]]
#SOURCE_TIME=1
SEED=as.numeric(snakemake@params[["seed"]])
#SEED=as.numeric(101)

# import true tracts
#all_tracts <- read.delim("~/Downloads/model_gf_0.3_gen_adm_50_tracts.tsv")
all_tracts <- read.delim(snakemake@input[[1]])


# import sources
#sources <- read.table("~/Downloads/sources_model_gf_0.3_source_50_gen_adm_50_source_time_1.txt", quote="\"") 
sources <- read.table(snakemake@input[[2]], quote="\"") 

# import targets
#targets <- read.table("~/Downloads/targets_model_gf_0.3_source_50_gen_adm_50_source_time_1.txt", quote="\"")
targets <- read.table(snakemake@input[[3]], quote="\"")

# extract unique ancestries
uniq_ancestry <- unique(sources$V2)
Anc0 <- uniq_ancestry[1]
Anc1 <- uniq_ancestry[2]


#################
#               #
#     FLARE     #
#               #
#################

options(scipen=999)

# import flare output (contains all samples)
#flare <- read.delim("~/Downloads/model_gf_0.3_source_50_gen_adm_50_source_time_1.txt", header=FALSE)
flare <- read.delim(snakemake@input[[4]], header=FALSE)

# replace chr in chrom column
flare$V1 <- str_replace(flare$V1, "chr", "")
flare$V1 <- as.numeric(flare$V1)

# remove unwanted column
flare_1 <- flare %>% select(-3)

# create dummy list
datalist = list()

# for each sample, split the column with the ancestry and ancestry probabilities
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


# change ancestry to pop_no with low probabilities
big_data$AN1 <- ifelse(big_data$PROB_HAPL1_ANC1<threshold & big_data$PROB_HAPL1_ANC2<threshold, 'pop_no', big_data$AN1)
big_data$AN2 <- ifelse(big_data$PROB_HAPL2_ANC1<threshold & big_data$PROB_HAPL2_ANC2<threshold, 'pop_no', big_data$AN2)


# create duplicated first row and replace with site=1
flare_hapl_dupl <- big_data %>%
  group_by(sample) %>%
  dplyr::slice(1L, row_number()) %>% #replicate the first row
  mutate(V2 = ifelse(row_number() == 1, 0, V2)) #replace the first snp to be 0


#### make ranges

# first haplotype
b1 <- flare_hapl_dupl %>%
  group_by(sample, rleid=with(rle(AN1), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1) 

flare_1_mod <- b1 %>%
  group_by(sample) %>%
  mutate(right=lead(V2)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(V2, "-", right, sep="")) %>% #make new range
  rename("left" = "V2")


# second haplotype
b2 <- flare_hapl_dupl %>%
  group_by(sample, rleid=with(rle(AN2), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1) 

flare_2_mod <- b2 %>%
  group_by(sample) %>%
  mutate(right=lead(V2)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(V2, "-", right, sep="")) %>% #make new range
  rename("left" = "V2") %>%
  mutate(left_2=left+sequence_length+1) %>% #create new right and left for plotting
  mutate(right_2=right+sequence_length+1) 

values = c("pop_b" = "chartreuse4", 
           "pop_c" = "cyan3",
           "pop_no" = "grey70")

# ggplot()+
#   geom_segment(data=flare_1_mod, aes(y = 0, yend = 0, x = left/1000000, xend = right/1000000, colour=AN1), linewidth = 15) +
#   geom_segment(data=flare_2_mod, aes(y = 0, yend = 0, x = left_2/1000000, xend = right_2/1000000, colour=AN2), linewidth = 15) +
#   theme_bw()+
#   theme(panel.spacing.y=unit(0.1, "lines"),
#         #panel.grid=element_blank(),
#         axis.title.x = element_text(margin=margin(t=5)),
#         axis.title.y = element_text(margin=margin(r=5)),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(color = "black"),
#         axis.ticks.x = element_line(linewidth = 0.3),
#         plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.text.y.left = element_text(size = 10, angle = 0),
#         legend.position = 'none',
#         strip.background =element_rect(fill="gray28"),
#         strip.text = element_text(colour = 'white')
#         #axis.line.x = element_line(size = 0.3)
#   ) +
#   scale_color_manual(name = "Ancestry", values = values) +
#   ylab("Sample") +
#   xlab("Genome position (Mb)") +
#   facet_grid(sample ~ .)



#################
#               #
#     MOSAIC    #
#               #
#################

library(MOSAIC)

la_results = snakemake@input[[5]]
#la_results = '~/Downloads/localanc_pop_mix.RData'

model_results = snakemake@input[[6]]
#model_results = '~/Downloads/pop_mix.RData'

mosaic_input_dir = snakemake@params[[3]]
#mosaic_input_dir = '~/Downloads/'

snp=read.table(snakemake@input[[7]], quote="\"", comment.char="")
#snp=read.table('~/Downloads/snpfile.1', quote="\"", comment.char="")

# load files
load(model_results)
load(la_results)

# localanc gives the local ancestry at each grid point
# get local ancestry probabilities at each SNP
local_pos=grid_to_pos(localanc, mosaic_input_dir, g.loc, chrnos)  #look here for code https://github.com/rwaples/lai-sim/blob/main/workflow/scripts/export_mosaic_results.R

local_pos_df <- reshape2::melt(local_pos[[1]])
colnames(local_pos_df) <- c("ancestry", "haplotype", "number", "prob")

# create dummy list
datalist_mosaic_hapl1 = list()
datalist_mosaic_hapl2 = list()

for(i in seq(1, max(local_pos_df$haplotype), 2)){
  
  mosaic_sample_hap1 <- local_pos_df %>%
    filter(haplotype %in% c(i)) %>% 
    pivot_wider(names_from = ancestry, values_from = prob, names_prefix="prob_anc") %>%
    mutate(AN1 = ifelse(prob_anc1>prob_anc2, Anc0, Anc1)) %>%
    mutate(AN1 = ifelse(prob_anc1<threshold & prob_anc2<threshold, 'pop_no', AN1)) %>%
    mutate(sample = paste("pop_mix_", (i+1)/2, sep = ""))
  mosaic_sample_hap1_snp <- cbind(mosaic_sample_hap1, snp[,4, drop=FALSE])
  
  datalist_mosaic_hapl1[[i]] <- mosaic_sample_hap1_snp
  
  mosaic_sample_hap2 <- local_pos_df %>%
    filter(haplotype %in% c(i+1)) %>% 
    pivot_wider(names_from = ancestry, values_from = prob, names_prefix="prob_anc") %>%
    mutate(AN2 = ifelse(prob_anc1>prob_anc2, Anc0, Anc1)) %>%
    mutate(AN2 = ifelse(prob_anc1<threshold & prob_anc2<threshold, 'pop_no', AN2)) %>%
    mutate(sample = paste("pop_mix_", (i+1)/2, sep = ""))
  mosaic_sample_hap2_snp <- cbind(mosaic_sample_hap2, snp[,4, drop=FALSE])
  
  datalist_mosaic_hapl2[[i]] <- mosaic_sample_hap2_snp
}

big_data_mosaic_hap1 = do.call(rbind, datalist_mosaic_hapl1) #combine all together hapl1
big_data_mosaic_hap2 = do.call(rbind, datalist_mosaic_hapl2) #combine all together hapl2


# create duplicated first row and replace with site=1
big_data_mosaic_hap1_dupl <- big_data_mosaic_hap1 %>%
  group_by(sample) %>%
  dplyr::slice(1L, row_number()) %>% #replicate the first row
  mutate(V4 = ifelse(row_number() == 1, 0, V4)) #replace the first snp to be 0


big_data_mosaic_hap2_dupl <- big_data_mosaic_hap2 %>%
  group_by(sample) %>%
  dplyr::slice(1L, row_number()) %>% #replicate the first row
  mutate(V4 = ifelse(row_number() == 1, 0, V4)) #replace the first snp to be 0


#### make ranges

# first haplotype
b1_mosaic <- big_data_mosaic_hap1_dupl %>%
  group_by(sample, rleid=with(rle(AN1), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1) 

mosaic_1_mod <- b1_mosaic %>%
  group_by(sample) %>%
  mutate(right=lead(V4)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(V4, "-", right, sep="")) %>% #make new range
  rename("left" = "V4")


# second haplotype
b2_mosaic <- big_data_mosaic_hap2_dupl %>%
  group_by(sample, rleid=with(rle(AN2), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1) 

mosaic_2_mod <- b2_mosaic %>%
  group_by(sample) %>%
  mutate(right=lead(V4)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(V4, "-", right, sep="")) %>% #make new range
  rename("left" = "V4") %>%
  mutate(left_2=left+sequence_length+1) %>% #create new right and left for plotting
  mutate(right_2=right+sequence_length+1) 

# ggplot()+
#   geom_segment(data=mosaic_1_mod, aes(y = 0, yend = 0, x = left/1000000, xend = right/1000000, colour=AN1), linewidth = 15) +
#   geom_segment(data=mosaic_2_mod, aes(y = 0, yend = 0, x = left_2/1000000, xend = right_2/1000000, colour=AN2), linewidth = 15) +
#   theme_bw()+
#   theme(panel.spacing.y=unit(0.1, "lines"),
#         #panel.grid=element_blank(),
#         axis.title.x = element_text(margin=margin(t=5)),
#         axis.title.y = element_text(margin=margin(r=5)),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(color = "black"),
#         axis.ticks.x = element_line(linewidth = 0.3),
#         plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.text.y.left = element_text(size = 10, angle = 0),
#         legend.position = 'none',
#         strip.background =element_rect(fill="gray28"),
#         strip.text = element_text(colour = 'white')
#         #axis.line.x = element_line(size = 0.3)
#   ) +
#   scale_color_manual(name = "Ancestry", values = values) +
#   ylab("Sample") +
#   xlab("Genome position (Mb)") +
#   facet_grid(sample ~ .)


#################
#               #
#     RFMIX     #
#               #
#################

# import rfmix output (contains all samples)
fb <- read.delim(snakemake@input[[8]], comment.char="#")
#fb <- read.delim("~/Downloads/model_gf_0.3_source_50_gen_adm_50_source_time_1.fb.tsv", comment.char="#")

# replace chr in chrom column
fb$chromosome <- str_replace(fb$chromosome, "chr", "")
fb$chromosome <- as.numeric(fb$chromosome)

# remove unwanted column
fb_1 <- fb %>% select(c(-1:-4))

# take 1 df per sample
a <- split.default(fb_1, rep(1:nrow(targets), each = 4))

# create dummy list
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

# assign ancestries based on probabilities, but mask with pop_no anything below our probability threshold
big_data_prob <- big_data_n %>% 
  mutate(AN1 = ifelse(PROB_HAPL1_ANC1>PROB_HAPL1_ANC2, Anc0, Anc1)) %>%
  mutate(AN1 = ifelse(PROB_HAPL1_ANC1<threshold & PROB_HAPL1_ANC2<threshold, 'pop_no', AN1)) %>%
  mutate(AN2 = ifelse(PROB_HAPL2_ANC1>PROB_HAPL2_ANC2, Anc0, Anc1)) %>%
  mutate(AN2 = ifelse(PROB_HAPL2_ANC1<threshold & PROB_HAPL2_ANC2<threshold, 'pop_no', AN2))


# create duplicated first row and replace with site=1
big_data_prob_dupl <- big_data_prob %>%
  group_by(sample) %>%
  dplyr::slice(1L, row_number()) %>% #replicate the first row
  mutate(pos = ifelse(row_number() == 1, 0, pos)) #replace the first snp to be 0


#### make ranges

# first haplotype
b1 <- big_data_prob_dupl %>%
  group_by(sample, rleid=with(rle(AN1), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1) 

rfmix_1_mod <- b1 %>%
  group_by(sample) %>%
  mutate(right=lead(pos)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(pos, "-", right, sep="")) %>% #make new range
  rename("left" = "pos")


# second haplotype
b2 <- big_data_prob_dupl %>%
  group_by(sample, rleid=with(rle(AN2), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1) 

rfmix_2_mod <- b2 %>%
  group_by(sample) %>%
  mutate(right=lead(pos)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(pos, "-", right, sep="")) %>% #make new range
  rename("left" = "pos") %>%
  mutate(left_2=left+sequence_length+1) %>% #create new right and left for plotting
  mutate(right_2=right+sequence_length+1) 


# ggplot()+
#   geom_segment(data=rfmix_1_mod, aes(y = 0, yend = 0, x = left/1000000, xend = right/1000000, colour=AN1), linewidth = 15) +
#   geom_segment(data=rfmix_2_mod, aes(y = 0, yend = 0, x = left_2/1000000, xend = right_2/1000000, colour=AN2), linewidth = 15) +
#   theme_bw()+
#   theme(panel.spacing.y=unit(0.1, "lines"),
#         #panel.grid=element_blank(),
#         axis.title.x = element_text(margin=margin(t=5)),
#         axis.title.y = element_text(margin=margin(r=5)),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(color = "black"),
#         axis.ticks.x = element_line(linewidth = 0.3),
#         plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.text.y.left = element_text(size = 10, angle = 0),
#         legend.position = 'none',
#         strip.background =element_rect(fill="gray28"),
#         strip.text = element_text(colour = 'white')
#         #axis.line.x = element_line(size = 0.3)
#   ) +
#   scale_color_manual(name = "Ancestry", values = values) +
#   ylab("Sample") +
#   xlab("Genome position (Mb)") +
#   facet_grid(sample ~ .)


##################
#                #
#     SIMPLAI    #
#                #
##################

##### input simplai output:
#simplai <- read.delim("~/Downloads/samples_model_gf_0.3_source_50_gen_adm_50_source_time_1_n_1000_m_1000_t_5_fromRec_withSingl.adm")
simplai <- read.delim(snakemake@input[[9]])

# remove first column from df
simplai_win <- simplai$beg
simplai_mod <- simplai[-1]

# getting number of columns in R 
cols <- ncol(simplai_mod) 

# extracting odd rows  
odd_cols <- seq_len(cols) %% 2 

# getting data from odd data frame 
data_mod_1 <- simplai_mod[, odd_cols == 1] 
data_mod_2 <- simplai_mod[, odd_cols == 0] 

# add sample names 
colnames(data_mod_1) <- targets$V1
colnames(data_mod_2) <- targets$V1

# add sites
simplai_hap1 <- cbind(simplai_win, data_mod_1)
simplai_hap2 <- cbind(simplai_win, data_mod_2)

# remove last row if the position is lower than the previous row, this is important since there is a bug in simplai, where sometimes the last bp or window is lower 
# than the previous one, and gives an error when making the ranges and testing for overlaps
if (simplai_hap1$simplai_win[nrow(simplai_hap1)] < simplai_hap1$simplai_win[nrow(simplai_hap1) - 1]) {
  # Remove the last row if the condition is met
  simplai_hap1 <- simplai_hap1[-nrow(simplai_hap1), ]
}

if (simplai_hap2$simplai_win[nrow(simplai_hap2)] < simplai_hap2$simplai_win[nrow(simplai_hap2) - 1]) {
  # Remove the last row if the condition is met
  simplai_hap2 <- simplai_hap2[-nrow(simplai_hap2), ]
}

# melt dataframe
simplai_hap1_tf <- reshape2::melt(simplai_hap1, id='simplai_win')
simplai_hap2_tf <- reshape2::melt(simplai_hap2, id='simplai_win')

# replace ancestries
simplai_hap1_tf$value[simplai_hap1_tf$value=='1'] <- Anc0
simplai_hap1_tf$value[simplai_hap1_tf$value=='2'] <- Anc1
simplai_hap2_tf$value[simplai_hap2_tf$value=='1'] <- Anc0
simplai_hap2_tf$value[simplai_hap2_tf$value=='2'] <- Anc1


# replace column names
colnames(simplai_hap1_tf) <- c('simplai_win', 'sample','AN1')
colnames(simplai_hap2_tf) <- c('simplai_win', 'sample','AN2')

# create ranges

# first haplotype
b1_simplai <- simplai_hap1_tf %>%
  group_by(sample, rleid=with(rle(AN1), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1)

simplai_1_mod <- b1_simplai %>%
  group_by(sample) %>%
  mutate(right=lead(simplai_win)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(simplai_win, "-", right, sep="")) %>% #make new range
  rename("left" = "simplai_win")

# second haplotype
b2_simplai <- simplai_hap2_tf %>%
  group_by(sample, rleid=with(rle(AN2), rep(seq_along(lengths), lengths))) %>% #group by sample, and also make ranges where ancestries change
  dplyr::slice(1) 

simplai_2_mod <- b2_simplai %>%
  group_by(sample) %>%
  mutate(right=lead(simplai_win)-1) %>%
  mutate(right = replace_na(right, sequence_length)) %>% #replace the last cell of each sample with the sequence length since NA was introduced in previous step
  mutate(range = paste(simplai_win, "-", right, sep="")) %>% #make new range
  rename("left" = "simplai_win") %>%
  mutate(left_2=left+sequence_length+1) %>% #create new right and left for plotting
  mutate(right_2=right+sequence_length+1) 


# ggplot()+
#   geom_segment(data=simplai_1_mod, aes(y = 0, yend = 0, x = left/1000000, xend = right/1000000, colour=AN1), linewidth = 15) +
#   geom_segment(data=simplai_2_mod, aes(y = 0, yend = 0, x = left_2/1000000, xend = right_2/1000000, colour=AN2), linewidth = 15) +
#   theme_bw()+
#   theme(panel.spacing.y=unit(0.1, "lines"),
#         #panel.grid=element_blank(),
#         axis.title.x = element_text(margin=margin(t=5)),
#         axis.title.y = element_text(margin=margin(r=5)),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(color = "black"),
#         axis.ticks.x = element_line(linewidth = 0.3),
#         plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.text.y.left = element_text(size = 10, angle = 0),
#         legend.position = 'none',
#         strip.background =element_rect(fill="gray28"),
#         strip.text = element_text(colour = 'white')
#         #axis.line.x = element_line(size = 0.3)
#   ) +
#   scale_color_manual(name = "Ancestry", values = values) +
#   ylab("Sample") +
#   xlab("Genome position (Mb)") +
#   facet_grid(sample ~ .)



#######################
#                     #
#     ALL TOGETHER    #
#                     #
#######################

# true tracts

# get first haplotype
all_tracts_hap1_join <- all_tracts %>%
  group_by(name) %>%
  filter(node_id==min(node_id))  %>%
  select('name','source_pop','left','right') %>%
  mutate(type='true')

# get second haplotype
all_tracts_hap2_join <- all_tracts %>%
  group_by(name) %>%
  filter(node_id==max(node_id)) %>% 
  mutate(left_2=left+sequence_length+1) %>%
  mutate(right_2=right+sequence_length+1) %>%
  select('name','source_pop','left_2','right_2') %>%
  mutate(type='true')


# flare tracts
flare_1_mod_join <- flare_1_mod %>% 
  select('sample','AN1','left','right') %>%
  rename("name" = "sample", "source_pop"="AN1") %>%
  mutate(type='flare')


flare_2_mod_join <- flare_2_mod %>% 
  select('sample','AN2','left_2','right_2') %>%
  rename("name" = "sample", "source_pop"="AN2") %>%
  mutate(type='flare')


# mosaic tracts
mosaic_1_mod_join <- mosaic_1_mod %>% 
  select('sample','AN1','left','right') %>%
  rename("name" = "sample", "source_pop"="AN1") %>%
  mutate(type='mosaic')


mosaic_2_mod_join <- mosaic_2_mod %>% 
  select('sample','AN2','left_2','right_2') %>%
  rename("name" = "sample", "source_pop"="AN2") %>%
  mutate(type='mosaic')

# rfmix tracts
rfmix_1_mod_join <- rfmix_1_mod %>% 
  select('sample','AN1','left','right') %>%
  rename("name" = "sample", "source_pop"="AN1") %>%
  mutate(type='rfmix')

rfmix_2_mod_join <- rfmix_2_mod %>% 
  select('sample','AN2','left_2','right_2') %>%
  rename("name" = "sample", "source_pop"="AN2") %>%
  mutate(type='rfmix')


# simplai tracts
simplai_1_mod_join <- simplai_1_mod %>% 
  select('sample','AN1','left','right') %>%
  rename("name" = "sample", "source_pop"="AN1") %>%
  mutate(type='simplai')


simplai_2_mod_join <- simplai_2_mod %>% 
  select('sample','AN2','left_2','right_2') %>%
  rename("name" = "sample", "source_pop"="AN2") %>%
  mutate(type='simplai')


# the columns I need are the ancestry one, the right and the left (or right_2 left_2 for the 2nd haplotype)
all_join <- rbind(all_tracts_hap1_join, flare_1_mod_join, mosaic_1_mod_join, rfmix_1_mod_join, simplai_1_mod_join)
all_join_2 <- rbind(all_tracts_hap2_join, flare_2_mod_join, mosaic_2_mod_join, rfmix_2_mod_join, simplai_2_mod_join)


values = c("pop_b" = "chartreuse4", 
           "pop_c" = "cyan3",
           "pop_no" = "grey70")


plot_tracts <- ggplot()+
  geom_segment(data=all_join %>% mutate(type=factor(type, levels=c("simplai","rfmix","mosaic","flare","true"))), aes(y = type, yend = type, x = left/1000000, xend = right/1000000, colour=source_pop), linewidth = 9) +
  geom_segment(data=all_join_2 %>% mutate(type=factor(type, levels=c("simplai","rfmix","mosaic","flare","true"))), aes(y = type, yend = type, x = left_2/1000000, xend = right_2/1000000, colour=source_pop), linewidth = 9) +
  theme_bw()+
  theme(panel.spacing.y=unit(0.1, "lines"),
        axis.title.x = element_text(margin=margin(t=5)),
        axis.title.y = element_text(margin=margin(r=5)),
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(size = 10, angle = 0),
        legend.position = 'none',
        strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white')
  ) +
  scale_color_manual(name = "Ancestry", values = values) +
  ylab("Sample") +
  xlab("Genome position (Mb)") +
  facet_grid(name~.)


# export plot
ggsave(snakemake@output[[1]], plot_tracts, width = 7, height = 17, bg='transparent')


##### EXPORT RANGES FOR EACH METHOD #####


# join both haplotypes into one file (the second haplotype is appended to the first coordinate wise)
flare_1_mod_join_export <- flare_1_mod_join %>% 
  mutate(range=paste(left,'-',right, sep="")) %>%
  select(name, source_pop, range, type)

flare_2_mod_join_export <- flare_2_mod_join %>% 
  mutate(range=paste(left_2,'-',right_2, sep="")) %>%
  select(name, source_pop, range, type)

mosaic_1_mod_join_export <- mosaic_1_mod_join %>%
  mutate(range=paste(left,'-',right, sep="")) %>%
  select(name, source_pop, range, type)

mosaic_2_mod_join_export <- mosaic_2_mod_join %>%
  mutate(range=paste(left_2,'-',right_2, sep="")) %>%
  select(name, source_pop, range, type)

rfmix_1_mod_join_export <- rfmix_1_mod_join %>%
  mutate(range=paste(left,'-',right, sep="")) %>%
  select(name, source_pop, range, type)

rfmix_2_mod_join_export <- rfmix_2_mod_join %>%
  mutate(range=paste(left_2,'-',right_2, sep="")) %>%
  select(name, source_pop, range, type)

simplai_1_mod_join_export <- simplai_1_mod_join %>%
  mutate(range=paste(left,'-',right, sep="")) %>%
  select(name, source_pop, range, type)

simplai_2_mod_join_export <- simplai_2_mod_join %>%
  mutate(range=paste(left_2,'-',right_2, sep="")) %>%
  select(name, source_pop, range, type)

# join haplotypes
flare_mod_join_export <- rbind(flare_1_mod_join_export, flare_2_mod_join_export)
mosaic_mod_join_export <- rbind(mosaic_1_mod_join_export, mosaic_2_mod_join_export)
rfmix_mod_join_export <- rbind(rfmix_1_mod_join_export, rfmix_2_mod_join_export)
simplai_mod_join_export <- rbind(simplai_1_mod_join_export, simplai_2_mod_join_export)


# add metadata to haplotype files
flare_mod_join_export1 <- flare_mod_join_export %>%
                   mutate(gf_rate=GF_RATE, source_pop_size=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIME, seed=SEED)

mosaic_mod_join_export1 <- mosaic_mod_join_export %>%
                   mutate(gf_rate=GF_RATE, source_pop_size=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIME, seed=SEED)

rfmix_mod_join_export1 <- rfmix_mod_join_export %>%
                   mutate(gf_rate=GF_RATE, source_pop_size=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIME, seed=SEED)

simplai_mod_join_export1 <- simplai_mod_join_export %>%
                   mutate(gf_rate=GF_RATE, source_pop_size=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIME, seed=SEED)


# export tracts
write.table(flare_mod_join_export1, file=snakemake@output[[2]], quote=FALSE, sep='\t', row.names=FALSE)
write.table(mosaic_mod_join_export1, file=snakemake@output[[3]], quote=FALSE, sep='\t', row.names=FALSE)
write.table(rfmix_mod_join_export1, file=snakemake@output[[4]], quote=FALSE, sep='\t', row.names=FALSE)
write.table(simplai_mod_join_export1, file=snakemake@output[[5]], quote=FALSE, sep='\t', row.names=FALSE)







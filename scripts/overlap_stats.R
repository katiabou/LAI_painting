library(GenomicRanges)
library(magrittr)
library(dplyr)
library(mltools)

sequence_length <- as.numeric(snakemake@params[["seq_length"]]) # megabases

GF_RATE=as.numeric(snakemake@params[["gf_rate"]])
#GF_RATE=as.numeric(0.3)
SOURCE_POP_SIZE=as.numeric(snakemake@params[["source_pop_size"]])
#SOURCE_POP_SIZE=as.numeric(50)
GEN_ADM=as.numeric(snakemake@params[["gen_adm"]])
#GEN_ADM=as.numeric(50)
#SOURCE_TIME=1
SOURCE_TIME=snakemake@params[["source_time"]]
#SEED=as.numeric(101)
SEED=as.numeric(snakemake@params[["seed"]])

# import true tracks:
#all_tracts <- read.delim("~/Downloads/model_gf_0.3_gen_adm_50_tracts.tsv")
all_tracts <- read.delim(snakemake@input[[1]])

# import inferred tracts:
flare_mod_join_export <- read.delim(snakemake@input[[2]])
mosaic_mod_join_export <- read.delim(snakemake@input[[3]])
rfmix_mod_join_export <- read.delim(snakemake@input[[4]])
simplai_mod_join_export <- read.delim(snakemake@input[[5]])


# fixing ranges for true tracts for both haplotypes and all samples
all_tracts_hap1_range <- all_tracts %>%
  group_by(name) %>%
  filter(node_id==min(node_id)) %>%
  mutate(range=paste(left,'-',right, sep="")) %>%
  filter(left != 0 | right != 0) %>%
  select(name, source_pop, range)

all_tracts_hap2_range <- all_tracts %>%
  group_by(name) %>%
  filter(node_id==max(node_id)) %>% 
  mutate(left_2=left+sequence_length+1) %>%
  mutate(right_2=right+sequence_length+1) %>%
  mutate(range=paste(left_2,'-',right_2, sep="")) %>%  #using left_2 and right_2 since it's probably easier to consider one continuous chromosome
  filter(left != 0 | right != 0) %>%
  select(name, source_pop, range)


# merge both haplotypes
all_tracts_range <- rbind(all_tracts_hap1_range, all_tracts_hap2_range)

# filter true tracts for the two source populations
true_popb <- all_tracts_range %>% filter(source_pop=='pop_b')
true_popc <- all_tracts_range %>% filter(source_pop=='pop_c')


# filter LAI tracts for the two source populations and unassigned population:
flare_tracts_pop_b <- flare_mod_join_export %>% filter(source_pop=='pop_b')
flare_tracts_pop_c <- flare_mod_join_export %>% filter(source_pop=='pop_c')
flare_tracts_pop_no <- flare_mod_join_export %>% filter(source_pop=='pop_no')

mosaic_tracts_pop_b <- mosaic_mod_join_export %>% filter(source_pop=='pop_b')
mosaic_tracts_pop_c <- mosaic_mod_join_export %>% filter(source_pop=='pop_c')
mosaic_tracts_pop_no <- mosaic_mod_join_export %>% filter(source_pop=='pop_no')

rfmix_tracts_pop_b <- rfmix_mod_join_export %>% filter(source_pop=='pop_b')
rfmix_tracts_pop_c <- rfmix_mod_join_export %>% filter(source_pop=='pop_c')
rfmix_tracts_pop_no <- rfmix_mod_join_export %>% filter(source_pop=='pop_no')

simplai_tracts_pop_b <- simplai_mod_join_export %>% filter(source_pop=='pop_b')
simplai_tracts_pop_c <- simplai_mod_join_export %>% filter(source_pop=='pop_c')
simplai_tracts_pop_no <- simplai_mod_join_export %>% filter(source_pop=='pop_no')


# use GRanges to create approprite format
get_hap1 <- function(segment){
  GRanges(
    seqnames=segment$name,
    ranges=segment$range,
    ID=segment$name
  )
}

# get true ranges for both ancestries
get_true_b <- get_hap1(true_popb)
get_true_c <- get_hap1(true_popc)

# get flare ranges for both ancestries and for no_pop ranges (grey)
get_flare_pop_no <- get_hap1(flare_tracts_pop_no)
get_flare_b <- get_hap1(flare_tracts_pop_b)
get_flare_c <- get_hap1(flare_tracts_pop_c)

# get mosaic ranges for both ancestries and for no_pop ranges (grey)
get_mosaic_pop_no <- get_hap1(mosaic_tracts_pop_no)
get_mosaic_b <- get_hap1(mosaic_tracts_pop_b)
get_mosaic_c <- get_hap1(mosaic_tracts_pop_c)

# get rfmix ranges for both ancestries and for no_pop ranges (grey)
get_rfmix_pop_no <- get_hap1(rfmix_tracts_pop_no)
get_rfmix_b <- get_hap1(rfmix_tracts_pop_b)
get_rfmix_c <- get_hap1(rfmix_tracts_pop_c)

# get simplai ranges for both ancestries and for no_pop ranges (grey)
get_simplai_pop_no <- get_hap1(simplai_tracts_pop_no)
get_simplai_b <- get_hap1(simplai_tracts_pop_b)
get_simplai_c <- get_hap1(simplai_tracts_pop_c)


# create vector with sample names
name <- unique(all_tracts$name) 


# function to use for LAI methods, specifying the b, c and grey tracts
overlap_stats <- function(query_b, subject_b, subject_c, subject_grey, method){
  
  # define empty list
  whole_list = list()
  
  # for loop to go over all samples for this run
  for (m in 1:length(name)){
    #for (m in 3:3){
    
    target_sample <- name[m]
    
    # true tracts -----------------------------------------------------------
    query_b_id <- query_b[query_b$ID == target_sample] %>% IRanges:::reduce()
    
    # true c tracts -----------------------------------------------------------
    query_c_id <- gaps(query_b_id)
    query_c_id <- query_c_id[strand(query_c_id) == "*", ]
    
    # gaps in output ----------------------------------------------------
    subject_grey_id <- subject_grey[subject_grey$ID == target_sample] %>% IRanges:::reduce()
    
    # ancestry B --------------------------------------------------------
    subject_b_id <- subject_b[subject_b$ID == target_sample] %>% IRanges:::reduce()
    
    # ancestry C --------------------------------------------------------
    subject_c_id <- subject_c[subject_c$ID == target_sample] %>% IRanges:::reduce()
    
    # mask out gaps in true tracts --------------------------------------------
    true_b_range <- GenomicRanges::subtract(query_b_id, subject_grey_id, ignore.strand=TRUE, minoverlap=1)
    true_b <- unlist(true_b_range)
    
    true_c_range <- GenomicRanges::subtract(query_c_id, subject_grey_id, ignore.strand=TRUE, minoverlap=1)
    true_c <- unlist(true_c_range)
    
    # compute accuracy metrics ------------------------------------------------
    tp_b <- GenomicRanges::intersect(true_b, subject_b_id)
    TP_b <- sum(width(tp_b))
    
    tp_c <- GenomicRanges::intersect(true_c, subject_c_id)
    TP_c <- sum(width(tp_c))
    
    fp_b <- GenomicRanges::intersect(true_c, subject_b_id)
    FP_b <- sum(width(fp_b))
    
    fp_c <- GenomicRanges::intersect(true_b, subject_c_id)
    FP_c <- sum(width(fp_c))
    
    tn_b <- GenomicRanges::intersect(true_c, subject_c_id)
    TN_b <- sum(width(tn_b))
    
    tn_c <- GenomicRanges::intersect(true_b, subject_b_id)
    TN_c <- sum(width(tn_c))
    
    fn_b <- GenomicRanges::intersect(true_b, subject_c_id)
    FN_b <- sum(width(fn_b))
    
    fn_c <- GenomicRanges::intersect(true_c, subject_b_id)
    FN_c <- sum(width(fn_c))
    
    # put into tables ------------------------------------------------
    all_stats_length_b <- data.frame(
      sample=target_sample,
      method=method,
      pop='pop_b',
      TP=TP_b,
      FP=FP_b,
      FN=FN_b,
      TN=TN_b)
    
    all_stats_length_c <- data.frame(
      sample=target_sample,
      method=method,
      pop='pop_c',
      TP=TP_c,
      FP=FP_c,
      FN=FN_c,
      TN=TN_c)
    
    # merge the two dfs
    all_stats_length <- rbind(all_stats_length_b, all_stats_length_c)
    
    whole_list[[m]] <- all_stats_length
    
  }
  
  big_list = do.call(rbind, whole_list) # combine all together 
  
  return(big_list) # return this from function
}

flare_stats <- overlap_stats(get_true_b, get_flare_b, get_flare_c, get_flare_pop_no, 'flare')
mosaic_stats <- overlap_stats(get_true_b, get_mosaic_b, get_mosaic_c, get_mosaic_pop_no, 'mosaic')
rfmix_stats <- overlap_stats(get_true_b, get_rfmix_b, get_rfmix_c, get_rfmix_pop_no, 'rfmix')
simplai_stats <- overlap_stats(get_true_b, get_simplai_b, get_simplai_c, get_simplai_pop_no, 'simplai')


# function for additional stats:
stats_all <- function(LAI_stats){
  LAI_stats_all <- LAI_stats %>%
    group_by(sample, pop) %>%
    mutate(recall = TP/(TP+FN)) %>%
    mutate(precision =  TP/(TP+FP)) %>%
    mutate(F1 = 2*(precision*recall/(precision+recall))) %>%
    mutate(mcc = mcc(TP=TP, FP=FP, TN=TN, FN=FN)) %>%
    mutate(nMCC = (mcc+1)/2) %>%
    mutate(accuracy = (TP + TN) / (TP + TN + FP + FN)) %>% # accuracy (worst value: 0; best value: 1)
    mutate(specificity = TN / (TN+FP)) %>% # Specificity 
    mutate(FDR = FP / (TP + FP)) %>% # False discovery rate
    mutate(gf_rate=GF_RATE, source_pop_size=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIME, seed=SEED)
  return(LAI_stats_all)
}

flare_stats_all <- stats_all(flare_stats)
mosaic_stats_all <- stats_all(mosaic_stats)
rfmix_stats_all <- stats_all(rfmix_stats)
simplai_stats_all <- stats_all(simplai_stats)

# export tables
write.table(flare_stats_all, file=snakemake@output[[1]], quote=FALSE, sep='\t', row.names=FALSE)
write.table(mosaic_stats_all, file=snakemake@output[[2]], quote=FALSE, sep='\t', row.names=FALSE)
write.table(rfmix_stats_all, file=snakemake@output[[3]], quote=FALSE, sep='\t', row.names=FALSE)
write.table(simplai_stats_all, file=snakemake@output[[4]], quote=FALSE, sep='\t', row.names=FALSE)



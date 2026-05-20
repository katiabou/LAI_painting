library(dplyr)
library(ggplot2)
library(slendr)
check_dependencies(python = TRUE, quit = TRUE) # dependencies must be present

init_env(quiet = TRUE)

# popgen parameters:
Ne <- 10000
generation_time <- 1 # units of years

# sequence simulation parameters:
sequence_length <- as.numeric(snakemake@params[["genome_size_bp"]]) # megabases
recombination_rate <- 1e-8 # using the human one (does not differ much between mammals) check this link here: http://book.bionumbers.org/what-is-the-rate-of-recombination/#:~:text=In%20humans%2C%20the%20average%20rate,scale%20inversely%20with%20genomic%20length. 
#mutation_rate <- 4e-09 # based on Skoglund et al., 2015
mutation_rate <- 1e-8 # commonly used for humans


gf_rate <- as.numeric(snakemake@params[["gf_rate"]])
source_pop_size_max <- 100 # set this number to max number of samples from each source population, will subset for smaller numbers using ts_sample below
gen_adm <- as.numeric(snakemake@params[["gen_adm"]])

# seed to use 
seed_slendr <- as.numeric(snakemake@params[["seed"]])

# we'll consider a model with two source populations (coming from the same ancestor since slendr needs somewhere to begin)
# the target population branches from one of the sources, and the other source later admixes with it. 
pop <- population("pop", time = 2000, N = Ne, remove = 1499)
pop_b <- population("pop_b", time = 1500, N = Ne, parent = pop)
pop_c <- population("pop_c", time = 1500, N = Ne, parent = pop)
pop_mix <- population("pop_mix", time = gen_adm+1, N = Ne, parent = pop_b)


# specify gene flow from one ancestral population
gf <- gene_flow(from = pop_c, to = pop_mix, start = gen_adm, end = gen_adm-1, rate = gf_rate)

model <- compile_model(
  populations = list(pop, pop_b, pop_c, pop_mix),
  gene_flow = gf,
  generation_time = generation_time
)

# plot the model
png(snakemake@output[[1]])
plot_model(model, order = c("pop_b", "pop_mix", "pop", "pop_c"), proportions=TRUE)
dev.off()


# define the sampling events (basically "at which time should how many diploid individuals be sampled from which population")
source_samples <- schedule_sampling(
  model, 
  times = round(gen_adm * c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)),
  list(pop_b, source_pop_size_max),
  list(pop_c, source_pop_size_max))

# schedule sampling for target at present time:
target_samples <- schedule_sampling(
  model, 
  times = 0,
  list(pop_mix, 10))

# simulate! the time depends on the amount of sequence and the number of recorded individuals)
# this takes 100 samples from each source population
ts <- msprime(model, sequence_length, recombination_rate, samples = rbind(source_samples,target_samples), random_seed=seed_slendr) %>%
  ts_mutate(mutation_rate, random_seed=seed_slendr)

# make lists with samples from each pop
samples <- ts_samples(ts)
samples_mix <- dplyr::filter(samples, pop == "pop_mix")$name
samples_b <- dplyr::filter(samples, pop == "pop_b")$name
samples_c <- dplyr::filter(samples, pop == "pop_c")$name

# export ts data
ts_vcf(ts, snakemake@output[[2]], chrom ="chr1")

# export metadata
meta_t <- ts_samples(ts)
write.table(meta_t, file=snakemake@output[[3]], quote=FALSE, sep='\t', row.names=FALSE)


####### GET TRACKS #######

# this gets tracks for pop_c since through the slendr model, it's the pop from which "gene flow" happens into the mix population, which seperates from the other ancestry (pop_b)
popc_tracts <- ts_tracts(ts, census = gen_adm, source = "pop_c")


# make tracts at beginning of chrom to assign to other ancestry (anc_b) in case it hasn't been assigned to anything
popc_tracts_rep <- function(poptracts, pop_1){
  blabla <- poptracts %>%
  group_by(name, node_id) %>% # group by name and node (so per haplotype per sample)
  dplyr::slice(1L, row_number()) %>% # replicate the first row for each sample/haplotype combo
  group_by_all() %>% # then group by everything
  filter(n()>1) %>% # filter for the duplicates
  distinct() %>% # extract distinct duplicate
  ungroup() %>% 
  mutate(left_2 = 0, right_2 = if_else(left==0, 0, left-1)) %>% #put into new df and then only assign new tract to other ancestry if they don't start with 0.
  select(-left, -right)
  
  blabla$length <- blabla$right_2-blabla$left_2
  blabla$source_pop <- pop_1
  names(blabla)[names(blabla) == 'left_2'] <- 'left'
  names(blabla)[names(blabla) == 'right_2'] <- 'right'
  return(blabla)
}

ttt <- popc_tracts_rep(popc_tracts, 'pop_b')


# make new coordinates for end of chromosomes to fill in second ancestry if no ancestry is present
end_tracts <- function(poptracts, pop_1){
  blabla2 <- poptracts %>% 
  group_by(name, node_id) %>%
  mutate(left_2 = right+1, right_2=lead(left)-1) 
  blabla2[is.na(blabla2)] <- sequence_length+1 #fix NAs with max mb size
  
  blabla3 <- blabla2 %>% select(-left, -right)
  blabla3$length <- blabla3$right_2-blabla3$left_2
  blabla3$source_pop <- pop_1
  
  names(blabla3)[names(blabla3) == 'left_2'] <- 'left'
  names(blabla3)[names(blabla3) == 'right_2'] <- 'right'
  
  return(blabla3)
}

fff <- end_tracts(popc_tracts,'pop_b')


# remove other columns from initial ancestry 1 coordinates
# remove_other_col <- function(poptracts){
#   blabla4 <- poptracts %>% 
#   group_by(name, node_id) %>%
#   mutate(left_2 = right+1, right_2=lead(left)-1)
#   blabla4[is.na(blabla4)] <- sequence_length
#   blabla5 <- blabla4 %>%select(-left_2, -right_2)
# }

# nnn <- remove_other_col(popc_tracts)


# merge all df together (missing tracts at the beginning, missing tracts at the end and the initial tracts)
#all_tracts <- rbind(ttt, fff, nnn)
all_tracts <- rbind(ttt, fff, popc_tracts)

# add info from run to track file
all_tracts_info <- all_tracts %>%
                   mutate(gen_adm=gen_adm, gf=gf_rate, seed=seed_slendr)
write.table(all_tracts_info, file=snakemake@output[[4]], quote=FALSE, sep='\t', row.names=FALSE)


####### FST of sources and targets #######

# make a function to extract the fst values from each time point for a given adm time
source_times <- unique(meta_t$time)

# fst between sources
get_fst <- function(sampling_time){
  pop_b_meta <- meta_t %>% filter(pop=='pop_b', time==sampling_time) %>% select(name)
  pop_c_meta <- meta_t %>% filter(pop=='pop_c', time==sampling_time) %>% select(name)
  fst <- ts_fst(ts, sample_sets = list(pop_b = pop_b_meta$name,
                                       pop_c = pop_c_meta$name))
  fin <- data.frame(
    fst=fst$Fst,
    source_time=sampling_time,
    gen_adm=gen_adm,
    seed=seed_slendr,
    comparison='source1_source2')
}

gg <- lapply(source_times, get_fst)
hh <- as.data.frame(do.call(rbind, gg))

# fst between sources and targets

# source 1
get_fst_1 <- function(sampling_time){
  pop_b_meta <- meta_t %>% filter(pop=='pop_b', time==sampling_time) %>% select(name)
  pop_mix_meta <- meta_t %>% filter(pop=='pop_mix', time==0) %>% select(name)
  fst <- ts_fst(ts, sample_sets = list(pop_b = pop_b_meta$name,
                                       pop_mix = pop_mix_meta$name))
  fin <- data.frame(
    fst=fst$Fst,
    source_time=sampling_time,
    gen_adm=gen_adm,
    seed=seed_slendr,
    comparison='source1_target')
}

gg1 <- lapply(source_times, get_fst_1)
hh1 <- as.data.frame(do.call(rbind, gg1))


# source 2
get_fst_2 <- function(sampling_time){
  pop_c_meta <- meta_t %>% filter(pop=='pop_c', time==sampling_time) %>% select(name)
  pop_mix_meta <- meta_t %>% filter(pop=='pop_mix', time==0) %>% select(name)
  fst <- ts_fst(ts, sample_sets = list(pop_c = pop_c_meta$name,
                                       pop_mix = pop_mix_meta$name))
  fin <- data.frame(
    fst=fst$Fst,
    source_time=sampling_time,
    gen_adm=gen_adm,
    seed=seed_slendr,
    comparison='source2_target')
}

gg2 <- lapply(source_times, get_fst_2)
hh2 <- as.data.frame(do.call(rbind, gg2))


# combine all together
all_fst <- rbind(hh, hh1, hh2)

# export 
#write.table(all_fst, file='~/Downloads/fst_sources_target.txt', quote=FALSE, sep='\t', row.names=FALSE)
write.table(all_fst, file=snakemake@output[[5]], quote=FALSE, sep='\t', row.names=FALSE)

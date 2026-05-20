library(dplyr)
library(gsubfn)

#meta <- read.delim("~/Downloads/model_gf_0.3_gen_adm_200_meta.tsv")
meta <- read.delim(snakemake@input[[1]])

source_size <- as.numeric(snakemake@params[["source_size"]]) #this would be from each source, and would be restricted to 4,8,12,16
gen_adm  <- as.numeric(snakemake@params[["gen_adm"]])

# names of source populations and target population
source_1 = snakemake@params[["source_1"]]
source_2 = snakemake@params[["source_2"]]
target = snakemake@params[["target"]]


#sampling scenarios
source_times_1 <- c(0.2,0.6,1.4,1.8) #sample from these specific time periods
source_times <- source_times_1*gen_adm


#I want to take 100 from each source pop 4 time periods before and after admixture
source_size_per_time <- source_size/4

#remove target sample and extract specific number per source and time period:
meta_source_sub <- meta %>% 
  filter(pop!=target) %>%
  filter(time %in% source_times) %>%
  group_by(pop, time) %>% 
  slice_head(n=source_size_per_time) %>%
  ungroup

#keep only target samples
meta_target <- meta %>% 
  filter(pop==target)


####  EXPORTING FORMATS FOR EACH LAI METHOD ####

# FLARE and RFMIX

#export target sample names (used for FLARE and RFMIX)
meta_target_name <- meta_target %>% 
                    select(name)

#export source sample names and population (used for FLARE and RFMIX)
meta_source_name_pop <- meta_source_sub %>% 
                        select(name, pop)


#export source sample names (used for RFMIX)
meta_source_name <- meta_source_sub %>% 
                    select(name)

#MOSAIC 

#source population 1 sample names
meta_source_name_1 <- meta_source_sub %>% 
                          filter(pop==source_1)%>%
                          select(name)

#source population 2 sample names
meta_source_name_2 <- meta_source_sub %>% 
                          filter(pop==source_2)%>%
                          select(name)


#source population 1 and sample names
meta_source_pop_name_pop_1 <- meta_source_sub %>% 
                              filter(pop==source_1)%>%
                              select(pop,name)

#source population 2 and sample names
meta_source_pop_name_pop_2 <- meta_source_sub %>% 
                              filter(pop==source_2)%>%
                              select(pop,name)

#target population and sample names
meta_target_pop_name <- meta_target %>% 
                        select(pop,name)

#this is needed to run MOSAIC
meta_target_source_pop_name <- rbind(meta_source_pop_name_pop_1, meta_source_pop_name_pop_2, meta_target_pop_name)


#SIMPLAI 

#source and target sample names
meta_target_source_name <- meta_target_source_pop_name %>%
                           select(name)



write.table(meta_target_name, file=snakemake@output[[1]], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
#write.table(meta_target_name, file='~/Downloads/test.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

write.table(meta_source_name_pop, file=snakemake@output[[2]], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
#write.table(meta_source_name_pop, file='~/Downloads/test.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

write.table(meta_source_name, file=snakemake@output[[3]], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
#write.table(meta_source_name, file='~/Downloads/test.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

write.table(meta_source_name_1, file=snakemake@output[[4]], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
#write.table(meta_source_name_1, file='~/Downloads/test.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

write.table(meta_source_name_2, file=snakemake@output[[5]], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
#write.table(meta_source_name_2, file='~/Downloads/test.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

write.table(meta_target_source_pop_name, file=snakemake@output[[6]], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
#write.table(meta_target_source_pop_name, file='~/Downloads/test.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

write.table(meta_target_source_name, file=snakemake@output[[7]], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
#write.table(meta_target_source_name, file='~/Downloads/test.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

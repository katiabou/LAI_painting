library(MOSAIC)

args <- commandArgs(trailingOnly = TRUE)

#import result files
la_results = args[1]
#la_results = '~/Downloads/mosaic_test/localanc_pop_mix_2way_1-10_1-1_420_204.312_0.99_100.RData'

model_results = args[2]
#model_results = '~/Downloads/mosaic_test/pop_mix_2way_1-10_1-1_420_204.312_0.99_100.RData'

#import run parameters
gf <- args[3]
#gf <- 0.2
source <- args[4]
#source <- 20
gen_adm <- args[5]
#gen_adm <- 50
sample_number <- as.numeric(args[6])

source_time <- args[7]

seed <- args[8]

#load files
load(model_results)
load(la_results)

#this is to get the adm dates:
test <- plot_coanccurves(acoancs,dr)
lamdas <- test$gens.matrix

date <- (lamdas[1,1]+lamdas[1,2]+lamdas[2,2])/3

#this is to get the adm proportions:
acoancs$ancprobs #take each entry and divide by 2

vec1 <- (acoancs$ancprobs[1,]/2)
vec2 <- (acoancs$ancprobs[2,]/2)

allvec1 <- sum(vec1)/sample_number
allvec2 <- sum(vec2)/sample_number

#export into table 
method <- 'MOSAIC'

mosaic_adm_time <- data.frame(V1=method, V2=gf, V3=source, V4=gen_adm, V5=source_time, V6=seed, V7=date)

mosaic_adm_prop <- data.frame(V1=method, V2=gf, V3=source, V4=gen_adm, V5=source_time, V6=seed, V7=allvec1, V8=allvec2)

#export
write.table(mosaic_adm_time, file=args[9], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)
write.table(mosaic_adm_prop, file=args[10], quote=FALSE, sep='\t', row.names=FALSE, col.names = FALSE)

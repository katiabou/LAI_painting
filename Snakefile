#this is for the real data analysis
import pandas as pd

##### load config file #####
configfile: 'config.yaml'

##### set wildcard constraints ##### 
wildcard_constraints:
    source="\d+",
    source_time=r"\d+(\.\d+)?",
    gf="\d+.\d+",
    gen_adm="\d+",
    prob="\d+.\d+",
    seed="\d+"

GF_RATE=[0.1, 0.3]
SOURCE_POP_SIZE=[2,4,6,8,10,12,14,16,20,40,50,80,100]
GEN_ADM=[50, 100, 200, 300, 400, 500]
SOURCE_TIMES=[0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2]
PROB=[0.0]
PROB_CUTOFF=[0.6, 0.7, 0.8, 0.9, 0.999]
SEQ_LENGTH=[300000000]
SEED=[101,102,103,104,105,106,107,108,109,110]
#SEED=config["seeds"] #this is to run specified values of the seed from the command line (or within the lai_painting.sh script called to run snakemake, in the form of --config seeds=${seeds})


##### Rules #####
include: "rules/slendr_simulations.smk"
include: "rules/LAI_tools.smk"
include: "rules/LAI_tracts_accuracy.smk"
include: "rules/LAI_tools_grouped_sources.smk"
include: "rules/LAI_tracts_accuracy_grouped_sources.smk"
include: "rules/dates_v4010.smk"
include: "rules/dates_v4010_grouped.smk"
include: "rules/adm_dating_proportion_comparison.smk"


rule all:
    input:
        expand( ## Run slendr simulations
             [
                #"output/track_length_slendr/tracts_gf_{gf}_histogram_length.png",
                #"output/fst_source_target/fst_gf_{gf}_sources.png",    
                #"output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt-PC1_PC2.png"
            ],
            gf=GF_RATE, 
            gen_adm=GEN_ADM, 
            source=SOURCE_POP_SIZE,
            source_time=SOURCE_TIMES,
            seed=SEED, 
        ),
        ## Compare inferred vs true for all LAI methods
        expand(
            [
                #"output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_nMCC_heatmap.png",
                #"output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_sensitivity_heatmap.png",
            ],
            gf=GF_RATE, 
            gen_adm=GEN_ADM, 
            source=SOURCE_POP_SIZE,
            source_time=SOURCE_TIMES,
            seed=SEED, 
            prob=PROB
        ),
        ## Get tract length distribution, posterior probability and admixture proportion plots
        expand( #SOS this should contain a subset of the wildcards!!!
            [
                # "output/compare_software/plots_all_LAI/model_gf_{gf}_source_subset_gen_adm_all_source_time_0_prob_{prob}-tract_lengths.png",
                # "output/compare_software/plots_all_LAI/model_gf_{gf}_source_sub_gen_adm_all_source_time_0_all_tracts_probcutoff_percentage.png",
                # "output/compare_software/plots/posterior_prob_gf_{gf}_pop_C_subset.png",
                # "output/compare_software/plots/adm_prop_raw_{gf}_pop_C_subset.png",
            ],
            gf=GF_RATE,
            seed=SEED, 
            source=[4,10,20,100], 
            gen_adm=GEN_ADM, 
            source_time=0, 
            prob=PROB
        ),
        ## Run DATES
        expand(
            [
                "output/compare_software/plots/all_simulated_adm_times_gf_{gf}_ribbon.png",
                "output/seed_{seed}/dates/plots/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_curve.png",
                "output/compare_software/dates/dates_gf_{gf}_nrmsd_all.png",
            ],
            gf=GF_RATE, 
            seed=SEED, 
            source=SOURCE_POP_SIZE, 
            gen_adm=GEN_ADM, 
            source_time=SOURCE_TIMES,
            prob=PROB
        ),
        expand( # Combo figure, keep this specific wildcards
            [
                #"output/compare_software/plots_all_LAI/model_gf_{gf}_seed_{seed}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_prob_{prob}-combo_figure.png",
            ],
            gf=[0.3],
            seed=[101],
            source=[10],
            gen_adm=[50,200],
            source_time=[0],
            prob=[0.0]
        ),
        ## Run LAI and DATES on grouped sources, be careful to specify specific source pop sizes
        expand(
            [
                #"output/compare_software/plots_all_LAI_grouped_sources/model_gf_{gf}_all_source_all_gen_adm_grouped_source_time_all_methods_prob_{prob}_nMCC_heatmap.png",
                #"output/seed_{seed}/dates/plots_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_curve.png",
                #"output/compare_software/plots_grouped_sources/all_simulated_adm_times_gf_{gf}_ribbon.png",
            ],
            gf=GF_RATE, 
            gen_adm=GEN_ADM, 
            source=[4,8,12,16],
            seed=SEED, 
            prob=PROB
        ),

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
#SEED=config["seeds"] #this is to run specified values of the seed from the command line (or within the bash.sh script called to run snakemake, in the form of --config seeds=${seeds})


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
        expand( #SOS this should contain a subset of the wildcards!!!
            [
               #"output/compare_software/plots_all_LAI/model_gf_{gf}_source_sub_gen_adm_all_source_time_0_all_tracts_probcutoff_percentage.png",
               #"output/compare_software/plots/posterior_prob_gf_{gf}_pop_C_subset.png",
               # "output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv",
               #"output/compare_software/plots/adm_prop_raw_{gf}_pop_C_subset.png"

            ],
            gf=GF_RATE,
            seed=SEED, 
            source=[4,10,20,100], 
            gen_adm=GEN_ADM, 
            source_time=0, 
            prob_cutoff=PROB_CUTOFF
        ),
        ### Run slendr simulations and then the LAI methods
        expand(
             [
                # "output/track_length_slendr/tracts_gf_{gf}_histogram_length.png",
                #"output/fst_source_target/fst_gf_{gf}_all.png",
                #"output/fst_source_target/fst_gf_{gf}_sources.png"    
                # "output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt-PC1_PC2.png"
        #         "output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv",
        #         "output/seed_{seed}/flare/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt",
        #         "output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/localanc_pop_mix.RData",
        #         "output/seed_{seed}/rfmix/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.fb.tsv",
        #         "output/seed_{seed}/simplai/output/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_n_2000_m_1000_t_5_fromRec_withSingl.adm",
        #         #"output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt-PC1_PC2.png", #RUN THESE LATER ON, MIGHT INCLUDE
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
                #"output/seed_{seed}/compare_software/plots/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_all_method_tracts_prob_{prob}.png",
                #"output/compare_software/plots_all_LAI/model_gf_{gf}_source_subset_gen_adm_all_source_time_0_prob_{prob}-tract_lengths.png",
                #"output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_nMCC_heatmap.png",
                #"output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_sensitivity_heatmap.png"
            ],
            gf=GF_RATE, 
            gen_adm=GEN_ADM, 
            source=SOURCE_POP_SIZE,
            source_time=SOURCE_TIMES,
            seed=SEED, 
            prob=PROB,
            # gf=GF_RATE, 
            # gen_adm=100, 
            # source=100,
            # source_time=0,
            # seed=101, 
            # prob=PROB
        ),
        ## Run LAI on grouped source run, specify specific source pop sizes, since I don't use all of them
        expand(
            [
                #"output/compare_software/plots_all_LAI_grouped_sources/model_gf_{gf}_all_source_all_gen_adm_grouped_source_time_all_methods_prob_{prob}_nMCC_heatmap.png",
                #"output/compare_software/plots_grouped_sources/all_simulated_adm_times_gf_{gf}.png",
                #"output/seed_{seed}/dates/plots_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_curve.png",
                #"output/seed_{seed}/dates/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-nrmsd_dates.txt"
                #"output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_adm_time.txt",
                #"output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_adm_time.txt"
            ],
            gf=GF_RATE, 
            gen_adm=GEN_ADM, 
            source=[4,8,12,16],
            seed=SEED, 
            prob=PROB
        ),
        expand(
             [
                #"output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_nMCC_heatmap.png",
        #         "output/compare_software/plots/poster_simulated_adm_times_gf_0.3.png",
        #         "output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_nMCC_heatmap-poster.png",
                  "output/compare_software/plots/all_simulated_adm_times_gf_{gf}.png",
        #         #"output/seed_{seed}/dates/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}-mean_dates.txt",
                  #"output/seed_{seed}/dates/plots/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_curve.png",
                  #"output/compare_software/dates/dates_gf_{gf}_nrmsd_all.png",
        #         #"output/seed_{seed}/dates/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}-nrmsd_dates.txt"
        #         #"output/seed_{seed}/dates/files/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/pop_mix:log",
            ],
            gf=GF_RATE, 
            seed=SEED, 
            source=SOURCE_POP_SIZE, 
            gen_adm=GEN_ADM, 
            source_time=SOURCE_TIMES,
            prob=PROB
        ),
       #grouped sources, SOS have to be careful with the input files, don't want them for all wildcard values
        # expand(
        #     [
        #         "output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_tracts_overlap_stats_prob_{prob}.tsv",
        #         "output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_tracts_overlap_stats_prob_{prob}.tsv",
        #         "output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_rfmix_tracts_overlap_stats_prob_{prob}.tsv",
        #         "output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_simplai_tracts_overlap_stats_prob_{prob}.tsv",
        #     ],
        # gf=GF_RATE, 
        # seed=SEED, 
        # source=[4,8,12,16], 
        # gen_adm=GEN_ADM,
        # prob=PROB,
        # ),
        # expand(
        #     [
        #         "output/empirical_data/humans/rfmix/plots/rfmix_allchrom_MAF_{maf}_recalibrated_INFO_{info}_gen_adm_{gen_adm_humans}_EBA1.png",
        #         "output/empirical_data/humans/mosaic/plots/mosaic_allchrom_MAF_{maf}_recalibrated_INFO_{info}_gen_adm_{gen_adm_humans}_EBA1.png",
        #         "output/empirical_data/humans/mosaic/plots_flat_rate/mosaic_allchrom_MAF_{maf}_recalibrated_INFO_{info}_gen_adm_{gen_adm_humans}_EBA1.png",
        #         "output/empirical_data/humans/mosaic/plots_no_gen/mosaic_allchrom_MAF_{maf}_recalibrated_INFO_{info}_EBA1.png",
        #         "output/empirical_data/humans/flare/plots/flare_allchrom_MAF_{maf}_recalibrated_INFO_{info}_gen_adm_{gen_adm_humans}_EBA1.png",
        #         "output/empirical_data/humans/dates/files/dates_allchrom_MAF_{maf}_recalibrated_INFO_{info}/output_dates.out",
        #         "output/empirical_data/humans/dates/files_flat_rm/dates_allchrom_MAF_{maf}_recalibrated_INFO_{info}/output_dates.log"
        #     ],
        #     gf=GF_RATE, 
        #     seed=SEED, 
        #     source=SOURCE_POP_SIZE, 
        #     gen_adm=GEN_ADM, 
        #     source_time=SOURCE_TIMES,
        #     prob=PROB,
        #     chrom_hum=CHROM_HUMAN,
        #     maf=config['maf'],
        #     info=config['info'],
        #     gen_adm_humans=GEN_ADM_HUMANS
        # ),
        # expand(
        #     [
        #         "output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_prob_{prob}.tsv",
        #     ],
        # gf=GF_RATE, 
        # seed=SEED, 
        # source=SOURCE_POP_SIZE, 
        # source_time=SOURCE_TIMES,
        # gen_adm=GEN_ADM,
        # prob=PROB,
        # ),
#grouped sources version 2, SOS have to be careful with the input files, don't want them for all wildcard values
        # expand(
        #     [
        #         "output/seed_{seed}/compare_software/files_grouped_sources_v2/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_v2_flare_tracts_overlap_stats_prob_{prob}.tsv",
        #         "output/seed_{seed}/compare_software/files_grouped_sources_v2/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_v2_mosaic_tracts_overlap_stats_prob_{prob}.tsv",
        #         "output/seed_{seed}/compare_software/files_grouped_sources_v2/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_v2_rfmix_tracts_overlap_stats_prob_{prob}.tsv",
        #         "output/seed_{seed}/compare_software/files_grouped_sources_v2/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_v2_simplai_tracts_overlap_stats_prob_{prob}.tsv",
        #         "output/compare_software/plots_all_LAI_grouped_sources_v2/model_gf_{gf}_all_source_all_gen_adm_grouped_source_time_all_methods_prob_{prob}_nMCC_heatmap.png"
        #     ],
        # gf=GF_RATE, 
        # seed=SEED, 
        # source=[4,8,12,16], 
        # gen_adm=GEN_ADM,
        # prob=PROB,
        # ),

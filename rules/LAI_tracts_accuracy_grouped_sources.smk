########################################
#  Plot LAI output for grouped sources #
########################################


rule extract_inferred_tracts_grouped_sources:
    """ 
    Extract and plot the true tracts and those inferred from the LAI tools
    """
    input:
        tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        meta_source_name_pop = 'output/seed_{seed}/files_grouped_sources/source_name_pop_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_target_name = 'output/seed_{seed}/files_grouped_sources/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        flare_txt = 'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        la_results = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/MOSAIC_RESULTS/localanc_pop_mix.RData',
        model_results = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/MOSAIC_RESULTS/pop_mix.RData',
        modern_imputed_snps = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/snpfile.1',
        rfmix_out = 'output/seed_{seed}/rfmix/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.fb.tsv',
        simplai_output_rec = 'output/seed_{seed}/simplai/output_grouped_sources/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_n_2000_m_1000_t_5_fromRec_withSingl.adm'
    output:
        all_tracts = 'output/seed_{seed}/compare_software/plots_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_all_method_tracts_prob_{prob}.png',
        flare_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_tracts_prob_{prob}.tsv',
        mosaic_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_tracts_prob_{prob}.tsv',
        rfmix_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_rfmix_tracts_prob_{prob}.tsv',
        simplai_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_simplai_tracts_prob_{prob}.tsv'
    params:
        seq_length = SEQ_LENGTH,
        prob_cutoff = '{prob}',
        dir_mosaic = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/',
        gf_rate = '{gf}',
        source_pop_size = '{source}',
        gen_adm = '{gen_adm}',
        seed = '{seed}'
    script:
        "../scripts/extract_inferred_tracts.R"

    
rule get_overlap_stats_grouped_sources:
    """ 
    Get stats from overlaps 
    """
    input:
        tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        flare_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_tracts_prob_{prob}.tsv',
        mosaic_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_tracts_prob_{prob}.tsv',
        rfmix_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_rfmix_tracts_prob_{prob}.tsv',
        simplai_tracts = 'output/seed_{seed}/compare_software/files_grouped_sources/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_simplai_tracts_prob_{prob}.tsv'
    output:
        flare_tracts_overlap = 'output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_tracts_overlap_stats_prob_{prob}.tsv',
        mosaic_tracts_overlap = 'output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_tracts_overlap_stats_prob_{prob}.tsv',
        rfmix_tracts_overlap = 'output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_rfmix_tracts_overlap_stats_prob_{prob}.tsv',
        simplai_tracts_overlap = 'output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_simplai_tracts_overlap_stats_prob_{prob}.tsv'
    params:
        gf_rate = '{gf}',
        source_pop_size = '{source}',
        gen_adm = '{gen_adm}',
        seq_length = SEQ_LENGTH,
        source_time = 'grouped_sources',
        seed = '{seed}'
    script:
        "../scripts/overlap_stats.R"


rule plot_all_grouped_sources:
    """ 
    Plot overlap stats such as nMCC across demographic scenarios and LAI methods
    Specifying source size to only include 4,8,12 and 16 
    """
    input:
        flare_tracts_overlap = expand('output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, allow_missing=True),
        mosaic_tracts_overlap = expand('output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, allow_missing=True),
        rfmix_tracts_overlap = expand('output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_rfmix_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, allow_missing=True),
        simplai_tracts_overlap = expand('output/seed_{seed}/compare_software/files_grouped_sources/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_simplai_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, allow_missing=True),
    output:
        heatmap_nmcc = 'output/compare_software/plots_all_LAI_grouped_sources/model_gf_{gf}_all_source_all_gen_adm_grouped_source_time_all_methods_prob_{prob}_nMCC_heatmap.png', 
        heatmap_sensitivity = 'output/compare_software/plots_all_LAI_grouped_sources/model_gf_{gf}_all_source_all_gen_adm_grouped_source_time_all_methods_prob_{prob}_sensitivity_heatmap.png',    
    params:
        prob = '{prob}',
        gf = '{gf}'
    script:
        "../scripts/all_plots_grouped_sources.R"




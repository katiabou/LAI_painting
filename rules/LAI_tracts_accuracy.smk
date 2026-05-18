#####################
#  Plot LAI output  #
#####################


rule extract_inferred_tracts:
    """ 
    Extract and plot the true tracts and those inferred from the LAI tools
    """
    input:
        tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        meta_source_name_pop = 'output/seed_{seed}/files/source_name_pop_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        meta_target_name = 'output/seed_{seed}/files/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        flare_txt = 'output/seed_{seed}/flare/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        la_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/localanc_pop_mix.RData',
        model_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/pop_mix.RData',
        modern_imputed_snps = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/snpfile.1',
        rfmix_out = 'output/seed_{seed}/rfmix/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.fb.tsv',
        simplai_output_rec = 'output/seed_{seed}/simplai/output/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_n_2000_m_1000_t_5_fromRec_withSingl.adm'
    output:
        all_tracts_plot = 'output/seed_{seed}/compare_software/plots/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_all_method_tracts_prob_{prob}.png',
        flare_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_prob_{prob}.tsv',
        mosaic_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_prob_{prob}.tsv',
        rfmix_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_prob_{prob}.tsv',
        simplai_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_prob_{prob}.tsv'
    params:
        seq_length = SEQ_LENGTH,
        prob_cutoff = '{prob}',
        dir_mosaic = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/',
        gf_rate = '{gf}',
        source_pop_size = '{source}',
        gen_adm = '{gen_adm}',
        source_time = '{source_time}',
        seed = '{seed}'
    script:
        "../scripts/extract_inferred_tracts.R"


rule plot_tract_lengths:
    """ 
    Plot inferred and true tract lengths for one source time point (0), and for four population sizes (4,10,20,100)
    """
    input:
        tracts = expand('output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv', seed=SEED, gen_adm=GEN_ADM, allow_missing=True),
        flare_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_0_flare_tracts_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, allow_missing=True),
        mosaic_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_0_mosaic_tracts_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, allow_missing=True),
        rfmix_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_0_rfmix_tracts_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, allow_missing=True),
        simplai_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_0_simplai_tracts_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, allow_missing=True),
    output:
        tract_lengths_plot = 'output/compare_software/plots_all_LAI/model_gf_{gf}_source_subset_gen_adm_all_source_time_0_prob_{prob}-tract_lengths.png',
    params:
        gf_rate = '{gf}',
        prob_cutoff = '{prob}',
        seq_length = SEQ_LENGTH,
        source_time = 0
    script:
        "../scripts/tracts_length_plots.R"



rule get_overlap_stats:
    """ 
    Get stats from overlaps 
    """
    input:
        tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        flare_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_prob_{prob}.tsv',
        mosaic_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_prob_{prob}.tsv',
        rfmix_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_prob_{prob}.tsv',
        simplai_tracts = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_prob_{prob}.tsv'
    output:
        flare_tracts_overlap = 'output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_overlap_stats_prob_{prob}.tsv',
        mosaic_tracts_overlap = 'output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_overlap_stats_prob_{prob}.tsv',
        rfmix_tracts_overlap = 'output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_overlap_stats_prob_{prob}.tsv',
        simplai_tracts_overlap = 'output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_overlap_stats_prob_{prob}.tsv'
    params:
        gf_rate = '{gf}',
        source_pop_size = '{source}',
        gen_adm = '{gen_adm}',
        seq_length = SEQ_LENGTH,
        source_time = '{source_time}',
        seed = '{seed}'
    script:
        "../scripts/overlap_stats.R"


rule plot_overlap_stats:
    """ 
    Plot overlap stats across demographic scenarios and LAI methods 
    """
    input:
        flare_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIMES, allow_missing=True),
        mosaic_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIMES, allow_missing=True),
        rfmix_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIMES, allow_missing=True),
        simplai_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIMES, allow_missing=True),
    output:
        heatmap_nmcc = 'output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_nMCC_heatmap.png',
        heatmap_sen = 'output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_sensitivity_heatmap.png',
        heatmap_spe = 'output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_specificity_heatmap.png',
        heatmap_pre = 'output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_precision_heatmap.png',
        heatmap_acc = 'output/compare_software/plots_all_LAI/model_gf_{gf}_all_source_all_gen_adm_all_source_time_all_methods_prob_{prob}_accuracy_heatmap.png',
    params:
        prob = '{prob}',
        gf = '{gf}'
    script:
        "../scripts/all_plots_LAI.R"


rule percentage_seq_pp:
    """
    Filter tracts based on different pp cutoffs
    """
    input:
        tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        meta_source_name_pop = 'output/seed_{seed}/files/source_name_pop_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        meta_target_name = 'output/seed_{seed}/files/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        flare_txt = 'output/seed_{seed}/flare/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        la_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/localanc_pop_mix.RData',
        model_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/pop_mix.RData',
        modern_imputed_snps = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/snpfile.1',
        rfmix_out = 'output/seed_{seed}/rfmix/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.fb.tsv',
        simplai_output_rec = 'output/seed_{seed}/simplai/output/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_n_2000_m_1000_t_5_fromRec_withSingl.adm'
    output:
        all_tracts_plot = 'output/seed_{seed}/compare_software/plots_pp_cutoff/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_all_method_tracts_probcutoff_{prob_cutoff}.png',
        flare_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_probcutoff_{prob_cutoff}.tsv',
        mosaic_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_probcutoff_{prob_cutoff}.tsv',
        rfmix_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_probcutoff_{prob_cutoff}.tsv',
        simplai_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_probcutoff_{prob_cutoff}.tsv'
    params:
        seq_length = SEQ_LENGTH,
        prob_cutoff = '{prob_cutoff}',
        dir_mosaic = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/',
        gf_rate = '{gf}',
        source_pop_size = '{source}',
        gen_adm = '{gen_adm}',
        source_time = '{source_time}',
        seed = '{seed}'
    script:
        "../scripts/extract_inferred_tracts.R"



rule get_overlap_stats_pp:
    """ 
    Get stats from overlaps 
    """
    input:
        tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        flare_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_probcutoff_{prob_cutoff}.tsv',
        mosaic_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_probcutoff_{prob_cutoff}.tsv',
        rfmix_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_probcutoff_{prob_cutoff}.tsv',
        simplai_tracts = 'output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_probcutoff_{prob_cutoff}.tsv'
    output:
        flare_tracts_overlap = 'output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv',
        mosaic_tracts_overlap = 'output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv',
        rfmix_tracts_overlap = 'output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv',
        simplai_tracts_overlap = 'output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv',
    params:
        gf_rate = '{gf}',
        source_pop_size = '{source}',
        gen_adm = '{gen_adm}',
        seq_length = SEQ_LENGTH,
        source_time = '{source_time}',
        seed = '{seed}'
    script:
        "../scripts/overlap_stats.R"


rule plot_percentage_seq_pp:
    """
    Plot the percentage of missing sequences after applying different pp cutoffs on each method (when applicable)
    """
    input:
        flare_tracts = expand('output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=0, prob_cutoff=PROB_CUTOFF, allow_missing=True),
        mosaic_tracts = expand('output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=0, prob_cutoff=PROB_CUTOFF, allow_missing=True),
        rfmix_tracts = expand('output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=0, prob_cutoff=PROB_CUTOFF, allow_missing=True),
        flare_tracts_overlap = expand('output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=0, prob_cutoff=PROB_CUTOFF, allow_missing=True),
        mosaic_tracts_overlap = expand('output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=0, prob_cutoff=PROB_CUTOFF, allow_missing=True),
        rfmix_tracts_overlap = expand('output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=0, prob_cutoff=PROB_CUTOFF, allow_missing=True),
    output:
        pp_percentage = 'output/compare_software/plots_all_LAI/model_gf_{gf}_source_sub_gen_adm_all_source_time_0_all_tracts_probcutoff_percentage.png'
    params:
        seq_length = SEQ_LENGTH,
        gf_rate = '{gf}',
    script:
        "../scripts/plot_pp_percentage_sequence.R"



rule plot_combo_examples:
    """
    Plot a selected combination of the above plots into main figure for 2 demographic scenarios
    """
    input:
        tracts1 = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        flare_tracts1 = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_prob_{prob}.tsv',
        mosaic_tracts1 = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_prob_{prob}.tsv',
        rfmix_tracts1 = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_prob_{prob}.tsv',
        simplai_tracts1 = 'output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_prob_{prob}.tsv',
        all_adm_prop = 'output/seed_{seed}/compare_software/files/posterior_probabilities/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}-all_simulated_post_prob.txt',
        dates_mix_expfit = 'output/seed_{seed}/dates/files/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/pop_mix:expfit.out',
        dates_mean = 'output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_dates_adm_time.txt',
        dates_nrmsd = 'output/seed_{seed}/dates/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}-nrmsd_dates.txt'
    output:
        combo_plot = 'output/compare_software/plots_all_LAI/model_gf_{gf}_seed_{seed}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_prob_{prob}-combo_figure.png'
    params:
        seq_length = SEQ_LENGTH,
        gf_rate = '{gf}',
        source_time = '{source_time}',
        source_pop_size = '{source}',
        gen_adm = '{gen_adm}',
        seed = '{seed}', 
        prob = '{prob}'
    script:
        "../scripts/combo_plot_scenarios.R"

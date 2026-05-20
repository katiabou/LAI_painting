#######################
#  Get adm dates plot #
#######################


rule flare_adm_time_prop:
    """
    Get admixture time and proportions from FLARE 
    """
    input:
        flare_model = 'output/seed_{seed}/flare/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.model'
    output:
        flare_adm_time = 'output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_adm_time.txt',
        flare_adm_prop = 'output/seed_{seed}/compare_software/files/admixture_proportions/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_adm_prop.txt'
    shell:
        '''
        date=$(grep -w T: {input.flare_model} -A 1 | sed -n '2p')

        prop=$(grep -w mu {input.flare_model} -A 1 | sed -n '2p')

        echo 'flare' {wildcards.gf} {wildcards.source} {wildcards.gen_adm} {wildcards.source_time} {wildcards.seed} $date > {output.flare_adm_time}

        echo 'flare' {wildcards.gf} {wildcards.source} {wildcards.gen_adm} {wildcards.source_time} {wildcards.seed} $prop > {output.flare_adm_prop}
        '''

rule mosaic_adm_time_prop:
    """ 
    Get admixture time and proportions from MOSAIC 
    """
    input:
        la_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/localanc_pop_mix.RData',
        model_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/pop_mix.RData',
    output:
        mosaic_adm_time = 'output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_adm_time.txt',
        mosaic_adm_prop = 'output/seed_{seed}/compare_software/files/admixture_proportions/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_adm_prop.txt'
    params:
        gf_rate = '{gf}',
        source_num = '{source}',
        gen_adm = '{gen_adm}',
        sample_number = 10,
        source_time = '{source_time}',
        seed = '{seed}'
    shell:
        '''
        Rscript scripts/mosaic_adm_time_prop.R \
        {input.la_results} \
        {input.model_results} \
        {params.gf_rate} \
        {params.source_num} \
        {params.gen_adm} \
        {params.sample_number} \
        {params.source_time} \
        {params.seed} \
        {output.mosaic_adm_time} \
        {output.mosaic_adm_prop}
        '''

rule rfmix_adm_prop:
    """
    Get admixture proportions from RFMix (averaging across all target samples) 
    """
    input:
        rfmix_out = 'output/seed_{seed}/rfmix/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.rfmix.Q',
    output:
        rfmix_adm_prop = 'output/seed_{seed}/compare_software/files/admixture_proportions/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_adm_prop.txt'
    shell:
        '''
        prop_b=$(cat {input.rfmix_out} | grep "pop_mix" | awk '{{ total += $2 }} END {{ print total/NR }}')
        prop_c=$(cat {input.rfmix_out} | grep "pop_mix" | awk '{{ total += $3 }} END {{ print total/NR }}')

        echo 'rfmix' {wildcards.gf} {wildcards.source} {wildcards.gen_adm} {wildcards.source_time} {wildcards.seed} $prop_b $prop_c > {output.rfmix_adm_prop}
        '''


rule plot_adm_times:
    """ 
    Plot all inferred adm time estimates from FLARE, MOSAIC and DATES
    """
    input:
        flare_adm_time = expand('output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_adm_time.txt', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIMES, allow_missing=True),
        mosaic_adm_time = expand('output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_adm_time.txt', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIMES, allow_missing=True),
        dates_mean = expand('output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_dates_adm_time.txt', seed=SEED, source=SOURCE_POP_SIZE, gen_adm=GEN_ADM, source_time=SOURCE_TIMES, allow_missing=True),
    output:
        all_adm_time_plot = 'output/compare_software/plots/all_simulated_adm_times_gf_{gf}.png',
        all_adm_time_plot_subset = 'output/compare_software/plots/all_simulated_adm_times_subset_gf_{gf}.png',
        all_adm_time = 'output/compare_software/files/all_simulated_adm_times_gf_{gf}.txt',
        all_adm_time_plot_ribbon = 'output/compare_software/plots/all_simulated_adm_times_gf_{gf}_ribbon.png',
        all_adm_time_plot_subset_ribbon = 'output/compare_software/plots/all_simulated_adm_times_subset_gf_{gf}_ribbon.png',
    params:
        geneflow = '{gf}'
    script:
        "../scripts/adm_times_all_plots.R"


rule posterior_probability:
    """ 
    Get posterior probabilities from each LAI method (when applicable)
    """
    input:
        meta_source_name_pop = 'output/seed_{seed}/files/source_name_pop_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        meta_target_name = 'output/seed_{seed}/files/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        flare_txt = 'output/seed_{seed}/flare/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.txt',
        la_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/localanc_pop_mix.RData',
        model_results = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/MOSAIC_RESULTS/pop_mix.RData',
        modern_imputed_snps = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/snpfile.1',
        rfmix_out = 'output/seed_{seed}/rfmix/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}.fb.tsv',
    output:
        all_adm_prop = 'output/seed_{seed}/compare_software/files/posterior_probabilities/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}-all_simulated_post_prob.txt'
    params:
        gf_rate = '{gf}',
        source_num = '{source}',
        gen_adm = '{gen_adm}',
        source_time = '{source_time}',
        seed = '{seed}',
        dir_mosaic = 'output/seed_{seed}/mosaic/output/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}/'
    script:
        "../scripts/posterior_probability_LAI_methods.R"



rule plot_posterior_probability:
    """ 
    Plot posterior probabilities from each LAI method, for a subset of scenarios
    Only for one seed
    """
    input:
        all_post_prob = expand('output/seed_{seed}/compare_software/files/posterior_probabilities/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}-all_simulated_post_prob.txt', seed=[101], source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], allow_missing=True),
    output:
        all_post_prob_plot1 = 'output/compare_software/plots/posterior_prob_gf_{gf}_pop_B_subset.png',
        all_post_prob_plot2 = 'output/compare_software/plots/posterior_prob_gf_{gf}_pop_C_subset.png'
    params:
        gf_rate = '{gf}',
    script:
        "../scripts/plot_posterior_probability_LAI_methods.R"



rule plot_adm_prop_raw:
    """ 
    Get admixture proportions from all methods by just summing up total tract length
    """
    input:
        flare_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        mosaic_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        rfmix_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        simplai_tracts = expand('output/seed_{seed}/compare_software/files/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        flare_tracts_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
        mosaic_tracts_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
        rfmix_tracts_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
        simplai_tracts_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/tracts/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
        flare_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        mosaic_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        rfmix_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        simplai_tracts_overlap = expand('output/seed_{seed}/compare_software/files/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_overlap_stats_prob_{prob}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob=[0.0], allow_missing=True),
        flare_tracts_overlap_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
        mosaic_tracts_overlap_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
        rfmix_tracts_overlap_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_rfmix_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
        simplai_tracts_overlap_pp = expand('output/seed_{seed}/compare_software/files_pp_cutoff/overlap/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_simplai_tracts_overlap_stats_probcutoff_{prob_cutoff}.tsv', seed=SEED, source=[4,10,20,100], gen_adm=GEN_ADM, source_time=[0], prob_cutoff=PROB_CUTOFF, allow_missing=True),
    output:
        adm_prop_raw_popb = 'output/compare_software/plots/adm_prop_raw_{gf}_pop_B_subset.png',
        adm_prop_raw_popc = 'output/compare_software/plots/adm_prop_raw_{gf}_pop_C_subset.png'
    params:
        seq_length = SEQ_LENGTH,
        geneflow = '{gf}'
    script:
        "../scripts/plot_adm_prop_raw_pp.R"


############################################
#  Get adm dates plot for grouped sources  #
############################################

rule flare_adm_time_prop_grouped:
    """
    Get admixture time and proportions from FLARE 
    """
    input:
        flare_model = 'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.model'
    output:
        flare_adm_time = 'output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_adm_time.txt',
        flare_adm_prop = 'output/seed_{seed}/compare_software/files_grouped_sources/admixture_proportions/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_adm_prop.txt'
    shell:
        '''
        date=$(grep -w T: {input.flare_model} -A 1 | sed -n '2p')

        prop=$(grep -w mu {input.flare_model} -A 1 | sed -n '2p')

        echo 'flare' {wildcards.gf} {wildcards.source} {wildcards.gen_adm} 'grouped_sources' {wildcards.seed} $date > {output.flare_adm_time}

        echo 'flare' {wildcards.gf} {wildcards.source} {wildcards.gen_adm} 'grouped_sources' {wildcards.seed} $prop > {output.flare_adm_prop}
        '''

rule mosaic_adm_time_prop_grouped:
    """ 
    Get admixture time and proportions from MOSAIC 
    """
    input:
        la_results = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/MOSAIC_RESULTS/localanc_pop_mix.RData',
        model_results = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/MOSAIC_RESULTS/pop_mix.RData',
    output:
        mosaic_adm_time = 'output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_adm_time.txt',
        mosaic_adm_prop = 'output/seed_{seed}/compare_software/files_grouped_sources/admixture_proportions/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_adm_prop.txt'
    params:
        gf_rate = '{gf}',
        source_num = '{source}',
        gen_adm = '{gen_adm}',
        sample_number = 10,
        source_time = 'grouped_sources',
        seed = '{seed}'
    shell:
        '''
        Rscript scripts/mosaic_adm_time_prop.R \
        {input.la_results} \
        {input.model_results} \
        {params.gf_rate} \
        {params.source_num} \
        {params.gen_adm} \
        {params.sample_number} \
        {params.source_time} \
        {params.seed} \
        {output.mosaic_adm_time} \
        {output.mosaic_adm_prop}
        '''

rule rfmix_adm_prop_grouped:
    """
    Get admixture proportions from RFMix (averaging across all target samples) 
    """
    input:
        rfmix_out = 'output/seed_{seed}/rfmix/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.rfmix.Q',
    output:
        rfmix_adm_prop = 'output/seed_{seed}/compare_software/files_grouped_sources/admixture_proportions/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_rfmix_adm_prop.txt'
    shell:
        '''
        prop_b=$(cat {input.rfmix_out} | grep "pop_mix" | awk '{{ total += $2 }} END {{ print total/NR }}')
        prop_c=$(cat {input.rfmix_out} | grep "pop_mix" | awk '{{ total += $3 }} END {{ print total/NR }}')

        echo 'rfmix' {wildcards.gf} {wildcards.source} {wildcards.gen_adm} 'grouped_sources' {wildcards.seed} $prop_b $prop_c > {output.rfmix_adm_prop}
        '''


rule plot_adm_times_group:
    """ 
    Plot all inferred adm time estimates from FLARE, MOSAIC and DATES
    """
    input:
        flare_adm_time_group = expand('output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_flare_adm_time.txt', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, allow_missing=True),
        mosaic_adm_time_group = expand('output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_mosaic_adm_time.txt', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, allow_missing=True),
        dates_mean_group = expand('output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_dates_adm_time.txt', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, allow_missing=True),
        flare_adm_time = expand('output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_flare_adm_time.txt', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, source_time=[0.2,0.6,1.4,1.8], allow_missing=True),
        mosaic_adm_time = expand('output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_mosaic_adm_time.txt', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, source_time=[0.2,0.6,1.4,1.8], allow_missing=True),
        dates_mean = expand('output/seed_{seed}/compare_software/files/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_source_time_{source_time}_dates_adm_time.txt', seed=SEED, source=[4,8,12,16], gen_adm=GEN_ADM, source_time=[0.2,0.6,1.4,1.8], allow_missing=True),
    output:
        all_adm_time_plot_group = 'output/compare_software/plots_grouped_sources/all_simulated_adm_times_gf_{gf}.png',
        all_adm_time_group = 'output/compare_software/files_grouped_sources/all_simulated_adm_times_gf_{gf}.txt',
        all_adm_time_plot_group_ribbon = 'output/compare_software/plots_grouped_sources/all_simulated_adm_times_gf_{gf}_ribbon.png',
    params:
        geneflow = '{gf}'
    script:
        "../scripts/adm_times_all_plots_grouped.R"


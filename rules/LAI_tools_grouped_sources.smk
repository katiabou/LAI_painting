#########################################
#  Run LAI methods with grouped sources #
#########################################


rule prep_sample_files_grouped_sources:
    """
    Prepare source and target files for flare and mosaic for time series analysis
    Be careful, in rule all, have to specify source=[4,8,12,16], since I'll only focus on these samples sizes
    """
    input:
        sim_vcf_meta = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_meta.tsv',
    output:
        meta_target_name = 'output/seed_{seed}/files_grouped_sources/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_source_name_pop = 'output/seed_{seed}/files_grouped_sources/source_name_pop_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_source_name = 'output/seed_{seed}/files_grouped_sources/source_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_source_name_1 = 'output/seed_{seed}/files_grouped_sources/source_name_1_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_source_name_2 = 'output/seed_{seed}/files_grouped_sources/source_name_2_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_target_source_pop_name = 'output/seed_{seed}/files_grouped_sources/target_source_pop_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_target_source_name = 'output/seed_{seed}/files_grouped_sources/target_source_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
    params:
        source_size = '{source}',
        gen_adm = '{gen_adm}',
        source_1 = 'pop_b',
        source_2 = 'pop_c',
        target = 'pop_mix',
    script:
        "../scripts/prep_sample_name_files_grouped_sources.R"


###############
#  Run flare  #
###############

rule flare_simulation_grouped_sources:
    """ 
    Run flare on simulated data for grouped model
    Have to specify --min-mac value lower than number of reference samples. Setting it to 1 to include all sites. 
    Leaving min-maf to default (0.005)
    """
    input:
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
        meta_target_name = 'output/seed_{seed}/files_grouped_sources/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        meta_source_name_pop = 'output/seed_{seed}/files_grouped_sources/source_name_pop_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        gen_map_flare = 'output/seed_{seed}/flare/files/model_gf_{gf}_gen_adm_{gen_adm}-genetic_map.txt',
    output:
        flare_out = 'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.anc.vcf.gz',
        flare_model = 'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.model'
    params:
        flare_path = config['flare_path'],
        prefix_out = 'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources',
        gen = '{gen_adm}'
    threads: 4
    log:
        'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.anc.vcf.gz.log'
    benchmark:
        'benchmarks/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.tsv'
    shell:
        '''
        java -Xmx500g -jar {params.flare_path} \
        ref={input.sim_vcf_filt} \
        ref-panel={input.meta_source_name_pop} \
        gt={input.sim_vcf_filt} \
        map={input.gen_map_flare} \
        gt-samples={input.meta_target_name} \
        gen={params.gen} \
        out={params.prefix_out} \
        min-mac=1 \
        probs=true \
        nthreads={threads}

        bcftools index -f {output.flare_out}
        '''

rule prep_file_flare_sim_grouped_sources:
    """ 
    Prepare flare file format for plotting
    """
    input:
        flare_out = 'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.anc.vcf.gz',
    output:
        flare_txt = 'output/seed_{seed}/flare/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
    shell:
        '''
        bcftools query -f '%CHROM\t%POS\t%FORMAT\n' \
        {input.flare_out} > {output.flare_txt}
        '''



################
#  Run mosaic  #
################

rule pop_source1_vcf_grouped_sources:
    """ 
    Extract selected source 1 population samples from vcf
    Make genotype file for MOSAIC
    """
    input:
        meta_source_name_1 = 'output/seed_{seed}/files_grouped_sources/source_name_1_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        pop_b_vcf = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_b.vcf.gz',
        genotype_file_pop_b = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_bgenofile.1',
    threads: 4
    shell:
        '''        
        bcftools view \
        {input.sim_vcf_filt} \
        -S {input.meta_source_name_1} \
        -m 2 -M 2 \
        --threads {threads} \
        -Oz -o {output.pop_b_vcf}

        bcftools index -f {output.pop_b_vcf}

        bcftools query -f '[%GT]\n' {output.pop_b_vcf} | sed 's/|//g' > {output.genotype_file_pop_b}
        '''

rule pop_source2_vcf_grouped_sources:
    """ 
    Extract selected source 2 population samples from vcf
    Make genotype file for MOSAIC
    """
    input:
        meta_source_name_2 = 'output/seed_{seed}/files_grouped_sources/source_name_2_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        pop_c_vcf = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_c.vcf.gz',
        genotype_file_pop_c = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_cgenofile.1',
    threads: 4
    shell:
        '''
        bcftools view \
        {input.sim_vcf_filt} \
        -S {input.meta_source_name_2} \
        -m 2 -M 2 \
        --threads {threads} \
        -Oz -o {output.pop_c_vcf}

        bcftools index -f {output.pop_c_vcf}

        bcftools query -f '[%GT]\n' {output.pop_c_vcf} | sed 's/|//g' > {output.genotype_file_pop_c}
        '''

rule pop_target_vcf_grouped_sources:
    """ 
    Extract selected target population samples from vcf
    Make genotype file for MOSAIC
    """
    input:
        meta_target_name = 'output/seed_{seed}/files_grouped_sources/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        pop_mix_vcf = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mix.vcf.gz',
        genotype_file_pop_mix = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mixgenofile.1',
    threads: 4
    shell:
        '''
        bcftools view \
        {input.sim_vcf_filt} \
        -S {input.meta_target_name} \
        -m 2 -M 2 \
        --threads {threads} \
        -Oz -o {output.pop_mix_vcf}

        bcftools index -f {output.pop_mix_vcf}

        bcftools query -f '[%GT]\n' {output.pop_mix_vcf} | sed 's/|//g' > {output.genotype_file_pop_mix}
        '''

rule snp_list_mosaic_grouped_sources:
    """ 
    Get list of snps with ref and alt 
    """
    input:
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        modern_imputed_snps = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/snpfile.1',
    shell:
        '''
        bcftools query -f '%CHROM\_%POS %CHROM 0 %POS %REF %ALT\n' {input.sim_vcf_filt} | sed 's/chr//g' > {output.modern_imputed_snps}
        '''

rule sample_list_mosaic_grouped_sources:
    """ 
    Copy list of samples into mosaic running directory (it's picky about where things are)
    """
    input:
        meta_target_source_pop_name = 'output/seed_{seed}/files_grouped_sources/target_source_pop_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
    output:
        samples_mosaic = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/sample.names',
    shell:
        '''
        cp {input.meta_target_source_pop_name} {output.samples_mosaic}
        '''

rule rate_files_mosaic_grouped_sources:
    """ 
    Transpose rate files for mosaic format, using number of sites and actual sites in snpfile
    """
    input:
        modern_imputed_snps = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/snpfile.1',
    output:
        temp_sites = temp('output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/temp_sites.1'),
        temp_rate = temp('output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/temp_rate.1'),
        rate_fie = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/rates.1',
    shell:
        '''
        SITES=$(wc -l < {input.modern_imputed_snps})
        echo ':sites:'$SITES > {output.temp_sites}

        cut -d$' ' -f 4  {input.modern_imputed_snps} | awk '{{print $1" "$1/1000000}}' | \
        awk '
        {{ 
            for (i=1; i<=NF; i++)  {{
                a[NR,i] = $i
            }}
        }}
        NF>p {{ p = NF }}
        END {{    
            for(j=1; j<=p; j++) {{
                str=a[1,j]
                for(i=2; i<=NR; i++){{
                    str=str" "a[i,j];
                }}
                print str
            }}
        }}' > {output.temp_rate}

        cat {output.temp_sites} {output.temp_rate} >> {output.rate_fie}
        '''


rule mosaic_grouped_sources:
    """ 
    Run MOSAIC, using the nophase option to not allow it to flip my phasing
    HAVE TO BE CAREFUL WITH THE NUMER OF SAMPLES SPECIFIED IN --number (equal to the number of targets sampled)
    Specifying the number of grid points based on Browning et al., 2023
    """
    input:
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
        genotype_file_pop_b = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_bgenofile.1',
        genotype_file_pop_c = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_cgenofile.1',
        genotype_file_pop_mix = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mixgenofile.1',
        modern_imputed_snps = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/snpfile.1',
        samples_mosaic = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/sample.names',
        rate_fie = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/rates.1',
    output:
        la_results = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/MOSAIC_RESULTS/localanc_pop_mix.RData',
        model_results = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/MOSAIC_RESULTS/pop_mix.RData',
    params:
        path_mosaic = config['mosaic_path'],
        dir_mosaic = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources',
        tmp_dir = '/tmp/seed_{seed}_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources',
        output_dir = 'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/MOSAIC_RESULTS',
        target = 'pop_mix',
        gens = '{gen_adm}'
    benchmark:
        'benchmarks/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.tsv'
    log:
        'output/seed_{seed}/mosaic/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/run_mosaic.log'
    threads: 6
    shell:
        '''
        SITES=$(bcftools view -H {input.sim_vcf_filt} | wc -l)
        GRID=$(echo "$SITES * 0.0012" | bc) 

        # make temporary folder
        mkdir {params.tmp_dir}

        Rscript {params.path_mosaic} \
        {params.target} \
        {params.dir_mosaic}/ \
        --number 10 \
        --ancestries 2 \
        --chromosomes 1:1 \
        --GpcM $GRID \
        --gens {params.gens} \
        --nophase \
        --output_dir {params.output_dir} \
        --fastfiles {params.tmp_dir} \
        --maxcores {threads} 2> {log}

        # try to fix the naming mess
		cd {params.output_dir}
		find localanc_{params.target}_*.RData -exec cp {{}} localanc_{params.target}.RData \;
		find {params.target}_*.RData -exec cp {{}} {params.target}.RData \;

        # remove tmp directory
        rm -r {params.tmp_dir}
        '''



###############
#  Run rfmix  #
###############

rule prep_target_vcf_rfmix_grouped_sources:
    """ 
    Extract target samples from simulated VCF
    """
    input:
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
        meta_target_name = 'output/seed_{seed}/files_grouped_sources/target_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
    output:
        target_vcf = 'output/seed_{seed}/rfmix/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_filt-targets.vcf.gz',
    threads: 4   
    shell:
        '''
        bcftools view \
        {input.sim_vcf_filt} \
        -S {input.meta_target_name} \
        -m 2 -M 2 \
        --threads {threads} \
        -Oz -o {output.target_vcf}
        
        bcftools index -f {output.target_vcf}
        '''

rule prep_source_vcf_rfmix_grouped_sources:
    """ 
    Extract source samples from simulated VCF
    """
    input:
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
        meta_source_name = 'output/seed_{seed}/files_grouped_sources/source_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
    output:
        source_vcf = 'output/seed_{seed}/rfmix/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_filt-source.vcf.gz',
    threads: 4   
    shell:
        '''
        bcftools view \
        {input.sim_vcf_filt} \
        -S {input.meta_source_name} \
        -m 2 -M 2 \
        --threads {threads} \
        -Oz -o {output.source_vcf}
        
        bcftools index -f {output.source_vcf}
        '''

rule run_rfmix_grouped_sources:
    """
    Look here: https://github.com/slowkoni/rfmix/blob/master/MANUAL.md 
    Have to check when I should use the --reanalyze-reference option. In the case a set of reference haplotypes may not be of "pure" ancestry and may themselves be somewhat admixed
    """
    input:
        gen_map_rfmix = 'output/seed_{seed}/rfmix/files/model_gf_{gf}_gen_adm_{gen_adm}-genetic_map.txt',
        target_vcf = 'output/seed_{seed}/rfmix/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_filt-targets.vcf.gz',
        source_vcf = 'output/seed_{seed}/rfmix/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_filt-source.vcf.gz',
        meta_source_name_pop = 'output/seed_{seed}/files_grouped_sources/source_name_pop_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
    output:
        rfmix_out = 'output/seed_{seed}/rfmix/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.fb.tsv',
    params:
        rfmix_path = config['rfmix_path'],
        prefix = 'output/seed_{seed}/rfmix/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources',
        rfmix_iterations = config['rfmix_iterations'],
        gen = '{gen_adm}'
    benchmark:
        'benchmarks/seed_{seed}/rfmix/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.tsv'
    threads: 6
    shell:
        '''
        {params.rfmix_path} \
        -f {input.target_vcf} \
        -r {input.source_vcf} \
        -m {input.meta_source_name_pop} \
        -g {input.gen_map_rfmix} \
        -o {params.prefix} \
        --chromosome=chr1 \
        --n-threads={threads} \
        -e {params.rfmix_iterations} \
        -G {params.gen}
        '''


#################
#  Run simplai  #
#################

rule extract_simplai_samples_grouped_sources:
    """
    Extract samples from VCF for simplai
    """
    input:
        meta_target_source_name = 'output/seed_{seed}/files_grouped_sources/target_source_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        simplai_vcf = 'output/seed_{seed}/simplai/files_grouped_sources/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.vcf.gz',
        simplai_gen = 'output/seed_{seed}/simplai/output_grouped_sources/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.gen',
    shell:
        '''
        bcftools view -S {input.meta_target_source_name} \
        {input.sim_vcf_filt} \
        -Oz -o {output.simplai_vcf}  

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {output.simplai_vcf} | sed "s/|/\t/g" | sed "s/chr//g" > {output.simplai_gen}
        '''

rule run_simplai_grouped_sources:
    """
    Be careful with the number of targets (ssa) and the name of the output file, has to match the chosen -n and -m parameters
    """
    input:
        simplai_gen = 'output/seed_{seed}/simplai/output_grouped_sources/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.gen',
    output:
        simplai_output_rec = 'output/seed_{seed}/simplai/output_grouped_sources/samples_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_n_2000_m_1000_t_5_fromRec_withSingl.adm'
    params:
        seq_length = SEQ_LENGTH,
        simplai_path = config['simplai_path'],
        target_pop_size = config['target_pop_size']
    benchmark:
        'benchmarks/seed_{seed}/simplai/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.tsv'
    shell:
        '''
        SOURCE={wildcards.source}
        TARGET={params.target_pop_size}
        DIPLOID=2
        HAPL=$((SOURCE * DIPLOID))
        TARGET_HAPL=$((TARGET * DIPLOID))

        {params.simplai_path} \
        -g {input.simplai_gen} \
        --ss1 $HAPL \
        --ss2 $HAPL \
        --ssa $TARGET_HAPL \
        -l {params.seq_length} \
        -s 1e6 \
        -i 500000 \
        -n 2000 \
        -m 1000 \
        -t 5
        '''


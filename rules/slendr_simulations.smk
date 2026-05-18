#######################################
#  Run slendr to simulate null models #
#######################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']

rule simulate_data:
    """
    Simulate genetic data with slendr
    """
    output:
        plot = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}.png',
        sim_vcf = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}.vcf.gz',
        sim_vcf_meta = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_meta.tsv',
        tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv',
        fst = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_fst.tsv',
    params:
        gf_rate='{gf}',
        gen_adm='{gen_adm}',
        genome_size_bp = SEQ_LENGTH,
        seed = '{seed}'
    script:
        "../scripts/slendr_null_size_new_model.R"


rule filter_simulated_vcf:
    """
    Filter for biallelic sites only and apply a MAF filter of 1%
    """
    input:
        sim_vcf = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}.vcf.gz',
    output:
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    shell:
        '''
        bcftools view \
        -m 2 -M 2 \
        -q 0.01:minor \
        {input.sim_vcf} \
        -Oz -o {output.sim_vcf_filt}

        bcftools index -f {output.sim_vcf_filt}
        '''

# rule add_info_tract_file:
#     """
#     Add to tract file details for each run
#     """
#     input:
#         tracts = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv'
#     output:
#         tracts_info = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts_info.tsv'
#     shell:
#         '''
#         awk -F'\t' '{{$(NF+1)={wildcards.gen_adm} FS {wildcards.gf} FS {wildcards.seed};}}1' OFS='\t' {input.tracts} > {output.tracts_info}
#         '''

rule tract_lengths:
    """
    Plot histogram of tract lengths
    """
    input:
        tracts = expand('output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_tracts.tsv', seed=SEED, gen_adm=GEN_ADM, allow_missing=True)
    output:
        tracts_hist = 'output/track_length_slendr/tracts_gf_{gf}_histogram_length.png',
        tracts_average = 'output/track_length_slendr/tracts_gf_{gf}_average_length.tsv',
    params:
        tracts_files=lambda wildcards, input: ','.join(input.tracts),
    shell:
        '''
        Rscript scripts/tract_length.R \
        {params.tracts_files} \
        {output.tracts_hist} \
        {output.tracts_average}
        '''

rule fst_plots:
    """
    Plot fst between sources and targets
    """
    input:
        fst = expand('output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_fst.tsv', seed=SEED, gen_adm=GEN_ADM, allow_missing=True)
    output:
        fst_all = 'output/fst_source_target/fst_gf_{gf}_all.png',
        fst_sources = 'output/fst_source_target/fst_gf_{gf}_sources.png',
    params:
        fst_files=lambda wildcards, input: ','.join(input.fst),
    shell:
        '''
        Rscript scripts/fst_plots.R \
        {params.fst_files} \
        {output.fst_all} \
        {output.fst_sources}
        '''

rule pca_source_target:
    """
    Create plink files
    """
    input:
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        expand('output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt'
    shell:
        '''
        plink \
        --vcf {input.sim_vcf_filt} \
        --make-bed \
        --biallelic-only strict \
        --allow-no-sex \
        --set-missing-var-ids @:#[b37] \
        --out {params.prefix} \
        --double-id
        '''

rule pcaone_source_target:
    """
    Run PCAone  
    """
    input:
        bed = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.bed',
        bim = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.bim',
        fam = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.fam',
    output:
        pcaone_eigenvals = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.eigvals',
        pcaone_eigenvecs = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.eigvecs'
    params:
        prefix = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt',
        prefix_output = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt'
    log:
        'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.eig.log'
    threads: 8
    shell:
        '''
        PCAone \
        --bfile {params.prefix} \
        --svd 2 \
        --maf 0.01 \
        -k 10 \
        --scale 0 \
        --threads {threads} \
        -o {params.prefix_output} 2> {log}
        '''

rule plot_pcaone_source_target:
    """
    Plot PCAone  
    """
    input:
        sim_vcf_meta = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_meta.tsv',
        pcaone_eigenvecs = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.eigvecs',
        pcaone_eigenvals = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt.eigvals',
    output:
        pcaone_p1p2 = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt-PC1_PC2.png',
        pcaone_p2p3 = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt-PC2_PC3.png',
        pcaone_p1p4 = 'output/seed_{seed}/slendr/pca/model_gf_{gf}_gen_adm_{gen_adm}_filt-PC1_PC4.png',
    script:
        "../scripts/pcaone_plot.R"

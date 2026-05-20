################
#  Run DATES   #
################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']

rule make_plink_dates_grouped:
    """
    Make plink files from simulated VCF
    """
    input:
        meta_target_source_pop_name = 'output/seed_{seed}/files_grouped_sources/target_source_pop_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        meta_target_source_pop_name_no_fid = temp('output/seed_{seed}/dates/files_grouped_sources/target_source_pop_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-no_fid.txt'),
        dates_plink = expand('output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources'
    log: 
        'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_plink.log'
    shell:
        '''
        awk -F"\t" '{{print $2 "\t" $2}}' {input.meta_target_source_pop_name} > {output.meta_target_source_pop_name_no_fid}

        plink \
        --vcf {input.sim_vcf_filt} \
        --keep {output.meta_target_source_pop_name_no_fid} \
        --make-bed \
        --out {params.prefix} \
        --double-id 2> {log}
        '''

rule prepare_convertf_dates_grouped:
    """
    Prepare for convertf
    """
    input:
        dates_plink_bim = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.bim',
        dates_plink_bed = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.bed',
        dates_plink_fam = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.fam',
        sim_vcf_filt = 'output/seed_{seed}/slendr/sim_data/model_gf_{gf}_gen_adm_{gen_adm}_filt.vcf.gz',
    output:
        dates_plink_bim_mod = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-mod.bim',
        dates_plink_fam_mod = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-mod.fam',
        convertf_parfile = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-convertf_parfile'
    params:
        prefix = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources'
    shell:
        '''
        awk 'BEGIN{{OFS="\t"}}$2=$1"_"$4' {input.dates_plink_bim} > {output.dates_plink_bim_mod}
        awk '{{$6=2 ; print ; }}' {input.dates_plink_fam} > {output.dates_plink_fam_mod}

        echo "                      
        genotypename: {input.dates_plink_bed}
        snpname:      {output.dates_plink_bim_mod}
        indivname:    {output.dates_plink_fam_mod}
        outputformat:    EIGENSTRAT                                     
        genooutfilename:   {params.prefix}.eigenstratgeno
        snpoutfilename:    {params.prefix}.snp
        indoutfilename:    {params.prefix}.ind
        familynames:       NO
        " >> {output.convertf_parfile}
        '''

rule run_convertf_dates_grouped:
    """
    Run convertf
    """
    input:
        convertf_parfile = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-convertf_parfile'
    output:
        dates_eigenstrat = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.eigenstratgeno',
        dates_snp = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.snp',
        dates_ind = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.ind',
    log: 
        'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_convertf.log'
    shell:
        '''
        convertf -p {input.convertf_parfile} 2> {log}
        '''

rule prepare_files_dates_grouped:
    """
    Prepare files for dates, adding populations to ind file and create admix file
    """
    input:
        dates_eigenstrat = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.eigenstratgeno',
        dates_snp = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.snp',
        dates_ind = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.ind',
        meta_target_source_pop_name = 'output/seed_{seed}/files_grouped_sources/target_source_pop_name_model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.txt',
    output:
        temp1 = temp('output/seed_{seed}/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-temp1.txt'),
        dates_ind_mod = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-mod.ind',
        admix = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-admix',
    params:
        prefix = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources'
    shell:
        '''
        awk '{{print $2 "\t" $1}}' {input.meta_target_source_pop_name} > {output.temp1}
        awk 'NR==FNR {{map[$1]=$2; next}} {{if ($1 in map) $3=map[$1]; print}}' {output.temp1} {input.dates_ind} > {output.dates_ind_mod}

        echo "pop_b   pop_c   pop_mix {params.prefix}" > {output.admix}
        '''

rule parfile_dates_grouped:
    """
    Prepare DATES parfile
    """
    input:
        dates_eigenstrat = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.eigenstratgeno',
        dates_snp = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources.snp',
        dates_ind_mod = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-mod.ind',
        admix = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-admix',
    output:
        dates_parfile = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-par.dates',
    params:
        prefix = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources'
    shell:
        '''
        echo "
        genotypename: {input.dates_eigenstrat}
        snpname:    {input.dates_snp}
        indivname:  {input.dates_ind_mod}
        admixlist:  {input.admix}   # This file contains the source and admixed populations to use for the analysis. Each line has the format: <source1> <source2> <testpopulation> <output_directory>, where source1 and source2 are the reference populations for the ancestral populations, testpop is the name of the admixed population and the output_directory is the name of the output directory.Output files are of the format output_directory/testpopulation.out
        binsize:    0.001   # in Morgans, range is from 0-1. Optimal binsize of 0.001 is recommended.
        maxdis:     1.0     # in Morgans, range is 0-1. For quicker runs, use max_distance < 1.0. However, for recent admixture,   ensure that max_distance is greater than the expected admixture LD blocks.
        seed:       77      # Random seed to ensure reproducibility of runs. 
        jackknife: NO      # if YES, program will run jackknife by dropping one chr in each run and estimate the mean and standard error across the 23 runs.
        qbin:       10      # discretization parameter on mesh size for the binned residuals. Higher qbin correlates loosely with higher accuracy and highly with longer run time.
        runfit:  YES        # run exponential fit using least squares on the output to infer the date of admixture?
        afffit:    YES      # use affine for the fit? 
        lovalfit:  0.45     # in centiMorgans, starting genetic distance.
        checkmap: NO        # ONLY TO USE FOR SIMULATIONS
        " >> {output.dates_parfile}
        '''

rule run_dates_grouped:
    """
    Run DATES
    Added the || true part at the end, since it was giving a segmentation fault error saying that -r parameter unknown (even though it is)
    """
    input:
        dates_parfile = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-par.dates',
    output:
        dates_log = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mix.out',
        dates_mix_log = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mix:log',
        dates_mix_expfit = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mix:expfit.out',
    shell:
        '''
        export PATH="/projects/racimolab/people/qcj125/programmes/dates/bin:$PATH"

        dates -p /maps/projects/racimolab/people/qcj125/MUNICH_FILES/MESO_DOGS/ANALYSES/SIMULATIONS_PAINTING_TEMP7/{input.dates_parfile} || true
        '''

rule get_estimate_dates_grouped:
    """
    Get mean estimate from DATES 
    """
    input:
        dates_mix_log = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mix:log'
    output:
        temp1 = temp('output/seed_{seed}/dates/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-temp1.txt'),
        temp2 = temp('output/seed_{seed}/dates/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-temp2.txt'),
        dates_mean = 'output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_dates_adm_time.txt'
    shell:
        '''
        echo 'DATES' {wildcards.gf} {wildcards.source} {wildcards.gen_adm} 'grouped_sources' {wildcards.seed} > {output.temp1}
        cat {input.dates_mix_log} | grep "mean (generations):" | awk -F" " '{{print $3}}' > {output.temp2}
        paste -d ' ' {output.temp1} {output.temp2} > {output.dates_mean}
        '''


rule get_nrmsd_dates_grouped:
    """
    Get nrmsd estimate from DATES 
    """
    input:
        dates_mix_log = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mix:log'
    output:
        temp1 = temp('output/seed_{seed}/dates/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-temp1_nrmsd.txt'),
        temp2 = temp('output/seed_{seed}/dates/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-temp2_nrmsd.txt'),
        dates_nrmsd = 'output/seed_{seed}/dates/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-nrmsd_dates.txt'
    shell:
        '''
        echo {wildcards.seed} {wildcards.gf} {wildcards.source} {wildcards.gen_adm} 'grouped_sources' > {output.temp1}
        cat {input.dates_mix_log} | grep "nrmsd" | awk -F" " '{{print $2}}' > {output.temp2}
        paste -d ' ' {output.temp1} {output.temp2} > {output.dates_nrmsd}
        '''

rule plot_curves_each_grouped:
    """
    Plot each run with it's own curve (this is normally produced by the jackknife step of DATES, but since we didn't use it we have to make it ourselves)
    """
    input:
        dates_mix_expfit = 'output/seed_{seed}/dates/files_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources/pop_mix:expfit.out',
        dates_mean = 'output/seed_{seed}/compare_software/files_grouped_sources/admixture_times/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_dates_adm_time.txt',
        dates_nrmsd = 'output/seed_{seed}/dates/output_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources-nrmsd_dates.txt'
    output:
        dates_mix_expfit_plot = 'output/seed_{seed}/dates/plots_grouped_sources/model_gf_{gf}_source_{source}_gen_adm_{gen_adm}_grouped_sources_curve.png',
    params:
        gf_rate='{gf}',
        gen_adm='{gen_adm}',
        source='{source}',
        source_time='grouped_sources',
        seed = '{seed}'
    script:
        "../scripts/dates_curves.R"

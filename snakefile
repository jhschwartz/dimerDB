from snakemake.utils import min_version
min_version('7.12')


subworkflow download_pdb:
    workdir:
        './'
    snakefile:
        'subworkflows/snakefile0_download_pdb'


subworkflow retrieve_all_homodimers:
    workdir:
        './'
    snakefile:
        'subworkflows/snakefile1_retrieve_all'


#subworkflow prune_homodimers:
#    workdir:
#        './'
#    snakefile:
#        'subworkflows/snakefile2_prune_dimers'
#
#
#subworkflow prune_seqs:
#    workdir:
#        './'
#    snakefile:
#        'subworkflows/snakefile3_prune_seqs'

rule all:
    input:
        download_pdb('intermediates/subflow0.done'),
        retrieve_all_homodimers('intermediates/subflow1.done')
        #prune_homodimers('intermediates/subflow2.done'),
        #prune_seqs('intermediates/subflow4.done'),
        #'intermediates/cleanup.done'


#rule cleanup:
#    input:
#        done = 'intermediates/cleanup.done'
#    shell:
#        '''
#        rm -rf intermediates/all_homodimer_tmp/;
#        touch {input.done};
#        '''

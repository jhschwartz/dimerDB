from snakemake.utils import min_version
min_version('7.12')

configfile: 'config.yaml'

subworkflow download_pdb:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile0']


subworkflow deduce_all_dimers:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile1']



subworkflow filter_homodimers:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile2.1']



subworkflow prune_homo_structures:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile3.1']


subworkflow prune_homo_seqs:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile4.1']


subworkflow filter_homodimers:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile2.2']


subworkflow prune_hetero_structures:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile3.2']


subworkflow prune_hetero_seqs:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile4.2']


rule all:
    input:
        download_pdb(config['snake_donefiles']['sub0_all_done']),
        deduce_all_dimers(config['snake_donefiles']['sub1_all_done']),
        filter_homodimers(config['snake_donefiles']['sub2_all_done'])
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

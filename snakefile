from snakemake.utils import min_version
min_version('7.12')

configfile: 'config.yaml'


import os
basepath = config['paths']['basepath']
if basepath != os.getcwd():
    raise RuntimeError('config does not match the current directory!')


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



subworkflow filter_prune_homodimers:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile2']


subworkflow prune_seqs:
    workdir:
        config['paths']['basepath']
    snakefile:
        config['workflow']['snakefile3']



rule all:
    input:
        download_pdb(config['snake_donefiles']['sub0_all_done']),
        deduce_all_dimers(config['snake_donefiles']['sub1_all_done'])
#        filter_prune_homodimers(config['snake_donefiles']['sub2_all_done']),
#        prune_seqs(config['snake_donefiles']['sub3_all_done']),
        #'intermediates/cleanup.done'


#rule cleanup:
#    input:
#        done = 'intermediates/cleanup.done'
#    shell:
#        '''
#        rm -rf intermediates/all_homodimer_tmp/;
#        touch {input.done};
#        '''

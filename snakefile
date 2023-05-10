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
        os.path.realpath('./new-subworkflows/0_local_pdb.smk')


rule all:
    input:
        download_pdb(config['workflow']['0_local_pdb_done'])

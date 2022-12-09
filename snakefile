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



rule all:
    input:
        download_pdb('intermediates/subflow0.done'),
        retrieve_all_homodimers('intermediates/subflow1.done')


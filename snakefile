from snakemake.utils import min_version
min_version('7.12')

configfile: 'config.yaml'

import os
basepath = config['paths']['basepath']
if basepath != os.getcwd():
    raise RuntimeError('config does not match the current directory!')


from datetime import datetime
start = datetime.utcnow()


# 0 -> output: intermediates/sub0.done
out_sub0 = config['subworkflow_done']['0_local_pdb']
subworkflow download_pdb:
    workdir:
        config['paths']['basepath']
    snakefile:
        os.path.realpath('workflow/0_local_pdb.smk')


# 1 -> output: intermediates/sub1.done
out_sub1 = config['subworkflow_done']['1_derive_dimers']
subworkflow derive_dimers:
    workdir:
        config['paths']['basepath']
    snakefile:
        os.path.realpath('workflow/1_derive_dimers.smk')


# 2 -> output: intermediates/sub2.done
out_sub2 = config['subworkflow_done']['2_seq_cluster']
subworkflow seq_cluster:
    workdir:
        config['paths']['basepath']
    snakefile:
        os.path.realpath('workflow/2_seq_cluster.smk')


# 3 -> output: intermediates/sub3.done
out_sub3 = config['subworkflow_done']['3_rm_structural_redundancy']
subworkflow rm_structural_redundancy:
    workdir:
        config['paths']['basepath']
    snakefile:
        os.path.realpath('workflow/3_rm_structural_redundancy.smk')



# 4 -> output: intermediates/sub4.done
out_sub4 = config['subworkflow_done']['4_prep_final_data']
subworkflow prep_final_data:
    workdir:
        config['paths']['basepath']
    snakefile:
        os.path.realpath('workflow/4_prep_final_data.smk')


localrules: all
rule all:
    input:
        download_pdb(               out_sub0 ),
        derive_dimers(              out_sub1 ),
        seq_cluster(                out_sub2 ),
        rm_structural_redundancy(   out_sub3 ),
        prep_final_data(            out_sub4 )
    run:
        print('PIPELINE FINISHED!')
        print('full pipeline start time:')
        shell(''' cat build_start.tmp ''')
        print('full pipeline start time:')
        print('full pipeline end time:')
        print(datetime.utcnow())
        
        # each outfile timestamps the end of each subworkflow, so we print
        shell('''
            cat {out_sub0};
            cat {out_sub1};
            cat {out_sub2};
            cat {out_sub3};
            cat {out_sub4};
              ''')



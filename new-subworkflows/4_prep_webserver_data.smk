import os
import sys

configfile: 'config.yaml'


sys.path.append('scripts')
import name_pdb


subworkflow_done = config['workflow']['4_prep_webserver_data_done']

files = config['paths']['pipeline_files']
webserver_data = config['paths']['webserver_data']

pdb_lib = config['paths']['pdb_lib'] # TODO put in config, match to sub0
chains_index = config['paths']['chains_index'] # TODO put in config, match to sub0




pdb_chains_tarfile_template = webserver_data+'/pdb/div/{div}.tar.gz'



### retrieve list of divs used in assemblies_index ###
###    and use for rule tar_pdb_chain_lib below    ###
chain_divs = set()
with open(chains_index, 'r') as f:
    for line in f:
        chain_name = line.strip()
        div = name_pdb.get_div(chain_name)
        chain_divs.add(div)

rule tar_pdb_chain_lib:
    output:
#        done = tar_chainfiles_done
    threads: 1
    resources:
        time = '12:00:00',
        memory = '10000'
    run:
        tar_pdb_chains(indexfile=input.index,
                       chains_lib=assemblies_lib,
                       target_dir=tar_pdb_parent_dir
                       )
        shell(''' touch {output.done} ''')







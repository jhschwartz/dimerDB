import os
import sys 
import pickle
import tarfile
from datetime import datetime

configfile: 'config.yaml'

sys.path.append('scripts')
from parallel_convert_split_cif import parallel_convert_split_rename_cifs, fill_empty_pdb
from tar_pdb_chains import tar_pdb_chains
from generate_rcsb_index import generate_rcsb_index

subworkflow_done = config['workflow']['0_local_pdb_done']

files = config['paths']['intermediates_dir']
tar_pdb_parent_dir = config['paths']['pdb_tar_lib']

assemblies_lib = f'{files}/sub0/no_sample/assemblies' # TODO move to config 
assemblies_download_done = f'{files}/sub0/no_sample/assemblies_download.done'
assemblies_ungzip_done = f'{files}/sub0/no_sample/assemblies_ungzip.done'
split_renamed_all_chains_done = f'{files}/sub0/no_sample/split+renamed.done'
tar_chainfiles_done = f'{files}/sub0/no_sample/chains_pdb_tar.done'
chains_index_file = f'{assemblies_lib}/index_pdb.txt' # TODO move to config for later use
empty_pdbs_list = f'{files}/sub0/no_sample/initially_empty_chain_pdbs'
empty_files_filled_done = f'{files}/sub0/no_sample/filled_init_empty_pdbs.done'


localrules: all
rule all:
    input: 
        done = empty_files_filled_done
    output:
        done = subworkflow_done
    run:
        with open(output.done, 'w') as f:
            f.write('Subworkflow 0 done at ')
            f.write(str(datetime.utcnow()))



localrules: download_assemblies
rule download_assemblies:
    output:
        done = assemblies_download_done
    shell:
        ''' 
        rsync -rlpt -v -z -q --delete --port=33444 \
        rsync.rcsb.org::ftp_data/assemblies/mmCIF/divided/ {assemblies_lib};
        touch {output.done};
        '''
        
        # The above should be unaffected by the deprecation
        # of ftp protocol by rcsb; while it is "ftp_data",
        # the protocol itself is rsync which will remain.




rule ungzip_assemblies:
    input:
        done = assemblies_download_done
    output:
        done = assemblies_ungzip_done
    threads: 32 
    resources:
        time = '12:00:00',
        mem_mb = '500000'
    shell:
        '''
        scripts/parallel_ungzip_all.sh {assemblies_lib} {threads} 1000;
        touch {output.done};
        '''
        






rule split_rename_chains:
    input:
        done = assemblies_ungzip_done
    output:
        done = split_renamed_all_chains_done
    threads: 32 
    resources:
       mem_mb = '500000',
       time = '24:00:00'
    run:
        exe = os.path.join('bin', 'USalign', 'cif2pdb')
        parallel_convert_split_rename_cifs(parent_dir=assemblies_lib,
                                           cif2pdb_exe=exe,
                                           threads=threads
                                           )
        shell(''' touch {output.done} ''')






rule index_chains_pdb:
    input:
        done = split_renamed_all_chains_done
    output:
        indexfile = chains_index_file,
        empty_files_txt = empty_pdbs_list
    threads: 1
    resources:
        time = '4:00:00'
    run:
        empty_files = generate_rcsb_index(rcsb_path=assemblies_lib, index_file=output.indexfile)
        with open(output.empty_files_txt, 'w') as f:
            for ef in empty_files:
                f.write(f'{ef}\n')




rule fill_empty_pdb_files:
    input:
        empty_files_txt = empty_pdbs_list
    output:
        done = empty_files_filled_done
    threads: 1
    resources:
        mem_mb = '15000'
    run:
        exe = os.path.join('bin', 'USalign', 'cif2pdb')
        with open(input.empty_files_txt, 'rb') as f:
            for line in f:
                empty_file = line.strip()
                fill_empty_pdb(filepath=empty_file, cif2pdb_exe=exe) 
        shell(''' touch {output.done} ''')
            


















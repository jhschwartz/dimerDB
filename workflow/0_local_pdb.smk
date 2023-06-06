import os
import sys 
import pickle
import tarfile
from datetime import datetime

configfile: 'config.yaml'

sys.path.append('scripts')
from parallel_convert_split_cif import parallel_convert_split_rename_cifs, fill_empty_pdb
from generate_rcsb_index import generate_rcsb_index

subworkflow_done = config['subworkflow_done']['0_local_pdb']

files = config['paths']['intermediates_dir']
lib = config['paths']['lib']

cif2pdb_exe = config['exe']['cif2pdb']

outfile = {
    'assemblies_dir': os.path.join(lib, 'rcsb'),
    'chains_index': os.path.join(lib, 'pdb_index.txt'),
    'empty_pdbs': os.path.join(files, 'sub_0', 'initially_empty_pdbs.txt'),
    'done': {
        'rsync': os.path.join(files, 'sub_0', 'assemblies_download.done'),
        'ungzip': os.path.join(files, 'sub_0', 'assemblies_ungzip.done'),
        'split': os.path.join(files, 'sub_0', 'split_rename.done'),
        'fill_empty': os.path.join(files, 'sub_0', 'filled_init_empty_pdbs.done')
    }
}


localrules: all
rule all:
    input: 
        done = outfile['done']['fill_empty']
    output:
        done = subworkflow_done
    run:
        with open(output.done, 'w') as f:
            f.write('Subworkflow 0 done at ')
            f.write(str(datetime.utcnow()))
            f.write('\n')



localrules: download_assemblies
rule download_assemblies:
    output:
        done = outfile['done']['rsync']
    params:
        rsync_target = outfile['assemblies_dir']
    shell:
        ''' 
        mkdir -p {params.rsync_target};
        rsync -rlpt -v -z -q --delete --update --port=33444 \
        rsync.rcsb.org::ftp_data/assemblies/mmCIF/divided/ {params.rsync_target};
        touch {output.done};
        '''
        
        # The above should be unaffected by the deprecation
        # of ftp protocol by rcsb; while it is "ftp_data",
        # the protocol itself is rsync which will remain.




rule ungzip_assemblies:
    input:
        done = outfile['done']['rsync']
    output:
        done = outfile['done']['ungzip']
    params:
        assemblies = outfile['assemblies_dir']
    threads: 32 
    resources:
        time = '12:00:00',
        mem_mb = '20000'
    shell:
        '''
        scripts/parallel_ungzip_all.sh {params.assemblies} {threads} 1000;
        touch {output.done};
        '''
        






rule split_rename_chains:
    input:
        done = outfile['done']['ungzip']
    output:
        done = outfile['done']['split']
    threads: 32 
    resources:
       mem_mb = '20000',
       time = '24:00:00'
    run:
        parallel_convert_split_rename_cifs(parent_dir=outfile['assemblies_dir'],
                                           cif2pdb_exe=cif2pdb_exe,
                                           threads=threads
                                           )
        shell(''' touch {output.done} ''')






rule index_chains_pdb:
    input:
        done = outfile['done']['split']
    output:
        indexfile = outfile['chains_index'],
        empty_files_txt = outfile['empty_pdbs']
    threads: 1
    resources:
        time = '4:00:00'
    run:
        empty_files = generate_rcsb_index(rcsb_path=outfile['assemblies_dir'], index_file=output.indexfile)
        with open(output.empty_files_txt, 'w') as f:
            for ef in empty_files:
                f.write(f'{ef}\n')




rule fill_empty_pdb_files:
    input:
        empty_files_txt = outfile['empty_pdbs']
    output:
        done = outfile['done']['fill_empty']
    threads: 1
    resources:
        mem_mb = '15000'
    run:
        with open(input.empty_files_txt, 'rb') as f:
            for line in f:
                empty_file = line.strip()
                fill_empty_pdb(filepath=empty_file, cif2pdb_exe=cif2pdb_exe) 
        shell(''' touch {output.done} ''')
            


















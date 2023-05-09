import os
import sys 
import pickle

configfile: 'config.yaml'

files = config['paths']['pipeline_files'] # TODO

targz_assemblies = f'{files}/no_sample/assemblies.tar.gz'
targz_cryo = f'{files}/no_sample/cryo.tar.gz'
assemblies_ungzip_done = f'{files}/no_sample/assemblies_ungzip.done'
cryo_ungzip_done = f'{files}/no_sample/cryo_ungzip.done'
assemblies_cryo_merged = f'{files}/no_sample/assemblies_cryo_merged.done'
split_renamed_all_chains_done = f'{files}/no_sample/split+renamed.done'




rule all:
    input: 
        split_renamed_all_chains_done
    run:
        raise NotImplementedError # TODO




rule download_assemblies:
    output:
        targz = targz_assemblies
    run:
        raise NotImplementedError # TODO



rule download_cryo:
    output:
        targz = targz_cryo
    run:
        raise NotImplementedError # TODO



rule ungzip_assemblies:
    input:
        targz = targz_assemblies
    output:
        assemblies_ungzip_done
    run:
        raise NotImplementedError # TODO



rule ungzip_cryo:
    input:
        targz = targz_cryo
    output:
        cryo_ungzip_done
    run:
        raise NotImplementedError # TODO



rule merge_assemblies_cryo:
    input:
        assemblies_ungzip_done,
        cryo_ungzip_done
    output:
        assemblies_cryo_merged
    run:
        raise NotImplementedError # TODO



rule split_rename_chains:
    output:
        split_renamed_all_chains_done
    run:
        raise NotImplementedError # TODO    




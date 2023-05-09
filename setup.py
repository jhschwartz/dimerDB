#!/usr/bin/python3
'''
This script writes a config.yaml based on your responses.
'''

raise NotImplementedError('this is a work in progress') # TODO

import os 

space = None
while space is None or not space in ['y','n']:
    space = input('Is there ???G available on the current disk? ')
if space == 'n':
    raise ValueError('You\'ll want to do this with more space available.')


needs_env = None
while needs_env is None or needs_env not in ['y', 'n']:
    needs_env = input('Install new conda env? (if first use of this pipeline, then yes) [y/n]: ')

if needs_env == 'y':
    with open('envneeded.maketmp','w') as _: pass
else:
    env_path = None
    while env_path is None or not os.path.exists(env_path): 
        env_path = input('Please enter the path to the conda environment: ')
    with open('envpath.maketmp', 'w') as f:
        f.write(env_path)

use_slurm = None
while use_slurm is None or use_slurm not in ['y','n']:
    use_slurm = input('Is this a system that uses Slurm? [y/n]: ')

if use_slurm:
    slurm_acct = None
    slurm_part = None
    node_cpu_count = None
    while slurm_acct is None:
        slurm_acct = input('What slurm account would you like to use? ')
    while slurm_part is None:
        slurm_part = input('What slurm queue/partition would you like to use? ')
    while node_cpu_count is None or not node_cpu_count.isdigit():
        node_cpu_count = input('How many cpus are on each slurm node? ')
with open('slurminfo.maketmp', 'w') as f:
    if not use_slurm:
        f.write('no_slurm')
    else:
        f.write(f'{slurm_acct}\t{slurm_part}\t{node_cpu_count}')









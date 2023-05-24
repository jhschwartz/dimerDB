from contextlib import contextmanager
import tempfile
import shutil 
import subprocess
import os
import glob

@contextmanager
def run_tmp_snakemake(snakefile, config, snakemake_exe, scripts, bin, rule=None, lib=None):
    init_dir = os.getcwd()
    try:
        td = tempfile.mkdtemp()
        tmp_snakefile = os.path.join(td, 'snakefile')

        # if no rule specified, then we are doing all, so let rule "all" run
        if not rule or rule == 'all':
            shutil.copy(snakefile, tmp_snakefile)
        else:
            with open(snakefile, 'r') as smk_in, \
                    open(tmp_snakefile, 'w') as smk_out:
                currently_reading_rule_all = False
                for line in smk_in:
                    
                    if 'BEGIN RULE ALL TARGETS' in line:
                        currently_reading_rule_all = True
                    elif 'END RULE ALL TARGETS' in line:
                        currently_reading_rule_all = False

                    if not currently_reading_rule_all:
                        smk_out.write(line)
                    else:
                        smk_out.write(f'# skipping rule all target...\n')

        if lib:
            shutil.copytree(lib, os.path.join(td, 'lib'))
        
        tmp_config = shutil.copy(config, os.path.join(td, 'config.yaml'))
        snakemake_exe = os.path.realpath(snakemake_exe)
        
        # copy scripts
        paths_scripts_depth_1 = glob.glob(f'{scripts}/*')
        paths_scripts_depth_1 = filter(lambda f: os.path.isfile(f), paths_scripts_depth_1)
        tmp_scripts = os.path.join(td, 'scripts')
        os.makedirs(tmp_scripts)
        for script in paths_scripts_depth_1:
            shutil.copy(script, tmp_scripts)

        # copy bin
        shutil.copytree(bin, os.path.join(td, 'bin'))
        
        os.chdir(td)
        rule_info = ''
        if rule:
            rule_info = f'-R {rule}'
        subprocess.run(f'{snakemake_exe} -s {tmp_snakefile} -c8 {rule_info} --printshellcmds', shell=True, check=True)
        yield tmp_snakefile
    finally:
        shutil.rmtree(td)
        os.chdir(init_dir)


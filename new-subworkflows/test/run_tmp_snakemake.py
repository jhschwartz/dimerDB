from contextlib import contextmanager
import tempfile
import shutil 
import subprocess
import os
import glob

@contextmanager
def run_tmp_snakemake(snakefile, config, snakemake_exe, scripts, rule=None):
    try:
        td = tempfile.mkdtemp()
        tmp_snakefile = shutil.copy(snakefile, td)
        tmp_config = shutil.copy(config, os.path.join(td, 'config.yaml'))
        snakemake_exe = os.path.realpath(snakemake_exe)
        
        # copy scripts
        paths_scripts_depth_1 = glob.glob(f'{scripts}/*')
        paths_scripts_depth_1 = filter(lambda f: os.path.isfile(f), paths_scripts_depth_1)
        tmp_scripts = os.path.join(td, 'scripts')
        os.makedirs(tmp_scripts)
        for script in paths_scripts_depth_1:
            shutil.copy(script, tmp_scripts)
        
        init_dir = os.getcwd()
        os.chdir(td)
        rule_info = ''
        if rule:
            rule_info = f'-R {rule}'
        subprocess.run(f'{snakemake_exe} -s {tmp_snakefile} -c1 {rule_info}', shell=True, check=True)
        yield tmp_snakefile
    finally:
        shutil.rmtree(td)
        os.chdir(init_dir)


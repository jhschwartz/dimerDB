import subprocess
import tempfile
import os


def run_mmseqs_cluster(infasta, outprefix, seq_id, mmseqs_exe, cores=1):
    '''
    This function wraps mmseqs' easy-cluster function,
    clustering the sequences of infasta to seq_id and
    writing the output files to outprefix.

    Output files are:
        {outprefix}_cluster.tsv
        {outprefix}_rep_seq.fasta,

    :param infasta: str, path to a fasta file of protein seqs
    :param outprefix: str, the output path prefix to write
                        the three output files
    :param seq_id: float, range [0, 100], the sequence identity
                    to which we are clustering the seqs

    :returns cluster_tsv: str, path to the resulting mmseqs file
                            which describes membership of input seqs
                            to clusters. Two columns, first is
                            cluster name and second is member seq name
    :returns reps_fasta: str, path to a fasta which includes only the
                            resulting seqs of clustering infasta to seq_id
    '''
    seq_id_scaled = float(seq_id) / 100

    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = f'{mmseqs_exe} easy-cluster {infasta} {outprefix} {tmpdir} --min-seq-id {seq_id_scaled} --cov-mode 0 --threads {cores}'
        with open(os.devnull, 'w') as NULL: # NULL to silence stdout
            subprocess.run(cmd, shell=True, check=True, stdout=NULL, stderr=subprocess.STDOUT)
        
        cluster_tsv = f'{outprefix}_cluster.tsv'
        reps_fasta = f'{outprefix}_rep_seq.fasta'
        
        all_fasta = f'{outprefix}_all_seqs.fasta'
        os.remove(all_fasta) # simply don't need this; it's identical to input

        return cluster_tsv, reps_fasta




def derive_mmseqs_reps_from_tsv(mmseqs_result_tsv):
    '''
    Simple function to read the resulting cluster_tsv from run_mmseqs_cluster
    and return the list of cluster names, which are representative sequence
    names of the clustered sequence set.

    :param mmseqs_result_tsv: str, path to an output file cluster_tsv from
                                run_mmseqs_cluster above
    
    :returns reps: list[str], list of sequence names which make up the
                    representatives of the resulting cluster sequences
    '''
    reps = set()
    with open(mmseqs_result_tsv, 'r') as f:
        for line in f:
            rep = line.split()[0]
            reps.add(rep)
    return sorted(list(reps))



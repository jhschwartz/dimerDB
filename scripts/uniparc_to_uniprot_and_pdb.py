import gzip
import yaml
from sortedcontainers import SortedSet, SortedDict



def _extract_chains(chain_str):
    #   'abcd:A; efghA:A; zxvd:B;' --> ['abcd_A', 'efgh_A', 'zxvd_B']
    chains_list = chain_str.split(';')
    chains_list_clean = [c.strip().replace(':','_') for c in chains_list]
    chains = set()
    for c in chains_list_clean:
        if c != '':
            chains.add(c)
    return chains


def _extract_ids(line):
    # read the ids of interest (uniprot, chains, and uniparc) from each line
    spl = line.split('\t')
    uniprot = spl[0] if spl[0] != '' else None
    uniparc = spl[10] if spl[10] != '' else None
    chains = _extract_chains(spl[5])
    return uniprot, chains, uniparc


def _sort_entries(entries):
    # convert each entry value to a list and sort it. We do this because yamls can't elegantly store sets.
    for uniparc in entries.keys():
        entries[uniparc]['uniprot'] = sorted(list(entries[uniparc]['uniprot']))
        entries[uniparc]['pdb'] = sorted(list(entries[uniparc]['pdb']))
    return entries


def _clear_entries_with_no_chains(entries):
    entries = {k:v for k,v in entries.items() if len(v['pdb']) > 0}
    return entries



def make_uniparc2others(infile, outfile):
    '''
    This function makes a dict, which is later saved in a yaml file, of single uniparc IDs matched to all 
    unicprot ACs and pdb chains that are that uniparc sequence. For example, for a protein with uniparc ID
    "UPI0000123567", there could be several matching uniprots like "XYZ123" and "ABC987", which each match
    pdb chains which could be, as a union, ['1abc_A', '1abc_B', '7cba_B']. The resulting entry would be:
        { 
            'UPI0000123567': {
                'uniprot': ['XYZ123', 'ABC987'],
                'pdb': ['1abc_A', '1abc_B', '7cba_B']
            }
        }
    The value of each entry, both of subkeys 'uniprot' and 'pdb', is a alphabetically sorted list without repeats.

    :param infile: str, the path of a idmapping_selected.tab.gz file from uniprot.org. This file is quite large (11G).
                    Each line is tab-delimited. We only care about the columns 0 (uniprotAC), 5 (PDB), and 10 (UniParcID).
                    The format of columns 0 and 10 (uniprot and uniparc) are simple.
                    The format of column 5 can contain zero, one, or many PDB chains in the format below:
                                "4H3B:B; 4H3B:D; 6DJL:B; 6DJL:C;" meaning "PDBCODE:CHAIN"
    
    :param outfile: str, the path of the output yaml file we are saving the dict in. Should end in .yml or .yaml.
    '''
    
    # init empty dict - this is ultimately what we will save
    entries = SortedDict() # SortedDict to make search must faster/finite

    # Yes, this is odd. But the id_mapping file is so big that we need to do memory management.
    # Remember that 1) we only care about uniparc IDs that have a matching PDB and that
    # 2) it is possible for a uniparc ID to be seen first without a PDB but later be seen
    # with a PDB but a different uniprot. To preserve uniprots, we will first create a list
    # of uniparc IDs that have a PDB, then reading from the beginning we will choose what to save.
    
    uniparcs_with_a_pdb = SortedSet() # SortedSet to make search must faster/finite
    with gzip.open(infile, 'rt') as f:
        while line := f.readline():
            uniprot, chains, uniparc = _extract_ids(line)
            if len(chains) > 0:
                uniparcs_with_a_pdb.add(uniparc)


    # open the g'zipped infile in read-text mode so that lines are read as strings 
    with gzip.open(infile, 'rt') as f:

        while line := f.readline():
           
            # pick out uniprot, pdbs, and uniparc of each line
            uniprot, chains, uniparc = _extract_ids(line)

            # if no structure for this uniparc, skip it
            if not uniparc in uniparcs_with_a_pdb:
                continue

            if uniparc in entries:
                entries[uniparc]['uniprot'].add(uniprot)
                entries[uniparc]['pdb'] = entries[uniparc]['pdb'].union(chains)     # does nothing if no chains
            else:
                entries[uniparc] = { 'uniprot': { uniprot }, 'pdb': set(chains) }

    
    entries = _clear_entries_with_no_chains(entries)
    entries = _sort_entries(entries)


    # write the dict to yaml
    with open(outfile, 'w') as f:
        yaml.dump(entries, f)


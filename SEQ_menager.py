from PDB_menager import AA_ATTRIBUTES

def read_fasta(filename):
    """Read fasta format file. Return string type seq."""
    fasta = ''
    with open(filename, 'r') as pdb:
        for line in pdb:
            if line[:1] != '>':
                fasta += line.rstrip('\n')
    return fasta


def read_conservation(clustal_aln, skylign):
    """Read matrix 21 cols x N rows, eg. skylign format file. Return conservation values for each position in query sequence."""
    conservation = []
    seq = ''
    try:
        with open(clustal_aln, 'r') as aln:
            for line in aln:
                if line.startswith('Query'):
                    seq += line[20:80]
    except it:
        seq = clustal_aln
    except:
        print("ERROR: Wrong format of clustal_aln! Provide clustal format file or string sequence matching to skylign file.")
    with open(skylign, 'r') as con:
        n = 0
        for line in con:
            if line[:1] != '#':
                tokens = line.split()[1:21]
                if seq[n] != '.':
                    conservation.append((seq[n], tokens[AA_ATTRIBUTES[seq[n]][0]]))
                n += 1
    return conservation

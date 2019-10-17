import re
from PDB_menager import AA_ATTRIBUTES

def read_fasta(filename):
    """Reads fasta format file. Return string type seq."""
    fasta = ''
    with open(filename, 'r') as pdb:
        for line in pdb:
            if line[:1] != '>':
                fasta += line.rstrip('\n')
    return fasta


def read_conservation(clustal_aln, skylign):
    """Reads matrix 21 cols x N rows, eg. skylign format file. Return conservation values for each position in query sequence."""
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


def read_LLPSDB(key):
    """Reads subset of LLPSDB database. @key is a subset name: DESIG, NATUR, C_ONE, C_TWO, C_MOR, P_DNA, P_RNA, P_PRT"""
    path="/Users/adawid/PROJECTS/databases/LLPSDB/"
    data=['/5_Uniprot_id', '/0_PID', '/14_DisProt_ID', '/7_Species', '/11_Sequence', '/12_IDR', '/13_LCR']
    vec=[]
    for n1, feature in enumerate(data):
        try:
            with open(path+key+feature, "r") as f:
                for n2, row in enumerate(f):
                    r=''
                    if n1 == 0:
                        vec.append([row.strip()])
                    elif n1 == 4:
                        for i in row.strip().split(';'):
                            if i.isalpha():
                                r += i
                        vec[n2].append(r)
                        vec[n2].append(len(r))
                    else:
                        vec[n2].append(row.strip())
        except:
            print("ERROR: File %s can't be open!"%str(path+key+feature))
    try:
        with open(path+key+"/1_P_name", "r") as names:
            for n, name in enumerate(names):
                LLPSDB[name.strip()] = vec[n]
    except:
        print("ERROR: File %s can't be open!"%str(path+key+'/1_P_name'))
    return LLPSDB


def aa_composition(sequence):######################
    """"""
    aa=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    counts={}
    length=len(sequence)
    for a in aa:
        counts[a] = sequence.count(a)/float(length)
        if a=='D' or a=='E' or a=='K' or a=='R' or a=='H':
            
        elif a=='N' or a=='Q' or a=='S' or a=='T' or a=='Y':
            
        if a=='F' or a=='W' or a=='H' or a=='Y':
            
    return counts


def search_motif(queries, sequence):
    """Searches of given as a list sequence motifs in given string sequence.
       * treats regular expressions in queries
       Dictionary output: key - motif, @1 - counts, @2+ - index in seq
    """
    counts={}
    for query in queries:
        s=re.compile(query)
        for q in s.findall(sequence):
            n=0
            if not q in counts.keys():
                counts[q] = [sequence.count(q)]
                for i in range(0, counts[q][0]):
                    position=sequence.find(q, n)
                    counts[q].append(position)
                    n=position+len(q)
    return counts


LLPSDB={}
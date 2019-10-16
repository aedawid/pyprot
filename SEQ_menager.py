import os
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
    r=''
    for n1, feature in enumerate(data):
        try:
            with open(path+key+feature, "r") as f:
                for n2, row in enumerate(f):
                    if n1 == 0:
                        vec.append([row.strip()])
                    elif n1 == 4:
                        for i in row.strip().split(';'):
                            if i.isalpha():
                                r += i+" "
                        vec[n2].append(r)
                        r=''
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

LLPSDB={}
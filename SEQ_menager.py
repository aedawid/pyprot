import re
import math
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
    path="/Users/adawid/Box/GROUP_data/LLPS-bioinformatics/databases/LLPSDB/"
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


def aa_composition(sequence):
    """Prepare statistics based on amino acid sequence."""
    counts={}
    length=len(sequence)
    for a in aa:
        ac = sequence.count(a)
        if ac != 0:
            counts[a] = ac/float(length)
        else:
            counts[a] = 0.0
    return counts


def aa_pattern(sequence):
    """Prepare statistics based on amino acid sequence."""
    hp=['V', 'L', 'I', 'M', 'F', 'C', 'W', 'A']
    pol=['S', 'D', 'N', 'T', 'E', 'Q', 'H', 'K', 'R', 'Y']
    pi=['F', 'Y', 'W', 'H', 'R', 'N', 'Q', 'D', 'E']
    ar=['F', 'Y', 'W', 'H']
    counts={}
    for key in ['hp', 'pol', 'pi', 'ch']:
        counts.setdefault(key, [])
        counts[key].append(sequence)
        counts[key].append(0.0)
    length=len(sequence)
    for a in aa:
        if a in hp:
            counts['hp'][0]=counts['hp'][0].replace(a, "1")
            counts['hp'][1]+=sequence.count(a)
        else:
            counts['hp'][0]=counts['hp'][0].replace(a, "0")
        if a in pol:
            counts['pol'][0]=counts['pol'][0].replace(a, "1")
            counts['pol'][1]+=sequence.count(a)
        else:
            counts['pol'][0]=counts['pol'][0].replace(a, "0")
        if a in pi:
            counts['pi'][1]+=sequence.count(a)
            if a in ar:
                counts['pi'][0]=counts['pi'][0].replace(a, "1")
            else:
                counts['pi'][0]=counts['pi'][0].replace(a, "2")
        else:
            counts['pi'][0]=counts['pi'][0].replace(a, "0")
        if a=='D' or a=='E':
            counts['ch'][0]=counts['ch'][0].replace(a, "1")
            counts['ch'][1]+=sequence.count(a)
        elif a=='K' or a=='R' or a=='H':
            counts['ch'][0]=counts['ch'][0].replace(a, "2")
            counts['ch'][1]+=sequence.count(a)
        else:
            counts['ch'][0]=counts['ch'][0].replace(a, "0")
    for key in ['hp', 'pol', 'pi', 'ch']:
        if counts[key][1] != 0.0:
            counts[key][1]/=float(length)
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


def sequence_charge_decoration(sequence):
    """Calculates SCD parameter of a given sequence.
       float output: SCD value
    """
    N=len(sequence)
    charge=charge_pattern(sequence)
    suma = 0.0
    for n,ch in enumerate(charge):
        for m in range(n+1,N):
            if ch != 0.0 and charge[m] != 0.0:
                suma += (charge[m]*ch*math.sqrt(m-n))
    return suma/N


def overall_charge_symmetry(sequence):
    """Calculates OCS parameter of a given sequence.
       float output: OCS value
    """
    N=len(sequence)
    blob=[5, 6]
    SUM=suma=0.0
    pos=(sequence.count('K') + sequence.count('R') + sequence.count('H'))/float(N)
    neg=(sequence.count('D') + sequence.count('E'))/float(N)
    FCR=pos+neg
    sigma=(pos-neg)*(pos-neg)/FCR
    for n in blob:
        m=n_seg=0
        suma=0.0
        length=n
        while m < N:
            if m+length >= N:
                length=N-m
            frag=sequence[m:m+length]
            pos=(frag.count('K') + frag.count('R') + frag.count('H'))/float(length)
            neg=(frag.count('D') + frag.count('E'))/float(length)
            if pos+neg == 0.0:
                sig_i=0.0
            else:
                sig_i=(pos-neg)*(pos-neg)/float(pos+neg)
            suma+=(sig_i-sigma)*(sig_i-sigma)
            n_seg+=1
            m+=length
        SUM+=suma/float(n_seg)
    return (FCR, SUM/float(len(blob)))

def charge_pattern(sequence):
    seq=list(sequence)
    charge=[0.0]*len(seq)
    for m,res in enumerate(seq):
        if res == 'K' or res == 'R' or res == 'H':
            charge[m] = 1.0
        elif res == 'D' or res == 'E':
            charge[m] = -1.0
    return charge


def write_fasta(seq_list):
    with open(seq_list, "r") as f:
        for row in f:
            tokens=row.split()
            out=open(tokens[0]+".fasta", 'a')
            out.write(">"+str(tokens[0])+"\n"+str(tokens[1])+"\n")
            out.close()


aa=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
LLPSDB={}
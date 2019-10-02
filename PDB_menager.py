    
class Pdb(object):

    def __init__(self, structure):

        self.structure = structure
        self.header = header
    
    
    def write_pdb(self, which = 1):

        n = which
        for model in self.structure:
            if n == which:
                print("MODEL%9s"%which)
                n += 1
            else:
                print("ENDMDL\nMODEL%9s"%n)
                n += 1
            for atom in model:
                print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s%6s%6s           %3s"%tuple(atom))
        print("ENDMDL")
    
    def write_model(self, which):

        print("MODEL%9s"%which)
        for atom in self.structure[which-1]:
            print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s%6s%6s           %3s"%tuple(atom))
        print("ENDMDL")
    
    def seq_from_struct(self):
        fasta = ''
        for atom in self.structure[0]:
            if atom[2] == ' CA ':
                fasta += AA_code(atom[3])
        return fasta
    
    def seq_from_header(self):
        fasta = ''
        for row in self.header['SEQRES']:
            tokens = row.split()
            for t in tokens:
                if len(t) == 3 and t.isalpha():
                    fasta += AA_code(t)
        return fasta
    
    def ss_for_struct(self):
        seq = self.seq_from_struct()
        ids = self.resi_ids()
        ch_r = self.chain_ranges()
        ss = 'C'*len(seq)
        s = list(ss)

        def parse_row(row, n, c):
            tokens = row.split()
            if tokens[2+n] == tokens[5+n]:
                for ch in ch_r:
                    if tokens[2+n] == ch[2]:
                        start = int(ch[0])
                        ss_first = start + int(tokens[3+n]) - int(ids[start])
                        ss_last = ss_first + int(tokens[6+n]) - int(tokens[3+n])
                        if seq[ss_first] == AA_code(tokens[1+n]):
                            for i in range(ss_first, ss_last + 1):
                                s[i] = c

        if len(seq) != len(ids):
            print("ERROR: len(seq) != len(ids)")
            print(len(seq), len(ids))
        else:
            for row in self.header['HELIX']:
                parse_row(row, 0, 'H')
            for row in self.header['SHEET']:
                parse_row(row, 1, 'E')

        ss = "".join(s)
        return ss
    
    def ss_for_seq(self, c = 'X'):
        ss = ''
        for i in self.resi_list():
            ss += i[3]
        return ss
    
    def compare_seqs(self):
        seq_full = list(self.seq_from_header())
        resi_list = self.resi_list()
        seq = ''
        if len(seq_full) != len(self.seq_from_struct()):
            for aa in range(0, len(seq_full)):
                if seq_full[aa] != resi_list[aa][2]:
                    print("ERROR: seq_from_struct + missing_resi differ from seq_from_header at position %s" %n)
                elif resi_list[aa][4] == 'm':
                    seq += '-'
                else:
                    seq += seq_full[aa]
        else:
            seq = "".join(seq_full)
        return seq
    
    def resi_list(self):
        seq_struct = list(self.seq_from_struct())
        ss_struct = list(self.ss_for_struct())
        r_ids = self.resi_ids()
        c_ids = self.chain_ids()
        missing = self.missing_resi()
        r_list = []
        for i in missing:
            r_list.append((int(i[2]), i[1], i[0], 'X', 'm'))
        for j in range(0, len(seq_struct)):
            r_list.append((int(r_ids[j]), c_ids[j], seq_struct[j], ss_struct[j], 's'))
        r_list = sorted(r_list, key = lambda tup: (tup[1], tup[0]))
        return r_list
    
    def missing_resi(self):
        missing = []
        for row in self.header['REMARK 465']:
            tokens = row.split()
            if len(tokens) == 3:
                missing.append((AA_code(tokens[0]), tokens[1], tokens[2]))
        return missing
    
    def __select(self, which):
        vec = []
        for atom in self.structure[0]:
            vec.append(atom[which])
        return vec
    
    def record_type(self):
        return self.__select(0)
    
    def atom_ids(self):
        return self.__select(1)
    
    def atom_names(self):
        return self.__select(2)
    
    def resi_names(self):
        return self.__select(3)
    
    def chain_ids(self):
        return self.__select(4)
    
    def chain_ranges(self):
        ch_ids = self.chain_ids()
        ranges = []
        if ch_ids[0] == ch_ids[len(ch_ids)-1]:
            ranges.append((0, len(ch_ids)-1))
        else:
            ch = ch_ids[0]
            first = 0
            n = 0
            for i in ch_ids:
                if ch != i:
                    ranges.append((first, n-1, ch))
                    first = n
                    ch = i
                elif n == len(ch_ids)-1:
                    ranges.append((first, n, ch))
                n += 1
        return ranges
    
    def resi_ids(self):
        return self.__select(5)
    
    def atom_x(self):
        return self.__select(6)
    
    def atom_y(self):
        return self.__select(7)
    
    def atom_z(self):
        return self.__select(8)
    
    def occupancy(self):
        return self.__select(9)
    
    def bfactor(self):
        return self.__select(10)
    
    def atom_type(self):
        return self.__select(11)
    
    def aa_feature(self, which):
        strctr = self.structure
        for atom in strctr[0]:
            atom[10] = aa_features[AA_code(atom[3])][which]
        return strctr


def read_pdb(filename, w_model = '0', w_chain = '0', w_atoms = [], alter = 'A'):
    
    def parse_line(line, model):

        atom = [line[:6], line[6:11], line[12:16], line[17:20],
        line[21], line[22:26], line[30:38], line[38:46],
        line[46:54], line[54:60], line[60:66], line[77:80]]
        if w_chain == '0':                         ##parse all chains
            if not len(w_atoms):                   ###parse all atoms
                model.append(atom)
            else:
                for at in w_atoms:                 ###parse atoms
                    if line[12:16] == at:
                        if line[16] == ' ' or line[16] == alter:
                            model.append(atom)
        elif line[21] == w_chain:                  ##parse single chain
            if not len(w_atoms):
                model.append(atom)
            else:
                for at in w_atoms:
                    if line[12:16] == at:
                        if line[16] == ' ' or line[16] == alter:
                            model.append(atom)
    
    def parse_header(line):
        for key in header:
            if line.startswith(key):
                header[key].append(line[11:80])
    
    model = []
    structure = []
    with open(filename, 'r') as pdb:
        if w_model == '0':                                 #parse all_models
            for line in pdb:
                if line[:4] == 'ATOM':# or line[:6] == "HETATM":
                    parse_line(line, model)
                elif line.startswith('ENDMDL'):
                    structure.append(model)
                    model = []
                else:
                    parse_header(line)
            if not len(structure):
                structure.append(model)
        else:                                              #parse single model
            is_ok = 'false'
            for line in pdb:
                if is_ok == 'true':
                    if line[:4] == 'ATOM':# or line[:6] == "HETATM":
                        parse_line(line, model)
                    elif line.startswith('ENDMDL'):
                        structure.append(model)
                        break
                elif line.startswith("MODEL%9s"%w_model):
                    is_ok = 'true'
                else:
                    parse_header(line)
    return structure


def pdb_to_fasta(filename):

    fasta = ''
    with open(filename, 'r') as pdb:
        for line in pdb:
            if line.startswith('ENDMDL'):
                break
            elif line[:4] == 'ATOM' or line[:6] == 'HETATM':
                if line[12:16] == ' CA ':
                    resid = AA_code(line[17:20])
                    fasta += resid
    return fasta


def read_fasta(filename):

    fasta = ''
    with open(filename, 'r') as pdb:
        for line in pdb:
            if line[:1] != '>':
                fasta += line.rstrip('\n')
    return fasta


def AA_code(resid):

    aa = {v: k for k, v in aa_names.items()}
    if len(resid) == 3:
        return aa[resid]
    else:
        return list(aa.keys())[list(aa.values()).index(resid)]


aa_names = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

##### 0 - id, 1 - weight, 2 - frequency, 3 - charge, 4 - polar, 5 - aromatic, 6 - hp_KD
aa_features = {
    'I': ( 0, 131.175, 5.49,  0.0, 0.0, 0.0,  4.5 ),
    'V': ( 1, 117.148, 6.73,  0.0, 0.0, 0.0,  4.2 ),
    'L': ( 2, 131.175, 9.68,  0.0, 0.0, 0.0,  3.8 ),
    'F': ( 3, 165.192, 3.87,  0.0, 0.0, 1.0,  2.8 ),
    'C': ( 4, 121.154, 1.38,  0.0, 0.0, 0.0,  2.5 ),
    'M': ( 5, 149.208, 2.32,  0.0, 0.0, 0.0,  1.9 ),
    'A': ( 6,  89.094, 8.76,  0.0, 0.0, 0.0,  1.8 ),
    'G': ( 7,  75.067, 7.03,  0.0, 0.0, 0.0, -0.4 ),
    'T': ( 8, 119.119, 5.53,  0.0, 1.0, 0.0, -0.7 ),
    'S': ( 9, 105.093, 7.14,  0.0, 1.0, 0.0, -0.8 ),
    'W': (10, 204.228, 1.25,  0.0, 0.0, 1.0, -0.9 ),
    'Y': (11, 181.191, 2.91,  0.0, 1.0, 1.0, -1.3 ),
    'P': (12, 115.132, 5.02,  0.0, 0.0, 0.0, -1.6 ),
    'H': (13, 155.156, 2.26,  0.5, 1.0, 1.0, -3.2 ),
    'E': (14, 147.131, 6.32, -1.0, 1.0, 0.0, -3.5 ),
    'Q': (15, 146.146, 3.90,  0.0, 1.0, 0.0, -3.5 ),
    'D': (16, 133.104, 5.49, -1.0, 1.0, 0.0, -3.5 ),
    'N': (17, 132.119, 3.93,  0.0, 1.0, 0.0, -3.5 ),
    'K': (18, 146.189, 5.19,  1.0, 1.0, 0.0, -3.9 ),
    'R': (19, 174.203, 5.78,  1.0, 1.0, 0.0, -4.5 )
}

header = {
    'HEADER'    : [],
    'TITLE'     : [],
    'COMPND'    : [],
    'SOURCE'    : [],
    'KEYWDS'    : [],
    'EXPDTA'    : [],
    'AUTHOR'    : [],
    'JRNL'      : [],
    'REMARK   2': [],
    'REMARK 800': [],
    'REMARK 465': [],
    'DBREF'     : [],
    'SEQRES'    : [],
    'HELIX'     : [],
    'SHEET'     : [],
    'CISPEP'    : []
}
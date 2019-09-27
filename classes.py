    
    
class Pdb(object):

    def __init__(self, structure):

        self.structure = structure
    
    
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
    
    def seq(self):
        fasta = ''
        for atom in self.structure[0]:
            if atom[2] == ' CA ':
                fasta += AA_code(atom[3])
        return fasta
    
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


class Chain(object):

    def __init__(self, seq, bond_length = 0.38):
        self.beads         = [AA(aa) for aa in seq]
        self.name          = [resid.name       for resid in self.beads]
        self.id            = [resid.id         for resid in self.beads]
        self.weight        = [resid.weight     for resid in self.beads]
        self.frequency     = [resid.frequency  for resid in self.beads]
        self.charge        = [resid.charge     for resid in self.beads]
        self.polarity      = [resid.polarity   for resid in self.beads]
        self.aromatic      = [resid.aromatic   for resid in self.beads]
        self.hp_KD         = [resid.hp_KD      for resid in self.beads]
        
        self.seq_length    = len(seq)
        self.nbonds        = len(seq)-1
        self.bond_length   = bond_length


class AA(object):

    def __init__(self, aa):
        if len(aa) == 3:
            aa = aa_names[resid]
        self.name = aa
        self.id = aa_features[aa][0]
        self.weight = aa_features[aa][1]
        self.frequency = aa_features[aa][2]
        self.charge = aa_features[aa][3]
        self.polarity = aa_features[aa][4]
        self.aromatic = aa_features[aa][5]
        self.hp_KD = aa_features[aa][6]

        self.x = 0.0
        self.y = 0.0
        self.z = 0.0


def read_pdb(filename, w_model = '0', w_chain = '0', w_atoms = []):
    
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
                        model.append(atom)
        elif line[21] == w_chain:                  ##parse single chain
            if not len(w_atoms):
                model.append(atom)
            else:
                for at in w_atoms:
                    if line[12:16] == at:
                        model.append(atom)
    
    model = []
    structure = []
    with open(filename, 'r') as pdb:
        if w_model == '0':                                 #parse all_models
            for line in pdb:
                if line[:4] == 'ATOM' or line[:6] == "HETATM":
                    parse_line(line, model)
                elif line.startswith('ENDMDL'):
                    structure.append(model)
                    model = []
        else:
            is_ok = 'false'
            for line in pdb:
                if is_ok == 'true':
                    if line[:4] == 'ATOM' or line[:6] == "HETATM":
                        parse_line(line, model)
                    elif line.startswith('ENDMDL'):
                        structure.append(model)
                        break
                elif line.startswith("MODEL%9s"%w_model):  #parse single model
                    is_ok = 'true'
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

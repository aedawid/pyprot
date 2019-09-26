class Pdb(object):

    def __init__(self, filename, w_model = '0', w_chain = '0', w_atoms = []):
        self.f  = filename
        self.m  = w_model
        self.ch = w_chain
        self.a  = w_atoms
#        read(filename, w_model, w_chain, w_atoms)
        
#        self.atom_ids    = 
#        self.atom_names  = 
#        self.resi_names  = 
#        self.chain_ids   = 
#        self.resi_ids    = 
#        self.atom_x      = 
#        self.atom_y      = 
#        self.atom_z      = 
#        self.bfactor     = 
#        self.atom_type   = 
        
        



    def read(self):
    
        def parse_line(line, model):

            atom = [line[:6], line[6:11], line[12:16], line[17:20],
            line[21], line[22:26], line[30:38], line[38:46],
            line[46:54], line[54:60], line[60:66], line[77:80]]
            if self.ch == '0':                  ##parse all chains
                if not len(self.a):             ###parse all atoms
                    model.append(atom)
                else:
                    for at in self.a:           ###parse atoms
                        if line[12:16] == at:
                            model.append(atom)
            elif line[21] == self.ch:           ##parse single chain
                if not len(self.a):
                    model.append(atom)
                else:
                    for at in self.a:
                        if line[12:16] == at:
                            model.append(atom)
    

        model = []
        structure = []
        with open(self.f, 'r') as pdb:
            if self.m == '0':                                      #parse all_models
                for line in pdb:
                    if line[:4] == 'ATOM' or line[:6] == "HETATM":
                        parse_line(line, model)
                        print(len(model))######################################
                    elif line.startswith('ENDMDL'):
                        structure.append(model)
                        model = []
            else:
                is_ok = 'false'
                for line in pdb:
                    if is_ok == 'true':
                        if line[:4] == 'ATOM' or line[:6] == "HETATM":
                            parse_line(line, model)
                            print(len(model))######################################
                        elif line.startswith('ENDMDL'):
                            structure.append(model)
                            break
                    elif line.startswith("MODEL%9s"%self.m):       #parse single model
                        is_ok = 'true'
        return structure
        
        

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

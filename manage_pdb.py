    


def AA_features(resid):

    if len(resid) == 3:
        resid = AA_code(resid)
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
    seq_keys   =   list(aa_features.keys())
    seq_values =   list(aa_features.values())

    aa_id      =   [value[0] for value in seq_values]
    weight     =   [value[1] for value in seq_values]
    frequency  =   [value[2] for value in seq_values]
    charge     =   [value[3] for value in seq_values]
    polarity   =   [value[4] for value in seq_values]
    aromatic   =   [value[5] for value in seq_values]
    hp_KD      =   [value[6] for value in seq_values]


def read_fasta(filename):

    fasta = ''
    with open(filename, 'r') as pdb:
        for line in pdb:
            if line[:1] != '>':
                fasta += line.rstrip('\n')
    return fasta


def AA_code(resid):

    aa_names = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    aa_names = {v: k for k, v in aa_names.items()}
    if len(resid) == 3:
        return aa_names[resid]
    else:
        return list(aa_names.keys())[list(aa_names.values()).index(resid)]


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

def read_pdb(filename, structure):

    model = []
    with open(filename, 'r') as pdb:
        for line in pdb:
            if line.startswith('ENDMDL'):
                structure.append(model)
                model = []
            elif line[:4] == 'ATOM' or line[:6] == "HETATM":
                atom = [line[:6], line[6:11], line[12:16], line[17:20],
                line[21], line[22:26], line[30:38], line[38:46],
                line[46:54], line[54:60], line[60:66], line[77:80]]
                model.append(atom)
    if not len(structure):
        structure.append(model)


def read_pdb_model(filename, w_model, w_chain, w_atoms, structure):

    model = []
    with open(filename, 'r') as pdb:
        if w_model == '0':					#parse all_models
            for line in pdb:
                if line[:4] == 'ATOM' or line[:6] == "HETATM":
                    __parse_line(line, w_chain, w_atoms, model)
                elif line.startswith('ENDMDL'):
                    structure.append(model)
                    model = []
        else:
            is_ok = 'false'
            for line in pdb:
                if is_ok == 'true':
                    if line[:4] == 'ATOM' or line[:6] == "HETATM":
                        __parse_line(line, w_chain, w_atoms, model)
                    elif line.startswith('ENDMDL'):
                        structure.append(model)
                        break
                elif line.startswith("MODEL%9s"%w_model):	#parse single model
                    is_ok = 'true'


def write_pdb(structure, which=1):

    n = which
    for model in structure:
        if n == which:
            print("MODEL%9s"%which)
            n += 1
        else:
            print("ENDMDL\nMODEL%9s"%n)
            n += 1
        for atom in model:
            print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s%6s%6s           %3s"%tuple(atom))
    print("ENDMDL")


def write_pdb_model(structure, which):

    print("MODEL%9s"%which)
    for atom in structure[which-1]:
        print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s%6s%6s           %3s"%tuple(atom))
    print("ENDMDL")


def __parse_line(line, w_chain, w_atoms, model):

    atom = [line[:6], line[6:11], line[12:16], line[17:20],
    line[21], line[22:26], line[30:38], line[38:46],
    line[46:54], line[54:60], line[60:66], line[77:80]]
    if w_chain == '0':			##parse all chains
        if not len(w_atoms):		###parse all atoms
            model.append(atom)
        else:
            for at in w_atoms:		###parse atoms
                if line[12:16] == at:
                    model.append(atom)
    elif line[21] == w_chain:		##parse single chain
        if not len(w_atoms):
            model.append(atom)
        else:
            for at in w_atoms:
                if line[12:16] == at:
                    model.append(atom)

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
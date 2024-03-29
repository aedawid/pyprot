class PDB(object):
    """Creates PDB object from loaded 'structure'.
       'structure' can be returned by function read_pdb(params).
       PDB object has a form of multilayer list, where:
       - the most external is 'structure', composed of:
       -- N 'models', composed of:
       --- M 'atoms', each atom has 12 string-type tokens with data
    """

    def __init__(self, structure):

        self.structure = structure
        """Returns structure object."""
        self.header = HEADER
        """Returns global variable that keeps selected info from PDB header."""
        ### ATOM SECTION
        self.record_type = self.__select(0)
        """Returns string list of record type ('ATOM' or 'HETATM') for all atoms. """
        self.atom_ids = self.__select(1)
        """Returns string list of atom ids for all atoms. """
        self.atom_names = self.__select(2)
        """Returns string list of atom names for all atoms."""
        self.resi_names = self.__select(3)
        """Returns string list of residue names for all atoms."""
        self.chain_ids = self.__select(4)
        """Returns string list of chain ids for all atoms."""
        self.resi_ids = self.__select(5)
        """Returns string list of residue ids for all atoms."""
        self.atom_x = self.__select(6)
        """Returns string list of x coordinate for all atoms."""
        self.atom_y = self.__select(7)
        """Returns string list of y coordinate for all atoms."""
        self.atom_z = self.__select(8)
        """Returns string list of z coordinate for all atoms."""
        self.occupancy = self.__select(9)
        """Returns string list of occupancy for all atoms."""
        self.bfactor = self.__select(10)
        """Returns string list of bfactor for all atoms."""
        self.atom_type = self.__select(11)
        """Returns string list of atom types for all atoms."""

        self.chains = self._chain_ranges()
        """Returns dictionary of chains in current structure: eg. ('A':(0, 50), 'B':(51:100))."""
        self.pdb_seq = self._seq_from_header()
        """Returns list of fasta strings for all chains in PDB header: eg. ('MSSGS', 'AGLSH')."""
        self.seq = self._seq_from_struct()
        """Returns list of fasta strings for chains in current structure: eg. ('MSSGS', 'AGLSH')."""
        self.pdb_chains = self._pdb_chains_ranges()
        """Returns dictionary of chains in PDB header: eg. ('A':(0, 50), 'B':(51:100))."""
        self.pdb_missing = self.missing_resi()
        """Returns list of missing residues in PDB header: eg. ()"""
        self.ss = self._ss_for_struct()
        """Returns list of strings of secondary structure in current structure: eg. ('CEEEC', 'CCHHH')."""
        self.residues = self.resi_list()
        """Returns list of consensus description of sequences in current structure."""
        self.missing = self.missing_resi(True)
        """Returns list of missing residues in current structure: eg. ()"""
        self.out_seq = self.outcome_seq()
        """Returns list of consensus fasta strings for chains in structure and header: eg. ('-SSGS', 'A--SH')."""
        self.pdb_ss = self._ss_for_seq()
        """Returns list of strings of secondary structure for seq from header: eg. ('CEEEC', 'CCHHH')."""
    
    ### LOADED STRUCTURE
    def write_pdb(self, which = 1):
        """Write 'structure' in pdb format."""
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
        """Write chosen 'model' in pdb format."""
        print("MODEL%9s"%which)
        for atom in self.structure[which-1]:
            print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s%6s%6s           %3s"%tuple(atom))
        print("ENDMDL")
    
    def _seq_from_struct(self):
        """Return AA sequences for chains from structure."""
        seq = []
        ch = self.structure[0][0][4]
        fasta = ''
        for atom in self.structure[0]:
            if atom[2] == ' CA ':
                if atom[4] == ch:
                    fasta += AA_code(atom[3])
                else:
                    seq.append(fasta)
                    ch = atom[4]
                    fasta = AA_code(atom[3])
        seq.append(fasta)
        return seq
    
    def _chain_ranges(self, per_resi=True):
        """Return list of ranges for chains; default per resi, else per atom."""
        if per_resi == True:
            ch_ids = self.__per_resi(self.chain_ids)
        else:
            ch_ids = self.chain_ids
        ranges = {}
        if ch_ids[0] == ch_ids[len(ch_ids)-1]:
            ranges[ch_ids[0]] = (0, len(ch_ids)-1, 0)
        else:
            n = 0
            ch = ch_ids[0]
            first = 0
            for num, i in enumerate(ch_ids):
                if ch != i:
                    ranges[ch] = (first, num-1, n)
                    first = num
                    ch = i
                    n += 1
                elif num == len(ch_ids)-1:
                    ranges[ch] = (first, num, n)
        return ranges
    
    def _ss_for_struct(self):
        """Return secondary structure assignment for residues from structure: eg. ('CHHHCC', 'CCEEEEECC')."""
        ss = []					#('cccccc', 'ccccccccc')
        def parse_row(row, n1, n2, n3, n4, c):
            if row[n1] == row[n2]:
                for num, ch in enumerate(ch_r):
                    if row[n1].upper() == ch:
                        for t, i in enumerate(range(ch_r[ch][0], ch_r[ch][1]+1)):
                            if int(ids[i]) == int(row[n3:n3+4]) and seq[num][t] == AA_code(row[n4:n4+3]):
                                for tt, k in enumerate(range(i, int(ch_r[ch][1])+1)):
                                    if int(ids[k]) == int(row[23:27]):
                                        s = list(ss[num])
                                        for ii in range(t, t+tt):
                                            s[ii] = c
                                        ss[num] = "".join(s)
                                        break
        seq = self.seq				#('ahagly', 'nanannana')
        ids = self.__per_resi(self.resi_ids)	#(33, 34, 35, 36, 37), N_resi, ids
        ch_r = self.chains			#('A':(0, 15), 'B':(16, 50)), idxs
        for chain in seq:
            ss.append('C'*len(chain))
        for row in self.header['HELIX']:
            parse_row(row, 9, 21, 11, 5, 'H')
        for row in self.header['SHEET']:
            parse_row(row, 11, 22, 12, 7, 'E')
        return ss
    
    ### PARSE HEADER
    def _seq_from_header(self):
        """Return full AA seq derived from PDB header (taking missing residues): eg. ('ahagly', 'nanannana')."""
        seq = []
        ch = self.header['SEQRES'][0].split()[0].upper()
        fasta = ''
        for row in self.header['SEQRES']:
            tokens = row.split()
            for t in tokens[2:]:
                if tokens[0].upper() != ch:
                    seq.append(fasta)
                    ch = tokens[0].upper()
                    if len(tokens) == 3:
                        if len(t.strip()) == 3:
                          fasta += AA_code(t)
                        else:
                          fasta += tokens[2].strip()
                    else:
                        fasta = AA_code(t)
                else:
                    if len(tokens) == 3:
                        if len(t.strip()) == 3:
                          fasta += AA_code(t)
                        else:
                          fasta += tokens[2].strip()
                    else:
                        fasta += AA_code(t)
        seq.append(fasta)
        return seq
    
    def _pdb_chains_ranges(self):
        """Returns ranges of residue_idxs for chains in full seq: eg. [(0, 100, 'A'), (101, 220, 'B')]."""
        seq = self.pdb_seq
        c = []
        s_range = {}
        first = 0
        for row in self.header['SEQRES']:
            if not row.split()[0].upper() in c:
                c.append(row.split()[0].upper())
        if len(c) == len(seq):
            for num, i in enumerate(c):
                s_range[i] = (first, first + len(seq[num])-1, num)
                first = first + len(seq[num])
        return s_range
    
    def _ss_for_seq(self):
        """Return secondary structure assignment for full seq, 'X' for unknown."""
        s = ''
        ss = []
        for chain in self.residues:
            for res in chain:
                s += res[3]
            ss.append(s)
        return ss
    
    def outcome_seq(self):
        """Return list of chain sequences, where missing residues are marked as '-'."""
        seq = []
        for num, ch in enumerate(self.chains.keys()):				#('A', 'B')
            resi_list = self.residues[num]
            s = list(self.pdb_seq[self.pdb_chains[ch][2]])
#            print("struct: ", self.seq[self.chains[ch][2]])##############################
#            print("seq   : ", "".join(s))###########################
#            print(len(self.seq[self.chains[ch][2]]), len(s), len(resi_list))#########################
            if len(s) != len(self.seq[self.chains[ch][2]]):
                for aa in range(0, len(s)):
                    if s[aa] != resi_list[aa][2]:
                        print("ERROR: seq_from_struct + missing_resi differ from seq_from_header at position %s" %aa)
                    if resi_list[aa][4] == 'm':
                        s[aa] = '-'
            seq.append("".join(s))
#        print("out_s : ", seq[0])#######################################
        return seq
    
    def resi_list(self):
        """Return list of tuples for full AA seq for structure composed of given chains:
           @1 - residue id
           @2 - chain code
           @3 - residue name (1 letter code)
           @4 - secondary structure assignment: H, E, C, X - unknown
           @5 - experimentally known coordinates: s - known, m - missing
        """
        r_l = []
        r_list = []
        r_ids = self.__per_resi(self.resi_ids)
        for ch in self.chains.keys():
            for i in self.pdb_missing:
                if i[1] == ch:
                    r_l.append((int(i[2]), i[1], i[0], 'X', 'm'))
            r_idx = int(r_ids[self.chains[ch][0]])
            for num, res in enumerate(range(self.chains[ch][0], self.chains[ch][1]+1)):
                if r_idx > int(r_ids[res]):
                    r_idx += 1
                else:
                    r_idx = int(r_ids[res])
                r_l.append((r_idx, ch, self.seq[self.chains[ch][2]][num], self.ss[self.chains[ch][2]][num], 's'))
            r_l = sorted(r_l, key = lambda tup: (tup[1], tup[0]))
            
            seq_h = list(self.pdb_seq[self.pdb_chains[ch][2]])
            if self.chains[ch][1]+1 - self.chains[ch][0] != len(seq_h):
                n = 0
                for num, aa in enumerate(seq_h):
                    if num >= len(r_l):
                        r_l.insert(num, (r_l[num-1][0]+1, ch, aa, 'X', 'm'))
                    elif aa != r_l[num][2]:
                        r_l.insert(num, ('nan', ch, aa, 'X', 'm'))
                        n += 1
                        if seq_h[num+1] == r_l[num+1][2]:###
                            for i in range(0, n):
                                tmp = list(r_l[num-i])
                                tmp[0] = r_l[num+1][0]-i-1
                                r_l[num-i] = tuple(tmp)
                            n = 0
                    elif num+1 < len(r_l)-1 and seq_h[num+1] != r_l[num+1][2]:###
                        if seq_h[num+1] == r_l[num][2] and seq_h[num+2] == r_l[num+1][2]:
                            r_l.insert(num, (r_l[num][0]-1, ch, aa, 'X', 'm'))
            r_list.append(r_l)
        return r_list
    
    def missing_resi(self, per_struct=False):
        """Return list of missing residues: resi name, chain, resi id."""
        missing = []
        if per_struct == False:
            for row in self.header['REMARK 465']:			#missing residues from header
                tokens = [row[5:8], row[9], row[11:16]]
                if tokens[2].strip().isdigit():
                    missing.append((AA_code(tokens[0]), tokens[1].upper(), tokens[2]))
            for row in self.header['REMARK 470']:			#missing CA atoms from header
                tokens = [row[5:8], row[9], row[10:15], row[17:21], row[22:26]]
                if tokens[2].strip().isdigit():
                    if tokens[3] == ' CA ' or tokens[4] == ' CA ':
                        missing.append((AA_code(tokens[0]), tokens[1].upper(), ''.join(c for c in tokens[2] if not c.isalpha())))
        else:
            for ch in self.residues:
                for res in ch:
                    if res[4] == 'm':
                        missing.append((res[2], res[1], res[0]))
        return missing
    
    def non_standard_resi(self):
        """Return list of non-standard AA in the structure."""
        r_ids = self.resi_ids
        r_ch = self.chain_ids
        resids=[]
        n=0
        for res in self.resi_names:
            if is_non_standard_AA(res):
                resids.append((res, AA_code(res), r_ids[n], r_ch[n]))
            n += 1
        return resids
    

    def aa_attribute(self, which):
        """Return a structure, where bfactor is replaced by chosen AA atribute:
        0 - id,
        1 - weight,
        2 - frequency,
        3 - charge,
        4 - polar,
        5 - aromatic,
        6 - hp_KD
        """
        strctr = self.structure
        for atom in strctr[0]:
            atom[10] = AA_ATTRIBUTES[AA_code(atom[3])][which]
        return strctr
    
    def binding_sites(self):
        """Return binding sites."""
        s = ''
        if len(self.header['REMARK 800']):
            for line in self.header['REMARK 800']:
                tokens = line.split(':')
                if tokens[0].startswith('SITE'):
                    if tokens[0] == 'SITE_IDENTIFIER':
                        s += tokens[1][1:4]+" -"
                    elif tokens[0] == 'SITE_DESCRIPTION':
                        s += tokens[1][:35]+","
        return s
    
    def binding_resi(self):
        """Return binding residues."""
        vec = []
        v = []
        if len(self.header['SITE']):
            site = self.header['SITE'][0].split()[0]
            for line in self.header['SITE']:
                if line[0:3] != site:
                    vec.append((site, v))
                    site = line[0:3]
                    v = []
                if line[7:10].isalpha() and line[7:10] != 'HOH':
                    v.append((line[7:10], line[12:16], line[11]))
                if line[18:21].isalpha() and line[18:21] != 'HOH':
                    v.append((line[18:21], line[23:27], line[22]))
                if line[29:32].isalpha() and line[29:32] != 'HOH':
                    v.append((line[29:32], line[34:38], line[33]))
                if line[40:43].isalpha() and line[40:43] != 'HOH':
                    v.append((line[40:43], line[45:49], line[44]))
            vec.append((site, v))
        return vec
    
    def ss_type(self):
        """Return type of protein due to its secondary structure: alpha, beta, albe."""
        ss = ''
        if self.is_alpha():
            ss += 'alpha,'
        if self.is_beta():
            ss += 'beta'
        return ss
    
    def is_alpha(self):
        """Return 'True' if structure contains helix."""
        return 'H' in list(self.pdb_ss)
    
    def is_beta(self):
        """Return 'True' if structure contains strand."""
        return 'E' in list(self.pdb_ss)
    
    def __select(self, which):
        """Return string list of chosen atom's attribute."""
        vec = []
        for atom in self.structure[0]:
            vec.append(atom[which])
        return vec
    
    def __per_resi(self, data):
        """Return data per resi."""
        a = self.atom_names
        vec = []
        for i in range(0, len(a)):
            if a[i] == ' CA ':
                vec.append(data[i])
        return vec


def read_pdb(filename, w_model = '0', w_chain = '0', w_atoms = [], alter = 'A'):
    """Read PDB file. Return sstructure as PDB object. Update HEADER variable.
    
       PARAMS:
       @filename - file in pdb format
       @w_model - choose single MODEL, eg. '1'; default '0' means 'all models'
       @w_chain - chose single CHAIN, eg. 'A'; default '0' means 'all chains'
       @w_atoms - chose multiple ATOMS, eg. [' CA ']; default [] means 'all atoms'
       @alter - choose alternative position, eg. 'A'; default 'A' means A variant
    """
    def parse_line(line, model):

        atom = [line[:6], line[6:11], line[12:16], line[17:20],
        line[21].upper(), line[22:26], line[30:38], line[38:46],
        line[46:54], line[54:60], line[60:66], line[77:80]]
        if w_chain == '0':                         ##parse all chains
            if not len(w_atoms):                   ###parse all atoms
                model.append(atom)
            else:
                for at in w_atoms:                 ###parse atoms
                    if line[12:16] == at:
                        model.append(atom)
        elif line[21].upper() == w_chain:          ##parse single chain
            if not len(w_atoms):
                model.append(atom)
            else:
                for at in w_atoms:
                    if line[12:16] == at:
                        model.append(atom)
    
    def parse_header(line):
        for key in HEADER:
            if line.startswith(key):
                if key == 'HET ':
                    HEADER[key].append(line[7:80])
                else:
                    HEADER[key].append(line[10:80])
    
    model = []
    structure = []
    with open(filename, 'r') as pdb:
        if w_model == '0':                                 #parse all_models
            for line in pdb:
                if line[:4] == 'ATOM' or line[:6] == "HETATM":
                    if line[16] == ' ' or line[16] == alter:
                        if line[17:20] in AA_MONOMERS.keys():
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
                    if line[:4] == 'ATOM' or line[:6] == "HETATM":
                        if line[16] == ' ' or line[16] == alter:
                            if line[17:20] in AA_MONOMERS.keys():
                                parse_line(line, model)
                    elif line.startswith('ENDMDL'):
                        structure.append(model)
                        break
                elif line.startswith("MODEL%9s"%w_model):
                    is_ok = 'true'
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    is_ok = 'true'
                    if line[17:20] in AA_MONOMERS.keys():
                        parse_line(line, model)
                else:
                    parse_header(line)
            if not len(structure):
                structure.append(model)
    return structure


def pdb_to_fasta(filename):
    """Convert pdb format file to string type seq. Equal to seq_from_struct()."""
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


def AA_code(resid):
    """Convert AA code: from 1-letter to 3-letter and vice versa.
       Return always one of 20 standard AA.
       Convert non-standard AA to standard one.
    """
    aa = {v: k for k, v in AA_CODES.items()}
    code = ''
    if len(resid) == 3:
        if resid in AA_CODES.values():
            code = aa[resid]
        elif resid in AA_MONOMERS.keys():
            code = aa[AA_MONOMERS[resid]]
        else:
            code = ''
#            print("The residue %s is unknown." %resid)
    elif len(resid) == 1:
        if resid in AA_CODES.keys():
            code = AA_CODES[resid]
        else:
            code = resid
    else:
        code = 'X'
    return code


def is_non_standard_AA(resid):
    """Return 'True' if resid is a non-standard known AA monomer."""
    if resid in AA_MONOMERS.keys():
        return not resid in AA_CODES.values()
    else:
        print("The residue %s is unknown." %resid)


HEADER = {
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
    'REMARK 470': [],
    'DBREF'     : [],
    'SEQRES'    : [],
    'HET '      : [],
    'HELIX'     : [],
    'SHEET'     : [],
    'CISPEP'    : [],
    'SITE'      : []
}

AA_CODES = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}


AA_ATTRIBUTES = {
    'A': ( 0,  89.094, 8.76,  0.0, 0.0, 0.0,  1.8 ),
    'C': ( 1, 121.154, 1.38,  0.0, 0.0, 0.0,  2.5 ),
    'D': ( 2, 133.104, 5.49, -1.0, 1.0, 0.0, -3.5 ),
    'E': ( 3, 147.131, 6.32, -1.0, 1.0, 0.0, -3.5 ),
    'F': ( 4, 165.192, 3.87,  0.0, 0.0, 1.0,  2.8 ),
    'G': ( 5,  75.067, 7.03,  0.0, 0.0, 0.0, -0.4 ),
    'H': ( 6, 155.156, 2.26,  0.5, 1.0, 1.0, -3.2 ),
    'I': ( 7, 131.175, 5.49,  0.0, 0.0, 0.0,  4.5 ),
    'K': ( 8, 146.189, 5.19,  1.0, 1.0, 0.0, -3.9 ),
    'L': ( 9, 131.175, 9.68,  0.0, 0.0, 0.0,  3.8 ),
    'M': (10, 149.208, 2.32,  0.0, 0.0, 0.0,  1.9 ),
    'N': (11, 132.119, 3.93,  0.0, 1.0, 0.0, -3.5 ),
    'P': (12, 115.132, 5.02,  0.0, 0.0, 0.0, -1.6 ),
    'Q': (13, 146.146, 3.90,  0.0, 1.0, 0.0, -3.5 ),
    'R': (14, 174.203, 5.78,  1.0, 1.0, 0.0, -4.5 ),
    'S': (15, 105.093, 7.14,  0.0, 1.0, 0.0, -0.8 ),
    'T': (16, 119.119, 5.53,  0.0, 1.0, 0.0, -0.7 ),
    'V': (17, 117.148, 6.73,  0.0, 0.0, 0.0,  4.2 ),
    'W': (18, 204.228, 1.25,  0.0, 0.0, 1.0, -0.9 ),
    'Y': (19, 181.191, 2.91,  0.0, 1.0, 1.0, -1.3 )
}

AA_MONOMERS = {
    '0CS': 'ALA',  # 0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
    '1AB': 'PRO',  # 1AB PRO  1,4-DIDEOXY-1,4-IMINO-D-ARABINITOL
    '1LU': 'LEU',  # 1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
    '1PA': 'PHE',  # 1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
    '1TQ': 'TRP',  # 1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
    '1TY': 'TYR',  # 1TY TYR
    '23F': 'PHE',  # 23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
    '23S': 'TRP',  # 23S TRP  MODIFIED TRYPTOPHAN
    '2BU': 'ALA',  # 2BU ADE
    '2ML': 'LEU',  # 2ML LEU  2-METHYLLEUCINE
    '2MR': 'ARG',  # 2MR ARG  N3, N4-DIMETHYLARGININE
    '2MT': 'PRO',  # 2MT PRO
    '2OP': 'ALA',  # 2OP ALA  (2S  2-HYDROXYPROPANAL
    '2TY': 'TYR',  # 2TY TYR
    '32S': 'TRP',  # 32S TRP  MODIFIED TRYPTOPHAN
    '32T': 'TRP',  # 32T TRP  MODIFIED TRYPTOPHAN
    '3AH': 'HIS',  # 3AH HIS
    '3MD': 'ASP',  # 3MD ASP  2S,3S-3-METHYLASPARTIC ACID
    '3TY': 'TYR',  # 3TY TYR  MODIFIED TYROSINE
    '4DP': 'TRP',  # 4DP TRP
    '4F3': 'ALA',  # 4F3 ALA  CYCLIZED
    '4FB': 'PRO',  # 4FB PRO  (4S)-4-FLUORO-L-PROLINE
    '4FW': 'TRP',  # 4FW TRP  4-FLUOROTRYPTOPHANE
    '4HT': 'TRP',  # 4HT TRP  4-HYDROXYTRYPTOPHAN
    '4IN': 'TRP',  # 4IN TRP  4-AMINO-L-TRYPTOPHAN
    '4PH': 'PHE',  # 4PH PHE  4-METHYL-L-PHENYLALANINE
    '5CS': 'CYS',  # 5CS CYS
    '6CL': 'LYS',  # 6CL LYS  6-CARBOXYLYSINE
    '6CW': 'TRP',  # 6CW TRP  6-CHLORO-L-TRYPTOPHAN
    'A0A': 'ASP',  # A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
    'AA4': 'ALA',  # AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
    'AAR': 'ARG',  # AAR ARG  ARGININEAMIDE
    'AB7': 'GLU',  # AB7 GLU  ALPHA-AMINOBUTYRIC ACID
    'ABA': 'ALA',  # ABA ALA  ALPHA-AMINOBUTYRIC ACID
    'ACB': 'ASP',  # ACB ASP  3-METHYL-ASPARTIC ACID
    'ACL': 'ARG',  # ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
    'ACY': 'GLY',  # ACY GLY  POST-TRANSLATIONAL MODIFICATION
    'AEI': 'THR',  # AEI THR  ACYLATED THR
    'AFA': 'ASN',  # AFA ASN  N-[7-METHYL-OCT-2,4-DIENOYL]ASPARAGINE
    'AGM': 'ARG',  # AGM ARG  4-METHYL-ARGININE
    'AGT': 'CYS',  # AGT CYS  AGMATINE-CYSTEINE ADDUCT
    'AHB': 'ASN',  # AHB ASN  BETA-HYDROXYASPARAGINE
    'AHO': 'ALA',  # AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
    'AHP': 'ALA',  # AHP ALA  2-AMINO-HEPTANOIC ACID
    'AIB': 'ALA',  # AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
    'AKL': 'ASP',  # AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
    'ALA': 'ALA',  # ALA ALA
    'ALC': 'ALA',  # ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
    'ALG': 'ARG',  # ALG ARG  GUANIDINOBUTYRYL GROUP
    'ALM': 'ALA',  # ALM ALA  1-METHYL-ALANINAL
    'ALN': 'ALA',  # ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
    'ALO': 'THR',  # ALO THR  ALLO-THREONINE
    'ALS': 'ALA',  # ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
    'ALT': 'ALA',  # ALT ALA  THIOALANINE
    'ALY': 'LYS',  # ALY LYS  N(6)-ACETYLLYSINE
    'AME': 'MET',  # AME MET  ACETYLATED METHIONINE
    'AP7': 'ALA',  # AP7 ADE
    'APH': 'ALA',  # APH ALA  P-AMIDINOPHENYL-3-ALANINE
    'API': 'LYS',  # API LYS  2,6-DIAMINOPIMELIC ACID
    'APK': 'LYS',  # APK LYS
    'AR2': 'ARG',  # AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
    'AR4': 'GLU',  # AR4 GLU
    'ARG': 'ARG',  # ARG ARG
    'ARM': 'ARG',  # ARM ARG  DEOXY-METHYL-ARGININE
    'ARO': 'ARG',  # ARO ARG  C-GAMMA-HYDROXY ARGININE
    'ASA': 'ASP',  # ASA ASP  ASPARTIC ALDEHYDE
    'ASB': 'ASP',  # ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
    'ASI': 'ASP',  # ASI ASP  L-ISO-ASPARTATE
    'ASK': 'ASP',  # ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
    'ASL': 'ASP',  # ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
    'ASN': 'ASN',  # ASN ASN
    'ASP': 'ASP',  # ASP ASP
    'AYA': 'ALA',  # AYA ALA  N-ACETYLALANINE
    'AYG': 'ALA',  # AYG ALA
    'AZK': 'LYS',  # AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
    'B2A': 'ALA',  # B2A ALA  ALANINE BORONIC ACID
    'B2F': 'PHE',  # B2F PHE  PHENYLALANINE BORONIC ACID
    'B2I': 'ILE',  # B2I ILE  ISOLEUCINE BORONIC ACID
    'B2V': 'VAL',  # B2V VAL  VALINE BORONIC ACID
    'B3A': 'ALA',  # B3A ALA  (3S)-3-AMINOBUTANOIC ACID
    'B3D': 'ASP',  # B3D ASP  3-AMINOPENTANEDIOIC ACID
    'B3E': 'GLU',  # B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
    'B3K': 'LYS',  # B3K LYS  (3S)-3,7-DIAMINOHEPTANOIC ACID
    'B3S': 'SER',  # B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
    'B3X': 'ASN',  # B3X ASN  (3S)-3,5-DIAMINO-5-OXOPENTANOIC ACID
    'B3Y': 'TYR',  # B3Y TYR
    'BAL': 'ALA',  # BAL ALA  BETA-ALANINE
    'BBC': 'CYS',  # BBC CYS
    'BCS': 'CYS',  # BCS CYS  BENZYLCYSTEINE
    'BCX': 'CYS',  # BCX CYS  BETA-3-CYSTEINE
    'BFD': 'ASP',  # BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
    'BG1': 'SER',  # BG1 SER
    'BHD': 'ASP',  # BHD ASP  BETA-HYDROXYASPARTIC ACID
    'BIF': 'PHE',  # BIF PHE
    'BLE': 'LEU',  # BLE LEU  LEUCINE BORONIC ACID
    'BLY': 'LYS',  # BLY LYS  LYSINE BORONIC ACID
    'BMT': 'THR',  # BMT THR
    'BNN': 'ALA',  # BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
    'BOR': 'ARG',  # BOR ARG
    'BPE': 'CYS',  # BPE CYS
    'BTR': 'TRP',  # BTR TRP  6-BROMO-TRYPTOPHAN
    'BUC': 'CYS',  # BUC CYS  S,S-BUTYLTHIOCYSTEINE
    'BUG': 'LEU',  # BUG LEU  TERT-LEUCYL AMINE
    'C12': 'ALA',  # C12 ALA
    'C1X': 'LYS',  # C1X LYS  MODIFIED LYSINE
    'C3Y': 'CYS',  # C3Y CYS  MODIFIED CYSTEINE
    'C5C': 'CYS',  # C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
    'C6C': 'CYS',  # C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
    'C99': 'ALA',  # C99 ALA
    'CAB': 'ALA',  # CAB ALA  4-CARBOXY-4-AMINOBUTANAL
    'CAF': 'CYS',  # CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
    'CAS': 'CYS',  # CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
    'CCS': 'CYS',  # CCS CYS  CARBOXYMETHYLATED CYSTEINE
    'CGU': 'GLU',  # CGU GLU  CARBOXYLATION OF THE CG ATOM
    'CH6': 'ALA',  # CH6 ALA
    'CH7': 'ALA',  # CH7 ALA
    'CHG': 'GLY',  # CHG GLY  CYCLOHEXYL GLYCINE
    'CHP': 'GLY',  # CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
    'CHS': 'PHE',  # CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
    'CIR': 'ARG',  # CIR ARG  CITRULLINE
    'CLB': 'ALA',  # CLB ALA
    'CLD': 'ALA',  # CLD ALA
    'CLE': 'LEU',  # CLE LEU  LEUCINE AMIDE
    'CLG': 'LYS',  # CLG LYS
    'CLH': 'LYS',  # CLH LYS
    'CLV': 'ALA',  # CLV ALA
    'CME': 'CYS',  # CME CYS  MODIFIED CYSTEINE
    'CML': 'CYS',  # CML CYS
    'CMT': 'CYS',  # CMT CYS  O-METHYLCYSTEINE
    'CQR': 'ALA',  # CQR ALA
    'CR2': 'ALA',  # CR2 ALA  POST-TRANSLATIONAL MODIFICATION
    'CR5': 'ALA',  # CR5 ALA
    'CR7': 'ALA',  # CR7 ALA
    'CR8': 'ALA',  # CR8 ALA
    'CRK': 'ALA',  # CRK ALA
    'CRO': 'THR',  # CRO THR  CYCLIZED
    'CRQ': 'TYR',  # CRQ TYR
    'CRW': 'ALA',  # CRW ALA
    'CRX': 'ALA',  # CRX ALA
    'CS1': 'CYS',  # CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
    'CS3': 'CYS',  # CS3 CYS
    'CS4': 'CYS',  # CS4 CYS
    'CSA': 'CYS',  # CSA CYS  S-ACETONYLCYSTEIN
    'CSB': 'CYS',  # CSB CYS  CYS BOUND TO LEAD ION
    'CSD': 'CYS',  # CSD CYS  3-SULFINOALANINE
    'CSE': 'CYS',  # CSE CYS  SELENOCYSTEINE
    'CSI': 'ALA',  # CSI ALA
    'CSO': 'CYS',  # CSO CYS  INE S-HYDROXYCYSTEINE
    'CSR': 'CYS',  # CSR CYS  S-ARSONOCYSTEINE
    'CSS': 'CYS',  # CSS CYS  1,3-THIAZOLE-4-CARBOXYLIC ACID
    'CSU': 'CYS',  # CSU CYS  CYSTEINE-S-SULFONIC ACID
    'CSW': 'CYS',  # CSW CYS  CYSTEINE-S-DIOXIDE
    'CSX': 'CYS',  # CSX CYS  OXOCYSTEINE
    'CSY': 'ALA',  # CSY ALA  MODIFIED TYROSINE COMPLEX
    'CSZ': 'CYS',  # CSZ CYS  S-SELANYL CYSTEINE
    'CTH': 'THR',  # CTH THR  4-CHLOROTHREONINE
    'CWR': 'ALA',  # CWR ALA
    'CXM': 'MET',  # CXM MET  N-CARBOXYMETHIONINE
    'CY0': 'CYS',  # CY0 CYS  MODIFIED CYSTEINE
    'CY1': 'CYS',  # CY1 CYS  ACETAMIDOMETHYLCYSTEINE
    'CY3': 'CYS',  # CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
    'CY4': 'CYS',  # CY4 CYS  S-BUTYRYL-CYSTEIN
    'CY7': 'CYS',  # CY7 CYS  MODIFIED CYSTEINE
    'CYD': 'CYS',  # CYD CYS
    'CYF': 'CYS',  # CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
    'CYG': 'CYS',  # CYG CYS
    'CYJ': 'LYS',  # CYJ LYS  MODIFIED LYSINE
    'CYQ': 'CYS',  # CYQ CYS
    'CYR': 'CYS',  # CYR CYS
    'CYS': 'CYS',  # CYS CYS
    'CZ2': 'CYS',  # CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
    'CZZ': 'CYS',  # CZZ CYS  THIARSAHYDROXY-CYSTEINE
    'DA2': 'ARG',  # DA2 ARG  MODIFIED ARGININE
    'DAB': 'ALA',  # DAB ALA  2,4-DIAMINOBUTYRIC ACID
    'DAH': 'PHE',  # DAH PHE  3,4-DIHYDROXYDAHNYLALANINE
    'DAL': 'ALA',  # DAL ALA  D-ALANINE
    'DAM': 'ALA',  # DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
    'DAR': 'ARG',  # DAR ARG  D-ARGININE
    'DAS': 'ASP',  # DAS ASP  D-ASPARTIC ACID
    'DBU': 'ALA',  # DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
    'DBY': 'TYR',  # DBY TYR  3,5 DIBROMOTYROSINE
    'DBZ': 'ALA',  # DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
    'DCL': 'LEU',  # DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
    'DCY': 'CYS',  # DCY CYS  D-CYSTEINE
    'DDE': 'HIS',  # DDE HIS
    'DGL': 'GLU',  # DGL GLU  D-GLU
    'DGN': 'GLN',  # DGN GLN  D-GLUTAMINE
    'DHA': 'ALA',  # DHA ALA  2-AMINO-ACRYLIC ACID
    'DHI': 'HIS',  # DHI HIS  D-HISTIDINE
    'DHL': 'SER',  # DHL SER  POST-TRANSLATIONAL MODIFICATION
    'DIL': 'ILE',  # DIL ILE  D-ISOLEUCINE
    'DIV': 'VAL',  # DIV VAL  D-ISOVALINE
    'DLE': 'LEU',  # DLE LEU  D-LEUCINE
    'DLS': 'LYS',  # DLS LYS  DI-ACETYL-LYSINE
    'DLY': 'LYS',  # DLY LYS  D-LYSINE
    'DMH': 'ASN',  # DMH ASN  N4,N4-DIMETHYL-ASPARAGINE
    'DMK': 'ASP',  # DMK ASP  DIMETHYL ASPARTIC ACID
    'DNE': 'LEU',  # DNE LEU  D-NORLEUCINE
    'DNG': 'LEU',  # DNG LEU  N-FORMYL-D-NORLEUCINE
    'DNL': 'LYS',  # DNL LYS  6-AMINO-HEXANAL
    'DNM': 'LEU',  # DNM LEU  D-N-METHYL NORLEUCINE
    'DPH': 'PHE',  # DPH PHE  DEAMINO-METHYL-PHENYLALANINE
    'DPL': 'PRO',  # DPL PRO  4-OXOPROLINE
    'DPN': 'PHE',  # DPN PHE  D-CONFIGURATION
    'DPP': 'ALA',  # DPP ALA  DIAMMINOPROPANOIC ACID
    'DPQ': 'TYR',  # DPQ TYR  TYROSINE DERIVATIVE
    'DPR': 'PRO',  # DPR PRO  D-PROLINE
    'DSE': 'SER',  # DSE SER  D-SERINE N-METHYLATED
    'DSG': 'ASN',  # DSG ASN  D-ASPARAGINE
    'DSN': 'SER',  # DSN SER  D-SERINE
    'DTH': 'THR',  # DTH THR  D-THREONINE
    'DTR': 'TRP',  # DTR TRP  D-TRYPTOPHAN
    'DTY': 'TYR',  # DTY TYR  D-TYROSINE
    'DVA': 'VAL',  # DVA VAL  D-VALINE
    'DYG': 'ALA',  # DYG ALA
    'DYS': 'CYS',  # DYS CYS
    'EFC': 'CYS',  # EFC CYS  S,S-(2-FLUOROETHYL)THIOCYSTEINE
    'ESB': 'TYR',  # ESB TYR
    'ESC': 'MET',  # ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
    'FCL': 'PHE',  # FCL PHE  3-CHLORO-L-PHENYLALANINE
    'FGL': 'ALA',  # FGL ALA  2-AMINOPROPANEDIOIC ACID
    'FGP': 'SER',  # FGP SER
    'FHL': 'LYS',  # FHL LYS  MODIFIED LYSINE
    'FLE': 'LEU',  # FLE LEU  FUROYL-LEUCINE
    'FLT': 'TYR',  # FLT TYR  FLUOROMALONYL TYROSINE
    'FME': 'MET',  # FME MET  FORMYL-METHIONINE
    'FOE': 'CYS',  # FOE CYS
    'FOG': 'PHE',  # FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
    'FOR': 'MET',  # FOR MET
    'FRF': 'PHE',  # FRF PHE  PHE FOLLOWED BY REDUCED PHE
    'FTR': 'TRP',  # FTR TRP  FLUOROTRYPTOPHANE
    'FTY': 'TYR',  # FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
    'GHG': 'GLN',  # GHG GLN  GAMMA-HYDROXY-GLUTAMINE
    'GHP': 'GLY',  # GHP GLY  4-HYDROXYPHENYLGLYCINE
    'GL3': 'GLY',  # GL3 GLY  POST-TRANSLATIONAL MODIFICATION
    'GLH': 'GLN',  # GLH GLN
    'GLN': 'GLN',  # GLN GLN
    'GLU': 'GLU',  # GLU GLU
    'GLY': 'GLY',  # GLY GLY
    'GLZ': 'GLY',  # GLZ GLY  AMINO-ACETALDEHYDE
    'GMA': 'GLU',  # GMA GLU  1-AMIDO-GLUTAMIC ACID
    'GMU': 'ALA',  # GMU 5MU
    'GPL': 'LYS',  # GPL LYS  LYSINE GUANOSINE-5'-MONOPHOSPHATE
    'GT9': 'CYS',  # GT9 CYS  SG ALKYLATED
    'GVL': 'SER',  # GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
    'GYC': 'CYS',  # GYC CYS
    'GYS': 'GLY',  # GYS GLY
    'H5M': 'PRO',  # H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
    'HHK': 'ALA',  # HHK ALA  (2S)-2,8-DIAMINOOCTANOIC ACID
    'HIA': 'HIS',  # HIA HIS  L-HISTIDINE AMIDE
    'HIC': 'HIS',  # HIC HIS  4-METHYL-HISTIDINE
    'HIP': 'HIS',  # HIP HIS  ND1-PHOSPHONOHISTIDINE
    'HIQ': 'HIS',  # HIQ HIS  MODIFIED HISTIDINE
    'HIS': 'HIS',  # HIS HIS
    'HLU': 'LEU',  # HLU LEU  BETA-HYDROXYLEUCINE
    'HMF': 'ALA',  # HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
    'HMR': 'ARG',  # HMR ARG  BETA-HOMOARGININE
    'HPE': 'PHE',  # HPE PHE  HOMOPHENYLALANINE
    'HPH': 'PHE',  # HPH PHE  PHENYLALANINOL GROUP
    'HPQ': 'PHE',  # HPQ PHE  HOMOPHENYLALANINYLMETHANE
    'HRG': 'ARG',  # HRG ARG  L-HOMOARGININE
    'HSE': 'SER',  # HSE SER  L-HOMOSERINE
    'HSL': 'SER',  # HSL SER  HOMOSERINE LACTONE
    'HSO': 'HIS',  # HSO HIS  HISTIDINOL
    'HTI': 'CYS',  # HTI CYS
    'HTR': 'TRP',  # HTR TRP  BETA-HYDROXYTRYPTOPHANE
    'HY3': 'PRO',  # HY3 PRO  3-HYDROXYPROLINE
    'HYP': 'PRO',  # HYP PRO  4-HYDROXYPROLINE
    'IAM': 'ALA',  # IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
    'IAS': 'ASP',  # IAS ASP  ASPARTYL GROUP
    'IGL': 'ALA',  # IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
    'IIL': 'ILE',  # IIL ILE  ISO-ISOLEUCINE
    'ILE': 'ILE',  # ILE ILE
    'ILG': 'GLU',  # ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
    'ILX': 'ILE',  # ILX ILE  4,5-DIHYDROXYISOLEUCINE
    'IML': 'ILE',  # IML ILE  N-METHYLATED
    'IPG': 'GLY',  # IPG GLY  N-ISOPROPYL GLYCINE
    'IT1': 'LYS',  # IT1 LYS
    'IYR': 'TYR',  # IYR TYR  3-IODO-TYROSINE
    'KCX': 'LYS',  # KCX LYS  CARBAMOYLATED LYSINE
    'KGC': 'LYS',  # KGC LYS
    'KOR': 'CYS',  # KOR CYS  MODIFIED CYSTEINE
    'KST': 'LYS',  # KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
    'KYN': 'ALA',  # KYN ALA  KYNURENINE
    'LA2': 'LYS',  # LA2 LYS
    'LAL': 'ALA',  # LAL ALA  N,N-DIMETHYL-L-ALANINE
    'LCK': 'LYS',  # LCK LYS
    'LCX': 'LYS',  # LCX LYS  CARBAMYLATED LYSINE
    'LDH': 'LYS',  # LDH LYS  N~6~-ETHYL-L-LYSINE
    'LED': 'LEU',  # LED LEU  POST-TRANSLATIONAL MODIFICATION
    'LEF': 'LEU',  # LEF LEU  2-5-FLUOROLEUCINE
    'LET': 'LYS',  # LET LYS  ODIFIED LYSINE
    'LEU': 'LEU',  # LEU LEU
    'LLP': 'LYS',  # LLP LYS
    'LLY': 'LYS',  # LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
    'LME': 'GLU',  # LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
    'LNT': 'LEU',  # LNT LEU
    'LPD': 'PRO',  # LPD PRO  L-PROLINAMIDE
    'LSO': 'LYS',  # LSO LYS  MODIFIED LYSINE
    'LYM': 'LYS',  # LYM LYS  DEOXY-METHYL-LYSINE
    'LYN': 'LYS',  # LYN LYS  2,6-DIAMINO-HEXANOIC ACID AMIDE
    'LYP': 'LYS',  # LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
    'LYR': 'LYS',  # LYR LYS  MODIFIED LYSINE
    'LYS': 'LYS',  # LYS LYS
    'LYX': 'LYS',  # LYX LYS  N''-(2-COENZYME A)-PROPANOYL-LYSINE
    'LYZ': 'LYS',  # LYZ LYS  5-HYDROXYLYSINE
    'M0H': 'CYS',  # M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
    'M2L': 'LYS',  # M2L LYS
    'M3L': 'LYS',  # M3L LYS  N-TRIMETHYLLYSINE
    'MAA': 'ALA',  # MAA ALA  N-METHYLALANINE
    'MAI': 'ARG',  # MAI ARG  DEOXO-METHYLARGININE
    'MBQ': 'TYR',  # MBQ TYR
    'MC1': 'SER',  # MC1 SER  METHICILLIN ACYL-SERINE
    'MCL': 'LYS',  # MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
    'MCS': 'CYS',  # MCS CYS  MALONYLCYSTEINE
    'MDO': 'ALA',  # MDO ALA
    'MEA': 'PHE',  # MEA PHE  N-METHYLPHENYLALANINE
    'MEG': 'GLU',  # MEG GLU  (2S,3R)-3-METHYL-GLUTAMIC ACID
    'MEN': 'ASN',  # MEN ASN  GAMMA METHYL ASPARAGINE
    'MET': 'MET',  # MET MET
    'MEU': 'GLY',  # MEU GLY  O-METHYL-GLYCINE
    'MFC': 'ALA',  # MFC ALA  CYCLIZED
    'MGG': 'ARG',  # MGG ARG  MODIFIED D-ARGININE
    'MGN': 'GLN',  # MGN GLN  2-METHYL-GLUTAMINE
    'MHL': 'LEU',  # MHL LEU  N-METHYLATED, HYDROXY
    'MHO': 'MET',  # MHO MET  POST-TRANSLATIONAL MODIFICATION
    'MHS': 'HIS',  # MHS HIS  1-N-METHYLHISTIDINE
    'MIS': 'SER',  # MIS SER  MODIFIED SERINE
    'MLE': 'LEU',  # MLE LEU  N-METHYLATED
    'MLL': 'LEU',  # MLL LEU  METHYL L-LEUCINATE
    'MLY': 'LYS',  # MLY LYS  METHYLATED LYSINE
    'MLZ': 'LYS',  # MLZ LYS  N-METHYL-LYSINE
    'MME': 'MET',  # MME MET  N-METHYL METHIONINE
    'MNL': 'LEU',  # MNL LEU  4,N-DIMETHYLNORLEUCINE
    'MNV': 'VAL',  # MNV VAL  N-METHYL-C-AMINO VALINE
    'MPQ': 'GLY',  # MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
    'MSA': 'GLY',  # MSA GLY  (2-S-METHYL) SARCOSINE
    'MSE': 'MET',  # MSE MET  ELENOMETHIONINE
    'MSO': 'MET',  # MSO MET  METHIONINE SULFOXIDE
    'MTY': 'PHE',  # MTY PHE  3-HYDROXYPHENYLALANINE
    'MVA': 'VAL',  # MVA VAL  N-METHYLATED
    'N10': 'SER',  # N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
    'NAL': 'ALA',  # NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
    'NAM': 'ALA',  # NAM ALA  NAM NAPTHYLAMINOALANINE
    'NBQ': 'TYR',  # NBQ TYR
    'NC1': 'SER',  # NC1 SER  NITROCEFIN ACYL-SERINE
    'NCB': 'ALA',  # NCB ALA  CHEMICAL MODIFICATION
    'NEP': 'HIS',  # NEP HIS  N1-PHOSPHONOHISTIDINE
    'NFA': 'PHE',  # NFA PHE  MODIFIED PHENYLALANINE
    'NIY': 'TYR',  # NIY TYR  META-NITRO-TYROSINE
    'NLE': 'LEU',  # NLE LEU  NORLEUCINE
    'NLN': 'LEU',  # NLN LEU  NORLEUCINE AMIDE
    'NLO': 'LEU',  # NLO LEU  O-METHYL-L-NORLEUCINE
    'NMC': 'GLY',  # NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
    'NMM': 'ARG',  # NMM ARG  MODIFIED ARGININE
    'NPH': 'CYS',  # NPH CYS
    'NRQ': 'ALA',  # NRQ ALA
    'NVA': 'VAL',  # NVA VAL  NORVALINE
    'NYC': 'ALA',  # NYC ALA
    'NYS': 'CYS',  # NYS CYS
    'NZH': 'HIS',  # NZH HIS
    'OAS': 'SER',  # OAS SER  O-ACETYLSERINE
    'OBS': 'LYS',  # OBS LYS  MODIFIED LYSINE
    'OCS': 'CYS',  # OCS CYS  CYSTEINE SULFONIC ACID
    'OCY': 'CYS',  # OCY CYS  HYDROXYETHYLCYSTEINE
    'OHI': 'HIS',  # OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
    'OHS': 'ASP',  # OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
    'OLT': 'THR',  # OLT THR  O-METHYL-L-THREONINE
    'OMT': 'MET',  # OMT MET  METHIONINE SULFONE
    'OPR': 'ARG',  # OPR ARG  C-(3-OXOPROPYL)ARGININE
    'ORN': 'ALA',  # ORN ALA  ORNITHINE
    'ORQ': 'ARG',  # ORQ ARG  N~5~-ACETYL-L-ORNITHINE
    'OSE': 'SER',  # OSE SER  O-SULFO-L-SERINE
    'OTY': 'TYR',  # OTY TYR
    'OXX': 'ASP',  # OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
    'P1L': 'CYS',  # P1L CYS  S-PALMITOYL CYSTEINE
    'P2Y': 'PRO',  # P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
    'PAQ': 'TYR',  # PAQ TYR  SEE REMARK 999
    'PAT': 'TRP',  # PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
    'PBB': 'CYS',  # PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
    'PBF': 'PHE',  # PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
    'PCA': 'PRO',  # PCA PRO  5-OXOPROLINE
    'PCS': 'PHE',  # PCS PHE  PHENYLALANYLMETHYLCHLORIDE
    'PEC': 'CYS',  # PEC CYS  S,S-PENTYLTHIOCYSTEINE
    'PF5': 'PHE',  # PF5 PHE  2,3,4,5,6-PENTAFLUORO-L-PHENYLALANINE
    'PFF': 'PHE',  # PFF PHE  4-FLUORO-L-PHENYLALANINE
    'PG1': 'SER',  # PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
    'PG9': 'GLY',  # PG9 GLY  D-PHENYLGLYCINE
    'PHA': 'PHE',  # PHA PHE  PHENYLALANINAL
    'PHD': 'ASP',  # PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
    'PHE': 'PHE',  # PHE PHE
    'PHI': 'PHE',  # PHI PHE  IODO-PHENYLALANINE
    'PHL': 'PHE',  # PHL PHE  L-PHENYLALANINOL
    'PHM': 'PHE',  # PHM PHE  PHENYLALANYLMETHANE
    'PIA': 'ALA',  # PIA ALA  FUSION OF ALA 65, TYR 66, GLY 67
    'PLE': 'LEU',  # PLE LEU  LEUCINE PHOSPHINIC ACID
    'PM3': 'PHE',  # PM3 PHE
    'POM': 'PRO',  # POM PRO  CIS-5-METHYL-4-OXOPROLINE
    'PPH': 'LEU',  # PPH LEU  PHENYLALANINE PHOSPHINIC ACID
    'PPN': 'PHE',  # PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
    'PR3': 'CYS',  # PR3 CYS  INE DTT-CYSTEINE
    'PRO': 'PRO',  # PRO PRO
    'PRQ': 'PHE',  # PRQ PHE  PHENYLALANINE
    'PRR': 'ALA',  # PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
    'PRS': 'PRO',  # PRS PRO  THIOPROLINE
    'PSA': 'PHE',  # PSA PHE
    'PSH': 'HIS',  # PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
    'PTH': 'TYR',  # PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
    'PTM': 'TYR',  # PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
    'PTR': 'TYR',  # PTR TYR  O-PHOSPHOTYROSINE
    'PYA': 'ALA',  # PYA ALA  3-(1,10-PHENANTHROL-2-YL)-L-ALANINE
    'PYC': 'ALA',  # PYC ALA  PYRROLE-2-CARBOXYLATE
    'PYR': 'SER',  # PYR SER  CHEMICALLY MODIFIED
    'PYT': 'ALA',  # PYT ALA  MODIFIED ALANINE
    'PYX': 'CYS',  # PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
    'R1A': 'CYS',  # R1A CYS
    'R1B': 'CYS',  # R1B CYS
    'R1F': 'CYS',  # R1F CYS
    'R7A': 'CYS',  # R7A CYS
    'RC7': 'ALA',  # RC7 ALA
    'RCY': 'CYS',  # RCY CYS
    'S1H': 'SER',  # S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
    'SAC': 'SER',  # SAC SER  N-ACETYL-SERINE
    'SAH': 'CYS',  # SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
    'SAR': 'GLY',  # SAR GLY  SARCOSINE
    'SBD': 'SER',  # SBD SER
    'SBG': 'SER',  # SBG SER  MODIFIED SERINE
    'SBL': 'SER',  # SBL SER
    'SC2': 'CYS',  # SC2 CYS  N-ACETYL-L-CYSTEINE
    'SCH': 'CYS',  # SCH CYS  S-METHYL THIOCYSTEINE GROUP
    'SCS': 'CYS',  # SCS CYS  MODIFIED CYSTEINE
    'SCY': 'CYS',  # SCY CYS  CETYLATED CYSTEINE
    'SDP': 'SER',  # SDP SER
    'SEB': 'SER',  # SEB SER  O-BENZYLSULFONYL-SERINE
    'SEC': 'ALA',  # SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
    'SEL': 'SER',  # SEL SER  2-AMINO-1,3-PROPANEDIOL
    'SEP': 'SER',  # SEP SER  E PHOSPHOSERINE
    'SER': 'SER',  # SER SER
    'SET': 'SER',  # SET SER  AMINOSERINE
    'SGB': 'SER',  # SGB SER  MODIFIED SERINE
    'SGR': 'SER',  # SGR SER  MODIFIED SERINE
    'SHC': 'CYS',  # SHC CYS  S-HEXYLCYSTEINE
    'SHP': 'GLY',  # SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
    'SIC': 'ALA',  # SIC ALA
    'SLZ': 'LYS',  # SLZ LYS  L-THIALYSINE
    'SMC': 'CYS',  # SMC CYS  POST-TRANSLATIONAL MODIFICATION
    'SME': 'MET',  # SME MET  METHIONINE SULFOXIDE
    'SMF': 'PHE',  # SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
    'SNC': 'CYS',  # SNC CYS  S-NITROSO CYSTEINE
    'SNN': 'ASP',  # SNN ASP  POST-TRANSLATIONAL MODIFICATION
    'SOC': 'CYS',  # SOC CYS  DIOXYSELENOCYSTEINE
    'SOY': 'SER',  # SOY SER  OXACILLOYL-ACYLATED SERINE
    'SUI': 'ALA',  # SUI ALA
    'SUN': 'SER',  # SUN SER  TABUN CONJUGATED SERINE
    'SVA': 'SER',  # SVA SER  SERINE VANADATE
    'SVV': 'SER',  # SVV SER  MODIFIED SERINE
    'SVX': 'SER',  # SVX SER  MODIFIED SERINE
    'SVY': 'SER',  # SVY SER  MODIFIED SERINE
    'SVZ': 'SER',  # SVZ SER  MODIFIED SERINE
    'SXE': 'SER',  # SXE SER  MODIFIED SERINE
    'TBG': 'GLY',  # TBG GLY  T-BUTYL GLYCINE
    'TBM': 'THR',  # TBM THR
    'TCQ': 'TYR',  # TCQ TYR  MODIFIED TYROSINE
    'TEE': 'CYS',  # TEE CYS  POST-TRANSLATIONAL MODIFICATION
    'TH5': 'THR',  # TH5 THR  O-ACETYL-L-THREONINE
    'THC': 'THR',  # THC THR  N-METHYLCARBONYLTHREONINE
    'THR': 'THR',  # THR THR
    'TIH': 'ALA',  # TIH ALA  BETA(2-THIENYL)ALANINE
    'TMD': 'THR',  # TMD THR  N-METHYLATED, EPSILON C ALKYLATED
    'TNB': 'CYS',  # TNB CYS  S-(2,3,6-TRINITROPHENYL)CYSTEINE
    'TOX': 'TRP',  # TOX TRP
    'TPL': 'TRP',  # TPL TRP  TRYTOPHANOL
    'TPO': 'THR',  # TPO THR  HOSPHOTHREONINE
    'TPQ': 'ALA',  # TPQ ALA  2,4,5-TRIHYDROXYPHENYLALANINE
    'TQQ': 'TRP',  # TQQ TRP
    'TRF': 'TRP',  # TRF TRP  N1-FORMYL-TRYPTOPHAN
    'TRN': 'TRP',  # TRN TRP  AZA-TRYPTOPHAN
    'TRO': 'TRP',  # TRO TRP  2-HYDROXY-TRYPTOPHAN
    'TRP': 'TRP',  # TRP TRP
    'TRQ': 'TRP',  # TRQ TRP
    'TRW': 'TRP',  # TRW TRP
    'TRX': 'TRP',  # TRX TRP  6-HYDROXYTRYPTOPHAN
    'TTQ': 'TRP',  # TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
    'TTS': 'TYR',  # TTS TYR
    'TY2': 'TYR',  # TY2 TYR  3-AMINO-L-TYROSINE
    'TY3': 'TYR',  # TY3 TYR  3-HYDROXY-L-TYROSINE
    'TYB': 'TYR',  # TYB TYR  TYROSINAL
    'TYC': 'TYR',  # TYC TYR  L-TYROSINAMIDE
    'TYI': 'TYR',  # TYI TYR  3,5-DIIODOTYROSINE
    'TYN': 'TYR',  # TYN TYR  ADDUCT AT HYDROXY GROUP
    'TYO': 'TYR',  # TYO TYR
    'TYQ': 'TYR',  # TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
    'TYR': 'TYR',  # TYR TYR
    'TYS': 'TYR',  # TYS TYR  INE SULPHONATED TYROSINE
    'TYT': 'TYR',  # TYT TYR
    'TYX': 'CYS',  # TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
    'TYY': 'TYR',  # TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
    'TYZ': 'ARG',  # TYZ ARG  PARA ACETAMIDO BENZOIC ACID
    'UMA': 'ALA',  # UMA ALA
    'VAD': 'VAL',  # VAD VAL  DEAMINOHYDROXYVALINE
    'VAF': 'VAL',  # VAF VAL  METHYLVALINE
    'VAL': 'VAL',  # VAL VAL
    'VDL': 'VAL',  # VDL VAL  (2R,3R)-2,3-DIAMINOBUTANOIC ACID
    'VLL': 'VAL',  # VLL VAL  (2S)-2,3-DIAMINOBUTANOIC ACID
    'HSD': 'HIS',  # up his
    'VME': 'VAL',  # VME VAL  O- METHYLVALINE
    'X9Q': 'ALA',  # X9Q ALA
    'XX1': 'LYS',  # XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
    'XXY': 'ALA',  # XXY ALA
    'XYG': 'ALA',  # XYG ALA
    'YCM': 'CYS',  # YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
    'YOF': 'TYR'}  # YOF TYR  3-FLUOROTYROSINE

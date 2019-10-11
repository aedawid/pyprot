import sys
from PDB_menager import PDB, read_pdb
from Structure import Structure, Atom
from SEQ_menager import read_conservation

# All-atom section
pdb = PDB(read_pdb('2JY6.pdb', '1', 'B', []))
structure = Structure(pdb.structure)

feature = read_conservation('clustal.aln', 'skylign.txt')

structure.place_seq_feature_in_bfcol(feature)
for atom in structure.structure:
    print(atom.__str__())

#print(len(structure.structure))
#print("CHAINS: ", structure.chains)
#print("RESIDUES: ", structure.residues)
#print("SEQ:", structure.sequence)
#print("SEQ_LEN:", structure.seq_len)

#United-atom section
#united = structure.create_united_representation(2, 2)

#print("U_CODES: ", structure.united_codes)

#for atom in united.structure:
#    print(atom.__str__())

#UN_ATTRIBUTES = united.united_attributes
#for key in UN_ATTRIBUTES:
#    print(key, UN_ATTRIBUTES[key])

#united.place_attribute_in_bfcol()
#for atom in united.structure:
#    print(atom.__str__())


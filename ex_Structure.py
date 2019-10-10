import sys
from PDB_menager import PDB, read_pdb
from Structure import Structure, Atom

pdb = PDB(read_pdb('2JY6.pdb', '1', 'A', []))
structure = Structure(pdb.structure)
#print(len(structure.structure))

#print("CHAINS: ", structure.chains)
#print("RESIDUES: ", structure.residues)
#print("SEQ:", structure.sequence)
#print("SEQ_LEN:", structure.seq_len)

united = structure.create_united_representation(2, 2)

#print("U_CODES: ", structure.united_codes)

for atom in united:
    print(atom.__str__())

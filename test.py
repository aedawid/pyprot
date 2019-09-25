#import manage_pdb
from manage_pdb import read_pdb, read_pdb_model, write_pdb, write_pdb_model

structure = []
#read_pdb('1D3Z_UB_NMR.pdb', structure)
atoms = [' N  ', ' CA ']
read_pdb_model('ubiquitin.pdb', '3', 'A', atoms, structure)
write_pdb(structure)
#write_pdb_model(structure, 11)
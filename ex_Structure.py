#import manage_pdb
from manage_pdb import read_pdb, read_pdb_model, write_pdb, write_pdb_model
from manage_pdb import pdb_to_fasta, AA_code, read_fasta
from classes import AA, Chain

#structure = []
#read_pdb('1D3Z_UB_NMR.pdb', structure)
#atoms = [' N  ', ' CA ']
#read_pdb_model('ubiquitin.pdb', '3', 'A', atoms, structure)
#write_pdb(structure)
#write_pdb_model(structure, 11)

#print(pdb_to_fasta('ubiquitin.pdb'))
#print(AA_code('A'))
#print(AA_code('ALA'))
#print(read_fasta('ubiquitin.fasta'))
seq = read_fasta('ubiquitin.fasta')
#i = AA('H')
#print(i.name, i.id, i.weight)
ch = Chain(seq)
#print(ch.name, ch.id, ch.weight)
b = ch.beads[4]
print(b.name, b.id, b.weight)
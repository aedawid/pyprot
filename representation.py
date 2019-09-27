from classes import AA, Chain, Pdb, read_pdb, read_fasta

s = Pdb(read_pdb('6MUN.pdb', '1', '0', []))
#ss = Pdb(s.aa_feature(3))
#ss.write_pdb()
#print(s.atom_ids())
print(s.seq())
#print(read_fasta('UBQLN2.fasta'))

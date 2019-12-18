import sys
from PDB_menager import PDB, read_pdb, AA_code


#print(AA_code('ALA'))
#print(AA_code('K'))
#print(AA_code('1AB'))
#print(AA_code('KKK'))


s = PDB(read_pdb('1YYB.pdb'))
print(s.pdb_seq[0])
#print(s.seq_ranges())
#print(s.missing_resi())
#print(s.outcome_seq())
#print(s.ss_for_seq())
#s.write_pdb()

#print(s.atom_ids)
#print(s.non_standard_resi())
#print(s.seq_from_struct())
#print(s.ss_for_struct())

#print(s.resi_names())
#print(s.resi_ids())
#print(len(s.seq_from_header()))
#print(len(s.seq_from_struct()))
#print(len(s.resi_ids()))
#print(len(s.chain_ids()))
#for k in s.header['SEQRES']:
#    print(k)


#print(s.outcome_seq())
#print(s.ss_for_seq())
#print(s.header['REMARK 465'])
#print(s.header['SEQRES'])
#print(len(s.header['HELIX']))
#print(s.header['HELIX'])
#print(len(s.header['SHEET']))
#print(s.header['SHEET'])

#print(s.chain_ids())
#print(s.chain_ranges())
#




#ss = Pdb(s.aa_feature(3))
#print(s.atom_ids())
#print(read_fasta('UBQLN2.fasta'))
##### 0 - id, 1 - weight, 2 - frequency, 3 - charge, 4 - polar, 5 - aromatic, 6 - hp_KD
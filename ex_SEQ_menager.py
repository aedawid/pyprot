import sys
from SEQ_menager import read_LLPSDB, search_motif, aa_composition


db = read_LLPSDB('P_PRT')
print(db['UBQLN2'][4])
seq=db['UBQLN2'][4]

#with open(query_list,'r') as f:
#    queries = [l.strip() for l in f]

#queries=['LQSQM', 'QQQL', '\wQQL', '\wQQ\w', 'AAAAA']
#print(search_motif(queries, seq))

composition=aa_composition(seq)
for key in composition:
    print(key, "%.3f" %composition[key])

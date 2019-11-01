import sys
from SEQ_menager import read_LLPSDB, search_motif, aa_composition

db_keys=['NATUR', 'C_ONE', 'C_TWO', 'C_MOR', 'P_DNA', 'P_RNA', 'P_PRT']#'DESIG', 
db={}
pnames=[]
for key in db_keys:
  db.clear()
  db = read_LLPSDB(key)
  for protein in db.keys():
    if protein in pnames:
        continue
    else:
        pnames.append(protein)
        seq=db[protein][4]
        length=db[protein][5]
        if len(seq) > 0:
            print("%20s" %protein, seq)






#  if "ADF3" in db.keys():
#    print(db['ADF3'][4])
#    break
#seq=db['UBQLN2'][4]

#with open(query_list,'r') as f:
#    queries = [l.strip() for l in f]

#queries=['LQSQM', 'QQQL', '\wQQL', '\wQQ\w', 'AAAAA']
#print(search_motif(queries, seq))

#composition=aa_composition(seq)
#for key in composition:
#    print(key, "%.3f" %composition[key])

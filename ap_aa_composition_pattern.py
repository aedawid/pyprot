import re
from SEQ_menager import read_LLPSDB, aa_composition, aa_pattern

db_keys=['NATUR', 'C_ONE', 'C_TWO', 'C_MOR', 'P_DNA', 'P_RNA', 'P_PRT']
pnames=[]
out1=open("aa_counts.txt", 'a')
out1.write("#            protein   "+"   A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y\n")

for db_key in db_keys:
    db = read_LLPSDB(db_key)
    GLOBAL_motifs={}
    out1=open("aa_counts.txt", 'a')
    out2=open("aa_hydrophobicity.txt", 'a')
    out3=open("aa_polarity.txt", 'a')
    out4=open("aa_aromaticity.txt", 'a')
    out5=open("aa_charge.txt", 'a')
    for protein in db.keys():
        if protein in pnames:
            continue
        else:
            pnames.append(protein)
            seq=db[protein][4]
            length=db[protein][5]
            if len(seq) > 0:
                aa_counts=aa_composition(seq)
                aa_compos=aa_pattern(seq)
                line=(str(protein).rjust(20))
                for key in aa_counts:
                    line+=(" "+str(round(aa_counts[key], 3)).rjust(6))
                out1.write(line+"\n")
                out2.write(str(protein).rjust(20)+" "+str(aa_compos['hp'][0])+" "+str(round(aa_compos['hp'][1], 3)).rjust(6)+"\n")
                out3.write(str(protein).rjust(20)+" "+str(aa_compos['pol'][0])+" "+str(round(aa_compos['pol'][1], 3)).rjust(6)+"\n")
                out4.write(str(protein).rjust(20)+" "+str(aa_compos['pi'][0])+" "+str(round(aa_compos['pi'][1], 3)).rjust(6)+" "+str(aa_compos['pi'][0].count("1")).rjust(6)+"\n")
                out5.write(str(protein).rjust(20)+" "+str(aa_compos['ch'][0])+" "+str(round(aa_compos['ch'][1], 3)).rjust(6)+" "+str(aa_compos['ch'][0].count("1")).rjust(6)+" "+str(aa_compos['ch'][0].count("2")).rjust(6)+"\n")
    out1.close()
    out2.close()
    out3.close()
    out4.close()
    out5.close()
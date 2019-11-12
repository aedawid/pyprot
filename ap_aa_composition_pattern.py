import re
from SEQ_menager import read_LLPSDB, aa_composition, aa_pattern, sequence_charge_decoration, overall_charge_symmetry

pnames=[]
out1=open("aa_counts.txt", 'a')
out1.write("#            protein   "+"   A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y\n")
out6=open("ch_decoration.txt", 'a')
out6.write("# protein_name".rjust(20)+" "+"SCD".rjust(9)+" "+"OCS".rjust(9)+" "+"FCR".rjust(9)+"\n")

with open("llps_sequences",'r') as f:
    for row in f:
        tokens=row.split()
        if len(tokens) > 0:
            out1=open("aa_counts.txt", 'a')
            out2=open("aa_hydrophobicity.txt", 'a')
            out3=open("aa_polarity.txt", 'a')
            out4=open("aa_aromaticity.txt", 'a')
            out5=open("aa_charge.txt", 'a')
            out6=open("ch_decoration.txt", 'a')
            seq=tokens[1]
            length=len(seq)
            if len(seq) > 0:
                aa_counts=aa_composition(seq)
                aa_compos=aa_pattern(seq)
                line=(str(tokens[0]).rjust(20))
                for key in aa_counts:
                    line+=(" "+str(round(aa_counts[key], 3)).rjust(6))
                ocs=overall_charge_symmetry(seq)
                out1.write(line+"\n")
                out2.write(str(tokens[0]).rjust(20)+" "+str(aa_compos['hp'][0])+" "+str(round(aa_compos['hp'][1], 3)).rjust(6)+"\n")
                out3.write(str(tokens[0]).rjust(20)+" "+str(aa_compos['pol'][0])+" "+str(round(aa_compos['pol'][1], 3)).rjust(6)+"\n")
                out4.write(str(tokens[0]).rjust(20)+" "+str(aa_compos['pi'][0])+" "+str(round(aa_compos['pi'][1], 3)).rjust(6)+" "+str(aa_compos['pi'][0].count("1")).rjust(6)+"\n")
                out5.write(str(tokens[0]).rjust(20)+" "+str(aa_compos['ch'][0])+" "+str(round(aa_compos['ch'][1], 3)).rjust(6)+" "+str(aa_compos['ch'][0].count("1")).rjust(6)+" "+str(aa_compos['ch'][0].count("2")).rjust(6)+"\n")
                out6.write(str(tokens[0]).rjust(20)+" "+str(round(sequence_charge_decoration(seq), 4)).ljust(6, '0').rjust(9)+" "+str(round(ocs[1], 4)).ljust(6, '0').rjust(9)+" "+str(round(ocs[0], 4)).ljust(6, '0').rjust(9)+"\n")
    out1.close()
    out2.close()
    out3.close()
    out4.close()
    out5.close()
    out6.close()




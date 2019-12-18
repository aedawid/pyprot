import sys
import os
import re
import argparse
sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')
from SEQ_menager import search_motif
#print("USAGE:\npython ap_parse_disorder.py -f file.glob")

def ss(filename):
  prefix=filename.split('_')[0]
  with open(filename,'r') as f:
    for row in f:
      tokens=row.strip().split()
      if len(tokens) == 2:
        protein=tokens[0]
        seq=tokens[1]
        length=len(seq)
        out0=open(prefix+"_counts.txt", 'a')
        out1=open(prefix+"_charged.txt", 'a')
        out2=open(prefix+"_lcr.txt", 'a')
###---Rich regions (a - excess of the AA or b - chraged residues) in 20-AA fragment
        n=0
        shift = 20
        new_motifs=[]
        aa=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        while n < len(seq)-shift-1:
          new_motifs.append(seq[n:n+shift])
          n+=1
        for frag in new_motifs:
          position=seq.find(frag)
          for a in aa:
            n=frag.count(a)
            out0.write(str(protein.rjust(10))+" "+a+" "+str(n).rjust(2)+"\n")
          nD=frag.count('D')
          nE=frag.count('E')
          nK=frag.count('K')
          nR=frag.count('R')
          nH=frag.count('H')
          NN=nD+nE
          NP=nK+nR+nH
          ALL=nD+nE+nK+nR+nH
          out1.write(str(protein.rjust(10))+" DEKRH "+str(ALL).rjust(2)+" "+str(NN).rjust(2)+" "+str(NP).rjust(2)+" "+str(nD).rjust(2)+" "+str(nE).rjust(2)+" "+str(nK).rjust(2)+" "+str(nR).rjust(2)+" "+str(nH).rjust(2)+"\n")
###---Low complexity when 20-aa fragment is composed of maximally 5 AA or at least 2 AA occure more than 4
          s=''.join(set(frag))
          ni=0
          nj=0
          ai=''
          aj=''
          for i in s:
            nk=frag.count(i)
            if nk >= ni:
              nj=ni
              aj=ai
              ni=nk
              ai=i
            elif nk > nj:
              nj=nk
              aj=i
          out2.write(str(protein.rjust(10))+" "+s.rjust(20)+" "+str(len(s)).rjust(2)+" "+ai+" "+str(ni).rjust(2)+" "+aj+" "+str(nj).rjust(2)+"\n")

      out0.close()
      out1.close()
      out2.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='wget',
        description=""" """,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=''
    )
    parser.add_argument(
        '-f', '--filename',
        help='input ss file',
        metavar='INPUT',
        dest='inputs',
        required=True
    )
    args = parser.parse_args()
    ss(args.inputs)



    
import sys
import os
import re
import argparse
sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')
from SEQ_menager import search_motif
#print("USAGE:\npython ap_parse_disorder.py -f file.glob")

def ss(filename):
  KNOWN_motifs={}
  GLOBAL_motifs={}
  with open(filename,'r') as f:
    for row in f:
      tokens=row.strip().split()
      if len(tokens) == 2:
        protein=tokens[0]
        seq=tokens[1]
        length=len(seq)
        out1=open("rich_regions.txt", 'a')
        out2=open("charged.txt", 'a')
        out3=open("lcr.txt", 'a')
        out4=open("known_motifs.txt", 'a')
        out5=open("new_motifs.txt", 'a')
###---Known motifs (like LARKS, RGG, etc.)
        known_motifs=['S\w\wE', 'S\w\wD', 'S\wR', 'G\w\wG\w\w', 'SYSGYS', 'SYSSYGQS', 'GYNGFG', 'STGGYG', 'GFGNFGTS', 'NKGAII']
        FOUND_motifs=search_motif(known_motifs, seq)
        for motif in FOUND_motifs:
          vec=[]
          vec.append(protein)
          vec.append(length)
          for i in FOUND_motifs[motif]:
            vec.append(i)
          if not motif in KNOWN_motifs.keys():
            KNOWN_motifs[motif]=[]
          KNOWN_motifs[motif].append(vec)
###---New motifs of length from 2 to 20
        new_motifs=[]
        for shift in range(2, 20):
          n=0
          while n < len(seq)-shift-1:
            new_motifs.append(seq[n:n+shift])
            n+=1
        FOUND_motifs=search_motif(new_motifs, seq)
#        FOUND_motifs={}
#        for motif in motifs:
#          if motifs[motif][0] > 6 and len(motif) < 4:
#            FOUND_motifs[motif] = motifs[motif]
#          elif motifs[motif][0] > 2 and len(motif) >= 4:
#            FOUND_motifs[motif] = motifs[motif]
        for motif in FOUND_motifs:
          if not motif in known_motifs:
            vec=[]
            vec.append(protein)
            vec.append(length)
            for i in FOUND_motifs[motif]:
              vec.append(i)
            if not motif in GLOBAL_motifs.keys():
              GLOBAL_motifs[motif]=[]
            GLOBAL_motifs[motif].append(vec)
###---Rich regions (a - excess of the AA or b - chraged residues) in 20-AA fragment
        n=0
        shift = 20
        new_motifs=[]
        aa=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        while n < len(seq)-shift-1:
          new_motifs.append(seq[n:n+shift])
          n+=5
        for frag in new_motifs:
          for a in aa:
            n=frag.count(a)
            if n >= 5:
              position=seq.find(frag)
              out1.write(str(protein.rjust(10))+" "+a+" "+str(n).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
          nD=frag.count('D')
          nE=frag.count('E')
          nK=frag.count('K')
          nR=frag.count('R')
          position=seq.find(frag)
          if nD+nE >=5:
            out2.write(str(protein.rjust(10))+"   DE "+str(nD+nE).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
          if nK+nR >=5:
            out2.write(str(protein.rjust(10))+"   KR "+str(nK+nR).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
          if nD+nE+nK+nR >=5:
            out2.write(str(protein.rjust(10))+" DEKR "+str(nD+nE+nK+nR).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
###---Low complexity when 20-aa fragment is composed of maximally 5 AA or at least 2 AA occure more than 4
          s=''.join(set(frag))
          if len(s) <= 5:
            out3.write(str(protein.rjust(10))+" 1 "+s.rjust(5)+" "+frag+str(position).rjust(4)+"\n")
          ni=0
          ss=''
          for i in s:
            if frag.count(i) >= 4:
              ni+=1
              ss+=i
          if ni >= 2:
            out3.write(str(protein.rjust(10))+" 2 "+ss.rjust(5)+" "+frag+str(position).rjust(4)+"\n")
###---Known and new motifs found several times in a single sequence or once (or more) in many sequences
      for typ in [KNOWN_motifs, GLOBAL_motifs]:
        for num, motif in enumerate(typ):
          for p in typ[motif]:
            string=''
            string+=str(p[0]).rjust(10)+" "+str(p[1]).rjust(4)+" "+str(p[2]).rjust(3)
            for i in range(3, len(p)):
              string+=" "+str(p[i])
            if typ == KNOWN_motifs:
              out4.write(str("MOTIF "+str(num).rjust(10)+" "+str(motif).rjust(20)+" "+string+"\n"))
            elif typ == GLOBAL_motifs:
              out5.write(str("MOTIF "+str(num).rjust(10)+" "+str(motif).rjust(20)+" "+string+"\n"))
      out1.close()
      out2.close()
      out3.close()
      out4.close()
      out5.close()

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



    
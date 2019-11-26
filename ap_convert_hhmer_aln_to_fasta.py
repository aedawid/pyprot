import sys
import os
import argparse
#sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

#print("USAGE:\npython ap_convert_hhmer_aln_to_fasta.py")

def ss(filename):
  prefix=str(filename).split("-")[0]
  aln={}
  with open(filename,'r') as f:
    for nn,row in enumerate(f):
      tokens=row.strip().split()
      if len(tokens) == 2:
        if not tokens[0].startswith("#"):
          if tokens[0] in aln.keys():
            aln[tokens[0]]+=str(tokens[1])
          else:
            aln[tokens[0]]=str(tokens[1])
  out=open(prefix+"."+str(len(aln))+".fasta.aln", 'a')
  for key in aln:
    out.write(str(aln[key]).replace(".", "-").upper()+"\n")
  out.close()

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

import sys
import os
import argparse
#sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

#print("USAGE:\npython ap_parse_disorder.py -f file.glob")

def ss(filename):
  prefix=str(filename).split(".")[0]
  typ=str(filename).split(".")[1]
  iup=anc=''
  with open(filename,'r') as f:
    for nn,row in enumerate(f):
      tokens=row.strip().split()
      if len(tokens):
        if not tokens[0].startswith("#"):
          if len(tokens) == 4 and len(tokens[1]) == 1:
            val=float(tokens[2])
            if val < 0.5:
              iup+='0'
            elif val < 0.75:
              iup+='1'
            else:
              iup+='2'
            val=float(tokens[3])
            if val < 0.5:
              anc+='0'
            elif val < 0.75:
              anc+='1'
            else:
              anc+='2'

  dis1=1.0-(iup.count('0')/float(len(iup)))
  dis2=1.0-(anc.count('0')/float(len(anc)))
  out=open(prefix+".DIS", 'a')
  out.write(str(prefix).rjust(10)+" "+typ.rjust(5)+" "+str(round(dis1,2)).ljust(4,'0').rjust(4)+"IUPRED2".rjust(8)+" "+str(iup)+"\n")
  out.write(str(prefix).rjust(10)+" "+typ.rjust(5)+" "+str(round(dis2,2)).ljust(4,'0').rjust(4)+"ANCHOR2".rjust(8)+" "+str(anc)+"\n")
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

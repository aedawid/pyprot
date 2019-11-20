import sys
import os
import argparse
#sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

#print("USAGE:\npython ap_parse_ss.py -f file.ss -c 2 -t PS")

def ss(filename):
  
  typ=str(filename).split(".")[0]
  seq=prefix=''
  fraq=0.0
  n='----'
  with open(filename,'r') as f:
    for nn,row in enumerate(f):
      tokens=row.strip().split()
      prefix=tokens[0]
      seq=tokens[1]
      fraq=round(float(tokens[2]),2)
      if len(tokens) > 3:
        n=str(round(float(tokens[3])/len(seq),2))
      
      out=open(prefix+".pattern", 'a')
      out.write(str(prefix).rjust(9)+" "+str(fraq).ljust(4, '0').rjust(4)+" "+n.ljust(4, '0').rjust(4)+str(typ).rjust(10)+" "+str(seq)+"\n")
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

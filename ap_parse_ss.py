import sys
import os
import argparse
#sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

#print("USAGE:\npython ap_parse_ss.py -f file.ss -c 2 -t PS")

def ss(filename, col=2, tool='SSSSSSSSS'):
  prefix=str(filename).split(".")[0]
  typ=str(filename).split(".")[1]
  ss=''
  with open(filename,'r') as f:
    for nn,row in enumerate(f):
      tokens=row.strip().split()
      if len(tokens):
        if not tokens[0].startswith("#"):
          n=len(tokens)
          if n == 6 or n == 11:		#ss3 or ss8 format file
            ss+=str(tokens[col])
  dis=ss.count('C')/float(len(ss))
  out=open(prefix+".SS", 'a')
  out.write(str(prefix).rjust(10)+" "+typ+" "+str(round(dis,2)).rjust(4)+str(tool).rjust(10)+" "+str(ss)+"\n")
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
    parser.add_argument(
        '-c', '--column',
        help='parsed column',
        metavar='COL',
        type=int,
        default='2',
        dest='col',
    )
    parser.add_argument(
        '-t', '--tool',
        help='used tool',
        metavar='TOOL',
        type=str,
        default='SS',
        dest='tool',
    )
    args = parser.parse_args()
    ss(args.inputs, args.col, args.tool)

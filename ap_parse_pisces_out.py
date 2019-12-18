import sys
import os
import argparse
sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

#print("USAGE:\npython ap_parse_pisces_out.py -f file.glob")

def ss(filename):
  prefix=str(filename).split(".")[0]
  with open(filename,'r') as f:
    for row in f:
      if row.startswith("SEQ_H:"):
        seq=row.strip().split(':')[2]
  out=open("pisces_sequences", 'a')
  out.write(str(prefix).rjust(6)+" "+str(seq)+"\n")
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

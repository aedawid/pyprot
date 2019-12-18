import sys
import os
import argparse
#sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

#print("USAGE:\npython ap_parse_disorder.py -f file.glob")

def ss(filename):
  prefix=str(filename).split(".")[0]
  suma=0.0
  data=[]
  with open(filename,'r') as f:
    for nn,row in enumerate(f):
      tokens=row.strip().split()
      data.append(tokens)
      if len(tokens) == 4:
        suma+=float(tokens[1])
  for i in data:
    out=open(prefix+".h", 'a')
    out.write(str(i[0]).rjust(6)+" "+str(round(float(i[1])/suma,2)).ljust(4, '0')+" "+str(round(float(i[2])/suma,2)).ljust(4, '0')+" "+str(round(float(i[3])/suma,2)).ljust(4, '0')+" "+i[1].rjust(2)+" "+i[2].rjust(2)+" "+i[3].rjust(2)+"\n")
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

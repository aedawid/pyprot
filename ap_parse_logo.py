import sys
import os
import argparse
sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')


def logo(filename):
  max_h=ave=max_v=0.0
  high=[]
  pattern=''
  with open(filename,'r') as f:
    prefix=str(filename).split(".")[0]
    for row in f:
        row=row.replace("(", "")
        row=row.replace(")", "")
        tokens=row.split()
        if row.startswith("max expected height"):
            max_h=float(row.split()[4])
        elif len(tokens) == 22:
            high.append(float(tokens[21]))
            if float(tokens[21]) > max_v:
               max_v=float(tokens[21])
            ave+=float(tokens[21])

  ave=ave/len(high)
  diff=(max_v-ave)/5.0
  print(max_h, max_v, ave, diff)#####################
  print(high)########################################
  for i in high:
    if i <= ave:
        pattern+="0"
    elif i <= ave+diff:
        pattern+="1"
    elif i <= ave+(2*diff):
        pattern+="2"
    elif i <= ave+(3*diff):
        pattern+="3"
    elif i <= ave+(4*diff):
        pattern+="4"
    else:
        pattern+="5"

  out=open("conservation.txt", 'a')
  out.write(str(prefix).rjust(20)+" "+str(pattern)+"\n")
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
        help='input logo file',
        metavar='LOGO',
        dest='logo',
        required=True
    )

    args = parser.parse_args()
    logo(args.logo)
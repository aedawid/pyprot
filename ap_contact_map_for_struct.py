import sys
import argparse

from PDB_menager import PDB, read_pdb
from Structure import Structure, Atom
from SEQ_menager import read_conservation

sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

def code(filename, model='1', chain='A', cutoff=8.0, output=None):

  pdb = PDB(read_pdb(filename, model, chain, [' CA ']))
  structure = Structure(pdb.structure)
  subset = structure.select_atoms('A')
  
#  for atom in subset:
#    print(atom.__str__())
  
  if output:
    outfile=open(output, 'a')
  else:
    outfile = sys.stdout

  for at1 in subset:
    for at2 in subset:
#      if at1.idx != at2.idx:
        idx = 0
        dist=at1.distance(at2)
        if dist <= cutoff:
          idx = 1
        outfile.write(str(at1.resi_id).rjust(4)+" "+str(at2.resi_id).rjust(4)+" "+str(round(dist,3)).rjust(7)+" "+str(idx).rjust(2)+"\n")
  outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='code',
        description="""Application for calculating contact maps for given structure.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=''
    )
    parser.add_argument(
        '-f', '--filename',
        help='input pdb file',
        metavar='PDB',
        dest='pdb',
        required=True
    )
    parser.add_argument(
        '-m', '--model-id',
        help='chosen model',
        metavar='MODEL',
        type=str,
        default='1',
        dest='model'
    )
    parser.add_argument(
        '-c', '--chain-id',
        help='chosen chain',
        metavar='CHAIN',
        type=str,
        default='0',
        dest='chain'
    )
    parser.add_argument(
        '-d', '--distance',
        help='contact distance cutoff',
        metavar='DIST',
        type=float,
        default=8.0,
        dest='dist'
    )
    parser.add_argument(
        '-o', '--output-file',
        help='save output file [default stdout]',
        default=None,
        metavar='OUTPUT',
        dest='out'
    )
    args = parser.parse_args()
    code(args.pdb, args.model, args.chain, args.dist, args.out)

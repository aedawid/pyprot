import os
import sys
import argparse

from classes import AA, Chain, Pdb, read_pdb, read_fasta

def wget(filename, output=None, model='1', chain='0'):

    s = Pdb(read_pdb(filename, model, chain, []))
    s.write_pdb()
    if output:
        outfile = open(output, 'w')
    else:
        outfile = sys.stdout


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='wget',
        description=""" """,
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
        '-o', '--output-pdb',
        help='save output pdb [default stdout]',
        default=None,
        metavar='OUTPUT',
        dest='out'
    )
    parser.add_argument(
        '-m', '--model-id',
        help='chosen model',
        metavar='MODEL',
        type=str,
        default='0',
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
    args = parser.parse_args()
    wget(args.pdb, args.model, args.chain)
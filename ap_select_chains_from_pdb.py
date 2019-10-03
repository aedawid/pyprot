import sys
import argparse
from PDB_menager import PDB, read_pdb


def ap_select_chains_from_pdb(filename, output=None, model='0', chain='0'):

    s = PDB(read_pdb(filename, model, chain, []))
    s.write_pdb()
    if output:
        outfile = open(output, 'w')
    else:
        outfile = sys.stdout
    
    outfile.write("\nMISSING,")
    for i in s.missing_resi():
        if i[1] == chain:
            st1=i[0]+" "+i[2]+" "+i[1]+","
            outfile.write(st1)
    
    outfile.write("\nNON_AA,")
    for i in s.non_standard_resi():
        if i[3] == chain:
            st2=i[0]+" "+i[2]+" "+i[3]+" "+i[1]+","
        outfile.write(st2)
    
    c_range = s.chain_resi(chain)
    st3=''
    st4=''
    st5=''
    for i in range(c_range[0], c_range[len(c_range)-1]):
        st3 += list(s.seq_from_header())[i]
        st4 += list(s.outcome_seq())[i]
        st5 += list(s.ss_for_seq())[i]
    outfile.write("\nSEQ_H,")
    outfile.write(st3)
    outfile.write("\nSEQ_O,")
    outfile.write(st4)
    outfile.write("\nSS   ,")
    outfile.write(st5)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ap_select_chains_from_pdb',
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
    ap_select_chains_from_pdb(args.pdb, args.out, args.model, args.chain)
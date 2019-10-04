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
    
    outfile.write(" MISS:")
    for i in s.missing_resi():
        if i[1] == chain:
            st1=i[0]+" "+i[2]+" "+i[1]+","
            outfile.write(st1)
    
    outfile.write("\n\nNON_A:")
    for i in s.non_standard_resi():
        if i[3] == chain:
            st2=i[0]+" "+i[2]+" "+i[3]+" "+i[1]+","
        outfile.write(st2)
    
    outfile.write("\n\nHETAT:")
    for ii in s.header["HET "]:
        i = ii.split()
        st = i[0]+" "+i[1]+" "+i[2]+" "+i[3]+","
        outfile.write(st)
    
    outfile.write("\n\nCIS_P:")
    for ii in s.header["CISPEP"]:
        i = ii.split()
        st = i[0]+" "+i[1]+" "+i[2]+" "+i[3]+" "+i[4]+" "+i[5]+","
        outfile.write(st)
    
    c_range = s.chain_resi(chain)
    st3 = s.seq_from_header()[c_range[0]:c_range[1]]
    st4 = s.outcome_seq()
    st5 = s.ss_for_seq()
    outfile.write("\n\nSEQ_H:")
    outfile.write(st3)
    outfile.write("\n\nSEQ_O:")
    outfile.write(st4)
    outfile.write("\n\n   SS:")
    outfile.write(st5)
    outfile.write("\n\n SIZE:")
    outfile.write(str(len(st4)))
    outfile.write("\n\n  RES:")
    outfile.write(s.header["REMARK   2"][1].split()[1])
    outfile.write("\n\n CATH:")
    outfile.write(s.ss_type())
    outfile.write("\n\nTITLE:")
    outfile.write(s.header["TITLE"][0])
    outfile.write("\n\nSITES:")
    outfile.write(s.binding_sites())
    outfile.write("\n\n SITE:")
    sites = s.binding_resi()
    for site in sites:
        outfile.write(site[0])
        for res in site[1]:
            outfile.write(","+res[0]+" "+res[1]+" "+res[2])
        outfile.write(";")

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
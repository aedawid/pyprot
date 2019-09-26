from classes import AA, Chain, Pdb

s = Pdb('ubiquitin.pdb', 1, 0, ' CA ')
ss = s.read()
for m in ss:
    for a in m:
        print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s%6s%6s           %3s"%tuple(a))
    print("ENDMDL")


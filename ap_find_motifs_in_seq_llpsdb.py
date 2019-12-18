import re
from SEQ_menager import read_LLPSDB, search_motif

db_keys=['DESIG', 'NATUR', 'C_ONE', 'C_TWO', 'C_MOR', 'P_DNA', 'P_RNA', 'P_PRT']
for db_key in db_keys:
    print(db_key)#####################
    db = read_LLPSDB(db_key)
    GLOBAL_motifs={}
    out1=open(db_key+"_rich_regions.txt", 'a')
    out2=open(db_key+"_charged.txt", 'a')
    out3=open(db_key+"_lcr.txt", 'a')
    output=open(db_key+"_motifs.txt", 'a')
    for protein in db.keys():
        seq=db[protein][4]
        length=db[protein][5]
        known_motifs=['S\w\wE', 'S\w\wD', 'S\wR', 'G\w\wG\w\w', 'SYSGYS', 'SYSSYGQS', 'GYNGFG', 'STGGYG', 'GFGNFGTS', 'NKGAII']
        FOUND_motifs=search_motif(known_motifs, seq)
        for motif in FOUND_motifs:
            vec=[]
            vec.append(protein)
            vec.append(length)
            for i in FOUND_motifs[motif]:
                vec.append(i)
            if not motif in GLOBAL_motifs.keys():
                GLOBAL_motifs[motif]=[]
            GLOBAL_motifs[motif].append(vec)
        
        new_motifs=[]
        for shift in range(2, 10):
            n=0
            while n < len(seq)-shift-1:
                new_motifs.append(seq[n:n+shift])
                n+=1
        motifs=search_motif(new_motifs, seq)
        FOUND_motifs={}
        for motif in motifs:
            if motifs[motif][0] > 6 and len(motif) < 4:
                FOUND_motifs[motif] = motifs[motif]
            elif motifs[motif][0] > 2 and len(motif) >= 4:
                FOUND_motifs[motif] = motifs[motif]
        for motif in FOUND_motifs:
            if not motif in known_motifs:
                vec=[]
                vec.append(protein)
                vec.append(length)
                for i in FOUND_motifs[motif]:
                    vec.append(i)
                if not motif in GLOBAL_motifs.keys():
                    GLOBAL_motifs[motif]=[]
                GLOBAL_motifs[motif].append(vec)
        
        n=0
        shift = 20
        new_motifs=[]
        aa=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        while n < len(seq)-shift-1:
            new_motifs.append(seq[n:n+shift])
            n+=5
        for frag in new_motifs:
            for a in aa:
                n=frag.count(a)
                if n >= 5:
                    position=seq.find(frag)
                    out1.write(str(protein.rjust(10))+" "+a+" "+str(n).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
            nD=frag.count('D')
            nE=frag.count('E')
            nK=frag.count('K')
            nR=frag.count('R')
            position=seq.find(frag)
            if nD+nE >=5:
                out2.write(str(protein.rjust(10))+"   DE "+str(nD+nE).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
            if nK+nR >=5:
                out2.write(str(protein.rjust(10))+"   KR "+str(nK+nR).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
            if nD+nE+nK+nR >=5:
                out2.write(str(protein.rjust(10))+" DEKR "+str(nD+nE+nK+nR).rjust(2)+" "+frag+str(position).rjust(4)+"\n")
            s=''.join(set(frag))
            if len(s) <= 5:
                out3.write(str(protein.rjust(10))+" 1 "+s.rjust(5)+" "+frag+str(position).rjust(4)+"\n")
            ni=0
            ss=''
            for i in s:
                if frag.count(i) >= 4:
                    ni+=1
                    ss+=i
            if ni >= 2:
                out3.write(str(protein.rjust(10))+" 2 "+ss.rjust(5)+" "+frag+str(position).rjust(4)+"\n")
    
    for num, motif in enumerate(GLOBAL_motifs):
        for p in GLOBAL_motifs[motif]:
            string=''
            string+=str(p[0]).rjust(10)+" "+str(p[1]).rjust(4)+" "+str(p[2]).rjust(3)
            for i in range(3, len(p)):
                string+=" "+str(p[i])
            output.write(str("MOTIF "+str(num).rjust(10)+" "+str(motif).rjust(10)+" "+string+"\n"))
    output.close()
    out1.close()
    out2.close()
    out3.close()
    
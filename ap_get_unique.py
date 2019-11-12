import sys

seq=[]
data={}
with open("p_data",'r') as f:
    for row in f:
        tokens=row.split()
        if len(tokens) > 0:
            if not tokens[2] in seq:
                seq.append(tokens[2])
                data[tokens[2]]=[[],[]]
            if not tokens[1] in data[tokens[2]][0]:
                data[tokens[2]][0].append(tokens[1])
            if not tokens[0] in data[tokens[2]][1]:
                data[tokens[2]][1].append(tokens[0])
for key in data.keys():
    print("%20s" %data[key][1], "%10s" %data[key][0], key)


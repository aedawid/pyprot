import sys

with open("PhaSePro",'r') as f:
    for row in f:
        tokens=row.split()
        if len(tokens) > 0:
            print("%20s" %tokens[0], "%10s" %tokens[1], tokens[2])

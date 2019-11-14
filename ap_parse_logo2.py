import sys
import os
import argparse
import operator
sys.path.append('/Users/adawid/PROJECTS/own/IDP_LLPS/UBQL/code/pyprot')

#print("USAGE:\npython ap_parse_logo.py -f file.logo")

def logo(filename, db="SP"):
  prefix=str(filename).split(".")[0]
  
  cons=seq=''
  nums=[]
  with open(str(prefix+".MSA"), 'r') as f1:
    for line in f1:
      msa=list(line)
  with open(str(prefix+".PROFILE"), 'r') as f2:
    for line in f2:
      profile=list(line)
  with open(str(prefix+".NUM"), 'r') as f3:
    for line in f3:
      nums.append(line.strip())
  
  new_profile=['-']*len(msa)
  for idx, t in nums:
    new_profile[int(t)]=profile[idx]
  for idx, t in new_profile:
    if t != '-':
      cons+=t
      if msa[idx] != '.':
        seq+=msa[idx]
      else:
        seq+='-'
    elif msa[idx] != '.':
      cons+='-'
      seq+=msa[idx]

  ave=max_v=0.0
  high=[]
  pattern_s=pattern_d=pattern_x=consensus=''
  with open(filename,'r') as f:
    for nn,row in enumerate(f):
      row=row.replace("(", "")
      row=row.replace(")", "")
      tokens=row.split()
      if len(tokens) == 22:
        n=0
        c=[0.0,'']
        amino=''
        counts={}
        for key in ['hp', 'pol', 'pi', 'ar', 'ch']:
          counts.setdefault(key, 0)
      ###-pattern-strength-###
        high.append(float(tokens[21]))
        if float(tokens[21]) > max_v:
           max_v=float(tokens[21])
        ave+=float(tokens[21])
      ###-pattern-diversity-###
        for num,aa in enumerate(tokens[1:20]):
          if float(aa) != 0.0:
            n+=1				#n - diversity at position [diversity level]
            amino+=AA[num+1]			#amino - aa above treshold at position [diversity type]
            if float(aa) > c[0]:
              c[0] = float(aa)			#c[0] - highest value at position [consensus]
              c[1] = AA[num+1]			#c[1] - highest aa at position [consensus]
        pattern_d+=diversity[n]
        consensus+=c[1]
        for a in amino:
          for key in attribute:
            if a in attribute[key]:
              counts[key]+=amino.count(a)
        if (max(counts.items(), key=operator.itemgetter(1))[1]) == 0:
          s='0'
        else:
          s=max(counts.items(), key=operator.itemgetter(1))[0]
        pattern_x+=feature[s]
        
        out2=open(prefix+".txt", 'a')
        out2.write(str(nn-2).rjust(4)+" "+str(amino).rjust(15)+" "+str(counts).rjust(30)+" "+str(s).rjust(5)+"\n")

  ave=ave/len(high)
  diff=(max_v-ave)/5.0
  for i in high:
    if i <= ave:
        pattern_s+="0"
    elif i <= ave+diff:
        pattern_s+="1"
    elif i <= ave+(2*diff):
        pattern_s+="2"
    elif i <= ave+(3*diff):
        pattern_s+="3"
    elif i <= ave+(4*diff):
        pattern_s+="4"
    else:
        pattern_s+="5"

  if len(cons) != len(pattern_s):
    strength=diversity=attribute=[]
    idx=0
    for i in list(cons):
      if i != '-':
        strength.append(list(pattern_s)[idx])
        diversity.append(list(pattern_d)[idx])
        attribute.append(list(pattern_x)[idx])
        idx+=1
      else:
        strength.append('-')
        diversity.append('-')
        attribute.append('-')
    pattern_s=strength.join()
    pattern_d=diversity.join()
    pattern_x=attribute.join()

  out=open(prefix+".conservation", 'a')
  out.write(str(prefix).rjust(10)+" "+str(db)+"  SEQUENCE "+str(seq)+"\n")
  out.write(str(prefix).rjust(10)+" "+str(db)+" CONSENSUS "+str(cons)+"\n")
  out.write(str(prefix).rjust(10)+" "+str(db)+"  STRENGTH "+str(pattern_s)+"\n")
  out.write(str(prefix).rjust(10)+" "+str(db)+" DIVERSITY "+str(pattern_d)+"\n")
  out.write(str(prefix).rjust(10)+" "+str(db)+" ATTRIBUTE "+str(pattern_x)+"\n")
  out.close()
  out2.close()

feature = {'hp': 'H', 'pol': 'P', 'pi': 'B', 'ar': 'A', 'ch': 'C', '0': '0'}

attribute = {
    'hp': ['V', 'L', 'I', 'M', 'C', 'A'],
   'pol': ['S', 'N', 'T', 'Q', 'Y'],
    'pi': ['R', 'N', 'Q', 'D', 'E'],
    'ar': ['F', 'Y', 'W', 'H'],
    'ch': ['K', 'R', 'H', 'D', 'E']
}

diversity = {
    1: '0', 2: '1', 3: '2', 4: '3', 5: '3',
    6: '4', 7: '4', 8: '4', 9: '4', 10: '4',
    11: '5', 12: '5', 13: '5', 14: '5', 15: '5',
    16: '5', 17: '5', 18: '5', 19: '5', 20: '5'
}

AA = {
    1: 'A', 2: 'C', 3: 'D', 4: 'E', 5: 'F',
    6: 'G', 7: 'H', 8: 'I', 9: 'K', 10: 'L',
    11: 'M', 12: 'N', 13: 'P', 14: 'Q', 15: 'R',
    16: 'S', 17: 'T', 18: 'V', 19: 'W', 20: 'Y'
}

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
    parser.add_argument(
        '-db', '--database',
        help='db used for MSA',
        metavar='DB',
        dest='db',
        required=True
    )

    args = parser.parse_args()
    logo(args.logo, args.db)
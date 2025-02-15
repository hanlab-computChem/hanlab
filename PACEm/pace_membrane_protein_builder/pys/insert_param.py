import sys
f=open(sys.argv[1],"r")

for rl in f:
  if "exclusions" in rl:
    break;

excl=[]

for rl in f:
  if "pairs" in rl:
    break;
  excl.append(rl[:-1])

pair=[]

for rl in f:
  if "impropers" in rl:
    break;
  pair.append(rl[:-1])

imp=[]

for rl in f:
  imp.append(rl[:-1])

f.close()

print ';Adding',len(excl),'excls',len(pair),'pairs',len(imp),'imps'

f=open(sys.argv[2],"r")

for rl in f:
  if 'pairs' in rl:
    break
  print rl[:-1]

print "[ exclusions ]"
for i in excl:
  print i

print rl[:-1]

for i in pair:
  print i

for rl in f:
  print rl[:-1]
  if "dihedrals" in rl:
    break

for i in imp:
  print i

for rl in f:
  print rl,
 

f.close()

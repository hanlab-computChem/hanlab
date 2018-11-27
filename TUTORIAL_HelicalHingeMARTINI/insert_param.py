import sys


fitp=open(sys.argv[1],'r')
fnewp=open(sys.argv[2],'r')
fname=sys.argv[1]

bbp=[]
cp=[]
bap=[]
bdpk=[]
bdpd=[]
for rl in fnewp:
    srl = rl.split()[:]
    if len(srl) == 7 and len(srl[-1]) == 13:
        bbp.append(rl)
    elif len(srl) == 6 and len(srl[-1]) == 13:
        cp.append(rl)
    elif len(srl) == 8 and len(srl[-1]) == 20:
        bap.append(rl)
    elif len(srl) == 10 and len(srl[-1]) == 27:
        if "(T)" not in rl:
            bdpk.append(rl)
        else: 
            bdpd.append(rl)

lines=[]
for rl in fitp:
    srl = rl.split()[:]
    if 'Sidechain bonds' in rl:
        lines += bbp   
    if len(srl) == 7 and len(srl[-1]) == 13:
        for i in range(len(bbp)):
            if srl[0:2] == bbp[i].split()[0:2]:
                rl = ";" + rl 
    if len(srl) == 6 and len(srl[-1]) == 13:
        for i in range(len(bbp)):
            if srl[0:2] == bbp[i].split()[0:2]:
                rl = ";" + rl 
        for i in range(len(cp)):
            if srl[0:2] == cp[i].split()[0:2]:
                rl = cp[i] 
    if len(srl) == 8 and len(srl[-1]) == 20:
        for i in range(len(bap)):
            if srl[0:3] == bap[i].split()[0:3]:
                rl = bap[i] 
    if len(srl) == 10 and len(srl[-1]) == 27:
        for i in range(len(bdpk)):
            if srl[0:4] == bdpk[i].split()[0:4]:
                rl = bdpk[i] 
    if 'Sidechain improper dihedrals' in rl:
        lines += bdpd
    lines.append(rl)

fw = open('new_'+fname,'w')

for i in range(len(lines)):
#    print lines[i]
    fw.write(lines[i])    


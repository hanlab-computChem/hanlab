import sys

atom=[]
backbone=[]
sidechainbone=[]
constraintsB=[]
constraintsS=[]
angleB=[]
angleS=[]
angleBS=[]
dihedralB=[]
dihedralS=[]
posit=[]

def itpRead(rl):
   srl = rl.split()[:]
   if len(srl) == 9 and srl[0] == srl[5] and str(srl[7]) == str(';'):
      srl[0]=int(srl[0])+atom_num
      srl[2]=int(srl[2])+resid_num
      srl[5]=int(srl[5])+atom_num
      atom.append(srl)
   elif len(srl) == 7 and int(srl[2]) == 1 and str(srl[5]) == str(';'):
      srl[0]=int(srl[0])+atom_num
      srl[1]=int(srl[1])+atom_num
      if len(srl[6]) == 13:
         backbone.append(srl)
      elif len(srl[6]) == 3:
         sidechainbone.append(srl)
   elif len(srl) == 6 and str(srl[4]) == str(';'):
      srl[0]=int(srl[0])+atom_num
      srl[1]=int(srl[1])+atom_num
      if len(srl[5]) == 13:
         constraintsB.append(srl)
      elif len(srl[5]) == 3:
         constraintsS.append(srl)
   elif len(srl) == 8 and int(srl[3]) == 2 and str(srl[6]) == str(';'):
      srl[0]=int(srl[0])+atom_num
      srl[1]=int(srl[1])+atom_num
      srl[2]=int(srl[2])+atom_num
      if len(srl[7]) == 20:
         angleB.append(srl)
      elif len(srl[7]) == 3:
         angleS.append(srl)
   elif len(srl) == 9 and str(srl[8]) == str('SBB'):
      srl[0]=int(srl[0])+atom_num
      srl[1]=int(srl[1])+atom_num
      srl[2]=int(srl[2])+atom_num
      angleBS.append(srl)
   elif len(srl) == 10 and int(srl[4]) == 1 and str(srl[8]) == str(';'):
      srl[0]=int(srl[0])+atom_num
      srl[1]=int(srl[1])+atom_num
      srl[2]=int(srl[2])+atom_num
      srl[3]=int(srl[3])+atom_num
      dihedralB.append(srl)
   elif len(srl) == 9 and str(srl[7]) == str(';'):
      srl[0]=int(srl[0])+atom_num
      srl[1]=int(srl[1])+atom_num
      srl[2]=int(srl[2])+atom_num
      srl[3]=int(srl[3])+atom_num
      dihedralS.append(srl)
   elif len(srl) == 5 and int(srl[1]) == 1 and str(srl[2]) == str('POSRES_FC'):
      srl[0]=int(srl[0])+atom_num
      posit.append(srl)
   return

file_num=sys.argv[1]
atom_num=0
resid_num=0
seq=''
ss=''
itpname=''
for i in range(int(file_num)):
    fitp=open(sys.argv[i+2])
    itpname+=sys.argv[i+2]+' '
    rl_num=0
    for rl in fitp:
       itpRead(rl)
       rl_num+=1
       if rl_num==5:
          seq+=rl.split()[-1]
       elif rl_num==7:
          ss+=rl.split()[-1]
    atom_num=int(atom[-1][0])
    resid_num=int(atom[-1][2])

atomline='%5s%6s%6s%6s%6s%6s%8s%2s %-2s'
bondline='%5s%6s%7s%10s%6s%2s %-14s'
bondCline='%5s%6s%7s%10s%2s %-14s'
angleline='%5s%6s%6s%7s%7s%7s%2s %-21s'
angleBSline='%5s%6s%6s%7s%7s%6s%2s%14s %-4s'
angleSSline='%5s%6s%6s%7s%7s%6s%2s %-4s'
dihedralline='%5s%6s%6s%6s%7s%8s%7s%6s%2s %-28s'
dihedralSline='%5s%6s%6s%6s%7s%7s%6s%2s %-4s'
positline='%7s%5s%13s%13s%13s'

f=open('Protein_total.itp','w')

f.write('; MARTINI (martini22) Coarse Grained topology file for %s\n'%itpname)

f.write('; Sequence:\n')
f.write('; %s\n'%seq)
f.write('; Secondary Structure:\n')
f.write('; %s\n'%ss)
f.write('\n')
f.write('[ moleculetype ]\n')
f.write('; Name         Exclusions\n')
f.write('Protein_A            1\n')
f.write('\n')

### writing atoms
f.write('[ atoms ]\n')
for i in range(len(atom)):
    line=atom[i]
    f.write(atomline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8]))
    f.write('\n')

### writing bonds
f.write('\n')
f.write('[ bonds ]\n')
f.write('; Backbone bonds\n')
for i in range(len(backbone)):
    line=backbone[i]
    f.write(bondline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6]))
    f.write('\n')
f.write('; Sidechain bonds\n')
for i in range(len(sidechainbone)):
    line=sidechainbone[i]
    f.write(bondline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6]))
    f.write('\n')
f.write('\n')
f.write('[ constraints ]\n')
for i in range(len(constraintsB)):
    line=constraintsB[i]
    f.write(bondCline%(line[0],line[1],line[2],line[3],line[4],line[5]))
    f.write('\n')
for i in range(len(constraintsS)):
    line=constraintsS[i]
    f.write(bondCline%(line[0],line[1],line[2],line[3],line[4],line[5]))
    f.write('\n')

### writing angles
f.write('\n')
f.write('[ angles ]\n')
f.write('; Backbone angles\n')
for i in range(len(angleB)):
    line=angleB[i]
    f.write(angleline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7]))
    f.write('\n')
f.write('; Backbone-sidechain angles\n')
for i in range(len(angleBS)):
    line=angleBS[i]
    f.write(angleBSline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8]))
    f.write('\n')
f.write('; Sidechain angles\n')
for i in range(len(angleS)):
    line=angleS[i]
    f.write(angleSSline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7]))
    f.write('\n')

### writing dihedrals
f.write('\n')
f.write('[ dihedrals ]\n')
f.write('; Backbone dihedrals\n')
for i in range(len(dihedralB)):
    line=dihedralB[i]
    f.write(dihedralline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9]))
    f.write('\n')
f.write('; Sidechain improper dihedrals\n')
for i in range(len(dihedralS)):
    line=dihedralS[i]
    f.write(dihedralSline%(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8]))
    f.write('\n')

### writing position restraint
f.write('\n')
f.write('#ifdef POSRES\n')
f.write('#ifndef POSRES_FC\n')
f.write('#define POSRES_FC 1000.00\n')
f.write('#endif\n')
f.write(' [ position_restraints ]\n')
for i in range(len(posit)):
    line=posit[i]
    f.write(positline%(line[0],line[1],line[2],line[3],line[4]))
    f.write('\n')
f.write('#endif')


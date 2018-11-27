import sys
import numpy as np
import MDAnalysis
import helanal

pdbaa = sys.argv[1]       ## All atom pdb file
itp = sys.argv[2]         ## standard MARTINI protein topology file
chain = sys.argv[3]       ## select a chain
hs = sys.argv[4]          ## transmembrane helix starting residue 
he = sys.argv[5]          ## transmembrane helix starting residue
    

###########################################
### obtain atom, residue, dssp from pdb ###
###########################################
pdbIndex = []
for rl in open(pdbaa):
    srl = rl.split()[:]
    if srl[0] == 'ATOM':
        if srl[2] == 'CA' and srl[4] == chain:
            pdbIndex.append(srl[1:])
pdbResin = [int(x[4]) for x in pdbIndex]                         

itpIndex = []
for rl in open(itp):
    srl = rl.split()[:]
    if len(srl) == 9 and srl[4] == 'BB':
            itpIndex.append(srl[:])
itpAtom = [int(x[0]) for x in itpIndex]                            ## itp BB atoms index
itpResit = [str(x[3]) for x in itpIndex]                           ## itp residue type index
itpResin = [int(x[2])+int(pdbResin[0])-1 for x in itpIndex]        ## pdb residue number 
Dsspt = [str(x[-1]) for x in itpIndex]                             ## secondary structure
resiDsspDict = dict(zip(itpResin, Dsspt))



###########################
### Secondary structure ###
###########################

##             F   E   H   1   2   3   T   S   C
## ssnum  = ( 13,  4,  2,  2,  2,  2,  6, 22,  0 )   
## .HTH.  using .272.   ##  one turn in disrupted helix
## .HTTH. using .2882.  ##  two turn in disrupted helix


def dsspRep(ss):
    dsspt=[ 'F', 'E', 'H', '1', '2', '3', 'T',  'S', 'C', 'T', 'T']   
    dsspn=[13.0, 4.0, 2.0, 2.0, 2.0, 2.0, 6.0, 22.0, 0.0, 7.0, 8.0]
    for index in range(len(dsspt)):
        if ss == dsspt[index]:
           ss = dsspn[index]
    return ss

pdbDssp = [float(dsspRep(x)) for x in Dsspt]

for i in range(1,len(pdbDssp)-1):
    if pdbDssp[i-1] == 2.0 and pdbDssp[i] == 6.0 and pdbDssp[i+1] == 2.0:
        pdbDssp[i] = 7.0
    elif pdbDssp[i-1] == 2.0 and pdbDssp[i] == 6.0 and pdbDssp[i+1] == 6.0 and pdbDssp[i+2] == 2.0:
        pdbDssp[i] = 8.0
        pdbDssp[i+1] = 8.0



#################################
### bending angle calculation ###
#################################

def bac(pdb,ss,resi):
    import MDAnalysis
    import helanal
    benddata=[]
    residSL=[]
    residEL=[]
    bendresi=[]
    for i in range(len(ss)-1):
        if ss[i] != 2.0 and ss[i+1] == 2.0:
            residSL.append(resi[i+1])
        elif ss[i] == 2.0 and ss[i+1] != 2.0:
            residEL.append(resi[i])
    for seg in range(len(residSL)):
        if residEL[seg]-residSL[seg] > 7:
            bendresi += [x for x in range(residSL[seg], residEL[seg]+1)]
            benddata += [0,0,0]
           # data = helanal.helanal_main(pdb, selection='segid %s and name CA and resnum %s-%s'%(chain, residSL[seg], residEL[seg]))
            data = helanal.helanal_main(pdb, selection='segid %s and (name CA or name N or name C or name O)'%chain, start=residSL[seg], end=residEL[seg]) 
            benddata += data
            benddata += [0,0,0]
        else:
           # print residSL[seg], residEL[seg]
            bendresi += [x for x in range(residSL[seg], residEL[seg]+1)]
            benddata += [0,]*(residEL[seg] - residSL[seg] + 1)
   # print bendresi, benddata
   # print len(bendresi), len(benddata)
    bendDict = dict(zip(bendresi, benddata))
    return bendDict
bendDict = bac(pdbaa,pdbDssp,itpResin)

###Get Bending angles
bendIndex =[]
for index in range(int(hs), int(he)+1):
    bendIndex.append(index)

fhinge = open('hinge.xvg','w')  ### generating hinge file
fhinge.write('Chain   Residues   Dssp   Bending Angles\n')
for index in bendIndex:
    if index in bendDict:
        fhinge.write('%3s%9d%10s%12.2f\n'%(chain,index,resiDsspDict[index],bendDict[index]))
    else:
        fhinge.write('%3s%9d%10s%12.2f\n'%(chain,index,resiDsspDict[index],0.00))



################################################
### bending angle To dihedral force constant ###
################################################

def df(benda,bendb):
    if benda > 0 and bendb > 0:
        dfa = (44.97/benda)**(1/0.342)     ### power function between bending angle and dihedral force constant
        dfb = (44.97/bendb)**(1/0.342)
        df = (dfa*dfb)**0.5
    elif benda > 0 and bendb == 0:
        dfa = (44.97/benda)**(1/0.342)
        df = dfa
    elif benda == 0 and bendb > 0:
        dfb = (44.97/bendb)**(1/0.342)
        df = dfb
    elif benda == 0 and bendb == 0:
        df = 497.31
    if df > 497.31:
        df = 497.31
    elif df < 5.0:
        df = 5.00
    return df



########################################
### write backbone bonded parameters ###
########################################

print '[ bonds ]'
print '; Backbone bonds'
###    1     3      1   0.35000  1250 ; THR(C)-ALA(C)
bondline='%5d %5d       1   %7.5f  %5d ; %3s(%1s)-%3s(%1s)'

for i in range(1,len(bendIndex)):
    for j in range(1,len(itpResin)):
        if bendIndex[i-1:i+1] == itpResin[j-1:j+1] and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 7.0:
           print bondline%(itpAtom[j-1],itpAtom[j],0.33000,12500,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j])
           print bondline%(itpAtom[j],itpAtom[j+1],0.31000,12500,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
        elif bendIndex[i-2:i] == itpResin[j-2:j] and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 8.0:
           print bondline%(itpAtom[j-1],itpAtom[j],0.33000,12500,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j])
           print bondline%(itpAtom[j],itpAtom[j+1],0.33000,12500,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print bondline%(itpAtom[j+1],itpAtom[j+2],0.31000,12500,itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
print ''
print '[ constraints]'
###   24    26      1   0.33000 ; GLU(T)-THR(1)
constraintline='%5d %5d       1   0.33000 ; %3s(%1s)-%3s(%1s)'
for i in range(1,len(bendIndex)):
    for j in range(1,len(itpResin)):
        if bendIndex[i-1:i+1] == itpResin[j-1:j+1] and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 7.0:
           print constraintline%(itpAtom[j-2],itpAtom[j-1],itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1])
           print constraintline%(itpAtom[j+1],itpAtom[j+2],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
        elif bendIndex[i-1:i+1] == itpResin[j-1:j+1] and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 8.0:
           print constraintline%(itpAtom[j-2],itpAtom[j-1],itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1])
           print constraintline%(itpAtom[j+2],itpAtom[j+3],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])

print ''
print '[ Angles ]'
print '; Backbone angles'
###   26    28    30      2     96   700 ; THR(1)-LEU(1)-TRP(1)
angleline='%5d %5d %5d      2%7.2f%7.2f ; %3s(%1s)-%3s(%1s)-%3s(%1s)'
for i in range(2,len(bendIndex)):
    for j in range(2,len(itpResin)):
        if bendIndex[i-2:i+1] == itpResin[j-2:j+1] and pdbDssp[j-2] == 2.0 and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 2.0:
           print angleline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],87.69,462.58,itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j])
        elif bendIndex[i-2:i+1] == itpResin[j-2:j+1] and pdbDssp[j-2] == 2.0 and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 7.0:
           print angleline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],85.56,101.37,itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[i],Dsspt[j])
           print angleline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],100.00,20.00,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print angleline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],96.53,288.58,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[i+2],Dsspt[i+2])
        elif bendIndex[i-2:i+1] == itpResin[j-2:j+1] and pdbDssp[j-2] == 2.0 and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 8.0 and pdbDssp[j+1] == 8.0:
           print angleline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],81.76, 77.42,itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j])
           print angleline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],79.93, 50.00,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print angleline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],100.00,20.00,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
           print angleline%(itpAtom[j+1],itpAtom[j+2],itpAtom[j+3],95.43,229.60,itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])

print ''
print '[ Dihedrals ]'
print '; Backbone dihedrals'
###   26    28    30    35      1   -120   400     1 ; THR(1)-LEU(1)-TRP(1)-LEU(1)
dihedralline='%5d %5d %5d %5d      1%8.2f%7.2f     %d ; %3s(%1s)-%3s(%1s)-%3s(%1s)-%3s(%1s)'
for i in range(3,len(bendIndex)):
    for j in range(3,len(itpResin)):
        if bendIndex[i-3:i+1] == itpResin[j-3:j+1] and pdbDssp[j-3] == 2.0 and pdbDssp[j-2] == 2.0 and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 2.0:
           print dihedralline%(itpAtom[j-3],itpAtom[j-2],itpAtom[j-1],itpAtom[j],-122.00,df(bendDict[bendIndex[i-2]],bendDict[bendIndex[i-1]]),1,itpResit[j-3],Dsspt[j-3],itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j])

        elif bendIndex[i-3:i+1] == itpResin[j-3:j+1] and pdbDssp[j-3] == 2 and pdbDssp[j-2] == 2 and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 7.0:
           print dihedralline%(itpAtom[j-3],itpAtom[j-2],itpAtom[j-1],itpAtom[j],-140.49,46.40,1,itpResit[j-3],Dsspt[j-3],itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j])

           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1], 172.89,20.03,1,itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1], -33.93, 7.40,2,itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1],-180.00, 7.88,4,itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1], 151.92,-3.52,8,itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])

           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2],  -2.08,-11.53,1,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2],   5.56, 18.26,2,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2], -81.13, -0.34,3,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2],  -7.72,  4.26,4,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2], -39.79,  4.06,5,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2],  26.43, -1.01,6,itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])

           print dihedralline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],itpAtom[j+3],-125.09, 50.00,1,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])

        elif bendIndex[i-3:i+1] == itpResin[j-3:j+1] and pdbDssp[j-3] == 2.0 and pdbDssp[j-2] == 2.0 and pdbDssp[j-1] == 2.0 and pdbDssp[j] == 8.0 and pdbDssp[j+1] == 8.0:
           print dihedralline%(itpAtom[j-3],itpAtom[j-2],itpAtom[j-1],itpAtom[j],-128.48,50.00,1,itpResit[j-3],Dsspt[i-3],itpResit[j-2],Dsspt[j-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j])

           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1], 115.72,-20.00,1,itpResit[j-2],Dsspt[i-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1], -54.24,19.99,2,itpResit[j-2],Dsspt[i-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1], -34.53, 5.03,3,itpResit[j-2],Dsspt[i-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1], -65.67,16.05,4,itpResit[j-2],Dsspt[i-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1],-117.98,-18.83,5,itpResit[j-2],Dsspt[i-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])
           print dihedralline%(itpAtom[j-2],itpAtom[j-1],itpAtom[j],itpAtom[j+1],-25.77,-9.65,6,itpResit[j-2],Dsspt[i-2],itpResit[j-1],Dsspt[j-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1])

           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2],-159.29,14.67,1,itpResit[j-1],Dsspt[i-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])
           print dihedralline%(itpAtom[j-1],itpAtom[j],itpAtom[j+1],itpAtom[j+2],-106.94,-5.77,3,itpResit[j-1],Dsspt[i-1],itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2])

           print dihedralline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],itpAtom[j+3], 31.70,-18.32,1,itpResit[j],Dsspt[i],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])
           print dihedralline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],itpAtom[j+3], 52.06, 14.41,2,itpResit[j],Dsspt[i],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])
           print dihedralline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],itpAtom[j+3],-13.82, -2.91,3,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])
           print dihedralline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],itpAtom[j+3], 96.38, 16.05,4,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])
           print dihedralline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],itpAtom[j+3], 20.33, 20.00,5,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])
           print dihedralline%(itpAtom[j],itpAtom[j+1],itpAtom[j+2],itpAtom[j+3],116.90, -8.76,6,itpResit[j],Dsspt[j],itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3])

           print dihedralline%(itpAtom[j+1],itpAtom[j+2],itpAtom[j+3],itpAtom[j+4],-110.73, 50.00,1,itpResit[j+1],Dsspt[j+1],itpResit[j+2],Dsspt[j+2],itpResit[j+3],Dsspt[j+3],itpResit[j+4],Dsspt[j+4])




import sys
import numpy as np

itp_s1 = sys.argv[1]       ## the topology file
itp_s2 = sys.argv[2]       ## the topology file


itpIndex_s1 = []
itpIndex_s2 = []
enbond_s1 = []
enbond_s2 = []

for rl in open(itp_s1):
    srl = rl.split()[:]
    if len(srl) == 9 and srl[4] == "BB":
        itpIndex_s1.append(srl[:])
    if len(srl) == 5 and srl[-1] == "RUBBER_FC*1.000000":
        enbond_s1.append(srl[:])

for rl in open(itp_s2):
    srl = rl.split()[:]
    if len(srl) == 9 and srl[4] == "BB":
        itpIndex_s2.append(srl[:])
    if len(srl) == 5 and srl[-1] == "RUBBER_FC*1.000000": 
        enbond_s2.append(srl[:])


coil=['F','E','T','S','C']


def enbondDssp(enbond):
    ss_s1 = ""
    ss_s2 = ""
    for i in range(len(itpIndex_s1)):
        if float(enbond[0]) <= float(itpIndex_s1[i][0]) <= float(enbond[1]):
           ss_s1 += itpIndex_s1[i][-1]
    for i in range(len(itpIndex_s2)):
        if float(enbond[0]) <= float(itpIndex_s2[i][0]) <= float(enbond[1]):
           ss_s2 += itpIndex_s2[i][-1]
    if (coil[0] in ss_s1 or coil[1] in ss_s1 or coil[2] in ss_s1 or coil[3] in ss_s1 or coil[4] in ss_s1) and (coil[0] in ss_s2 or coil[1] in ss_s2 or coil[2] in ss_s2 or coil[3] in ss_s2 or coil[4] in ss_s2):
        return enbond

enData = []
for i in range(len(enbond_s1)):
    for j in range(len(enbond_s2)):
        diff = min([float(enbond_s1[i][3]),float(enbond_s2[j][3])])*0.1
        if enbond_s1[i][0] == enbond_s2[j][0] and enbond_s1[i][1] == enbond_s2[j][1] and -diff <= float(enbond_s1[i][3])-float(enbond_s2[j][3]) <= diff:
            if enbondDssp(enbond_s1[i]) != None:
                enData.append(enbondDssp(enbond_s1[i]))

print '#ifndef NO_RUBBER_BANDS'
print '#ifndef RUBBER_FC'
print '#define RUBBER_FC 500.000000'
print '#endif'
for i in range(len(enData)):
    print '%5s%6s%7s%10s%19s'%(enData[i][0], enData[i][1],enData[i][2],enData[i][3],enData[i][4])
print '#endif'




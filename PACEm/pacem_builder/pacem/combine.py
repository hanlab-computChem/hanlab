import sys
import math


pace_pdb = open(sys.argv[1],'r')    ### pace.pdb 
cg_pdb   = open(sys.argv[2],'r')    ### cg.pdb
index    = open(sys.argv[3],'r')    ### pace_residue.txt
pace_itp = sys.argv[4]              ### pace.itp
cg_itp   = sys.argv[5]              ### cg.itp


def readPDB(pdb):
    dataPDB=[]
    for rl in pdb:
        if rl[:4] == "ATOM":
           dataPDB.append(rl[:-1])
    return dataPDB

def get_resPDB(resnum, dataPDB, resCN):
    resPDB=[] 
    for rl in dataPDB:
        res = int(rl[22:27])
        res_atom_type = str(rl[12:15])
        if resnum == res and resCN == "res_noCN":           
           resPDB.append(rl)
        if resnum == res and resCN == "res_NTF":
           if res_atom_type != " N " and res_atom_type != " H ":
              resPDB.append(rl)
        if resnum == res and resCN == "res_CTF":
           if res_atom_type != " C " and res_atom_type != " O ":
              resPDB.append(rl)
    return resPDB

#  pace replace
resPACE  = []
resNTF   = []
resCTF   = []
num = 0
for rl in index:
    srl = rl.split()[:]
    r1 = int(srl[0])
    r2 = int(srl[1])
    if num == 0:
       resStart = int(r1)
       resEnd   = int(r2)
    else:
       for i in range(r1, r2 + 1):
           resPACE.append(i)
       resNTF.append(r1)
       resCTF.append(r2)
    num += 1 
#print(resPACE)

resPACE_TER_noStart_noEnd = []
for res_id in resNTF:
    if res_id != resStart and res_id != resEnd:
       resPACE_TER_noStart_noEnd.append(res_id)
for res_id in resCTF:
    if res_id != resStart and res_id != resEnd:
       resPACE_TER_noStart_noEnd.append(res_id)
#print(resPACE_TER_noStart_noEnd)

pacePDB = readPDB(pace_pdb)
cgPDB   = readPDB(cg_pdb)
fwpdb   = open("pacem.pdb",'w')
new_pdb_backbone = []
new_pdb          = []
atomNum      = 1
atomComb     = []
atomCombCG   = []
atomCombPACE = []
for resnum in range(resStart, resEnd+1):
    aaType = "CG" 
    if resnum in resPACE:
       aaType = "PACE"
       if resnum == resStart or resnum == resEnd:
          resPDB = get_resPDB(resnum, pacePDB, "res_noCN")
       else:
          if resnum not in resCTF and resnum not in resNTF:
             resPDB = get_resPDB(resnum, pacePDB, "res_noCN") 
          elif resnum in resCTF:
             resPDB = get_resPDB(resnum, pacePDB, "res_CTF")
          elif resnum in resNTF:
             resPDB = get_resPDB(resnum, pacePDB, "res_NTF")
    else:
       resPDB = get_resPDB(resnum, cgPDB, "res_noCN")
    for line in resPDB:
       fwpdb.write("%s %6d %s\n"%(line[:4],atomNum,line[12:]))
       if line[12:17] == " CA  " or line[12:17] == "  BB ":
          new_pdb_backbone.append("%s %6d %s\n"%(line[:4],atomNum,line[12:]))
       new_pdb.append("%s %6d %s\n"%(line[:4],atomNum,line[12:]))
       atomNum_o = int(line[4:12])
       atomComb.append([aaType, atomNum, atomNum_o])
       if aaType == "CG":
          atomCombCG.append([aaType, atomNum, atomNum_o])
       elif aaType == "PACE":
          atomCombPACE.append([aaType, atomNum, atomNum_o])
       atomNum += 1

pace_atoms_o = [atomCombPACE[i][2] for i in range(len(atomCombPACE))]
pace_atoms_n = [atomCombPACE[i][1] for i in range(len(atomCombPACE))]
cg_atoms_o   = [atomCombCG[i][2] for i in range(len(atomCombCG))]
cg_atoms_n   = [atomCombCG[i][1] for i in range(len(atomCombCG))]


### Read ITP files
pace_itp_type = ['[ atoms ]','[ bonds ]','[ exclusions ]','[ pairs ]','[ angles ]','[ dihedrals ]','[ end ]']
cg_itp_type   = ['[ atoms ]','[ bonds ]','[ constraints ]','[ angles ]','[ dihedrals ]', '#ifdef POSRES']
def extract_content(filename, itp_type):
    itp_content = []    
    for i in range(len(itp_type)-1):
        start_line = itp_type[i]
        end_line   = itp_type[i+1]
        content    = []
        #content.append(";; " + start_line)
        is_extracting = False
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line == start_line:
                   is_extracting = True
                elif line == end_line:
                   is_extracting = False
                elif is_extracting:
                   content.append(line)
        itp_content.append(content)
    return itp_content
pace_itp_content = extract_content(pace_itp, pace_itp_type)
cg_itp_content   = extract_content(cg_itp, cg_itp_type)

atom_line      = "%7s %5s %5s %5s %5s %7s %9s  %s  %s"
bond_line      = "%7s %7s   %s "
angle_line     = "%7s %7s %7s   %s "
dihedral_line  = "%7s %7s %7s %7s   %s "
def get_new_itp(itp_content, itp_type, model):
    if model == "PACE":
       atoms_o = pace_atoms_o
       atoms_n = pace_atoms_n
    elif model == "CG":
       atoms_o = cg_atoms_o
       atoms_n = cg_atoms_n
    new_itp_content = []
    new_itp_content.append(itp_type)
    new_itp_content.append("; %s %s"%(model,itp_type.split()[1]))
    if itp_type == "[ atoms ]":
       for line in itp_content:
           srl = line.split()[:]
           if len(line) != 0:
              if str(line[0]) != ";":
                 a1 = int(srl[0])
                 for i in range(len(atoms_o)):
                     if a1 == atoms_o[i] :
                        new_itp_content.append(atom_line%(atoms_n[i], srl[1], srl[2], srl[3], srl[4], atoms_n[i], srl[6], srl[7], ' '.join(srl[8:])))
              else:
                 new_itp_content.append(line)

    elif itp_type == "[ bonds ]" or itp_type == "[ constraints ]" or itp_type == "[ exclusions ]" or itp_type == "[ pairs ]":
       for line in itp_content:
           srl = line.split()[:]
           if len(line) != 0:
              if str(line[0]) not in [";", "#"]:
                 a1 = int(srl[0])
                 a2 = int(srl[1])
                 for i in range(len(atoms_o)):
                     for j in range(len(atoms_o)):
                         if a1 == atoms_o[i] and a2 == atoms_o[j]:
                            new_itp_content.append(bond_line%(atoms_n[i], atoms_n[j], ' '.join(srl[2:])))
              else:
                 new_itp_content.append(line)

    elif itp_type == "[ angles ]":
       for line in itp_content:
           srl = line.split()[:]
           if len(line) != 0:
              if str(line[0]) not in [";", "#"]:
                 a1 = int(srl[0])
                 a2 = int(srl[1])
                 a3 = int(srl[2])
                 for i in range(len(atoms_o)):
                     if a1 == atoms_o[i]:
                        for j in range(len(atoms_o)):
                            if a2 == atoms_o[j]:
                               for k in range(len(atoms_o)):
                                   #if a1 == atoms_o[i] and a2 == atoms_o[j] and a3 == atoms_o[k]:
                                   if a3 == atoms_o[k]:
                                      new_itp_content.append(angle_line%(atoms_n[i], atoms_n[j], atoms_n[k], ' '.join(srl[3:])))
              else:
                 new_itp_content.append(line)
    elif itp_type == "[ dihedrals ]":
       for line in itp_content:
           srl = line.split()[:]
           if len(line) != 0:
              if str(line[0]) not in [";", "#"]:
                 a1 = int(srl[0])
                 a2 = int(srl[1])
                 a3 = int(srl[2])
                 a4 = int(srl[3])
                 for i in range(len(atoms_o)):                   
                     if a1 == atoms_o[i]: 
                        for j in range(len(atoms_o)):                   
                            if a2 == atoms_o[j]: 
                               for k in range(len(atoms_o)):                   
                                   if a3 == atoms_o[k]: 
                                      for l in range(len(atoms_o)):                   
                                          if a4 == atoms_o[l]: 
                                             new_itp_content.append(dihedral_line%(atoms_n[i], atoms_n[j], atoms_n[k],atoms_n[l], ' '.join(srl[4:])))
              else:
                 new_itp_content.append(line)       
    return new_itp_content



### writing pacem.itp 
###
print("[ moleculetype ]")
print("; Name         Exclusions")
print("PROA_P            3\n")

print("[ atoms ]")
print(";   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB")
atoms_PACE = get_new_itp(pace_itp_content[0],"[ atoms ]","PACE")
atoms_CG   = get_new_itp(cg_itp_content[0],"[ atoms ]","CG")
for rl in new_pdb:
    atom_ID = int(rl[5:11])
    for rl_cg in atoms_CG:
      if rl_cg[0] not in ["[",";","#"]:
        atom_ID_CG = int(rl_cg.split()[0])
        if atom_ID == atom_ID_CG:
           print("%s"%rl_cg)
    for rl_pace in atoms_PACE:
      if rl_pace[0] not in ["[",";","#"]:
        atom_ID_pace = int(rl_pace.split()[0])
        if atom_ID == atom_ID_pace:
           print("%s"%rl_pace)
print

### writing [ bonds ] [ angles ] [ dihedrals ] and so on
col = 1
for itp_type in pace_itp_type[1:-1]:
    new_line = get_new_itp(pace_itp_content[col],itp_type,"PACE")
    for line in new_line:
        print line
    col += 1
    print
col = 1
for itp_type in cg_itp_type[1:-1]:
    new_line = get_new_itp(cg_itp_content[col],itp_type,"CG")
    for line in new_line:
        print line
    col += 1
    print


### write CG exclusion 1-2 1-3 LJ interaction
bonds_CG         = get_new_itp(cg_itp_content[1],"[ bonds ]","CG")
constraints_CG   = get_new_itp(cg_itp_content[2],"[ constraints ]","CG")
bonds_data       = bonds_CG + constraints_CG
bonds_data_noEN  = []

atom_pair_ID      = []
atom_pair_ID_type = []
for rl in atoms_CG:
    if rl[0] not in [";","#","["]:
       srl           = rl.split()[:]
       atom_pair_ID.append(int(srl[0]))
       atom_pair_ID_type.append(str(srl[1]))
atom_ID_to_type = dict(zip(atom_pair_ID, atom_pair_ID_type))

atom_pair_pace_ID      = []
atom_pair_pace_ID_type = []
for rl in atoms_PACE:
    if rl[0] not in [";","#","["]:
       srl           = rl.split()[:]
       atom_pair_pace_ID.append(int(srl[0]))
       atom_pair_pace_ID_type.append(str(srl[1]))
atom_ID_to_type_pace = dict(zip(atom_pair_pace_ID, atom_pair_pace_ID_type))
 
### martini bonds 
for rl in bonds_data:
    if rl[0] not in [";","#","["]:
       srl = rl.split()[:]
       if srl[-1] != "RUBBER_FC*1.000000":
          bonds_data_noEN.append([int(srl[0]),int(srl[1])])
          bonds_data_noEN.append([int(srl[1]),int(srl[0])])

### PACE bonds
bonds_PACE        = get_new_itp(pace_itp_content[1],"[ bonds ]","PACE")
bonds_data_PACE = []
for rl in bonds_PACE:
    if rl[0] not in [";","#","["]:
       srl = rl.split()[:]
       bonds_data_PACE.append([int(srl[0]),int(srl[1])])       

#### lj potential
pairs_lj_cg = []
for rl in open('ffPACE_1.3vdW.itp','r'):
    srl=rl.split()[:]
    if len(srl) >= 7 and srl[5] == ';':
       pairs_lj_cg.append(srl)
def generate_pairs_lj(at1, at2):
    for pair in pairs_lj_cg:
       if (at1 == pair[0] and at2 == pair[1]) or (at2 == pair[0] and at1 == pair[1]):
          content ='1  %s  %s  ; %s  %s  %s '%(pair[3], pair[4],at1, at2, ' '.join(pair[6:]))
    return content


pairs_CG = []
print "[ pairs ]"
print "; CG 1-2 1-3"
for atom_pair_ID0 in atom_pair_ID:
    atom_pair_ID0_type = atom_ID_to_type[int(atom_pair_ID0)]
    for i in range(len(bonds_data_noEN)):
           atom_pair_ID1 = bonds_data_noEN[i][0]    
           atom_pair_ID2 = bonds_data_noEN[i][1]    
           if atom_pair_ID0 == atom_pair_ID1:
              for j in range(len(bonds_data_noEN)):
                  atom_pair_ID3 = bonds_data_noEN[j][0]    
                  atom_pair_ID4 = bonds_data_noEN[j][1]    
                  atom_pair_ID4_type = atom_ID_to_type[int(atom_pair_ID4)]
                  atom_pair1 = [atom_pair_ID0, atom_pair_ID4]
                  atom_pair2 = [atom_pair_ID4, atom_pair_ID0]
                  if atom_pair_ID2 == atom_pair_ID3 and atom_pair_ID0 != atom_pair_ID4 and atom_pair1 not in pairs_CG and atom_pair2 not in pairs_CG and atom_pair1 not in bonds_data_noEN and atom_pair2 not in bonds_data_noEN:
                     print("%7s %7s %s"%(atom_pair_ID0, atom_pair_ID4, generate_pairs_lj(atom_pair_ID0_type, atom_pair_ID4_type)))
                     pairs_CG.append(atom_pair1)

                     for k in range(len(bonds_data_noEN)):
                         atom_pair_ID5 = bonds_data_noEN[k][0]    
                         atom_pair_ID6 = bonds_data_noEN[k][1]    
                         atom_pair_ID6_type = atom_ID_to_type[int(atom_pair_ID6)]
                         atom_pair3 = [atom_pair_ID0, atom_pair_ID6]
                         atom_pair4 = [atom_pair_ID6, atom_pair_ID0]
                         if atom_pair_ID4 == atom_pair_ID5 and atom_pair_ID0 != atom_pair_ID6 and atom_pair3 not in pairs_CG and atom_pair4 not in pairs_CG and atom_pair3 not in bonds_data_noEN and atom_pair4 not in bonds_data_noEN:
                            print("%7s %7s %s"%(atom_pair_ID0, atom_pair_ID6, generate_pairs_lj(atom_pair_ID0_type, atom_pair_ID6_type)))
                            pairs_CG.append(atom_pair3)


### print joint residue bonded parameters


def find_pair(index, bonds_type, exclude_index):
    pair_atoms = []
    for i in range(len(bonds_type)):
         bonds_line = bonds_type[i]
         atom_p1 = int(bonds_line[0])
         atom_p2 = int(bonds_line[1])
         if atom_p1 in index and atom_p2 not in index and atom_p2 not in pair_atoms and atom_p2 not in exclude_index:
            pair_atoms.append(atom_p2)
         if atom_p2 in index and atom_p1 not in index and atom_p1 not in pair_atoms and atom_p1 not in exclude_index:
            pair_atoms.append(atom_p1)
    return pair_atoms           

print 
print "; interface residue bonded parameters"
def bond_interaction(ln1, ln2):
    x1 = float(ln1[30:38])
    y1 = float(ln1[38:46])
    z1 = float(ln1[46:54])
    x2 = float(ln2[30:38])
    y2 = float(ln2[38:46])
    z2 = float(ln2[46:54])
    bond = 0.1*((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
    return bond
def angle_interaction(ln1, ln2, ln3):
    x1 = float(ln1[30:38])
    y1 = float(ln1[38:46])
    z1 = float(ln1[46:54])
    x2 = float(ln2[30:38])
    y2 = float(ln2[38:46])
    z2 = float(ln2[46:54])
    x3 = float(ln3[30:38])
    y3 = float(ln3[38:46])
    z3 = float(ln3[46:54])
    AB = [x2 - x1, y2 - y1, z2 - z1]
    BC = [x3 - x2, y3 - y2, z3 - z2]
    dot_product = AB[0] * BC[0] + AB[1] * BC[1] + AB[2] * BC[2]
    magnitude_AB = math.sqrt(AB[0]**2 + AB[1]**2 + AB[2]**2)
    magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2 + BC[2]**2)
    angle_rad = math.pi - math.acos(dot_product / (magnitude_AB * magnitude_BC))
    angle_deg = math.degrees(angle_rad)
    return angle_deg


for res_id in resNTF:
  if res_id != resStart and res_id != resEnd:
    for line in new_pdb_backbone:
        res_id_pdb  = int(line[22:27])
        atom_id_pdb = int(line[5:11])
        if res_id_pdb == res_id - 2:
           ln1 = line
           atomID1 = atom_id_pdb
        if res_id_pdb == res_id - 1:
           ln2 = line
           atomID2 = atom_id_pdb
        if res_id_pdb == res_id:
           ln3 = line
           atomID3 = atom_id_pdb
        if res_id_pdb == res_id + 1:
           ln4 = line
           atomID4 = atom_id_pdb
        #if res_id_pdb == res_id + 2:
        #   ln0 = line
        #   atomID0 = atom_id_pdb
    print("[ bonds ]")
    print("%7s %7s  1  %7.5f  150000 ; %s"%(atomID2, atomID3, bond_interaction(ln2, ln3), res_id))
    print("#ifndef NO_RUBBER_BANDS")
    print("#ifndef RUBBER_FC")
    print("#define RUBBER_FC 500.000000")
    print("#endif ")
    for line in new_pdb_backbone:
        res_id_pdb  = int(line[22:27])
        atom_id_pdb = int(line[5:11])
        res_atom_type = str(line[12:17]) 
        if res_atom_type == " CA  " and abs(res_id_pdb - (res_id - 2)) > 2:
           ln5 = line
           en1 = bond_interaction(ln1, ln5)
           if en1 <= 0.9 and en1 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID1, en1, res_id - 2))
        if res_atom_type == " CA  " and abs(res_id_pdb - (res_id - 1)) > 2:
           ln5 = line
           en2 = bond_interaction(ln2, ln5)
           if en2 <= 0.9 and en2 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID2, en2, res_id - 1))        
        if res_atom_type == "  BB " and abs(res_id_pdb - res_id) > 2:
           ln5 = line
           en3 = bond_interaction(ln3, ln5)
           if en3 <= 0.9 and en3 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID3, en3, res_id))
        if res_atom_type == "  BB " and abs(res_id_pdb - (res_id + 1)) > 2:
           ln5 = line
           en4 = bond_interaction(ln4, ln5)
           if en4 <= 0.9 and en4 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID4, en4, res_id + 1))
       # if res_atom_type == "  BB " and abs(res_id_pdb - (res_id + 2)) > 3:
       #    ln5 = line
       #    en0 = bond_interaction(ln0, ln5)
       #    if en0 <= 0.9 and en0 >= 0.5:
       #       print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID0, en0, res_id + 2))
    print("#endif ")
 
    ### interface pairs 
    print("[ pairs ]")
    l3 = [atomID2]  ## CG
    l4 = [atomID3]  ## PACE
    l5 = find_pair(l4, bonds_data_PACE,[])
    l6 = find_pair(l5, bonds_data_PACE,l4)
    l2 = find_pair(l3, bonds_data_noEN,[])
    l1 = find_pair(l2, bonds_data_noEN,l3)
    at3 = atom_ID_to_type[int(atomID2)]
    at4 = atom_ID_to_type_pace[int(atomID3)]
    if len(l1) != 0:
       for a1 in l1:
           at1 = atom_ID_to_type[int(a1)]
           print("%7s %7s %s"%(a1, atomID3, generate_pairs_lj(at1, at4)))
    if len(l2) != 0:
       for a2 in l2:
           at2 = atom_ID_to_type[int(a2)]
           print("%7s %7s %s"%(a2, atomID3, generate_pairs_lj(at2, at4)))
           if len(l5) != 0:
              for a5 in l5:
                  at5 = atom_ID_to_type_pace[int(a5)]
                  print("%7s %7s %s"%(a2, a5, generate_pairs_lj(at2, at5)))
    if len(l5) != 0:
       for a5 in l5:
           at5 = atom_ID_to_type_pace[int(a5)]
           print("%7s %7s %s"%(atomID2, a5, generate_pairs_lj(at3, at5)))
    if len(l6) != 0:
       for a6 in l6:
           at6 = atom_ID_to_type_pace[int(a6)]
           print("%7s %7s %s"%(atomID2, a6, generate_pairs_lj(at3, at6)))

    ## backbone angles
    print("[ angles ]")
    print("%7s %7s %7s  2  %10.5f   40 ; %s"%(atomID1, atomID2, atomID3, angle_interaction(ln1, ln2, ln3), res_id))
    print("%7s %7s %7s  2  %10.5f   40 ; %s"%(atomID2, atomID3, atomID4, angle_interaction(ln2, ln3, ln4), res_id))

print
for res_id in resCTF:
  if res_id != resStart and res_id != resEnd:
    for line in new_pdb_backbone:
        res_id_pdb  = int(line[22:27])
        atom_id_pdb = int(line[5:11])
        #if res_id_pdb == res_id - 2:
        #   ln0 = line
        #   atomID0 = atom_id_pdb
        if res_id_pdb == res_id - 1:
           ln1 = line
           atomID1 = atom_id_pdb
        if res_id_pdb == res_id:
           ln2 = line
           atomID2 = atom_id_pdb
        if res_id_pdb == res_id + 1:
           ln3 = line
           atomID3 = atom_id_pdb
        if res_id_pdb == res_id + 2:
           ln4 = line
           atomID4 = atom_id_pdb
    print("[ bonds ]")
    print("%7s %7s  1  %7.5f  150000 ; %s"%(atomID2, atomID3, bond_interaction(ln2, ln3), res_id))
    print("#ifndef NO_RUBBER_BANDS")
    print("#ifndef RUBBER_FC")
    print("#define RUBBER_FC 500.000000")
    print("#endif ")
    for line in new_pdb_backbone:
        res_id_pdb  = int(line[22:27])
        atom_id_pdb = int(line[5:11])
        res_atom_type = str(line[12:17]) 
        #if res_atom_type == "  BB " and abs(res_id_pdb - (res_id-2)) > 3:
        #   ln5 = line
        #   en0 = bond_interaction(ln0, ln5)
        #   if en0 <= 0.9 and en0 >= 0.5:
        #      print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID0, en0, res_id-2))        
        if res_atom_type == "  BB " and abs(res_id_pdb - (res_id-1)) > 2:
           ln5 = line
           en1 = bond_interaction(ln1, ln5)
           if en1 <= 0.9 and en1 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID1, en1, res_id - 1))
        if res_atom_type == "  BB " and abs(res_id_pdb - res_id) > 2:
           ln5 = line
           en2 = bond_interaction(ln2, ln5)
           if en2 <= 0.9 and en2 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID2, en2, res_id))
        if res_atom_type == " CA  " and abs(res_id_pdb - (res_id-1)) > 2:
           ln5 = line
           en3 = bond_interaction(ln3, ln5)
           if en3 <= 0.9 and en3 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID3, en3, res_id + 1))
        if res_atom_type == " CA  " and abs(res_id_pdb - res_id) > 2:
           ln5 = line
           en4 = bond_interaction(ln4, ln5)
           if en4 <= 0.9 and en4 >= 0.5:
              print("%7s %7s   1  %7.5f  RUBBER_FC*1.000000 ; %s"%(atom_id_pdb, atomID4, en4, res_id + 2))
    print("#endif ")

    ### interface pairs 
    print("[ pairs ]")
    l3 = [atomID2]  ## PACE
    l4 = [atomID3]  ## CG
    l5 = find_pair(l4, bonds_data_noEN,[])
    l6 = find_pair(l5, bonds_data_noEN,l4)
    l2 = find_pair(l3, bonds_data_PACE,[])
    l1 = find_pair(l2, bonds_data_PACE,l3)
    at3 = atom_ID_to_type_pace[int(atomID2)]
    at4 = atom_ID_to_type[int(atomID3)]
    if len(l1) != 0:
       for a1 in l1:
           at1 = atom_ID_to_type_pace[int(a1)]
           print("%7s %7s %s"%(a1, atomID3, generate_pairs_lj(at1, at4)))
    if len(l2) != 0:
       for a2 in l2:
           at2 = atom_ID_to_type_pace[int(a2)]
           print("%7s %7s %s"%(a2, atomID3, generate_pairs_lj(at2, at4)))
           if len(l5) != 0:
              for a5 in l5:
                  at5 = atom_ID_to_type[int(a5)]
                  print("%7s %7s %s"%(a2, a5, generate_pairs_lj(at2, at5)))
    if len(l5) != 0:
       for a5 in l5:
           at5 = atom_ID_to_type[int(a5)]
           print("%7s %7s %s"%(atomID2, a5, generate_pairs_lj(at3, at5)))
    if len(l6) != 0:
       for a6 in l6:
           at6 = atom_ID_to_type[int(a6)]
           print("%7s %7s %s"%(atomID2, a6, generate_pairs_lj(at3, at6)))
    ### backbone angle
    print("[ angles ]")
    print("%7s %7s %7s  2  %10.5f   40 ; %s"%(atomID1, atomID2, atomID3, angle_interaction(ln1, ln2, ln3), res_id))
    print("%7s %7s %7s  2  %10.5f   40 ; %s"%(atomID2, atomID3, atomID4, angle_interaction(ln2, ln3, ln4), res_id))  

print
### postion restraint
print("#ifdef POSRES")
print("#ifndef POSRES_FC")
print("#define POSRES_FC 1000.00")
print("#endif")
print(" [ position_restraints ]")

posre_line = "%7s   1   POSRES_FC    POSRES_FC    POSRES_FC"
for line in new_pdb_backbone:
    atom_ID_backbone = int(line[5:11])
    print(posre_line%(atom_ID_backbone))
print("#endif")



import sys
import numpy as np

for rl in open(sys.argv[1],'r'):
  ss=rl
helix=["H","G"]
turn=["T"]
_h_res=[]
for i in range(len(ss)-4):
    s0=ss[i]
    s1=ss[i+1]
    s2=ss[i+2]
    s3=ss[i+3]
    if s0 in helix and s1 in helix and s2 in helix and s3 in helix:
       _h_res.append(i+1)  

if len(_h_res) != 0:
    print "; helical residue ID: ", " ".join([str(item) for item in _h_res])
else:
    print "; helical residue ID: NAN"
    

_t1_res=[]
for i in range(len(ss)-5):
    s0=ss[i]
    s1=ss[i+1]
    s2=ss[i+2]
    s3=ss[i+3]
    s4=ss[i+4]
    if s0 in helix and s1 in helix and s2 in helix and s3 in turn and s4 in helix:
       _t1_res.append(i+1)
       _t1_res.append(i+2)
       _t1_res.append(i+3)
       _t1_res.append(i+4)


_t2_res=[]
for i in range(len(ss)-6):
    s0=ss[i]
    s1=ss[i+1]
    s2=ss[i+2]
    s3=ss[i+3]
    s4=ss[i+4]
    s5=ss[i+5]
    if s0 in helix and s1 in helix and s2 in helix and s3 in turn and s4 in turn and s5 in helix:
       _t2_res.append(i+1)
       _t2_res.append(i+2)
       _t2_res.append(i+3)
       _t2_res.append(i+4)
       _t2_res.append(i+5)


_atom=[]
_resid=[]
for rl in open(sys.argv[2],'r'):
   srl=rl.split()[:]
   col=len(srl)
   #print col
   if col==11: 
      if srl[9]=="qtot":
         _atom.append(int(srl[0]))
         _resid.append(int(srl[2]))   
_a2r_dict=dict(zip(_atom,_resid))


HB_index = sys.argv[3]

if HB_index == "1.25":
    _h_pair="%s    %s  1  1.404606E-02    2.684245E-06 ;  %s\n"    ## 1.25
elif HB_index == "1.5":
    _h_pair="%s    %s  1  1.685528E-02	  3.221094E-06 ;  %s\n"   ## 1.5
elif HB_index == "2" or HB_index == "2.0":
    _h_pair="%s    %s  1  2.247370E-02    4.294792E-06 ;  %s\n"   ## 2.0 
else:
    _h_pair="%s    %s  1  1.123685E-02    2.147396E-06 ;  %s\n"   ## 1.0 


for rl in open(sys.argv[2],'r'):
   srl=rl.split()[:]
   col=len(srl)
   if col==7:
      if srl[6]=="Oi-Ni+3" or srl[6]=="Oi-Ni+4":
         _atom1=int(srl[0])
         _atom2=int(srl[1])
         _res1=_a2r_dict[_atom1]
         _res2=_a2r_dict[_atom2]
         #if _res1 in _h_res and _res2 in _h_res:
         if _res1 in _h_res or _res1 in _t1_res or _res1 in _t2_res:
             rl=_h_pair%(_atom1,_atom2,srl[6])
   print rl,




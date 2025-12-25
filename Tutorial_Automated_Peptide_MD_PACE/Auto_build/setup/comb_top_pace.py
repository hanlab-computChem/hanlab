import sys

ff_nm='./pace-asm.ff/forcefield.itp'
ff_other ='''
'''
ff_water_nm='./pace-asm.ff/cgWater.itp'
ff_ion_nm='./pace-asm.ff/ions.itp'

# main #

n_comp = int((len(sys.argv)-1)/2)

print( '; # components ',n_comp)

nm=[]
cnt=[]

for _i in range(n_comp):

 _nm_c = sys.argv[_i*2+1]
 _cnt_c = int(sys.argv[_i*2+2])

 print( '; comp name',_nm_c,'mol count', _cnt_c)

 nm.append(_nm_c)
 cnt.append(_cnt_c)



print( f'#include "{ff_nm}"')
print(ff_other)
print() 

mol_nm=[]
_suf=0

for _i in nm:
  print( '#include "%s.itp"'%_i)
 
  f_wrt = open(_i+'.itp','w')

  _st = False
  _end = False
  _find_mol = True
  
  for rl in open(_i+'.top','r'):
   if 'moleculetype' in rl: _st =True
   if _st:
    
    if _find_mol and rl[0]!=';' and rl[0]!='[':
     srl=rl[:-1].split()
     if len(srl)>=2:
#      if srl[0] in mol_nm:
#       mol_nm.append(srl[0]+'_'+chr(ord('A')+_suf))
#       _suf+=1
#      else:
#       mol_nm.append(srl[0])
      mol_nm.append(_i)
      f_wrt.write('%s    %s\n'%( mol_nm[-1],srl[1]))
      _find_mol = False
      continue

    f_wrt.write(rl)
    if _end and 'endif' in rl: break
    if 'POSRES' in rl: _end = True

  f_wrt.close()





print()
print( f'#include "{ff_water_nm}"')
print( f'#include "{ff_ion_nm}"')


print()
print( '[ system ]')
print( 'Co-assembly system')
print()
print( '[ molecules ]')
print( ';	Compound		#mols')

for _i in range(n_comp):
 print( '	%s			%d'%(mol_nm[_i], cnt[_i]))
print()


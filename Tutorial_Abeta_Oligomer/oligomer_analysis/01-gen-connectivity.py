import sys
sys.path.append('../basic_code/')
import numpy
import MDAnalysis
import geom_opt as geom
import networkutil
import time

## Usage ##
'''
python2 01-gen-connectivity.py
'''

#MAIN
top_nm = 'data/ab40.pdb'
trj_nm = 'data/ab40-demo.xtc' 

_num_ctt = 10 
fo_ctt = open('data/ab40-demo-cont%d.xvg'%(_num_ctt),'w')

system = MDAnalysis.Universe(top_nm,trj_nm)
_all = system.select_atoms('all')

n_mono = 100
# contact parameters
clash_cut = 4.5 # unit A
cell_size = 9.0 # cell size
clash_cut2 = clash_cut*clash_cut

nbbatom_per_mono = len(_all)/n_mono
fo_ctt.write('Number of backbone atoms in a monomer ' + str(nbbatom_per_mono)+'\n')
print nbbatom_per_mono


i = 0
for fr in system.trajectory:
 box = system._trajectory.ts._unitcell[:3]
 i+= 1
 if i ==1:
  cl = geom.link_cell(box, clash_cut+1.0)
 else:
  cl.clear_atoms()
  cl.update_box(box)

 _num = numpy.zeros((n_mono,n_mono),dtype=numpy.int)
 _all_cor = _all.positions

 try:
   cl.load_atoms(_all_cor)
 except IndexError:
   continue

 for _i,_j in enumerate(_all_cor):
  _m_id = _i/nbbatom_per_mono
  for _k in cl.contact(_j, _all_cor, clash_cut2, PBC_corr=False):
   _num[_m_id,_k/nbbatom_per_mono] += 1

 ############## For Contact 10 ###################
 _contact_list = [[] for _i in range(n_mono)]
 _idx = numpy.where(_num >= _num_ctt)   

 for _ee,_rr in enumerate(_idx[0]):
  _contact_list[_rr] += [_idx[1][_ee]]

 fo_ctt.write('===time ' + str(fr.time) +' fr id '+str(fr.frame)+'\n')
 _res_contact=[]
 for _k,_i in enumerate(_contact_list):
  _res_contact.append(networkutil.sort_unq(_i))
  fo_ctt.write(str(_k)+':')
  for _kk in  _res_contact[-1]:
   fo_ctt.write(str(_kk)+' ')
  fo_ctt.write('\n')

fo_ctt.close()

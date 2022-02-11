import sys
sys.path.append('/share1/home/caoy/pyscript/')
import numpy as np
import os
import analysis_topology as anaTop 
import cal_cont_cluster as cal_cont


## Usage ##
'''
python2 04-l-cycle.py  
'''
_frame = 450 ###ps
_nchain = 100
_start = 5050000
_end = 5100000
_time_fr10 , _time_p=_start,_start


_f10 = 'data/ab40-demo-cont10.xvg'
_f10_p = 'data/ab40-demo-cont10-olig.xvg'

fconn = open('data/ab40-demo-clu_mem-clu_conn.xvg','w')

_fr_id10,_fr_id = 0,0
while _time_fr10+_frame <_end:
  ### smaller than the end time of this file
  _cal_ctt = cal_cont.cal_cont('','',100)
  _time_fr10, _conn,_fr_id10 = _cal_ctt.load_cluster_noself(_f10, _fr_id10)
  _time_p, clu_size, clu_mem, _fr_id = _cal_ctt.load_patch(_f10_p, _fr_id)
  fconn.write('=== time %d ===\n'%(_time_fr10))
  
  if _time_fr10 != _time_p:
    print '%d did not match with %d'%(_time_fr10,_time_p)
    sys.exit()

  for _index_clu,_size in enumerate(clu_size):
    _mem = clu_mem[_index_clu]

    if len(_mem) >2:
      fconn.write('cluster member: ')

      _mem_conn = []
      for _member in _mem:
        fconn.write(' %d '%_member)
        _mem_conn.append(_conn[_member])
      fconn.write('\n')
      fconn.write(str(_mem_conn)+'\n')

fconn.close()

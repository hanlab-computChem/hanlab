import sys
sys.path.append('../basic_code/')
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

_fcycle = 'data/ab40-demo-cont10-olig-cycle.xvg'
_fchain = 'data/ab40-demo-cont10-olig-chain.xvg'

_fr_id10,_fr_id = 0,0
with open(_fchain,'w') as fchain:
  with open(_fcycle,'w') as fcycle:
    while _time_fr10+_frame <_end:

      ### smaller than the end time of this file
      _cal_ctt = cal_cont.cal_cont('','',100)
      _time_fr10, _conn,_fr_id10 = _cal_ctt.load_cluster_noself(_f10, _fr_id10)
      _time_p, clu_size, clu_mem, _fr_id = _cal_ctt.load_patch(_f10_p, _fr_id)
      fcycle.write('=== time %d\n'%(_time_fr10))
      fchain.write('=== time %d\n'%(_time_p))

      if _time_fr10 != _time_p:
        print '%d did not match with %d'%(_time_fr10,_time_p)
        sys.exit()
      for _index_clu,_size in enumerate(clu_size):
        _mem = clu_mem[_index_clu]
        _mem_conn = []
        for _member in _mem:
          _rela_conn = []
          for _neighbor in _conn[_member]:
            _rela_conn.append(_mem.index(_neighbor)) 
          _mem_conn.append(_rela_conn)
        cycles=[]
        anaTop.cycle_in_graph(_mem_conn, cycles)

        _long_short_chain = anaTop.long_short_path(_mem_conn)

        for _e,_cyl in enumerate(cycles):
          if _cyl !=[]:
            if len(_cyl) == 1:
              fcycle.write('cycle: ')
              _r = _cyl[0]
              for _m in _r:
                fcycle.write('%d '%_mem[_m])
              fcycle.write('\n')

            elif len(_cyl) >1:
              for _l,_m in enumerate(_cyl):
                fcycle.write('cycle: ')
                for _r,_s in enumerate(_m):
                  fcycle.write('%d '%_mem[_s])
                fcycle.write('\n')

        if _size<2:continue
        fchain.write('%d Mer : '%_size)
        for _cid in _long_short_chain:
          fchain.write('%d '%_mem[_cid])
        fchain.write('\n')
 

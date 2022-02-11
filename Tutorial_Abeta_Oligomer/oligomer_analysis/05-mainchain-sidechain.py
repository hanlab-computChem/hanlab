import sys
sys.path.append('../basic_code/')
import numpy as np
import os
import analysis_topology as anaTop 
import cal_cont_cluster as cal_cont

## Usage ##
'''
python2 05-mainchain-sidechain.py
'''

_frame = 450 ###ps
_start = 5050000
_end = 5100000

_largest_clu = 40
n_res = 40

_time = _start 
_time_fr10 , _time_p=_start,_start

_f10 = 'data/ab40-demo-cont10.xvg'
_f10_p = 'data/ab40-demo-cont10-olig.xvg'
_fcycle = 'data/ab40-demo-cont10-olig-cycle.xvg'
_fchain = 'data/ab40-demo-cont10-olig-chain.xvg'

_fmain_side = 'data/ab40-demo-cont10-olig-main-side.xvg'


_fr_p,_fr_id10,_fr_id,_fr_ctt = 0,0,0,0
_frame_id = 0

with open(_fmain_side,'w') as fo:

  while _time_fr10+_frame <_end:

    _cal_ctt = cal_cont.cal_cont('','', 100)

    _time_ctt, _conn,_fr_ctt = _cal_ctt.load_cluster_noself(_f10, _fr_ctt)
    _time_p, clu_size, clu_mem, _fr_p = _cal_ctt.load_patch(_f10_p, _fr_p)
    _time_fr10, clu_cycle,_fr_id10 = _cal_ctt.load_cycle(_fcycle, _fr_id10)
    _time_fr, clu_chain,_fr_id = _cal_ctt.load_main_chain(_fchain,_fr_id)

    if _time_fr10 != _time_fr:
      print '%d did not match with %d'%(_time_fr10,_time_fr)
      sys.exit()
  
    fo.write('===  time %d  ===\n'%_time_fr10)

    _extra_chain = [i for i in clu_chain]    
    _main_chain,_side_chain = [],[]
    #print 'clu_cycle',clu_cycle

    if clu_cycle !=[]:
      for _cycle in clu_cycle:  
        for _idx,_ref in enumerate(clu_mem): ## _idx :which oligomer a cycle belongs to
          if _cycle[0] in _ref:
            break

        _repeat = list(set(_cycle) & set(clu_chain[_idx]))
        if len(_repeat) > 0.5*len(_cycle):
          _extra_chain[_idx] += _cycle

      for _e,_chain in enumerate(_extra_chain):  
        _main_chain.append(np.unique(_chain))
        _side_chain.append(list(set(clu_mem[_e])-set(_main_chain[_e])))

        fo.write('%d mainchain : '%_e)
        for _l in _main_chain[_e]:
          fo.write('%d '%_l)
        fo.write('\n')

        fo.write('%d sidechain : '%_e)
        for _m in _side_chain[_e]:
          fo.write('%d '%_m)
        fo.write('\n')
    else:
      for _e,_chain in enumerate(clu_chain):  ## for each mainchain
        _sidechain,_main = [],[_c for _c in _chain]
        fo.write('%d mainchain : '%_e)
        for _l in _chain:
          fo.write('%d '%_l)

        _tmp = list(set(clu_mem[_e])-set(_chain))
        if _tmp !=[]:  ### partical not in cycle and also not in mainchain
          for _m in _tmp:
            if list(set(_conn[_m]) &set(_chain)) ==[]:
              _sidechain.append(_m)
              continue  ### contact with mainchain particle or not
            else:
              _main.append(_m)
              fo.write('%d '%_m)
        fo.write('\n')

        fo.write('%d sidechain : '%_e)
        for _r in _sidechain:
          fo.write('%d '%_r)
        fo.write('\n')

 

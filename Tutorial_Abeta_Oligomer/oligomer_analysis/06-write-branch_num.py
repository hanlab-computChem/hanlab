import numpy as np
import sys
sys.path.append('../basic_code/')
import analysis_topology as anaTop
import networkutil

## Usage ##
'''
python2 06-write-branch_num.py 
'''


def load_frame(f):

 _t = f.readline()

 if len(_t)==0:
  return (-1,[])

 _txt=[]
 if '===' in _t:
  _time = int(_t[:-1].split()[2])
 else:
  _time = -1

 while not '===' in _t:
  _txt.append(_t)
  _t=f.readline()
  if len(_t)==0: break
  if '===' in _t:
   _time = int(_t[:-1].split()[2])

 return (_time, _txt)

def _match(mem, l):

 _re=[]
 for _i in l:
  if _i in mem: _re.append(mem[_i])
 return _re

class oligo_top:

 def __init__(self, mem = "", conn = "", _time = 0):
  #print 'mem',mem
  _srl_m = mem.split()
  self.mem_id = np.array([int(_k) for _k in _srl_m])

  self.mem = dict([(int(_j), _i) for _i,_j in enumerate(_srl_m)])
  #print 'self.mem',self.mem

  _srl_c = ((conn.replace('], [','|')).replace('[[','')).replace(']]','').split('|')

  self.n = len(_srl_m)  
  self.conn=[]
  self.ne = 0

  for _i in _srl_c:
   _t = (_i.replace(',',' ')).split()
   _tc = []
   for _j in _t:
    _tc.append(self.mem[int(_j)])
   self.conn.append(_tc)
   self.ne+=len(_tc)

  self.time = _time

 def load_long_short_path(self, _path):
  _re = _match(self.mem, _path)
  if len(_path)!=len(_re):
   print 'Wrong path info', _path
   print 'for oligomer with members', self.mem.keys()
   print 'at time', self.time
   sys.exit()

  self.path = _re 
  self.path_mem = self.mem_id[_re]
  #print 'self.path_mem',self.path_mem

 def load_cycle(self, _all_cycle):
  _c = [[] for _i in xrange(self.n+1)]
  
  for _ic in _all_cycle:
   _re = _match(self.mem, _ic)
   if len(_re)==0: continue
   if len(_re)!=len(_ic):
    print 'Wrong cycle info', _ic
    print 'for oligomer with members', self.mem.keys()
    print 'at time', self.time
    sys.exit()
   
   _c[len(_re)].append(_re)

  self.cycles= _c
  #print 'self.cycles',self.cycles

 def cycle_cnt(self):
  _t=0
  for _i in  self.cycles: _t+=len(_i)
  return _t

 def cycle_size(self, ring_size):

  for _j,_i in enumerate(self.cycles):
   ring_size[_j]+= len(_i)

 def n_branches(self):
  _mark = self.mark_nm 
  _conn=[]
  for _k,_i in enumerate(self.conn):
   if _mark[_k]==1:
    _conn.append([])
    continue
   _t=[]
   for _j in _i:
    if _mark[_j]==1: continue
    _t.append(_j)
   _conn.append(_t)

  _re = networkutil.strong_connected_components(_conn)

#  print 'all connect',self.conn
#  print 'long short path', self.path
#  print 'cycles', self.cycles
#  print 'n in m marked',_mark
#  print 'connect excluding marked',_conn
#  print 'connected comp',_re
#  print
  return len(_re) - _mark.sum()

 def node_in_main(self):
  mark = np.zeros(self.n,dtype=np.int32)
  mark_p = np.zeros(self.n, dtype=np.int32)
  mark_c = np.zeros(self.n, dtype=np.int32)
  mark[self.path]=1
  mark_p[self.path]=1
  _cnt = mark.sum()
  for _i in xrange(self.n+1):
   _c = self.cycles[_i]
   for _j in _c:
    mark_c[:]=0
    mark_c[_j]=1
    _nj = len(_j)
#    print _j, mark_c, mark_p, mark
#    print mark_p*mark_c
    if (mark_p*mark_c).sum()*2>_nj: 
#     print _j, 'on chain'
     mark[_j]=1

#  print self.cycles
#  print self.path
#  print np.where(mark>0)
  self.mark_nm = mark
  return mark.sum()    

 def mainchain(self):   ### show the cluster member of node in mainchain
  mark = np.zeros(self.n,dtype=np.int32)
  mark_p = np.zeros(self.n, dtype=np.int32)
  mark_c = np.zeros(self.n, dtype=np.int32)
  mark[self.path]=1
  mark_p[self.path]=1
  _cnt = mark.sum()
  for _i in xrange(self.n+1):
   _c = self.cycles[_i]
   for _j in _c:
    mark_c[:]=0
    mark_c[_j]=1
    _nj = len(_j)
    if (mark_p*mark_c).sum()*2>_nj:
#     print _j, 'on chain'
     mark[_j]=1

#  print np.where(mark>0)
  self.main_chain = self.mem_id[np.where(mark>0)[0]]
  return

 def branchitem(self):
   #for i in self.mem_id:
   #  if i not self.main_chain:
   #    _out.append()
   self.branch_item = list(set(self.mem_id) - set(self.main_chain))
   self.branch_item.sort()

 def thickeritem(self):
   self.thicker_item = list(set(self.main_chain) -set(self.path_mem))
   self.thicker_item.sort()

 def e_to_v(self):
  return float(self.ne)/float(self.n)  

# main #
_f10 = 'data/ab40-demo-clu_mem-clu_conn.xvg'
_fcycle = 'data/ab40-demo-cont10-olig-cycle.xvg'
_fchain = 'data/ab40-demo-cont10-olig-chain.xvg'
_fbranch = 'data/ab40-demo-cont10-olig-branch-num.xvg'


# largest oligomer to be analyzed
n = 35
_frame = 450 ###ps
_nchain = 100
_start = 5050000  ## demo
_end = 5100000    ## demo

_time_list = []


_f =  open(_f10 , 'r')
_f_cycle = open(_fcycle, 'r')
_f_path = open(_fchain, 'r')

fo = open(_fbranch, 'w')

_show_time = 0
_t , _txt = load_frame(_f)
_t_cycle , _txt_cycle = load_frame(_f_cycle)
_t_path , _txt_path = load_frame(_f_path)

_time_list.append(_t)
if _t==_t_cycle and _t==_t_path:
 pass
else:
 print 'Error! mismatched time frames between',fnm,fnm_cycle,fnm_path
 print 'at',_t,_t_cycle, _t_path
 sys.exit()

while _t >=0:
 _time = _t

 _t, _txt = load_frame(_f)
 _t_cycle , _txt_cycle = load_frame(_f_cycle)
 _t_path , _txt_path = load_frame(_f_path)

 _time_list.append(_t)
 print _t

 if _t==_t_cycle and _t==_t_path:
  pass
 else:
  print 'Error! mismatched time frames between',fnm,fnm_cycle,fnm_path
  print 'at',_t,_t_cycle, _t_path
  sys.exit()

 _all_cycle = [[int(__j) for __j in (__i.replace('cycle:','')).split()] for __i in  _txt_cycle]


 if _txt ==[]: continue
 if _time< _start or _time> _end: continue

 _i_path = 0 
 fo.write('===time %d ===\n'%(_time_list[-2]))

 for _i in xrange(0,len(_txt),2):
  _o = oligo_top(mem = _txt[_i].replace('cluster member:',''), \
       	 conn = _txt[_i+1], _time = _time)
  _o.load_long_short_path([int(__t) for __t in _txt_path[_i_path].split()[3:]])
  _i_path+=1
  _o.load_cycle(_all_cycle)

  _nt = len(_o.mem)
  if _nt>n: continue

  _nm = _o.node_in_main()
  _o.mainchain()

  _n_brch = _o.n_branches()
  _o.branchitem()
  _o.thickeritem()

  fo.write('branch_num : %d\n'%_n_brch)

_f.close()
_f_cycle.close()
_f_path.close()

fo.close()





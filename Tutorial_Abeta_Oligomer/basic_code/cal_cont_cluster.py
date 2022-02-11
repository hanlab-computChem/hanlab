import os
import sys
sys.path.append('/share/home/hanw/pyscript/ff_op')
import numpy
import MDAnalysis
import geom_opt as geom
import networkutil
import time

class cal_cont(object):
  def __init__(self, top_nm, trj_nm, n_mono):
    self.top_nm = top_nm
    self.trj_nm = trj_nm
    self.n_mono = n_mono

    # contact parameters
    self.clash_cut = 4.5 # unit A
    self.cell_size = 9.0 # cell size
    self.clash_cut2 = self.clash_cut* self.clash_cut

    self._nstep = 100
    self._tot= 5000000 ## ps
    self._part = self._tot//self._nstep

  def find_patch(self,_conn):
    _re = networkutil.strong_connected_components(_conn)
    _rl = ''
    for _i in _re:
      _rl +='%d '%len(_i)+str(_i)+'\n'
    return _rl

  def final_cluster(self, fi):
    _load_frame=False
    _num, _out = 0 ,''
    with open(fi,'r') as f:
      _lines = f.readlines()

      #for rl in open(fi,'r'):
      for rl in _lines:
        if '===' in rl:
          _conn =[]
          _num = 0
          _load_frame=True
          continue
        if _load_frame:
          _num += 1
          _conn.append([int(_i) for _i in (rl[:-1].split(':')[1]).split()[0:]])
        if _num == self.n_mono:
          _out += self.find_patch(_conn)
    
    return _out

  def load_side_main_chain(self, _file,_line_start):
    _id = _line_start
    _load = True
    _f = open(_file)
    _lines = _f.readlines()
    _num_chain = 0
    nclu,nchain = 35,35
    mainchain, sidechain =numpy.ones([nclu,nchain],dtype=numpy.int)*-1,numpy.ones([nclu,nchain],dtype=numpy.int)*-1

    _time_line = _lines[_id]
    if 'time' in _time_line:
      if _id != len(_lines)-1:
        _time_fr = int(float(_time_line.split()[2]))
        _load = True
        _id += 1
      else:
        _time_fr = int(float(_time_line.split()[2]))
        _load = False

    _sym_m, _sym_s = -1,-1

    while _load:
      if _id == len(_lines):break
      _line = _lines[_id]
      if _line == '':
        continue
      elif 'time' in _line:
        #_id += 1
        _load = False
        break
      else:
        _out = []
      
        sl = _line.split()
        for sll in sl[3:]:
          _out.append(int(sll))
        _id += 1

        if 'mainchain' in _line:
          _sym_m += 1
          for _e,_o in enumerate(_out):
            mainchain[_sym_m,_e] = _o

        elif 'sidechain' in _line:
          _sym_s += 1
          for _e1,_o1 in enumerate(_out):
            sidechain[_sym_s,_e1] = _o1

    return _time_fr, sidechain, mainchain,_id

  def load_side_main_chain_100(self, _file,_line_start):  ### compared with load_side_main_chain(), the nclu is larger
    _id = _line_start
    _load = True
    _f = open(_file)
    _lines = _f.readlines()
    _num_chain = 0
    nclu,nchain = 100,35
    mainchain, sidechain =numpy.ones([nclu,nchain],dtype=numpy.int)*-1,numpy.ones([nclu,nchain],dtype=numpy.int)*-1

    _time_line = _lines[_id]
    if 'time' in _time_line:
      if _id != len(_lines)-1:
        _time_fr = int(float(_time_line.split()[2]))
        _load = True
        _id += 1
      else:
        _time_fr = int(float(_time_line.split()[2]))
        _load = False

    _sym_m, _sym_s = -1,-1

    while _load:
      if _id == len(_lines):break
      _line = _lines[_id]
      if _line == '':
        continue
      elif 'time' in _line:
        #_id += 1
        _load = False
        break
      else:
        _out = []

        sl = _line.split()
        for sll in sl[3:]:
          _out.append(int(sll))
        _id += 1

        if 'mainchain' in _line:
          _sym_m += 1
          for _e,_o in enumerate(_out):
            mainchain[_sym_m,_e] = _o

        elif 'sidechain' in _line:
          _sym_s += 1
          for _e1,_o1 in enumerate(_out):
            sidechain[_sym_s,_e1] = _o1

    return _time_fr, sidechain, mainchain,_id

  def load_main_chain(self, _file,_line_start):
    _id = _line_start
    _load = True
    _f = open(_file)
    _lines = _f.readlines()
    _num_chain = 0
    clu_all =[]
 
    _time_line = _lines[_id]
    if 'time' in _time_line:
      if _id != len(_lines)-1:
        _time_fr = int(float(_time_line.split()[2]))
        _load = True
        _id += 1
      else:
        _time_fr = int(float(_time_line.split()[2]))
        _load = False
    
    while _load:
      if _id == len(_lines):break
      _line = _lines[_id]
      if _line == '':
        continue
      elif 'time' in _line:
        #_id += 1
        _load = False
        break
      else:
        #=== time 15000007
        #cycle: 46 7 47
        _out = []
        sl = _line.split()
        for sll in sl[3:]:
          _out.append(int(sll))
        _id += 1
        clu_all.append(_out)

    return _time_fr,clu_all,_id

  def load_cycle(self, _file,_line_start):
    _id = _line_start
    _load = True
    _f = open(_file)
    _lines = _f.readlines()
    _num_chain = 0
    clu_all =[]
    
    _time_line = _lines[_id]
    if 'time' in _time_line:
      if _id != len(_lines)-1:
        _time_fr = int(float(_time_line.split()[2]))
        _load = True
        _id += 1
      else:
        _time_fr = int(float(_time_line.split()[2]))
        _load = False
        
    while _load :
      if _id == len(_lines):break

      _line = _lines[_id]
      if _line == '':
        continue
      elif 'time' in _line:
        #_id += 1
        _load = False
        break
      else:
        #=== time 15000007
        #cycle: 46 7 47
        _out = []
        sl = _line.split()
        for sll in sl[1:]:
          _out.append(int(sll))
        _id += 1
        clu_all.append(_out)

    return _time_fr,clu_all,_id

  def load_patch(self,_file,_line_start):   ### frame_num 9135 1 [4]
    _id = _line_start
    _load = True
    _f = open(_file)
    _lines = _f.readlines()
    _num_chain = 0
    clu_size, clu_mem =[],[] 
    
    while _load :
      _line = _lines[_id]
      if _line == '':
        continue
      elif 'time' in _line:
        _time_fr = int(float(_line.split()[1]))
        _first = True
        _id += 1
        continue
      else:
        _num_chain += int(_line.split()[0])
        if _num_chain == self.n_mono:
          _load = False

      _out = []
      sl = _line.replace('[','').replace(']','').replace(',','').split()
      clu_size.append(int(sl[0]))

      for _i in sl[1:]:
        _out.append(int(_i))        
      clu_mem.append(_out)
      _id += 1
    return _time_fr,clu_size, clu_mem,_id


  def load_cluster_noself_new(self,_file, _line_start):
    _id = _line_start
    _load = True
    _conn_noself = []
    _f = open(_file)
    _lines = _f.readlines()

    while _load :
      if _id == len(_lines):break   #### Here is different with load_cluster_noself
      else:
        _line = _lines[_id]
        #print '_line,_id',_line,_id
        if _line == '':
          continue
        elif 'Number' in _line:
          _id += 1
          continue
        elif 'time' in _line:
          _time_fr = int(float(_line.split()[1]))
          _id += 1
          continue
        elif int(_line.split(':')[0]) == self.n_mono-1:
          _load = False
        _d2 = []
        for _j in (_line.split(':')[1]).split():
          if int(_j) !=int(_line.split(':')[0]):
            _d2.append(int(_j))
        _conn_noself.append(_d2)
        _id += 1

    _f.close()
    return _time_fr, _conn_noself,_id

  def load_cluster_noself(self,_file, _line_start):
    _id = _line_start
    _load = True
    _conn_noself = []
    _f = open(_file)
    _lines = _f.readlines()

    while _load :
      _line = _lines[_id]
      if _line == '':
        continue
      elif 'Number' in _line:
        _id += 1
        continue
      elif 'time' in _line:
        _time_fr = int(float(_line.split()[1]))
        _id += 1
        continue
      elif int(_line.split(':')[0]) == self.n_mono-1:
        _load = False

      _d2 = []
      for _j in (_line.split(':')[1]).split():
        if int(_j) !=int(_line.split(':')[0]):
          _d2.append(int(_j))
      _conn_noself.append(_d2)
      _id += 1

    _f.close()
    return _time_fr, _conn_noself,_id


  def load_edge_oligo_label(self,_file, cluster_label,_line_start):
    _id = _line_start
    _load = True
    with open(_file) as _f:
      _lines = _f.readlines()
      _time_num = 0

      while _load :
        if _id == len(_lines):break
        _line = _lines[_id]
        #print _line
        if _line == '':
          continue
        if 'edge_atom1' in _line:
          _id += 1
          continue
        elif 'time' in _line:
          if _time_num ==0:
            _time_fr = int(float(_line.split()[2]))
            _time_num += 1
            _id += 1
            continue
          elif _time_num ==1:
            _load = False
            continue
        _oligo_id =  int(_line.split()[3])
        _clu_id = int(_line.split()[5])
        cluster_label[_oligo_id, _clu_id] +=1

        _id += 1

    #print '_line, _time_fr, _id',_line, _time_fr, _id 
    return _time_fr, cluster_label,_id

  def load_edge_oligo_label_diff(self,_file, _line_start):
    _id = _line_start
    _load = True
    _edge_nm, _clu_id,_oligomer_size = [],[],[]

    with open(_file) as _f:
      _lines = _f.readlines()
      _time_num = 0

      while _load :
        if _id == len(_lines):break
        _line = _lines[_id]
        #print _line
        if _line == '':
          continue
        if 'edge_atom1' in _line:
          _id += 1
          continue
        elif 'time' in _line:
          if _time_num ==0:
            _time_fr = int(float(_line.split()[2]))
            _time_num += 1
            _id += 1
            continue
          elif _time_num ==1:
            _load = False
            continue
        _oligomer_size.append(int(_line.split()[3]))
        _clu_id.append(int(_line.split()[5]))
        _edge_nm.append([int(i1) for i1 in _line.split(':')[0].split()])
        _id += 1
    #print '_line, _time_fr, _id',_line, _time_fr, _id 
    return _time_fr, _edge_nm, _clu_id,_oligomer_size,_id



  def load_edge_label(self,_file,_line_start):
    _id = _line_start
    _load = True
    _edge_label = []
    _edge = []

    with open(_file) as _f:
      _lines = _f.readlines()
      _time_num = 0

      while _load :
        if _id == len(_lines):break
        _line = _lines[_id]
        #print _line
        if _line == '':
          continue
        if 'edge_atom1' in _line:
          _id += 1
          continue
        elif 'time' in _line:
          if _time_num ==0:
            _time_fr = int(float(_line.split()[2]))
            _time_num += 1
            _id += 1
            continue
          elif _time_num ==1:
            _load = False
            continue
        _edge_label.append(int(_line.split()[5]))
        _new = []
        for _ele in _line.split()[:2]:
          _new.append(int(_ele))
        _edge.append(_new)

        _id += 1
    #print '_id',_id
     
    return _time_fr, _edge, _edge_label,_id

  def load_cluster(self,_file, _line_start):
    _id = _line_start
    _load = True
    _d1,_d2 = [],[]
    _f = open(_file)
    _lines = _f.readlines()

    while _load :
      if _id == len(_lines):break
      else:
        _line = _lines[_id]
        #print '_line,_id',_line,_id
        if _line == '':
          continue
        elif 'Num' in _line:
          _id += 1
          continue
        elif 'time' in _line:
          _time_fr = int(float(_line.split()[1]))
          _id += 1
          continue
        elif int(_line.split(':')[0]) == self.n_mono-1:
          _load = False
          
        for _j in (_line.split(':')[1]).split():
          _d1.append(int(_line.split(':')[0]))
          _d2.append(int(_j))
        _id += 1

    _f.close()
    _index = tuple([numpy.array(_d1),numpy.array(_d2)])
    return _time_fr, _index,_id


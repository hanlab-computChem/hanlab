import sys
sys.path.append('../basic_code/')
import numpy as np
import analysis_topology as anaTop
import networkutil
import geom_opt as geom
import sys

def print_itf_condense(__i, _itf,to_old,typ,typ_char):
 if _itf==[]: return
 print typ,typ_char,to_old[__i],':',
 for _a in _itf:
  print to_old[_a],
 print
  
def print_itf(_im, _itf, to_old, itf_typ):

 for _i in _itf:
  print 'interface mc',to_old[_im],':',to_old[_i],'sc','type',itf_typ


def cyc_cnt_ana(_count_allseg, _ring_size_allseg, _i_seg, _a,_b,n):
 _tot_o_cnt = _count_allseg[_i_seg,_a:_b+1].sum()
 if _tot_o_cnt>0:

  print "==cycle count for oligomer size", _a,"to",_b
  
  for _ci in xrange(2,n+1):
   _tot_cnt = _ring_size_allseg[_i_seg,_a:_b+1 , _ci].sum()
   print 'size',_ci,float(_tot_cnt)/float(_tot_o_cnt)

  print '==cycle ana end'


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

def load_map(_contactmap, _txt):

 for _rl in _txt:
  _srl = _rl[:-1].split()
  _contactmap[int(_srl[0]), int(_srl[1])] = int(_srl[3])
  _contactmap[int(_srl[1]), int(_srl[0])] = int(_srl[3])

def load_monoconf(_monoconf, _txt):
 for _rl in _txt:
  _srl=_rl[:-1].split()
  _monoconf[int(_srl[0])] = int(_srl[2])

def load_xcom(_xcom, _txt):
 _i = 0
 for _rl in _txt:
  _srl=_rl[:-1].split()
  for _j, _k in enumerate(_srl):
   _xcom[_i, _j] = float(_k)
  _i+=1


def print_distr(txt, _data, singleline = False, normalize=True,\
		showsum = False, print2d =False, hidezero = False, print3d=False):
 _s = _data.sum()

 if _s<=0: return
 if not normalize: _s = 1

 if print3d:
  print txt,':'
  _nx, _ny, _nz = _data.shape

  for _ix in xrange(_nx):
   for _iy in xrange(_ny):
    if _data[_ix, _iy].sum()==0 and hidezero: continue
    
    print _ix, _iy, _data[_ix, _iy].sum(),
    print float(_data[_ix, _iy].sum())/float(_s),':',
    for _iz in xrange(_nz):
     print float(_data[_ix,_iy,_iz])/float(_data[_ix,_iy].sum()),
    print
  print
  return



 if print2d:
  print txt,':'
  _nx, _ny = _data.shape

  for _ix in xrange(_nx):
   for _iy in xrange(_ny):
    if _data[_ix, _iy]==0 and hidezero: continue
    print _ix, _iy, _data[_ix, _iy],
    print float(_data[_ix, _iy])/float(_s)
  print
  return


 if singleline:
  print txt,':',
  if showsum: print _s,
 else:
  print txt,':'
 for _i in xrange(len(_data)):
  if singleline:
   print float(_data[_i])/float(_s),
  else:
   print _i, float(_data[_i])/float(_s)

 print

def _match(mem, l):

 _re=[]
 for _i in l:
  if _i in mem: _re.append(mem[_i])
 return _re

class oligo_top:

 def __init__(self, mem = "", conn = "", _time = 0):

  _srl_m = mem.split()

# old-->new
  self.mem = dict([(int(_j), _i) for _i,_j in enumerate(_srl_m)])

# new --> old
  self.to_old = [int(_j) for _j in _srl_m]

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
#  print mem
#  print conn
#  print self.time
#  print self.mem
#  print self.conn
#  print

#  self.cycles = []
#  ana_top.cycle_in_graph(self.conn, self.cycles)
#  self.path = ana_top.long_short_path(self.conn)  
#  print self.cycles
#  print self.path


 def substruct_typ(self, _substruct_typ, anaBranch=False,\
		   _brch_pt_typ=None, _monoconf=None,\
		   _brch_pt_typ_conf=None, printIdx=False,\
		   countNode=False, _brch_node_typ=None):
  _n = len(self.conn)
  _r_edge = np.zeros((_n,_n), dtype=np.int32)
  _in_r = np.zeros(_n, dtype=np.int32)
  _r_edge_list = [[[] for _j in xrange(_n)] for _i in xrange(_n)]


  #consider only triangle ring here
  #this hypothesis dont work
#  _triangle = self.cycles[3]

  _triangle=[]
  for _rr in self.cycles:
   for _rrr in _rr: _triangle.append(_rrr)

#  if len(_triangle)>len(self.cycles[3]):
#   print 'more polygon than triangle',len(_triangle), len(self.cycles[3])

  _nr = len(_triangle)
  for __i,_ri in enumerate(_triangle):
   
    for _i in xrange(len(_ri)):
     _a = _ri[_i]
     _in_r[_a]=1
     if _i+1>=len(_ri):
      _b = _ri[0]
     else:
      _b = _ri[_i+1]

     _r_edge_list[_a][_b].append(__i)
     _r_edge_list[_b][_a].append(__i)

     _r_edge[_a,_b]+=1
     _r_edge[_b,_a]+=1


  _conn = [[] for _i in xrange(_nr)]

  for _a in xrange(_n):
   for _b in xrange(_n):
    if _a<_b:
     _e = _r_edge_list[_a][_b]
     _nsh = len(_e)
     if _nsh>1:
      for __j in xrange(_nsh-1):
        _conn[_e[__j]].append(_e[__j+1])
        _conn[_e[__j+1]].append(_e[__j])


  _conn_s = [networkutil.sort_unq(_i) for _i in _conn]

#  print _conn_s
#  if _conn_s==[]: _conn_s=[[]]
  if _conn_s==[]:
   _re = [] 
  else:
   _re = networkutil.strong_connected_components(_conn_s)

#  print _re

#  for _i in _re:
#   for _j in _i:
#    for _k in _triangle[_j]:
#     print self.to_old[_k],
#   print

###
# identify main chain from longest shortest path
# for branch analysis
####
  _mark_mc = np.zeros(_n,dtype=np.int32)
  _mark_ele_in_mc = np.zeros(len(_re), dtype=np.int32)
  _path=self.path
  _mark_path= np.zeros(_n, dtype=np.int32)
  for _i in _path: _mark_path[_i]=1

  if anaBranch:

   
   for __ir, _i in enumerate(_re):
    __t=[]
    for _j in _i:
     for _k in _triangle[_j]:
      __t.append(_k)
    ___t = networkutil.sort_unq(__t)

    _on_path = _mark_path[___t].sum()
    if (len(___t)==1 and _on_path==1) or _on_path>1:
     _mark_ele_in_mc[__ir]=1
     _mark_mc[___t]=1
   _mark_mc[_path]=1
   
#  print 'mc',
#  for _i in xrange(_n):
#   if _mark_mc[_i]==1:
#    print self.to_old[_i],
#  print   
#   
#  for __ir, _i in enumerate(_re):
#   if _mark_ele_in_mc[__ir]==1:
#    for _j in _i:
#     for _k in _triangle[_j]:
#      print self.to_old[_k],
#    print


#search for merged or single ring
  _n_brch = 0
#this label make sure when doing branch analysis no single sc nodes
#are analyzed more than once
  _label = np.zeros(_n, dtype=np.int32)
#this label distinguish nodes in an ring system from those outside 
  _mark = np.zeros(_n, dtype=np.int32)
  for __ir,_i in enumerate(_re):
   _mark[:]=0 

   for _j in _i:
    for _k in _triangle[_j]:
     _mark[_k] =1
#     print self.to_old[_k],
#   print

   _r_s = _mark.sum()


   _c_l = 0
   _c_r = 0

   for _im in xrange(_n):
    if anaBranch:
     _conf = _monoconf[self.to_old[_im]]
    if _mark[_im]==1:
     _c = self.conn[_im]

     _n_r = 0
     _n_l = 0

# shared edge
     _n_sh = 0
# ring edge
     _n_re = 0




    
     for _il in _c:
      if _mark[_il]==0:
       _ve = _r_edge[_im, _il]
 
       if _ve>1:
        _n_sh+=1


       if _ve==0:
        _n_l+=1
       else:
        _n_re+=1
        _n_r+=1
        _n_r+=(1- _ve)




     _n_r/=2

     _c_l+= _n_l
     _c_r+= _n_r
     
     if anaBranch:
     
      if _mark_ele_in_mc[__ir]==1:
# sidechain
       _n_r_mc = 0
       _n_l_mc = 0
       _n_r_sc = 0
       _n_l_sc = 0
# shared edge
       _n_sh_mc = 0
       _n_sh_sc = 0
# ring edge
       _n_re_mc = 0
       _n_re_sc = 0


       for _il in _c:
        if _mark[_il]==1:
         _ve = _r_edge[_im, _il]

         if _ve>1:
          _n_sh_mc+=1


         if _ve==0:
          _n_l_mc+=1
         else:
          _n_re_mc+=1
          _n_r_mc+=1
          _n_r_mc+=(1- _ve)



       _n_r_mc/=2

       _itf_l=[]
       _itf_r=[]

       for _il in _c:
        if _mark[_il]==0 and _mark_mc[_il]==0 and _label[_il]==0:
         #print self.to_old[_im],self.to_old[_il]
         _label[_il]=1
         _ve = _r_edge[_im, _il]

         if _ve>1:
          _n_sh_sc+=1


         if _ve==0:
          _n_l_sc+=1
          _itf_l.append(_il)
         else:
          _itf_r.append(_il)
          _n_re_sc+=1
          _n_r_sc+=1
          _n_r_sc+=(1- _ve)




       _n_r_sc/=2
       
        #  mc |>o-  sc
       if _n_l_sc>0 and _n_sh_mc==0:
        _n_brch+= _n_l_sc
        _brch_pt_typ[2]+=_n_l_sc
        _brch_pt_typ_conf[_conf, 2]+= _n_l_sc
        if countNode: _brch_node_typ[2]+=1
 
#        print 'mc |>o-  sc',self.to_old[_im]
        if printIdx: print_itf(_im, _itf_l, self.to_old, 2)

       #   mc |>0<| or E| sc
       if _n_r_sc>0 and _n_sh_mc==0:
        _n_brch+= _n_r_sc
        _brch_pt_typ[3]+=_n_r_sc
        _brch_pt_typ_conf[_conf, 3]+= _n_r_sc
        if countNode: _brch_node_typ[3]+=1
#        print 'mc |>0<| or E| sc',self.to_old[_im]
        if printIdx: print_itf(_im, _itf_r, self.to_old, 3)

        #  mc |3o-  sc
       if _n_l_sc>0 and _n_sh_mc>0:
        _n_brch+= _n_l_sc
        _brch_pt_typ[4]+=_n_l_sc
        _brch_pt_typ_conf[_conf, 4]+= _n_l_sc
        if countNode: _brch_node_typ[4]+=1
#        print 'mc |3o-  sc', self.to_old[_im]
        if printIdx: print_itf(_im, _itf_l, self.to_old, 4)

       #   mc |3o<| or E| sc
       if _n_r_sc>0 and _n_sh_mc>0:
        _n_brch+= _n_r_sc
        _brch_pt_typ[5]+=_n_r_sc
        _brch_pt_typ_conf[_conf, 5]+= _n_r_sc
        if countNode: _brch_node_typ[5]+=1
#        print 'mc |3o<| or E| sc', self.to_old[_im]
        if printIdx: print_itf(_im, _itf_r, self.to_old, 5)


#   print 'num ring connections',_c_r,'num lin connections', _c_l

   if _r_s<3 or _r_s>10 : continue
   if _c_l+_c_r > 10: continue
   _substruct_typ[_r_s, _c_l+_c_r]+=1 
#   print 'size',_r_s,'conn',_c_l+_c_r

  for _in in xrange(_n):
   if _in_r[_in]==0:
    _c = self.conn[_in]
    _n_l = len(_c)
#    print 'node',self.to_old[_in],'conn',_n_l
    _substruct_typ[1, _n_l]+=1

  if anaBranch:
   for _im in _path:
    _conf = _monoconf[self.to_old[_im]]
    _c = self.conn[_im]
# sidechain
    _n_r_sc = 0
    _n_l_sc = 0
# shared edge
    _n_sh_sc = 0
# ring edge
    _n_re_sc = 0

    _itf_l=[]
    _itf_r=[]

    for _il in _c:
       if _mark_mc[_il]==0 and _label[_il]==0:
         #print self.to_old[_im],self.to_old[_il]
         _label[_il]=1
         _ve = _r_edge[_im, _il]

         if _ve>1:
          _n_sh_sc+=1


         if _ve==0:
          _n_l_sc+=1
          _itf_l.append(_il)
         else:
          _itf_r.append(_il)
          _n_re_sc+=1
          _n_r_sc+=1
          _n_r_sc+=(1- _ve)

    _n_r_sc/=2

    #   mc >o-  sc
    if _n_l_sc>0:
     _n_brch+= _n_l_sc
     _brch_pt_typ[0]+=_n_l_sc
     _brch_pt_typ_conf[_conf, 0]+= _n_l_sc
     if countNode: _brch_node_typ[0]+=1
#     print 'mc >o-  sc',self.to_old[_im]
     if printIdx: print_itf(_im, _itf_l, self.to_old, 0)

    #  mc >o<| or E| sc
    if _n_r_sc>0:
     _n_brch+= _n_r_sc
     _brch_pt_typ[1]+=_n_r_sc
     _brch_pt_typ_conf[_conf, 1]+= _n_r_sc
     if countNode: _brch_node_typ[1]+=1
#     print 'mc >o<| or E| sc',self.to_old[_im]
     if printIdx: print_itf(_im, _itf_r, self.to_old, 1)

#  print 'branch cnt', _n_brch
  if anaBranch: return _n_brch


 def coord_typ(self, _monoconf, _coord_typ, _coord_typ_conf, printIdx=False):
  _n = len(self.conn)
  _r_edge = np.zeros((_n,_n), dtype=np.int32)
  _in_r = np.zeros(_n, dtype=np.int32)

  
  for _rs in  self.cycles:
   for _ri in _rs:
    
    for _i in xrange(len(_ri)):
     _a = _ri[_i]
     _in_r[_a]=1
     if _i+1>=len(_ri):
      _b = _ri[0]
     else:
      _b = _ri[_i+1]

     _r_edge[_a,_b]+=1
     _r_edge[_b,_a]+=1
     

  for __i, __c in enumerate(self.conn):
   _conf = _monoconf[self.to_old[__i]]
   if _in_r[__i]==0:
   # type 0, non ring with 2 connectivity  -o-
    if len(__c)==2:
     _coord_typ[0]+=1
     _coord_typ_conf[_conf,0]+=1
     if printIdx: print_itf_condense(__i, __c, self.to_old, 0,'l')
   #type 1 , -O-
   #          |
    if len(__c)>2:
     _coord_typ[1]+=1
     _coord_typ_conf[_conf,1]+=1
     if printIdx: print_itf_condense(__i, __c,self.to_old,1,'l')
   else:
    _rc = 0
    _re = 0
# ring bond
    _n_r = 0
#linear bond
    _n_l = 0
# shared edge
    _n_sh = 0
# ring edge
    _n_re = 0

    _itf_l=[]
    _itf_r=[]
    for __k in __c:
     _rc+= _in_r[__k]
     _re+= _r_edge[__i, __k]

     _ve = _r_edge[__i, __k]
     if _ve>1:
        _n_sh+=1
     if _ve==0:
        _n_l+=1
        _itf_l.append(__k)
     else:
        _itf_r.append(__k)
        _n_re+=1
        _n_r+=1
        _n_r+=(1- _ve)

    _n_r/=2

#    print 'node',self.to_old[__i],
    #  type 2: |>o
    if len(__c)==2 and _rc==2:
     _coord_typ[2]+=1
     _coord_typ_conf[_conf, 2]+=1
     if printIdx: print_itf_condense(__i, _itf_l,self.to_old,2,'l')
     if printIdx: print_itf_condense(__i, _itf_r,self.to_old,2,'r')
#     print 'type |>o',
    # type3 : |>o-
    if _n_re==2 and len(__c)>2:
     _coord_typ[3]+=1
     _coord_typ_conf[_conf, 3]+=1
     if printIdx: print_itf_condense(__i, _itf_l,self.to_old,3,'l')
     if printIdx: print_itf_condense(__i, _itf_r,self.to_old,3,'r')

#     print 'type |>o-',

    #type 4:     o
    #           <|>
    if _n_r==1 and _n_sh>0 and _n_l==0:
     _coord_typ[4]+=1
     _coord_typ_conf[_conf, 4]+=1
     if printIdx: print_itf_condense(__i, _itf_l,self.to_old,4,'l')
     if printIdx: print_itf_condense(__i, _itf_r,self.to_old,4,'r')

#     print 'type oE',
    #type5:        |
    #              o
    #             <|>
    if _n_r==1 and _n_sh>0 and _n_l>0:
     _coord_typ[5]+=1
     _coord_typ_conf[_conf, 5]+=1
     if printIdx: print_itf_condense(__i, _itf_l,self.to_old,5,'l')
     if printIdx: print_itf_condense(__i, _itf_r,self.to_old,5,'r')

#     print 'type -oE',
    #type6:      |>o<|
    if _n_r>1:
     _coord_typ[6]+=1
     _coord_typ_conf[_conf, 6]+=1
     if printIdx: print_itf_condense(__i, _itf_l,self.to_old,6,'l')
     if printIdx: print_itf_condense(__i, _itf_r,self.to_old,6,'r')

#     print 'type |>o<|',

#    print
 
 
  


    


 def coord_ring_geom(self, _xcom, _ang_ring_distr, _dih_ring_distr):
  for _rsize in xrange(3,min(7, len(self.cycles))):
 
   _a_distr = _ang_ring_distr[_rsize]
   _d_distr = _dih_ring_distr[_rsize]

   _cyc = self.cycles[_rsize]
#   print _cyc

   for _ring in _cyc:
    _old_r = [self.to_old[_i] for _i in _ring]
    _c_r = _old_r + _old_r
#    print _ring
#    print _c_r  
    for _k in xrange(_rsize): 
     _a = geom.angle(_xcom[_c_r[_k]], _xcom[_c_r[_k+1]], _xcom[_c_r[_k+2]])
     _idx = int(_a/10.0)
     _a_distr[_idx]+=1
#     print _c_r[_k], _c_r[_k+1], _c_r[_k+2]

    if _rsize<4: continue
   
    for _k in xrange(_rsize):
     _a = geom.dihedral(_xcom[_c_r[_k]], _xcom[_c_r[_k+1]], _xcom[_c_r[_k+2]], _xcom[_c_r[_k+3]])
     _idx = int((_a + 180.0)/10.0)
     _d_distr[_idx]+=1


 def coord_mainpath_geom(self, _xcom, _ang_distr, _dih_distr):
  path = [self.to_old[_i] for _i in self.path]
#  print path
  _lpath = len(path)
  for _i in xrange(_lpath-1):
   if _i+2 < _lpath:
    _a = geom.angle(_xcom[path[_i]], _xcom[path[_i+1]], _xcom[path[_i+2]])
    _idx = int(_a/10.0)
    _ang_distr[_idx]+=1

   if _i+3 < _lpath:
    _a = geom.dihedral(_xcom[path[_i]], _xcom[path[_i+1]], _xcom[path[_i+2]], _xcom[path[_i+3]])
    _idx = int((_a + 180.0)/10.0)
    _dih_distr[_idx]+=1


 def mainpath_geom(self, _xcom, _l_sq, _r_sq, _l_sq_cnt, _r_sq_cnt):
  path = [self.to_old[_i] for _i in self.path]
  _len = len(path)

  _r = geom.bond(_xcom[path[0]], _xcom[path[-1]])

  _r_sq[_len]+= _r*_r
  _r_sq_cnt[_len]+=1

  for __i in xrange(_len-1):
   _l = geom.bond(_xcom[path[__i]], _xcom[path[__i+1]])
   _l_sq[_len] += _l*_l
   _l_sq_cnt[_len]+=1

 def coord_D_geom(self, _xcom, _monoconf,\
                    _ang_dih_D_I_distr, _ang_dih_D_II_distr,\
                    _ang_dih_D_conf_I_distr, _ang_dih_D_conf_II_distr,\
                    _ang_D_distr, _ang_D_conf_distr):

   _n = len(self.conn)

   for __i in xrange(_n):
    _c_id = self.to_old[__i]
    _conn = self.conn[__i]

    _D = len(_conn)

    _conf = _monoconf[_c_id]
    
#    print __i, _c_id
#    print _conn

    _to_sort_list = [[self.to_old[__j], geom.bond(_xcom[_c_id], _xcom[self.to_old[__j]])]\
		     for __j in _conn]

#    print _to_sort_list

    _to_sort_list.sort(key = lambda x:x[1])

#    print _to_sort_list
    _nid = [_x[0] for _x in _to_sort_list]


    if _D>1:
     #  b-a-c
     _a  = geom.angle(_xcom[_nid[0]], _xcom[_c_id], _xcom[_nid[1]])
     _idx = int(_a/10.0)
     _ang_D_distr[_D,_idx]+=1
     _ang_D_conf_distr[_D, _conf, _idx]+=1
     _a_1 = _a
     if _D>2:
     #  b-a-c
     #    |
     #    d
      _a = geom.angle(_xcom[_nid[0]], _xcom[_c_id], _xcom[_nid[2]])

      if _a < _a_1:
       _a = _a_1
       _d = geom.dihedral(_xcom[_nid[1]], _xcom[_c_id], _xcom[_nid[0]], _xcom[_nid[2]])
      else:
       _d = geom.dihedral(_xcom[_nid[2]], _xcom[_c_id], _xcom[_nid[0]], _xcom[_nid[1]])
      _idxa = int(_a/15.0)
      _idxd = int((_d + 180.0)/15.0)
      _ang_dih_D_I_distr[_D, _idxa,_idxd]+=1
      _ang_dih_D_conf_I_distr[_D, _conf, _idxa,_idxd]+=1
      
      if _D>3:
      #     e
      #  b--a--c
      #     d

       _a = geom.angle(_xcom[_nid[0]], _xcom[_c_id], _xcom[_nid[3]])
       _d = geom.dihedral(_xcom[_nid[0]], _xcom[_c_id], _xcom[_nid[3]], _xcom[_nid[1]])
       _idxa = int(_a/15.0)
       _idxd = int((_d + 180.0)/15.0)
       _ang_dih_D_II_distr[_D, _idxa,_idxd]+=1
       _ang_dih_D_conf_II_distr[_D, _conf, _idxa,_idxd]+=1


     




 def load_long_short_path(self, _path):
  _re = _match(self.mem, _path)
  if len(_path)!=len(_re):
   print 'Wrong path info', _path
   print 'for oligomer with members', self.mem.keys()
   print 'at time', self.time
   sys.exit()

  self.path = _re
#  print self.path
#  for _i in self.path: print self.to_old[_i],
#  print

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
#  print self.cycles

 def cycle_cnt(self):
  _t=0
  for _i in  self.cycles: _t+=len(_i)
  return _t

 def cycle_size(self, ring_size, countOnce=False):

  for _j,_i in enumerate(self.cycles):
   if countOnce:
     if len(_i)>0:
      ring_size[_j]+= 1
   else:
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

 def node_in_main(self, countRing=False):
  mark = np.zeros(self.n,dtype=np.int32)
  mark_p = np.zeros(self.n, dtype=np.int32)
  mark_c = np.zeros(self.n, dtype=np.int32)
  mark[self.path]=1
  mark_p[self.path]=1
  _ring_cnt = 0
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
     _ring_cnt+=1
#  print self.cycles
#  print self.path
#  print np.where(mark>0)
  self.mark_nm = mark
  if countRing:
   return _ring_cnt
  else:
   return mark.sum()    

 def e_to_v(self):
  return float(self.ne)/float(self.n)  

 def edge_conf_type_ana(self,_contactmap,_edge_all,  _edge_in_ring, _edge_out_ring,\
                    _edge_in_nodemain, _edge_in_path,\
                    _edge_between_mc_sc, _edge_between_path_nm,\
                     _monoconf, _conf_in_ring, _conf_out_ring,\
                     _conf_in_nodemain, _conf_in_path, _conf_in_sc,\
                     _conf_in_path_ter,\
                     _conf_pair_in_path,\
		     _conf_pair_between_mc_sc, _conf_pair_between_mc_nm,\
		     showBond=False):

  _to_old = self.to_old
  mark_p = np.zeros(self.n, dtype=np.int32)
  mark_ring = np.zeros(self.n, dtype=np.int32)
  mark_c = np.zeros((self.n, self.n), dtype=np.int32)
  mark_p[self.path]=1
  mark_nm = self.mark_nm

  for _i in xrange(self.n+1):
   _c = self.cycles[_i]
   for _j in _c:
    mark_ring[_j]=1
    _f = _j[:-1]
    _s = _j[1:]
    _f_2 = _j[-1]
    _s_2 = _j[0]
    mark_c[_f, _s]=1
    mark_c[_s, _f]=1
    mark_c[_f_2, _s_2]=1
    mark_c[_s_2, _f_2]=1


#  print 'conn', self.conn
#  print 'to old', _to_old
#  print 'path', self.path
#  print 'clcyles', self.cycles

#  print 'mark p',mark_p
#  print 'mark nm', mark_nm
#  print 'mark c'
#  print mark_c

  for _i in xrange(len(self.conn)):
   _conf_i = _monoconf[_to_old[_i]]
   if _conf_i<0:
    print 'Undefined conf for mono', _to_old[_i]
    sys.exit()

   if mark_ring[_i]==1:
    _conf_in_ring[_conf_i]+=1
   else:
    _conf_out_ring[_conf_i]+=1

   if mark_nm[_i]==1:
    _conf_in_nodemain[_conf_i]+=1
   else:
    _conf_in_sc[_conf_i]+=1

   if mark_p[_i]==1:
    _conf_in_path[_conf_i]+=1

   if _i == self.path[0] or _i == self.path[-1]:
    _conf_in_path_ter[_conf_i]+=1

  for _i, _j in enumerate(self.conn):
   _conf_i = _monoconf[_to_old[_i]]
   for _k in _j:
    if _i<_k:
     _ct = _contactmap[_to_old[_i],_to_old[_k]]
     _conf_k = _monoconf[_to_old[_k]]

     if _ct<0:
      print '-1 type found between monomers',_to_old[_i], _to_old[_k]
      sys.exit()
     _edge_all[_ct] +=1
   
#  edge in/out ring
     if mark_c[_i,_k]==1:
      _edge_in_ring[_ct]+=1
     else:
      _edge_out_ring[_ct]+=1
#  edge in main path
     if mark_p[_i]==1 and mark_p[_k]==1:
      _edge_in_path[_ct]+=1
      _conf_pair_in_path[min(_conf_i, _conf_k),max(_conf_i, _conf_k), _ct]+=1
# edge in main chain
     if mark_nm[_i]==1 and mark_nm[_k]==1:
      _edge_in_nodemain[_ct]+=1

# edge between sidechain and mainchain
     if mark_nm[_i] + mark_nm[_k]==1:
      _edge_between_mc_sc[_ct]+=1      
      if mark_nm[_i]==1:
       _mc_c = _conf_i
       _sc_c = _conf_k
       _mc_i = _to_old[_i]
       _sc_i = _to_old[_k]
      else:
       _mc_c = _conf_k
       _sc_c = _conf_i
       _mc_i = _to_old[_k]
       _sc_i = _to_old[_i]

      _conf_pair_between_mc_sc[_mc_c, _sc_c, _ct]+=1
      if showBond:
       print 'mc', _mc_i,'to sc',_sc_i

# edge between main path and other partis in node main
     if mark_p[_i]+mark_p[_k]==1 and mark_nm[_i]+mark_nm[_k]==2:
      _edge_between_path_nm[_ct]+=1

      if mark_p[_i]==1:
       _mc_c = _conf_i
       _nm_c = _conf_k
       _mc_i = _to_old[_i]
       _nm_i = _to_old[_k]
      else:
       _mc_c = _conf_k
       _nm_c = _conf_i
       _mc_i = _to_old[_k]
       _nm_i = _to_old[_i]

      _conf_pair_between_mc_nm[_mc_c, _nm_c, _ct]+=1
      if showBond:
       print 'mc', _mc_i,'to nm',_nm_i


#  print _edge_all, 'all'
#  print _edge_in_ring, 'in ring'
#  print _edge_out_ring, 'not in ring'
#  print _edge_in_path, 'main path'
#  print _edge_in_nodemain, 'node main'
#  print _edge_between_mc_sc, 'between mc and sc'
#  print _edge_between_path_nm, 'between main path and other part of node main'


 def bond_type_between_conf_pair(self, _contactmap, _monoconf, \
                                   _bond_distr_by_conf_pair):
  _to_old = self.to_old


  for _i, _j in enumerate(self.conn):
   _conf_i  = _monoconf[_to_old[_i]]

   if _conf_i<0:
    print 'Undefined conf for mono', _to_old[_i]
    sys.exit()


   for _k in _j:

    _conf_k  = _monoconf[_to_old[_k]]

    if _conf_k<0:
     print 'Undefined conf for mono', _to_old[_k]
     sys.exit()


#    if _i<_k:
    _ct = _contactmap[_to_old[_i],_to_old[_k]]

    if _ct<0:
      print '-1 type found between monomers',_to_old[_i], _to_old[_k]
      sys.exit()


    _bond_distr_by_conf_pair[_conf_i, _conf_k, _ct]+=1


 def edge_corr_path(self, _contactmap,_edge_pair_path):

  _path = self.path
  _to_old = self.to_old

  if len(_path)<=2: return

  for _i in xrange(2,len(_path)):
   
   _ti = _contactmap[_to_old[_path[_i]], _to_old[_path[_i-1]]]
   _tj = _contactmap[_to_old[_path[_i-2]], _to_old[_path[_i-1]]]
   if _ti<0:
    print 'non exist contact between',_to_old[_path[_i]],'and',_to_old[_path[_i-1]]

   if _tj<0:
    print 'non exist contact between',_to_old[_path[_i-2]],'and',_to_old[_path[_i-1]]

   _edge_pair_path[min(_ti,_tj),max(_ti,_tj)]+=1


 def bond_type_multivalency(self, _contavtmap, _cor_num_distr_by_type,_count_node_by_type):
  _conn = self.conn
  _to_old = self.to_old

  for _k,_i in enumerate(_conn):
   _l = len(_i)
   _tk = _to_old[_k]
   for _j in _i:
    _tj = _to_old[_j]
    _ct = _contactmap[_tk, _tj]
    if _ct<0:
     print 'Error! contact between',_tk,'and',_tj,'does not exist'
     sys.exit()

    _cor_num_distr_by_type[_ct,_l]+=1
    _count_node_by_type[_ct]+=1
   
 def monoconf_edge_couple(self, _contactmap, _monoconf, _cor_num_distr_by_monoconf,\
                            _cnt_bondtype_by_monoconf):

  _conn = self.conn
  _to_old = self.to_old

  for _k,_i in enumerate(_conn):

   _l = len(_i)
   _tk = _to_old[_k]

   _confk = _monoconf[_tk]
   if _confk<0:
    print 'Error: undefined conf state for mono', _tk
    sys.exit()

   _cor_num_distr_by_monoconf[_confk, _l]+=1

   for _j in _i:
    _tj = _to_old[_j]
    _ct = _contactmap[_tk, _tj]

    _confj = _monoconf[_tj]
    if _confj<0:
     print 'Error: undefined conf state for mono', _tj
     sys.exit()

    if _ct<0:
     print 'Error! contact between',_tk,'and',_tj,'does not exist'
     sys.exit()

    _cnt_bondtype_by_monoconf[_confk, _ct]+=1

#end of oligomer obj


# main #
#name of files of connectivity for analysis
# this is for ab40
#fnm_list = ['ab40-traj1-clu-info-re/ab40-traj1-clu_mem-clu_conn-larger-than-trimer.xvg',\
#		'ab40-traj2-clu-info-re/ab40-traj2-clu_mem-clu_conn-larger-than-trimer.xvg']
#fnm_cycle_list = ['ab40-traj1-clu-info-re/ab40-traj1-clu_cycle.xvg',\
#		  'ab40-traj2-clu-info-re/ab40-traj2-clu_cycle.xvg']
#fnm_path_list = ['ab40-traj1-clu-info-re/ab40-traj1-clu_main_chain.xvg',\
#		 'ab40-traj2-clu-info-re/ab40-traj2-clu_main_chain.xvg']
#fnm_contact_list = ['ab40-traj1-clu-info-re/ab40-traj1-contact-map-label-k5-b5000017-e14961982-real-right-rename.xvg',\
#		    'ab40-traj2-clu-info-re/ab40-traj2-contact-map-label-k5-b5000017-e14961982-real-right-rename.xvg']
#fnm_monoconf_list = ['ab40-intramolecular-cluster/ab40-traj1-intramolecular-ctt-map-ave300fr-b5000017-e14961982-each_fr.xvg',\
#		     'ab40-intramolecular-cluster/ab40-traj2-intramolecular-ctt-map-ave300fr-b5000017-e14961982-each_fr.xvg']
#fnm_xcom_list = ['COM/ab40-traj1-com-cor.xvg','COM/ab40-traj2-com-cor.xvg']
#exclu_list = [5666490,14017612]
#test_list= [10000507]
# largest oligomer to be analyzed
#n=23


# this is for ab42
#fnm_list = ['ab42-traj1-clu-info-re/ab42-traj1-clu_mem-clu_conn-larger-than-trimer.xvg',\
#               'ab42-traj2-clu-info-re/ab42-traj2-clu_mem-clu_conn-larger-than-trimer.xvg']
#fnm_cycle_list = ['ab42-traj1-clu-info-re/ab42-traj1-clu_cycle.xvg',\
#                  'ab42-traj2-clu-info-re/ab42-traj2-clu_cycle.xvg']
#fnm_path_list = ['ab42-traj1-clu-info-re/ab42-traj1-clu_main_chain.xvg',\
#                 'ab42-traj2-clu-info-re/ab42-traj2-clu_main_chain.xvg']
#fnm_contact_list =['ab42-traj1-clu-info-re/ab42-traj1-contact-map-label-k5-b5000017-e14961982-real-right-rename.xvg',\
#		   'ab42-traj2-clu-info-re/ab42-traj2-contact-map-label-k5-b5000017-e14961982-real-right-rename.xvg']
#fnm_monoconf_list = ['ab42-intramolecular-cluster/ab42-traj1-intramolecular-ctt-map-ave300fr-b5000017-e14961982-each_fr.xvg',\
#                     'ab42-intramolecular-cluster/ab42-traj2-intramolecular-ctt-map-ave300fr-b5000017-e14961982-each_fr.xvg']
#fnm_xcom_list = ['COM/ab42-traj1-com-cor.xvg','COM/ab42-traj2-com-cor.xvg']
#exclu_list = [8257297, 12844552, 13764914]
# largest oligomer to be analyzed
#n=25


#explanation
# fnm_list: list of files contain connectivity info of all oligomers in a frame using adsolute ids of chains
# fnm_cycle_list: info on all possible cycles in networks in a frame
# fnm_path: info on chains on the longest shortest paths for each oligomers in a frame
# fnm_contact: info on contact state between those chains that form physical contact. contact state is defined by clustering contact maps of interfacial contacts
# fnm_monoconf: conf state of each chain in a frame. conf state is defined by using clustering based on intrachain contact map
# fnm_xcom: coordinate of COM of each chain in a frame
# exclu_list excludes frames (by time) that leads to unstable numerical results and crash calculation in rare cases  
# as an example, a.xvg-f.xvg contain all infro regarding a single frame needed for all the analysis conducted with this code
# infor more than one frames can be included in these files and be handled by this code
# but all info in these files should be synchronized.

#for debugging
fnm_list=['a.xvg']
fnm_cycle_list = ['b.xvg']
fnm_path_list =  ['c.xvg']
fnm_contact_list = ['d.xvg']
fnm_monoconf_list = ['e.xvg']
fnm_xcom_list = ['f.xvg']
exclu_list=[]
n = 23


#here parameters for ab42

#for contact map
_n_mono = 100
_n_cluster = 10
_n_cluster_mono = 10
_contactmap = np.ones((_n_mono,_n_mono),dtype=np.int32)
_monoconf = np.ones(_n_mono, dtype=np.int32)
_xcom = np.zeros((_n_mono,3),dtype=np.float32)


# length of first part of trj dropped 
t_cut = 5000000
t_mid = 10000000
t_end = 15500000
n_seg = 4

#output option
_to_showtime = False
_to_showdetail = True

_edge_all_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_edge_in_ring_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_edge_out_ring_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_edge_in_nodemain_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_edge_in_path_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_edge_between_mc_sc_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_edge_between_path_nm_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_edge_pair_path_allseg = np.zeros((n_seg,_n_cluster, _n_cluster), dtype=np.int32)
_cor_num_distr_by_type_allseg = np.zeros((n_seg,_n_cluster, n+1), dtype=np.int32)
_count_node_by_type_allseg = np.zeros((n_seg,_n_cluster), dtype=np.int32)
_cor_num_distr_by_monoconf_allseg = np.zeros((n_seg,_n_cluster_mono, n+1), dtype=np.int32)
_cnt_bondtype_by_monoconf_allseg = np.zeros((n_seg,_n_cluster_mono, _n_cluster), dtype=np.int32)
_conf_in_ring_allseg = np.zeros((n_seg,_n_cluster_mono), dtype=np.int32)
_conf_out_ring_allseg = np.zeros((n_seg,_n_cluster_mono), dtype=np.int32)
_conf_in_nodemain_allseg = np.zeros((n_seg,_n_cluster_mono), dtype=np.int32)
_conf_in_path_allseg = np.zeros((n_seg,_n_cluster_mono), dtype=np.int32)
_conf_in_sc_allseg = np.zeros((n_seg,_n_cluster_mono), dtype=np.int32)
_bond_type_by_conf_pair_allseg = np.zeros((n_seg,_n_cluster_mono, _n_cluster_mono, _n_cluster), dtype=np.int32)

_conf_in_path_ter_allseg = np.zeros((n_seg,_n_cluster_mono), dtype=np.int32)
_conf_pair_in_path_allseg = np.zeros((n_seg,_n_cluster_mono, _n_cluster_mono, _n_cluster), dtype=np.int32)
_conf_pair_between_mc_sc_allseg = np.zeros((n_seg,_n_cluster_mono, _n_cluster_mono, _n_cluster), dtype=np.int32)
_conf_pair_between_mc_nm_allseg = np.zeros((n_seg,_n_cluster_mono, _n_cluster_mono, _n_cluster), dtype=np.int32)


_count_allseg = np.zeros((n_seg,n+1), dtype=np.int32)
_tot_cycle_allseg = np.zeros((n_seg,n+1), dtype = np.int32)
_len_path_allseg = np.zeros((n_seg,n+1), dtype = np.int32)
_len_distr_allseg = np.zeros((n_seg,n+1, n+1), dtype = np.int32)
_n_branch_allseg = np.zeros((n_seg,n+1, n+1), dtype = np.int32)
_extra_node_main_allseg = np.zeros((n_seg,n+1, n+1), dtype = np.int32)
_ring_in_main_allseg = np.zeros((n_seg,n+1, n+1), dtype = np.int32)
_node_sc_allseg = np.zeros((n_seg,n+1, n+1), dtype = np.int32)

_node_main_distr_allseg = np.zeros((n_seg,n+1, n+1), dtype = np.int32)
_node_main_allseg = np.zeros((n_seg,n+1), dtype=np.int32)
_e_to_v_allseg = np.zeros((n_seg,n+1), dtype=np.float64)
_ring_size_allseg = np.zeros((n_seg,n+1,n+1), dtype = np.int32)
_e_to_v_det_allseg = np.zeros((n_seg,n+1,n+1), dtype=np.float64)

#polymer

_l_sq_allseg = np.zeros((n_seg, _n_mono), dtype = np.float64)
_r_sq_allseg =  np.zeros((n_seg,  _n_mono), dtype = np.float64)
_l_sq_cnt_allseg = np.zeros((n_seg, _n_mono), dtype = np.int64) 
_r_sq_cnt_allseg = np.zeros((n_seg,  _n_mono), dtype = np.int64)

# geom distribution 15 deg per bin
_angle_dih_distr_allseg = np.zeros((n_seg, 12, 24), dtype = np.int32)
_angle_dih_D_I_distr_allseg = np.zeros((n_seg, 15, 12, 24), dtype = np.int32)
_angle_dih_D_II_distr_allseg = np.zeros((n_seg, 15, 12, 24), dtype = np.int32)
_angle_dih_D_conf_I_distr_allseg = np.zeros((n_seg, 15, _n_cluster, 12, 24), dtype = np.int32)
_angle_dih_D_conf_II_distr_allseg = np.zeros((n_seg, 15, _n_cluster,12, 24), dtype = np.int32)

# 10 deg
_angle_distr_allseg = np.zeros((n_seg, 18), dtype = np.int32)
_dih_distr_allseg = np.zeros((n_seg, 36), dtype = np.int32)
_angle_ring_distr_allseg = np.zeros((n_seg, 7,18), dtype = np.int32)
_dih_ring_distr_allseg= np.zeros((n_seg, 7,36), dtype = np.int32)

_angle_D_distr_allseg = np.zeros((n_seg, 15, 18), dtype = np.int32)
_angle_D_conf_distr_allseg = np.zeros((n_seg, 15, _n_cluster,18), dtype = np.int32)

#coord type
_coord_typ_allseg = np.zeros((n_seg, 7), dtype = np.int32)
_coord_typ_conf_allseg = np.zeros((n_seg, _n_cluster, 7), dtype = np.int32)

_brch_pt_typ_allseg = np.zeros((n_seg, 6), dtype = np.int32)
_brch_pt_typ_conf_allseg = np.zeros((n_seg, _n_cluster, 6), dtype = np.int32)
_brch_cnt_allseg = np.zeros((n_seg,n+1), dtype=np.int32)
_brch_node_typ_allseg = np.zeros((n_seg, 6), dtype = np.int32)

#substract type
_substruct_typ_allseg = np.zeros((n_seg,11 , 11), dtype=np.int32)


i_seg = -1 
# -1

for fnm, fnm_cycle, fnm_path, fnm_contact, fnm_monoconf, fnm_xcom  in zip(fnm_list,\
				     fnm_cycle_list,\
				     fnm_path_list,\
				     fnm_contact_list,\
				     fnm_monoconf_list,\
				     fnm_xcom_list):
 _f =  open(fnm , 'r')
 _f_cycle = open(fnm_cycle, 'r')
 _f_path = open(fnm_path, 'r')
 _f_contact = open(fnm_contact,'r')
 _f_monoconf = open(fnm_monoconf, 'r')
 _f_xcom = open(fnm_xcom, 'r')



 print 'openning', fnm, fnm_cycle, fnm_path, fnm_contact, fnm_monoconf, fnm_xcom 
 _pass_mid = False
 i_seg+=1
 _count = _count_allseg[i_seg]
 _tot_cycle = _tot_cycle_allseg[i_seg]
 _len_path = _len_path_allseg[i_seg]
 _len_distr = _len_distr_allseg[i_seg]
 _n_branch = _n_branch_allseg[i_seg]
 _extra_node_main = _extra_node_main_allseg[i_seg]
 _node_main = _node_main_allseg[i_seg]
 _node_main_distr = _node_main_distr_allseg[i_seg]
 _e_to_v = _e_to_v_allseg[i_seg]
 _ring_size = _ring_size_allseg[i_seg]
 _e_to_v_det = _e_to_v_det_allseg[i_seg] 
 _ring_in_main = _ring_in_main_allseg[i_seg]
 _node_sc = _node_sc_allseg[i_seg]

 _edge_all = _edge_all_allseg[i_seg]
 _edge_in_ring = _edge_in_ring_allseg[i_seg]
 _edge_out_ring = _edge_out_ring_allseg[i_seg]
 _edge_in_nodemain = _edge_in_nodemain_allseg[i_seg]
 _edge_in_path = _edge_in_path_allseg[i_seg]
 _edge_between_mc_sc = _edge_between_mc_sc_allseg[i_seg]
 _edge_between_path_nm = _edge_between_path_nm_allseg[i_seg]
 _edge_pair_path = _edge_pair_path_allseg[i_seg]
 _cor_num_distr_by_type = _cor_num_distr_by_type_allseg[i_seg]
 _count_node_by_type = _count_node_by_type_allseg[i_seg]
 _cor_num_distr_by_monoconf = _cor_num_distr_by_monoconf_allseg[i_seg]
 _cnt_bondtype_by_monoconf = _cnt_bondtype_by_monoconf_allseg[i_seg]
 _conf_in_ring = _conf_in_ring_allseg[i_seg]
 _conf_out_ring = _conf_out_ring_allseg[i_seg]
 _conf_in_nodemain = _conf_in_nodemain_allseg[i_seg]
 _conf_in_path = _conf_in_path_allseg[i_seg]
 _conf_in_sc = _conf_in_sc_allseg[i_seg]
 _bond_type_by_conf_pair = _bond_type_by_conf_pair_allseg[i_seg]
 _conf_in_path_ter = _conf_in_path_ter_allseg[i_seg]
 _conf_pair_in_path = _conf_pair_in_path_allseg[i_seg]
 _conf_pair_between_mc_sc = _conf_pair_between_mc_sc_allseg[i_seg]
 _conf_pair_between_mc_nm = _conf_pair_between_mc_nm_allseg[i_seg]

 _angle_dih_distr = _angle_dih_distr_allseg[i_seg]
 _angle_distr = _angle_distr_allseg[i_seg]
 _dih_distr = _dih_distr_allseg[i_seg]
 _angle_ring_distr = _angle_ring_distr_allseg[i_seg]
 _dih_ring_distr = _dih_ring_distr_allseg[i_seg] 

 _angle_dih_D_I_distr = _angle_dih_D_I_distr_allseg[i_seg] 
 _angle_dih_D_II_distr = _angle_dih_D_II_distr_allseg[i_seg]
 _angle_dih_D_conf_I_distr = _angle_dih_D_conf_I_distr_allseg[i_seg]
 _angle_dih_D_conf_II_distr = _angle_dih_D_conf_II_distr_allseg[i_seg]
 _angle_D_distr = _angle_D_distr_allseg[i_seg]
 _angle_D_conf_distr = _angle_D_conf_distr_allseg[i_seg]

 _l_sq = _l_sq_allseg[i_seg]
 _r_sq = _r_sq_allseg[i_seg]
 _l_sq_cnt = _l_sq_cnt_allseg[i_seg]
 _r_sq_cnt = _r_sq_cnt_allseg[i_seg]


 _coord_typ = _coord_typ_allseg[i_seg]
 _coord_typ_conf = _coord_typ_conf_allseg[i_seg]

 _brch_pt_typ = _brch_pt_typ_allseg[i_seg]
 _brch_pt_typ_conf = _brch_pt_typ_conf_allseg[i_seg]
 _brch_cnt = _brch_cnt_allseg[i_seg]
 _brch_node_typ = _brch_node_typ_allseg[i_seg]


 _substruct_typ = _substruct_typ_allseg[i_seg]


 _show_time = 0
 _t , _txt = load_frame(_f)
 _t_cycle , _txt_cycle = load_frame(_f_cycle)
 _t_path , _txt_path = load_frame(_f_path)
 _t_contact, _txt_contact = load_frame(_f_contact)
 _t_monoconf, _txt_monoconf = load_frame(_f_monoconf)
 _t_xcom, _txt_xcom = load_frame(_f_xcom)

 if _t==_t_cycle and _t==_t_path and _t==_t_contact and  _t==_t_monoconf and _t==_t_xcom:
  pass
 else:
  print 'Error! mismatched time frames between',fnm,fnm_cycle,fnm_path, fnm_contact, fnm_monoconf,fnm_xcom
  print 'at',_t,_t_cycle, _t_path, _t_contact, _t_monoconf, _t_xcom
  sys.exit()

# _all_cycle = [[int(__j) for __j in (__i.replace('cycle:','')).split()] for __i in  _txt_cycle]


# print _t,_txt  
# print _t_path,_txt_path
# print _t_cycle, _txt_cycle
# print _t_contact, _txt_contact

# sys.exit()
 
 while _t >=0:
  _time = _t
  if _to_showdetail:
   print '=== time', _time,'==='

  if _time>_show_time and _to_showtime:
   print 'loading at t=', _show_time
   _show_time+=500000

  _t, _txt = load_frame(_f)


  _t_cycle , _txt_cycle = load_frame(_f_cycle)
  _t_path , _txt_path = load_frame(_f_path)
  _t_contact, _txt_contact = load_frame(_f_contact)
  _t_monoconf, _txt_monoconf = load_frame(_f_monoconf)
  _t_xcom, _txt_xcom = load_frame(_f_xcom)

  if _t==_t_cycle and _t==_t_path and _t==_t_contact and _t==_t_monoconf and _t==_t_xcom:
   pass
  else:
   print 'Error! mismatched time frames between',fnm,fnm_cycle,fnm_path, fnm_contact, fnm_monoconf, fnm_xcom
   print 'at',_t,_t_cycle, _t_path, _t_contact,_t_monoconf, _t_xcom
   sys.exit()

  if _time in exclu_list: continue
#  if not _time in test_list: continue
#  print _t,_txt  
#  print _t_path,_txt_path
#  print _t_cycle, _txt_cycle
#  print _t_contact, _txt_contact

#  sys.exit()




  _all_cycle = [[int(__j) for __j in (__i.replace('cycle:','')).split()] for __i in  _txt_cycle]

  _contactmap[:,:]=-1
  load_map(_contactmap, _txt_contact)
  _monoconf[:] = -1
  load_monoconf(_monoconf, _txt_monoconf)
  load_xcom(_xcom, _txt_xcom)
#  print _xcom

  if _txt ==[]: continue
  if _time< t_cut or _time> t_end: continue
  if (not _pass_mid) and _time> t_mid:
   i_seg+=1
   _count = _count_allseg[i_seg]
   _tot_cycle = _tot_cycle_allseg[i_seg]
   _len_path = _len_path_allseg[i_seg]
   _len_distr = _len_distr_allseg[i_seg]
   _n_branch = _n_branch_allseg[i_seg]
   _node_main = _node_main_allseg[i_seg]
   _node_main_distr = _node_main_distr_allseg[i_seg]
   _extra_node_main = _extra_node_main_allseg[i_seg]
   _e_to_v = _e_to_v_allseg[i_seg]
   _e_to_v_det = _e_to_v_det_allseg[i_seg]
   _ring_size = _ring_size_allseg[i_seg]
   _ring_in_main = _ring_in_main_allseg[i_seg]
   _node_sc = _node_sc_allseg[i_seg]

   _edge_all = _edge_all_allseg[i_seg] 
   _edge_in_ring = _edge_in_ring_allseg[i_seg] 
   _edge_out_ring = _edge_out_ring_allseg[i_seg] 
   _edge_in_nodemain = _edge_in_nodemain_allseg[i_seg] 
   _edge_in_path = _edge_in_path_allseg[i_seg] 
   _edge_between_mc_sc = _edge_between_mc_sc_allseg[i_seg] 
   _edge_between_path_nm = _edge_between_path_nm_allseg[i_seg]
   _edge_pair_path = _edge_pair_path_allseg[i_seg]
   _cor_num_distr_by_type = _cor_num_distr_by_type_allseg[i_seg] 
   _count_node_by_type = _count_node_by_type_allseg[i_seg] 
   _cor_num_distr_by_monoconf = _cor_num_distr_by_monoconf_allseg[i_seg]
   _cnt_bondtype_by_monoconf = _cnt_bondtype_by_monoconf_allseg[i_seg] 
   _conf_in_ring = _conf_in_ring_allseg[i_seg] 
   _conf_out_ring = _conf_out_ring_allseg[i_seg] 
   _conf_in_nodemain = _conf_in_nodemain_allseg[i_seg] 
   _conf_in_path = _conf_in_path_allseg[i_seg] 
   _conf_in_sc = _conf_in_sc_allseg[i_seg] 
   _bond_type_by_conf_pair = _bond_type_by_conf_pair_allseg[i_seg]
   _conf_in_path_ter = _conf_in_path_ter_allseg[i_seg] 
   _conf_pair_in_path = _conf_pair_in_path_allseg[i_seg] 
   _conf_pair_between_mc_sc = _conf_pair_between_mc_sc_allseg[i_seg] 
   _conf_pair_between_mc_nm = _conf_pair_between_mc_nm_allseg[i_seg]

   _angle_dih_distr = _angle_dih_distr_allseg[i_seg]
   _angle_distr = _angle_distr_allseg[i_seg]
   _dih_distr = _dih_distr_allseg[i_seg]
   _angle_ring_distr = _angle_ring_distr_allseg[i_seg]
   _dih_ring_distr = _dih_ring_distr_allseg[i_seg]    


   _angle_dih_D_I_distr = _angle_dih_D_I_distr_allseg[i_seg]       
   _angle_dih_D_II_distr = _angle_dih_D_II_distr_allseg[i_seg]
   _angle_dih_D_conf_I_distr = _angle_dih_D_conf_I_distr_allseg[i_seg]
   _angle_dih_D_conf_II_distr = _angle_dih_D_conf_II_distr_allseg[i_seg]
   _angle_D_distr = _angle_D_distr_allseg[i_seg]
   _angle_D_conf_distr = _angle_D_conf_distr_allseg[i_seg]

   _coord_typ = _coord_typ_allseg[i_seg]
   _coord_typ_conf = _coord_typ_conf_allseg[i_seg]

   _brch_pt_typ = _brch_pt_typ_allseg[i_seg]
   _brch_pt_typ_conf = _brch_pt_typ_conf_allseg[i_seg]
   _brch_cnt = _brch_cnt_allseg[i_seg]
   _brch_node_typ = _brch_node_typ_allseg[i_seg]

   _substruct_typ = _substruct_typ_allseg[i_seg]

   _l_sq = _l_sq_allseg[i_seg]
   _r_sq = _r_sq_allseg[i_seg]
   _l_sq_cnt = _l_sq_cnt_allseg[i_seg]
   _r_sq_cnt = _r_sq_cnt_allseg[i_seg]


   _pass_mid = True

  _i_path = 0 
  for _i in xrange(0,len(_txt),2):
   _o = oligo_top(mem = _txt[_i].replace('cluster member:',''), \
		 conn = _txt[_i+1], _time = _time)
#   print _txt_path[_i_path]
   _o.load_long_short_path([int(__t) for __t in _txt_path[_i_path].split()[3:]])
   _i_path+=1
   _o.load_cycle(_all_cycle)
 #doing staistics  

   _nt = len(_o.mem)
   if _nt>n: continue
   _count[_nt]+=1
   _tot_cycle[_nt]+= _o.cycle_cnt()
   _len_path[_nt] += len(_o.path)
   _len_distr[_nt, len(_o.path)]+=1
   _nm = _o.node_in_main()
#   _n_brch = _o.n_branches()
#   _n_branch[_nt, len(_o.path)] += _n_brch
#   _n_branch[_nt, _nm] += _n_brch
   _e_to_v_det[_nt, len(_o.path)]+= float(_nm)/float(len(_o.path)) - 1.
   _extra_node_main[_nt, len(_o.path)]+=(_nm - len(_o.path))
   _ring_in_main[_nt, len(_o.path)]+= _o.node_in_main(countRing=True)
   _node_sc[_nt, len(_o.path)]+= (_nt - _nm)
   _node_main[_nt]+= _nm
   _node_main_distr[_nt, _nm]+= 1
   _e_to_v[_nt]+= _o.e_to_v()
   if _nt >=4 :
    pass
    _o.cycle_size(_ring_size[_nt], countOnce=True)
# substructure element analysis
#    _o.substruct_typ( _substruct_typ)

# branch structure analysis

#    _old_cnt = _brch_node_typ[0]


#    _brch_cnt[_nt]+=_o.substruct_typ(_substruct_typ, anaBranch=True,\
#                   _brch_pt_typ=_brch_pt_typ, _monoconf=_monoconf,\
#                   _brch_pt_typ_conf=_brch_pt_typ_conf, printIdx=False,\
#		   countNode=True, _brch_node_typ=_brch_node_typ)

#    _old_cnt= _brch_node_typ[0]- _old_cnt

#    _old_cnt_a = _coord_typ[1]

#    _o.coord_typ(_monoconf, _coord_typ, _coord_typ_conf, printIdx=False)

#    _old_cnt_a = _coord_typ[1] - _old_cnt_a
#    if _old_cnt > _old_cnt_a:
#     print 'Warning: branch typ 0 cnt greater than out-of-ring multivalent',\
#		_old_cnt, _old_cnt_a


#     sys.exit()
## coupling between mono conf and bond type

   # we can adjust range of oligomer order for the analysis
#main path
#   if len(_o.path)>3:
#    pass
#    _o.mainpath_geom(_xcom, _l_sq, _r_sq, _l_sq_cnt, _r_sq_cnt)


#   if _nt>=4 and _nt<=10:
#   if _nt > 10:
#   if _nt>10:  
   if _nt>=4:
    pass
   #coordination geometry analysis
#    _o.coord_mainpath_geom(_xcom, _angle_distr, _dih_distr)
#    _o.coord_ring_geom(_xcom, _angle_ring_distr, _dih_ring_distr)

#    _o.coord_D_geom(_xcom, _monoconf,\
#		    _angle_dih_D_I_distr, _angle_dih_D_II_distr,\
#		    _angle_dih_D_conf_I_distr, _angle_dih_D_conf_II_distr,\
#		    _angle_D_distr, _angle_D_conf_distr)



   # bond type between conf types
#    _o.bond_type_between_conf_pair(_contactmap, _monoconf, \
#				   _bond_type_by_conf_pair)

#    _o.monoconf_edge_couple(_contactmap, _monoconf, _cor_num_distr_by_monoconf,\
#			    _cnt_bondtype_by_monoconf)

##various bond type analysis
##conf analysis
#    _o.edge_conf_type_ana(_contactmap,_edge_all, _edge_in_ring, _edge_out_ring,\
#		    _edge_in_nodemain, _edge_in_path,\
#		    _edge_between_mc_sc, _edge_between_path_nm,\
#		     _monoconf, _conf_in_ring, _conf_out_ring,\
#		     _conf_in_nodemain, _conf_in_path, _conf_in_sc,\
#		     _conf_in_path_ter,\
#		     _conf_pair_in_path,_conf_pair_between_mc_sc,_conf_pair_between_mc_nm,\
#		     showBond = _to_showdetail)

#   _o.edge_corr_path(_contactmap,_edge_pair_path)

#   if _nt>3:
#    _o.bond_type_multivalency(_contactmap, _cor_num_distr_by_type,_count_node_by_type)


 
 _f.close()
 _f_cycle.close()
 _f_path.close()
 _f_contact.close()
 _f_monoconf.close()

print '===overall analysis==='

for _i_seg in xrange(4):
 print 'segment', _i_seg


# analyze substructure to prove linear topology
# if _substruct_typ_allseg[_i_seg,1,:].sum()>0:
#
#  
#  print "=== linear substructure  num of conn  frequency"
#  for _i in xrange(1,11):
#   print _i, _substruct_typ_allseg[_i_seg,1,_i]
#
#  print '==end'
# #ring substruct
# if _substruct_typ_allseg[_i_seg, 3:,:].sum()>0:
#  print "=== ring substructure ring_sez num_conn frequency"
#  for _i in xrange(3,11):
#   for _j in xrange(1,11):
#    print _i,_j, _substruct_typ_allseg[_i_seg,_i,_j]
#  print '=== end'





# coord geom
# if _angle_distr_allseg[_i_seg].sum()>0:
#  print 'mainpath ang distr:' 
#  for _i in xrange(18):
#    print (_i+0.5)*10.0,\
#     float(_angle_distr_allseg[_i_seg, _i])/float(_angle_distr_allseg[_i_seg].sum())
#  print '==='
#
# if _dih_distr_allseg[_i_seg].sum()>0:
#  print 'mainpath dih distr:'
#  for _i in xrange(36):
#    print (_i+0.5)*10.0-180.0,\
#    float( _dih_distr_allseg[_i_seg, _i])/float(_dih_distr_allseg[_i_seg].sum())
#
#  print '==='

# corrdingation in ring
# for _rsize in xrange(3,7):
#  if _angle_ring_distr_allseg[_i_seg, _rsize].sum()>0:
#   print 'ring',_rsize,'ang distr:'
#   for _i in xrange(18):
#     print (_i+0.5)*10.0, _angle_ring_distr_allseg[_i_seg,_rsize, _i]
#   print '==='
#
#  if _dih_ring_distr_allseg[_i_seg, _rsize].sum()>0:
#   print 'ring', _rsize, 'dih distr:'
#   for _i in xrange(36):
#     print (_i+0.5)*10.0-180.0, _dih_ring_distr_allseg[_i_seg,_rsize, _i]
#
#   print '==='

### main path polymer property


# for _i in xrange(_n_mono):
#  if _l_sq_cnt_allseg[_i_seg,_i]>0 and _r_sq_cnt_allseg[_i_seg, _i]>0:
#   print 'len=', _i,
#   _r_2 = _r_sq_allseg[_i_seg, _i]/float(_r_sq_cnt_allseg[_i_seg,_i])
#   _l_2 = _l_sq_allseg[_i_seg, _i]/float(_l_sq_cnt_allseg[_i_seg, _i])
#
#   print '<r^2>=',_r_2,'<l^2>=', _l_2, 'C=', _r_2/_l_2/float(_i-1)



##coord type
# _s = _coord_typ_allseg[_i_seg].sum()
# if _s>0:
#  print '===coord by type'
#  for _ci in xrange(7):
#   print 'type',_ci,_coord_typ_allseg[_i_seg, _ci],\
#		float(_coord_typ_allseg[_i_seg, _ci])/float(_s)
#  for _ci in xrange(7):
#   
#    _cs = _coord_typ_conf_allseg[_i_seg, :, _ci].sum()
#    if _cs>0:
#     print '===coord by conf for type', _ci
#
#     for __i in xrange(_n_cluster):
#      print 'conf',__i, _coord_typ_conf_allseg[_i_seg, __i, _ci],\
#		float(_coord_typ_conf_allseg[_i_seg, __i, _ci])/float(_cs)


##branch typ
# _s = _brch_pt_typ_allseg[_i_seg].sum()
# if _s>0:
#  print '===branch by type'
#  for _ci in xrange(6):
#   print 'type',_ci,_brch_pt_typ_allseg[_i_seg, _ci],\
#               float(_brch_pt_typ_allseg[_i_seg, _ci])/float(_s),\
#	       'node cnt', _brch_node_typ_allseg[_i_seg, _ci]
#  for _ci in xrange(6):
#   
#    _cs = _brch_pt_typ_conf_allseg[_i_seg, :, _ci].sum()
#    if _cs>0:
#     print '===branch by conf for type', _ci
#
#     for __i in xrange(_n_cluster):
#      print 'conf',__i, _brch_pt_typ_conf_allseg[_i_seg, __i, _ci],\
#               float(_brch_pt_typ_conf_allseg[_i_seg, __i, _ci])/float(_cs)
# print '===ave branch count by oligomer size'
# for _in in xrange(3,n+1):
#  if _count_allseg[_i_seg,_in]>0:
#   print _in, float(_brch_cnt_allseg[_i_seg,_in])/float(_count_allseg[_i_seg,_in])
# print '===end branch analysis'
#****



###coord geom by D
#***

# for _D in [2,3,4]:
#  if _angle_D_distr_allseg[_i_seg,_D].sum()>0:
#   print '===angle b-a-c by D =', _D
#   for _i in xrange(18):
#     print (_i+0.5)*10.0, _angle_D_distr_allseg[_i_seg,_D, _i]
#  for _c in xrange(_n_cluster):
#   if _angle_D_conf_distr_allseg[_i_seg,_D,_c].sum()>0:
#    print '===angle b-a-c by D =', _D,'conf =',_c
#    for _i in xrange(18):
#      print (_i+0.5)*10.0, _angle_D_conf_distr_allseg[_i_seg,_D,_c, _i]
#
#  if _angle_dih_D_I_distr_allseg[_i_seg,_D].sum()>0:
#   print '$$$angle and dih b-a-d by D =', _D
#   for _i in xrange(12):
#    for _j in xrange(24):
#     print (_i+0.5)*15.0,(_j+0.5)*15.0-180.0, _angle_dih_D_I_distr_allseg[_i_seg,_D, _i,_j]
#   print '$$$end'
#
#  for _c in xrange(_n_cluster):
#   if _angle_dih_D_conf_I_distr_allseg[_i_seg,_D,_c].sum()>0:
#    print '$$$angle and dih b-a-d by D =', _D,'conf = ',_c
#    for _i in xrange(12):
#     for _j in xrange(24):
#      print (_i+0.5)*15.0,(_j+0.5)*15.0-180.0, _angle_dih_D_conf_I_distr_allseg[_i_seg,_D,_c, _i,_j]
#    print '$$$end'
#
#  if _angle_dih_D_II_distr_allseg[_i_seg,_D].sum()>0:
#   print '$$$angle and dih b-a-e by D =', _D
#   for _i in xrange(12):
#    for _j in xrange(24):
#     print (_i+0.5)*15.0,(_j+0.5)*15.0-180.0, _angle_dih_D_II_distr_allseg[_i_seg,_D, _i,_j]
#   print '$$$end'
#
#  for _c in xrange(_n_cluster):
#   if _angle_dih_D_conf_II_distr_allseg[_i_seg,_D,_c].sum()>0:
#    print '$$$angle and dih b-a-e by D =', _D,'conf = ',_c
#    for _i in xrange(12):
#     for _j in xrange(24):
#      print (_i+0.5)*15.0,(_j+0.5)*15.0-180.0, _angle_dih_D_conf_II_distr_allseg[_i_seg,_D,_c, _i,_j]
#    print '$$$end'
#*****end


#cycle analysis
# cyc_cnt_ana(_count_allseg, _ring_size_allseg, _i_seg, 4,10,n)
# cyc_cnt_ana(_count_allseg, _ring_size_allseg, _i_seg, 11,n,n)
 cyc_cnt_ana(_count_allseg, _ring_size_allseg, _i_seg, 4,n,n)



 for _i in xrange(3,n+1):
  if _count_allseg[_i_seg,_i]>0:
#  print _i,float(_tot_cycle[_i])/float(_count[_i]),\
#	float(_len_path[_i])/float(_count[_i]),
#  print float(_node_main[_i])/float(_count[_i]),\
#	_e_to_v[_i]/float(_count[_i])

   print _i ,  _count_allseg[_i_seg,_i],\
     'ave_len',float(_len_path_allseg[_i_seg,_i])/float(_count_allseg[_i_seg,_i]),
   for _j in xrange(2,n+1):
    pass
#    print float(_ring_size_allseg[_i_seg, _i,_j])/float(_count_allseg[_i_seg, _i]),
#    print float(_len_distr_allseg[_i_seg,_i,_j])/float(_count_allseg[_i_seg, _i]),
#    print float(_node_main_distr_allseg[_i_seg,_i,_j])/float(_count_allseg[_i_seg, _i]),

   

#for branch count for a given chain len of a given oligomer size
#    if _len_distr_allseg[_i_seg, _i, _j]>0:
#     print float(_n_branch_allseg[_i_seg, _i, _j])/float(_len_distr_allseg[_i_seg, _i, _j]),
#    else:
#     print '0.0',

#for addtional monomer in mainchain for a given chain len of a given oligomer size
#    if _len_distr_allseg[_i_seg, _i, _j]>0:
#     print float(_extra_node_main_allseg[_i_seg, _i, _j])/float(_len_distr_allseg[_i_seg, _i, _j]),
#    else:
#     print '0.0',

# ring count in main chain for a given chain length and oligomer size
#    if _len_distr_allseg[_i_seg, _i, _j]>0:
#     print float(_ring_in_main_allseg[_i_seg, _i, _j])/float(_len_distr_allseg[_i_seg, _i, _j]),
#    else:
#     print '0.0',

# node in sc for a given chain length and a given oligomer size
#    if _len_distr_allseg[_i_seg, _i, _j]>0:
#     print float(_node_sc_allseg[_i_seg, _i, _j])/float(_len_distr_allseg[_i_seg, _i, _j]),
#    else:
#     print '0.0',


#for branch count for a given main chain noce count of a given oligomer size
#    if _node_main_distr_allseg[_i_seg, _i, _j]>0:
#     print float(_n_branch_allseg[_i_seg, _i, _j])/float(_node_main_distr_allseg[_i_seg, _i, _j]),
#    else:
#     print '0.0',

# e to v
#    if _len_distr_allseg[_i_seg, _i, _j]>0:
#     print float(_e_to_v_det_allseg[_i_seg, _i, _j])/float(_len_distr_allseg[_i_seg, _i, _j]),
#    else:
#     print '-0.25',



   print

# ouptut the results of bond type analysis
print_list_bond = [['all bond', _edge_all_allseg],['bond in ring',_edge_in_ring_allseg ],\
	      ['bond out ring', _edge_out_ring_allseg],['bond in path', _edge_in_path_allseg],\
	      ['bond in nodemain', _edge_in_nodemain_allseg],['bond between mc and sc',_edge_between_mc_sc_allseg],\
	      ['bond between path and nm', _edge_between_path_nm_allseg],\
	      ['conf in ring', _conf_in_ring_allseg],\
	      ['conf out ring', _conf_out_ring_allseg],\
	      ['conf in path', _conf_in_path_allseg],\
	      ['conf in sc', _conf_in_sc_allseg],\
	      ['conf in nodemain', _conf_in_nodemain_allseg],\
	      ['conf in path terminus', _conf_in_path_ter_allseg]]


print '==========***=========='
print 'bond/conf type analysis'
print '======================='
for _ii in xrange(n_seg):
 print 'segment', _ii

# print all types of statiscal results from print_list
# for _j in print_list_bond:
#  print_distr(_j[0], _j[1][_ii], normalize=False)

#for bond type
# for _jj in xrange(_n_cluster):
#  _s = _count_node_by_type_allseg[_ii, _jj]
#  if _s<=0: continue
#  print 'bond type',_jj,":",
#  for _kk in _cor_num_distr_by_type_allseg[_ii, _jj]:
#   print float(_kk)/float(_s),
#  print

# print_distr('conf pair in path', _conf_pair_in_path_allseg[_ii], hidezero = True,\
#	      print3d = True)

# print_distr('conf pair between mc and sc',\
#	      _conf_pair_between_mc_sc_allseg[_ii],  hidezero = True,\
#	      print3d=True)

# print_distr('conf pair between path and nm',\
#              _conf_pair_between_mc_nm_allseg[_ii],  hidezero = True,\
#              print3d=True)


# for mono conf type
# for _jj in xrange(_n_cluster_mono):
#  _s = _cor_num_distr_by_monoconf_allseg[_ii, _jj].sum()
#  if _s<=0: continue
#  print 'conf type',_jj,'tot frame',_s ,": coordination number distribution (from 0)",
#  for _kk in _cor_num_distr_by_monoconf_allseg[_ii, _jj]:
#   print float(_kk)/float(_s),
#  print

# for mono conf bond coupling
# for _jj in xrange(_n_cluster_mono):
#  _s = _cnt_bondtype_by_monoconf_allseg[_ii, _jj].sum()
#  if _s<=0: continue
#  print 'conf type',_jj,": involved bond type distribution (from type 0)",
#  for _kk in _cnt_bondtype_by_monoconf_allseg[_ii, _jj]:
#   print float(_kk)/float(_s),
#  print

# bond type by conf pair
# for _iii in xrange(_n_cluster_mono):
#  for _jjj in xrange(_n_cluster_mono):
#   print_distr("conf type %d %d"%(_iii, _jjj),\
#	       _bond_type_by_conf_pair_allseg[_ii, _iii, _jjj],\
#	       singleline=True, showsum=True)

# output results of edge correlation in main path

# if _edge_pair_path_allseg[_ii].sum()>0:
#  print 'edge pair correlation in path'
#  _s = _edge_pair_path_allseg[_ii].sum()
#  for _i in xrange(_n_cluster):
#   for _j in xrange(_n_cluster):
#    if _edge_pair_path_allseg[_ii, _i, _j]>0:
#     print _i , _j, float(_edge_pair_path_allseg[_ii, _i, _j])/float(_s)

print '===========***==========='
print 'end of bond/conf analysis'
print '========================='


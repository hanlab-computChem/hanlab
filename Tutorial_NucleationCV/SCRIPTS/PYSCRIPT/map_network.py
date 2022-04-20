import sys
import numpy as np
import networkutil
import math

# base time 0.1 us
# base conf time
_tc = 0.45 
# base diff time at 1mM
_td = 1.08


def add_nb(_c, _idx, _n,_i,_j):
 if _n<0 or _n>18 or _i<0 or _i>17 or _j<0 or _j>17: return

 if _idx[_n,_i,_j]>=0:
  _c.append(_idx[_n,_i,_j])


def _get_v(_i, _j, _node, _data, _conc):

 _dn = _node[_j][0] - _node[_i][0]
 _di = _node[_j][1] - _node[_i][1]
 _dj = _node[_j][2] - _node[_i][2]

 if _dn==0 and _di==0 and _dj==0:
  print 'Error: transition between',_node[_i],'and',_node[_j],'seems not right'
  sys.exit()

 lnPi = -_data[_node[_i][0],_node[_i][1],_node[_i][2]]
 lnPj = -_data[_node[_j][0],_node[_j][1],_node[_j][2]]


 if _dn==0:
#conf change
  return math.log(1./_tc)

 elif _dn>0:
  if lnPi>lnPj:
   return math.log(1./_td*_conc) + lnPi-lnPj
  else:
   return math.log(1./_td*_conc)
 else:
  if lnPi>lnPj:
   return math.log(1./_td*_conc)
  else:
   return math.log(1./_td*_conc) + lnPj - lnPi




# states with energy above this cutoff will not be considered
_ene_cut = 50.0


_data = np.ones((19,18,18),dtype = np.float64)*200.0
_idx = np.ones((19,18,18),dtype = np.int32)*(-1)

_s_id = 0
_all_cnt =0

_conc = float(sys.argv[2])

print 'Info: loading file', sys.argv[1]
print 'Info: Concentration', _conc,'mM'

for _rl in open(sys.argv[1], 'r'):

 srl = _rl[:-1].split()

 _all_cnt +=1
 if float(srl[3])< _ene_cut:
  

  _data[int(srl[0]), int(float(srl[1])+0.01), int(srl[2])] = float(srl[3])
  if _idx[int(srl[0]), int(float(srl[1])+0.01), int(srl[2])]!=-1:
   print 'Error: repeated data point found'
   sys.exit()

  _idx[int(srl[0]), int(float(srl[1])+0.01), int(srl[2])] = _s_id
  _s_id+=1


print 'Info: ',_all_cnt,'data points loaded and',_s_id,'data points used'


_node=[]
_conn=[]

__id = 0
for _n in xrange(19):
 for _i in xrange(18):
  for _j in xrange(18):

   if _idx[_n,_i,_j]>=0:

    if _idx[_n, _i, _j]!=__id:
      print 'Error: idx ordering is incorrect.'
      sys.exit()

    _node.append([_n, _i, _j])
    _c=[]
    add_nb(_c, _idx, _n-1,_i,_j)
    add_nb(_c, _idx, _n+1,_i,_j)
    add_nb(_c, _idx, _n,_i-1,_j)
    add_nb(_c, _idx, _n,_i+1,_j)
    add_nb(_c, _idx, _n,_i,_j-1)
    add_nb(_c, _idx, _n,_i,_j+1)

    _conn.append(_c)
    __id+=1

_re = networkutil.strong_connected_components(_conn)

if len(_re[0])<_s_id:
 print 'Info: largest connected componet contain',len(_re[0]),'data points'
 print 'Info: the unconnected data will be discarded'

 _idx[:,:,:]=-1

 for _j,_i in enumerate(_re[0]):
  _idx[_node[_i][0], _node[_i][1], _node[_i][2]] = _j

 _node=[]
 _conn=[]

 __id = 0
 for _n in xrange(19):
  for _i in xrange(18):
   for _j in xrange(18):

    if _idx[_n,_i,_j]>=0:

     if _idx[_n, _i, _j]!=__id:
       print 'Error: idx ordering is incorrect.'
       sys.exit()

     _node.append([_n, _i, _j])
     _c=[]
     add_nb(_c, _idx, _n-1,_i,_j)
     add_nb(_c, _idx, _n+1,_i,_j)
     add_nb(_c, _idx, _n,_i-1,_j)
     add_nb(_c, _idx, _n,_i+1,_j)
     add_nb(_c, _idx, _n,_i,_j-1)
     add_nb(_c, _idx, _n,_i,_j+1)

     _conn.append(_c)
     __id+=1

 _re = networkutil.strong_connected_components(_conn)

 print 'Info: now network contains',len(_re[0]),'data points'





# network file used for anaNetwork_t -in
f= open('network.xvg','w')
# base rate file used for anaNetwork_t -baseRate
f_br= open('br.xvg','w')
# coord of data point
f_crd = open('crd.xvg','w')
# FE
f_fe = open('fe.xvg','w')

for _in in xrange(len(_conn)):


 _s1 = '%d : '%_in
 for _j in _conn[_in]:
  _s1+=' %d'%_j

  _s2 = '%d %d '%(_in, _j)

  _v = _get_v(_in, _j, _node, _data, _conc)

  f_br.write('%s %f  [%d %d %d] fe %f ==>[%d %d %d] fe %f\n'%(_s2, _v,\
   _node[_in][0],_node[_in][1],_node[_in][2], _data[_node[_in][0],_node[_in][1],_node[_in][2]],
   _node[_j][0],_node[_j][1],_node[_j][2], _data[_node[_j][0],_node[_j][1],_node[_j][2]]))
  

 f_fe.write('%d : %f\n'%(_in, _data[_node[_in][0],_node[_in][1],_node[_in][2]]))
 f.write('%s\n'%_s1)

 f_crd.write('%d : %d %d %d\n'%(_in, \
   _node[_in][0],_node[_in][1],_node[_in][2]))


f_fe.close()
f.close()
f_br.close()
f_crd.close()

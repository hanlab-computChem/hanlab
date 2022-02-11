import numpy
import sys

import scipy
import scipy.sparse.linalg
import scipy.optimize

def print_max(rk):
 print 'Info : Max residue ', rk

#binary search for an element in an array
#return its index, return -1 if not found
def get_index(v, _list):
 if _list==[]: return -1
 l = 0
 r = len(_list)-1
 m = (l+r)/2

 if v==_list[m]: return m

 while l<=r:
  m = (l+r)/2
  if v==_list[m]: return m
  if v<_list[m]: r = m-1
  if v>_list[m]: l=m+1

 return -1


def load_transition_matrix (fnm, fe, abs_2_id, dt = 0.01, NotCRS=False, RateMatrix=False):


 if RateMatrix:
  print 'Info: Since outputing rate matrix dt will be always 1.0'
  print 'Info: Results will not be outputed in a sparse matrix.'
  dt = 1.0
  NotCRS = True

  
 
 f=open(fnm,'r')
# all id are relative ids
 conn=[]
 k=0
 for rl in f:
  srl = rl[:-1].split()
  if len(srl)==0: continue
  conn.append(sort_unq([abs_2_id[int(i)] for i in srl[2:]]+[k]))
  k+=1
 f.close()

 if len(conn)!=len(fe):
  print 'Error : connectivity data dont match free energy data'
  sys.exit()

# gen transition matrix

 T=[]
 for i,ic in enumerate(conn):
  _idx = -1
  _sum = 0.0
  _t=[]
  for k,j in enumerate(ic):
   if j==i: 
    _idx = k
    _t.append(0.0)
    continue

   # i-->j transition
   fe_i = fe[i]
   fe_j = fe[j]
   if fe_i>fe_j:
    _t_ij = dt
   else:
    _t_ij = numpy.exp(fe_i - fe_j)*dt
   
   _sum+=_t_ij
  
   _t.append(_t_ij)
  
  if _idx<0:
   print 'Error: elements ii for',i,'not found'
   sys.exit()

  if _sum>=1.0 and not RateMatrix:
   print 'Error: elements ii',i,' is non-positive, try smaller dt'
   sys.exit()

  if RateMatrix:
   _t[_idx] = - _sum
  else:  
   _t[_idx] = 1.0- _sum

  T.append(_t)
#i,j relative cluster ids
# Tij  transition probability of i->j
# p(t)^T . T(dt) =  p(t+dt)^T


 _iptr = 0
 _data = []
 _ind = []
 _indptr=[0]

 for ic, it in zip(conn, T):
  if len(ic)!=len(it):
   print 'Error : unmatched row when converting transition matrix'
   sys.exit()

  _data+=it
  _ind+=ic
  _iptr+=len(ic)
  _indptr.append(_iptr)

 #print _data
 #print _ind
 #print _indptr

 if NotCRS: return (conn, T)


 return scipy.sparse.csr_matrix(\
(numpy.array(_data, numpy.float64),\
 numpy.array(_ind, numpy.int32),\
 numpy.array( _indptr, numpy.int32))\
 , shape=(len(conn), len(conn)) )


def construct_sparse_from_conn (conn, T):
 _iptr = 0
 _data = []
 _ind = []
 _indptr=[0]

 for ic, it in zip(conn, T):
  if len(ic)!=len(it):
   print 'Error : unmatched row when converting transition matrix'
   sys.exit()

  _data+=it
  _ind+=ic
  _iptr+=len(ic)
  _indptr.append(_iptr)

 #print _data
 #print _ind
 #print _indptr

 # sparse matrix accounts for i->j transition

 return scipy.sparse.csr_matrix(\
(numpy.array(_data, numpy.float64),\
 numpy.array(_ind, numpy.int32),\
 numpy.array( _indptr, numpy.int32))\
 , shape=(len(conn), len(conn)) )



def get_dist(eigv, v):
 """
 return array of dist between list of vector and a vector
 """
 

 _sum = numpy.zeros([len(eigv)],numpy.float64)
 for i in range(len(v)):
  _dx = eigv[:,i] - v[i]
  _sum = _sum + _dx * _dx

 return numpy.sqrt(_sum + 1e-35)

def f_ax_b(x, m, s):
 v = numpy.dot(m, x) -s
 return numpy.dot(v,v)

def pcca(eigv, usedist=False):
#  print eigv
  n, m = eigv.shape
  m+=1
  print 'Info : PCCA clustering ',n,'microstates into ',m,'macrostates'
  if m<2:
   print 'error : at least 2 macrostates are needed'
   sys.exit()

  _w = numpy.zeros([n,m-1], numpy.float32)
   
  _sum = numpy.zeros([n],numpy.float64)
# find most ditant 2
  print 'Info : look for two most distant points...'
  sys.stdout.flush()
  maxdist = -0.1
  maxpair = (-1,-1)
  for k in range(n):
   _sum[:]=0.0
   _p_c = eigv[k,:]
  # print _p_c
   for i in range(m-1):
    _dx =  eigv[:, i] - _p_c[i]  
    _sum= _sum + _dx*_dx 
  #  print _dx
  #  print _sum
   _mi = numpy.argmax(_sum) 
  # print k, _mi, _sum[_mi]
   if _sum[_mi]>maxdist:
    maxdist = _sum[_mi]
    maxpair = (k, _mi)

   if k%1000==0:
    print 'Info :',k,'nodes exmained'
    sys.stdout.flush()

  mark=numpy.zeros([n], numpy.bool)
# mark states as macrostate
  mark[maxpair[0]] =True 
  mark[maxpair[1]] =True 

  tomacrostat =get_dist(eigv, eigv[maxpair[0]])+\
                get_dist(eigv, eigv[maxpair[1]]) 

 
  #tolist=[maxpair[0], maxpair[1]]

  print 'Info : Looking for more nodes..'

  for i in range(2,m,1):
   _which  = numpy.where(mark==False)  
   _mi = numpy.argmax(tomacrostat[_which])    
   _midx = _which[0][_mi]
   mark[_midx]=True
   tomacrostat = tomacrostat + get_dist(eigv, eigv[_midx]) 

#now m points are located
  
  _macro = numpy.where(mark)
#  group= [[i] for i in _macro[0]]
  group = [[] for i in _macro[0]]

  print len(_macro[0]),'macrostates found'
  if len(_macro[0])!=m:
   print 'error : unexpected macrostate count'
   sys.exit()


  if usedist:
   print 'Info : use minimum dist to determine the memberships'
   _macrocen = eigv[_macro[0],:]
   for i in range(n):
    if not mark[i]:
     _dist = get_dist(_macrocen, eigv[i,:])
     _belongto = numpy.argmin(_dist)
     group[_belongto].append(i)
   return group


  print 'Info : construct simplex distance'
  print 'Info : construct macrostate matrix with dimension',m,'by',m


  # [ M1-Mm M2-Mm M3-Mm ... Mm-1 - Mm ] phi = M-Mm
# [M1 .. Mm] phi =M
#0<= M1<=1. Sum_(Mm) = 1


  _mmatr = numpy.zeros([m-1,m], numpy.float64)
  _phival = numpy.ones([m], numpy.float64)*1.0/m
  _tosolve = numpy.zeros([m-1,n], numpy.float64)
#  _Mm = eigv[_macro[0][m-1], :]

  for i in range(m):
   _mmatr[:, i] = eigv[_macro[0][i],:] 

#  print _mmatr

  # get phi for each micro state by solving the equations
# slow
#  for i in range(n):
#   if not mark[i]:
#    _phi = numpy.linalg.solve(_mmatr, eigv[i, :] - _Mm) 
#    _phival[:m-1] = _phi[:] 
#    _phival[m-1] = 1.0 - _phi.sum()

# fill solver matrix
  for i in range(n):
   _tosolve[:, i] = eigv[i, :] 

#  _phi = numpy.linalg.solve(_mmatr, _tosolve)

# solve the linear eq with restraint

#  res = scipy.optimize.minimize(\
#	lambda x: numpy.dot(numpy.dot(_mmatr, x)- _tosolve), \
#	numpy.dot(_mmatr, x)- _tosolve), _

  for i in range(n):
#   if not mark[i]:
#     _phival[:m] = _phi[:, i] 


    cons = ({'type':'eq','fun':lambda x:x.sum()-1.0})
    res = scipy.optimize.minimize(\
        f_ax_b, _phival,args=(_mmatr, _tosolve[:,i]),\
	 method ='SLSQP', bounds=[(0.0,1.0)]*m,\
	constraints = cons, options = {'maxiter':10000})
#numpy.dot(numpy.dot(_mmatr, x)- _tosolve[:,i]), \
#        numpy.dot(_mmatr, x)- _tosolve[:,i]), [1.0/m]*m,\
#
#	method='SLSQP', constraints = cons,\
#	options={'disp': False})
#	bounds=(0.0,1.0))

    _phival[:]=res['x']

#    _phival[m-1] = 1.0 - _phi[:, i].sum()





#  p1*M1 + p2*M2 + ... +pm*Mm = M    sum(p) = 1
# i belongs to j if |pj| is max
#    print _phival 
    _belongto = numpy.argmax(numpy.absolute(_phival)) 
#    print _belongto
    group[_belongto].append(i)
    #print _phival

    if i%100==0:
     print 'Info :',i,'states exmained'
     sys.stdout.flush()

  
  return group
 
def logSum(log_terms):
   """Compute the log of a sum of terms whose logarithms are provided.

   REQUIRED ARGUMENTS
      log_terms is the array (possibly multidimensional) containing the logs of the terms to be summed.

   RETURN VALUES
      log_sum is the log of the sum of the terms.

   """

   # compute the maximum argument
   max_log_term = log_terms.max()

   # compute the reduced terms
   terms = numpy.exp(log_terms - max_log_term)

   # compute the log sum
   log_sum = numpy.log( terms.sum() ) + max_log_term

   # return the log sum
   return log_sum


def load_list(fnm):
 f=open(fnm, 'r')
 _bind_idx=[]
 for rl in f:
  _bind_idx+=[int(i) for i in rl[:-1].split()]

 f.close()

 return _bind_idx

def write_list (fnm, _l):
 f=open(fnm, 'w')
 for k,i in enumerate(_l):
  f.write('%d '%i)
  if (k+1)%8==0:
   f.write('\n')

 f.write('\n')
 f.close()


def remove_idx(nonfib_idx, rm_list):
 rm_list.sort()

 _rmidx = rm_list
 j = nonfib_idx

 k=0
 l=0
 _new_c = []

 while k<len(_rmidx) and l<len(j):
  if _rmidx[k] < j[l]:
   k+=1
  elif _rmidx[k] > j[l]:
   _new_c.append(j[l])
   l+=1
  else:
   k+=1
   l+=1

 while l<len(j):
  _new_c.append(j[l])
  l+=1

 return _new_c


def sort_unq(_conn):
 if _conn==[]: return []
 _tt = _conn[:]
# print ii, _tt
 _tt.sort()
 _tc=[_tt[0]]
 for j in range(1,len(_tt),1):
  if _tt[j-1]!=_tt[j]:
   _tc.append(_tt[j])

 return _tc


def load_rate(fnm):


 rate = None
 f=open(fnm,'r')
 for rl in f:
  if 'rate' in rl:
   rate = float(rl[:-1].split()[3])

 f.close()

 if rate is None:
  print 'Wrong rate file', fnm
  sys.exit()

 return rate

def load_clu_prop_single(fnm, col=2):
 f= open(fnm,'r')
 _cluidx=[]
 _data = []



 for rl in f:
  if rl[0] in ['#','@']: continue
  if 'Info' in rl: continue
  srl=rl[:-1].split()
  if ']' in rl: continue 
  _cluidx.append(int(srl[0]))
  _data.append(float(srl[col]) )
 f.close()


 abs_2_id = dict([(i,j) for j,i in enumerate(_cluidx)])

 return (_cluidx, _data, abs_2_id)


def load_clu_prop_multiple(fnm, col=[2]):
 f= open(fnm,'r')
 _cluidx=[]
 _data = []



 for rl in f:
  if rl[0] in ['#','@']: continue
  if 'Info' in rl: continue

  srl=rl[:-1].split()
  _cluidx.append(int(srl[0]))
  _data.append([float(srl[ii]) for ii in col] )
 f.close()


 abs_2_id = dict([(i,j) for j,i in enumerate(_cluidx)])

 return (_cluidx, _data, abs_2_id)


def overlap_idx(nonfib_idx, rm_list, to_sort=True):

 if to_sort:
  rm_list.sort()
  nonfib_idx.sort()

 _rmidx = rm_list
 j = nonfib_idx

 k=0
 l=0
 _new_c = []

 while k<len(_rmidx) and l<len(j):
  if _rmidx[k] < j[l]:
   k+=1
  elif _rmidx[k] > j[l]:
   #_new_c.append(j[l])
   l+=1
  else:
   _new_c.append(j[l])
   k+=1
   l+=1

 #while l<len(j):
 # _new_c.append(j[l])
 # l+=1

 return _new_c

def load_paths(pathnm):

 onpath=[]
 paths = []
 p_fold = []
#loading paths
 f = open(pathnm, 'r')

 for rl in f:
  if 'Info : state' in rl:
   srl=rl[:-1].split()
   p_fold.append( float(srl[5]))

  if 'Info : onpath' in rl:
   srl=rl[:-1].split()
   onpath.append( float(srl[4]))


  if 'Info : Flow' in rl:
    srl=rl[:-1].split()
    _p=[]
    k = srl[6].replace('-->',' ').split()
    for j in k:
     _p.append(int(j))

    paths.append([float(srl[4]), _p])

 f.close()


 print len(paths),'paths loaded'

 return (paths, p_fold, onpath)



def cal_log_trans_rate (fnm, fe, abs_2_id):



 f=open(fnm,'r')
# all id are relative ids
 conn=[]
 for rl in f:
  srl = rl[:-1].split()
  if len(srl)==0: 
   print 'Error: The following node is disconnected from the other nodes'
   print rl,
   sys.exit()

 # print rl,
  conn.append(sort_unq([abs_2_id[int(i)] for i in srl[2:]]))
 f.close()

 if len(conn)!=len(fe):
  print 'Error : connectivity data dont match free energy data'
  sys.exit()


# gen transition matrix

 T=[]
 for i,ic in enumerate(conn):
  _t=[]
  for k,j in enumerate(ic):
   if j==i: 
    print 'Warning: Node',i,'contains ii element'
    continue
   # i-->j transition
   fe_i = fe[i]
   fe_j = fe[j]
   if fe_i>fe_j:
    _t_ij = 0.0
   else:
    _t_ij = fe_i - fe_j
   
  
   _t.append(_t_ij)
  
 

  T.append(_t)


 return (conn, T)

# all relative id
def cal_committor(conn, T, fe, _s, _t, gpfold, _tol = 1e-6, ntrial = 10000):

# to cal Ai->j = ki->j / sum(l!=i)(ki->l)
# need to first calculate cal sum(...)
# note that ii elements are not included in T

 N = len(T)
 if len(conn)!= N:
  print 'Error : connectivity and nodes do not match'
  sys.exit()

 _sumk= numpy.zeros(N , dtype= numpy.float64)

 for i in range(N):
  _sumk[i] = logSum(numpy.array(T[i], dtype = numpy.float64))

 print 'Info : Comput Ai->j ...'
# now calculate Ai->j
 _A=[]
 for i in range(N):
   _A.append( numpy.exp( numpy.array(T[i], dtype = numpy.float64) - _sumk[i] ))
     
# calculate Bi = - sum(j!=i, j in t, i not in s or t)(Ai->j)

 s = sort_unq(_s[:])
 t = sort_unq(_t[:])
 s_t = sort_unq(_s[:] + _t[:])

 print 'Info :  Comput Bi ...'

 B = numpy.zeros(N-len(s_t), dtype = numpy.float64)
 guessed_pfold = numpy.zeros(N-len(s_t), dtype = numpy.float64)
 Is_Nst = numpy.zeros(N, dtype = numpy.bool)
 Idx_Nst = numpy.zeros(N- len(s_t), dtype = numpy.int32)
 OtN_idx = numpy.zeros(N, dtype = numpy.int32) - 1
 i = 0 
 
 for k in range(N):
  _idxA = get_index(k, s_t)
  if _idxA < 0 :
   Is_Nst[k] = True
   Idx_Nst[i] = k   # store original id of states not in s or t
   OtN_idx[k] = i   # get new id from original id
   guessed_pfold[i] = gpfold[k]   # guessed pfold value of state k stored in new id i
   _Ak = _A[k]
   _Ck = conn[k]
   _tosum = [0.0]
   for l,v in zip(_Ck, _Ak):
    if get_index(l, t)>=0:
      _tosum.append(v) 

   _tosum.sort()
   B[i] = -1.0 * sum(_tosum)
   i+=1

# construct Matrix A, P and B 
 print 'Info : Construcut Aux Matrix with dimension ', N-len(s_t)

 conn_con=[]
 T_con= []

 for ik in Idx_Nst:
  _Cik = conn[ik]
  _Aik = _A[ik]
  _pl = [[OtN_idx[ik], -1.0]] + \
     [[OtN_idx[_ic], _ia] for _ic, _ia in zip(_Cik, _Aik) if Is_Nst[_ic]]

  _pl.sort(key = lambda x : x[0])
  
  conn_con.append([kk[0] for kk in _pl])
  T_con.append([kk[1] for kk in _pl])



# construct sparse matrix
#i,j relative cluster ids
# Tij  transition probability of i->j
# p(t)^T . T(dt) =  p(t+dt)^T


 _iptr = 0
 _data = []
 _ind = []
 _indptr=[0]

 for ic, it in zip(conn_con, T_con):
  if len(ic)!=len(it):
   print 'Error : unmatched row when converting transition matrix'
   sys.exit()

  _data+=it
  _ind+=ic
  _iptr+=len(ic)
  _indptr.append(_iptr)

 #print _data
 #print _ind
 #print _indptr



 _Am =  scipy.sparse.csr_matrix(\
(numpy.array(_data, numpy.float64),\
 numpy.array(_ind, numpy.int32),\
 numpy.array( _indptr, numpy.int32))\
 , shape=(len(conn_con), len(conn_con)) )

 print 'Info : guessed pfold', guessed_pfold

 print 'Info :  Calculate commitors Pfold...'
 print 'Info : Initial Opt ...'
 sys.stdout.flush()
 x, info = scipy.sparse.linalg.gmres(_Am, B,x0=guessed_pfold, callback = print_max, maxiter = N/50 )

 print 'Info : Improved opt ...'
 x, info = scipy.sparse.linalg.lgmres(_Am, B, x0 =x, callback = print_max, maxiter = ntrial, tol = _tol )


 if info>0:
  print 'Info : converged'
 elif info==0:
  print 'Info : finished but not converged'
 else:
  print 'Info : Error in solving sparse systems'
  sys.exit()



 p_fold = numpy.zeros(N, dtype = numpy.float64)
 
 p_fold[s] = 0.0
 p_fold[t] = 1.0

 for kk,v in zip(Idx_Nst, x):
  p_fold[kk] = v

 return p_fold

###check strongly connected components in a directed graph
def strong_connected_components(_conn): 
 _n_node = len(_conn)
 mark =[-1]* _n_node
 left_order = [-1]* _n_node
 right_order = [-1]* _n_node

 ##first round of depth-first search to populate _left_order and _right_order for each node

 label = 0
 for i in range(_n_node):
  if mark[i]>=0: continue
  label = depth_first_search_nonrecur(i,_conn, mark, _label = label,\
    _left_order = left_order, _right_order = right_order )

# print left_order
# print right_order
 _sort_right_order = [[i,j] for i,j in enumerate(right_order)]
 _sort_right_order.sort(key= lambda(x): x[1])
 _sort_right_order.reverse()

 mark=[-1]*_n_node

#reverse _conn
 _conn_reverse = directed_graph_reverse(_conn)
# print _conn_reverse

 label = 0
 for i in range(_n_node):
  i_sort = _sort_right_order[i][0]
  if mark[i_sort]>=0: continue
  depth_first_search_nonrecur(i_sort,_conn_reverse, mark, _label = label)
#  print mark
  label+=1

 if min(mark)<0:
  print 'Error: in search for strongly connected components, results do not look right.'
  sys.exit()

 re=[[] for i in range(max(mark)+1)]
 for i in range(_n_node):
  re[mark[i]].append(i)

 re.sort(key=lambda(x):len(x))
 re.reverse()
 return re

def depth_first_search_nonrecur(i,_conn, _mark, _label = 1,\
    _left_order = None, _right_order = None ):

 _mark[i] = _label
# print _mark

 if not _left_order is None:
  _left_order[i] = _label
  _label+=1

 _stack=[i]
 _stack_ptr=[0]
 while len(_stack)>0:

  _ci = _stack[-1]
  _ptri = _stack_ptr[-1]
  if _ptri>=len(_conn[_ci]):
# all children searched
   if not _left_order is None:
    _right_order[_ci] = _label
    _label+=1
   _stack.pop()
   _stack_ptr.pop()
   continue

  i_next = _conn[_ci][_ptri]
  _stack_ptr[-1]+=1
  if _mark[i_next]<0:
   _mark[i_next]=_label
   _stack.append(i_next)
   _stack_ptr.append(0)
   if not _left_order is None:
    _left_order[i_next] = _label
    _label+=1

 return _label

  


def depth_first_search(i,_conn, _mark, _label = 1,\
    _left_order = None, _right_order = None ):
 
 _mark[i] = _label
# print _mark

 if not _left_order is None:
  _left_order[i] = _label
  _label+=1

 for k in _conn[i]:
  if _mark[k]>=0: continue
  if not _left_order is None:
   _label = depth_first_search(k, _conn, _mark, _label, _left_order, _right_order) 
  else:
   depth_first_search(k, _conn, _mark, _label)

 if not _left_order is None:
  _right_order[i] = _label
  _label+=1

 return _label
  

def directed_graph_reverse(_conn):
 _n_node = len(_conn)

 _conn_re = [[] for i in range(_n_node)]

 for i in range(_n_node):
  for j in _conn[i]:
   _conn_re[j].append(i)

 for i in range(_n_node):
  _conn_re[i] = sort_unq(_conn_re[i])

 return _conn_re

def find_connection_between_str_comp(_conn, _groups):

## assuming groups are arranged according to numbers of elements contained
 _n_node = len(_conn)
 _bl_group = [-1]*_n_node
 for i,j in enumerate(_groups):
  for k in j:
   _bl_group[k]=i

 _grp_conn=[[] for i in _groups]
 _grp_conn_detail = [[] for i in _groups]

 for i,k in enumerate(_conn):
  for j in k:
   if _bl_group[i]!=_bl_group[j]: # from two groups
    _gi = _bl_group[i]
    _gj = _bl_group[j]
    
    _idx = get_index(_gj, _grp_conn[_gi])
    if _idx<0:  
     _ins_id = 0
     for __i in _grp_conn[_gi]:
      if _gj<__i: break
      _ins_id+=1
     
      
     _grp_conn[_gi].insert(_ins_id,_gj)
     _grp_conn_detail[_gi].insert(_ins_id,[(i,j)])
          
 
    else: 
     _grp_conn_detail[_gi][_idx].append((i,j))

 return (_grp_conn, _grp_conn_detail )

def merge_network(_conna, _connb):
 if len(_conna)!=len(_connb):
  print "Error: two networks to be merged have different numbers of nodes"
  sys.exit()

 _conn_m = [[] for i in _conna]
 for i in range(len(_conna)):
  _conn_m[i]=sort_unq( _conna[i]+_connb[i])
   
 return _conn_m

def remove_node(_ni, _conn, _conn_v, is_list=False, throwaway =False ):
 if throwaway and not is_list:
  print 'Error: throwaway needs is_list=True'
  sys.exit()
 if len(_conn)!=len(_conn_v):
  print 'Error: connectivities do not match'
  sys.exit()

 if is_list:
  _i_list = sort_unq(_ni)

 _re_conn=[[] for i in _conn]
 _re_conn_v=[[] for i in _conn]
 for i in range(len(_conn)):
   if is_list:
    if get_index(i, _i_list)>=0: continue
   elif i == _ni:
    continue
   for j,k in zip(_conn[i],_conn_v[i]):
    if is_list:
     if get_index(j, _i_list)>=0: continue
    elif j == _ni:
     continue

    _re_conn[i].append(j)
    _re_conn_v[i].append(k)

 if throwaway:
  _new = [i for j,i in enumerate(_re_conn) if get_index(j, _i_list)<0]
  _new_v = [i for j,i in enumerate(_re_conn_v) if get_index(j, _i_list)<0]
  n_clu = len(_conn)
  _keep_list = remove_idx(range(n_clu), _i_list)

  _abs_to_rela = dict([(j,i) for i,j in enumerate(_keep_list)])

  for i in range(len(_new)):
   _new[i]=[_abs_to_rela[k] for k in _new[i]]

  return (_new, _new_v)
 else:
  return (_re_conn, _re_conn_v)


def load_transition_network(f_nm, addSelf=False):
######
##load transition count matrix, not nessarity undirected##
#######
##return a disctionary
##key of the dictionary :
##conn: connectivity
##cluidx: abs idx for each node
##conn_v: transition probability/count for each connection
##conn_back: backward connectivity


 _cluidx=[]
 _conn_v=[]
 _conn=[]
 _conn_back=[]

 read_idx = True
 for rl in open(f_nm,'r'):

  srl=rl[:-1].split()
# if len(srl)==0: continue

  if read_idx:
   _cluidx.append(int(srl[0])) ## append adsolute id of cluster
   read_idx = False
   _conn.append([int(i) for i in srl[1:]])
  else:
   read_idx = True
   _conn_v.append([float(i) for i in srl])

## sanity check
 if len(_cluidx)!=len(_conn) or len(_cluidx)!=len(_conn_v):
  print 'Error: Numbers of clusters look not right.'
  sys.exit()

 n_clu = len(_conn)

 for k,i,j in zip(range(n_clu), _conn, _conn_v):
  if len(i)!=len(j):
   print 'Error: cluster',_cluidx[k],' has different numbers of connectivity.'
   print i
   print j
   sys.exit()


## add missing self transition nodes
 if addSelf:
  for i in range(n_clu):
   _conni = _conn[i]
   if get_index(i,_conni)<0:
    _idx = 0
    while _idx< len(_conni):
     if i< _conni[_idx]: break
     _idx+=1
    _conni.insert(_idx,i)
    _conn_v[i].insert(_idx, 0.0)

## add back back transition connectivity
 _conn_back = [[] for i in _conn]
 for i in range(n_clu):
  for j in _conn[i]:
   _conn_back[j].append(i)

 for i in range(n_clu):
  _conn_back[i]=sort_unq(_conn_back[i])

 return {"conn":_conn,"cluidx":_cluidx,\
	"conn_v":_conn_v,"conn_back":_conn_back}





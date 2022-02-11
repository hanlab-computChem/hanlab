import numpy as np

def cycle_in_graph(_conn, cycles):
# input will be connectivity array
# and an empty array [] to pass the results
# in cycles, cycles[n] will be a list of cycles with size n
# each element of this list is a list of indices of nodes in a cycle.
# all indices start with 0
 
 n=len(_conn)
 _cycles=[]
 for _i in xrange(n+1): 
  _cycles.append([])
  cycles.append([])
# _connect = np.zeros((n,n),dtype=np.int16)
# for _i in xrange(n):
#  for _j in conn[_i]:
#   _connect[_i,_j] =1
#   _connect[_j,_i] =1
 

 _cap = n
 _mark= np.ones(_cap, dtype = np.int16 )*(-1)
 _stack= np.zeros(_cap, dtype = np.int32)
 _stack_ptr = np.zeros(_cap, dtype = np.int32)

 i=0
  
 _mark[i] = 0
# print _mark

 
 _stack[0]= i
 _stack_ptr[0]=0
 len_stack=1
 while len_stack>0:
#  print _stack[:len_stack]
#  print _stack_ptr[:len_stack]
#  print
  
  _ci = _stack[len_stack-1]
  _ptri = _stack_ptr[len_stack-1]

  if _ptri>=len(_conn[_ci]):
# all children searched
#   if not _left_order is None:
#    _right_order[_ci] = _label
#    _label+=1
#   _stack.pop()
#   _stack_ptr.pop()
   _mark[_ci] = -1
   len_stack-=1
   continue

  i_next = _conn[_ci][_ptri]
  _stack_ptr[len_stack-1]+=1


  if _mark[i_next]==-1:
   _mark[i_next]= len_stack
   _stack[len_stack]= i_next
   _stack_ptr[len_stack] = 0
   len_stack+=1


 # i_next connected with any nodes on the sampled path except for _ci
 # to form a cycle?
   if len_stack>2:
    _conn_next = _conn[i_next]    
    _t1 = _mark[_ci]
    _t2 = _mark[i_next]
    _mark[_ci]= -1
    _mt = _mark[_conn_next]
    
    _mark[_ci]= _t1
    _mark[i_next] = _t2
    
    _mx = np.max(_mt)

    if _mx>=0:
#     print 'found'
     _cycle_size = len_stack  - _mx
     _cycles[_cycle_size].append(np.copy(_stack[_mx:len_stack]))
     _mark[i_next] = -1
     len_stack-=1


#   if not _left_order is None:
#    _left_order[i_next] = _label
#    _label+=1

 
#remove redundancy

 _all_s = []
 _sort_cycles = [[] for i in xrange(n+1)]
 for _i in xrange(n+1):
  _s_c = _sort_cycles[_i]
  _t_c = _cycles[_i]
  _c = cycles[_i]

  for _j in _t_c:
#   print _j,
   _s = np.sort(_j)
   _same = False
   for _k in _s_c:
    if np.array_equal(_k, _s):
     _same = True
     break
   if not _same:
    _s_c.append(_s)

    _all_same = False

    
    for _k in _all_s:
     # print _k
      n_k = len(_k)
      xy = np.intersect1d(_k,_s)
      if n_k==len(xy):
       _all_same = True
       break
  
    if not _all_same:
     _all_s.append(_s)

     _c.append(_j)
  #   print '*',

 #  print 

def long_short_path(conn):
 #print i
 n=len(conn)
 
 cycles = [[] for _i in xrange(n)]

 _connect = np.zeros((n,n),dtype=np.int16)
 for _i in xrange(n):
  for _j in conn[_i]:
   _connect[_i,_j] =1
   _connect[_j,_i] =1
  _connect[_i,_i]=0

 _cap = 2*n
 mark= np.ones(n, dtype = np.int16 )
 que= np.zeros(_cap, dtype = np.int32)
 que_back = np.zeros(_cap, dtype = np.int32)
 que_level= np.zeros(_cap, dtype = np.int32)
 _path = np.zeros(_cap, dtype = np.int32)
 m_len = -1
 for i in xrange(n):
# finding all pathways starting from i and 
# search for possible cycles formed along the pathways

  mark[:]=1
  que[:]=0
  que_back[:]=-1
  que_level[:]=0


  len_que = 0  
  que[len_que] = i
  que_level[len_que] = 0
  len_que+=1
  mark[i]=0 

  head=0
#  mark[i]=1

 #print i,que
  while head<len_que:
   i_c= que[head]
   i_lev=que_level[head]
   #print head,':',i_c,'|',

#   print que[head:len_que]
#   print que_back[head:len_que]
#   print que_level[head:len_que]
#   print

   _con_i = np.where(_connect[i_c]*mark>0)

   _con_i_idx = _con_i[0]
   _n_con = len(_con_i_idx)


   if _n_con>0:
    que[len_que:len_que + _n_con] = _con_i_idx[:]
    que_level[len_que:len_que + _n_con] = i_lev +1
    que_back[len_que:len_que + _n_con] = head
    len_que+= _n_con
    mark[_con_i_idx] = 0
 
    if i_lev + 1 > m_len:
     m_len = i_lev + 1
     __k=0
     _path[__k] = que[len_que-1]  
     _l = que_back[len_que-1]
     while _l>=0:
       __k+=1
       _path[__k] = que[_l]
       _l = que_back[_l]

   head+=1

 return np.copy(_path[:m_len+1])


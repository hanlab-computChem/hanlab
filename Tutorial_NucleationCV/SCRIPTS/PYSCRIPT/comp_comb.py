import numpy as np
import random
import sys

def _select_hl(n):
 h = random.randint(0,n-1)
 l = random.randint(0,n-1)
 while h==l:
  l = random.randint(0,n-1)

 return (h,l)


def _reshuffle_comp(n, contmat, _frA, _frB, verbose=False):
# n # of identical monomers
# _contact_A/B n * n * nres * nres 4D array storing adjacency matrix of 
# each pair of monomers
 _contact_A = contmat[_frA]

 _contact_B = contmat[_frB]

 _count_A = np.zeros((n,n), dtype= np.int32 )
 _count_B = np.zeros((n,n), dtype= np.int32 ) 
 _shared  = np.zeros((n,n), dtype= np.int32 )
 _all     = np.zeros((n,n), dtype= np.int32 )

 _shared_h_p = np.zeros(n, dtype = np.int32)
 _shared_l_p = np.zeros(n, dtype = np.int32)
 _all_h_p = np.zeros(n, dtype = np.int32)
 _all_l_p = np.zeros(n, dtype = np.int32)

 Q_shared=0
 Q_all =0
 for _i in range(n):
  for _j in range(n):
   _count_A[_i,_j] = _contact_A[_i,_j,:,:].sum()
   _count_B[_i,_j] = _contact_B[_i,_j,:,:].sum()
   _sh = (_contact_A[_i,_j,:,:]*_contact_B[_i,_j,:,:]).sum()
   _shared[_i,_j] = _sh
   _all[_i, _j] = _count_A[_i,_j] + _count_B[_i,_j] - _sh
   Q_shared += _sh
   Q_all += _all[_i,_j]

#   print _i,_j
#   print _contact_A[_i,_j,:,:]
#   print _contact_B[_i,_j,:,:]
#   print _count_A[_i,_j], _count_B[_i,_j]
#   print _all[_i,_j], _shared[_i,_j]
#   print

 _idx = np.array(range(n), dtype= np.int32)
 _idx_p = np.array(range(n), dtype= np.int32)
 _fail_in_a_row =0

# print 'all',Q_all,'shared',Q_shared
# MC cycle
 cnt=0
 Q_zero = 0
 while _fail_in_a_row < 30*n:
  cnt+=1
  h, l = _select_hl(n)
  _idx_h = _idx[h]
  _idx_l = _idx[l]

#  print _fail_in_a_row,':'
#  print _shared
#  print _all
#  print h , l
#  print _idx_h, _idx_l

  Q_shared_p = Q_shared - 2*_shared[_idx_h,:].sum() - 2* _shared[_idx_l,:].sum()\
		+ _shared[_idx_h, _idx_h] + _shared[_idx_h, _idx_l]\
                + _shared[_idx_l, _idx_h] + _shared[_idx_l, _idx_l]

  Q_all_p = Q_all - 2*_all[_idx_h,:].sum() - 2* _all[_idx_l,:].sum()\
                + _all[_idx_h, _idx_h] + _all[_idx_h, _idx_l]\
                + _all[_idx_l, _idx_h] + _all[_idx_l, _idx_l]

#  print Q_shared_p, Q_all_p
#  sys.exit()
 
  _idx_p[:] = _idx[:]

  _t = _idx_p[h]
  _idx_p[h] = _idx_p[l]
  _idx_p[l] = _t

#  print _idx
#  print _idx_p
#  sys.exit()

  _idx_h_p = _idx_p[h]
  _idx_l_p = _idx_p[l]

#  print _count_A
#  print _count_B 
#  print

  _count_A_h = _count_A[h ,:]
  _count_A_l = _count_A[l, :]
  _count_B_h = _count_B[_idx_h_p, :]
  _count_B_l = _count_B[_idx_l_p, :]
  _count_B_h_resh = _count_B_h[_idx_p]
  _count_B_l_resh = _count_B_l[_idx_p]

#  print _count_A_h
#  print _count_A_l
#  print(_idx_p)
#  print(_count_B_h)
#  print(_count_B_l)
#  print(_count_B_h_resh)
#  print(_count_B_l_resh)
#  sys.exit()

  _eva_idx_h = np.where((_count_A_h + _count_B_h_resh)>0)
  _eva_idx_l = np.where((_count_A_l + _count_B_l_resh)>0)

#  print _eva_idx_h
#  print _eva_idx_l
#  sys.exit()

  _contact_A_h = _contact_A[h, _eva_idx_h[0],:,:]
  _contact_B_for_h = _contact_B[_idx_h_p, _idx_p[_eva_idx_h], :,:]
  _contact_A_l = _contact_A[l, _eva_idx_l[0], :,:]
  _contact_B_for_l = _contact_B[_idx_l_p, _idx_p[_eva_idx_l], :,:]

#  print _contact_B[_idx_h_p,:,:,:]
#  print _contact_B_for_h

#  sys.exit()


  _shared_h_p[:]=0
  _all_h_p[:]=0
#  print _contact_A_h.shape
#  print _contact_B_for_h.shape
#  print _contact_A.shape
#  print _contact_B.shape 
  for _i, _j in enumerate(_eva_idx_h[0]):
#   print _i, _j, _idx_p[_j]
#   print _contact_A_h[_i]
#   print _contact_B_for_h[_i]

   _sh = (_contact_A_h[_i]*_contact_B_for_h[_i]).sum()
   _shared_h_p[_idx_p[_j]] = _sh
   _all_h_p[_idx_p[_j]] = _contact_A_h[_i].sum() + _contact_B_for_h[_i].sum() - _sh

#  sys.exit()
  
  _shared_l_p[:]=0
  _all_l_p[:]=0
  for _i, _j in enumerate(_eva_idx_l[0]):
   _sh = (_contact_A_l[_i]*_contact_B_for_l[_i]).sum()
   _shared_l_p[_idx_p[_j]] = _sh
   _all_l_p[_idx_p[_j]] = _contact_A_l[_i].sum() + _contact_B_for_l[_i].sum() - _sh

  Q_shared_p+= (2*_shared_h_p.sum() + 2*_shared_l_p.sum()\
		- _shared_h_p[_idx_h_p] - _shared_h_p[_idx_l_p]\
		- _shared_l_p[_idx_h_p] - _shared_l_p[_idx_l_p])

  Q_all_p+= (2*_all_h_p.sum() + 2*_all_l_p.sum()\
                - _all_h_p[_idx_h_p] - _all_h_p[_idx_l_p]\
                - _all_l_p[_idx_h_p] - _all_l_p[_idx_l_p])

  if float(Q_all_p) == 0.0:
    print('_frA:',_frA)
  if float(Q_shared_p)/float(Q_all_p)> float(Q_shared)/float(Q_all):
    if verbose:
     print('A',h,':','B',_idx_h,'and A',l,': B',_idx_l,"<-->",)
     print('A',h,':','B',_idx_h_p,'and A',l,': B',_idx_l_p,'with Qrel from',)
     print(float(Q_shared)/float(Q_all),'to', float(Q_shared_p)/float(Q_all_p))
 

    _idx[h] = _idx_l
    _idx[l] = _idx_h
    Q_shared = Q_shared_p
    Q_all    = Q_all_p

    _shared[_idx_h_p, :] = _shared_h_p[:]
    _shared[:, _idx_h_p] = _shared_h_p[:]
    _shared[_idx_l_p, :] = _shared_l_p[:]
    _shared[:, _idx_l_p] = _shared_l_p[:]

    _all[_idx_h_p, :] = _all_h_p[:]
    _all[:, _idx_h_p] = _all_h_p[:]
    _all[_idx_l_p, :] = _all_l_p[:]
    _all[:, _idx_l_p] = _all_l_p[:]

    _fail_in_a_row = 0

   
  else:
    _fail_in_a_row +=1
  #else:
   #Q_zero += 1
   #if Q_zero >= 30*n:
    #_fail_in_a_row = 30*n
    #float(Q_all) == 0.0 

 print('Converged after ',cnt,'MC cycles')
 print('Re:')
 print('Best Q', float(Q_shared)/float(Q_all))
 print('Order', _idx)
 #if float(Q_all) == 0.0:
 #  return 1
 #else:
 return 1-(float(Q_shared)/float(Q_all))

'''
#main 

#generate random tetramer

n = 18
nres = 6

_contact_A = np.zeros((n,n,nres,nres),dtype=np.int32)
_contact_B = np.zeros((n,n,nres,nres),dtype=np.int32)

for _i in range(n):
 for _j in range(_i,n):
  _contact_B[_i,_j, : , :] = np.random.randint(2, size=(nres,nres))
  if _i!=_j:
   _contact_B[_j,_i,:,:] = _contact_B[_i,_j,:,:]

# generating A from B by reversing monomer oder

_a = np.arange(n)
np.random.shuffle(_a)
print(_a)

for _i in range(n):
 for _j in range(n):
  _contact_A[_i, _j,:,:] = _contact_B[_a[_i], _a[_j],:,:]

  #print(_i,_j)
  #print(_contact_A[_i, _j,:,:])
  #print(_contact_B[_a[_i], _a[_j],:,:])

_reshuffle_comp(n, _contact_A, _contact_B)
'''





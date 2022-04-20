import sys
import numpy

_n = numpy.zeros(740,dtype=numpy.int32)

C1 = float(sys.argv[3]) # monomer conc in mM
n_lim = 18   # size of cluster from which explicit attachment and dettachment frequency are not available
Ce = 0.36 # critical conc mM
n = int(sys.argv[2]) # 2<=n
x = Ce/C1 # x=g_inf/f_inf

if x>=1.0 or n<2:
 print 'C1 cannot be smaller than Ce=',Ce,'mM'
 print 'n cannot be smaller than',2


for rl in open("crd_n_beta_m.xvg","r"):
 srl = rl[:-1].split()
 _n[int(srl[0])] = int(srl[2])

_t = numpy.zeros(n_lim+1, dtype=numpy.float64)
_cnt = numpy.zeros(n_lim+1, dtype = numpy.int32)
_fwd = numpy.zeros(n_lim+1, dtype=numpy.float64)
_fwd_cnt = numpy.zeros(n_lim+1, dtype = numpy.int32)
_bk = numpy.zeros(n_lim+1, dtype=numpy.float64)
_bk_cnt = numpy.zeros(n_lim+1, dtype = numpy.int32)
_g_f = numpy.zeros(n_lim+1, dtype=numpy.float64)

_load_trj =False
_cur_n = -1
_cur_st = 0.
_cur_end = 0.
_new_trj = False
_new_seg = True
_trj_cnt = 0
_n_last = -1
_t_last = -1.0
_stay=0.0
_ave_t=0.0
for rl in open(sys.argv[1],"r"):

 if 'Trj' in rl:
  if _new_trj:
   _t[_cur_n]+=(_cur_end - _cur_st)
   if _cur_end > _cur_st:
    _cnt[_cur_n]+= 1

  _new_trj = True

  _load_trj = True

  _new_seg = True
  _trj_cnt+=1
  _ave_t+=_cur_end
  if _trj_cnt>1000: break
  continue


 if _load_trj:
  srl = rl[:-1].split()
  _cur_end = float(srl[0])
  _c_n = _n[int(srl[1])]

  if _new_seg:
   _cur_st = _cur_end
   _cur_n = _c_n
   _new_seg = False
   _n_last = _c_n
   _t_last = _cur_end
   _stay = 0.0

  else:

   if _cur_n!= _c_n:
    if _stay>0.0:
     _t[_cur_n]+=_stay
     _cnt[_cur_n]+=1

    _cur_st = _cur_end
    _cur_n = _c_n

    
    

   _stay = _cur_end - _cur_st

   if _c_n > _n_last:
    _fwd[_n_last]+=(_cur_end - _t_last)
    _fwd_cnt[_n_last]+=1
   elif _c_n< _n_last:
    _bk[_n_last]+=(_cur_end - _t_last)
    _bk_cnt[_n_last]+=1

   _n_last = _c_n
   _t_last = _cur_end

print 'trj count:',_trj_cnt
print 'mfpt:',_ave_t/_trj_cnt
#print _t
#print _cnt
print 'mfpt decomposition'
print 'n		g_n		f_n		g_n/f_n'
#_a_t=0.0
for i in xrange(2,n_lim,1):
 _g_f[i] = float(_bk_cnt[i])/float(_fwd_cnt[i]) 
 print '%d		%f		%f		%f'\
   %(i,float(_fwd_cnt[i])/_trj_cnt,\
       float(_bk_cnt[i])/_trj_cnt,\
       _g_f[i])

#print _a_t

_a = 1.0
for _i in xrange(n-1,1,-1):
 
 if _i<n_lim:
  _a=_a*_g_f[_i] + 1.0
 else:
  _a=_a*x + 1.0

_b = 1.0
_c = 1.0
for _i in xrange(n_lim-1,1,-1):

  _b=_b*_g_f[_i] + 1.0
  _c*= _g_f[_i]

print 'p(%d)'%n,_a/(_b+_c*x/(1.0-x))



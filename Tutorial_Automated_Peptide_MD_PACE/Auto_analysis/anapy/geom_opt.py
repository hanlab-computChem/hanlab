import math
import sys
import numpy

def xp_box(xp, box):
    _xp = xp + 0.0
    _dx = _xp - _xp[0]
    for d in range(3):
        _idx = numpy.where(_dx[:, d] > 0.5 * box[d])
        _xp[_idx, d] -= box[d]
        _idx = numpy.where(_dx[:, d] < -0.5 * box[d])
        _xp[_idx, d] += box[d]
    return _xp

def comass(xp):
    n_f = 1. / float(len(xp))
    _com = numpy.zeros(3)
    for i in range(3):
        _com[i] = xp[:, i].sum()
    return _com * n_f

def asphericty(FMOC_position):
    _tensor = numpy.dot(FMOC_position.T, FMOC_position)

    _tensor*=(1.0/len(FMOC_position))

    w, v = numpy.linalg.eig(_tensor)

    return 1.5*(w[0]**4+w[1]**4+w[2]**4)/(w[0]**2+w[1]**2+w[2]**2)**2-0.5


def p_ax(FMOC_position):
    _tensor = numpy.dot(FMOC_position.T, FMOC_position)

    _tensor*=(1.0/len(FMOC_position))

    w, v = numpy.linalg.eig(_tensor)
    return (w,v)

def cal_dih(_idx,xp,edge_len):

 xp_b = xp_box(xp[_idx,:], edge_len)
 return dihedral(xp_b[0], xp_b[1],\
                 xp_b[2], xp_b[3])


def load_dih_val(dih_val, dih_idx, xp, edge_len):

 for _idx in dih_idx:
  dih_val.append(cal_dih(_idx, xp, edge_len))


def cal_ang(_idx,xp,edge_len):

 xp_b = xp_box(xp[_idx,:], edge_len)
 return angle(xp_b[0], xp_b[1],\
                 xp_b[2])


def load_ang_val(dih_val, dih_idx, xp, edge_len):

 for _idx in dih_idx:
  dih_val.append(cal_ang(_idx, xp, edge_len))


def _dcopy( xp_src):
 _xp=[]
 for i in xp_src:
  _xp.append(i[:])
 return _xp

def  _MC(_xp, _idx):
 cen = centroid([_xp[i] for i in _idx])
 for i in _xp:
  i[0]-=cen[0]
  i[1]-=cen[1]
  i[2]-=cen[2]

def dist(x,y):
 d=x[0]-y[0]
 d2=d*d
 d=x[1]-y[1]
 d2+=d*d
 d=x[2]-y[2]
 d2+=d*d
 return d2**0.5

def dist2(x,y):
 d=x[0]-y[0]
 d2=d*d
 d=x[1]-y[1]
 d2+=d*d
 d=x[2]-y[2]
 d2+=d*d
 return d2

def dist2_box(x,y,b):
 d= abs(x[0]-y[0])
 if d>0.5*b[0]: d=b[0]-d
 r2=d*d
 d= abs(x[1]-y[1])
 if d>0.5*b[1]: d=b[1]-d
 r2+=d*d
 d= abs(x[2]-y[2])
 if d>0.5*b[2]: d=b[2]-d
 r2+=d*d

 return r2


def draw_line(x,y,mol=0):
    print( "graphics %d cylinder {"%mol+\
    "{0:.2f} {1:.2f} {2:.2f}".format\
    (x[0],x[1], x[2])+\
    "} {"+ "{0:.2f} {1:.2f} {2:.2f}"\
    .format(y[0],y[1],y[2])+"} radius 0.5 resolution 20 filled yes")


def calxyz( c1, c2, c3,  b,  a,  d):
    ci=[0.,0.,0.]
    r=[0.,0.,0.]
    h=[0.,0.,0.]
    r21=[0.,0.,0.]
    r23=[0.,0.,0.]
    rb=[0.,0.,0.]
    rb21=[0.,0.,0.]
    rc=[0.,0.,0.]
   


  
    PI=3.141592654
    ca=a;
    cb=b;
    cd=d;
    ca=ca/180.*PI;
    cd=cd/180.*PI;
             
    VECSUB(r21, c1, c2);
    VECSUB(r23, c3, c2);
             
    VECOUT(r, r21, r23);

    VECOUT(h, r21, r);

             

    rv=VECLEN(r);
    hv=VECLEN(h);
    r21v=VECLEN(r21);


    VECELONG(r,r,1./rv);
             
    VECELONG(h,h,1./hv);

    VECELONG(r21,r21,1./r21v);


    rbv=cb*math.sin(ca);

    VECELONG(rb21, r21, math.cos(PI-ca)*cb); 
             
    VECELONG(rc, h, math.cos(PI-math.fabs(cd)) );
    VECELONG(rb, r, math.sin(math.fabs(cd))*_fsgn(cd) );
    VECADD (rb, rb, rc);
    VECELONG (rb, rb, rbv);


    VECADD( ci, rb, rb21);
    VECADD( ci, ci, c1);

    return ci        


def bond(ci, c1):


 r14=[0.0,0.0,0.0]
   
 VECSUB(r14,ci,c1)
   
 return VECLEN(r14)


def angle(ci,c1,c2):

   r14=[0.0,0.0,0.0]
   r21=[0.0,0.0,0.0]
   

   VECSUB(r14,ci,c1);
   r14v=VECLEN(r14);

   VECSUB(r21,c1,c2);
   r21v=VECLEN(r21);
      
   cs=VECDOT(r14,r21)/r14v/r21v;


    
   ca=math.atan(math.sqrt(1-cs*cs)/math.fabs(cs));
   PI=3.141592654
   if cs>0.0: ca=PI-ca
   return ca/PI*180

def dihedral(ci,c1,c2,c3):


    r=[0.0,0.0,0.0]
    r14=[0.0,0.0,0.0]
    r21=[0.0,0.0,0.0]
    r23=[0.0,0.0,0.0]
    ra=[0.0,0.0,0.0]
    rb=[0.0,0.0,0.0]
    rc=[0.0,0.0,0.0]
    rd=[0.0,0.0,0.0]
    

     
    VECSUB(r14, ci,c1)
  
     


    VECSUB(r21,c1,c2)
      

    r21v=VECLEN(r21)

    VECSUB(r23,c3,c2)
      
    VECOUT(r,r21,r23)

       

    cs=VECDOT(r21,r14)/r21v/r21v
    VECELONG(ra,r21,cs)

      
     
      

    cs=VECDOT(r21,r23)/r21v/r21v
    VECELONG(rc,r21,cs)

    VECSUB(rb,r14,ra)

    VECSUB(rd,r23,rc)
       
    rdv=VECLEN(rd)
    rbv=VECLEN(rb)

       

    cs=VECDOT(rb,rd)/rdv/rbv



    cd=math.atan(math.sqrt(math.fabs(1-cs*cs))/math.fabs(cs));

    PI=3.141592654
    if cs<0: cd=PI-cd
    cd*=_fsgn(VECDOT(r14,r))
    return cd/PI*180

#add a new H atom to ci if c1 c2 c3
#connecting to ci
def add1H_pyra(ci, c1, c2, c3):


  ci1=[0.,0.,0.]
  ci2=[0.,0.,0.]
  ci3=[0.,0.,0.]
  hyg=[0.,0.,0.]

  VECSUB(ci1,c1,ci)
  VECSUB(ci2,c2,ci)
  VECSUB(ci3,c3,ci)
  c1v=VECLEN(ci1)
  c2v=VECLEN(ci2)
  c3v=VECLEN(ci3)

  VECELONG(ci1,ci1,1.0/c1v)
  VECELONG(ci2,ci2, 1.0/c2v)
  VECELONG(ci3,ci3, 1.0/c3v)

  VECADD(ci1,ci1,ci2)
  VECADD(ci1,ci1,ci3)
  c1v=VECLEN(ci1)
  VECELONG (ci1,ci1,1.0/c1v)
  VECSUB(hyg,ci,ci1)


  return hyg






def _fsgn(a):
 if a>0: return 1
 if a<0: return -1
 return 0
  
def VECADD(c,a,b):
 c[0]=a[0]+b[0]
 c[1]=a[1]+b[1]
 c[2]=a[2]+b[2]

def VECSUB(c,a,b):
 c[0]=a[0]-b[0]
 c[1]=a[1]-b[1]
 c[2]=a[2]-b[2]

def VECLEN(a):     
 return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])**0.5

#define VECDOT(a,b)   ((a)[XX]*(b)[XX]+(a)[YY]*(b)[YY]+(a)[ZZ]*(b)[ZZ])

def VECDOT(a,b):
 return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def VECOUT(c,a,b):
  c[0]=a[1]*b[2]-a[2]*b[1]
  c[1]=a[2]*b[0]-a[0]*b[2]
  c[2]=a[0]*b[1]-a[1]*b[0]

def VECELONG(c,a,v):
  c[0]=a[0]*(v)
  c[1]=a[1]*(v)
  c[2]=a[2]*(v)

#define PRVEC(a)     printf("%f %f %f\n", (&a)[XX], (&a)[YY], (&a)[ZZ])
#define VECCPY(c,a)    (c)[XX]=(a)[XX];(c)[YY]=(a)[YY];(c)[ZZ]=(a)[ZZ]
#define ZMAX_DIH(a,b)  ((a)==NULL?(b):( (*(a))+(b) ))      

# rotate xp in place about z axis by deg
def rotate_in_z(xp, deg):
 deg_pi = deg/180.0*3.141592654
 for i in range(len(xp)):     
  _x = math.cos(deg_pi)*xp[i][0]  - math.sin(deg_pi)*xp[i][1] 
  _y = math.sin(deg_pi)*xp[i][0] + math.cos(deg_pi)*xp[i][1]
  xp[i][0] = _x
  xp[i][1] = _y
 
def rotate_in_y(xp, deg):
 deg_pi = deg/180.0*3.141592654
 for i in range(len(xp)):     
  _x = math.cos(deg_pi)*xp[i][0]  + math.sin(deg_pi)*xp[i][2] 
  _z = - math.sin(deg_pi)*xp[i][0] + math.cos(deg_pi)*xp[i][2]
  xp[i][0] = _x
  xp[i][2] = _z


def rotate_in_x(xp, deg):
 deg_pi = deg/180.0*3.141592654
 for i in range(len(xp)):     
  _y = math.cos(deg_pi)*xp[i][1]  - math.sin(deg_pi)*xp[i][2] 
  _z = math.sin(deg_pi)*xp[i][1] + math.cos(deg_pi)*xp[i][2]
  xp[i][1] = _y
  xp[i][2] = _z


def find_max_min(xp):
 n=len(xp)
 minb=[0.,0.,0.]
 maxb=[0.,0.,0.]
 first=True
 for x in xp:
  if first:
   first=False
   for dim,xi in enumerate(x):
    minb[dim]=xi
    maxb[dim]=xi
  else:
   for dim, xi in enumerate(x):
    if minb[dim]>xi: minb[dim]=xi
    if maxb[dim]<xi: maxb[dim]=xi
 return maxb, minb


def make_large_sbox(sbox, sbox_size, sbox_rep, pf='S'):
 na=[]
 reslist = sbox.reslist
 segid=1
 resid=1
 atoms=sbox.atoms
 for ix in range(sbox_rep[0]):
  for iy in range(sbox_rep[1]):
   for iz in range(sbox_rep[2]):
    bx = ix * sbox_size[0]
    by = iy *sbox_size[1]
    bz = iz * sbox_size[2]
    for i in reslist:
     for j in range(i[0],i[0]+i[1]):
      na.append([99999, atoms[j][1], atoms[j][2], atoms[j][3], resid,\
                 atoms[j][5]+bx, atoms[j][6]+by, atoms[j][7]+bz,
                 pf+str(segid)])
     resid+=1
     if resid>999:
      resid=1
      segid+=1
 return na

class link_cell:
 def __init__(self, box, cut):
 
  self.box_size=[]
  self.cell_dim=[]
  self.cell_size=[]
  self.box_size_org=[]
  for i in box:
   dim=int(i/cut)
   self.box_size.append(i)
   self.box_size_org.append(i)
   self.cell_dim.append(dim)
   self.cell_size.append(i/float(dim))
  
  self.cell_link=[]
  self.cell_list=[]
  for xi in range(self.cell_dim[0]):
   for yi in range(self.cell_dim[1]):
    for zi in range(self.cell_dim[2]):
     self.cell_list.append([])
     self.cell_link.append(make_cell_link(xi, yi, zi, self.cell_dim))

 def update_box(self,box):
# note box should not change too much
  for __i in range(3):
   if abs(1.0-box[__i]/self.box_size_org[__i])>0.10:
    print( 'Error! The shape of box changes too much!')
    sys.exit()

    dim=self.cell_dim[__i]
    self.box_size[__i] = box[__i]
    
    self.cell_size.append(box[__i]/float(dim))
  

    
 def load_atoms(self,xp, in_box = False):
  c_s=self.cell_size
  b_s=self.box_size
  c_d = self.cell_dim
# put x in [0,b) in one D with periodicity of b
# if __x<0.0: return __x + b
# if __x>=b: return __x-b

  _xp = xp+0.0 

  if not in_box:  
   for __i in range(3):
    _idx = numpy.where(_xp[:,__i]<0.0)
    _xp[_idx,__i]+=b_s[__i]
    _idx_l = numpy.where(_xp[:,__i]>=b_s[__i])
    _xp[_idx_l,__i]-=b_s[__i]


  for __i in range(3):
   _xp[:,__i]*=(1.0/c_s[__i])

  x_i = _xp.astype(int)

  if not in_box:
   for __i in range(3):
    _idx = numpy.where(x_i[:,__i]<0)
    x_i[:,__i]+=c_d[__i]
    _idx_l = numpy.where(x_i[:,__i]>=c_d[__i])
    x_i[:,__i]-=c_d[__i] 


# return x*cell_dy*cell_dz+\
#        y*cell_dz + z

  c_id_list = x_i[:,0]*c_d[1]*c_d[2] + x_i[:,1]*c_d[2] + x_i[:,2]

  _mark = numpy.zeros(c_d[0]*c_d[1]*c_d[2], dtype=numpy.int32)

  for i, c_id in enumerate(c_id_list):
   
#   xi = int(x_in_box(x[0],b_s[0])/c_s[0])
#   yi = int(x_in_box(x[1],b_s[1])/c_s[1])
#   zi = int(x_in_box(x[2],b_s[2])/c_s[2])
   
#   c_id = xyz2c_id(_xi[0],_xi[1],_xi[2], self.cell_dim)
   #print x, xi, yi, zi, c_id
#   if _mark[c_id]==0:
#    self.cell_list[c_id]=[]
#    _mark[c_id]=1
   self.cell_list[c_id].append(i)

#  _mark = numpy.zeros(c_d[0]*c_d[1]*c_d[2], dtype=numpy.int32)
#  _mark[c_id_list]=1

#  _o_idx = numpy.array(range(len(c_id_list)), dtype=numpy.int32)

#  _idx_exst = numpy.where(_mark==1)

#  print c_d[0]*c_d[1]*c_d[2],len(_idx_exst[0]), _idx_exst[0]
#  for __i in _idx_exst[0]:
#   _in_cell_idx = numpy.where(c_id_list==__i)
#   self.cell_list[__i]= _o_idx[_in_cell_idx]
#   for __j in self.cell_link[__i]:
#    _mark[__j]=1

#  _idx_exst_more = numpy.where(_mark==1)  

#  print c_d[0]*c_d[1]*c_d[2],len(_idx_exst_more[0]), _idx_exst[0]  
  
#  __c=[]
#  for _i in _idx_exst_more[0]:
#   _t_c=[]
#   for _j in self.cell_link[_i]:
#    _t_c.append(self.cell_list[_j])

#   self.cell_list[_i]= numpy.concatenate(_t_c)
#   __c.append(_t_c)
   
   
#  self.cell_list = __c

  __c=[]
  for _i in range(self.cell_dim[0]*self.cell_dim[1]*self.cell_dim[2]):
   _t_c=[]
   for _j in self.cell_link[_i]:
#    if _mark[_j]==1:
     _t_c+=self.cell_list[_j][:]
   __c.append(numpy.array(_t_c, dtype=numpy.int32))

  self.cell_list = __c



 def clear_atoms(self):
  __n= len(self.cell_list)
  self.cell_list=[[] for i in range(__n)]

 def contact(self, x, xp, clash_cut2, PBC_corr = True):

  re=[]
  c_s=self.cell_size
  b_s=numpy.array(self.box_size, dtype=numpy.float32)
  half_edge = 0.5* b_s
  
#  x = numpy.array(_x, dtype=numpy.float32)
#  xp = numpy.array(_xp, dtype=numpy.float32)

  xi = int(x_in_box(x[0],b_s[0])/c_s[0])
  yi = int(x_in_box(x[1],b_s[1])/c_s[1])
  zi = int(x_in_box(x[2],b_s[2])/c_s[2])


  c_id = xyz2c_id(xi,yi,zi, self.cell_dim)
  c_list= self.cell_list[c_id]
  _xp_cell = xp[c_list, :] 

#  print x,xi, yi, zi, c_id, self.cell_link[c_id]

# append atom list from neiboring cells
#  for n_l in self.cell_link[c_id]:
#   c_list+= self.cell_list[n_l]

#  print len(c_list)

# fatser routine using numpy

  _dx = x - _xp_cell

  if PBC_corr:
   for __i in range(3):
    _idx_gt_hb = numpy.where(_dx[:,__i]>half_edge[__i])
    _dx[_idx_gt_hb,__i]-=b_s[__i]
    _idx_gt_nhb = numpy.where(_dx[:,__i]<-half_edge[__i])
    _dx[_idx_gt_nhb,__i]+=b_s[__i]

  _dx2 = _dx*_dx
 
  _r2 = _dx2[:,0]+_dx2[:,1]+_dx2[:,2]

  _idx_within_rc = numpy.where(_r2<=clash_cut2)

### original implementation
#  for idx in c_list:
#   x_r = xp[idx,:]
#   d2= dist2_box(x, x_r, b_s)
#   if d2<clash_cut2:
#    re.append(idx)

### original implementation


#  return re  
#  print _idx_within_rc
  return c_list[_idx_within_rc]


def xyz2c_id(x,y,z, cell_dim):
 cell_dx=cell_dim[0]
 cell_dy=cell_dim[1]
 cell_dz=cell_dim[2]
 if x>= cell_dx: x-=cell_dx
 if x<0: x+=cell_dx
 if y>= cell_dy: y-=cell_dy
 if y<0: y+=cell_dy
 if z>= cell_dz: z-=cell_dz
 if z<0: z+=cell_dz

 #print x,y,z
 return x*cell_dy*cell_dz+\
        y*cell_dz + z

def x_in_box(x, b):
# put x in [0,b) in one D with periodicity of b

 __x = x
 if __x<0.0: return __x + b
 if __x>=b: return __x-b
 return __x

def make_in_box(_xp, b_s):
   for __i in range(3):
    _idx = numpy.where(_xp[:,__i]<0.0)
    _xp[_idx,__i]+=b_s[__i]
    _idx_l = numpy.where(_xp[:,__i]>=b_s[__i])
    _xp[_idx_l,__i]-=b_s[__i]


def make_cell_link(xi, yi, zi, cell_dim):
 t_list=[]
 t_list.append(xyz2c_id(xi,yi,zi, cell_dim))
 for i in [-1,0,1]:
  for j in [-1,0,1]:
   for k in [-1,0,1]:
    if i==0 and j==0 and k==0: continue
    c_id= xyz2c_id(xi+i,yi+j, zi+k, cell_dim)
    t_list.append(c_id)
 return t_list


def centroid(xp_list):
 re=[0.,0.,0.]
 for i in xp_list:
  VECADD(re,re,i)
 VECELONG(re,re,1./float(len(xp_list)))
 return re


#make whole of polymer by fixpoint th atom
def makewhole(xp, fixpoint, box):
 if fixpoint >= len(xp) or fixpoint<0:
  print( 'wrong idx')
  sys.exit()

 for i in range(fixpoint+1, len(xp)):
  pbc_move(xp[i], xp[i-1], box)


 for i in range(fixpoint-1, -1, -1):
  pbc_move(xp[i], xp[i+1], box)


#mov x close to ref according to box
def pbc_move(x, ref, box):
 for i in range(3):
  if x[i]-ref[i]>0.5*box[i]:
   x[i]-=box[i]
   
  if x[i]-ref[i]<-0.5*box[i]:
   x[i]+=box[i]

#convert coeff in RB dih to fourier seires
def ur2dihe(c):
 k=[]
 k.append(-c[1]-5./8.*c[5]-3./4.*c[3])
 k.append(0.5*(c[2]+c[4]))
 k.append(-5./16.*c[5]-c[3]/4.)
 k.append(c[4]/8.)
 k.append(-c[5]/16.)
 return k

#convert c6 c12 to emin rmin

def Coff2EminRmin(C6, C12):
 if C12==0.0 or C6==0.0: return (0.0,0.0)
 epsilon = C6*C6/C12/4.
 rmin = (2.*C12/C6)**(1/6.)
 return (epsilon, rmin)

def EminRmin2Coff(emin, rmin):
 return (2.*emin*(rmin**6), emin*(rmin**12))

def Coff2sigmaeps(C6, C12):
 if C12==0.0 or C6==0.0: return (0.0,0.0)
 return ((C12/C6)**(1./6.), C6**2./C12/4.)

def Ecos(phi, param):
 return param[0]*(1.+math.cos((param[1]*phi-param[2])/180.*3.141592654)) 

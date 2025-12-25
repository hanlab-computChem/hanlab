import sys
import MDAnalysis as mda
import numpy as np
from networkutil import strong_connected_components

def frag_frag_distance (_xp1, xp2, _mark, edge_len=[1e+6,1e+6,1e+6], _cut = 4.5,\
                         onlyContact=True, singleSite = False, dualBoundary=False,\
                         _cut_sec= 4.5, returnDx=False):
# cal contact matrix for coord of fragment 1 (xp1) and fragment 2 (xp2)
# _mark: flag to mark which atom of frag 2 will be investigated when checking if it is in vicinity of frag 1
# onlyContact: ture only contact matrix false: all distance matrix needed for energy evalulation  

# available option by flags
# onlyContact = True and singleSite =True: input xp1 is a verctori for reference point xp1 not a list of vectors
#                   return adjacency vector, distance square vector, a list of dx for other vectors in xp2 
# onlyContact = False and singleSite = True: return adjecency vectorm distance square vector, distance vector,
#                                             r^6 vector, r^12 vector
# onlyContact = True and singeSite = False and dualBoundary = False: xp1 is a list of vector
#                                   if returnDx=True: return adjacency matrix, distance square matrix, matrix of dx vectors
#                                    else: adjacency matrix
# onlyContact = False and singeSite = False: xp1 is a list of vector
#                                          return adjacency matrix, distance suqare matix, distance matrix, r^6 matrix, r^12 matrix
# onlyContact = True and singeSite = False and dualBoundary = True: xp1 is a list of vector
#                                   if returnDx=true: return adjacency matrix for cut1, adjacency matrix for cut2
#                                                            distance square matrix for cut1, matrix of dx vectors for cut1
#                                           else: adjacency matrix for cut1, adjacency matrix for cut2 

 if singleSite:
  xp1 = np.zeros((1,3))
  xp1[0,:] = _xp1[:]
  _n_x1 = 1
 else:
  xp1=_xp1
  _n_x1 = len(xp1)

 _n_x2 = len(xp2)
 _cut2 = _cut*_cut
 _cut2_sec = _cut_sec*_cut_sec 

 _r2 = np.zeros(_n_x2, dtype=np.float64)

 if not onlyContact: 
  _con_stat_1 = np.zeros((_n_x1, _n_x2),dtype=np.float64)
  _con_stat_6 = np.zeros((_n_x1, _n_x2),dtype=np.float64)
  _con_stat_12 = np.zeros((_n_x1, _n_x2),dtype=np.float64)


 _con_stat2 = np.zeros((_n_x1, _n_x2),dtype=np.float64)

 _con_stat = np.zeros((_n_x1, _n_x2),dtype=np.float64)

 if dualBoundary:
  _con_stat_sec = np.zeros((_n_x1, _n_x2),dtype=np.float64)


 if returnDx:
  _dx_all = np.zeros((_n_x1,_n_x2,3))

 for _j in range(_n_x1):

    if returnDx:
     _dx_from_i  = _dx_all[_j]
     _dx_from_i[:,:]  = xp2 -xp1[_j,:]
    else:
     _dx_from_i  =  xp2 -xp1[_j,:]

    _r2[:]=0.0

    for _d in range(3):
      _dx_t = _dx_from_i[:,_d]

      _idx_gt_hb = np.where(_dx_t>\
                                0.5*edge_len[_d])
      _dx_t[_idx_gt_hb]-= edge_len[_d]
      _idx_gt_hb = np.where(_dx_t<-0.5*\
                                edge_len[_d])
      _dx_t[_idx_gt_hb]+= edge_len[_d]

      _dx_from_i[:,_d] = _dx_t[:]

      _r2+= (_dx_t*_dx_t)

    if not onlyContact:
     _cs_1 = _con_stat_1[_j,:]
     _cs_6 = _con_stat_6[_j,:]
     _cs_12 = _con_stat_12[_j,:]



    _cs2 = _con_stat2[_j,:]
    _cs = _con_stat[_j,:]

    if dualBoundary:
     _cs_sec = _con_stat_sec[_j,:]

#    _idx_within_rc = np.where(np.absolute(_r2-(_cut2-0.2))<=0.2)

    _idx_within_rc = np.where(_r2<=_cut2)
    
    if dualBoundary:
     _idx_within_rc_sec = np.where(_r2<=_cut2_sec)

    _r2_cut = _r2[_idx_within_rc]+ 1e-8
    if not onlyContact:
     _r2_i = np.reciprocal(_r2_cut)
     _r1_i = np.sqrt(_r2_i)
#    _r1 = np.reciprocal(_r1_i)
     _r6_i = _r2_i*_r2_i*_r2_i
     _r12_i = _r6_i*_r6_i

     _cs_1[_idx_within_rc]+=(_r1_i*_mark[_idx_within_rc])
     _cs_6[_idx_within_rc]+=(_r6_i*_mark[_idx_within_rc])
     _cs_12[_idx_within_rc]+=(_r12_i*_mark[_idx_within_rc])



    _cs[_idx_within_rc]+=(_mark[_idx_within_rc])
# calculate r2 for all atom pairs rather than only for those within cutoff
    _cs2+=(_r2*_mark)
    #_cs2[_idx_within_rc]+=(_r2_cut*_mark[_idx_within_rc])

    if dualBoundary:
     _cs_sec[_idx_within_rc_sec]+=(_mark[_idx_within_rc_sec])

 if onlyContact:
  if singleSite:
   return _con_stat[0],_con_stat2[0],_dx_from_i
  else:
   if dualBoundary:
    if returnDx:
     return (_con_stat, _con_stat_sec, _con_stat2, _dx_all)
    else:
     return (_con_stat, _con_stat_sec)
   else:
    if returnDx:
     return (_con_stat,_con_stat2, _dx_all)
    else:
     return _con_stat
 else:
  if singleSite:
   return (_con_stat[0],_con_stat2[0], _con_stat_1[0], _con_stat_6[0],_con_stat_12[0])
  else:
   return (_con_stat,_con_stat2, _con_stat_1, _con_stat_6,_con_stat_12)



#r_cut = 15.0

def wrap_box(pep_com, box, returnClu=False, r_cut = 15.0):

#  com_p = np.mod(pep_com , box)
  com_p = pep_com
  n = len(com_p)
  com_fin = com_p*0.0


  contact_pair,_r2,_,_,_ = frag_frag_distance(com_p, com_p, np.ones(n),edge_len=box,
                                   _cut = r_cut,onlyContact=False)
  np.fill_diagonal(contact_pair, 0)
  conn = []
  for i in range(n):
      conn.append([k for k,j in enumerate(contact_pair[i]) if j>0])

  clusters = strong_connected_components(conn)
  n_clu = len(clusters)
  #print(len(clusters[0]), end=' ')
  clu_cen_fin = np.zeros((n_clu, 3))
  for _i in range(n_clu):
    clu_idx = clusters[_i]
    _r2_i = _r2[clu_idx, :][:,clu_idx]
#   cal the diameter of the cluster _i about the center of _j th monomer in the cluster
#   the monomer leading to min diameters is closest to the COM of the cluster
    _diameter = np.zeros(len(clu_idx))
    for _j in range(len(clu_idx)):
       _diameter[_j] = max(_r2_i[_j])
# monomer for cluster center
    cen_idx= clu_idx[np.argmin(_diameter)]
#    print(cen_idx)
#    clu_cen_last = comass(xp_box(com_p[clu_idx], box))
#    clu_cen_new = comass(xp_box(com_p[clu_idx], box, refPos = clu_cen_last))
#    while np.linalg.norm(clu_cen_new - clu_cen_last) > 0.1:
#      clu_cen_last[:] = clu_cen_new[:]
#      clu_cen_new[:] = comass(xp_box(com_p[clu_idx], box, refPos = clu_cen_last))
    com_fin[clu_idx,:] = xp_box(com_p[clu_idx], box, refPos = com_p[cen_idx])
  if returnClu:
      return clusters[0]
  else:
      return com_fin








def xp_box(xp, box, refPos = None):
# use the first atom as the ref
 _xp = xp + 0.0
 if refPos is None:
  _xp0 = _xp[0]
 else:
  _xp0 = refPos

 _dx = _xp - _xp0
 for d in range(3):
  _idx = np.where(_dx[:,d]>0.5*box[d])
  _xp[_idx, d]-= box[d]
  _idx = np.where(_dx[:,d]<-0.5*box[d])
  _xp[_idx, d]+= box[d]

 return _xp

def comass(xp):
 n_f = 1./float(len(xp))

 _com = np.zeros(3)


 for i in range(3):
  _com[i] = xp[:,i].sum()

 return _com*n_f


if __name__=='__main__':
  import argparse
  def parse_args():
        parser = argparse.ArgumentParser(description='Process input files and parameters.')
        parser.add_argument('--pdb_file', type=str, help='structural file (pdb)')
        parser.add_argument('--npep', type=int, help='Number of peptides')
        parser.add_argument('--input_file', type=str,  help='trajectory file (xtc)')
        parser.add_argument('--output_file', type=str, help='output trjactory (xtc)')
        parser.add_argument('--jump_step', type=int,default=10,  help='Jump step (default: 10)')
        return parser.parse_args()


# Load the universe with topology and trajectory
  args = parse_args()
  u = mda.Universe(args.pdb_file, args.input_file)
  npep = args.npep
  natom = u.atoms.n_atoms
  assert natom%npep==0, f'Error: total atom count ({natom}) is not multiple of peptide count ({npep})'
  n_pep_atom = natom // npep

  pep_com = np.zeros((npep,3))
# Define the continuous frame range (e.g., frames 100 to 200, inclusive)

# Create a writer for the output XTC file
  with mda.Writer(args.output_file, u.atoms.n_atoms) as W:
    # Iterate over the continuous range of frames
    for ts in u.trajectory[::args.jump_step]:
        # Modify coordinates (e.g., shift all atoms by 5 Å along x-axis)
        xp = ts.positions
        xp_new = xp+0.0
        edge_len = ts.dimensions[:3]
        for pep in range(npep):
          xp_new[pep*n_pep_atom:(pep+1)*n_pep_atom]=\
               xp_box(xp[pep*n_pep_atom:(pep+1)*n_pep_atom], edge_len)
          pep_com[pep] = comass(xp_new[pep*n_pep_atom:(pep+1)*n_pep_atom])
       
        com_fin = wrap_box(pep_com, edge_len) 

        for pep in range(npep):
          xp_new[pep*n_pep_atom:(pep+1)*n_pep_atom]+=(com_fin[pep] - pep_com[pep]) 

        u.atoms.positions*=0.0
        u.atoms.positions+=xp_new

   
        # Write the modified frame to the new trajectory
        W.write(u.atoms)



import numpy as np
import sys
from atom_interact_ene import frag_frag_distance
import networkutil
import glob

DEBUG = False

def pdb_line(rl, autoAssign=False, ass_idx=0):
    atominfo = []
    if (rl[0:4] != 'ATOM'):
        return 'None'
    if not autoAssign:
        if not is_idx(rl[6:11]):
            atominfo.append(99999)
        else:
            atominfo.append(int(rl[6:11]))
    else:
        atominfo.append(ass_idx)
    atominfo.append(rl[12:17].strip())
    atominfo.append(rl[17:21].strip())
    atominfo.append(rl[21:22])
    atominfo.append(int(rl[22:26]))
    atominfo.append(float(rl[30:38]))
    atominfo.append(float(rl[38:46]))
    atominfo.append(float(rl[46:54]))
    if len(rl) >= 75:
        kl = len(rl)
        if kl > 76: kl = 76
        atominfo.append(rl[72:kl].strip())
    else:
        atominfo.append('')
    return atominfo

def is_idx(term):
    if term == '': return False
    for i in term:
        if not i in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ' ']:
            return False
    return True

def count_atom(fPDB):
    _n=0
    for rl in open(fPDB, 'r'):
        _aif = pdb_line(rl)
        if _aif == 'None': continue
        _n+=1
    return _n

def get_list(fPDB, anm, n_pep_atom):
    _aid = 0
    _l = []
    for rl in open(fPDB, 'r'):
        _aif = pdb_line(rl)
        if _aif == 'None': continue
        if _aif[1] == anm:
            _l.append(_aid)
        _aid += 1
        if _aid >= n_pep_atom: break
    return np.array(_l, dtype=np.int32)

def xp_box(xp, box):
    _xp = xp + 0.0
    _dx = _xp - _xp[0]
    for d in range(3):
        _idx = np.where(_dx[:, d] > 0.5 * box[d])
        _xp[_idx, d] -= box[d]
        _idx = np.where(_dx[:, d] < -0.5 * box[d])
        _xp[_idx, d] += box[d]
    return _xp

def comass(xp):
    n_f = 1. / float(len(xp))
    _com = np.zeros(3)
    for i in range(3):
        _com[i] = xp[:, i].sum()
    return _com * n_f

def check_amyloid_contact(xp, edge_len, _n_mono):
    _n = len(xp)
    if _n % _n_mono != 0:
        print('Error! there are', _n, 'sites but each monomer has', _n_mono, 'sites.')
        sys.exit()
    _mark = np.zeros(_n, dtype=np.int32)
    _conn = [[] for i in range(_n)]
    for _i in range(_n):
        _mark[:] = 1
        _m_id = int(_i / _n_mono)
        _mark[_m_id * _n_mono:(_m_id + 1) * _n_mono] = 0
        _c, _r2, _dx = frag_frag_distance(xp[_i], xp, _mark, edge_len=edge_len, _cut=5.5, singleSite=True)
        _idx_rc = np.where(_c > 0.0)[0]
        _n_rc = len(_idx_rc)
        _r_1_rc = np.reciprocal(np.sqrt(_r2[_idx_rc]))
        _dx_rc = _dx[_idx_rc]
        _r_vecA = np.zeros((_n_rc * _n_rc, 3))
        _r_vecB = np.zeros((_n_rc * _n_rc, 3))
        _rA = np.zeros(_n_rc * _n_rc)
        _rB = np.zeros(_n_rc * _n_rc)
        _idx_A = np.zeros(_n_rc * _n_rc, dtype=np.int32)
        _idx_B = np.zeros(_n_rc * _n_rc, dtype=np.int32)
        for _j in range(_n_rc):
            _rA[_j * _n_rc:(_j + 1) * _n_rc] = _r_1_rc[_j]
            _rB[_j * _n_rc:(_j + 1) * _n_rc] = _r_1_rc
            _r_vecA[_j * _n_rc:(_j + 1) * _n_rc, :] = _dx_rc[_j]
            _r_vecB[_j * _n_rc:(_j + 1) * _n_rc, :] = _dx_rc
            _idx_A[_j * _n_rc:(_j + 1) * _n_rc] = _idx_rc[_j]
            _idx_B[_j * _n_rc:(_j + 1) * _n_rc] = _idx_rc
        _m_idx_A = (_idx_A / _n_mono).astype(np.int32)
        _m_idx_B = (_idx_B / _n_mono).astype(np.int32)
        _dot = _r_vecA * _r_vecB
        _cos = (_dot[:, 0] + _dot[:, 1] + _dot[:, 2]) * _rA * _rB
        _st = np.where(_cos < -0.9)[0]
        for _k in _st:
            if _m_idx_A[_k] == _m_idx_B[_k]: continue
            _conn[_i].append(_idx_A[_k])
            _conn[_idx_A[_k]].append(_i)
            _conn[_i].append(_idx_B[_k])
            _conn[_idx_B[_k]].append(_i)
    for _i in range(_n): _conn[_i] = networkutil.sort_unq(_conn[_i])
    _Re = networkutil.strong_connected_components(_conn)
    mark_re = np.zeros(_n, dtype=np.int32)
    for __i in _Re:
        if len(__i) > 1:
            mark_re[__i] = 1
    return mark_re

def in_fib(_clu):
    _n = sum([len(j) for j in _clu])
    _mark = np.zeros(_n, dtype=np.int32)
    for i in _clu:
        if len(i) >= 3:
            _mark[i] = 1
    return _mark

def iden_fib_frag(_clu, _mark_fib, _conn):
    for i in _clu:
        if len(i) >= 3:
            _m = _mark_fib * 0
            _m[i] = 1
            _idx = np.where(_m * _mark_fib > 0)[0]
            if len(_idx) >= 2:
                for j in range(len(_idx) - 1):
                    _conn[_idx[j]].append(_idx[j + 1])
                    _conn[_idx[j + 1]].append(_idx[j])

def iden_bundle(_fib_frag, xp_grp, box):
    n_frag = len(_fib_frag)
    if n_frag == 1: return [[0]]
    n_grp = len(xp_grp)
    _nb_vector_pos = [[] for i in range(n_frag)]
    _nb_vector_dir = [[] for i in range(n_frag)]
    for _i in range(n_frag):
        _pos = _nb_vector_pos[_i]
        _dir = _nb_vector_dir[_i]
        _idx = _fib_frag[_i]
        _n = len(_idx)
        _mark = np.ones(_n, dtype=np.int32)
        _I_m = np.identity(_n)
        if DEBUG:
            print(_idx)
        for _j in range(n_grp):
            xp = xp_box(xp_grp[_j][_idx], box)
            _nb_mat_f = frag_frag_distance(xp, xp, _mark, edge_len=box, _cut=6.5)
            _nb_mat = (_nb_mat_f - _I_m + 1e-5).astype(np.int32)
            _conn = [[] for _l in range(_n)]
            _n_nb = np.zeros(_n, dtype=np.int32)
            for _ii in range(_n):
                _nbidx = np.where(_nb_mat[_ii, :] > 0)[0]
                _nn = len(_nbidx)
                if _nn > 2:
                    _nb_mat[_ii, :] = 0
                    _nb_mat[:, _ii] = 0
            for _ii in range(_n):
                _nbidx = np.where(_nb_mat[_ii, :] > 0)[0]
                _n_nb[_ii] = len(_nbidx)
                for _jj in _nbidx:
                    _conn[_ii].append(_jj)
            _t_pos = np.zeros((_n, 3))
            _t_dir = np.zeros((_n, 3))
            _t_cnt = 0
            _mark_s = np.zeros(_n, dtype=np.int32)
            for _kk in range(_n):
                if _mark_s[_kk] == 0 and _n_nb[_kk] == 1:
                    _st_idx = _kk
                    _seq = [_st_idx]
                    _mark_s[_st_idx] = 1
                    _c_idx = _conn[_st_idx][0]
                    _seq.append(_c_idx)
                    _mark_s[_c_idx] = 1
                    while _n_nb[_c_idx] > 1:
                        _a = _conn[_c_idx][0]
                        _b = _conn[_c_idx][1]
                        if _mark_s[_a] == 0:
                            _c_idx = _a
                        elif _mark_s[_b] == 0:
                            _c_idx = _b
                        else:
                            print('Error! wrong aro group order in column.')
                            print(_seq)
                            print(_conn)
                            sys.exit()
                        _seq.append(_c_idx)
                        _mark_s[_c_idx] = 1
                    _n_seq = len(_seq)
                    if DEBUG:
                        print(_seq)
                    if _n_seq >= 2:
                        xp_1 = xp[_seq[:_n_seq - 1]]
                        xp_2 = xp[_seq[1:_n_seq]]
                        xp_pos = (xp_1 + xp_2) * 0.5
                        dxp = xp_2 - xp_1
                        dxp2 = dxp * dxp
                        r_1 = np.reciprocal(np.sqrt(dxp2[:, 0] + dxp2[:, 1] + dxp2[:, 2]))
                        dxp[:, 0] *= r_1
                        dxp[:, 1] *= r_1
                        dxp[:, 2] *= r_1
                        _t_pos[_t_cnt:_t_cnt + _n_seq - 1, :] = xp_pos
                        _t_dir[_t_cnt:_t_cnt + _n_seq - 1, :] = dxp
                        _t_cnt += _n_seq - 1
                if _t_cnt > 0:
                    _pos.append(_t_pos[:_t_cnt])
                    _dir.append(_t_dir[:_t_cnt])
                if DEBUG:
                    print(_t_cnt)
    _conn_f = [[] for _i in range(n_frag)]
    for _i in range(n_frag - 1):
        for _j in range(1, n_frag, 1):
            _cnt = 0
            _pos_i = _nb_vector_pos[_i]
            _dir_i = _nb_vector_dir[_i]
            _pos_j = _nb_vector_pos[_j]
            _dir_j = _nb_vector_dir[_j]
            for _k in range(len(_pos_i)):
                for _l in range(len(_pos_j)):
                    _cnt += _bundle_contact(_pos_i[_k], _dir_i[_k], _pos_j[_l], _dir_j[_l], box)
                    if _cnt >= 2: break
                if _cnt >= 2: break
            if _cnt >= 2:
                _conn_f[_i].append(_j)
                _conn_f[_j].append(_i)
    return networkutil.strong_connected_components(_conn_f)

def _bundle_contact(_pos_i, _dir_i, _pos_j, _dir_j, box):
    _cnt = 0
    _n_i = len(_pos_i)
    _n_j = len(_pos_j)
    _mark = np.ones(_n_j, dtype=np.int32)
    for _ii in range(_n_i):
        _con, _dummy1, _dummy2 = frag_frag_distance(_pos_i[_ii], _pos_j, _mark, edge_len=box, _cut=6.5, singleSite=True)
        _cos = np.abs(_dir_i[_ii, 0] * _dir_j[:, 0] + _dir_i[_ii, 1] * _dir_j[:, 1] + _dir_i[_ii, 2] * _dir_j[:, 2])
        _bundle_con = np.where(_con * _cos > 0.85)[0]
        _cnt += len(_bundle_con)
    return _cnt

if __name__ == '__main__':
    import argparse
    import MDAnalysis

    def parse_args():
        parser = argparse.ArgumentParser(description='Process input files and parameters.')
        parser.add_argument('--xtc', type=str, help='xtc file')
        parser.add_argument('--num_peptides', type=int, help='Number of peptides')
#        parser.add_argument('--num_atoms_per_peptide', type=int, help='Number of atoms per peptide')
        parser.add_argument('--b_time', type=float,default=0.0,  help='starting time for analysis (default: 0 ns)')
        parser.add_argument('--e_time', type=float,default=1E+8,  help='end time for analysis (default: 1E+8 ns)')
        parser.add_argument('--jump_step', type=int, default=10, help='Jump step for trajectory analysis (default: 10)')
        parser.add_argument('--out_directory', type=str, help='Directory containing output files')
        parser.add_argument('--atomnm', type=str,default = 'C', help='Atoms used for calc Shb (Default: C)')
        parser.add_argument('--suffix', type=str,default="", help='Suffix of output filename')
        parser.add_argument('--pdb_file', type=str, help='PDB')
        return parser.parse_args()
    
    args = parse_args()
#    input_dir = args.directory
    num_peptides = args.num_peptides
#    num_atoms_per_peptide = args.num_atoms_per_peptide
    st_time = args.b_time
    end_time = args.e_time
    jump_step = args.jump_step
    atomnm_calc = args.atomnm
    pdbin = args.pdb_file
    suffix = args.suffix
    output_dir = args.out_directory

    file_cores = []
    for file_name in glob.glob(f'{args.xtc}'):
        file_cores.append(file_name.split('/')[-1])
        # only try first xtc file, if the others than the first are ignored
        break
    print(file_cores)

    for core_name in file_cores:
        pdb_file = pdbin
        atom_count = count_atom(pdb_file)
        assert atom_count%num_peptides==0,f'Error: {pdb_file} contains {atom_count} atoms which is not multiple of peptide count ({num_peptides}).'
        num_atoms_per_peptide = atom_count // num_peptides
        xtc_file = f'{args.xtc}'
        c_atom_indices = get_list(pdb_file, atomnm_calc, num_atoms_per_peptide)
        o_atom_indices = get_list(pdb_file, 'O', num_atoms_per_peptide)

        if len(c_atom_indices) != len(o_atom_indices):
            c_atom_indices = c_atom_indices[:-1]

        num_peptide_bonds = len(c_atom_indices)

        c_atom_indices_full = np.zeros(num_peptides * num_peptide_bonds, dtype=np.int32)
        o_atom_indices_full = np.zeros(num_peptides * num_peptide_bonds, dtype=np.int32)

        for i in range(num_peptides):
            c_atom_indices_full[i * num_peptide_bonds:(i + 1) * num_peptide_bonds] = num_atoms_per_peptide * i + c_atom_indices
            o_atom_indices_full[i * num_peptide_bonds:(i + 1) * num_peptide_bonds] = num_atoms_per_peptide * i + o_atom_indices

#        print(c_atom_indices_full)

        with open(f'{output_dir}/SHB{suffix}.txt', 'w') as output_file:
          universe = MDAnalysis.Universe(pdb_file, xtc_file)
                    
#            ts = universe.trajectory[-1]
          for ts in universe.trajectory[::jump_step]:
            _time = universe.trajectory.time/1000.0
            if ts.frame%100==0: print(_time)
            if _time < st_time: continue
            if _time > end_time: break
            xp = ts.positions
            edge_len = ts.dimensions[:3]

            _clu_C = check_amyloid_contact(xp[c_atom_indices_full], edge_len, num_peptide_bonds)
            _clu_O = check_amyloid_contact(xp[o_atom_indices_full], edge_len, num_peptide_bonds)

            # only use C atom in peptide bond
            mark_CO = np.where(_clu_C*_clu_O>0)[0]
#            mark_CO = np.where(_clu_C>0)[0]

            output_file.write('time %f '%(_time))
            for __ii in c_atom_indices_full[mark_CO]:
                output_file.write('%d '%__ii)
            output_file.write(f'{len(mark_CO)/len(c_atom_indices_full)}\n')

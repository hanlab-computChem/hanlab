import numpy as np
import sys
from atom_interact_ene import frag_frag_distance
import glob
from wrapbox import wrap_box

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

def asphericty(FMOC_position):
    _tensor = np.dot(FMOC_position.T, FMOC_position)

    _tensor*=(1.0/len(FMOC_position))

    w, v = np.linalg.eig(_tensor)

    return 1.5*(w[0]**4+w[1]**4+w[2]**4)/(w[0]**2+w[1]**2+w[2]**2)**2-0.5

def aspectR(FMOC_position):
    _tensor = np.dot(FMOC_position.T, FMOC_position)

    _tensor*=(1.0/len(FMOC_position))

    w, v = np.linalg.eig(_tensor)
    _ma = max(w)
    _md = (w[0]*w[1]*w[2]/_ma)**0.5

    return _ma/_md

#    return 1.5*(w[0]**4+w[1]**4+w[2]**4)/(w[0]**2+w[1]**2+w[2]**2)**2-0.5




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
        parser.add_argument('--atomnm', type=str,default = 'CA', help='Atoms used for calc Shb (Default: CA)')
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

        ca_atom_indices = get_list(pdb_file, atomnm_calc, num_atoms_per_peptide)
#        o_atom_indices = get_list(pdb_file, 'O', num_atoms_per_peptide)

#        if len(c_atom_indices) != len(o_atom_indices):
#            c_atom_indices = c_atom_indices[:-1]

        num_peptide_bonds = len(ca_atom_indices)

        ca_atom_indices_full = np.zeros(num_peptides * num_peptide_bonds, dtype=np.int32)
#        o_atom_indices_full = np.zeros(num_peptides * num_peptide_bonds, dtype=np.int32)

        for i in range(num_peptides):
            ca_atom_indices_full[i * num_peptide_bonds:(i + 1) * num_peptide_bonds] = num_atoms_per_peptide * i + ca_atom_indices
#            o_atom_indices_full[i * num_peptide_bonds:(i + 1) * num_peptide_bonds] = num_atoms_per_peptide * i + o_atom_indices
#        print(ca_atom_indices_full)
#        print(c_atom_indices_full)

        with open(f'{output_dir}/Asphere{suffix}.txt', 'w') as output_file:
          universe = MDAnalysis.Universe(pdb_file, xtc_file)
                    
#            ts = universe.trajectory[-1]
          for ts in universe.trajectory[::jump_step]:
            _time = universe.trajectory.time/1000.0
            if _time < st_time: continue
            if _time > end_time: break
            xp = ts.positions
            edge_len = ts.dimensions[:3]
            xp_CA = xp[ca_atom_indices_full]
            cluidx = wrap_box(xp_CA, edge_len, returnClu=True)
            xp_max = xp_CA[cluidx]
            com = comass(xp_max)
            dx = xp_max - com            
            _asp = asphericty(dx)
            _asr = aspectR(dx)
#            _clu_C = check_amyloid_contact(xp[c_atom_indices_full], edge_len, num_peptide_bonds)
#            _clu_O = check_amyloid_contact(xp[o_atom_indices_full], edge_len, num_peptide_bonds)

            # only use C atom in peptide bond
#            mark_CO = np.where(_clu_C*_clu_O>0)[0]
#            mark_CO = np.where(_clu_C>0)[0]

            output_file.write('time %f %d %f %f\n'%(_time,len(cluidx), _asp, _asr))

#            for __ii in c_atom_indices_full[mark_CO]:
#                output_file.write('%d '%__ii)
#            output_file.write(f'{len(mark_CO)/len(c_atom_indices_full)}\n')

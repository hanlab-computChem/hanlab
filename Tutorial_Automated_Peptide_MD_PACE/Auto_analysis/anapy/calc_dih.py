import numpy as np
import sys
from atom_interact_ene import frag_frag_distance
import networkutil
import glob
from geom_opt import dihedral, xp_box, comass

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

dihnm_base=['phi','psi','theta']
dih_anm=[['-C','N','CA','C'],['N','CA','C','+N'],['SC','CA','+CA','+SC']]

def get_dih_list(fPDB, n_pep_atom):
    _dihnm=[]
    _dihs=[]
    _aid = 0
    _res_prv = ''
    _l = []
    res_atm_list=[]
    resnm=[]
    for rl in open(fPDB, 'r'):
        _aif = pdb_line(rl)
        if _aif == 'None': continue
        _rid = _aif[4]
        _rnm = _aif[2]
        _res = f'{_rnm}{_rid}'

        if _res!=_res_prv:
          if _res_prv!='':
            resnm.append(_res_prv)
            res_atm_list.append(dict([(i,j) for i,j in zip(_r_a_nm, _r_a_id)]))
          _r_a_id =[]
          _r_a_nm = []
          _res_prv = _res

        _r_a_id.append(_aid)
        # H2 in c-cap considered as CA for dih calculation
        _aif_t = 'CA' if _aif[1]=='H2' else _aif[1]
        assert not _aif_t in _r_a_nm, f'Error: duplicate atom {_aif_t} in {_rnm}{_rid}'
        _r_a_nm.append(_aif_t)

            
        _aid += 1
        if _aid >= n_pep_atom: break
    resnm.append(_res_prv)
    res_atm_list.append(dict([(i,j) for i,j in zip(_r_a_nm, _r_a_id)]))

    for rid in range(len(resnm)):
      for _di,_dj in zip(dihnm_base, dih_anm):
        dih_idlist=[]
        succ=True
        for _dk in _dj:
          _rid_t = rid
          if _dk[0]=='-':
            _rid_t-=1
            _anm_t = _dk[1:]
          elif _dk[0]=='+':
            _rid_t+=1
            _anm_t = _dk[1:]
          else:
            _anm_t = _dk
          if _rid_t<0 or _rid_t>= len(resnm):
            succ=False
            break
          _r_dict = res_atm_list[_rid_t]
          _a_list_t=[]
          if _anm_t=='SC':
            for _l in _r_dict.keys():
                      # this list will increase for new bb types
              if _l in ['C','N','O','CA','H','HN','H1','H2']: continue
              _a_list_t.append(_r_dict[_l])
            if _a_list_t==[]:
              succ=False
              break
          elif not _anm_t in _r_dict:
            succ=False
            break
          else:
            _a_list_t.append(_r_dict[_anm_t])
          dih_idlist.append(_a_list_t)
        if succ:
          _dihnm.append(f'{_di}{rid}')
          _dihs.append(dih_idlist)
#    print(resnm)
#    print(res_atm_list)
#    print(_dihnm)
#    print(_dihs)
#    exit()
    return (_dihnm, _dihs)

def out_dih(dihnm, dihs, xp_i, output_file, must_wo=None, must_w=None ):

   for i,j in zip(dihnm, dihs):
     if not must_wo is None:
       if must_wo in i: continue
     if not must_w is None:
       if not must_w in i: continue
     c0 = comass(xp_i[j[0]])
     c1 = comass(xp_i[j[1]])
     c2 = comass(xp_i[j[2]])
     c3 = comass(xp_i[j[3]])
     output_file.write(f',{int(dihedral(c0,c1,c2,c3))}')

if __name__ == '__main__':
    import argparse
    import MDAnalysis

    def parse_args():
        parser = argparse.ArgumentParser(description='Process input files and parameters.')
        parser.add_argument('--directory', type=str, help='Directory containing input files')
        parser.add_argument('--num_peptides', type=int, help='Number of peptides')
#        parser.add_argument('--num_atoms_per_peptide', type=int, help='Number of atoms per peptide')
        parser.add_argument('--b_time', type=float,default=0.0,  help='starting time for analysis (default: 0 ns)')
        parser.add_argument('--e_time', type=float,default=1E+8,  help='end time for analysis (default: 1E+8 ns)')
        parser.add_argument('--dt', type=float,default=-1.0,  help='end time for analysis (default: no dt, if specified, jump_step will be overidden)')
        parser.add_argument('--jump_step', type=int, default=10, help='Jump step for trajectory analysis (default: 10)')
        parser.add_argument('--out_directory', type=str,default = 'None', help='Directory containing output files')
#        parser.add_argument('--atomnm', type=str,default = 'C', help='Atoms used for calc Shb (Default: C)')
        parser.add_argument('--suffix', type=str,default="", help='Suffix of output filename')
        parser.add_argument('--pdb_file', type=str,default="", help='Use specified PDB instead of searching from the input directory')
        parser.add_argument('--xtc_file', type=str,default="", help='Use specified XTC instead of searching from the input directory')
        return parser.parse_args()
    
    args = parse_args()
    input_dir = args.directory
    num_peptides = args.num_peptides
#    num_atoms_per_peptide = args.num_atoms_per_peptide
    st_time = args.b_time
    end_time = args.e_time
    jump_step = args.jump_step
    __dt = args.dt
    if __dt>=0.0: jump_step = 1
#    atomnm_calc = args.atomnm
    pdbin = args.pdb_file
    suffix = args.suffix
    output_dir = input_dir if args.out_directory=='None' else args.out_directory
    xtcin = args.xtc_file

# load index for dih calculation

    pdb_file = glob.glob(f'{input_dir}/*.pdb')[0] if pdbin=="" else pdbin
    atom_count = count_atom(pdb_file)
    assert atom_count%num_peptides==0,f'Error: {pdb_file} contains {atom_count} atoms which is not multiple of peptide count ({num_peptides}).'
    num_atoms_per_peptide = atom_count // num_peptides
    dihnm, dihs = get_dih_list(pdb_file, num_atoms_per_peptide)
#    exit()
    file_cores = []
    for file_name in glob.glob(f'{input_dir}/*.xtc'):
        file_cores.append(file_name.split('/')[-1])
        # only try first xtc file, if the others than the first are ignored
        break
    print(file_cores)

    for core_name in file_cores:
        xtc_file =xtcin if xtcin!='' else f'{input_dir}/{core_name}'




#        for i in range(num_peptides):
#            c_atom_indices_full[i * num_peptide_bonds:(i + 1) * num_peptide_bonds] = num_atoms_per_peptide * i + c_atom_indices
#            o_atom_indices_full[i * num_peptide_bonds:(i + 1) * num_peptide_bonds] = num_atoms_per_peptide * i + o_atom_indices

#        print(c_atom_indices_full)

        with open(f'{output_dir}/DIH{suffix}.txt', 'w') as output_file:
          universe = MDAnalysis.Universe(pdb_file, xtc_file)
                    
#            ts = universe.trajectory[-1]
          output_file.write('time,pid')
          for ii in dihnm:
            if not 'theta' in ii: output_file.write(f',{ii}')
          for ii in dihnm: 
            if 'theta' in ii: output_file.write(f',{ii}')
          output_file.write('\n')
          _prv_time=-100.0
          for ts in universe.trajectory[::jump_step]:
            _time = universe.trajectory.time/1000.0
            if _time < st_time: continue
            if _time > end_time: break
            if __dt >=0:
              if _time< _prv_time: continue
              _prv_time= _time + __dt
            xp = ts.positions
            edge_len = ts.dimensions[:3]

            for pi in range(num_peptides):
              output_file.write(f'{universe.trajectory.time},P{pi}')
              xp_i = xp_box(xp[pi * num_atoms_per_peptide:(pi + 1) * num_atoms_per_peptide], edge_len)
              out_dih(dihnm, dihs, xp_i, output_file, must_wo='theta' )
              out_dih(dihnm, dihs, xp_i, output_file, must_w='theta' )
              output_file.write('\n')


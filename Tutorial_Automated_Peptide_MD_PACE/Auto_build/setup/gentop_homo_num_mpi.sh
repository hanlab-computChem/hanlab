#!/bin/bash

# pansheng
#source /opt/gmx2022.6.cuda/bin/GMXRC
# main cluster
source /share/scripts/gromacs-2023_cuda-12.1_plumed-2.9.0.env

python comb_top_pace.py "$1" 80 > "$1"_80.top

# pansheng use gmx_mpi
mpirun -np 1  gmx_mpi insert-molecules -ci "$1".pdb -nmol 80 -radius 0.5 -box 10 10 10 -o "$1"_80.gro


mpirun -np 1  gmx_mpi solvate -cp "$1"_80 -cs  pace-asm.ff/cg216water -p "$1"_80.top -o "$1"_sov

mpirun -np 1  gmx_mpi grompp -v -f em -c "$1"_sov -p "$1"_80.top -o em -maxwarn 1 || exit 1 

mpirun -np 1  gmx_mpi genion -s em -p "$1"_80.top -o "$1"_ion -neutral<<EOF
13
EOF
mpirun -np 1  gmx_mpi grompp -v -f em -c "$1"_ion -p "$1"_80.top -o em || exit 1
mpirun -np 1  gmx_mpi mdrun   -ntomp 4  -v -deffnm em
mpirun -np 1  gmx_mpi grompp -v -f pr -c em -p "$1"_80.top -o pr_"$1"_80 -maxwarn 10 || exit 1
mpirun -np 1  gmx_mpi mdrun   -ntomp 4 -v -deffnm pr_"$1"_80
rm *.tpr
mpirun -np 1  gmx_mpi grompp -v -f prod -c pr_"$1"_80.gro -p "$1"_80.top -o "$1".tpr -maxwarn 2 || exit 1

rm '#'*

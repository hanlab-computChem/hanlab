#!/bin/bash

source /share/scripts/gromacs-4.5.4.env
editconf -f trimer.gro -o box.gro -c -box 9 9 9 || exit 
genbox -cp box.gro -cs cg216water.gro -p pace.top -vdwd 0.235 -o sov.gro  || exit
grompp -v -f em-posres.mdp -c sov.gro -p pace.top -o sov.tpr || exit 
genion -s sov.tpr -o ion.gro -conc 0.15 -neutral -pname NA -nname CL -p pace.top << EOF
13
EOF

grompp -v -f em-posres.mdp -c ion.gro -p pace.top -o em-posres.tpr || exit

source /share/scripts/gromacs-5.0.6-avx2.env
mpirun mdrun_mpi -s em-posres.tpr -c em-posres.gro || exit

source /share/scripts/gromacs-4.5.4.env
export GMXLIB=
grompp -v -f em.mdp -c em-posres.gro -p pace.top -o em.tpr || exit

source /share/scripts/gromacs-5.0.6-avx2.env
mpirun mdrun_mpi -s em.tpr -c em.gro || exit

source /share/scripts/gromacs-4.5.4.env
grompp -v -f nvt.mdp -c em.gro -p pace.top -o nvt.tpr || exit

source /share/scripts/gromacs-5.0.6-avx2.env
mpirun mdrun_mpi -s nvt.tpr -c nvt.gro || exit

source /share/scripts/gromacs-4.5.4.env
grompp -v -f npt-no.mdp -c nvt.gro -p pace.top -o npt-no.tpr || exit

source /share/scripts/gromacs-5.0.6-avx2.env
mpirun mdrun_mpi -s npt-no.tpr -c npt-no.gro || exit

source /share/scripts/gromacs-4.5.4.env
grompp -v -f full300.mdp -c npt-no.gro -p pace.top -o md.tpr || exit




rm \#*

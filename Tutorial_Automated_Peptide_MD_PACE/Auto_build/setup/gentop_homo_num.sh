#!/bin/bash

# choose gmx or gmx_mpi depending on GMX executable installed
gmx_mpi=gmx

python comb_top_pace.py "$1" 80 > "$1"_80.top

# pansheng use gmx_mpi
$gmx_mpi insert-molecules -ci "$1".pdb -nmol 80 -radius 0.5 -box 10 10 10 -o "$1"_80.gro


$gmx_mpi solvate -cp "$1"_80 -cs  pace-asm.ff/cg216water -p "$1"_80.top -o "$1"_sov

$gmx_mpi grompp -v -f em -c "$1"_sov -p "$1"_80.top -o em -maxwarn 1 || exit 1 

# 13 is assumed to represent SOL group for genion
# in some GMX versions, this is not always true. Please change this number accordingly.

$gmx_mpi genion -s em -p "$1"_80.top -o "$1"_ion -neutral<<EOF
13
EOF
$gmx_mpi grompp -v -f em -c "$1"_ion -p "$1"_80.top -o em || exit 1
$gmx_mpi mdrun   -ntomp 4  -v -deffnm em
$gmx_mpi grompp -v -f pr -c em -p "$1"_80.top -o pr_"$1"_80 -maxwarn 10 || exit 1
$gmx_mpi mdrun   -ntomp 4 -v -deffnm pr_"$1"_80
rm *.tpr
$gmx_mpi grompp -v -f prod -c pr_"$1"_80.gro -p "$1"_80.top -o "$1".tpr -maxwarn 2 || exit 1

rm '#'*

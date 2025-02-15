source /share/scripts/gromacs-2018.4-plumed-2.5.0.env

numomp=6

gmx_mpi grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_charmm2gmx.pdb -r step3_charmm2gmx.pdb -p system.top -n index.ndx
gmx_mpi mdrun -deffnm step4.0_minimization -ntomp $numomp

# step4.1 minimization
gmx_mpi grompp -f step4.1_minimization.mdp -o step4.1_minimization.tpr -c step4.0_minimization.gro -r step3_charmm2gmx.pdb -p system.top -n index.ndx
gmx_mpi mdrun -deffnm step4.1_minimization -ntomp $numomp

# step4.2 equlibiration
gmx_mpi grompp -f step4.2_equilibration.mdp -o step4.2_equilibration.tpr -c step4.1_minimization.gro -r step3_charmm2gmx.pdb -p system.top -n index.ndx -maxwarn 1
gmx_mpi mdrun  -deffnm step4.2_equilibration -ntomp $numomp

# md run
gmx_mpi grompp -f full_300.mdp  -c step4.2_equilibration.gro -o md.tpr -p system.top -n index.ndx -maxwarn 1
gmx_mpi mdrun -deffnm md -cpi md.cpt -maxh 2 -ntomp $numomp





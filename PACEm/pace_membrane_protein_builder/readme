
### For details guidance on the parameters and usage of the PACE force field, please refer to 'PACE for Gromacs' section found at https://github.com/hanlab-pkusz/hanlab.
### Please refer to the specific procedure outlined in the 'user-guide-gromacs.pdf' document found in the 'PACE for Gromacs' package for installing gromacs v3 and the PACE force field. Afterwards, configure the gromacs v3 operating environment on the server.


### Step 1
### Generating PACE pdb and basic files.
pdb="your membrane protein name"
export GMXLIB=/share/home/XX/gro3/share/gromacs/top
~/gro3/bin/pdb2gmx -f $pdb -o pace.pdb -ignh -missing << EOF
11    # Note: Option 11 signifies the selection of the ffPACE1.5 force field. This force field, ffPACE1.5, is typically utilized for soluble protein simulations.
EOF


### Step 2
### Generate topology files for PACE.
atom_num=`tail -n 3 pace.pdb | head -n 1 | awk '{print $2}'`
resid_num=`tail -n 3 pace.pdb | head -n 1 | awk '{print $6}'`
echo $atom_num
echo $resid_num
./pys/genPairPACE13 $atom_num $resid_num pace.pdb 1 > pace.patch
python pys/insert_param.py pace.patch topol.top > topol-pace.top


### Step 3
### Enhancing backbone hydrogen bonding interaction of transmembrane helices
./pys/dssp pace.pdb | sed -e '1,25d' | cut -b 17 | awk -F "+" '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] ;print ""}}' | sed -e "s/ /N/g" > protein.dssp   # Obtain secondary structure data of membrane proteins using DSSP software. 
python pys/helix_pair_dssp.py protein.dssp topol-pace.top 1.5 > topol-pace-HB.itp  
# The input parameter "1.5" means that the helical hydrogen bond interaction is enhanced by 1.5 times. To enhance it by 2 times, input "2".
# If four consecutive amino acids are identified as a helical structure (H or G), the hydrogen bond interaction between the first and fourth amino acids will be strengthened.
# Note: Do not enhance hydrogen bonds in helical structures located in extramembrane regions. This can be controlled by manually adjusting the "protein.dssp" file. For example, if a helical structure exists in an extramembrane region, you can manually change the helical structure marker "H" to "N" for that region to ensure no enhancement is applied.


### Step 4
### Build a Martini coarse-grained membrane and insert pace.pdb into the membrane. 
python pys/insane.py  -o cg_membrane.gro -p cg.top -pbc square -box 8,8,9 -l POPC:1 -u POPC:1 -sol W  -dm 0
gmx_mpi editconf -f pace.pdb -o pace-pbc.gro -center 4 4 4.5 -box 8 8 9
gmx_mpi solvate -cp pace-pbc.gro -cs cg_membrane.gro -o system.gro -p topol-pace-HB.top

gmx_mpi grompp -f mdp/01_min.mdp -c system.gro -p topol-pace-HB.top -o sol.tpr -maxwarn 1
gmx_mpi genion -s sol.tpr  -p topol-pace-HB.top -pname NA -nname CL -neutral -conc 0.15 -o solvated.pdb  # solvation


### Step 5
### Running MD simulation
gmx_mpi grompp -f mdp/01_min.mdp -c solvated.pdb -p topol-pace-HB.top -o 01_min.tpr -maxwarn 1
gmx_mpi mdrun -deffnm  01_min

gmx_mpi grompp -f mdp/02_min.mdp -c 01_min.gro -p topol-pace-HB.top -o 02_min.tpr -maxwarn 1
gmx_mpi mdrun -deffnm  02_min

gmx_mpi grompp -f mdp/01_eq.mdp -c 02_min.gro -p topol-pace-HB.top -o 01_eq.tpr -r 02_min.gro -maxwarn 1 
gmx_mpi mdrun -deffnm  01_eq

gmx_mpi grompp -f mdp/02_eq.mdp -c 01_eq.gro -p topol-pace-HB.top -o 02_eq.tpr -r 01_eq.gro -maxwarn 1 
gmx_mpi mdrun -deffnm  02_eq

gmx_mpi grompp -f mdp/md.mdp -c 02_eq.gro -p topol-pace-HB.top -o md.tpr -r 02_eq.gro -maxwarn 1 
gmx_mpi mdrun -deffnm  md




rm cg.pdb
rm pace.pdb
rm pacem.pdb
rm cg.itp
rm pace.itp
pdbname='3gb1.pdb'   #### input pdb name


######################################
### generating PACE pdb and topology files
### For guidance on the parameters and usage of the PACE force field, please refer to 'PACE for Gromacs' section found at https://github.com/hanlab-pkusz/hanlab.
### Please refer to the specific procedure outlined in the 'user-guide-gromacs.pdf' document found in the 'PACE for Gromacs' package for installing gromacs v3 and the PACE force field. Afterwards, configure the gromacs v3 operating environment on the server.
export GMXLIB=/share/home/XX/gro3/share/gromacs/top
~/gro3/bin/pdb2gmx -f $pdbname -o pace.pdb -ignh -missing << EOF
11
EOF
### Option 11 signifies the selection of the ffPACE1.3vdW force field. This force field, ffPACE1.3vdW, is typically utilized for soluble protein simulations.
### If you want to run membrane protein simulations, please select the ffPACE1.5 force field to obtain the topology file for this protein.


atom_num=`tail -n 3 pace.pdb | head -n 1 | awk '{print $2}'`
resid_num=`tail -n 3 pace.pdb | head -n 1 | awk '{print $6}'`
echo $atom_num
echo $resid_num
./genPairPACE13 $atom_num $resid_num pace.pdb 1 > pace.patch
python insert_param.py pace.patch topol.top > topol-pace.top
sed -e '2,16d' topol-pace.top | tac | sed 1,25d | tac > pace.itp
echo "[ end ]" >> pace.itp
rm -rf \#*


##########################################
### generating MARTINI pdb and topology files
python martini.py -f pace.pdb -o cg.top -x cg.pdb -dssp ./dssp -p backbone -ff elnedyn22
cp Protein_P.itp cg.itp
### For guidance on using the MARTINI force field, please refer to http://www.cgmartini.nl.
 

###########################################
### combining UA and CG to generate PACEm pdb and topology files
################
#index.txt
#1 56        ### "1" means pace.pdb start residue ID,   "56" means pace.pdb end residue ID
#1 20        ### UA part 1 from residue 1 to 20
#42 56       ### UA part 2 from residue 42 to 56
################
python combine.py pace.pdb cg.pdb index.txt pace.itp cg.itp > PROA_P.itp
cp PROA_P.itp ../

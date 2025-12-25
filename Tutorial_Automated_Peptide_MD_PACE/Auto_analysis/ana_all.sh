#!/bin/bash

WORKDIR=$(pwd)
cd $WORKDIR
cd ../

###################################
###MAJOR parameters for analysis###
###################################

# length of peptide sequence, e.g., 4 for tetrapeptides
seqLen=3

# frame frequency for analysis
jump_step=10

# select correct directories where the results for each sequence 
# directory where the simulation results are stored
# This is for demonstration
RTDIR=$WORKDIR/example
# for acutal application, please use the following statement
#RTDIR=PATH_TO_SIM_RESULT

# directory where the analysis results are stored
# This is for demonstration
TGTDIR=$WORKDIR/example/Result
# for acutal application, please use the following statement
#TGTDIR=PATH_TO_ANA_RESULT


# name of GMX executable usually gmx or gmx_mpi
gmx_mpi=gmx

###########END OF parameter SESSION

anapyDIR=$WORKDIR/anapy
seqLenP=$((seqLen + 1))
seqMid=$(( (seqLenP / 2) + 1 ))

cd $WORKDIR
PUR='\033[0;35m'
GRE='\033[0;32m'
CYN='\033[0;36m'
GREY='\033[0;37m'
RED='\033[0;31m'
YEL='\033[0;33m'
cle='\033[0m'
for dir in $RTDIR/"$1"*/; do
    # Check if it's a directory
    dir_bn=$(basename $dir)
    if [ -d "$dir" ] && [ ${#dir_bn} -eq $seqLen ]; then
        # Print the name of the subdirectory
        echo "Processing directory: $dir"
        
        # You can add your commands here to perform actions on each subdirectory
        # For example, you could change into the directory and run a command:
        cd $dir
        # check if analysis has been done previously
        printf "for seq ${dir_bn}: "


        if [ -f "${dir_bn}.tpr" ]; then
          printf "${dir_bn}.tpr ${GRE}found${cle}; "
        else
          printf "${dir_bn}.tpr ${RED}missing${cle}\n"
          continue
        fi


        if ls *.xtc >/dev/null 2>&1; then
          printf "xtc_files ${GRE}found${cle}; "
        else
          printf "xtc_files ${RED}missing${cle}\n"
          continue
        fi


        #doing analysis
        cd $TGTDIR
        if [ -d $dir_bn ]; then
          printf ""
        else
          mkdir $dir_bn
        fi
        cd $dir_bn
        printf "${CYN}Analysis${cle}:  "
# create critical xtc file
        xtcfile=$(ls "$dir"/*.xtc 2>/dev/null | head -n 1)
        pdbfile=""
# check if any pdb file can be found, in case no pdb found
# create one
#      There will always be a pdbfile xxx.pdb in target dir
        pdbfile_A=$TGTDIR/"$dir_bn"/"$dir_bn".pdb
        if [ -f $pdbfile_A ]; then
          pdbfile=$pdbfile_A
          printf "$pdbfile ${GRE}used${cle}; "
         else
          printf "pdb ${YEL}not found${cle}; generating new... "
          pdbfile=$TGTDIR/"$dir_bn"/"$dir_bn".pdb
          $gmx_mpi trjconv -f $xtcfile  -s $dir/"$dir_bn".tpr  -o $pdbfile -pbc mol  -dump 0 <<EOF
1
EOF

          if [ -f $pdbfile ]; then
            printf "$pdbfile ${GRE}created${cle}; "
          else
            printf "pdb ${RED}missing${cle};\n"
            continue
          fi
         fi





        npep=$(grep 'CA ' $pdbfile | grep 'ACE' | wc -l | awk '{print $1}')
        printf "use $pdbfile which contains ${CYN}$npep${cle} peptides; "
        outfile=trj_200ns_wrap.xtc
        if [ -f $outfile ]; then
          printf "$outfile ${GRE}found${cle}; "
        else
          printf "$outfile ${YEL}extracting${cle} "
          $gmx_mpi trjconv -f $xtcfile  -s $dir/"$dir_bn".tpr  -o trj_200ns.xtc -dt 0.2  -b 190.0 -e 200.0 -tu ns<<EOF
1
EOF
          python $anapyDIR/wrapbox.py --pdb_file $pdbfile --npep $npep --input_file trj_200ns.xtc --output_file $outfile
          if [ -f $outfile ]; then
            printf "$outfile ${GRE}done${cle}; "
          else
            printf "$outfile ${RED}failed${cle}; "
          fi
        fi

        outfile=trj_2000ns_wrap.xtc
        if [ -f $outfile ]; then
          printf "$outfile ${GRE}found${cle}; "
        else
          printf "$outfile ${YEL}extracting${cle} "
          $gmx_mpi trjconv -f $xtcfile  -s $dir/"$dir_bn".tpr  -o trj_2000ns.xtc -dt 0.2 -b 1990.0 -e 2000.0 -tu ns<<EOF
1
EOF
          python $anapyDIR/wrapbox.py --pdb_file $pdbfile --npep $npep --input_file trj_2000ns.xtc --output_file $outfile
          if [ -f $outfile ]; then
            printf "$outfile ${GRE}done${cle}; "
          else
            # trj_2000ns_wrap.xtc must be obtained for analysis to proceed
            printf "$outfile ${RED}failed${cle}; this sequence is IGNORED \n"
            continue 
          fi
        fi


# extract final frame
        outfile="$dir_bn"_final.pdb
        if [ -f $outfile ]; then
          printf "$outfile ${GRE}found${cle}; "
        elif [ -f trj_2000ns_wrap.xtc  ]; then
          printf "$outfile ${YEL}extracting${cle} "
          $gmx_mpi trjconv -f trj_2000ns_wrap.xtc  -s $pdbfile  -o $outfile -dump 2000.0 -tu ns<<EOF
0
EOF
          if [ -f $outfile ]; then
            printf "$outfile ${GRE}done${cle}; "
          else
            printf "$outfile ${RED}failed${cle}; "
          fi
         else
          printf "trj_2000ns_wrap.xtc ${RED}missing${cle};\n"
          continue
        fi

# analyze dihedral
        outfile=DIH_agg.txt
        if [ -f $outfile ]; then
          printf "$outfile ${GRE}found${cle}; "
        else
          printf "$outfile ${YEL}analyzing${cle} "
        python $anapyDIR/calc_dih.py --suffix _agg --directory $dir --num_peptides $npep --dt 1.0 --out_directory . --xtc_file trj_2000ns.xtc --pdb_file $pdbfile
          if [ -f $outfile ]; then
            printf "$outfile ${GRE}done${cle}; "
          else
            printf "$outfile ${RED}failed${cle}; "
          fi
        fi

# analyze morphology
        outfile=Asphere_final.txt
        if [ -f $outfile ]; then
          printf "$outfile ${GRE}found${cle}; "
        elif [ -f trj_2000ns_wrap.xtc  ]; then
          printf "$outfile ${YEL}extracting${cle} "
#add here
        python $anapyDIR/calc_asph.py --xtc trj_2000ns_wrap.xtc --jump_step 1 --pdb_file $pdbfile \
        --num_peptides $npep --b_time 1990.0 --e_time 2000.0 --suffix _final --out_directory .
          if [ -f $outfile ]; then
            printf "$outfile ${GRE}done${cle}; "
          else
            printf "$outfile ${RED}failed${cle}; "
          fi
         else
          printf "trj_2000ns_wrap.xtc ${RED}missing${cle};\n"
          continue
        fi

# calculate ALT
        outfile=SHB_final_strict.txt
        if [ -f $outfile ]; then
          printf "$outfile ${GRE}found${cle}; "
        else
          printf "$outfile ${YEL}analyzing${cle} "
        python $anapyDIR/calc_shb_strict_xtc.py --xtc trj_2000ns_wrap.xtc --jump_step 1 --pdb_file $pdbfile \
        --num_peptides $npep  --suffix _final_strict --out_directory .
          if [ -f $outfile ]; then
            printf "$outfile ${GRE}done${cle}; "
          else
            printf "$outfile ${RED}failed${cle}; "
          fi
        fi

# AP calculation
        # AP analysis
        outfile=AP.txt
        if [ -f $outfile ]; then
          printf "$outfile ${GRE}found${cle}; "
        else
          printf "$outfile ${YEL}analyzing${cle} "
          xtcfile=$(ls "$dir"/*.xtc 2>/dev/null | head -n 1)
          $gmx_mpi sasa -s $dir/"$dir_bn".tpr -f $xtcfile -o _sa0.xvg -b 0.0 -e 0.0 -tu ns <<EOF
1
EOF
          $gmx_mpi sasa -s $dir/"$dir_bn".tpr -f trj_200ns_wrap.xtc -o _sa_inte.xvg  -tu ns <<EOF
1       
EOF
          $gmx_mpi sasa -s $dir/"$dir_bn".tpr -f trj_2000ns_wrap.xtc -o _sa1.xvg  -tu ns <<EOF
1
EOF

          if [ -f "_sa0.xvg" ] && [ -f "_sa_inte.xvg" ]; then
            printf "APanalysis ${GRE}done${cle}; "

            av0=$(grep -v '\#' _sa0.xvg | grep -v '\@' | awk '{sum += $2; count++} END {print sum/count}')
            av1=$(grep -v '\#' _sa1.xvg | grep -v '\@' | awk '{sum += $2; count++} END {print sum/count}')
            avint=$(grep -v '\#' _sa_inte.xvg | grep -v '\@' | awk '{sum += $2; count++} END {print sum/count}')

# av0/av1
            result=$(echo "scale=3; $av0 / $av1" | bc)
            resultA=$(echo "scale=3; $av0 / $avint" | bc)
            echo "$dir_bn  :  t=0 $av0 t=200ns $avint t=2000ns $av1 AP(200ns) $resultA AP(2000ns) $result" >> $outfile


          else
            printf "APanalysis ${RED}failed${cle}; "
          fi
        fi


        rm "#"*
        printf "\n"
    fi


done

#!/bin/bash


WORKDIR=$(pwd)
cd $WORKDIR
cd ../
RTDIR=$(pwd)

cd $WORKDIR

bash gentop_homo_num.sh "$1"  || exit 1


mv "$1".tpr  runtpr/

rm *.gro
rm *.tpr
rm *.itp
rm *.edr
rm *.trr
rm *.xtc
rm *.pdb
rm *.log
rm *.top
rm *.cpt


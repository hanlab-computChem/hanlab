#!/bin/bash


WORKDIR=$(pwd)
cd $WORKDIR
cd ../
RTDIR=$(pwd)

cd $WORKDIR

bash gentop_homo_num.sh "$1"  || exit 1

mkdir $RTDIR/$2
mkdir $RTDIR/$2/$1

mv "$1".tpr  $RTDIR/$2/$1/
cp $RTDIR/$2/$1/"$1".tpr $RTDIR/$2/$1/prod.tpr
cp run_pace.pbs $RTDIR/$2/$1/

cd $RTDIR/$2/$1/

qsub run_pace.pbs


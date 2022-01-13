step1: 

Download https://github.com/plumed/plumed2/releases/tag/v2.3.7/plumed-2.3.7.tgz

tar zxvf plumed-2.3.7.tgz -C /directory/


step2: 

mv SCRIPTS/ADO.cpp /directory/plumed-2.3.7/src/colvar/

step3: 

PLUMED Installization

./configure --prefix=/directory/plumed2.3/

make

make install

PLUMED Installation:

▶︎ step1: 

Download https://github.com/plumed/plumed2/releases/tag/v2.3.7/plumed-2.3.7.tgz

• tar zxvf plumed-2.3.7.tgz -C /directory/


▶︎ step2: 

mv SCRIPTS/PLUMED/ADO.cpp /directory/plumed-2.3.7/src/colvar/

▶︎ step3: 

PLUMED Installization

• ./configure --prefix=/directory/plumed2.3/

• make

• make install




*fftw version: 3.3.8 https://fftw.org/pub/fftw/

*GROMACS version: 5.1.4 https://manual.gromacs.org/documentation/5.1.4/download.html




Example of Metadynamics:

• Files of gro, top, itp, mdp and script of building simulation are in Folder SETUP

• Files of production run are in Folder TOPO

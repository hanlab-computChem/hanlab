import sys
sys.path.append('../basic_code/')
import networkutil

## Usage ##
'''
python 02-get-oligomer.py data/ab40-demo-cont10.xvg > data/ab40-demo-cont10-olig.xvg
'''

def find_patch(_conn):
 _re = networkutil.strong_connected_components(_conn)
 for _i in _re: print len(_i),_i

#main#
_load_frame=False
_num = 0

for rl in open(sys.argv[1],'r'):

 if '===' in rl:

  #if _load_frame: find_patch(_conn)
  print rl,
  _conn =[]
  _num = 0 
  _load_frame=True
  continue

 if _load_frame:
  _num += 1
  _conn.append([int(_i) for _i in (rl[:-1].split(':')[1]).split()[0:]])

 if _num == 100: find_patch(_conn)

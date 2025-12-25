import numpy
import sys

def sort_unq(_conn):
 if _conn==[]: return []
 _tt = _conn[:]
# print ii, _tt
 _tt.sort()
 _tc=[_tt[0]]
 for j in range(1,len(_tt),1):
  if _tt[j-1]!=_tt[j]:
   _tc.append(_tt[j])

 return _tc

def strong_connected_components(_conn): 
 _n_node = len(_conn)
 mark =[-1]* _n_node
 left_order = [-1]* _n_node
 right_order = [-1]* _n_node

 ##first round of depth-first search to populate _left_order and _right_order for each node

 label = 0
 for i in range(_n_node):
  if mark[i]>=0: continue
  label = depth_first_search_nonrecur(i,_conn, mark, _label = label,\
    _left_order = left_order, _right_order = right_order )

# print left_order
# print right_order
 _sort_right_order = [[i,j] for i,j in enumerate(right_order)]
 _sort_right_order.sort(key= lambda x: x[1])
 _sort_right_order.reverse()

 mark=[-1]*_n_node

#reverse _conn
 _conn_reverse = directed_graph_reverse(_conn)
# print _conn_reverse

 label = 0
 for i in range(_n_node):
  i_sort = _sort_right_order[i][0]
  if mark[i_sort]>=0: continue
  depth_first_search_nonrecur(i_sort,_conn_reverse, mark, _label = label)
#  print mark
  label+=1

 if min(mark)<0:
  print( 'Error: in search for strongly connected components, results do not look right.')
  sys.exit()

 re=[[] for i in range(max(mark)+1)]
 for i in range(_n_node):
  re[mark[i]].append(i)

 re.sort(key=lambda x:len(x))
 re.reverse()
 return re

def depth_first_search_nonrecur(i,_conn, _mark, _label = 1,\
    _left_order = None, _right_order = None ):

 _mark[i] = _label
# print _mark

 if not _left_order is None:
  _left_order[i] = _label
  _label+=1

 _stack=[i]
 _stack_ptr=[0]
 while len(_stack)>0:

  _ci = _stack[-1]
  _ptri = _stack_ptr[-1]
  if _ptri>=len(_conn[_ci]):
# all children searched
   if not _left_order is None:
    _right_order[_ci] = _label
    _label+=1
   _stack.pop()
   _stack_ptr.pop()
   continue

  i_next = _conn[_ci][_ptri]
  _stack_ptr[-1]+=1
  if _mark[i_next]<0:
   _mark[i_next]=_label
   _stack.append(i_next)
   _stack_ptr.append(0)
   if not _left_order is None:
    _left_order[i_next] = _label
    _label+=1

 return _label

  


  

def directed_graph_reverse(_conn):
 _n_node = len(_conn)

 _conn_re = [[] for i in range(_n_node)]

 for i in range(_n_node):
  for j in _conn[i]:
   _conn_re[j].append(i)

 for i in range(_n_node):
  _conn_re[i] = sort_unq(_conn_re[i])

 return _conn_re


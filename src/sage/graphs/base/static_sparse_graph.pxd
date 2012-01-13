ctypedef unsigned short ushort

include "../../misc/bitset_pxd.pxi"

ctypedef struct short_digraph_s:
   ushort n
   ushort * edges
   ushort ** neighbors

ctypedef short_digraph_s short_digraph[1]

cdef int init_short_digraph(short_digraph g, G) except -1
cdef void free_short_digraph(short_digraph g)

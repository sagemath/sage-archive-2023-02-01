"""
Exterior algebras backend
"""

from sage.data_structures.bitset cimport FrozenBitset
from sage.rings.integer cimport Integer

cdef inline Integer bitset_to_int(FrozenBitset X)
cdef inline FrozenBitset int_to_bitset(Integer n)
cdef inline unsigned long degree(FrozenBitset X)
cdef inline FrozenBitset leading_supp(f)
cpdef tuple get_leading_supports(tuple I)

# Grobner basis functions
cdef inline build_monomial(supp, E)
cdef inline partial_S_poly(f, g, E, int side)
cdef inline set preprocessing(list P, list G, E, int side)
cdef inline list reduction(list P, list G, E, int side)
#cpdef tuple compute_groebner(tuple I, int side)


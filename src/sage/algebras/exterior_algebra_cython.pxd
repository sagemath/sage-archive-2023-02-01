"""
Exterior algebras backend
"""

from sage.data_structures.bitset cimport FrozenBitset
from sage.rings.integer cimport Integer

cdef Integer bitset_to_int(FrozenBitset X)
cdef FrozenBitset int_to_bitset(Integer n)
cdef unsigned long degree(FrozenBitset X)
cdef FrozenBitset leading_supp(f)
cpdef tuple get_leading_supports(tuple I)

# Grobner basis functions
cdef build_monomial(supp, E)
cdef partial_S_poly(f, g, E, int side)
cdef set preprocessing(list P, list G, E, int side)
cdef list reduction(list P, list G, E, int side)


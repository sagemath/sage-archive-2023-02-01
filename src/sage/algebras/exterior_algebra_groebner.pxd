"""
Exterior algebras Gröbner bases
"""

from sage.data_structures.bitset cimport FrozenBitset
from sage.rings.integer cimport Integer
from sage.algebras.clifford_algebra_element cimport CliffordAlgebraElement
from sage.structure.parent cimport Parent

cdef Integer bitset_to_int(FrozenBitset X)
cdef FrozenBitset int_to_bitset(Integer n)
cdef long degree(FrozenBitset X)
cdef FrozenBitset leading_supp(CliffordAlgebraElement f)
cpdef tuple get_leading_supports(tuple I)

# Grobner basis functions
cdef build_monomial(Parent E, FrozenBitset supp)
cdef partial_S_poly(CliffordAlgebraElement f, CliffordAlgebraElement g, Parent E, int side)
cdef set preprocessing(list P, list G, Parent E, int side)
cdef list reduction(list P, list G, Parent E, int side)


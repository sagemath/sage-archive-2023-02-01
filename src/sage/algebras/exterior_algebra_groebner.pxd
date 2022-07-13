"""
Exterior algebras Gröbner bases
"""

from sage.data_structures.bitset cimport FrozenBitset
from sage.rings.integer cimport Integer
from sage.algebras.clifford_algebra_element cimport CliffordAlgebraElement
from sage.structure.parent cimport Parent
from sage.structure.element cimport MonoidElement

cdef long degree(FrozenBitset X)
cdef CliffordAlgebraElement build_monomial(Parent E, FrozenBitset supp)

# Grobner basis functions
cdef class GroebnerStrategy:
    cdef Parent E  # the exterior algebra
    cdef int side
    cdef MonoidElement ideal

    cdef inline FrozenBitset leading_supp(self, CliffordAlgebraElement f)
    cdef inline partial_S_poly_left(self, CliffordAlgebraElement f, CliffordAlgebraElement g)
    cdef inline partial_S_poly_right(self, CliffordAlgebraElement f, CliffordAlgebraElement g)
    cdef set preprocessing(self, list P, list G)
    cdef list reduction(self, list P, list G)

    # These are the methods that determine the ordering of the monomials.
    # Override these for other orderings.
    # TODO: Make them abstract methods that must be implemented in subclasses
    cdef inline Integer bitset_to_int(self, FrozenBitset X)
    cdef inline FrozenBitset int_to_bitset(self, Integer n)


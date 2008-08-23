include "../../libs/polybori/decl.pxi"

from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens
from sage.rings.polynomial.multi_polynomial_ring_generic cimport \
                                                MPolynomialRing_generic
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.structure.element cimport MonoidElement

cdef class BooleanPolynomialRing(MPolynomialRing_generic):
    cdef PBRing _pbring
    cdef Py_ssize_t* pbind
    cdef public _monom_monoid
    cdef public object __interface

cdef class BooleanPolynomial(MPolynomial):
    cdef PBPoly _pbpoly

cdef class BooleSet:
    cdef BooleanPolynomialRing _ring
    cdef PBSet _pbset

cdef class CCuddNavigator:
    cdef PBNavigator _pbnav

cdef class DD:
    cdef PBDD _pbdd

cdef class BooleanMonomial(MonoidElement):
    cdef PBMonom _pbmonom
    cdef BooleanPolynomialRing _ring

cdef class BooleanMonomialVariableIterator:
    cdef object parent
    cdef BooleanPolynomialRing _ring
    cdef BooleanMonomial obj
    cdef PBMonomVarIter _iter
    cdef PBMonomVarIter _end

cdef class BooleanMonomialIterator:
    cdef BooleanMonomial obj
    cdef PBMonomIter _iter
    cdef PBMonomIter _end

cdef class BooleanPolynomialIterator:
    cdef BooleanPolynomial obj
    cdef PBPolyIter _iter
    cdef PBPolyIter _end

cdef class BooleSetIterator:
    cdef object _parent
    cdef BooleanPolynomialRing _ring
    cdef PBSetIter _iter
    cdef PBSetIter _end
    cdef BooleSet obj

cdef class GroebnerStrategy:
    cdef GBStrategy _strat
    cdef BooleanPolynomialRing _parent

cdef class BooleanPolynomialVector:
    cdef PBPolyVector _vec
    cdef BooleanPolynomialRing _parent

cdef class BooleanPolynomialVectorIterator:
    cdef BooleanPolynomialVector obj
    cdef BooleanPolynomialRing _parent
    cdef PBPolyVectorIter _iter
    cdef PBPolyVectorIter _end

cdef class VariableBlock_base:
    cdef BooleanPolynomialRing _ring
    cdef int size
    cdef int start_index
    cdef int offset
    cdef public object __name__

cdef class BooleVariable:
    cdef PBVar _pbvar

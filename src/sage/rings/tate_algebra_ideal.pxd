from sage.rings.tate_algebra_element cimport TateAlgebraTerm
from sage.rings.tate_algebra_element cimport TateAlgebraElement

cdef _groebner_basis_buchberger(I, prec, bint integral_first)
#cdef _groebner_basis_F5(I, prec)

cdef TateAlgebraElement regular_reduce(sgb, TateAlgebraTerm s, TateAlgebraElement v, verbose, stopval)
cdef TateAlgebraElement reduce(gb, TateAlgebraElement v, verbose, stopval)

from sage.misc.sagex_ds import BinaryTree
from sage.misc.sagex_ds cimport BinaryTree
include "../../ext/interrupt.pxi"

cdef enum:
    pdUNIVAR
    pdVAR
    pdCOEFF
    pdSQR
    pdADD
    pdMUL
    pdPOW
    pdABC


cdef class poly_dag
cdef class generic_pd

cdef class CompiledPolynomialFunction:
    cdef generic_pd _dag
    cdef object _coeffs

    cdef object _parse_structure(CompiledPolynomialFunction)
    cdef poly_dag _get_gap(CompiledPolynomialFunction, BinaryTree, int)
    cdef void _fill_gaps_quick(CompiledPolynomialFunction, BinaryTree)
    cdef object eval(CompiledPolynomialFunction, object)

cdef class poly_dag:
    cdef generic_pd dag
    cdef int label
    cdef object val(poly_dag, object, object)
#    cdef generic_pd get_dag(poly_dag)


cdef class generic_pd:
    cdef object value
    cdef int refs, hits
    cdef void eval(self, x, coeffs)
    cdef void accelerate(generic_pd)

cdef class var_pd(generic_pd):
    cdef int index

cdef class univar_pd(generic_pd):
    pass

cdef class coeff_pd(generic_pd):
    cdef int index


cdef class unary_pd(generic_pd):
    cdef generic_pd operand
    cdef poly_dag operand_p
    cdef void accelerate(unary_pd)

cdef class sqr_pd(unary_pd):
    pass

cdef class pow_pd(unary_pd):
    cdef object exponent


cdef class binary_pd(generic_pd):
    cdef generic_pd left, right
    cdef poly_dag left_p, right_p
    cdef void accelerate(binary_pd)

cdef class add_pd(binary_pd):
    pass

cdef class mul_pd(binary_pd):
    pass

cdef class abc_pd(binary_pd):
    cdef int index







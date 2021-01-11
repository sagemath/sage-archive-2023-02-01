include "sage/libs/linkages/padics/lazy/flint.pxi"
include "lazy_template.pxi"

cdef inline type lazy_class_zero = pAdicLazyElement_zero
cdef inline type lazy_class_one = pAdicLazyElement_one
cdef inline type lazy_class_value = pAdicLazyElement_value
cdef inline type lazy_class_random = pAdicLazyElement_random
cdef inline type lazy_class_shift = pAdicLazyElement_shift
cdef inline type lazy_class_add = pAdicLazyElement_add
cdef inline type lazy_class_mul = pAdicLazyElement_mul
cdef inline type lazy_class_muldigit = pAdicLazyElement_muldigit
cdef inline type lazy_class_div = pAdicLazyElement_div
cdef inline type lazy_class_sqrt = pAdicLazyElement_sqrt
cdef inline type lazy_class_teichmuller = pAdicLazyElement_teichmuller
cdef inline type lazy_class_selfref = pAdicLazyElement_selfref

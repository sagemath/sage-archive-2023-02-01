include "sage/libs/linkages/padics/lazy/flint.pxi"
include "lazy_template.pxi"

cdef type lazy_class_zero = pAdicLazyElement_zero
cdef type lazy_class_one = pAdicLazyElement_one
cdef type lazy_class_value = pAdicLazyElement_value
cdef type lazy_class_random = pAdicLazyElement_random
cdef type lazy_class_shift = pAdicLazyElement_shift
cdef type lazy_class_add = pAdicLazyElement_add
cdef type lazy_class_mul = pAdicLazyElement_mul
cdef type lazy_class_muldigit = pAdicLazyElement_muldigit
cdef type lazy_class_div = pAdicLazyElement_div
cdef type lazy_class_sqrt = pAdicLazyElement_sqrt
cdef type lazy_class_teichmuller = pAdicLazyElement_teichmuller
cdef type lazy_class_selfref = pAdicLazyElement_selfref

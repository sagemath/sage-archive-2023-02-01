include "../libs/ginac/decl.pxi"

from sage.structure.element cimport CommutativeRingElement


cdef class Expression(CommutativeRingElement):
    cdef GEx _gobj
    cdef Expression coerce_in(self, z)
    cpdef object _eval_self(self, R)
    cpdef bint is_polynomial(self, var)
    cpdef bint is_relational(self)
    cpdef object pyobject(self)
    cpdef Expression _subs_expr(self, expr)

cdef Expression new_Expression_from_GEx(GEx juice)

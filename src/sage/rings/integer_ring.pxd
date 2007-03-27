include "../libs/ntl/decl.pxi"

from ring cimport PrincipalIdealDomain
from integer cimport Integer

cdef class IntegerRing_class(PrincipalIdealDomain):
    cdef Integer _coerce_ZZ(self, ntl_c_ZZ *z)

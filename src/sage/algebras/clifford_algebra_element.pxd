"""
Clifford algebra elements
"""

from sage.modules.with_basis.indexed_element cimport IndexedFreeModuleElement

cdef class CliffordAlgebraElement(IndexedFreeModuleElement):
    pass

cdef class ExteriorAlgebraElement(CliffordAlgebraElement):
    pass


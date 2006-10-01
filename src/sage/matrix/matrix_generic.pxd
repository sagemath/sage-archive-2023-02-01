from sage.structure.element cimport ModuleElement
from sage.structure.mutability_pyx cimport Mutability

cdef class Matrix(ModuleElement):
    cdef object _mutability
    cdef public object _parent
    cdef object __nrows
    cdef object __ncols
    cdef object __dict
    cdef object __determinant
    cdef object __sparse_columns
    cdef object __sparse_rows


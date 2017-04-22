from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense

cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    cpdef _hermite_form_euclidean(self, transformation=*)
    cpdef _reversed_hermite_form_euclidean(self, transformation=*)
    cpdef _weak_popov_form(self, transformation=*, shifts=*)

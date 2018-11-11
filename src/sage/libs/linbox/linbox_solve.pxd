# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: language = c++

from sage.modules.vector_modn_sparse cimport c_vector_modint

cdef class Linbox_modn_sparse:
    cdef c_vector_modint *rows
    cdef size_t nrows
    cdef size_t ncols
    cdef unsigned int modulus

    cdef set(self, int modulus, size_t nrows, size_t ncols, c_vector_modint *rows)
    cdef void solve(self, c_vector_modint **x, c_vector_modint *b, int method)

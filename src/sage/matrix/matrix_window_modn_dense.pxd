from matrix_window cimport MatrixWindow

cdef extern from "../ext/multi_modular.h":
    ctypedef unsigned long mod_int
    mod_int MOD_INT_MAX

cdef class MatrixWindow_modn_dense(MatrixWindow):
    pass
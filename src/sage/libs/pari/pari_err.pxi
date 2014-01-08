# Let Cython know that pari_err.h includes interrupt.h
cdef extern from 'interrupt.h':
    pass

from sage.libs.pari.handle_error cimport _pari_check_warning

cdef extern from 'sage/libs/pari/pari_err.h':
    int pari_catch_sig_on() except 0
    int pari_catch_sig_str(char *) except 0
    void pari_catch_sig_off()

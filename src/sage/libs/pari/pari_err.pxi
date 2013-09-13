# We don't need anything from here, as we have PARI-specific signal
# handling.  We still include this such that Cython detects the
# dependency on interrupt.h for recompiling gen.pyx.
cdef extern from 'interrupt.h':
    pass

cdef extern from 'sage/libs/pari/pari_err.h':
    int pari_catch_sig_on() except 0
    int pari_catch_sig_str(char *) except 0
    void pari_catch_sig_off()

from sage.libs.pari.gen cimport _pari_trap

################################################################
# WARNING: Only modify the version of this file in the
# sage/ext/ directory.  It overwrites all other copies
# of this file during the SAGE build.
################################################################

cdef extern from 'interrupt.h':
    void sig_on()
    void sig_str(char*)
    void sig_off()

    # These provide backwards compatibility with sage-4.6 and earlier
    int _sig_on
    void _sig_str(char*)
    int _sig_off

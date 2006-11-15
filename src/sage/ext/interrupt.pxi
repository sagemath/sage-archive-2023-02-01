################################################################
# WARNING: Only modify the version of this file in the
# sage/ext/ directory.  It overwrites all other copies
# of this file during the SAGE build.
################################################################

cdef extern from 'interrupt.h':
    int _sig_on, _sig_off, _sig_check
    void _sig_str(char*)
    void setup_signal_handler()

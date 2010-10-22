# vim:ft=pyrex:

cdef extern from 'pari_err.h':
    int _pari_catch
    int _pari_endcatch
    int sig_on "_pari_sig_on" () except 0
    int sig_str "_pari_sig_str" (char *) except 0
    void sig_off "_pari_sig_off" ()

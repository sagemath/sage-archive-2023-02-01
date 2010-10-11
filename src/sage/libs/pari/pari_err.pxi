# vim:ft=pyrex:

cdef extern from 'pari_err.h':
    int _pari_catch
    int _pari_endcatch
    void sig_on "_pari_sig_on" ()
    void sig_off "_pari_sig_off" ()
    void sig_str "_pari_sig_str" (char *)

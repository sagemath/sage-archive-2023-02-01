
cdef extern from 'pari_err.h':
    int _pari_catch
    int _pari_endcatch
    int _sig_on "_pari_sig_on"
    int _sig_off "_pari_sig_off"
    void _sig_str "_pari_sig_str" (char *)

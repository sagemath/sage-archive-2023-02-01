cdef make_ZZ(ZZ* x):
    cdef ntl_ZZ y
    _sig_off
    y = ntl_ZZ(_INIT)
    #y.x = <ZZ*> x
    y.set(<ZZ*> x)
    return y


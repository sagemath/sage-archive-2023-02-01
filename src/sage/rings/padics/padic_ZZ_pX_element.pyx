

cdef class pAdicZZpXElement(pAdicExtElement):
    def __init__(self, parent):
        self.prime_pow = parent.prime_pow
        pAdicExtElement.__init__(self, parent)

    def _ntl_rep(self):
        raise NotImplementedError

    cdef int _set_from_ZZ_pX_c(self, ZZ_pX_c poly, ntl_ZZ_pContext_class ctx) except -1:
        raise NotImplementedError

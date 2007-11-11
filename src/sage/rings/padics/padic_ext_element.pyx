cdef class pAdicExtElement(pAdicGenericElement):
    cdef ext_p_list(self, bint pos):
        raise NotImplementedError

    def _ext_p_list(self, pos):
        return self.ext_p_list(pos)

    cdef int _set_from_mpz(self, mpz_t x) except -1:
        raise NotImplementedError

    cdef int _set_from_mpq(self, mpq_t x) except -1:
        raise NotImplementedError

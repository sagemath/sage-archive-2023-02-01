from sage.libs.gmp.all cimport mpz_t

cdef mpz_stirling_s2(mpz_t s, unsigned long n, unsigned long k)

cdef list from_word(list w, list base_set)

cdef class PerfectMatchingsIterator(object):
    cdef Py_ssize_t n
    cdef bint* available
    cdef list current
    cdef bint at_end

    cpdef list next(self)

cdef list convert(Py_ssize_t* f, Py_ssize_t n)


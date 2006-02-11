

cdef class export:
    cdef int c_p1_normalize_int(self, int N, int u, int v,
                                int* uu, int* vv, int* ss,
                                int compute_s) except -1

    cdef int c_p1_normalize_llong(self, int N, int u, int v,
                                     int* uu, int* vv, int* ss,
                                     int compute_s) except -1



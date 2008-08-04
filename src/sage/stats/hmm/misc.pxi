#############################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version at your option.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

cdef extern from "string.h":
    char *strcpy(char *s1, char *s2)
    void* memcpy(void* dst, void* src, size_t len)

cdef double* to_double_array(v) except NULL:
    """
    Transform a Python list of floats to a C array of doubles.  The caller is
    responsible for deallocating the resulting memory.

    INPUT:
        v -- a list of objects coercible to floats
    OUTPUT:
        a newly allocated C array of doubles
    """
    cdef double x
    cdef double* w = <double*> safe_malloc(sizeof(double)*len(v))
    cdef Py_ssize_t i = 0
    for x in v:
        w[i] = x
        i += 1
    return w

cdef int* to_int_array(v) except NULL:
    """
    Transform a Python list of ints to a C array of ints.  The caller is
    responsible for deallocating the resulting memory.

    INPUT:
        v -- a list of objects coercible to ints
    OUTPUT:
        a newly allocated C array of ints
    """
    cdef int x
    cdef int* w = <int*> safe_malloc(sizeof(int)*len(v))
    cdef Py_ssize_t i = 0
    for x in v:
        w[i] = x
        i += 1
    return w

cdef void* safe_malloc(int bytes) except NULL:
    """
    malloc the given bytes of memory and check that the malloc
    succeeds -- if not raise a MemoryError.

    INPUT:
        bytes -- a nonnegatie integer

    OUTPUT:
        void pointer or raise a MemoryError.
    """
    cdef void* t = sage_malloc(bytes)
    if not t:
        raise MemoryError, "error allocating memory for Hidden Markov Model"
    return t


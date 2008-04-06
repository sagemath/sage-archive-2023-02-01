cdef extern from *:
    ctypedef struct PyObject:
        Py_ssize_t ob_refcnt

def ref_count(o):
    """
    Return the number of references to \code{o}.

    EXAMPLES:
      sage: from sage.misc.refcount import ref_count
      sage: ref_count(33)
      1
      sage: L = [-7] * 5
      sage: ref_count(L[0])
      6
    """
    return (<PyObject *>o).ob_refcnt

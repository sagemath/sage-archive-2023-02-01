cdef class WeakValueDictionary(dict):
    cdef __weakref__
    cdef callback
    cdef int _guard_level
    cdef list _pending_removals

    cdef int _enter_iter(self) except -1
    cdef int _exit_iter(self) except -1

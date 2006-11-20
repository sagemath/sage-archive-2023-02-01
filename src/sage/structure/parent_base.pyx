cdef class ParentWithBase(parent.Parent):
    def __init__(self, base):
        self._base = base

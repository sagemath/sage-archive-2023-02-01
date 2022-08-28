from sage.misc.classcall_metaclass cimport ClasscallMetaclass


cdef class ConstructorBaseclassMetaclass(ClasscallMetaclass):

    def __cinit__(self, *args, **opts):
        r"""
        TESTS::
        """
        if '__constructor__' in self.__dict__:
            def constructor(cls, *a, **o):
                return self.__constructor__(*a, **o)
            self.classcall = constructor
        else:
            self.classcall = None

        self.classcontains = getattr(self, "__classcontains__", None)
        self.classget = getattr(self, "__classget__", None)




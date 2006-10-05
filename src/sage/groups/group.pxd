cimport sage.structure.gens
import  sage.structure.gens

cdef class Group(sage.structure.gens.Generators):
    pass

cdef class AbelianGroup(Group):
    pass

cdef class FiniteGroup(Group):
    pass

cdef class AlgebraicGroup(Group):
    pass


import gens
cimport gens

cdef class Group(gens.Generators):
    pass

cdef class AbelianGroup(Group):
    pass

cdef class FiniteGroup(Group):
    pass

cdef class AlgebraicGroup(Group):
    pass


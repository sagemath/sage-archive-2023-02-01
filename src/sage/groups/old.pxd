cimport sage.structure.parent_gens

cdef class Group(sage.structure.parent_gens.ParentWithGens):
    pass

cdef class AbelianGroup(Group):
    pass

cdef class FiniteGroup(Group):
    pass

cdef class AlgebraicGroup(Group):
    pass


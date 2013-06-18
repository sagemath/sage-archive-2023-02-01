from sage.structure.parent cimport Parent
from sage.categories.morphism cimport Morphism
from sage.structure.element cimport AlgebraElement, Element, ModuleElement, RingElement

cdef class WrapperParent_model0(Parent):
    cdef Parent R

cdef class WrapperParent_model1(Parent):
    cdef Parent R

cdef class WrapperElement(AlgebraElement):
    cdef Element val


# It is important to keep this line here, basically to trick Pyrex.
# If you remove this line then other modules that cimport element
# from other directories will fail.

cimport sage.structure.sage_object
from sage.structure.parent cimport Parent

cimport sage_object
import  sage_object

cdef class Element(sage_object.SageObject):
    cdef Parent _parent
    cdef _richcmp_c_impl(left, Element right, int op)
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef public _richcmp(self, right, int op)
    cdef public _cmp(self, right)
    cdef _set_parent_c(self, Parent parent)
    cdef base_extend_c(self, Parent R)       # do *NOT* override, but OK to call directly
    cdef base_extend_c_impl(self, Parent R)  # OK to override, but do NOT call
    cdef _rich_to_bool(self, int op, int r)

    cpdef _act_on_(self, x, bint self_on_left)
    cpdef _acted_upon_(self, x, bint self_on_left)

cdef class ElementWithCachedMethod(Element):
    cdef public dict __cached_methods

cdef class ModuleElement(Element)       # forward declaration

cdef class RingElement(ModuleElement)   # forward declaration

cdef class ModuleElement(Element):

    cpdef ModuleElement _add_(self, ModuleElement right)
    cpdef ModuleElement _sub_(self, ModuleElement right)
    cpdef ModuleElement _neg_(self)
    # self._rmul_(x) is x * self
    cpdef ModuleElement _lmul_(self, RingElement right)
    # self._lmul_(x) is self * x, to abide with Python conventions.
    cpdef ModuleElement _rmul_(self, RingElement left)

    cdef ModuleElement _mul_long(self, long n)

    # Inplace operations, override, do *NOT* call directly
    cpdef ModuleElement _iadd_(self, ModuleElement right)
    cpdef ModuleElement _isub_(self, ModuleElement right)
    cpdef ModuleElement _ilmul_(self, RingElement right)

    # Coerce x to the base ring of self and return the result.
    cdef RingElement coerce_to_base_ring(self, x)

cdef class MonoidElement(Element):
    cpdef MonoidElement _mul_(self, MonoidElement right)

cdef class MultiplicativeGroupElement(MonoidElement):
    cpdef MultiplicativeGroupElement _div_(self, MultiplicativeGroupElement right)


cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(ModuleElement):
    cpdef RingElement _mul_(self, RingElement right)
    cpdef RingElement _div_(self, RingElement right)

    # Inplace operations, override, do *NOT* call directly
    cpdef RingElement _imul_(self, RingElement right)
    cpdef RingElement _idiv_(self, RingElement right)

    cdef RingElement _add_long(self, long n)

cdef class CommutativeRingElement(RingElement):
    pass

cdef class IntegralDomainElement(CommutativeRingElement):
    pass

cdef class DedekindDomainElement(IntegralDomainElement):
    pass

cdef class PrincipalIdealDomainElement(DedekindDomainElement):
    pass

cdef class EuclideanDomainElement(PrincipalIdealDomainElement):
    pass

cdef class FieldElement(CommutativeRingElement):
    pass

cdef class AlgebraElement(RingElement):
    pass

cdef class CommutativeAlgebraElement(CommutativeRingElement):
    pass

cdef class CommutativeAlgebra(AlgebraElement):
    pass

cdef class InfinityElement(RingElement):
    pass


cdef class Vector(ModuleElement):
    cdef Py_ssize_t _degree

    # Returns the dot product, using the simple metric $e_i \cdot e_j = \delta_{ij}$.
    cpdef Element _dot_product_(Vector left, Vector right)  # override, call if parents the same

    cpdef Vector _pairwise_product_(Vector left, Vector right) # override, call if parents the same

    cdef bint is_sparse_c(self)
    cdef bint is_dense_c(self)


cdef class Matrix(ModuleElement):
    # All matrix classes must be written in Cython
    cdef Py_ssize_t _nrows
    cdef Py_ssize_t _ncols

    cdef Vector _vector_times_matrix_(matrix_right, Vector vector_left)    # OK to override, AND call directly
    cdef Vector _matrix_times_vector_(matrix_left, Vector vector_right)    # OK to override, AND call directly
    cdef Matrix _matrix_times_matrix_(left, Matrix right)                  # OK to override, AND call directly

    cdef bint is_sparse_c(self)
    cdef bint is_dense_c(self)




cdef class CoercionModel:
    cpdef canonical_coercion(self, x, y)
    cpdef bin_op(self, x, y, op)


cdef generic_power_c(a, nn, one)

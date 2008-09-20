from sage.rings.ring cimport FiniteField
from sage.rings.ring cimport Ring
from sage.structure.element cimport FiniteFieldElement, Element, RingElement, ModuleElement

from sage.structure.parent  cimport Parent


cdef extern from "givaro/givconfig.h":
    pass

cdef extern from "givaro/givrandom.h":
    ctypedef struct GivRandom "GivRandom":
        pass

    GivRandom GivRandomSeeded  "GivRandom"(unsigned long seed)

cdef extern from "givaro/givgfq.h":
    ctypedef struct intvec "std::vector<unsigned int>":
        void (* push_back)(int elem)

    ctypedef struct constintvec "const std::vector<unsigned int>"

    intvec intvec_factory "std::vector<unsigned int>"(int len)

cdef extern from "givaro/givgfq.h":

    ctypedef struct GivaroGfq "GFqDom<int>":
        #attributes
        unsigned int one
        unsigned int zero

        # methods
        int (* mul)(int r, int a, int b)
        int (* add)(int r, int a, int b)
        int (* sub)(int r, int a, int b)
        int (* div)(int r, int a, int b)
        int (* inv)(int r, int x)
        int (* neg)(int r, int x)
        int (* mulin)(int a, int b)
        unsigned int (* characteristic)()
        unsigned int (* cardinality)()
        int (* exponent)()
        int (* random)(GivRandom gen, int res)
        int (* initi "init")(int res, int e)
        int (* initd "init")(int res, double e)
        int (* axpyin)(int r, int a, int x)
        int (* sage_generator)() # SAGE specific method, not found upstream
        int (* convert)(int r, int p)
        int (* read)(int r, int p)
        int (* axpy)(int r, int a, int b, int c)
        int (* axmy)(int r, int a, int b, int c)
        int (* amxy)(int r, int a, int b, int c)
        bint (* isZero)(int e)
        bint (* isOne)(int e)
        bint (* isunit)(int e)

    GivaroGfq *gfq_factorypk "new GFqDom<int>" (unsigned int p, unsigned int k)
    # SAGE specific method, not found upstream
    GivaroGfq *gfq_factorypkp "new GFqDom<int>" (unsigned int p, unsigned int k, intvec poly)
    GivaroGfq *gfq_factorycopy "new GFqDom<int>"(GivaroGfq orig)
    GivaroGfq  gfq_deref "*"(GivaroGfq *orig)
    void delete "delete "(void *o)
    int gfq_element_factory "GFqDom<int>::Element"()


cdef class FiniteField_givaro(FiniteField):
    cdef GivaroGfq *objectptr # C++ object
    cdef object _polynomial
    cdef object _polynomial_ring
    cdef object _prime_subfield
    cdef object _array
    cdef object _is_conway
    cdef object _hash
    cdef int repr
    cdef gen_array(FiniteField_givaro self)
    cdef order_c(FiniteField_givaro self)
    cdef _coerce_c_impl(self, x)
    cdef prime_subfield_C(FiniteField_givaro self)

cdef class FiniteField_givaro_iterator:
    cdef int iterator
    cdef FiniteField_givaro _parent

cdef class FiniteField_givaroElement(FiniteFieldElement):
    cdef int element
    cdef object __multiplicative_order
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef RingElement _mul_c_impl(self, RingElement right)
    cdef RingElement _div_c_impl(self, RingElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef FiniteField_givaroElement _new_c(self, int value)




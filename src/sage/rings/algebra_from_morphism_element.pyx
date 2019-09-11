from sage.structure.element cimport CommutativeAlgebraElement
from sage.structure.element cimport Element


cdef class AlgebraFMElement(CommutativeAlgebraElement):
    r"""
    Generic class for elements lying in ring extensions

    AUTHOR:

    - Xavier Caruso (2016)
    """
    def __init__(self, parent, x, *args, **kwds):
        from sage.rings.algebra_from_morphism import AlgebraFromMorphism
        if not isinstance(parent, AlgebraFromMorphism):
            raise TypeError("%s is not an instance of AlgebraFromMorphism" % parent)
        if isinstance(x, AlgebraFMElement):
            x = x._backend()
        try:
            parentx = x.parent()
            if parent.base().has_coerce_map_from(parentx):
                x = parent.base().coerce_map_from(parentx)(x)
                x = parent.defining_morphism()(x)
        except AttributeError:
            pass
        Element.__init__(self, parent)
        ring = parent._backend()
        self._element = ring(x, *args, **kwds)

    def _repr_(self):
        r"""
        Return a string representation of this element

        By default, it is the same as the strict representation
        of this element viewed as an element of the underlying
        ring.

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = L.gen()
            sage: E(x)._repr_()
            'z4'
            sage: x._repr_()
            'z4'
        """
        return str(self._element)

    def _latex_(self):
        r"""
        Return a latex representation of this element

        By default, it is the same as the latex representation
        of this element viewed as an element of the underlying
        ring.

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = L.gen()
            sage: E(x)._latex_()
            'z_{4}'
            sage: x._latex_()
            'z_{4}'
        """
        from sage.misc.latex import latex
        return str(latex(self._element))

    def _backend(self):
        return self._element

    cpdef _add_(self,other):
        r"""
        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.random_element()
            sage: y = E.random_element()
            sage: (x+y).parent() is E
            True
        """
        return self.__class__(self._parent, self._element + other._backend())

    cpdef _neg_(self):
        return self.__class__(self._parent, -self._element)

    cpdef _sub_(self,other):
        r"""
        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.random_element()
            sage: y = E.random_element()
            sage: (x-y).parent() is E
            True
        """
        return self.__class__(self._parent, self._element - other._backend())

    cpdef _mul_(self,other):
        r"""
        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.random_element()
            sage: y = E.random_element()
            sage: (x*y).parent() is E
            True
        """
        return self.__class__(self._parent, self._element * other._backend())

    #cpdef _div_(self,other)
    #cpdef _floordiv_(self,other)

    def additive_order(self):
        r"""
        Return the additive order of this element

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.gen()
            sage: x.additive_order()
            5
        """
        return self._element.additive_order()

    def multiplicative_order(self):
        r"""
        Return the multiplicite order of this element

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.gen()
            sage: x.multiplicative_order()
            624
        """
        return self._element.multiplicative_order()

    def is_unit(self):
        r"""
        Return ``True`` if this element is a unit in the ring
        of the parent extension, ``False`` otherwise

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = RingExtension(A,QQ)
            sage: E(4).is_unit()
            True
            sage: E(x).is_unit()
            False
        """
        return self._element.is_unit()

    def is_nilpotent(self):
        r"""
        Return ``True`` if this element is nilpotent in the ring
        of the parent extension, ``False`` otherwise

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = RingExtension(A,QQ)
            sage: E(0).is_nilpotent()
            True
            sage: E(x).is_nilpotent()
            False
        """
        return self._element.is_nilpotent()

    def is_prime(self):
        r"""
        Return ``True`` if this element is a prime element 
        in the ring of the parent extension, ``False`` otherwise

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = RingExtension(A,QQ)
            sage: E(x^2+1).is_prime()
            True
            sage: E(x^2-1).is_prime()
            False
        """
        return self._element.is_prime()


cdef class RingExtensionWithBasisElement(AlgebraFMElement):
    def _repr_(self):
        parent = self._parent
        basis = parent.basis()
        _, _, j = parent.vector_space()
        coeffs = j(self)
        s = ""
        for i in range(len(basis)):
            if coeffs[i].is_zero(): continue
            c = coeffs[i]
            sign = 1
            if (-c)._is_atomic():
                c = -c
                sign = -sign
            b = basis[i]
            sign = 1
            if (-b)._is_atomic():
                b = -b
                sign = -sign
            if sign == 1:
                s += " + "
            else:
                s += " - "
            if c != 1:
                if c._is_atomic():
                    s += "%s" % c
                else:
                    s += "(%s)" % c
                if b != 1: s += "*"
            if b != 1:
                if b._is_atomic():
                    s += "%s" % b
                else:
                    s += "(%s)" % b
            if b == 1 and c == 1:
                s += "1"
        if s == "": return "0"
        return s[3:]

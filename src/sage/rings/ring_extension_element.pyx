#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.cpython.getattr cimport AttributeErrorMessage

from sage.structure.element cimport CommutativeAlgebraElement
from sage.structure.element cimport Element
from sage.rings.integer_ring import ZZ
from sage.categories.fields import Fields
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.rings.ring_extension cimport RingExtension_class, RingExtensionWithGen
from sage.rings.ring_extension_morphism cimport MapRelativeFieldToVectorSpace
from sage.rings.ring_extension_conversion cimport backend_parent, backend_element
from sage.rings.ring_extension_conversion import to_backend, from_backend


# Classes
#########

cdef class RingExtensionElement(CommutativeAlgebraElement):
    r"""
    Generic class for elements lying in ring extensions

    AUTHOR:

    - Xavier Caruso (2016)
    """
    def __init__(self, RingExtension_class parent, x, *args, **kwds):
        if not isinstance(parent, RingExtension_class):
            raise TypeError("%s is not a ring extension" % parent)
        x = backend_element(x)
        try:
            parentx = x.parent()
            if parent._base.has_coerce_map_from(parentx):
                x = parent._base.coerce_map_from(parentx)(x)
                x = parent._backend_defining_morphism(x)
        except AttributeError:
            pass
        Element.__init__(self, parent)
        ring = parent._backend
        self._backend = ring(x, *args, **kwds)

    def __reduce__(self):
        return self._parent, (self._backend,)

    def __getattr__(self, name):
        method = None
        if (<RingExtension_class>self._parent)._import_methods and hasattr(self._backend, name):
            method = getattr(self._backend, name)
        if not callable(method):
            raise AttributeError(AttributeErrorMessage(self, name))
        def wrapper(*args, **kwargs):
            output = method(*to_backend(args), **to_backend(kwargs))
            return from_backend(output, self._parent)
        wrapper.__doc__ = method.__doc__
        return wrapper

    def __dir__(self):
        d = dir(self.__class__)
        if not (<RingExtension_class>self._parent)._import_methods:
            return d
        for name in dir(self._backend):
            try:
                attribute = getattr(self._backend, name)
                if callable(attribute):
                    d.append(name)
            except:
                pass
        return sorted(set(d))

    def __hash__(self):
        return hash(self._backend)

    def __repr__(self):
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
        cdef RingExtension_class parent = self._parent
        print_as = parent._print_elements_as
        if print_as is not None:
            return str(print_as(self._backend))
        return self._repr_()

    def _repr_(self):
        return str(self._backend)

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
        return str(latex(self._backend))

    cpdef _richcmp_(left, right, int op):
        return left._backend._richcmp_(backend_element(right), op)

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
        return self._parent(self._backend + (<RingExtensionElement>other)._backend)

    cpdef _neg_(self):
        return self._parent(-self._backend)

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
        return self._parent(self._backend - (<RingExtensionElement>other)._backend)

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
        return self._parent(self._backend * (<RingExtensionElement>other)._backend)

    cpdef _div_(self,other):
        return self._parent.fraction_field()(self._backend / (<RingExtensionElement>other)._backend)

    cpdef _floordiv_(self, other):
        return self._parent(self._backend // (<RingExtensionElement>other)._backend)

    def sqrt(self, extend=True, all=False, name=None):
        return self._parent(self._backend.sqrt(extend, all, name))

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
        return self._backend.additive_order()

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
        return self._backend.multiplicative_order()

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
        return self._backend.is_unit()

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
        return self._backend.is_nilpotent()

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
        return self._backend.is_prime()

    def is_square(self, root=False):
        is_sq = self._backend.is_square()
        sq = None
        if root and is_sq:
            sq = self.sqrt()
        if root:
            return is_sq, sq
        else:
            return is_sq


# Fraction fields
#################

cdef class RingExtensionFractionFieldElement(RingExtensionElement):
    def _repr_(self):
        num = self.numerator()
        denom = self.denominator()
        if denom == 1:
            sd = ""
        elif denom == -1:
            num = -num
            sd = ""
        elif denom._is_atomic():
            sd = "/%s" % denom
        elif (-denom)._is_atomic():
            sd = "/%s" % (-denom)
            num = -num
        if num._is_atomic():
            return "%s%s" % (num, sd)
        else:
            return "(%s)%s" % (num, sd)

    def numerator(self):
        ring = (<RingExtensionFractionField>self._parent)._ring
        try:
            num = self._backend.numerator()
        except AttributeError:
            num = self._backend * self._backend.denominator()
        return ring(num)

    def denominator(self):
        parent = (<RingExtension_class>self._parent)._base
        denom = self._backend.denominator()
        return parent(denom)


# Finite free extensions
########################

cdef class RingExtensionWithBasisElement(RingExtensionElement):
    def __hash__(self):
        return hash(self._backend)

    @cached_method
    def _repr_(self):
        names = (<RingExtensionWithBasis>self._parent)._basis_names
        coeffs = self.vector()
        s = ""
        for i in range(len(names)):
            if coeffs[i].is_zero(): continue
            c = coeffs[i]
            sign = 1
            ss = ""
            if c == 1:
                pass
            elif c == -1:
                sign = -1
            else:
                atomic = c._is_atomic()
                if not atomic and (-c)._is_atomic():
                    c = -c
                    sign = -sign
                    atomic = True
                sc = str(c)
                if atomic:
                    ss += sc
                else:
                    ss += "(" + sc + ")"
                if names[i] != "": ss += "*"
            if ss != "" and ss[0] == "-":
                ss = ss[1:]
                sign *= -1
            if s == "":
                if sign == -1: s = "-"
            else:
                s += " + " if sign == 1 else " - "
            ss += names[i]
            if ss == "": ss += "1"
            s += ss
        if s == "": return "0"
        if s[0] == "(" and s[-1] == ")":
            s = s[1:-1]
        return s

    def vector(self, base=None):
        base = (<RingExtension_class>self._parent)._check_base(base)
        return self._vector(base)

    cdef _vector(self, CommutativeRing base):
        _, _, j = (<RingExtensionWithBasis>self._parent)._free_module(base, map=True)
        return j(self)

    def polynomial(self, base=None, var='x'):
        if base is None:
            base = (<RingExtension_class>self._parent)._base
        degrees = [ ]
        b = self._parent
        degree = 1
        while b is not base:
            if not isinstance(b, RingExtensionWithGen):
                raise NotImplementedError
            reldeg = b.relative_degree()
            degree *= reldeg
            degrees.append(reldeg)
            if b is b.base_ring():
                raise ValueError("not defined over (%s)" % (self, base))
            b = b.base_ring()
        coeffs = { }
        v = self.vector(base)
        S = PolynomialRing(base, len(degrees), names=var)
        for i in range(degree):
            ii = ZZ(i)
            exponents = [ ]
            for j in range(len(degrees)):
                ii, exponent = ii.quo_rem(degrees[j])
                exponents.append(exponent)
            coeffs[tuple(exponents)] = v[i]
        return S(coeffs)

    def matrix(self, base=None):
        cdef RingExtension_class parent = self._parent
        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        return self._matrix(base)

    cdef _matrix(self, CommutativeRing base):
        from sage.matrix.matrix_space import MatrixSpace
        cdef RingExtensionWithBasis parent = self._parent
        _, _, j = parent._free_module(base, map=True)
        x = self._backend
        M = [ j(x * (<RingExtensionElement>b)._backend) for b in parent._basis_over(base) ]
        return MatrixSpace(base, len(M))(M)

    def trace(self, base=None):
        cdef RingExtension_class parent = self._parent
        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        return self._trace(base)

    cdef _trace(self, CommutativeRing base):
        cdef RingExtensionWithBasis parent = self._parent
        cdef CommutativeRing b
        if base is parent:
            return self
        b = parent._base
        t = self._matrix(b).trace()
        if base is b:
            return t
        return (<RingExtensionWithBasisElement>t)._trace(base)

    def norm(self, base=None):
        cdef RingExtension_class parent = self._parent
        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        return self._norm(base)

    cdef _norm(self, CommutativeRing base):
        cdef RingExtensionWithBasis parent = self._parent
        cdef CommutativeRing b
        if base is parent:
            return self
        b = parent._base
        n = self._matrix(b).determinant()
        if base is b:
            return n
        return (<RingExtensionWithBasisElement>n)._norm(base)

    def charpoly(self, base=None, var='x'):
        return self.matrix(base).charpoly(var)

    cpdef minpoly(self, base=None, var='x'):
        from sage.modules.free_module import VectorSpace
        cdef RingExtensionWithBasis parent = self._parent
        cdef MapRelativeFieldToVectorSpace j

        base = parent._check_base(base)
        if not (parent._is_finite_over(base) and parent._is_free_over(base)):
            raise ValueError("the extension is not finite free")
        if not base in Fields():
            raise NotImplementedError("minpoly is only implemented when the base is a field")
        K = backend_parent(base)
        degree = parent._degree_over(base)
        _, _, j = parent._free_module(base, map=True)
        V = VectorSpace(K, degree)
        vector = [K(1)] + (degree-1)*[K(0)]
        vectors = [vector]
        W = V.span(vectors)
        elt = self
        while True:
            vector = V(j.backend_coefficients(elt))
            if vector in W: break
            vectors.append(vector)
            W += V.span([vector])
            elt *= self
        W = V.span_of_basis(vectors)
        coeffs = [ -c for c in W.coordinate_vector(vector) ] + [K(1)]
        return PolynomialRing(base, name=var)(coeffs)

"""
Base class for finite field elements.

AUTHORS::

- David Roe (2010-1-14) -- factored out of sage.structure.element

"""
include "../../ext/stdsage.pxi"

from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.integer import Integer

def is_FiniteFieldElement(x):
    """
    Returns if x is a finite field element.

    EXAMPLE::

        sage: from sage.rings.finite_rings.element_ext_pari import is_FiniteFieldElement
        sage: is_FiniteFieldElement(1)
        False
        sage: is_FiniteFieldElement(IntegerRing())
        False
        sage: is_FiniteFieldElement(GF(5)(2))
        True
    """
    from sage.rings.finite_rings.finite_field_base import is_FiniteField
    return isinstance(x, Element) and is_FiniteField(x.parent())

cdef class FiniteFieldElement(FieldElement):

    def _im_gens_(self, codomain, im_gens):
        """
        Used for applying homomorphisms of finite fields.

        EXAMPLES::

            sage: k.<a> = FiniteField(73^2, 'a')
            sage: K.<b> = FiniteField(73^4, 'b')
            sage: phi = k.hom([ b^(73*73+1) ])
            sage: phi(0)
            0
            sage: phi(a)
            7*b^3 + 13*b^2 + 65*b + 71

            sage: phi(a+3)
            7*b^3 + 13*b^2 + 65*b + 1
        """
        ## NOTE: see the note in sage/rings/number_field_element.pyx,
        ## in the comments for _im_gens_ there -- something analogous
        ## applies here.
        return codomain(self.polynomial()(im_gens[0]))

    def minpoly(self,var='x'):
        """
        Returns the minimal polynomial of this element
        (over the corresponding prime subfield).

        EXAMPLES::

            sage: k.<a> = FiniteField(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: b=a**20;p=b.charpoly("x");p
            x^2 + 15*x + 4
            sage: factor(p)
            (x + 17)^2
            sage: b.minpoly('x')
            x + 17
        """
        p=self.charpoly(var);
        for q in p.factor():
            if q[0](self)==0:
                return q[0]
        # This shouldn't be reached, but you never know!
        raise ArithmeticError("Could not find the minimal polynomial")

        ## We have two names for the same method
        ## for compatibility with sage.matrix
    def minimal_polynomial(self,var='x'):
        """
        Returns the minimal polynomial of this element
        (over the corresponding prime subfield).

        EXAMPLES::

            sage: k.<a> = FiniteField(3^4)
            sage: parent(a)
            Finite Field in a of size 3^4
            sage: b=a**20;p=charpoly(b,"y");p
            y^4 + 2*y^2 + 1
            sage: factor(p)
            (y^2 + 1)^2
            sage: b.minimal_polynomial('y')
            y^2 + 1
        """
        return self.minpoly(var)

    def vector(self, reverse=False):
        r"""
        See :meth:`_vector_`.

        EXAMPLE::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 + 1
            sage: e.vector() # random-ish error message
            doctest:1: DeprecationWarning:The function vector is replaced by _vector_.
            (1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        from sage.misc.misc import deprecation
        deprecation("The function vector is replaced by _vector_.")
        return self._vector_()

    def _vector_(self, reverse=False):
        """
        Return a vector in self.parent().vector_space() matching
        self. The most significant bit is to the right.

        INPUT::

        - ``reverse`` -- reverse the order of the bits
          from little endian to big endian.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: k(v)
            a^2 + 1

            sage: k.<a> = GF(3^16)
            sage: e = 2*a^2 + 1
            sage: v = vector(e)
            sage: v
            (1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: k(v)
            2*a^2 + 1

        You can also compute the vector in the other order::

            sage: e._vector_(reverse=True)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1)
        """
        #vector(foo) might pass in ZZ
        if PY_TYPE_CHECK(reverse, Parent):
            raise TypeError, "Base field is fixed to prime subfield."

        k = self.parent()

        v = self.polynomial().list()

        ret = [v[i] for i in range(len(v))]

        for i in range(k.degree() - len(ret)):
            ret.append(0)

        if reverse:
            ret = list(reversed(ret))
        return k.vector_space()(ret)

    def matrix(self, reverse=False):
        r"""
        See :meth:`_matrix_`.

        EXAMPLE::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 + 1
            sage: e.matrix() # random-ish error message
            doctest:1: DeprecationWarning:The function matrix is replaced by _matrix_.
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
            [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
            [1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0]
            [0 1 0 1 0 0 0 0 0 0 0 0 0 0 1 1]
            [0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 1]
            [0 0 0 1 0 1 0 0 0 0 0 0 0 0 1 0]
            [0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1]
        """
        from sage.misc.misc import deprecation
        deprecation("The function matrix is replaced by _matrix_.")
        return self._matrix_()


    def _matrix_(self, reverse=False):
        """
        Return the matrix of right multiplication by the element on
        the power basis `1, x, x^2, \ldots, x^{d-1}` for the field
        extension.  Thus the \emph{rows} of this matrix give the images
        of each of the `x^i`.

        INPUT:

        - ``reverse`` - if True act on vectors in reversed order

        EXAMPLE::

            sage: k.<a> = GF(2^4)
            sage: a._vector_(reverse=True), a._matrix_(reverse=True) * a._vector_(reverse=True)
            ((0, 0, 1, 0), (0, 1, 0, 0))
            sage: vector(a), matrix(a) * vector(a)
            ((0, 1, 0, 0), (0, 0, 1, 0))
        """
        import sage.matrix.matrix_space

        K = self.parent()
        a = K.gen()
        x = K(1)
        d = K.degree()

        columns = []

        if not reverse:
            l = xrange(d)
        else:
            l = reversed(range(d))

        for i in l:
            columns.append( (self * x)._vector_() )
            x *= a

        k = K.base_ring()
        M = sage.matrix.matrix_space.MatrixSpace(k, d)

        if reverse:
            return M(columns)
        else:
            return M(columns).transpose()
    def _latex_(self):
        r"""
        Return the latex representation of self, which is just the
        latex representation of the polynomial representation of self.

        EXAMPLES::

            sage: k.<b> = GF(5^2); k
            Finite Field in b of size 5^2
            sage: b._latex_()
            'b'
            sage: (b^2+1)._latex_()
            'b + 4'
        """
        if self.parent().degree()>1:
            return self.polynomial()._latex_()
        else:
            return str(self)

    def _pari_init_(self, var=None):
        """
        Return string that when evaluated in PARI defines this element.

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b._pari_init_()
            'Mod(b, Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'
            sage: (2*b+3)._pari_init_()
            'Mod(2*b + 3, Mod(1, 5)*b^2 + Mod(4, 5)*b + Mod(2, 5))'
        """
        g = self.parent()._finite_field_ext_pari_modulus_as_str()
        f = str(self.polynomial())
        s = 'Mod(%s, %s)'%(f, g)
        if var is None:
            return s
        return s.replace(self.parent().variable_name(), var)

    def charpoly(self, var='x', algorithm='matrix'):
        """
        Return the characteristic polynomial of self as a polynomial with given variable.

        INPUT:

        - ``var`` - string (default: 'x')

        - ``algorithm`` - string (default: 'matrix')

          - 'matrix' - return the charpoly computed from the matrix of
            left multiplication by self

          - 'pari' -- use pari's charpoly routine on polymods, which
            is not very good except in small cases

        The result is not cached.

        EXAMPLES::

            sage: k.<a> = GF(19^2)
            sage: parent(a)
            Finite Field in a of size 19^2
            sage: a.charpoly('X')
            X^2 + 18*X + 2
            sage: a^2 + 18*a + 2
            0
            sage: a.charpoly('X', algorithm='pari')
            X^2 + 18*X + 2
        """
        if algorithm == 'matrix':
            return self._matrix_().charpoly(var)
        elif algorithm == 'pari':
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self.parent().prime_subfield(), var)
            return R(self._pari_().charpoly('x').lift())
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm

    def norm(self):
        """
        Return the norm of self down to the prime subfield.

        This is the product of the Galois conjugates of self.

        EXAMPLES::

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: b.norm()
            2
            sage: b.charpoly('t')
            t^2 + 4*t + 2

        Next we consider a cubic extension::

            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: a.norm()
            2
            sage: a.charpoly('t')
            t^3 + 3*t + 3
            sage: a * a^5 * (a^25)
            2
        """
        f = self.charpoly('x')
        n = f[0]
        if f.degree() % 2 != 0:
            return -n
        else:
            return n

    def trace(self):
        """
        Return the trace of this element, which is the sum of the
        Galois conjugates.

        EXAMPLES::

            sage: S.<a> = GF(5^3); S
            Finite Field in a of size 5^3
            sage: a.trace()
            0
            sage: a.charpoly('t')
            t^3 + 3*t + 3
            sage: a + a^5 + a^25
            0
            sage: z = a^2 + a + 1
            sage: z.trace()
            2
            sage: z.charpoly('t')
            t^3 + 3*t^2 + 2*t + 2
            sage: z + z^5 + z^25
            2
        """
        return self.parent().prime_subfield()(self._pari_().trace().lift())

    def multiplicative_order(self):
        """
        Return the multiplicative order of this field element.

        """
        import sage.rings.arith

        if self.is_zero():
            raise ArithmeticError, "Multiplicative order of 0 not defined."
        n = self._parent.order() - 1
        order = 1
        for p, e in sage.rings.arith.factor(n):
            # Determine the power of p that divides the order.
            a = self**(n/(p**e))
            while a != 1:
                order = order * p
                a = a**p

        return order

    def nth_root(self, int n, extend = False, all = False):
        r"""
        Returns an nth root of self.

        INPUT:

        - ``n`` - integer >= 1 (must fit in C int type)

        - ``extend`` - bool (default: True); if True, return an nth
          root in an extension ring, if necessary. Otherwise, raise a
          ValueError if the root is not in the base ring.  Warning:
          this option is not implemented!

        - ``all`` - bool (default: False); if True, return all nth
          roots of self, instead of just one.

        OUTPUT:

        If self has an nth root, returns one (if all == False) or a
        list of all of them (if all == True).  Otherwise, raises a
        ValueError (if extend = False) or a NotImplementedError (if
        extend = True).

        .. warning::

           The 'extend' option is not implemented (yet).

        AUTHOR:

        - David Roe (2007-10-3)

        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.integer import Integer
        if extend:
            raise NotImplementedError
        R = PolynomialRing(self.parent(), "x")
        f = R([-self] + [self.parent()(Integer(0))] * (n - 1) + [self.parent()(1)])
        L = f.roots()
        if all:
            return [x[0] for x in L]
        else:
            if len(L) == 0:
                raise ValueError, "no nth root"
            else:
                return L[0][0]

    def additive_order(self):
        """
        Return the additive order of this finite field element.

        EXAMPLES::

            sage: k.<a> = FiniteField(2^12, 'a')
            sage: b = a^3 + a + 1
            sage: b.additive_order()
            2
            sage: k(0).additive_order()
            1
        """
        if self.is_zero():
            from sage.rings.integer import Integer
            return Integer(1)
        return self.parent().characteristic()

    def pth_power(self, int k = 1):
        """
        Return the `(p^k)^{th}` power of self, where `p` is the
        characteristic of the field.

        INPUT:

        - ``k`` - integer (default: 1, must fit in C int type)

        Note that if `k` is negative, then this computes the appropriate root.

        EXAMPLES::

            sage: F.<a> = GF(29^2)
            sage: z = a^2 + 5*a + 1
            sage: z.pth_power()
            19*a + 20
            sage: z.pth_power(10)
            10*a + 28
            sage: z.pth_power(-10) == z
            True
            sage: F.<b> = GF(2^12)
            sage: y = b^3 + b + 1
            sage: y == (y.pth_power(-3))^(2^3)
            True
            sage: y.pth_power(2)
            b^7 + b^6 + b^5 + b^4 + b^3 + b
        """
        p = self.additive_order()
        n = self.parent().degree()
        return self**(p**(k % n))

    frobenius = pth_power

    def pth_root(self, int k = 1):
        """
        Return the `(p^k)^{th}` root of self, where `p` is the characteristic
        of the field.

        INPUT:

        - ``k`` - integer (default: 1, must fit in C int type)

        Note that if `k` is negative, then this computes the appropriate power.

        EXAMPLES::

            sage: F.<b> = GF(2^12)
            sage: y = b^3 + b + 1
            sage: y == (y.pth_root(3))^(2^3)
            True
            sage: y.pth_root(2)
            b^11 + b^10 + b^9 + b^7 + b^5 + b^4 + b^2 + b
        """
        return self.pth_power(-k)

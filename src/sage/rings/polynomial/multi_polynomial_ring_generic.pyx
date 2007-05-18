include '../../ext/stdsage.pxi'

from sage.structure.parent_gens cimport ParentWithGens
import sage.misc.latex
import multi_polynomial_ideal


def is_MPolynomialRing(x):
    return bool(PY_TYPE_CHECK(x, MPolynomialRing_generic))

cdef class MPolynomialRing_generic(sage.rings.ring.CommutativeRing):
    def __init__(self, base_ring, n, names, order):
        if not isinstance(base_ring, sage.rings.ring.CommutativeRing):
            raise TypeError, "Base ring must be a commutative ring."
        n = int(n)
        if n < 0:
            raise ValueError, "Multivariate Polynomial Rings must " + \
                  "have more than 0 variables."
        self.__ngens = n
        self.__term_order = order
        self._has_singular = False #cannot convert to Singular by default
        ParentWithGens.__init__(self, base_ring, names)

    def is_integral_domain(self):
        return self.base_ring().is_integral_domain()

    cdef _coerce_c_impl(self, x):
        """
        Return the canonical coercion of x to this multivariate
        polynomial ring, if one is defined, or raise a TypeError.

        The rings that canonically coerce to this polynomial ring are:
            * this ring itself
            * polynomial rings in the same variables over any base ring that
              canonically coerces to the base ring of this ring
            * any ring that canonically coerces to the base ring of this
              polynomial ring.

        """
        try:
            P = x.parent()
            # polynomial rings in the same variable over the any base that coerces in:
            if is_MPolynomialRing(P):
                if P.variable_names() == self.variable_names():
                    if self.has_coerce_map_from(P.base_ring()):
                        return self(x)
                    else:
                        raise TypeError, "no natural map between bases of polynomial rings"
                else:
                    raise TypeError, "polynomial ring has a different indeterminate names"

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this polynomial ring.
        return self._coerce_try(x, [self.base_ring()])

    def __richcmp__(left, right, int op):
        return (<ParentWithGens>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
        if not is_MPolynomialRing(right):
            return cmp(type(left),type(right))
        else:
            return cmp((left.base_ring(), left.__ngens, left.variable_names(), left.__term_order),
                       (right.base_ring(), (<MPolynomialRing_generic>right).__ngens, right.variable_names(), (<MPolynomialRing_generic>right).__term_order))

    def __contains__(self, x):
        """
        This definition of containment does not involve a natural
        inclusion from rings with less variables into rings with more.
        """
        try:
            return x.parent() == self
        except AttributeError:
            return False

    def _repr_(self):
        return "Polynomial Ring in %s over %s"%(", ".join(self.variable_names()), self.base_ring())

    def _latex_(self):
        vars = str(self.latex_variable_names()).replace('\n','').replace("'",'')
        return "%s[%s]"%(sage.misc.latex.latex(self.base_ring()), vars[1:-1])


    def _ideal_class_(self):
        return multi_polynomial_ideal.MPolynomialIdeal

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            # all that is needed is that elements of the base ring
            # of the polynomial ring canonically coerce into codomain.
            # Since poly rings are free, any image of the gen
            # determines a homomorphism
            codomain._coerce_(self.base_ring()(1))
        except TypeError:
            return False
        return True

    def _magma_(self, magma=None):
        """
        Used in converting this ring to the corresponding ring in MAGMA.

        EXAMPLES:
            sage: R.<y,z,w> = PolynomialRing(QQ,3)
            sage: magma(R) # optional
            Polynomial ring of rank 3 over Rational Field
            Graded Reverse Lexicographical Order
            Variables: y, z, w

            sage: magma(PolynomialRing(GF(7),4, 'x')) #optional
            Polynomial ring of rank 4 over GF(7)
            Graded Reverse Lexicographical Order
            Variables: x0, x1, x2, x3

            sage: magma(PolynomialRing(GF(49,'a'),10, 'x')) #optional
            Polynomial ring of rank 10 over GF(7^2)
            Graded Reverse Lexicographical Order
            Variables: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9

            sage: magma(PolynomialRing(ZZ['a,b,c'],3, 'x')) #optional
            Polynomial ring of rank 3 over Polynomial ring of rank 3 over Integer Ring
            Graded Reverse Lexicographical Order
            Variables: x0, x1, x2
        """
        if magma == None:
            import sage.interfaces.magma
            magma = sage.interfaces.magma.magma

        try:
            if self.__magma is None:
                raise AttributeError
            m = self.__magma
            m._check_valid()
            if not m.parent() is magma:
                raise ValueError
            return m
        except (AttributeError,ValueError):
            B = magma(self.base_ring())
            R = magma('PolynomialRing(%s, %s, %s)'%(B.name(), self.ngens(),self.term_order().magma_str()))
            R.assign_names(self.variable_names())
            self.__magma = R
            return R

    def var_dict(self):
        """
        Return dictionary of paris varname:var of the variables
        of this multivariate polynomial ring.
        """
        return dict([(str(g),g) for g in self.gens()])


    def is_finite(self):
        if self.ngens() == 0:
            return self.base_ring().is_finite()
        return False

    def is_field(self):
        """
        Return True if this multivariate polynomial ring is a field, i.e.,
        it is a ring in 0 generators over a field.
        """
        if self.ngens() == 0:
            return self.base_ring().is_field()
        return False

    def term_order(self):
        return self.__term_order

    def characteristic(self):
        """
        Return the characteristic of this polynomial ring.

        EXAMPLES:
            sage: R = MPolynomialRing(QQ, 'x', 3)
            sage: R.characteristic()
            0
            sage: R = MPolynomialRing(GF(7),'x', 20)
            sage: R.characteristic()
            7
        """
        return self.base_ring().characteristic()

    def gen(self, n=0):
        if n < 0 or n >= self.__ngens:
            raise ValueError, "Generator not defined."
        return self._gens[int(n)]

    def gens(self):
        return self._gens

    def krull_dimension(self):
        return self.base_ring().krull_dimension() + self.ngens()

    def ngens(self):
        return self.__ngens

    def _monomial_order_function(self):
        raise NotImplementedError

    def latex_variable_names(self):
        """
        Returns the list of variable names suitable for latex output.

        All '_SOMETHING' substrings are replaced by '_{SOMETHING}' recursively
        so that subscripts of subscripts work.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ,'x',12).objgens()
            sage: x
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)
            sage: print R.latex_variable_names ()
            ['x_{0}', 'x_{1}', 'x_{2}', 'x_{3}', 'x_{4}', 'x_{5}', 'x_{6}', 'x_{7}', 'x_{8}', 'x_{9}', 'x_{10}', 'x_{11}']
            sage: f = x[0]^3 + 15/3 * x[1]^10
            sage: print latex(f)
            5 x_{1}^{10} + x_{0}^{3}
        """
        if self._latex_names is not None:
            return self._latex_names
        names = []
        for g in self.variable_names():
            i = len(g)-1
            while i >= 0 and g[i].isdigit():
                i -= 1
            if i < len(g)-1:
                g = '%s_{%s}'%(g[:i+1], g[i+1:])
            names.append(g)
        self._latex_names = names
        return names

    def __reduce__(self):
        """
        """

        base_ring = self.base_ring()
        n = self.ngens()
        names = self.variable_names()
        order = self.term_order()

        return unpickle_MPolynomialRing_generic_v1,(base_ring, n, names, order)


####################
# Leave *all* old versions!

def unpickle_MPolynomialRing_generic_v1(base_ring, n, names, order):
    from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing
    return MPolynomialRing(base_ring, n, names=names, order=order)


def unpickle_MPolynomialRing_generic(base_ring, n, names, order):
    from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing

    return MPolynomialRing(base_ring, n, names=names, order=order)

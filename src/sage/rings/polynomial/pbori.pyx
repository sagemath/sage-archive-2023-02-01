r"""
Boolean Polynomials.

Elements of the quotient ring

    $$\F_2[x_1,...,x_n]/<x_1^2+x_1,...,x_n^2+x_n>.$$

are called boolean polynomials. Boolean polynomials arise naturally in
cryptography, coding theory, formal logic, chip design and other
areas. This implementation is a thin wrapper around the \PolyBoRi
library by Michael Brickenstein and Alexander Dreyer.

``Boolean polynomials can be modelled in a rather simple way, with
both coefficients and degree per variable lying in $\{0, 1\}$. The
ring of Boolean polynomials is, however, not a polynomial ring, but
rather the quotient ring of the polynomial ring over the field with
two elements modulo the field equations $x^2=x$ for each
variable $x$. Therefore, the usual polynomial data structures seem not
to be appropriate for fast Groebner basis computations.  We introduce
a specialised data structure for Boolean polynomials based on
zero-suppressed binary decision diagrams (ZDDs), which is capable of
handling these polynomials more efficiently with respect to memory
consumption and also computational speed.  Furthermore, we concentrate
on high-level algorithmic aspects, taking into account the new data
structures as well as structural properties of Boolean
polynomials.'' -- [BD07]

AUTHORS:
    -- Michael Brickenstein: \PolyBoRi author
    -- Alexander Dreyer: \PolyBoRi author
    -- Burcin Erocal <burcin@erocal.org>: main \SAGE wrapper author
    -- Martin Albrecht <malb@informatik.uni-bremen.de>: some contributions to the \SAGE wrapper

EXAMPLES:

Consider the ideal

  $$<ab + cd + 1, ace + de, abe + ce, bc + cde + 1>.$$

First, we compute the lexicographical Groebner basis in the polynomial
ring $$R = \F_2[a,b,c,d,e].$$

    sage: P.<a,b,c,d,e> = PolynomialRing(GF(2), 5, order='lex')
    sage: I1 = ideal([a*b + c*d + 1, a*c*e + d*e, a*b*e + c*e, b*c + c*d*e + 1])
    sage: for f in I1.groebner_basis():
    ...     f
    d^4*e^2 + d^4*e + d^3*e + d^2*e^2 + d^2*e + d*e + e
    c*e + d^3*e^2 + d^3*e + d^2*e^2 + d*e
    b*e + d*e^2 + d*e + e
    b*c + d^3*e^2 + d^3*e + d^2*e^2 + d*e + e + 1
    a + c^2*d + c + d^2*e

If one wants to solve this system over the algebraic closure of $\F_2$
then this Groebner basis was the one to consider. If one wants
solutions over $\F_2$ only then one adds the field polynomials to the
ideal to force the solutions in $\F_2$.

    sage: J = I1 + sage.rings.ideal.FieldIdeal(P)
    sage: for f in J.groebner_basis():
    ...     f
    e
    d^2 + d
    c + 1
    b + 1
    a + d + 1

So the solutions over $\F_2$ are $\{e=0, d=1, c=1, b=1, a=0\}$ and $\{e=0,
d=0, c=1, b=1, a=1\}$.

We can express the restriction to $\F_2$ by considering the quotient
ring. If $I$ is an ideal in $\F[x_1, ..., x_n]$ then the ideals in the
quotient ring $\F[x_1, ..., x_n]/I$ are in one-to-one correspondence
with the ideals of $\F[x_0, ..., x_n]$ containing $I$ (that is, the
ideals $J$ satisfying $I \subset J \subset P$).

    sage: Q = P.quotient( sage.rings.ideal.FieldIdeal(P) )
    sage: I2 = ideal([Q(f) for f in I1.gens()])
    sage: for f in I2.groebner_basis():
    ...     f
    ebar
    cbar + 1
    bbar + 1
    abar + dbar + 1

This quotient ring is exactly what \PolyBoRi handles well.

    sage: B.<a,b,c,d,e> = BooleanPolynomialRing(5, order='lex')
    sage: I2 = ideal([B(f) for f in I1.gens()])
    sage: for f in I2.groebner_basis():
    ...     f
    a + d + 1
    b + 1
    c + 1
    e

Note that $d^2 + d$ is not representable in $B == Q$. Also note, that
\PolyBoRi cannot play out its strength in such small examples,
i.e. working in the polynomial ring might be faster for small examples
like this.

\subsection{Implementation specific notes}

\PolyBoRi comes with a Python wrapper. However this wrapper does not
match \SAGE's style and is written using Boost. Thus \SAGE's wrapper
is a reimplementation of Python bindings to \PolyBoRi's C++
library. This interface is written in Cython like all of \SAGE's C/C++
library interfaces. An interface in \PolyBoRi style is also provided
which is effectively a reimplementation of the official Boost wrapper
in Cython. This means that some functionality of the official wrapper
might be missing from this wrapper and this wrapper might have bugs
not present in the offical Python interface.

\subsection{Access to the original \PolyBoRi interface}

The re-implementation \PolyBoRi's native wrapper is available to the user too:

    sage: from polybori import *
    sage: declare_ring([Block('x',2),Block('y',3)],globals())
    Boolean PolynomialRing in x(0), x(1), y(0), y(1), y(2)
    sage: r
    Boolean PolynomialRing in x(0), x(1), y(0), y(1), y(2)

    sage: [Variable(i) for i in xrange(r.ngens())]
    [x(0), x(1), y(0), y(1), y(2)]

For details on this interface see:

   \url{http://polybori.sourceforge.net/doc/tutorial/tutorial.html}.


REFERENCES:
    [BD07] Michael Brickenstein, Alexander Dreyer; '\PolyBoRi: A
    Groebner basis framework for Boolean polynomials';
    http://www.itwm.fraunhofer.de/zentral/download/berichte/bericht122.pdf
"""

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
include '../../libs/polybori/decl.pxi'

from sage.rings.integer import Integer
from sage.rings.finite_field import FiniteField as GF

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.term_order import TermOrder

from sage.structure.element cimport Element
from sage.structure.element cimport RingElement
from sage.structure.element cimport ModuleElement

from sage.structure.parent cimport Parent
from sage.structure.sequence import Sequence

from sage.monoids.monoid import Monoid_class

order_dict= {"lp":      lp,
             "dlex":    dlex,
             "dp_asc":  dp_asc,
             "block_dlex":   block_dlex,
             "block_dp_asc": block_dp_asc,
             }

order_mapping = {'lp':   lp,
                 'Dp':   dlex,
                 'dp':   dp_asc}

cdef BooleanPolynomialRing cur_ring

cdef class BooleanPolynomialRing(MPolynomialRing_generic):
    def __init__(self, n, names, order='lex'):
        """
        Construct a boolean polynomial ring with the following
        parameters:

        INPUT:
            n -- number of variables (an integer > 1)
            names -- names of ring variables, may be a string of
                     list/tuple
            order -- term order (default: lex)

        EXAMPLES:
            sage: R.<x, y, z> = BooleanPolynomialRing(3)
            sage: R
            Boolean PolynomialRing in x, y, z

            sage: p = x*y + x*z + y*z
            sage: x*p
            x*y*z + x*y + x*z

            sage: R.term_order()
            Lexicographic term order

            sage: R = BooleanPolynomialRing(5,'x',order='deglex(3),deglex(2)')
            sage: R.term_order()
            deglex(3),deglex(2) term order

            sage: R = BooleanPolynomialRing(3,'x',order='degrevlex')
            sage: R.term_order()
            Degree reverse lexicographic term order

        TESTS:
            sage: P.<x,y> = BooleanPolynomialRing(2,order='degrevlex')
            sage: x > y
            True

            sage: P.<x0, x1, x2, x3> = BooleanPolynomialRing(4,order='degrevlex(2),degrevlex(2)')
            sage: x0 > x1
            True
            sage: x2 > x3
            True
        """
        cdef Py_ssize_t i, j, bstart, bsize
        try:
            n = int(n)
        except TypeError, msg:
            raise TypeError, "Number of variables must be an integer"

        if n < 1:
            raise ValueError, "Number of variables must be greater than 1."

        self.pbind = <Py_ssize_t*>sage_malloc(n*sizeof(Py_ssize_t))
        cdef char *_n

        order = TermOrder(order, n)

        try:
            pb_order_code = order_mapping[order.blocks[0][0]]
        except KeyError:
            raise ValueError, "Only lex, deglex, degrevlex orders are supported."

        if len(order.blocks) > 1:
            if pb_order_code is lp:
                raise ValueError, "Only deglex and degrevlex are supported for block orders."
            elif pb_order_code is dlex:
                pb_order_code = block_dlex
            elif pb_order_code is dp_asc:
                pb_order_code = block_dp_asc
            for i in range(1, len(order.blocks)):
                if order.blocks[0][0] != order.blocks[i][0]:
                    raise ValueError, "Each block must have the same order type (deglex or degrevlex) for block orderings."

        if (pb_order_code is dlex) or (pb_order_code is lp) or \
                (pb_order_code is block_dlex):
            for i from 0 <= i < n:
                self.pbind[i] = i
        elif pb_order_code is dp_asc:
            for i from 0 <= i < n:
                self.pbind[i] = n - i -1
        else:
            # pb_order_code is block_dp_asc:
            bstart = 0
            for i from 0 <= i < len(order.blocks):
                bsize = order.blocks[i][1]
                for j from 0 <= j < bsize:
                    self.pbind[bstart + j] = bstart + bsize - j -1
                bstart += bsize

        PBRing_construct(&self._pbring, n, pb_order_code)

        MPolynomialRing_generic.__init__(self, GF(2), n, names, order)

        counter = 0
        for i in range(len(order.blocks)-1):
            counter += order.blocks[i][1]
            pb_append_block(counter)

        for i from 0 <= i < n:
            _n = self._names[self.pbind[i]]
            pb_set_variable_name(i, _n)

        self._zero_element = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._zero_element)._pbpoly, 0)
        self._one_element  = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._one_element)._pbpoly, 1)

        self._monom_monoid = BooleanMonomialMonoid(self)

        global cur_ring
        cur_ring = self

    def __dealloc__(self):
        sage_free(self.pbind)
        PBRing_destruct(&self._pbring)

    def __reduce__(self):
        """
        EXAMPLE:
            sage: P.<a,b> = BooleanPolynomialRing(2)
            sage: loads(dumps(P)) == P # indirect doctest
            True
        """
        n = self.ngens()
        names = self.variable_names()
        order = self.term_order()
        return unpickle_BooleanPolynomialRing,(n, names, order)

    def ngens(self):
        """
        Returns the number of variables in self.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: P.ngens()
            2

            sage: P = BooleanPolynomialRing(1000, 'x')
            sage: P.ngens()
            1000
        """
        return self._pbring.nVariables()

    def gen(self, int i=0):
        """
        Returns the i-th generator of self.

        INPUT:
            i -- an integer

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.gen()
            x
            sage: P.gen(2)
            z

        TESTS:
            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='dp')
            sage: P.gen(0)
            x
        """
        if i < 0 or i >= self._pbring.nVariables():
            raise ValueError, "Generator not defined."
        return new_BP_from_DD(self, self._pbring.variable(self.pbind[i]))

    def gens(self):
        """
        Return the tuple of variables in this ring.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.gens()
            (x, y, z)

            sage: P = BooleanPolynomialRing(10,'x')
            sage: P.gens()
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)

        TESTS:
            sage: P.<x,y,z> = BooleanPolynomialRing(3,order='degrevlex')
            sage: P.gens()
            (x, y, z)
        """
        return tuple([new_BP_from_DD(self,
            self._pbring.variable(self.pbind[i])) \
                for i from 0<= i < self.__ngens])

    def _repr_(self):
        """
        EXAMPLE:
            sage: P.<x, y> = BooleanPolynomialRing(2)
            sage: P # indirect doctest
            Boolean PolynomialRing in x, y
        """
        gens = ", ".join(map(str,self.gens()))
        return "Boolean PolynomialRing in %s"%(gens)

    cdef _coerce_c_impl(self, other):
        """
        Canonical conversion of elements from other objects to self.

        EXAMPLES:

        Coerce elements of self.

            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: p = x*y + x
            sage: P._coerce_(p)
            x*y + x

        Coerce from monomials over the same ring.

            sage: P._coerce_(p.lm())
            x*y

        Coerce from a different BooleanPolynomialRing.

            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: R = BooleanPolynomialRing(2,'y,x')
            sage: p = R._coerce_(x+y+x*y+1)
            sage: p.parent()
            Boolean PolynomialRing in y, x
            sage: p
            y*x + y + x + 1

        Coerce from polynomials over the integers.

            sage: P = BooleanPolynomialRing(3,'x,y,z')
            sage: R.<z,x,y> = ZZ['z,x,y']
            sage: t = x^2*z+5*y^3
            sage: p = P._coerce_(t)
            sage: p.parent()
            Boolean PolynomialRing in x, y, z
            sage: p
            x*z + y

        Coerce from integers.
            sage: P = BooleanPolynomialRing(3,'x,y,z')
            sage: p = P._coerce_(1)
            sage: p.is_one()
            True
            sage: p = P._coerce_(6)
            sage: p.is_zero()
            True

        Coerce from GF(2).
            sage: P = BooleanPolynomialRing(3,'x,y,z')
            sage: F = GF(2)
            sage: p = P._coerce_(F.zero_element())
            sage: p.is_zero()
            True
            sage: p = P._coerce_(F.one_element())
            sage: p.is_one()
            True

        Coerce from BooleanMonomials over a different BooleanPolynomialRing.
            sage: R.<y,x> = BooleanPolynomialRing(2)
            sage: M = R._monom_monoid
            sage: P = BooleanPolynomialRing(3,'x,y,z')
            sage: t = P._coerce_(M(x*y))
            sage: t
            x*y
            sage: t.parent()
            Boolean PolynomialRing in x, y, z

        TESTS:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: R = BooleanPolynomialRing(1,'y')
            sage: p = R._coerce_(x+y+x*y+1)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from <type 'sage.rings.polynomial.pbori.BooleanPolynomial'> to Boolean PolynomialRing in y

            sage: P = BooleanPolynomialRing(2,'x,y')
            sage: R.<z,x,y> = ZZ['z,x,y']
            sage: t = x^2*z+5*y^3
            sage: p = P._coerce_(t)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'> to Boolean PolynomialRing in x, y

        Test coercion from a ring that compares equal.
            sage: P = BooleanPolynomialRing(2,'x,y')
            sage: R.<x,y> = BooleanPolynomialRing(2)
            sage: P == R
            True
            sage: P(x)
            x
        """
        cdef BooleanPolynomial p
        # we check for other PolyBoRi types first since this conversion
        # is used by the PolyBoRi python code often
        if PY_TYPE_CHECK(other, DD):
            return new_BP_from_DD(self, (<DD>other)._pbdd)
        elif PY_TYPE_CHECK(other, BooleSet):
            return new_BP_from_PBSet(self, (<BooleSet>other)._pbset)

        elif PY_TYPE_CHECK(other, int) or PY_TYPE_CHECK(other, Integer):
            if other %2:
                return self._one_element
            else:
                return self._zero_element
        elif PY_TYPE_CHECK(other, BooleanMonomial):
            if (<BooleanMonomial>other).ring is self:
                p = new_BP_from_PBMonom(self, (<BooleanMonomial>other)._pbmonom)
                return p
            elif (<BooleanMonomial>other)._parent.ngens() <= \
                    self._pbring.nVariables():
                try:
                    var_mapping = get_var_mapping(self, other.parent())
                except NameError, msg:
                    raise ValueError, "cannot coerce monomial %s to %s: %s"%(other,self,msg)
                p = self._one_element
                for i in other.iterindex():
                    p *= var_mapping[i]
                return p
            else:
                raise ValueError, "cannot coerce monomial %s to %s: %s"%(other,self,msg)
        elif PY_TYPE_CHECK(other,BooleanPolynomial) and \
            ((<BooleanPolynomialRing>(<BooleanPolynomial>other)\
            ._parent)._pbring.nVariables() <= self._pbring.nVariables()):
                    try:
                        var_mapping = get_var_mapping(self, other.parent())
                    except NameError, msg:
                        raise ValueError, "cannot coerce polynomial %s to %s: %s"%(other,self,msg)
                    p = self._zero_element
                    for monom in other:
                        new_monom = self._monom_monoid._one_element
                        for i in monom.iterindex():
                            new_monom *= var_mapping[i]
                        p += new_monom
                    return p
        elif (PY_TYPE_CHECK(other, MPolynomial) or \
                PY_TYPE_CHECK(other, Polynomial)) and \
                self.base_ring().has_coerce_map_from(other.base_ring()) and \
                (other.parent().ngens() <= self._pbring.nVariables()):
                    try:
                        var_mapping = get_var_mapping(self, other.parent())
                    except NameError, msg:
                        raise ValueError, "cannot coerce polynomial %s to %s: %s"%(other,self,msg)
                    p = self._zero_element
                    exponents = other.exponents()
                    coefs = other.coefficients()
                    for i in range(len(coefs)):
                        if self._base._coerce_(coefs[i]).is_one():
                            m = self._monom_monoid._one_element
                            for j in range(len(exponents[i])):
                                if exponents[i][j] > 0:
                                    m *= var_mapping[j]
                            p += m
                    return p
        elif PY_TYPE_CHECK(other, Element) and \
                self.base_ring().has_coerce_map_from(other.parent()):
                    if self.base_ring()(other).is_zero():
                        return self._zero_element
                    else:
                        return self._one_element
        else:
            raise TypeError, "cannot coerce from %s to %s" % \
                    (type(other), str(self))

    def __call__(self, other):
        """
        Convert elements of other objects to this boolean polynomial
        ring.

        EXAMPLE:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: P(5)
            1

            sage: P(x+y)
            x + y

            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: R = BooleanPolynomialRing(1,'y')
            sage: p = R(y); p
            y
            sage: p.parent()
            Boolean PolynomialRing in y

            sage: P = BooleanPolynomialRing(2,'x,y')
            sage: R.<z,x,y> = ZZ['z,x,y']
            sage: t = x^2*y + 5*y^3
            sage: p = P(t); p
            x*y + y
            sage: p.parent()
            Boolean PolynomialRing in x, y

        TESTS:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: R = BooleanPolynomialRing(1,'y')
            sage: p = R(x+y+x*y+1)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert polynomial x*y + x + y + 1 to Boolean PolynomialRing in y: name x not defined

            sage: P = BooleanPolynomialRing(2,'x,y')
            sage: R.<z,x,y> = ZZ['z,x,y']
            sage: t = x^2*z+5*y^3
            sage: p = P(t)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert polynomial z*x^2 + 5*y^3 to Boolean PolynomialRing in x, y: name z not defined
        """
        cdef int i

        try:
            return self._coerce_c(other)
        except NameError, msg:
            raise NameError, msg
        except TypeError:
            pass

        if PY_TYPE_CHECK(other, BooleanMonomial) and \
            ((<BooleanMonomial>other)._pbmonom.deg() <= self._pbring.nVariables()):
                try:
                    var_mapping = get_var_mapping(self, other)
                except NameError, msg:
                    raise ValueError, "cannot convert monomial %s to %s: %s"%(other,self,msg)
                p = self._one_element
                for i in other.iterindex():
                    p *= var_mapping[i]
                return p
        elif PY_TYPE_CHECK(other,BooleanPolynomial) and \
                ((<BooleanPolynomial>other)._pbpoly.nUsedVariables() <= \
                self._pbring.nVariables()):
                    try:
                        var_mapping = get_var_mapping(self, other)
                    except NameError, msg:
                        raise ValueError, "cannot convert polynomial %s to %s: %s"%(other,self,msg)
                    p = self._zero_element
                    for monom in other:
                        new_monom = self._monom_monoid._one_element
                        for i in monom.iterindex():
                            new_monom *= var_mapping[i]
                        p += new_monom
                    return p
        elif (PY_TYPE_CHECK(other, MPolynomial) or \
                PY_TYPE_CHECK(other, Polynomial)) and \
                self.base_ring().has_coerce_map_from(other.base_ring()):
                    try:
                        var_mapping = get_var_mapping(self, other)
                    except NameError, msg:
                        raise ValueError, "cannot convert polynomial %s to %s: %s"%(other,self,msg)
                    p = self._zero_element
                    exponents = other.exponents()
                    coefs = other.coefficients()
                    for i in range(len(coefs)):
                        if self._base._coerce_(coefs[i]).is_one():
                            m = self._monom_monoid._one_element
                            for j in range(len(exponents[i])):
                                if exponents[i][j] > 0:
                                    m *= var_mapping[j]
                            p += m
                    return p

        elif PY_TYPE_CHECK(other, str):
            other = other.replace("^","**")
            p = self(eval(other, self.gens_dict(), {}))
            return p

        try:
            i = int(other)
        except:
            raise TypeError, "cannot convert %s to BooleanPolynomial"%(type(other))

        i = i % 2
        if i:
            return self._one_element
        else:
            return self._zero_element

    def __richcmp__(left, right, int op):
        """
        BooleanPolynomialRings are equal if they have the same
            - number of variables
            - variable names
            - order type

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: R.<x,y> = BooleanPolynomialRing(2)
            sage: P == R
            True

            sage: Q.<x,z> = BooleanPolynomialRing(2)
            sage: P == Q
            False

            sage: S.<x,y> = BooleanPolynomialRing(2, order='deglex')
            sage: P == S
            False
        """
        return (<Parent>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
        r"""
        See \code{self.__richcmp__}
        """
        if PY_TYPE_CHECK(right, BooleanPolynomialRing):
            return cmp( (type(left), map(str, left.gens()), left.term_order()),
                    (type(right), map(str, right.gens()), right.term_order()))
        else:
            return -1

    def __hash__(self):
        """
        Return a hash of this boolean polynomial ring.

        EXAMPLE:
            sage: P.<a,b,c,d> = BooleanPolynomialRing(4, order='lex')
            sage: P
            Boolean PolynomialRing in a, b, c, d
            sage: {P:1} # indirect doctest
            {Boolean PolynomialRing in a, b, c, d: 1}
        """
        return hash(str(self))

    def ideal(self, *gens, **kwds):
        """
        Create an ideal in this ring.

        INPUT:
            gens -- list or tuple of generators
            coerce -- bool (default: True) automatically coerce the
                      given polynomials to this ring to form the ideal

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.ideal(x+y)
            Ideal (x + y) of Boolean PolynomialRing in x, y, z

            sage: P.ideal(x*y, y*z)
            Ideal (x*y, y*z) of Boolean PolynomialRing in x, y, z

            sage: P.ideal([x+y, z])
            Ideal (x + y, z) of Boolean PolynomialRing in x, y, z
        """
        from sage.misc.flatten import flatten
        coerce = kwds.get('coerce', True)
        gens = flatten(gens)
        return BooleanPolynomialIdeal(self, gens, coerce)

    def random_element(self, degree=2, terms=5, choose_degree=True,
                       vars_set=None, seed=None):
        """
        Return a random boolean polynomial. Generated polynomial has
        the given number of terms, and at most given degree.

        INPUT:
            degree -- maximum degree (default: 2)
            terms -- number of terms (default: 5)
            choose_degree -- choose degree of monomials randomly first, rather
                             than monomials uniformly random
            vars_set -- list of integer indicies of generators of self to use
                        in the generated polynomial

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.random_element(degree=3, terms=4) # random output
            x*y + x*z + z + 1

            sage: P.random_element(degree=1, terms=2) # random output
            x + 1

        TESTS:
            sage: P.random_element(degree=4)
            Traceback (most recent call last):
            ...
            ValueError: Given degree should be less than or equal to number of variables (3)

            sage: t = P.random_element(degree=1, terms=5)
            Traceback (most recent call last):
            ...
            ValueError: Cannot generate random polynomial with 5 terms and maximum degree 1 using 3 variables

            sage: t = P.random_element(degree=2,terms=5,vars_set=(0,1))
            Traceback (most recent call last):
            ...
            ValueError: Cannot generate random polynomial with 5 terms using 2 variables

        """
        from sage.rings.integer import Integer
        from sage.rings.arith import binomial

        if not vars_set:
            vars_set=range(self.ngens())
        nvars = len(vars_set)

        if degree > nvars:
            raise ValueError, "Given degree should be less than or equal to number of variables (%s)"%(nvars)

        if Integer(terms-1).bits() > nvars:
            raise ValueError, "Cannot generate random polynomial with %s terms using %s variables"%(terms, nvars)

        tot_terms=0
        monom_counts = []
        for i from 0 <= i <= degree:
            tot_terms += binomial(nvars,i)
            monom_counts.append(tot_terms)

        if terms > tot_terms:
            raise ValueError, "Cannot generate random polynomial with %s terms and maximum degree %s using %s variables"%(terms, degree, nvars)

        p = self._zero_element
        while len(p) < terms:
            p=self(p.set().union(\
                self._random_uniform_rec(degree, monom_counts, vars_set, choose_degree, terms-len(p))\
                .set()))
        return p

    def _random_uniform_rec(self, degree, monom_counts, vars_set, dfirst, l):
        r"""
        Recursively generate a random polynomial in in this ring,
        using the variables from \code{vars_set}.

        INPUT:
            degree -- maximum degree
            monom_counts -- a list containing total number of
                            monomials up to given degree
            vars_set -- list of variable indicies to use in the
                        generated polynomial
            dfirst -- if \code{True} choose degree first, otherwise
                      choose the monomial uniformly
            l -- number of monomials to generate

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P._random_uniform_rec(2, [1, 3, 4], (0,1), True, 2) # random
            y + 1
        """
        from sage.rings.integer import Integer
        from sage.rings.integer_ring import ZZ
        if l == 0:
            return self._zero_element
        if l == 1:
            if dfirst:
                return self._random_monomial_dfirst(degree, vars_set)
            else:
                return self._random_monomial_uniform(monom_counts, vars_set)

        return self._random_uniform_rec(degree, monom_counts,
                    vars_set, dfirst, l//2) + \
               self._random_uniform_rec(degree, monom_counts,
                    vars_set, dfirst, l - l//2)

    def _random_monomial_uniform(self, monom_counts, vars_set):
        r"""
        Choose a random monomial uniformly from set of monomials in
        the variables indexed by \code{vars_set} in self.

        INPUT:
            monom_counts -- list of number of monomials up to given
                            degree
            vars_set -- list of variable indicies to use in the
                        generated monomial

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P._random_monomial_uniform([1, 3, 4], (0,1)) # random output
            x
        """
        from sage.rings.integer_ring import ZZ
        from sage.combinat.choose_nk import from_rank

        t = ZZ.random_element(0,monom_counts[-1])
        if t == 0:
            return self._one_element
        i = 1
        while t >= monom_counts[i]:
            i+=1
        mind = t-monom_counts[i-1]
        var_inds = from_rank(mind,len(vars_set),i)
        M = self._monom_monoid
        m = M._one_element
        for j in var_inds:
            m*=M.gen(vars_set[j])
        return self(m)

    def _random_monomial_dfirst(self, degree, vars_set):
        r"""
        Choose a random monomial using variables indexed in
        \code{vars_set} up to given \code{degree}. The degree of the
        monomial, $d$, is chosen uniformly in the interval [0,degree]
        first, then the monomial is generated by selecting a random
        sample of size $d$ from \code{vars_set}.

        INPUT:
            degree -- maximum degree
            vars_set -- list of variable indicies of self

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P._random_monomial_dfirst(3, (0,1,2)) # random output
            x*y
        """
        from sage.rings.integer_ring import ZZ
        from random import sample
        d = ZZ.random_element(0,degree+1)
        vars = sample(vars_set, d)
        M = self._monom_monoid
        m = M._one_element
        for j in vars:
            m*=M.gen(j)
        return self(m)

###
#
# Methods for compatibility with PolyBoRi
#
###

    def _change_ordering(self, int order):
        r"""
        Change the ordering of this boolean polynomial ring. Do NOT
        call this method, unless you know very well what you are
        doing.

        INPUT:
            order -- an integer (0 <= order <= 4)

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3,order='deglex')
            sage: y*z > x
            True

        Now we call the internal method and change the ordering to 'lex':

            sage: B._change_ordering(0)
            sage: y*z > x
            False

        However, this change is not -- and should not be -- picked up
        by the public interface.

            sage: B.term_order()
            Degree lexicographic term order

        WARNING: Do not use this method. It is provided for
        compatibility reasons with \PolyBoRi but parents are supposed
        to be immutable in \Sage.
        """
        if order < 0 or order > 4:
            raise ValueError, "order value %s is not supported"%(order)
        pbenv_changeOrdering(<ordercodes>order)


    def _set_variable_name(self, i, s):
        r"""
        Set variable name of i-th variable to s.

        This function is used by \PolyBoRi python functions.

        INPUT:
            i -- index of variable
            s -- new variable name

        EXAMPLES:
            sage: P.<x0,x1> = BooleanPolynomialRing(2)
            sage: P
            Boolean PolynomialRing in x0, x1

            sage: P._set_variable_name(0, 't')
            sage: P
            Boolean PolynomialRing in t, x1

        WARNING: Do not use this method. It is provided for
        compatibility reasons with \PolyBoRi but parents are supposed
        to be immutable in \Sage.
        """
        self._pbring.activate()
        pb_set_variable_name(i, s)
        t = list(self._names)
        t[i] = s
        self._names = tuple(t)

def get_var_mapping(ring, other):
    """
    Return a variable mapping between variables of \var{other} and
    \var{ring]. When other is a parent object, the mapping defines
    images for all variables of other. If it is an element, only
    variables occuring in other are mapped.

    Raises \code{NameError} if no such mapping is possible.

    EXAMPLES:
        sage: P.<x,y,z> = BooleanPolynomialRing(3)
        sage: R.<z,y> = QQ[]
        sage: sage.rings.polynomial.pbori.get_var_mapping(P,R)
        [z, y]
        sage: sage.rings.polynomial.pbori.get_var_mapping(P, z^2)
        [z, None]

        sage: R.<z,x> = BooleanPolynomialRing(2)
        sage: sage.rings.polynomial.pbori.get_var_mapping(P,R)
        [z, x]
        sage: sage.rings.polynomial.pbori.get_var_mapping(P, x^2)
        [None, x]
    """
    my_names = list(ring._names) # we need .index(.)
    if PY_TYPE_CHECK(other, ParentWithGens):
        vars = range(other.ngens())
        ovar_names = other._names
    else:
        ovar_names = other.parent().variable_names()
        if PY_TYPE_CHECK(other, BooleanPolynomial):
            vars = other.vars().iterindex()
        elif PY_TYPE_CHECK(other, BooleanMonomial):
            vars = other.iterindex()
        else:
            t = other.variables()
            ovar_names = list(ovar_names)
            vars = [ovar_names.index(str(var)) for var in t]
    var_mapping = [None] * len(ovar_names)
    for i in vars:
        try:
            ind = my_names.index(ovar_names[i])
        except ValueError:
            # variable name not found in list of our variables
            # raise an exception and bail out
            raise NameError, "name %s not defined"%(ovar_names[i])
        var_mapping[i] = ring.gen(ind)
    return var_mapping

class BooleanMonomialMonoid(Monoid_class):
    def __init__(self, BooleanPolynomialRing polring):
        """
        Construct a boolean monomial monoid given a boolean polynomial
        ring.

        This object provides a parent for boolean monomials.

        INPUT:
            polring -- the polynomial ring our monomials lie in

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M
            MonomialMonoid of Boolean PolynomialRing in x, y

            sage: M.gens()
            (x, y)
            sage: type(M.gen(0))
            <type 'sage.rings.polynomial.pbori.BooleanMonomial'>
        """
        cdef BooleanMonomial m
        self._ring = polring
        ParentWithGens.__init__(self, GF(2), polring._names)

        m = new_BM(self, polring)
        PBMonom_construct(&m._pbmonom)
        self._one_element = m

    def _repr_(self):
        """
         sage: P.<x,y> = BooleanPolynomialRing(2)
         sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
         sage: M # indirect doctest
         MonomialMonoid of Boolean PolynomialRing in x, y
        """
        return "MonomialMonoid of %s" % (str(self._ring))

    def __hash__(self):
        """
        Return a hash for this monoid.

        EXAMPLE:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: {M:1} # indirect doctest
            {MonomialMonoid of Boolean PolynomialRing in x, y: 1}
        """
        return hash(str(self))

    def ngens(self):
        """
        Returns the number of variables in this monoid.

        EXAMPLES:
            sage: P = BooleanPolynomialRing(100, 'x')
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M.ngens()
            100
        """
        return self._ring.ngens()

    def gen(self, int i=0):
        """
        Return the i-th generator of self.

        INPUT:
            i -- an integer

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M.gen(0)
            x
            sage: M.gen(2)
            z

            sage: P = BooleanPolynomialRing(1000, 'x')
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M.gen(50)
            x50
        """
        if i < 0 or i >= self.ngens():
            raise ValueError, "Generator not defined."

        return new_BM_from_DD(self, (<BooleanPolynomialRing>self._ring),
                (<BooleanPolynomialRing>self._ring)._pbring.variable(i))

    def gens(self):
        """
        Return the tuple of generators of this monoid.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M.gens()
            (x, y, z)
        """
        return tuple([new_BM_from_DD(self, (<BooleanPolynomialRing>self._ring),
            (<BooleanPolynomialRing>self._ring)._pbring.variable(i)) \
                for i in xrange(self.ngens())])

    def _coerce_impl(self, other):
        """
        Canonical conversion of elements from other objects to this
        monoid.

        EXAMPLES:

        Coerce elements of self.

            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: x_monom = M(x); x_monom
            x
            sage: M._coerce_(x_monom) # indirect doctest
            x

        Coerce elements from BooleanMonomialMonoids where the generators of
        self include the generators of the other monoid.
            sage: R.<z,y> = BooleanPolynomialRing(2)
            sage: N = sage.rings.polynomial.pbori.BooleanMonomialMonoid(R)
            sage: m = M._coerce_(N(y*z)); m
            y*z
            sage: m.parent() is M
            True

        TESTS:
            sage: R.<t,y> = BooleanPolynomialRing(2)
            sage: N = sage.rings.polynomial.pbori.BooleanMonomialMonoid(R)
            sage: m = M._coerce_(N(y)); m
            Traceback (most recent call last):
            ...
            ValueError: cannot coerce monomial y to MonomialMonoid of Boolean PolynomialRing in x, y, z: name t not defined

            sage: R.<t,x,y,z> = BooleanPolynomialRing(4)
            sage: N = sage.rings.polynomial.pbori.BooleanMonomialMonoid(R)
            sage: m = M._coerce_(N(x*y*z)); m
            Traceback (most recent call last):
            ...
            TypeError: coercion from <type 'sage.rings.polynomial.pbori.BooleanMonomial'> to MonomialMonoid of Boolean PolynomialRing in x, y, z not implemented
        """
        if PY_TYPE_CHECK(other, BooleanMonomial) and \
            ((<BooleanMonomial>other)._parent.ngens() <= \
            (<BooleanPolynomialRing>self._ring)._pbring.nVariables()):
                try:
                    var_mapping = get_var_mapping(self, other.parent())
                except NameError, msg:
                    raise ValueError, "cannot coerce monomial %s to %s: %s"%(other,self,msg)
                m = self._one_element
                for i in other.iterindex():
                    m *= var_mapping[i]
                return m
        raise TypeError, "coercion from %s to %s not implemented" % \
            (type(other), str(self))

    def __call__(self, other = None):
        r"""
        Convert elements of other objects to elements of this monoid.

        INPUT:
            other -- element to convert,
                     - if \code{None} a BooleanMonomial representing 1
                       is returned
                     - only \code{BooleanPolynomials} with the same
                       parent ring as \code{self} which have a single
                       monomial is converted

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: x_monom = M(x); x_monom
            x

            sage: M(x*y)
            x*y

            sage: M(x+y)
            Traceback (most recent call last):
            ...
            TypeError: cannot convert to BooleanMonomialMonoid

        Convert elements of self.

            sage: M(x_monom)
            x

        Convert from other BooleanPolynomialRings.

            sage: R.<z,x> = BooleanPolynomialRing(2)
            sage: t = M(z); t
            z
            sage: t.parent() is M
            True

        Convert BooleanMonomials over other BooleanPolynomialRings.

            sage: N = sage.rings.polynomial.pbori.BooleanMonomialMonoid(R)
            sage: t = M(N(x*z)); t
            x*z
            sage: t.parent() is M
            True
        """
        cdef BooleanMonomial m
        cdef PBMonom t

        # this is needed for the PolyBoRi python code
        if other is None:
            return self._one_element

        try:
            return self._coerce_(other)
        except ValueError:
            pass
        except TypeError:
            pass

        if PY_TYPE_CHECK(other, BooleanPolynomial) and \
            (<BooleanPolynomial>other)._pbpoly.isSingleton():
                if (<BooleanPolynomial>other)._parent is self._ring:
                    return new_BM_from_PBMonom(self,
                            (<BooleanPolynomialRing>self._ring),
                            (<BooleanPolynomial>other)._pbpoly.lead())
                elif ((<BooleanPolynomial>other)._pbpoly.nUsedVariables() <= \
                    (<BooleanPolynomialRing>self._ring)._pbring.nVariables()):
                        try:
                            var_mapping = get_var_mapping(self, other)
                        except NameError, msg:
                            raise ValueError, "cannot convert polynomial %s to %s: %s"%(other,self,msg)
                        t = (<BooleanPolynomial>other)._pbpoly.lead()

                        m = self._one_element
                        for i in new_BMI_from_PBMonomIter(t, t.begin()):
                            m*= var_mapping[i]
                        return m
                else:
                    raise ValueError, "cannot convert polynomial %s to %s: %s"%(other,self,msg)

        elif PY_TYPE_CHECK(other, BooleanMonomial) and \
            ((<BooleanMonomial>other)._pbmonom.deg() <= \
            (<BooleanPolynomialRing>self._ring)._pbring.nVariables()):
                try:
                    var_mapping = get_var_mapping(self, other)
                except NameError, msg:
                    raise ValueError, "cannot convert monomial %s to %s: %s"%(other,self,msg)
                m = self._one_element
                for i in other:
                    m *= var_mapping[i]
                return m
        elif PY_TYPE_CHECK(other, Element) and \
                self.base_ring().has_coerce_map_from(other.parent()) and \
                        self.base_ring()(other).is_one():
                            return self._one_element
        elif isinstance(other, (int, long)) and other % 2:
            return self._one_element

        raise TypeError, "cannot convert to BooleanMonomialMonoid"

cdef class BooleanMonomial(MonoidElement):
    def __init__(self, parent):
        r"""
        Construct a boolean monomial object.

        INPUT:
            parent -- parent monoid this element lies in

        EXAMPLE:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: sage.rings.polynomial.pbori.BooleanMonomial(M)
            1

        NOTE: Use the \code{__call__} method of
        \code{BooleanMonomialMonoid} and not this constructor to
        construct these objects.
        """
        PBMonom_construct(&self._pbmonom)
        _parent = <ParentWithBase>parent
        self.ring = parent._ring

    def __dealloc__(self):
        PBMonom_destruct(&self._pbmonom)

    def __richcmp__(left, right, int op):
        """
        Compare BooleanMonomial objects.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M(x) < M(y)
            False

            sage: M(x) > M(y)
            True

            sage: M(x) == M(x)
            True

            sage: M(x) <= M(x)
            True

            sage: M(x) >= M(x)
            True
        """
        # boilerplate code from sage.structure.parent
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        cdef comparecodes res
        res = left._pbmonom.compare((<BooleanMonomial>right)._pbmonom)
        return res

    def _repr_(self):
        """
        Return a string representing self.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M(x*y) # indirect doctest
            x*y

            sage: R.<t,u> = BooleanPolynomialRing(2)
            sage: M(x*y)
            x*y
        """
        return PBMonom_to_str(&self._pbmonom)

    def _eval(self, d):
        """
        Evaluate this monomial.

        INPUT:
            d -- dictionary with integer indicies

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: m = P._monom_monoid(x*y)
            sage: m._eval({0:y,1:z})
            y*z
        """
        res = 1
        for i in self.iterindex():
            if d.has_key(i):
                res *= d[i]
            else:
                res *= self._parent.gen(i)
        return res


    def __call__(self, *args, **kwds):
        """
        Evaluate this monomial.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y
            sage: m = f.lm()
            sage: m(B(0),B(1))
            0
            sage: m(x=B(1))
            y
        """
        P = self.parent()
        if args and kwds:
            raise ValueError, "Using keywords and regular arguments not supported."
        if args:
            d = {}
            if len(args) > self._parent.ngens():
                raise ValueError, "Number of arguments is greater than the number of variables of parent ring."
            for i in range(len(args)):
                d[i] = args[i]
        elif kwds:
            d = list(self._parent.gens())
            gd = dict(zip(self._parent.variable_names(),range(len(d))))
            for var,val in kwds.iteritems():
                d[gd[var]] = val
        res = self._parent._one_element
        for var in self.iterindex():
            res *= d[var]
        return res

    def __hash__(self):
        """
        Return a hash of this monomial.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y
            sage: m = f.lm()
            sage: {m:1} #indirect doctest
            {x*y: 1}
        """
        return self._pbmonom.hash()

    def stableHash(self):
        return self._pbmonom.stableHash()

    def index(self):
        """
        Return the variable index of the first variable in this
        monomial.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y
            sage: m = f.lm()
            sage: m.index()
            0
        """
        return self._pbmonom.firstIndex()

    def deg(BooleanMonomial self):
        """
        Return total degree of this monomial.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M(x*y).deg()
            2

            sage: M(x*x*y*z).deg()
            3
        """
        return self._pbmonom.deg()

    def divisors(self):
        """
        Return a set of boolean monomials with all divisors of this
        monomial.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y
            sage: m = f.lm()
            sage: m.divisors()
            {{x,y}, {x}, {y}, {}}
        """
        return new_BS_from_PBSet(self._pbmonom.divisors(), self.ring)

    def multiples(self, BooleanMonomial rhs):
        r"""
        Return a set of boolean monomials with all multiples of this
        monomial up the the bound \var{rhs}.

        INPUT:
            rhs -- a boolean monomial

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x
            sage: m = f.lm()
            sage: g = x*y*z
            sage: n = g.lm()
            sage: m.multiples(n)
            {{x,y,z}, {x,y}, {x,z}, {x}}
            sage: n.multiples(m)
            {{x,y,z}}

        NOTE: The returned set always contains \code{self} even if the
        bound \var{rhs} is smaller than \code{self}.
        """
        return new_BS_from_PBSet(self._pbmonom.multiples(rhs._pbmonom),
                self.ring)

    def reducibleBy(self, BooleanMonomial rhs):
        r"""
        Return \code{True} if \code{self} is reducible by \var{rhs}.

        INPUT:
            rhs -- a boolean monomial

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y
            sage: m = f.lm()
            sage: m.reducibleBy((x*y).lm())
            True
            sage: m.reducibleBy((x*z).lm())
            False
        """
        return self._pbmonom.reducibleBy(rhs._pbmonom)

    def set(self):
        r"""
        Return a boolean set of variables in this monomials.

        EXAMPLE:
           sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y
            sage: m = f.lm()
            sage: m.set()
            {{x,y}}
        """
        return new_BS_from_PBSet(self._pbmonom.set(), self.ring)

    def __len__(BooleanMonomial self):
        """
        Return number of variables in this monomial. This is
        equivalent to the total degree for BooleanMonomials.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: len(M(x*y))
            2
        """
        return self._pbmonom.deg()

    def __iter__(self):
        """
        Return an iterator over the variables in this monomial.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: list(M(x*z)) # indirect doctest
            [x, z]
        """
        return new_BMVI_from_BooleanMonomial(self)

    def iterindex(self):
        """
        Return an iterator over the indicies of the variables in self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: list(M(x*z).iterindex())
            [0, 2]
        """
        return new_BMI_from_PBMonomIter(self._pbmonom, self._pbmonom.begin())

    cdef MonoidElement _mul_c_impl(left, MonoidElement right):
        """
        Multiply self with another boolean monomial.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: x = M(x); xy = M(x*y); z=M(z)
            sage: x*x
            x

            sage: xy*y
            x*y

            sage: xy*z
            x*y*z

            sage: x*1   # todo: not implemented
            x
            sage: 1*x   # todo: not implemented

            sage: x*int(1)
            x
            sage: int(1)*x
            x
        """
        cdef BooleanMonomial m = new_BM_from_PBMonom(\
                (<BooleanMonomial>left)._parent,
                (<BooleanMonomial>left).ring,
                (<BooleanMonomial>left)._pbmonom)
        m._pbmonom.imul( (<BooleanMonomial>right)._pbmonom )
        return m

    def __add__(left, right):
        """
        Addition operator. Returns a boolean polynomial.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: x = M(x); xy = M(x*y)
            sage: x + xy
            x*y + x

            sage: x+0
            x
            sage: 0+x   # todo: not implemented
            x

            sage: x+1
            x + 1
            sage: 1 + x     # todo: not implemented
            x + 1
        """
        # Using canonical coercion is not possible for this case.
        # The coercion model cannot find the common parent
        # BooleanPolynomialRing for argument types BooleanMonomial and Integer.
        # This is a common case we should handle as in the examples above.
        # Since we know the result will be a BooleanPolynomial, we let
        # BooleanPolynomial handle the coercion.
        cdef BooleanPolynomial res
        cdef BooleanMonomial monom
        if PY_TYPE_CHECK(left, BooleanMonomial):
            monom = left
            other = right
        elif PY_TYPE_CHECK(right, BooleanMonomial):
            monom = right
            other = left
        else:
            raise TypeError, "BooleanMonomial.__add__ called with not supported types %s and %s" %(type(right),type(left))

        res = new_BP_from_PBMonom(monom.ring, monom._pbmonom)
        return res.__iadd__(monom.ring._coerce_c(other))

###
#
# Various internal constructors for boolean polynomials from various
# other formats.
#
###

cdef inline BooleanMonomial new_BM(parent, BooleanPolynomialRing ring):
    cdef BooleanMonomial m
    m = <BooleanMonomial>PY_NEW(BooleanMonomial)
    m._parent = parent
    m.ring = ring
    return m

cdef inline BooleanMonomial new_BM_from_PBMonom(parent,
        BooleanPolynomialRing ring, PBMonom juice):
    cdef BooleanMonomial m = new_BM(parent, ring)
    PBMonom_construct_pbmonom(&m._pbmonom,juice)
    return m

cdef inline BooleanMonomial new_BM_from_PBVar(parent,
        BooleanPolynomialRing ring, PBVar juice):
    cdef BooleanMonomial m = new_BM(parent, ring)
    PBMonom_construct_pbvar(&m._pbmonom,juice)
    return m

cdef inline BooleanMonomial new_BM_from_DD(parent,
        BooleanPolynomialRing ring, PBDD juice):
    cdef BooleanMonomial m = new_BM(parent, ring)
    PBMonom_construct_dd(&m._pbmonom,juice)
    return m

cdef class BooleanMonomialVariableIterator:
    def __iter__(self):
        """
        Return an iterator over the variables of a boolean monomial.

        EXAMPLE:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y + z + 1
            sage: for m in f: list(m)# indirect doctest
            [x, y]
            [z]
            []
        """
        return self

    def __dealloc__(self):
        PBMonomVarIter_destruct(&self.iter)

    def __next__(self):
        """
        EXAMPLE:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y + z + 1
            sage: m = f.lm()
            sage: iter(m).next()
            x
        """
        cdef PBVar val
        if self.iter.equal(self.obj.variableEnd()):
            raise StopIteration
        val = self.iter.value()
        self.iter.next()
        return new_BM_from_PBVar(self.parent, self.ring, val)

cdef inline BooleanMonomialVariableIterator new_BMVI_from_BooleanMonomial(\
                            BooleanMonomial monom):
    """
    Construct a new BooleanMonomialIterator
    """
    cdef BooleanMonomialVariableIterator m
    m = <BooleanMonomialVariableIterator>PY_NEW(BooleanMonomialVariableIterator)
    m.parent = monom._parent
    m.ring = monom.ring
    m.obj = monom._pbmonom
    m.iter = m.obj.variableBegin()
    return m

cdef class BooleanMonomialIterator:
    """
    An iterator over the variable indices of a monomial.
    """
    def __iter__(self):
        """
        EXAMPLE:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y + z + 1
            sage: for m in f: list(m.iterindex())# indirect doctest
            [0, 1]
            [2]
            []
        """
        return self

    def __dealloc__(self):
        PBMonomIter_destruct(&self._iter)

    def __next__(self):
        """
        EXAMPLE:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y + z + 1
            sage: m = f.lm()
            sage: m.iterindex().next()
            0
        """
        cdef int val
        if self._iter.hash() == self._obj.end().hash():
            raise StopIteration
        val = self._iter.value()
        self._iter.next()
        return val

cdef inline BooleanMonomialIterator new_BMI_from_PBMonomIter(\
                            PBMonom parent, PBMonomIter juice):
    """
    Construct a new BooleanMonomialIterator
    """
    cdef BooleanMonomialIterator m
    m = <BooleanMonomialIterator>PY_NEW(BooleanMonomialIterator)
    m._iter = juice
    m._obj = parent
    return m

cdef class BooleanPolynomial(MPolynomial):
    def __init__(self, parent):
        r"""
        Construct a boolean polynomial object in the given boolean
        polynomial ring.

        INPUT:
            parent -- a boolean polynomial ring

        TEST:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: sage.rings.polynomial.pbori.BooleanPolynomial(B)
           0

        NOTE: Do not use this method to construct boolean polynomials,
        but use the approriate \code{__call__} method in the parent.
        """
        PBPoly_construct(&self._pbpoly)
        self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBPoly_destruct(&self._pbpoly)

    def _repr_(self):
        """
        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: repr(a+b+z^2+1) # indirect doctest
           'a + b + z + 1'
        """
        return PBPoly_to_str(&self._pbpoly)

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        """
        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: f = a*z + b + 1
           sage: g = b + z
           sage: f + g
           a*z + z + 1
        """
        cdef BooleanPolynomial p = new_BP_from_PBPoly(\
                (<BooleanPolynomial>left)._parent, (<BooleanPolynomial>left)._pbpoly)
        p._pbpoly.iadd( (<BooleanPolynomial>right)._pbpoly )
        return p

    cdef ModuleElement _sub_c_impl(left, ModuleElement right):
        """
        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: f = a*z + b + 1
           sage: g = b + z
           sage: f - g
           a*z + z + 1
        """
        return left._add_c_impl(right)

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        """
        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: k = B.base_ring()
           sage: f = a*z + b + 1
           sage: f*k(1)
           a*z + b + 1
        """
        if left:
            return new_BP_from_PBPoly(left._parent, self._pbpoly)
        else:
            return 0

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        """
        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: k = B.base_ring()
           sage: f = a*z + b + 1
           sage: k(0)*f
           0
        """
        return self._rmul_c_impl(right)

    cdef RingElement _mul_c_impl(left, RingElement right):
        """
        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: f = a*z + b + 1
           sage: g = b + z
           sage: f * g
           a*b*z + a*z + b*z + z
        """
        cdef BooleanPolynomial p = new_BP_from_PBPoly(\
                (<BooleanPolynomial>left)._parent, (<BooleanPolynomial>left)._pbpoly)
        p._pbpoly.imul( (<BooleanPolynomial>right)._pbpoly )
        return p

    def is_equal(self, BooleanPolynomial right):
        """
        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: f = a*z + b + 1
           sage: g = b + z
           sage: f.is_equal(g)
           False

           sage: f.is_equal( (f + 1) - 1 )
           True
        """
        return self._pbpoly.is_equal(right._pbpoly)

    def __richcmp__(left, right, int op):
        """
        Compare left and right and return -1, 0, 1 for <, ==, and > respectively.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: x < x+y
            True

            sage: y*z < x
            True

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: y*z < x
            False

            sage: P(0) == 0
            True
        """
        #boilerplate from sage.structure.element
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        # see __richcmp__
        cdef int res
        from itertools import izip
        for lm, rm in izip(left, right):
            res = cmp(lm, rm)
            if res != 0:
                return res
        return cmp(len(left),len(right))

    def __iter__(self):
        r"""
        Return an iterator over the monomials of \code{self}, in the order of
        the parent ring.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: p = x + z + x*y + y*z + x*y*z
            sage: list(iter(p))
            [x*y*z, x*y, x, y*z, z]

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: p = x + z + x*y + y*z + x*y*z
            sage: list(iter(p))
            [x*y*z, x*y, y*z, x, z]

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='degrevlex')
            sage: p = x + z + x*y + y*z + x*y*z
            sage: list(iter(p))
            [z*y*x, y*x, z*y, x, z]

        TESTS:
            sage: R = BooleanPolynomialRing(1,'y')
            sage: list(iter(y))
            [y]
            sage: R
            Boolean PolynomialRing in y
        """
        return new_BPI_from_PBPolyIter(self, self._pbpoly.orderedBegin())

    def __pow__(BooleanPolynomial self, int exp, ignored):
        r"""
        Return \code{self^(exp)}.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: p = x + y
            sage: p^0
            1

            sage: p^1
            x + y

            sage: p^5
            x + y

            sage: p^-1
            Traceback (most recent call last):
            ...
            NotImplementedError: Negative exponents for non constant boolean polynomials not implemented.

            sage: z = P(0)
            sage: z^0
            1

            sage: z^1
            0
        """
        if exp > 0:
            return self
        elif exp == 0:
            return self._parent._one_element
        elif self._pbpoly.isOne():
            return self
        elif self._pbpoly.isZero():
            raise ZeroDivisionError
        else:
            raise NotImplementedError, "Negative exponents for non constant boolean polynomials not implemented."

    def __neg__(BooleanPolynomial self):
        """
        Return -self.

        EXAMPLE:
           sage: B.<a,b,z> = BooleanPolynomialRing(3)
           sage: f = a*z + b + 1
           sage: -f
           a*z + b + 1
        """
        return self

    def total_degree(BooleanPolynomial self):
        """
        Return the total degree of self.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: (x+y).total_degree()
            1

            sage: P(1).total_degree()
            0

            sage: (x*y + x + y + 1).total_degree()
            2
        """
        return self._pbpoly.deg()

    def degree(self):
        """
        Return the total degree of self.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: (x+y).degree()
            1

            sage: P(1).degree()
            0

            sage: (x*y + x + y + 1).degree()
            2
        """
        return self._pbpoly.deg()

    def lm(BooleanPolynomial self):
        r"""
        Return the leading monomial of boolean polynomial, with
        respect to the order of parent ring.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x+y+y*z).lm()
            x

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: (x+y+y*z).lm()
            y*z
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._parent,
                self._pbpoly.lead())

    def lt(BooleanPolynomial self):
        """
        Return the leading term of this boolean polynomial, with
        respect to the order of the parent ring.

        Note that for boolean polynomials this is equivalent to
        returning leading monomials.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x+y+y*z).lt()
            x

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: (x+y+y*z).lt()
            y*z
        """
        return self.lm()

    def is_zero(BooleanPolynomial self):
        r"""
        Check if \code{self} is zero.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: P(0).is_zero()
            True

            sage: x.is_zero()
            False

            sage: P(1).is_zero()
            False
        """
        return self._pbpoly.isZero()

    def __nonzero__(self):
        r"""
        Check if \code{self} is zero.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: bool(P(0))
            False

            sage: bool(x)
            True

            sage: bool(P(1))
            True
        """
        return not self._pbpoly.isZero()

    def is_one(BooleanPolynomial self):
        """
        Check if self is 1.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: P(1).is_one()
            True

            sage: P.one_element().is_one()
            True

            sage: x.is_one()
            False

            sage: P(0).is_one()
            False
        """
        return self._pbpoly.isOne()

    def is_unit(BooleanPolynomial self):
        """
        Check if self is invertible in the parent ring.

        Note that this condition is equivalent to being 1 for boolean polynomials.

        EXAMPLE:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: P.one_element().is_unit()
            True

            sage: x.is_unit()
            False
        """
        return self._pbpoly.isOne()

    def is_constant(BooleanPolynomial self):
        """
        Check if self is constant.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: P(1).is_constant()
            True

            sage: P(0).is_constant()
            True

            sage: x.is_constant()
            False

            sage: (x*y).is_constant()
            False
        """
        return self._pbpoly.isConstant()

    def lm_degree(BooleanPolynomial self):
        """
        Returns the total degree of the leading monomial of self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: p = x + y*z
            sage: p.lm_degree()
            1

            sage: P.<x,y,z> = BooleanPolynomialRing(3,order='deglex')
            sage: p = x + y*z
            sage: p.lm_degree()
            2

            sage: P(0).lm_degree()
            0
        """
        return self._pbpoly.lmDeg()

    def vars(self):
        r"""
        Return a boolean monomial with all the variables appearing in
        \code{self}.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x + y).vars()
            x*y

            sage: (x*y + z).vars()
            x*y*z

            sage: P.zero_element().vars()
            1

            sage: P.one_element().vars()
            1

        TESTS:
            sage: R = BooleanPolynomialRing(1, 'y')
            sage: y.vars()
            y
            sage: R
            Boolean PolynomialRing in y
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid,
                self._parent, self._pbpoly.usedVariables())

    def variables(self):
        r"""
        Return a list of all variables appearing in \code{self}.

        EXAMPLE:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x + y).variables()
            [x, y]

            sage: (x*y + z).variables()
            [x, y, z]

            sage: P.zero_element().variables()
            []

            sage: P.one_element().variables()
            [1]
        """
        P = self.parent()
        o = P.one_element()
        if self is o or self == o:
            return [o]
        return list(self.vars())

    def nvariables(self):
        """
        Return the number of variables used to form this boolean
        polynomial.

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4)
            sage: f = a*b*c + 1
            sage: f.nvariables()
            3
        """
        return self._pbpoly.nUsedVariables()

    def monomials(self):
        r"""
        Return a list of monomials appearing in \code{self} ordered
        largest to smallest.

        EXAMPLE:
            sage: P.<a,b,c> = BooleanPolynomialRing(3,order='lex')
            sage: f = a + c*b
            sage: f.monomials()
            [a, b*c]

            sage: P.<a,b,c> = BooleanPolynomialRing(3,order='degrevlex')
            sage: f = a + c*b
            sage: f.monomials()
            [c*b, a]
        """
        return list(self)

    def monomial_coefficient(self, mon):
        r"""
        Return the coefficient of the monomial \var{mon} in
        \code{self}, where \var{mon} must have the same parent as
        \code{self}.

        INPUT:
            mon -- a monomial

        EXAMPLE:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: x.monomial_coefficient(x)
            1
            sage: x.monomial_coefficient(y)
            0
            sage: R.<x,y,z,a,b,c>=BooleanPolynomialRing(6)
            sage: f=(1-x)*(1+y); f
            x*y + x + y + 1

            sage: f.monomial_coefficient(1)
            1

            sage: f.monomial_coefficient(0)
            0
        """
        cdef BooleanPolynomialRing B = <BooleanPolynomialRing>self._parent
        k = B._base
        mon = B._coerce_c(mon)
        if mon in set(self.set()):
            return k._one_element
        else:
            return k._zero_element

    def __hash__(self):
        r"""
        Return hash for \code{self}.

        EXAMPLE:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: {x:1} # indirect doctest
            {x: 1}
        """
        return self._pbpoly.hash()

    def __len__(self):
        """
        Return number of monomials in self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: len(x + y)
            2

            sage: len(P.one_element())
            1

            sage: len(x*y + y + z + x*z)
            4

            sage: len(P.zero_element())
            0
        """
        return self._pbpoly.length()

    def __call__(self, *args, **kwds):
        """
        Evaluate this boolean polynomials.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y + z + 1
            sage: f(0,1,1)
            0
            sage: f(z,y,x)
            x + y*z + 1
            sage: f(x=z)
            y*z + z + 1

            sage: P.<a,b,c> = PolynomialRing(QQ)
            sage: f(a,b,c)
            a*b + c + 1
            sage: f(x=a,y=b,z=1)
            a*b + 2

        Evaluation of polynomials can be used fully symbolic:

            sage: f(x=var('a'),y=var('b'),z=var('c'))
            c + a*b + 1
            sage: f(var('a'),var('b'),1)
            a*b
        """
        P = self._parent
        cdef int N = P.ngens()
        if args and kwds:
            raise ValueError, "Using keywords and regular arguments not supported."
        if args:
            d = {}
            if len(args) != N:
                raise ValueError, "Number of arguments is different from the number of variables of parent ring."
            for i in range(N):
                arg = args[i]
                try:
                    arg = P(arg)
                    if arg.constant():
                        # TODO: We should collect those and reduce once only
                        self = ll_red_nf(self, (P.gen(i) + arg).set())
                    else:
                        d[i] = arg
                except TypeError:
                    d[i] = arg
            if not len(d):
                return self
        elif kwds:
            d = dict(zip(range(P.ngens()), P.gens()))
            gd = dict(zip(P.variable_names(),range(P.ngens())))
            for var,val in kwds.iteritems():
                d[gd[var]] = val

        res = 0
        for m in self:
            res += m._eval(d)
        return res

    def subs(self, in_dict=None, **kwds):
        r"""
        Fixes some given variables in a given boolean polynomial and
        returns the changed boolean polynomials. The polynomial itself
        is not affected. The variable,value pairs for fixing are to be
        provided as dictionary of the form \code{\{variable:value\}}
        or named parameters (see examples below).

        INPUT:
            in_dict -- (optional) dict with variable:value pairs
            **kwds -- names parameters

        EXAMPLE:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: f = x*y + z + y*z + 1
            sage: f.subs(x=1)
            y*z + y + z + 1
            sage: f.subs(x=0)
            y*z + z + 1

            sage: f.subs(x=y)
            y*z + y + z + 1

            sage: f.subs({x:1},y=1)
            0
            sage: f.subs(y=1)
            x + 1
            sage: f.subs(y=1,z=1)
            x + 1
            sage: f.subs(z=1)
            x*y + y
            sage: f.subs({'x':1},y=1)
            0

            This method can work fully symbolic:

            sage: f.subs(x=var('a'),y=var('b'),z=var('c'))
            b*c + c + a*b + 1

            sage: f.subs({'x':var('a'),'y':var('b'),'z':var('c')})
            b*c + c + a*b + 1
        """
        P = self._parent

        fixed = {}
        if in_dict is not None:
            for var,val in in_dict.iteritems():
                if PY_TYPE_CHECK(var, basestring):
                    var = P(var)
                elif var.parent() is not P:
                    var = P(var)
                try:
                    v = P(val)
                    if v.constant():
                        self = ll_red_nf(self, (var + v).set())
                    else:
                        fixed[var.lm().index()] = val
                except TypeError:
                    fixed[var.lm().index()] = val
        if kwds:
            gdict = P._monom_monoid.gens_dict()

        for var,val in kwds.iteritems():
            var = gdict[var]
            try:
                v =  P(val)
                if v.constant():
                    self = ll_red_nf(self, (var + v).set())
                else:
                    fixed[var.index()] = val
            except TypeError:
                fixed[var.index()] = val

        if not len(fixed):
            return self
        res = 0
        for m in self:
            res += m._eval(fixed)
        return res

    def __reduce__(self):
        """
        EXAMPLE:
            sage: P.<a,b> = BooleanPolynomialRing(2)
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_BooleanPolynomial, (self._parent, PBPoly_to_str(&self._pbpoly))

    def set(self):
        r"""
        Return a \code{BooleSet} with all monomials apprearing in this
        polynomial.

        EXAMPLE:
            sage: B.<a,b,z> = BooleanPolynomialRing(3)
            sage: (a*b+z+1).set()
            {{a,b}, {z}, {}}
        """
        return new_BS_from_PBSet(self._pbpoly.set(), self._parent)

    def deg(self):
        """
        Return the total degree of self.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: (x+y).degree()
            1

            sage: P(1).degree()
            0

            sage: (x*y + x + y + 1).degree()
            2
        """
        return self._pbpoly.deg()

    def elength(self):
        return self._pbpoly.eliminationLength()

    def lead(self):
        r"""
        Return the leading monomial of boolean polynomial, with
        respect to to the order of parent ring.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x+y+y*z).lead()
            x

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: (x+y+y*z).lead()
            y*z
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._parent,
                                   self._pbpoly.lead())

    def lexLead(self):
        r"""
        Return the leading monomial of boolean polynomial, with
        respect to the lexicographical term ordering..

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x+y+y*z).lexLead()
            x

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: (x+y+y*z).lexLead()
            x
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._parent,
                                                self._pbpoly.lexLead())
    def lexLmDeg(self):
        """
        Return degree of leading monomial with respect to the
        lexicographical ordering.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3,order='lex')
            sage: f = x + y*z
            sage: f
            x + y*z
            sage: f.lexLmDeg()
            1

            sage: B.<x,y,z> = BooleanPolynomialRing(3,order='deglex')
            sage: f = x + y*z
            sage: f
            y*z + x
            sage: f.lexLmDeg()
            1
        """
        return self._pbpoly.lexLmDeg()

    def lmDeg(self):
        """
        Return the degree of the leading monomial with respect to the
        lexicographical orderings.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3,order='lex')
            sage: f = x + y*z
            sage: f
            x + y*z
            sage: f.lmDeg()
            1

            sage: B.<x,y,z> = BooleanPolynomialRing(3,order='deglex')
            sage: f = x + y*z
            sage: f
            y*z + x
            sage: f.lmDeg()
            2
        """
        return self._pbpoly.lmDeg()

    def constant(self):
        r"""
        Return \code{True} if this element is constant.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: x.constant()
            False

            sage: B(1).constant()
            True
        """
        return self._pbpoly.isConstant()

    def isZero(self):
        r"""
        Return \code{True} if this element is zero.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: x.isZero()
            False
            sage: B(0).isZero()
            True
            sage: B(1).isZero()
            False
        """
        return self._pbpoly.isZero()

    def isOne(self):
        r"""
        Return \code{True} if this element is one.

        EXAMPLE:
            sage: B.<x,y,z> = BooleanPolynomialRing(3)
            sage: x.isOne()
            False
            sage: B(0).isOne()
            False
            sage: B(1).isOne()
            True
        """
        return self._pbpoly.isOne()

    def navigation(self):
        return new_CN_from_PBNavigator(self._pbpoly.navigation())

    def mapEveryXToXPlusOne(self):
        r"""
        Map every variable \var{x_i} in this polynomial to $x_i + 1$.

        EXAMPLE:
            sage: B.<a,b,z> = BooleanPolynomialRing(3)
            sage: f = a*b + z + 1; f
            a*b + z + 1
            sage: f.mapEveryXToXPlusOne()
            a*b + a + b + z + 1
            sage: f(a+1,b+1,z+1)
            a*b + a + b + z + 1
        """
        return new_BP_from_PBPoly(self._parent,
                pb_map_every_x_to_x_plus_one(self._pbpoly))

    def lmDivisors(self):
        r"""
        Return a \code{BooleSet} of all divisors of the leading
        monomial.

        EXAMPLE:
            sage: B.<a,b,z> = BooleanPolynomialRing(3)
            sage: f = a*b + z + 1
            sage: f.lmDivisors()
            {{a,b}, {a}, {b}, {}}
        """
        return new_BS_from_PBSet(self._pbpoly.lmDivisors(), self._parent)

    def firstTerm(self):
        r"""
        Return the first term with respect to the lexicographical term
        ordering.

        EXAMPLE:
            sage: B.<a,b,z> = BooleanPolynomialRing(3,order='lex')
            sage: f = b*z + a + 1
            sage: f.firstTerm()
            a
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._parent,
                self._pbpoly.firstTerm())

    def reducibleBy(self, BooleanPolynomial rhs):
        r"""
        Return \code{True} if this boolean polynomial is reducible by
        the polynomial \var{rhs}.

        INPUT:
            rhs -- a boolean polynomial

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4,order='degrevlex')
            sage: f = (a*b + 1)*(c + 1)
            sage: f.reducibleBy(d)
            False
            sage: f.reducibleBy(c)
            True
            sage: f.reducibleBy(c + 1)
            True

        """
        return self._pbpoly.reducibleBy(rhs._pbpoly)

    def nNodes(self):
        return self._pbpoly.nNodes()

    def nVars(self):
        """
        Return the number of variables used to form this boolean
        polynomial.

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4)
            sage: f = a*b*c + 1
            sage: f.nVars()
            3
        """
        return self._pbpoly.nUsedVariables()

    def totalDegree(self):
        """
        Return total degree of this boolean polynomial.

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4)
            sage: f = a*b*c + 1
            sage: f.totalDegree()
            3
        """
        return self._pbpoly.totalDeg()

    def gradedPart(self, int deg):
        r"""
        Return graded part of this boolean polynomial of degree
        \var{deg}.

        INPUT:
            deg -- a degree

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4)
            sage: f = a*b*c + c*d + a*b + 1
            sage: f.gradedPart(2)
            a*b + c*d

            sage: f.gradedPart(0)
            1

        TESTS:
            sage: f.gradedPart(-1)
            0
        """
        return new_BP_from_PBPoly(self._parent,
                self._pbpoly.gradedPart(deg))

    def hasConstantPart(self):
        r"""
        Return \code{True} if this boolean polynomial has a constant
        part, i.e. if $1$ is a term.

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4)
            sage: f = a*b*c + c*d + a*b + 1
            sage: f.hasConstantPart()
            True

            sage: f = a*b*c + c*d + a*b
            sage: f.hasConstantPart()
            False
        """
        return self._pbpoly.hasConstantPart()

    def zeroesIn(self, BooleSet s):
        return new_BS_from_PBSet(pb_zeroes(self._pbpoly, s._pbset), self._parent)

    def spoly(self, BooleanPolynomial rhs):
        r"""
        Return the S-Polynomial of this boolean polynomial and the
        other boolean polynomial \var{rhs}.

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4)
            sage: f = a*b*c + c*d + a*b + 1
            sage: g = c*d + b
            sage: f.spoly(g)
            a*b + a*c*d + c*d + 1
        """
        return new_BP_from_PBPoly(self._parent,
                pb_spoly(self._pbpoly, rhs._pbpoly))

    def stableHash(self):
        return self._pbpoly.stableHash()

    def ring(self):
        """
        Return the parent of this boolean polynomial.

        EXAMPLE:
            sage: B.<a,b,c,d> = BooleanPolynomialRing(4)
            sage: a.ring() is B
            True
        """
        return self._parent

cdef class BooleanPolynomialIterator:
    def __iter__(self):
        return self

    def __dealloc__(self):
        PBPolyIter_destruct(&self._iter)

    def __next__(self):
        cdef PBMonom val
        if self._iter.equal(self._obj._pbpoly.orderedEnd()):
            raise StopIteration
        val = self._iter.value()
        self._iter.next()
        return new_BM_from_PBMonom(self._obj._parent._monom_monoid,
                self._obj._parent, val)

cdef inline BooleanPolynomialIterator new_BPI_from_PBPolyIter(\
                            BooleanPolynomial parent, PBPolyIter juice):
    """
    Construct a new BooleanMonomialIterator
    """
    cdef BooleanPolynomialIterator m
    m = <BooleanPolynomialIterator>PY_NEW(BooleanPolynomialIterator)
    m._iter = juice
    m._obj = parent
    return m

class BooleanPolynomialIdeal(MPolynomialIdeal):
    def __init__(self, ring, gens=[], coerce=True):
        r"""
        Construct an ideal in the boolean polynomial ring.

        INPUT:
            ring -- the ring this ideal is defined in
            gens -- a list of generators
            coerce -- coerce all elements to the ring \var{ring} (default: \code{True})

        EXAMPLES:
            sage: P.<x0, x1, x2, x3> = BooleanPolynomialRing(4)
            sage: I = P.ideal(x0*x1*x2*x3 + x0*x1*x3 + x0*x1 + x0*x2 + x0)
            sage: I
            Ideal (x0*x1*x2*x3 + x0*x1*x3 + x0*x1 + x0*x2 + x0) of Boolean PolynomialRing in x0, x1, x2, x3
        """
        MPolynomialIdeal.__init__(self, ring, gens, coerce)

    def groebner_basis(self, **kwds):
        r"""
        Return a Groebner basis of this ideal.

        INPUT:
            heuristic -- Turn off heuristic by setting \code{heuristic=False}
                         (default: True)
            lazy  --  (default: True)
            red_tail  --  use tail reduction (default: True)
            redsb  --  return reduced Groebner basis (default: True)
            minsb  --  (default: True)
            invert -- setting \code{invert=True} input and output get
                      a transformation $x+1$ for each variable $x$,
                      which shouldn't effect the calculated GB, but
                      the algorithm.
            prot  --  show protocol (default: False)
            full_prot  --  show full protocol (default: False)
            faugere -- use a variant of Faugere's F4 (default: False)
            aes  --  input is AES system (default: False)
            coding  --  input is coding theory system (default: False)
            ll  --  (default: False)
            llfirst  --  (default: False)
            llfirstonthefly  --  (default: False)
            gauss_on_linear_first  --  (default: True)
            linearAlgebraInLastBlock  --  (default: True)
            max_growth  --  (default: 2.0)
            exchange  --  (default: True)
            selection_size  --  (default: 1000)
            implementation  -- either 'Python' or anything else (default: 'Python')
            deg_bound  --  (default: 1000000000000)
            recursion  --  (default: False)
            implications  --  (default: False)
            step_factor  --  (default: 1)

        EXAMPLES:
            sage: P.<x0, x1, x2, x3> = BooleanPolynomialRing(4)
            sage: I = P.ideal(x0*x1*x2*x3 + x0*x1*x3 + x0*x1 + x0*x2 + x0)
            sage: I.groebner_basis()
            [x0*x1 + x0*x2 + x0, x0*x2*x3 + x0*x3]
        """
        from polybori.gbcore import groebner_basis
        _sig_on
        gb = groebner_basis(self.gens(), **kwds)
        _sig_off
        return Sequence(gb, self.ring(), check=False, immutable=True)

##
#
# Various internal constructors for boolean polynomials from various
# other formats.
#
##

cdef inline BooleanPolynomial new_BP(BooleanPolynomialRing parent):
    cdef BooleanPolynomial p
    p = <BooleanPolynomial>PY_NEW(BooleanPolynomial)
    p._parent = parent
    return p

cdef inline BooleanPolynomial new_BP_from_DD(BooleanPolynomialRing parent, PBDD juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_dd(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBPoly(BooleanPolynomialRing parent, PBPoly juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbpoly(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBMonom(BooleanPolynomialRing parent, PBMonom juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbmonom(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBSet(BooleanPolynomialRing parent, PBSet juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbset(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_int(BooleanPolynomialRing parent, int juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_int(&p._pbpoly,juice)
    return p

cdef class DD:
    def __call__(self):
        return self

    def __dealloc__(self):
        PBDD_destruct(&self._pbdd)

    def empty(self):
        return self._pbdd.emptiness()

    def navigation(self):
        return new_CN_from_PBNavigator(self._pbdd.navigation())

    def subset0(self, idx):
        return new_DD_from_PBDD(self._pbdd.subset0(idx))

    def subset1(self, idx):
        return new_DD_from_PBDD(self._pbdd.subset1(idx))

    def union(self, rhs):
        return new_DD_from_PBDD(self._pbdd.unite((<DD>rhs)._pbdd))

cdef inline DD new_DD_from_PBDD(PBDD juice):
    """
    Construct a new DD
    """
    cdef DD d
    d = <DD>PY_NEW(DD)
    d._pbdd = juice
    return d

cdef class BooleSet:
    def __init__(self, param=None, ring=None):
        if PY_TYPE_CHECK(param, CCuddNavigator):
            if ring is None:
                raise TypeError, "BooleSet constructor requires parent ring argument"
            self.ring = ring
            PBSet_construct_pbnav(&self._pbset, (<CCuddNavigator>param)._pbnav, (<BooleanPolynomialRing>ring)._pbring)
        elif PY_TYPE_CHECK(param, BooleSet):
            PBSet_construct_pbset(&self._pbset, (<BooleSet>param)._pbset)
            self.ring = (<BooleSet>param).ring
        else:
            PBSet_construct(&self._pbset)
            global cur_ring
            self.ring = cur_ring

    def __dealloc__(self):
        PBSet_destruct(&self._pbset)

    def __repr__(self):
        return PBSet_to_str(&self._pbset)

    def set(self):
        return self

    def empty(self):
        return self._pbset.emptiness()

    def navigation(self):
        return new_CN_from_PBNavigator(self._pbset.navigation())

    def cartesianProduct(self, BooleSet rhs):
        return new_BS_from_PBSet(
                self._pbset.cartesianProduct((<BooleSet>rhs)._pbset), self.ring)

    def diff(self, rhs):
        cdef PBSet s
        if PY_TYPE_CHECK(rhs, BooleSet):
            s = (<BooleSet>rhs)._pbset
        elif PY_TYPE_CHECK(rhs, BooleanPolynomial):
            s = (<BooleanPolynomial>rhs)._pbpoly.set()
        else:
            raise TypeError, "Argument 'rhs' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(rhs))
        return new_BS_from_PBSet(self._pbset.diff(s), self.ring)


    def change(self, ind):
        return new_BS_from_PBSet(self._pbset.change(ind), self.ring)

    def vars(self):
        return new_BM_from_PBMonom(self.ring._monom_monoid, self.ring,
                                            self._pbset.usedVariables())

    def nNodes(self):
        return self._pbset.nNodes()

    def nSupport(self):
        return self._pbset.nSupport()

    def union(self, rhs):
        cdef PBSet s
        if PY_TYPE_CHECK(rhs, BooleSet):
            s = (<BooleSet>rhs)._pbset
        elif PY_TYPE_CHECK(rhs, BooleanPolynomial):
            s = (<BooleanPolynomial>rhs)._pbpoly.set()
        else:
            raise TypeError, "Argument 'rhs' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(rhs))
        return new_BS_from_PBSet(self._pbset.unite(s), self.ring)

    def __iter__(self):
        """
        Create an iterator over elements of self.

        EXAMPLES:
            sage: P.<x, y> = BooleanPolynomialRing(2)
            sage: f = x*y+x+y+1; s = f.set()
            sage: list(s)
            [x*y, x, y, 1]
        """
        return new_BSI_from_PBSetIter(self)

    def __len__(self):
        return self._pbset.length()

    def __hash__(self):
        return self._pbset.hash()

    def __mod__(self, BooleSet vs):
        return mod_mon_set(self, vs)

    def __contains__(self, BooleanMonomial m):
        return self._pbset.owns(m._pbmonom)

    def stableHash(self):
        return self._pbset.stableHash()

    def divide(self, BooleanMonomial rhs):
        return new_BS_from_PBSet(self._pbset.divide(rhs._pbmonom), self.ring)

    def subset0(self, int i):
        return new_BS_from_PBSet(self._pbset.subset0(i), self.ring)

    def subset1(self, int i):
        return new_BS_from_PBSet(self._pbset.subset1(i), self.ring)

    def includeDivisors(self):
        return new_BS_from_PBSet(pb_include_divisors(self._pbset), self.ring)

    def minimalElements(self):
        return new_BS_from_PBSet(pb_minimal_elements(self._pbset), self.ring)


cdef inline BooleSet new_BS_from_PBSet(PBSet juice, BooleanPolynomialRing ring):
    """
    Construct a new BooleSet
    """
    cdef BooleSet s
    s = <BooleSet>PY_NEW(BooleSet)
    s._pbset = juice
    s.ring = ring
    return s

cdef class BooleSetIterator:
    def __iter__(self):
        return self

    def __dealloc__(self):
        PBSetIter_destruct(&self._iter)

    def __next__(self):
        cdef PBMonom val
        if self._iter.equal(self._obj.end()):
            raise StopIteration
        val = self._iter.value()
        self._iter.next()
        return new_BM_from_PBMonom(self._parent, self.ring, val)

cdef inline BooleSetIterator new_BSI_from_PBSetIter(BooleSet s):
    """
    Construct a new BooleSetIterator
    """
    cdef BooleSetIterator m
    m = <BooleSetIterator>PY_NEW(BooleSetIterator)
    m.ring = s.ring
    m._parent = m.ring._monom_monoid
    m._obj = s._pbset
    PBSetIter_construct(&m._iter)
    m._iter = s._pbset.begin()
    return m

cdef class CCuddNavigator:
    def __call__(self):
        return self

    def __dealloc__(self):
        PBNavigator_destruct(&self._pbnav)

    def value(self):
        return self._pbnav.value()

    def elseBranch(self):
        return new_CN_from_PBNavigator(self._pbnav.elseBranch())

    def thenBranch(self):
        return new_CN_from_PBNavigator(self._pbnav.thenBranch())

    def constant(self):
        return self._pbnav.isConstant()

    def terminalOne(self):
        return self._pbnav.isTerminated()

cdef class BooleanPolynomialVector:
    def __init__(self):
        # This is used by PolyBoRi python code
        PBPolyVector_construct(&self._vec)
        global cur_ring
        self._parent = cur_ring

    def __dealloc__(self):
        PBPolyVector_destruct(&self._vec)

    def __iter__(self):
        return new_BPVI_from_PBPolyVectorIter(self)

    def __len__(self):
        return self._vec.size()

    def __getitem__(self, ind):
        #FIXME: no index checking
        return new_BP_from_PBPoly(self._parent, self._vec.get(ind))

    def append(self, el):
        cdef PBPoly p
        if PY_TYPE_CHECK(el, BooleanPolynomial):
            p = (<BooleanPolynomial>el)._pbpoly
        elif PY_TYPE_CHECK(el, BooleanMonomial):
            PBPoly_construct_pbmonom(&p, (<BooleanMonomial>el)._pbmonom)
        else:
            raise TypeError, "Argument 'el' has incorrect type (expected BooleanPolynomial or BooleanMonomial, got %s)"%(type(el))
        self._vec.push_back(p)

cdef inline BooleanPolynomialVector new_BPV_from_PBPolyVector(\
        BooleanPolynomialRing parent, PBPolyVector juice):
    cdef BooleanPolynomialVector m
    m = <BooleanPolynomialVector>PY_NEW(BooleanPolynomialVector)
    m._vec = juice
    m._parent = parent
    return m

cdef class BooleanPolynomialVectorIterator:
    def __dealloc__(self):
        PBPolyVectorIter_destruct(&self._iter)

    def __iter__(self):
        return self

    def __next__(self):
        cdef PBPoly val
        if PBPolyVectorIter_equal(self._iter, self._obj.end()):
            raise StopIteration

        val = self._iter.value()
        self._iter.next()
        return new_BP_from_PBPoly(self._parent, val)

cdef inline BooleanPolynomialVectorIterator new_BPVI_from_PBPolyVectorIter(\
        BooleanPolynomialVector vec):
    """
    Construct a new BooleanPolynomialVectorIterator
    """
    cdef BooleanPolynomialVectorIterator m
    m = <BooleanPolynomialVectorIterator>PY_NEW(BooleanPolynomialVectorIterator)
    m._parent = vec._parent
    m._obj = vec._vec
    m._iter = m._obj.begin()
    return m


cdef class GroebnerStrategy:
    def __init__(self, param = None):
        if PY_TYPE_CHECK(param, GroebnerStrategy):
            GBStrategy_construct_gbstrategy(&self._strat,
                    (<GroebnerStrategy>param)._strat)
            self._parent = (<GroebnerStrategy>param)._parent
        else:
            GBStrategy_construct(&self._strat)
            global cur_ring
            self._parent = cur_ring

    def __dealloc__(self):
        GBStrategy_destruct(&self._strat)

    def addGeneratorDelayed(self, BooleanPolynomial p):
        if p._pbpoly.isZero():
            raise ValueError, "zero generators not allowed."
        self._strat.addGeneratorDelayed(p._pbpoly)

    def addGenerator(self, BooleanPolynomial p, bint is_impl=False):
        if p._pbpoly.isZero():
            raise ValueError, "zero generators not allowed."
        if self._strat.leadingTerms.owns(p._pbpoly.lead()):
            raise ValueError, "strategy already contains a polynomial with same lead"
        return self._strat.addGenerator(p._pbpoly, is_impl)

    def addAsYouWish(self, BooleanPolynomial p):
        if p._pbpoly.isZero():
            raise ValueError, "zero generators not allowed."
        self._strat.addAsYouWish(p._pbpoly)

    def implications(self, ind):
        implications(self._strat, ind)

    def cleanTopByChainCriterion(self):
        self._strat.cleanTopByChainCriterion()

    def symmGB_F2(self):
        self._strat.symmGB_F2()

    def containsOne(self):
        return self._strat.containsOne()

    def faugereStepDense(self, BooleanPolynomialVector v):
        return new_BPV_from_PBPolyVector(self._parent,
                                    self._strat.faugereStepDense(v._vec))

    def minimalize(self):
        return new_BPV_from_PBPolyVector(self._parent, self._strat.minimalize())

    def minimalizeAndTailReduce(self):
        return new_BPV_from_PBPolyVector(self._parent,
                                    self._strat.minimalizeAndTailReduce())

    def npairs(self):
        return self._strat.npairs()

    def topSugar(self):
        return pairs_top_sugar(self._strat)

    def someSpolysInNextDegree(self, n):
        return new_BPV_from_PBPolyVector(self._parent,
                someNextDegreeSpolys(self._strat, n))

    def allSpolysInNextDegree(self):
        return new_BPV_from_PBPolyVector(self._parent,
                nextDegreeSpolys(self._strat))

    def smallSpolysInNextDegree(self, double f, int n):
        return new_BPV_from_PBPolyVector(self._parent,
                small_next_degree_spolys(self._strat, f, n))

    def llReduceAll(self):
        self._strat.llReduceAll()

    def nextSpoly(self):
        return new_BP_from_PBPoly(self._parent, self._strat.nextSpoly())

    def allGenerators(self):
        return new_BPV_from_PBPolyVector(self._parent,
                self._strat.allGenerators())

    def suggestPluginVariable(self):
        return self._strat.suggestPluginVariable()

    def variableHasValue(self, int idx):
        return self._strat.variableHasValue(idx)

    def nf(self, BooleanPolynomial p):
        return new_BP_from_PBPoly(self._parent, self._strat.nf(p._pbpoly))

    def select(self, BooleanMonomial m):
        return pb_select1(self._strat, m._pbmonom)

    def __len__(self):
        return self._strat.nGenerators()

    def __getitem__(self, int i):
        cdef PBPoly t
        if (i < 0) or (i >= self._strat.nGenerators()):
            raise IndexError
        return new_BP_from_PBPoly(self._parent, GB_get_ith_gen(self._strat, i))

    def __getattr__(self, name):
        if name is 'monomials':
            return new_BS_from_PBSet(self._strat.monomials, self._parent)
        if name is 'llReductor':
            return new_BS_from_PBSet(self._strat.llReductor, self._parent)
        if name is 'optLazy':
            return self._strat.optLazy
        if name is 'optExchange':
            return self._strat.optExchange
        if name is 'optAllowRecursion':
            return self._strat.optAllowRecursion
        if name is 'enabledLog':
            return self._strat.enabledLog
        if name is 'optLL':
            return self._strat.optLL
        if name is 'optLinearAlgebraInLastBlock':
            return self._strat.optLinearAlgebraInLastBlock
        if name is 'redByReduced':
            return self._strat.reduceByTailReduced
        if name is 'reductionSteps':
            return self._strat.reductionSteps
        if name is 'normalForms':
            return self._strat.normalForms
        if name is 'currentDegree':
            return self._strat.currentDegree
        if name is 'chainCriterions':
            return self._strat.chainCriterions
        if name is 'variableChainCriterions':
            return self._strat.variableChainCriterions
        if name is 'easyProductCriterions':
            return self._strat.easyProductCriterions
        if name is 'extendedProductCriterions':
            return self._strat.extendedProductCriterions
        if name is 'averageLength':
            return self._strat.averageLength
        if name is 'optRedTail':
            return self._strat.optRedTail
        if name is 'optDelayNonMinimals':
            return self._strat.optDelayNonMinimals
        if name is 'optBrutalReductions':
            return self._strat.optBrutalReductions
        if name is 'optRedTailDegGrowth':
            return self._strat.optRedTailDegGrowth
        if name is 'optStepBounded':
            return self._strat.optStepBounded
        if name is 'optRedTailInLastBlock':
            return self._strat.optRedTailInLastBlock,
        if name is 'minimalLeadingTerms':
            return new_BS_from_PBSet(self._strat.minimalLeadingTerms,
                    self._parent)
        if name is 'leadingTerms':
            return new_BS_from_PBSet(self._strat.leadingTerms, self._parent)
        raise AttributeError, name

    def __setattr__(self, name, val):
        if name is 'enabledLog':
            self._strat.enabledLog = val
        elif name == 'optRedTail':
            self._strat.optRedTail = val
        elif name is 'optLazy':
            self._strat.optLazy = val
        elif name is 'optLL':
            self._strat.optLL = val
        elif name is 'optDelayNonMinimals':
            self._strat.optDelayNonMinimals = val
        elif name is 'optBrutalReductions':
            self._strat.optBrutalReductions = val
        elif name is 'optExchange':
            self._strat.optExchange = val
        elif name is 'optAllowRecursion':
            self._strat.optAllowRecursion = val
        elif name is 'optRedTailDegGrowth':
            self._strat.optRedTailDegGrowth = val
        elif name is 'optStepBounded':
            self._strat.optStepBounded = val
        elif name is 'optLinearAlgebraInLastBlock':
            self._strat.optLinearAlgebraInLastBlock = val
        elif name is 'optRedTailInLastBlock':
            self._strat.optRedTailInLastBlock = val
        elif name is 'redByReduced':
            self._strat.reduceByTailReduced = val
        else:
            raise AttributeError, name

cdef inline CCuddNavigator new_CN_from_PBNavigator(PBNavigator juice):
    """
    Construct a new CCuddNavigator
    """
    cdef CCuddNavigator n
    n = <CCuddNavigator>PY_NEW(CCuddNavigator)
    n._pbnav = juice
    return n

cdef class VariableBlock_base:
    def __init__(self, size, start_index, offset):
        self.size = size
        self.start_index = start_index
        self.offset = offset

cdef class VariableBlockTrue(VariableBlock_base):
    def __init__(self, size, start_index, offset):
        global cur_ring
        self.ring = cur_ring
        VariableBlock_base.__init__(self, size, start_index, offset)

    def __call__(self, int i):
        #FIXME: no index checking
        cdef PBVar v
        PBVar_construct_int(&v, self.offset+self.start_index+self.size-1-i)
        return new_BM_from_PBVar(self.ring._monom_monoid, self.ring, v)

cdef class VariableBlockFalse(VariableBlock_base):
    def __init__(self, size, start_index, offset):
        global cur_ring
        self.ring = cur_ring
        VariableBlock_base.__init__(self, size, start_index, offset)

    def __call__(self, int i):
        #FIXME: no index checking
        cdef PBVar v
        PBVar_construct_int(&v, i-self.start_index+self.offset)
        return new_BM_from_PBVar(self.ring._monom_monoid, self.ring, v)

def VariableBlock(size, start_index, offset, reverse):
    if reverse:
        return VariableBlockTrue(size, start_index, offset)
    else:
        return VariableBlockFalse(size, start_index, offset)


def add_up_polynomials(BooleanPolynomialVector v):
    return new_BP_from_PBPoly(v._parent, pb_add_up_polynomials(v._vec))

def nf3(GroebnerStrategy s, BooleanPolynomial p, BooleanMonomial m):
    return new_BP_from_PBPoly(s._parent,
            pb_nf3(s._strat, p._pbpoly, m._pbmonom))

def red_tail(GroebnerStrategy s, BooleanPolynomial p):
    return new_BP_from_PBPoly(p._parent, pb_red_tail(s._strat, p._pbpoly))

def map_every_x_to_x_plus_one(BooleanPolynomial p):
    return new_BP_from_PBPoly(p._parent,
            pb_map_every_x_to_x_plus_one(p._pbpoly))

def zeroes(pol, BooleSet s):
    cdef PBPoly p
    if PY_TYPE_CHECK(pol, BooleanPolynomial):
        p = (<BooleanPolynomial>pol)._pbpoly
    elif PY_TYPE_CHECK(pol, BooleanMonomial):
        PBPoly_construct_pbmonom(&p, (<BooleanMonomial>pol)._pbmonom)
    else:
        raise TypeError, "Argument 'p' has incorrect type (expected BooleanPolynomial or BooleanMonomial, got %s)"%(type(pol))
    return new_BS_from_PBSet(pb_zeroes(p, s._pbset), s.ring)

def interpolate(zero, one):
    cdef PBSet z, o
    cdef BooleanPolynomialRing ring
    if PY_TYPE_CHECK(zero, BooleSet):
        z = (<BooleSet>zero)._pbset
        ring = (<BooleSet>zero).ring
    elif PY_TYPE_CHECK(zero, BooleanPolynomial):
        z = (<BooleanPolynomial>zero)._pbpoly.set()
        ring = (<BooleanPolynomial>zero)._parent
    else:
        raise TypeError, "Argument 'zero' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(zero))
    if PY_TYPE_CHECK(one, BooleSet):
        o = (<BooleSet>one)._pbset
    elif PY_TYPE_CHECK(one, BooleanPolynomial):
        o = (<BooleanPolynomial>one)._pbpoly.set()
    else:
        raise TypeError, "Argument 'one' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(one))
    return new_BP_from_PBPoly(ring, pb_interpolate(z, o))

def interpolate_smallest_lex(zero, one):
    cdef PBSet z, o
    cdef BooleanPolynomialRing ring
    if PY_TYPE_CHECK(zero, BooleSet):
        z = (<BooleSet>zero)._pbset
        ring = (<BooleSet>zero).ring
    elif PY_TYPE_CHECK(zero, BooleanPolynomial):
        z = (<BooleanPolynomial>zero)._pbpoly.set()
        ring = (<BooleanPolynomial>zero)._parent
    else:
        raise TypeError, "Argument 'zero' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(zero))
    if PY_TYPE_CHECK(one, BooleSet):
        o = (<BooleSet>one)._pbset
    elif PY_TYPE_CHECK(one, BooleanPolynomial):
        o = (<BooleanPolynomial>one)._pbpoly.set()
    else:
        raise TypeError, "Argument 'one' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(one))
    return new_BP_from_PBPoly(ring, pb_interpolate_smallest_lex(z, o))

def contained_vars(BooleSet m):
    return new_BS_from_PBSet(pb_contained_variables_cudd_style(m._pbset),
            m.ring)

def mod_var_set(BooleSet a, BooleSet v):
    return new_BS_from_PBSet(pb_mod_var_set(a._pbset, v._pbset), a.ring)

def mult_fact_sim_C(BooleanPolynomialVector v):
    return new_BP_from_PBPoly(v._parent, pb_mult_fast_sim(v._vec))

def recursively_insert(CCuddNavigator n, int ind, BooleSet m):
    cdef PBSet b
    b = pb_recursively_insert((<CCuddNavigator>n)._pbnav, ind,
                                                (<BooleSet>m)._pbset)
    return new_BS_from_PBSet(b, m.ring)

def ll_red_nf(p, BooleSet reductors):
    cdef PBPoly t
    cdef PBPoly res
    cdef BooleanPolynomialRing parent
    if PY_TYPE_CHECK(p, BooleSet):
        PBPoly_construct_pbset(&t, (<BooleSet>p)._pbset)
        parent = (<BooleSet>p).ring
    elif PY_TYPE_CHECK(p, BooleanPolynomial):
        t = (<BooleanPolynomial>p)._pbpoly
        parent = (<BooleanPolynomial>p)._parent
    else:
        raise TypeError, "Argument 'p' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(p))

    res = pb_ll_red_nf(t, reductors._pbset)

    return new_BP_from_PBPoly(parent, res)

def ll_red_nf_noredsb(BooleanPolynomial p, BooleSet reductors):
    cdef PBPoly t
    t = pb_ll_red_nf_noredsb(p._pbpoly, reductors._pbset)
    return new_BP_from_PBPoly(p._parent, t)

def mod_mon_set(BooleSet as, BooleSet vs):
    cdef PBSet b
    b = pb_mod_mon_set((<BooleSet>as)._pbset, (<BooleSet>vs)._pbset)
    return new_BS_from_PBSet(b, as.ring)

def get_order_code():
    return pbenv_getOrderCode()

def change_ordering(order):
    global cur_ring
    cur_ring._change_ordering(order)

def parallel_reduce(BooleanPolynomialVector inp, GroebnerStrategy strat, \
                                    int average_steps, double delay_f):
    return new_BPV_from_PBPolyVector(inp._parent, \
        pb_parallel_reduce(inp._vec, strat._strat, average_steps, delay_f))

def have_degree_order():
    return pbenv_isDegreeOrder()

def set_variable_name( i, s):
    global cur_ring
    cur_ring._set_variable_name(i,s)

def append_ring_block(i):
    pb_append_block(i)

def if_then_else(int ind, a, b):
    cdef PBSet a_set, b_set
    cdef PBSet res
    cdef BooleanPolynomialRing ring
    if PY_TYPE_CHECK(b, BooleSet):
        b_set = (<BooleSet>b)._pbset
        ring = (<BooleSet>b).ring
    elif PY_TYPE_CHECK(b, BooleanPolynomial):
        b_set = (<BooleanPolynomial>b)._pbpoly.set()
        ring = (<BooleanPolynomial>b)._parent
    else:
        raise TypeError, "Argument 'b' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(b))

    if PY_TYPE_CHECK(a, BooleSet):
        a_set = (<BooleSet>a)._pbset
    elif PY_TYPE_CHECK(a, BooleanPolynomial):
        a_set = (<BooleanPolynomial>a)._pbpoly.set()
    else:
        raise TypeError, "Argument 'a' has incorrect type (expected BooleSet or BooleanPolynomial, got %s)"%(type(a))

    if ind >= a_set.navigation().value() or ind >= b_set.navigation().value():
        raise IndexError, "value of ind must be less than the values of roots of the branches."
    PBSet_construct_indsetset(&res, ind, a_set.navigation(),
            b_set.navigation(), ring._pbring);
    return new_BS_from_PBSet(res, ring)

def top_index(s):
    if PY_TYPE_CHECK(s, BooleSet):
        return (<BooleSet>s)._pbset.navigation().value()
    elif PY_TYPE_CHECK(s, BooleanMonomial):
        return (<BooleanMonomial>s)._pbmonom.firstIndex()
    elif PY_TYPE_CHECK(s, BooleanPolynomial):
        return (<BooleanPolynomial>s)._pbpoly.navigation().value()
    else:
        raise TypeError, "Argument 's' has incorrect type (expected BooleSet, BooleanMonomial or BooleanPolynomial, got %s)"%(type(s))

def get_cring():
    """
    Return the currently active global ring, this is only relevant for
    the native \PolyBoRi interface.

    sage: from polybori import *
    sage: R = declare_ring([Block('x',2),Block('y',3)],globals())
    sage: Q = get_cring(); Q
    Boolean PolynomialRing in x(0), x(1), y(0), y(1), y(2)
    sage: R is Q
    True
    """
    global cur_ring
    return cur_ring

def set_cring(BooleanPolynomialRing R):
    """
    Set the currently active global ring, this is only relevant for the native
    \PolyBoRi interface.

    sage: from polybori import *
    sage: declare_ring([Block('x',2),Block('y',3)],globals())
    Boolean PolynomialRing in x(0), x(1), y(0), y(1), y(2)
    sage: R = get_cring(); R
    Boolean PolynomialRing in x(0), x(1), y(0), y(1), y(2)

    sage: declare_ring([Block('x',2),Block('y',2)],globals())
    Boolean PolynomialRing in x(0), x(1), y(0), y(1)

    sage: get_cring()
    Boolean PolynomialRing in x(0), x(1), y(0), y(1)

    sage: set_cring(R)
    sage: get_cring()
    Boolean PolynomialRing in x(0), x(1), y(0), y(1), y(2)
    """
    global cur_ring
    cur_ring = R

def unpickle_BooleanPolynomial(ring, string):
    """
    Unpickle boolean polynomials

    EXAMPLE:
        sage: T = TermOrder('deglex',2)+TermOrder('deglex',2)
        sage: P.<a,b,c,d> = BooleanPolynomialRing(4,order=T)
        sage: loads(dumps(a+b)) == a+b # indirect doctest
        True
    """
    return ring(eval(string,ring.gens_dict()))

def unpickle_BooleanPolynomialRing(n, names, order):
    """
    Unpickle boolean polynomial rings.

    EXAMPLE:
        sage: T = TermOrder('deglex',2)+TermOrder('deglex',2)
        sage: P.<a,b,c,d> = BooleanPolynomialRing(4,order=T)
        sage: loads(dumps(P)) == P  # indirect doctest
        True
    """
    return BooleanPolynomialRing(n, names=names, order=order)

###
#
# M4RI, this needs to be revisited as soon as PolyBoRi uses the shared
# library version of M4RI.
#
###

cdef int M4RI_init = 0

cdef init_M4RI():
    global M4RI_init
    if M4RI_init is int(0):
        buildAllCodes()
        setupPackingMasks()
        M4RI_init = 1

def free_m4ri():
    destroyAllCodes()

init_M4RI()

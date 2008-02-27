r"""
Boolean Polynomials via PolyBoRi.

We call boolean polynomials elements of the quotient ring

    $F_2[x_1,...,x_n]/<x_1^2+x_1,...,x_n^2+x_n>$.

AUTHOR:
    -- Burcin Erocal <burcin@erocal.org>
    -- Martin Albrecht <malb@informatik.uni-bremen.de>

REFERENCES:
    Michael Brickenstein, Alexander Dreyer; 'POLYBORI: A Groebner basis
    framework for Boolean polynomials';
    http://www.itwm.fraunhofer.de/zentral/download/berichte/bericht122.pdf

"""

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
include '../../libs/polybori/decl.pxi'

from sage.structure.element cimport Element
from sage.structure.element cimport RingElement
from sage.structure.element cimport ModuleElement

from sage.structure.parent cimport Parent

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.integer import Integer

from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.finite_field import GF
from sage.monoids.monoid import Monoid_class

from sage.structure.sequence import Sequence

order_dict= {"lp":      lp,
             "dlex":    dlex,
             "dp_asc":  dp_asc,
             "block_dlex":   block_dlex,
             "block_dp_asc": block_dp_asc,
             }

order_mapping = {'lp':   lp,
                 'Dp':   dlex,
                 'dp':   dp_asc}

cdef class BooleanPolynomialRing(MPolynomialRing_generic):
    """
    Boolean Polynomial Ring.
    """
    def __init__(self, n, names, order='lex'):
        """
        Construct a BooleanPolynomialRing with the following parameters:

        INPUT:
            n -- number of variables (an integer > 1)
            names -- names of ring variables, may be a string of list/tuple
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
        """
        try:
            n = int(n)
        except TypeError, msg:
            raise TypeError, "Number of variables must be an integer"

        if n < 1:
            raise ValueError, "Number of variables must be greater than 1."

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

        PBRing_construct(&self._pbring, n, pb_order_code)

        MPolynomialRing_generic.__init__(self, GF(2), n, names, order)

        counter = 0
        for i in range(len(order.blocks)-1):
            counter += order.blocks[i][1]
            pb_append_ring_block(counter)

        for i in range(self.ngens()):
            _n = self._names[i]
            self._pbring.setRingVariableName(i,_n)

        self._zero_element = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._zero_element)._pbpoly, 0)
        self._one_element  = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._one_element)._pbpoly, 1)

        self._monom_monoid = BooleanMonomialMonoid(self)

        set_cring(self)

    def __dealloc__(self):
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

    def gen(self, int n=0):
        """
        Returns the n-th generator of self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.gen()
            x
            sage: P.gen(2)
            z
        """
        return new_BP_from_DD(self, self._pbring.variable(n))

    def gens(self):
        """
        Return the tuple of variables in self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: P.gens()
            (x, y, z)

            sage: P = BooleanPolynomialRing(10,'x')
            sage: P.gens()
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
        """
        return tuple([new_BP_from_DD(self, self._pbring.variable(i)) \
                for i in xrange(self.ngens())])

    def _repr_(self):
        """
        EXAMPLE:
            sage: P.<x, y> = BooleanPolynomialRing(2)
            sage: P
            Boolean PolynomialRing in x, y
        """
        self._pbring.activate()
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
        if PY_TYPE_CHECK(other, int) or PY_TYPE_CHECK(other, Integer):
            if other %2:
                return self._one_element
            else:
                return self._zero_element
        elif PY_TYPE_CHECK(other, BooleanMonomial):
            if (<BooleanMonomial>other)._parent._ring is self:
                p = new_BP_from_PBMonom(self, (<BooleanMonomial>other)._pbmonom)
                return p
            elif (<BooleanMonomial>other)._parent.ngens() <= \
                    self._pbring.nVariables():
                try:
                    var_mapping = get_var_mapping(self, other.parent())
                except NameError, msg:
                    raise ValueError, "cannot coerce monomial %s to %s: %s"%(other,self,msg)
                p = self._one_element
                for i in other:
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
                        for i in monom:
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
        Convert elements of other objects to self.

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

        # we check for other PolyBoRi types first since this conversion
        # is used by the PolyBoRi python code often
        if PY_TYPE_CHECK(other, DD):
            return new_BP_from_DD(self, (<DD>other)._pbdd)
        elif PY_TYPE_CHECK(other, BooleSet):
            return new_BP_from_PBSet(self, (<BooleSet>other)._pbset)

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
                for i in other:
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
                        for i in monom:
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
        return (<Parent>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
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
        if PY_TYPE_CHECK(right, BooleanPolynomialRing):
            return cmp( (type(left), map(str, left.gens()), left.term_order()),
                    (type(right), map(str, right.gens()), right.term_order()))
        else:
            return -1

    def __hash__(self):
        """
        Return a hash of self.

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
        Create an ideal of this ring.

        INPUT:
            gens -- list or tuple of generators
            coerce -- bool (default: True) automatically coerce the given polynomials to this ring to form the ideal

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
        Return a random boolean polynomial. Generated polynomial has the given number of terms, and at most given degree.

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
        """
        Recursively generate a random polynomial in self, using the variables from \code{vars_set}.

        INPUT:
            degree -- maximum degree
            monom_counts -- a list containing total number of monomials up to given degree
            vars_set -- list of variable indicies to use in the generated polynomial
            dfirst -- if \code{True} choose degree first, otherwise choose the monomial uniformly
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
        """
        Choose a random monomial uniformly from set of monomials in the variables indexed by \code{vars_set} in self.

        INPUT:
            monom_counts -- list of number of monomials up to given degree
            vars_set -- list of variable indicies to use in the generated monomial

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
        """
        Choose a random monomial using variables indexed in \code{vars_set} up to given \code{degree}. The degree of the monomial, $d$, is chosen uniformly in the interval [0,degree] first, then the monomial is generated by selecting a random sample of size $d$ from \code{vars_set}.

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

def get_var_mapping(ring, other):
    """
    Return a variable mapping between variables of other and ring. When other is a parent object, the mapping defines images for all variables of other. If it is an element, only variables occuring in other are mapped.

    Raises NameError if no such mapping is possible.
    """
    my_names = list(ring._names) # we need .index(.)
    if PY_TYPE_CHECK(other, ParentWithGens):
        vars = range(other.ngens())
        ovar_names = other._names
    else:
        ovar_names = other.parent().variable_names()
        if PY_TYPE_CHECK(other, BooleanPolynomial):
            vars = other.vars()
        elif PY_TYPE_CHECK(other, BooleanMonomial):
            vars = other
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
    """
    Monomial Monoid of a Boolean Polynomial Ring

    Provides a parent object for BooleanMonomials.
    """
    def __init__(self, BooleanPolynomialRing polring):
        """
        Construct a BooleanMonomialMonoid given a BooleanPolynomialRing

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

        m = new_BM(self)
        PBMonom_construct(&m._pbmonom)
        self._one_element = m

    def _repr_(self):
        return "MonomialMonoid of %s" % (str(self._ring))

    def __hash__(self):
        """
        Return a hash for self.
        """
        return hash(str(self))

    def ngens(self):
        """
        Returns the number of variables in self.

        EXAMPLES:
            sage: P = BooleanPolynomialRing(100, 'x')
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M.ngens()
            100
        """
        return self._ring.ngens()

    def gen(self, int n=0):
        """
        Return the n-th generator of self.

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
        return new_BM_from_DD(self,
                (<BooleanPolynomialRing>self._ring)._pbring.variable(n))

    def gens(self):
        """
        Return the tuple of generators of self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M.gens()
            (x, y, z)
        """
        return tuple([new_BM_from_DD(self,
            (<BooleanPolynomialRing>self._ring)._pbring.variable(i)) \
                for i in xrange(self.ngens())])

    def _coerce_impl(self, other):
        """
        Canonical conversion of elements from other objects to self.

        EXAMPLES:

        Coerce elements of self.

            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: x_monom = M(x); x_monom
            x
            sage: M._coerce_(x_monom)
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
                for i in other:
                    m *= var_mapping[i]
                return m
        raise TypeError, "coercion from %s to %s not implemented" % \
            (type(other), str(self))

    def __call__(self, other = None):
        """
        Convert elements of other objects to elements of self.

        INPUT:
            other -- element to convert,
                     if None a BooleanMonomial representing 1 is returned
                     only BooleanPolynomials with the same parent ring as self which have a single monomial is converted

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
    """
    BooleanMonomial
    """
    def __init__(self, parent):
        """
        Construct a BooleanMonomial object.

        INPUT:
            parent -- parent monoid this element lies in

        """
        PBMonom_construct(&self._pbmonom)
        _parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBMonom_destruct(&self._pbmonom)

    def __richcmp__(left, right, int op):
        # boilerplate code from sage.structure.parent
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
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
        cdef comparecodes res
        res = left._pbmonom.compare((<BooleanMonomial>right)._pbmonom)
        return res

    def _repr_(self):
        """
        Return a string representing self.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M(x*y)
            x*y

            sage: R.<t,u> = BooleanPolynomialRing(2)
            sage: M(x*y)
            x*y
        """
        (<BooleanPolynomialRing>self._parent._ring)._pbring.activate()
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
        for i in self:
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
            sage: m = f.monomials()[0]
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
        for var in self:
            res *= d[var]
        return res

    def __hash__(self):
        """
        Return a hash of self.
        """
        return self._pbmonom.hash()

    def deg(BooleanMonomial self):
        """
        Return total degree of self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: M(x*y).deg()
            2

            sage: M(x*x*y*z).deg()
            3
        """
        return self._pbmonom.deg()

    def __len__(BooleanMonomial self):
        """
        Return number of variables in self. This is equivalent to the total degree for BooleanMonomials.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: len(M(x*y))
            2
        """
        return self._pbmonom.deg()

    def __iter__(self):
        """
        Return an iterator over the indicies of the variables in self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = sage.rings.polynomial.pbori.BooleanMonomialMonoid(P)
            sage: list(iter(M(x*z)))
            [0, 2]
        """
        return new_BMI_from_PBMonomIter(self._pbmonom, self._pbmonom.begin())

    cdef MonoidElement _mul_c_impl(left, MonoidElement right):
        """
        Multiply self with another BooleanMonomial.

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
                (<BooleanMonomial>left)._parent, (<BooleanMonomial>left)._pbmonom)
        m._pbmonom.imul( (<BooleanMonomial>right)._pbmonom )
        return m

    def __add__(left, right):
        """
        Addition operator. Returns a BooleanPolynomial.

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

        res = new_BP_from_PBMonom(monom._parent._ring, monom._pbmonom)
        return res.__iadd__(other)


cdef inline BooleanMonomial new_BM(parent):
    """
    Construct a new BooleanMonomial
    """
    cdef BooleanMonomial m
    m = <BooleanMonomial>PY_NEW(BooleanMonomial)
    m._parent = parent
    return m

cdef inline BooleanMonomial new_BM_from_PBMonom(parent, PBMonom juice):
    cdef BooleanMonomial m = new_BM(parent)
    PBMonom_construct_pbmonom(&m._pbmonom,juice)
    return m

cdef inline BooleanMonomial new_BM_from_PBVar(parent, PBVar juice):
    cdef BooleanMonomial m = new_BM(parent)
    PBMonom_construct_pbvar(&m._pbmonom,juice)
    return m

cdef inline BooleanMonomial new_BM_from_DD(parent, PBDD juice):
    cdef BooleanMonomial m = new_BM(parent)
    PBMonom_construct_dd(&m._pbmonom,juice)
    return m

cdef class BooleanMonomialIterator:
    def __iter__(self):
        return self

    def __dealloc__(self):
        PBMonomIter_destruct(&self._iter)

    def __next__(self):
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
    """
    BooleanPolynomial
    """
    def __init__(self, parent):
        """
        Construct a BooleanPolynomial object in the given BooleanPolynomialRing.
        """
        PBPoly_construct(&self._pbpoly)
        self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBPoly_destruct(&self._pbpoly)

    def _repr_(self):
        """
        Return a string representing self.
        """
        (<BooleanPolynomialRing>self._parent)._pbring.activate()
        return PBPoly_to_str(&self._pbpoly)

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        cdef BooleanPolynomial p = new_BP_from_PBPoly(\
                (<BooleanPolynomial>left)._parent, (<BooleanPolynomial>left)._pbpoly)
        p._pbpoly.iadd( (<BooleanPolynomial>right)._pbpoly )
        return p

    cdef ModuleElement _sub_c_impl(left, ModuleElement right):
        return left._add_c_impl(right)

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        if left:
            return new_BP_from_PBPoly(left._parent, self._pbpoly)
        else:
            return 0

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        return self._rmul_c_impl(right)

    cdef RingElement _mul_c_impl(left, RingElement right):
        cdef BooleanPolynomial p = new_BP_from_PBPoly(\
                (<BooleanPolynomial>left)._parent, (<BooleanPolynomial>left)._pbpoly)
        p._pbpoly.imul( (<BooleanPolynomial>right)._pbpoly )
        return p

    def is_equal(self, BooleanPolynomial right):
        return self._pbpoly.is_equal(right._pbpoly)

    def __richcmp__(left, right, int op):
        #boilerplate from sage.structure.element
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
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
        cdef int res
        from itertools import izip
        for lm, rm in izip(left, right):
            res = cmp(lm, rm)
            if res != 0:
                return res
        return cmp(len(left),len(right))

    def __iter__(self):
        """
        Return an iterator over the monomials of self, in the order of the parent ring.

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
            [x*y*z, y*z, x*y, z, x]

        TESTS:
            sage: R = BooleanPolynomialRing(1,'y')
            sage: list(iter(y))
            [y]
            sage: R
            Boolean PolynomialRing in y
        """
        (<BooleanPolynomialRing>self._parent)._pbring.activate()
        return new_BPI_from_PBPolyIter(self, self._pbpoly.orderedBegin())

    def __pow__(BooleanPolynomial self, int exp, ignored):
        """
        Return self^(exp).

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
        return self._pbpoly.deg()

    def lm(BooleanPolynomial self):
        """
        Return the leading monomial of self, according to the order of parent ring.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x+y+y*z).lm()
            x

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: (x+y+y*z).lm()
            y*z
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._pbpoly.lead())

    def lm_lex(BooleanPolynomial self):
        """
        Return the leading monomial of self in lexicographic order.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: (x+y+y*z).lm_lex()
            x

            sage: P.<x,y,z> = BooleanPolynomialRing(3, order='deglex')
            sage: (x+y+y*z).lm_lex()
            x
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid,
                                                self._pbpoly.lexLead())

    def lt(BooleanPolynomial self):
        """
        Return the leading term of self, according to the order of the parent ring.

        Note that for boolean polynomials this is equivalent to returning leading monomials.

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
        """
        Check if self is zero.

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
        """
        Return a BooleanMonomial with all the variables appearing in self.

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
        (<BooleanPolynomialRing>self._parent)._pbring.activate()
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._pbpoly.usedVariables())

    def elimination_length(self):
        return self._pbpoly.eliminationLength()

    def variables(self):
        """
        Return a list of all variables appearing in self.

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
        V = self.vars()
        return [P.gen(i) for i in V]

    def monomials(self):
        """
        Return a list of monomials appearing in self ordered largest
        to smallest.

        EXAMPLE:
            sage: P.<a,b,c> = BooleanPolynomialRing(3,order='lex')
            sage: f = a + c*b
            sage: f.monomials()
            [a, b*c]

            sage: P.<a,b,c> = BooleanPolynomialRing(3,order='degrevlex')
            sage: f = a + c*b
            sage: f.monomials()
            [b*c, a]
        """
        return list(self)

    def __hash__(self):
        """
        Return hash for self.

        EXAMPLE:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: {x:1} # indirect doctest
            {x: 1}
        """
        (<BooleanPolynomialRing>self._parent)._pbring.activate()
        return hash(PBPoly_to_str(&self._pbpoly))

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
                        fixed[list(list(var)[0])[0]] = val
                except TypeError:
                    fixed[list(list(var)[0])[0]] = val

        for var,val in kwds.iteritems():
            var = P(var)
            try:
                v =  P(val)
                if v.constant():
                    self = ll_red_nf(self, (var + v).set())
                else:
                    fixed[list(list(var)[0])[0]] = val
            except TypeError:
                fixed[list(list(var)[0])[0]] = val

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
        (<BooleanPolynomialRing>self._parent)._pbring.activate()
        return unpickle_BooleanPolynomial, (self._parent, PBPoly_to_str(&self._pbpoly))

    def set(self):
        return new_BS_from_PBSet(self._pbpoly.set())

    def deg(self):
        return self._pbpoly.deg()

    def elength(self):
        return self._pbpoly.eliminationLength()

    def lead(self):
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._pbpoly.lead())

    def lexLead(self):
        return new_BM_from_PBMonom(self._parent._monom_monoid,
                                                self._pbpoly.lexLead())
    def lexLmDeg(self):
        """
        Return degree of leading monomial w.r.t to lex ordering
        """
        return self._pbpoly.lexLmDeg()

    def constant(self):
        return self._pbpoly.isConstant()

    def lmDeg(self):
        return self._pbpoly.lmDeg()

    def isZero(self):
        return self._pbpoly.isZero()

    def isOne(self):
        return self._pbpoly.isOne()

    def navigation(self):
        return new_CN_from_PBNavigator(self._pbpoly.navigation())

    def mapEveryXToXPlusOne(self):
        return new_BP_from_PBPoly(self._parent, map_every_x_to_x_plus_one(self._pbpoly))

cdef class BooleanPolynomialIterator:
    def __iter__(self):
        return self

    def __dealloc__(self):
        PBPolyIter_destruct(&self._iter)

    def __next__(self):
        (<BooleanPolynomialRing>self._obj._parent)._pbring.activate()
        cdef PBMonom val
        if self._iter.equal(self._obj._pbpoly.orderedEnd()):
            raise StopIteration
        val = self._iter.value()
        self._iter.next()
        return new_BM_from_PBMonom(self._obj._parent._monom_monoid, val)

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
    """
    BooleanPolynomialIdeal
    """
    def __init__(self, ring, gens=[], coerce=True):
        """
        Construct a BooleanPolynomialIdeal object.

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
            faugere -- use Faugere's F4 (default: False)
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

cdef inline BooleanPolynomial new_BP(BooleanPolynomialRing parent):
    """
    Construct a new BooleanPolynomial
    """
    cdef BooleanPolynomial p
    p = <BooleanPolynomial>PY_NEW(BooleanPolynomial)
    p._parent = parent
    return p

cdef inline BooleanPolynomial new_BP_from_DD(BooleanPolynomialRing parent,
        PBDD juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_dd(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBPoly(BooleanPolynomialRing parent,
        PBPoly juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbpoly(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBMonom(BooleanPolynomialRing parent,
        PBMonom juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbmonom(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBSet(BooleanPolynomialRing parent,
        PBSet juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbset(&p._pbpoly,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_int(BooleanPolynomialRing parent,
        int juice):
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
    def __init__(self, param = None):
        if PY_TYPE_CHECK(param, DD):
            PBSet_construct_dd(&self._pbset, (<DD>param)._pbdd)
        elif PY_TYPE_CHECK(param, CCuddNavigator):
            PBSet_construct_pbnav(&self._pbset, (<CCuddNavigator>param)._pbnav)
        elif PY_TYPE_CHECK(param, BooleSet):
            PBSet_construct_pbset(&self._pbset, (<BooleSet>param)._pbset)
        else:
            PBSet_construct(&self._pbset)

    def __dealloc__(self):
        PBSet_destruct(&self._pbset)

    def __call__(self):
        return self

    def _repr_(self):
        return PBSet_to_str(&self._pbset)

    def empty(self):
        return self._pbset.emptiness()

    def navigation(self):
        return new_CN_from_PBNavigator(self._pbset.navigation())

    def cartesianProduct(self, rhs):
        return new_BS_from_PBSet(self._pbset.cartesianProduct((<BooleSet>rhs)._pbset))

    def diff(self, rhs):
        return new_BS_from_PBSet(self._pbset.diff((<BooleSet>rhs)._pbset))

    def change(self, ind):
        return new_BS_from_PBSet(self._pbset.change(ind))

    def usedVariables(self):
        return new_BM_from_PBMonom(get_cring()._monom_monoid,
                                            self._pbset.usedVariables())

    def if_then_else(self, int ind, BooleSet a, BooleSet b):
        cdef PBSet res
        if ind >= a.navigation().value() or ind >= b.navigation().value():
            raise IndexError, "value of ind must be less than the values of roots of the branches."
        PBSet_construct_indsetset(&res, ind, a._pbset.navigation(),
                b._pbset.navigation());
        return new_BS_from_PBSet(res)

    def union(self, rhs):
        return new_BS_from_PBSet(self._pbset.unite((<BooleSet>rhs)._pbset))

    def __iter__(self):
        return new_BSI_from_PBSetIter(self, get_cring())

    def subset0(self, i):
        return new_BS_from_PBSet(self._pbset.subset0(int(i)))

    def subset1(self, i):
        return new_BS_from_PBSet(self._pbset.subset1(int(i)))




cdef inline BooleSet new_BS_from_PBSet(PBSet juice):
    """
    Construct a new BooleSet
    """
    cdef BooleSet s
    s = <BooleSet>PY_NEW(BooleSet)
    s._pbset = juice
    return s

cdef inline BooleSet new_BS_from_PBDD(PBDD juice):
    """
    Construct a new BooleSet
    """
    cdef BooleSet s
    s = <BooleSet>PY_NEW(BooleSet)
    PBSet_construct_dd(&s._pbset, juice)
    return s

cdef class BooleSetIterator:
    def __iter__(self):
        return self

    def __dealloc__(self):
        PBSetIter_destruct(&self._iter)

    def __next__(self):
        cdef PBMonom val
        if self._iter.equal(self._obj._pbset.end()):
            raise StopIteration
        val = self._iter.value()
        self._iter.next()
        return new_BM_from_PBMonom(self._ring._monom_monoid, val)

cdef inline BooleSetIterator new_BSI_from_PBSetIter(\
                BooleSet parent, BooleanPolynomialRing cring):
    """
    Construct a new BooleSetIterator
    """
    cdef BooleSetIterator m
    m = <BooleSetIterator>PY_NEW(BooleSetIterator)
    m._obj = parent
    m._ring = cring
    m._iter = parent._pbset.begin()
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
        PBPolyVector_construct(&self._vec)
        self._parent = get_cring()

    def __dealloc__(self):
        PBPolyVector_destruct(&self._vec)

    def __iter__(self):
        #return new_BPVI_from_PBPolyVectorIter(self._parent, self._vec,
        #    self._vec.begin())
        return BooleanPolynomialVectorIterator(self._parent, self)

    def __len__(self):
        return self._vec.size()

    def __getitem__(self, ind):
        #FIXME: no index checking
        return new_BP_from_PBPoly(self._parent, self._vec.get(ind))

    def append(self, BooleanPolynomial poly):
        self._vec.push_back(poly._pbpoly)

cdef inline BooleanPolynomialVector new_BPV_from_PBPolyVector(\
        BooleanPolynomialRing parent, PBPolyVector juice):
    cdef BooleanPolynomialVector m
    m = <BooleanPolynomialVector>PY_NEW(BooleanPolynomialVector)
    m._vec = juice
    m._parent = parent
    return m

cdef class BooleanPolynomialVectorIterator:
    def __init__(self, BooleanPolynomialRing parent,
                                    BooleanPolynomialVector vector):
        self._parent = parent
        self._obj = vector._vec
        self._iter = self._obj.begin()

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
        BooleanPolynomialRing parent, PBPolyVector juice, PBPolyVectorIter itr):
    """
    Construct a new BooleanMonomialIterator
    """
    cdef BooleanPolynomialVectorIterator m
    m = <BooleanPolynomialVectorIterator>PY_NEW(BooleanPolynomialVectorIterator)
    m._obj = juice
    m._iter = itr
    m._parent = parent
    return m


cdef class GroebnerStrategy:
    def __init__(self, param = None):
        self._parent = get_cring()
        if PY_TYPE_CHECK(param, GroebnerStrategy):
            GBStrategy_construct_gbstrategy(&self._strat,
                    (<GroebnerStrategy>param)._strat)
        else:
            GBStrategy_construct(&self._strat)

    def __dealloc__(self):
        GBStrategy_destruct(&self._strat)

    def addGeneratorDelayed(self, BooleanPolynomial p):
        if p.isZero():
            raise ValueError, "zero generators not allowed."
        self._strat.addGeneratorDelayed(p._pbpoly)

    def addGenerator(self, BooleanPolynomial p, bint is_impl=False):
        if p.isZero():
            raise ValueError, "zero generators not allowed."
        if self._strat.leadingTerms.owns(p._pbpoly.lead()):
            raise ValueError, "strategy already contains a polynomial with same lead"
        return self._strat.addGenerator(p._pbpoly, is_impl)

    def addAsYouWish(self, BooleanPolynomial p):
        if p.isZero():
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
        cdef PBPolyVector v = someNextDegreeSpolys(self._strat, n)
        return new_BPV_from_PBPolyVector(self._parent, v)

    def __len__(self):
        return self._strat.nGenerators()

    def __getitem__(self, int i):
        #FIXME: no index checking
        return new_BP_from_PBPoly(self._parent, GB_get_ith_gen(self._strat, i))

    def __getattr__(self, name):
        if name is 'enabledLog':
            return self._strat.enabledLog
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
        if name is 'optLazy':
            return self._strat.optLazy
        if name is 'optLL':
            return self._strat.optLL
        if name is 'optDelayNonMinimals':
            return self._strat.optDelayNonMinimals
        if name is 'optBrutalReductions':
            return self._strat.optBrutalReductions
        if name is 'optExchange':
            return self._strat.optExchange
        if name is 'optAllowRecursion':
            return self._strat.optAllowRecursion
        if name is 'optRedTailDegGrowth':
            return self._strat.optRedTailDegGrowth
        if name is 'optStepBounded':
            return self._strat.optStepBounded
        if name is 'optLinearAlgebraInLastBlock':
            return self._strat.optLinearAlgebraInLastBlock
        if name is 'optRedTailInLastBlock':
            return self._strat.optRedTailInLastBlock,
        if name is 'redByReduced':
            return self._strat.reduceByTailReduced
        if name is 'monomials':
            return new_BS_from_PBSet(self._strat.monomials)
        if name is 'minimalLeadingTerms':
            return new_BS_from_PBSet(self._strat.minimalLeadingTerms)
        if name is 'leadingTerms':
            return new_BS_from_PBSet(self._strat.leadingTerms)
        if name is 'llReductor':
            return new_BS_from_PBSet(self._strat.llReductor)
        else:
            raise AttributeError, name

    def __setattr__(self, name, val):
        if name is 'enabledLog':
            self._strat.enabledLog = val
        elif name is 'reductionSteps':
            self._strat.reductionSteps = val
        elif name is 'normalForms':
            self._strat.normalForms = val
        elif name is 'currentDegree':
            self._strat.currentDegree = val
        elif name is 'averageLength':
            self._strat.averageLength = val
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

cdef class BooleVariable:
    def index(self):
        return self._pbvar.index()

    def is_equal(self, BooleVariable other):
        return self._pbvar.is_equal(other._pbvar)

cdef inline BooleVariable new_BV_from_PBVar(PBVar juice):
    """
    Construct a new BooleVariable
    """
    cdef BooleVariable n
    n = <BooleVariable>PY_NEW(BooleVariable)
    n._pbvar = juice
    return n

cdef inline BooleVariable new_BV_from_int(int juice):
    """
    Construct a new BooleVariable
    """
    cdef BooleVariable n
    n = <BooleVariable>PY_NEW(BooleVariable)
    PBVar_construct_int(&n._pbvar, juice)
    return n

cdef class VariableBlock_base:
    def __init__(self, size, start_index, offset):
        self.size = size
        self.start_index = start_index
        self.offset = offset

cdef class VariableBlockTrue(VariableBlock_base):
    def __init__(self, size, start_index, offset):
        VariableBlock_base.__init__(self, size, start_index, offset)

    def __call__(self, int i):
        #FIXME: no index checking
        cdef PBVar v
        PBVar_construct_int(&v, self.offset+self.start_index+self.size-1-i)
        return new_BM_from_PBVar(get_cring()._monom_monoid, v)

cdef class VariableBlockFalse(VariableBlock_base):
    def __init__(self, size, start_index, offset):
        VariableBlock_base.__init__(self, size, start_index, offset)

    def __call__(self, int i):
        #FIXME: no index checking
        cdef PBVar v
        PBVar_construct_int(&v, i-self.start_index+self.offset)
        return new_BM_from_PBVar(get_cring()._monom_monoid, v)

def VariableBlock(size, start_index, offset, reverse):
    if reverse:
        return VariableBlockTrue(size, start_index, offset)
    else:
        return VariableBlockFalse(size, start_index, offset)

cdef int M4RI_init = 0

def init_M4RI():
    global M4RI_init
    if M4RI_init is int(0):
        buildAllCodes()
        setupPackingMasks()
        M4RI_init = 1

def free_m4ri():
    destroyAllCodes()

def recursively_insert(CCuddNavigator n, int ind, CCuddNavigator m):
    cdef PBSet b
    b = pb_recursively_insert((<CCuddNavigator>n)._pbnav, ind,
                                                (<CCuddNavigator>m)._pbnav)
    return new_BS_from_PBSet(b)

def ll_red_nf(BooleanPolynomial p, BooleSet reductors):
    cdef PBPoly t
    t = pb_ll_red_nf(p._pbpoly, reductors._pbset)
    return new_BP_from_PBPoly(p._parent, t)

def ll_red_nf_noredsb(BooleanPolynomial p, BooleSet reductors):
    cdef PBPoly t
    t = pb_ll_red_nf_noredsb(p._pbpoly, reductors._pbset)
    return new_BP_from_PBPoly(p._parent, t)

def mod_mon_set(BooleSet as, BooleSet vs):
    cdef PBSet b
    b = pb_mod_mon_set((<BooleSet>as)._pbset, (<BooleSet>vs)._pbset)
    return new_BS_from_PBSet(b)

def get_order_code():
    R = get_cring()
    return (<BooleanPolynomialRing>R)._pbring.getOrderCode()

def change_ordering(order):
    global cur_ring
    pb_change_ordering(order)
    cur_ring._pbring = get_current_ring()

def parallel_reduce(BooleanPolynomialVector inp, GroebnerStrategy strat, \
                                    int average_steps, double delay_f):
    return new_BPV_from_PBPolyVector(inp._parent, \
        pb_parallel_reduce(inp._vec, strat._strat, average_steps, delay_f))

def have_degree_order():
    global cur_ring
    return cur_ring._pbring.isDegreeOrder()

def set_variable_name( i, s):
    pb_set_variable_name(i, s)

def append_ring_block(i):
    pb_append_ring_block(i)

cdef BooleanPolynomialRing cur_ring
ring_callbacks = []

def get_cring():
    global cur_ring
    return cur_ring

def set_cring(BooleanPolynomialRing R):
    global cur_ring, ring_callbacks
    cur_ring = R
    for f in ring_callbacks:
        f()

def set_ring_callback(func):
    global ring_callbacks
    ring_callbacks.append(func)

def unpickle_BooleanPolynomial(ring, string):
    """
    Unpickle BooleanPolynomial

    EXAMPLE:
        sage: T = TermOrder('deglex',2)+TermOrder('deglex',2)
        sage: P.<a,b,c,d> = BooleanPolynomialRing(4,order=T)
        sage: loads(dumps(a+b)) == a+b
        True
    """
    return ring(eval(string,ring.gens_dict()))

def unpickle_BooleanPolynomialRing(n, names, order):
    """
    Unpickle BooleanPolynomialRing

    EXAMPLE:
        sage: T = TermOrder('deglex',2)+TermOrder('deglex',2)
        sage: P.<a,b,c,d> = BooleanPolynomialRing(4,order=T)
        sage: loads(dumps(P)) == P
        True
    """
    return BooleanPolynomialRing(n, names=names, order=order)

init_M4RI()

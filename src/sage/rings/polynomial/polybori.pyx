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

from sage.rings.integer import Integer
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.finite_field import GF
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
            n = Integer(n)
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
                    raise ValueError, "Each block should have the same order type (deglex or degrevlex) for block orderings."

        PBRing_construct(&self._R, n, pb_order_code)

        MPolynomialRing_generic.__init__(self, GF(2), n, names, order)

        counter = 0
        for i in range(len(order.blocks)-1):
            counter += order.blocks[i][1]
            pb_append_ring_block(counter)

        for i in range(self.ngens()):
            _n = self._names[i]
            self._R.setRingVariableName(i,_n)

        self._zero_element = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._zero_element)._P, 0)
        self._one_element  = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._one_element)._P, 1)

        self._monom_monoid = BooleanMonomialMonoid(self)

        set_cring(self)

    def __dealloc__(self):
        PBRing_destruct(&self._R)

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
        return self._R.nVariables()

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
        return new_BP_from_DD(self, self._R.variable(n))

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
        return tuple([new_BP_from_DD(self, self._R.variable(i)) \
                for i in xrange(self.ngens())])

    def _repr_(self):
        """
        EXAMPLE:
            sage: P.<x, y> = BooleanPolynomialRing(2)
            sage: P
            Boolean PolynomialRing in x, y
        """
        self._R.activate()
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
        """
        cdef BooleanPolynomial p
        if PY_TYPE_CHECK(other, BooleanMonomial) and \
                (<BooleanMonomial>other)._parent._ring is self:
            p = new_BP_from_PBMonom(self, (<BooleanMonomial>other)._M)
            return p
        else:
            raise TypeError, "coercion from %s to %s not implemented" % \
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
        """
        cdef int i
        cdef BooleanPolynomial p
        try:
            return self._coerce_c(other)
        except TypeError:
            pass

        if PY_TYPE_CHECK(other, DD):
            p = new_BP_from_DD(self, (<DD>other)._D)
        elif PY_TYPE_CHECK(other, BooleSet):
            p = new_BP_from_PBSet(self, (<BooleSet>other)._S)
        else:
            i = int(other)
            i = i % 2
            p = new_BP_from_int(self,i)
        return p

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
            sage: M = BooleanMonomialMonoid(P)
            sage: M
            MonomialMonoid of Boolean PolynomialRing in x, y

            sage: M.gens()
            (x, y)
            sage: type(M.gen(0))
            <type 'sage.rings.polynomial.polybori.BooleanMonomial'>
        """
        self._ring = polring

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
            sage: M = BooleanMonomialMonoid(P)
            sage: M.ngens()
            100
        """
        return self._ring.ngens()

    def gen(self, int n=0):
        """
        Return the n-th generator of self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: M.gen(0)
            x
            sage: M.gen(2)
            z

            sage: P = BooleanPolynomialRing(1000, 'x')
            sage: M = BooleanMonomialMonoid(P)
            sage: M.gen(50)
            x50
        """
        return new_BM_from_DD(self,
                (<BooleanPolynomialRing>self._ring)._R.variable(n))

    def gens(self):
        """
        Return the tuple of generators of self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: M.gens()
            (x, y, z)
        """
        return tuple([new_BM_from_DD(self,
            (<BooleanPolynomialRing>self._ring)._R.variable(i)) \
                for i in xrange(self.ngens())])

    def __call__(self, other = None):
        """
        Convert elements of other objects to elements of self.

        INPUT:
            other -- element to convert,
                     if None a BooleanMonomial representing 1 is returned
                     only BooleanPolynomials with the same parent ring as self which have a single monomial is converted

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: M(x+y)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce to BooleanMonomialMonoid

            sage: M(x*y)
            x*y
        """
        cdef BooleanMonomial m
        if other is None:
            m = new_BM(self)
            PBMonom_construct(&m._M)
            return m
        if PY_TYPE_CHECK(other, BooleanPolynomial) and \
                        (<BooleanPolynomial>other)._parent is self._ring and \
                        (<BooleanPolynomial>other)._P.isSingleton():
            m = new_BM_from_PBMonom(self,
                    (<BooleanPolynomial>other)._P.lead())
            return m
        else:
            raise TypeError, "cannot coerce to BooleanMonomialMonoid"

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
        PBMonom_construct(&self._M)
        _parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBMonom_destruct(&self._M)

    def __richcmp__(left, right, int op):
        # boilerplate code from sage.structure.parent
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Compare BooleanMonomial objects.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
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
        res = left._M.compare((<BooleanMonomial>right)._M)
        return res

    def _repr_(self):
        """
        Return a string representing self.

        EXAMPLES:
            sage: P.<x,y> = BooleanPolynomialRing(2)
            sage: M = BooleanMonomialMonoid(P)
            sage: M(x*y)
            x*y

            sage: R.<t,u> = BooleanPolynomialRing(2)
            sage: M(x*y)
            x*y
        """
        (<BooleanPolynomialRing>self._parent._ring)._R.activate()
        return PBMonom_to_str(&self._M)

    def __hash__(self):
        """
        Return a hash of self.
        """
        return self._M.hash()

    def deg(BooleanMonomial self):
        """
        Return total degree of self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: M(x*y).deg()
            2

            sage: M(x*x*y*z).deg()
            3
        """
        return self._M.deg()

    def __len__(BooleanMonomial self):
        """
        Return number of variables in self. This is equivalent to the total degree for BooleanMonomials.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: len(M(x*y))
            2
        """
        return self._M.deg()

    def __iter__(self):
        """
        Return an iterator over the indicies of the variables in self.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: list(iter(M(x*z)))
            [0, 2]
        """
        return new_BMI_from_PBMonomIter(self._M, self._M.begin())

    cdef MonoidElement _mul_c_impl(left, MonoidElement right):
        """
        Multiply self with another BooleanMonomial.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: x = M(x); xy = M(x*y); z=M(z)
            sage: x*x
            x

            sage: xy*y
            x*y

            sage: xy*z
            x*y*z
        """
        cdef BooleanMonomial m = new_BM_from_PBMonom(\
                (<BooleanMonomial>left)._parent, (<BooleanMonomial>left)._M)
        m._M.imul( (<BooleanMonomial>right)._M )
        return m

    def __add__(left, right):
        """
        Addition operator. Returns a BooleanPolynomial.

        EXAMPLES:
            sage: P.<x,y,z> = BooleanPolynomialRing(3)
            sage: M = BooleanMonomialMonoid(P)
            sage: x = M(x); xy = M(x*y)
            sage: x + xy
            x*y + x

            sage: x+0
            x

            sage: x+1
            x + 1
        """
        cdef BooleanPolynomial res
        if PY_TYPE_CHECK(left, BooleanMonomial):
            monom = left
            other = right
        elif PY_TYPE_CHECK(right, BooleanMonomial):
            monom = right
            other = left
        else:
            raise TypeError, "BooleanMonomial.__add__ called with not supported types %s and %s" %(type(right),type(left))

        res = new_BP_from_PBMonom( (<BooleanMonomial>monom)._parent._ring, \
                                        (<BooleanMonomial>monom)._M)
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
    PBMonom_construct_pbmonom(&m._M,juice)
    return m

cdef inline BooleanMonomial new_BM_from_PBVar(parent, PBVar juice):
    cdef BooleanMonomial m = new_BM(parent)
    PBMonom_construct_pbvar(&m._M,juice)
    return m

cdef inline BooleanMonomial new_BM_from_DD(parent, PBDD juice):
    cdef BooleanMonomial m = new_BM(parent)
    PBMonom_construct_dd(&m._M,juice)
    return m

cdef class BooleanMonomialIterator:
    def __iter__(self):
        return self

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
        PBPoly_construct(&self._P)
        self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBPoly_destruct(&self._P)

    def _repr_(self):
        """
        Return a string representing self.
        """
        (<BooleanPolynomialRing>self._parent)._R.activate()
        return PBPoly_to_str(&self._P)

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        cdef BooleanPolynomial p = new_BP_from_PBPoly(\
                (<BooleanPolynomial>left)._parent, (<BooleanPolynomial>left)._P)
        p._P.iadd( (<BooleanPolynomial>right)._P )
        return p

    cdef ModuleElement _sub_c_impl(left, ModuleElement right):
        return left._add_c_impl(right)

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        if left:
            return new_BP_from_PBPoly(left._parent, self._P)
        else:
            return 0

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        return self._rmul_c_impl(right)

    cdef RingElement _mul_c_impl(left, RingElement right):
        cdef BooleanPolynomial p = new_BP_from_PBPoly(\
                (<BooleanPolynomial>left)._parent, (<BooleanPolynomial>left)._P)
        p._P.imul( (<BooleanPolynomial>right)._P )
        return p

    def is_equal(self, BooleanPolynomial right):
        return self._P.is_equal(right._P)

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
        """
        return new_BPI_from_PBPolyIter(self, self._P.orderedBegin())

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
        elif self._P.isOne():
            return self
        elif self._P.isZero():
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
        return self._P.deg()

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
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._P.lead())

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
                                                self._P.lexLead())

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
        return self._P.isZero()

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
        return self._P.isOne()

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
        return self._P.isOne()

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
        return self._P.isConstant()

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
        return self._P.lmDeg()

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
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._P.usedVariables())

    def elimination_length(self):
        return self._P.eliminationLength()

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
        return self._P.length()

    def __getattr__(self, name):
        # this provides compatibility with the boost-python wrappers
        # of PolyBoRi and prevents us from polluting the sage namespace
        # i.e., these don't show up in tab completion lists
        if name == 'set':
            return new_BS_from_PBSet(self._P.set())
        elif name == 'deg':
            return self.total_degree
        elif name == 'elength':
            return self.elimination_length
        elif name == 'lead':
            return self.lm
        elif name == 'lexLead':
            return self.lm_lex
        elif name == 'constant':
            return self.is_constant
        elif name == 'lmDeg':
            return self.lm_degree
        elif name == 'isZero':
            return self.is_zero
        elif name == 'isOne':
            return self.is_one
        elif name == 'navigation':
            return new_CN_from_PBNavigator(self._P.navigation())
        else:
            raise AttributeError, name

cdef class BooleanPolynomialIterator:
    def __iter__(self):
        return self

    def __next__(self):
        cdef PBMonom val
        if self._iter.equal(self._obj._P.orderedEnd()):
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
    def __init__(self, ring, gens=[], coerce=True):
        MPolynomialIdeal.__init__(self, ring, gens, coerce)

    def groebner_basis(self):
        return groebner_basis_c_impl(self.ring(), self.gens())

cdef groebner_basis_c_impl(BooleanPolynomialRing R, g):
    cdef int i
    cdef PBPoly t
    cdef BooleanPolynomial p, r
    cdef PBPolyVector vec
    cdef GBStrategy strat

    GBStrategy_construct(&strat)
    for p in g:
        strat.addGeneratorDelayed(p._P)
    strat.symmGB_F2()
    vec = strat.minimalize()
    lvec = vec.size()
    res = []
    for i from 0 <= i < lvec:
        r = new_BP_from_PBPoly(R, vec.get(i) )
        res.append(r)
    return res

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
    PBPoly_construct_dd(&p._P,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBPoly(BooleanPolynomialRing parent,
        PBPoly juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbpoly(&p._P,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBMonom(BooleanPolynomialRing parent,
        PBMonom juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbmonom(&p._P,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBSet(BooleanPolynomialRing parent,
        PBSet juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbset(&p._P,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_int(BooleanPolynomialRing parent,
        int juice):
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_int(&p._P,juice)
    return p

cdef class DD:
    def __call__(self):
        return self

    def __dealloc__(self):
        PBDD_destruct(&self._D)

    def empty(self):
        return self._D.emptiness()

    def navigation(self):
        return new_CN_from_PBNavigator(self._D.navigation())

    def subset0(self, idx):
        return new_DD_from_PBDD(self._D.subset0(idx))

    def subset1(self, idx):
        return new_DD_from_PBDD(self._D.subset1(idx))

    def union(self, rhs):
        return new_DD_from_PBDD(self._D.unite((<DD>rhs)._D))

cdef inline DD new_DD_from_PBDD(PBDD juice):
    """
    Construct a new DD
    """
    cdef DD d
    d = <DD>PY_NEW(DD)
    d._D = juice
    return d

cdef class BooleSet:
    def __init__(self, param = None):
        if PY_TYPE_CHECK(param, DD):
            PBSet_construct_dd(&self._S, (<DD>param)._D)
        elif PY_TYPE_CHECK(param, CCuddNavigator):
            PBSet_construct_pbnav(&self._S, (<CCuddNavigator>param)._N)
        elif PY_TYPE_CHECK(param, BooleSet):
            PBSet_construct_pbset(&self._S, (<BooleSet>param)._S)
        else:
            PBSet_construct(&self._S)

    def __dealloc__(self):
        PBSet_destruct(&self._S)

    def __call__(self):
        return self

    def _repr_(self):
        return PBSet_to_str(&self._S)

    def empty(self):
        return self._S.emptiness()

    def navigation(self):
        return new_CN_from_PBNavigator(self._S.navigation())

    def cartesianProduct(self, rhs):
        return new_BS_from_PBSet(self._S.cartesianProduct((<BooleSet>rhs)._S))

    def diff(self, rhs):
        return new_BS_from_PBSet(self._S.diff((<BooleSet>rhs)._S))

    def change(self, ind):
        return new_BS_from_PBSet(self._S.change(ind))

    def usedVariables(self):
        return new_BM_from_PBMonom(get_cring()._monom_monoid,
                                            self._S.usedVariables())

    def if_then_else(self, int ind, BooleSet a, BooleSet b):
        cdef PBSet res
        if ind >= a.navigation().value() or ind >= b.navigation().value():
            raise IndexError, "value of ind must be less than the values of roots of the branches."
        PBSet_construct_indsetset(&res, ind, a._S.navigation(),
                b._S.navigation());
        return new_BS_from_PBSet(res)


    def __iter__(self):
        return new_BSI_from_PBSetIter(self, get_cring())


cdef inline BooleSet new_BS_from_PBSet(PBSet juice):
    """
    Construct a new BooleSet
    """
    cdef BooleSet s
    s = <BooleSet>PY_NEW(BooleSet)
    s._S = juice
    return s

cdef inline BooleSet new_BS_from_PBDD(PBDD juice):
    """
    Construct a new BooleSet
    """
    cdef BooleSet s
    s = <BooleSet>PY_NEW(BooleSet)
    PBSet_construct_dd(&s._S, juice)
    return s

cdef class BooleSetIterator:
    def __iter__(self):
        return self

    def __next__(self):
        cdef PBMonom val
        if self._iter.equal(self._obj._S.end()):
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
    m._iter = parent._S.begin()
    return m

cdef class CCuddNavigator:
    def __call__(self):
        return self

    def __dealloc__(self):
        PBNavigator_destruct(&self._N)

    def value(self):
        return self._N.value()

    def elseBranch(self):
        return new_CN_from_PBNavigator(self._N.elseBranch())

    def thenBranch(self):
        return new_CN_from_PBNavigator(self._N.thenBranch())

    def constant(self):
        return self._N.isConstant()

    def terminalOne(self):
        return self._N.isTerminated()

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
        self._vec.push_back(poly._P)

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
            GBStrategy_construct_gbstrategy(&self._S,
                    (<GroebnerStrategy>param)._S)
        else:
            GBStrategy_construct(&self._S)

    def __dealloc__(self):
        GBStrategy_destruct(&self._S)

    def addGeneratorDelayed(self, BooleanPolynomial p):
        if p.isZero():
            raise ValueError, "zero generators not allowed."
        self._S.addGeneratorDelayed(p._P)

    def addGenerator(self, BooleanPolynomial p, bint is_impl=False):
        if p.isZero():
            raise ValueError, "zero generators not allowed."
        if self._S.leadingTerms.owns(p._P.lead()):
            raise ValueError, "strategy already contains a polynomial with same lead"
        return self._S.addGenerator(p._P, is_impl)

    def addAsYouWish(self, BooleanPolynomial p):
        if p.isZero():
            raise ValueError, "zero generators not allowed."
        self._S.addAsYouWish(p._P)

    def implications(self, ind):
        implications(self._S, ind)

    def cleanTopByChainCriterion(self):
        self._S.cleanTopByChainCriterion()

    def symmGB_F2(self):
        self._S.symmGB_F2()

    def containsOne(self):
        return self._S.containsOne()

    def faugereStepDense(self, BooleanPolynomialVector v):
        return new_BPV_from_PBPolyVector(self._parent,
                                    self._S.faugereStepDense(v._vec))

    def minimalize(self):
        return new_BPV_from_PBPolyVector(self._parent, self._S.minimalize())

    def minimalizeAndTailReduce(self):
        return new_BPV_from_PBPolyVector(self._parent,
                                    self._S.minimalizeAndTailReduce())

    def npairs(self):
        return self._S.npairs()

    def topSugar(self):
        return pairs_top_sugar(self._S)

    def someSpolysInNextDegree(self, n):
        cdef PBPolyVector v = someNextDegreeSpolys(self._S, n)
        return new_BPV_from_PBPolyVector(self._parent, v)

    def __len__(self):
        return self._S.nGenerators()

    def __getitem__(self, int i):
        #FIXME: no index checking
        return new_BP_from_PBPoly(self._parent, GB_get_ith_gen(self._S, i))

    def __getattr__(self, name):
        if name is 'enabledLog':
            return self._S.enabledLog
        if name is 'reductionSteps':
            return self._S.reductionSteps
        if name is 'normalForms':
            return self._S.normalForms
        if name is 'currentDegree':
            return self._S.currentDegree
        if name is 'chainCriterions':
            return self._S.chainCriterions
        if name is 'variableChainCriterions':
            return self._S.variableChainCriterions
        if name is 'easyProductCriterions':
            return self._S.easyProductCriterions
        if name is 'extendedProductCriterions':
            return self._S.extendedProductCriterions
        if name is 'averageLength':
            return self._S.averageLength
        if name is 'optRedTail':
            return self._S.optRedTail
        if name is 'optLazy':
            return self._S.optLazy
        if name is 'optLL':
            return self._S.optLL
        if name is 'optDelayNonMinimals':
            return self._S.optDelayNonMinimals
        if name is 'optBrutalReductions':
            return self._S.optBrutalReductions
        if name is 'optExchange':
            return self._S.optExchange
        if name is 'optAllowRecursion':
            return self._S.optAllowRecursion
        if name is 'optRedTailDegGrowth':
            return self._S.optRedTailDegGrowth
        if name is 'optStepBounded':
            return self._S.optStepBounded
        if name is 'optLinearAlgebraInLastBlock':
            return self._S.optLinearAlgebraInLastBlock
        if name is 'optRedTailInLastBlock':
            return self._S.optRedTailInLastBlock,
        if name is 'redByReduced':
            return self._S.reduceByTailReduced
        if name is 'monomials':
            return new_BS_from_PBSet(self._S.monomials)
        if name is 'minimalLeadingTerms':
            return new_BS_from_PBSet(self._S.minimalLeadingTerms)
        if name is 'leadingTerms':
            return new_BS_from_PBSet(self._S.leadingTerms)
        if name is 'llReductor':
            return new_BS_from_PBSet(self._S.llReductor)
        else:
            raise AttributeError, name

    def __setattr__(self, name, val):
        if name is 'enabledLog':
            self._S.enabledLog = val
        elif name is 'reductionSteps':
            self._S.reductionSteps = val
        elif name is 'normalForms':
            self._S.normalForms = val
        elif name is 'currentDegree':
            self._S.currentDegree = val
        elif name is 'averageLength':
            self._S.averageLength = val
        elif name == 'optRedTail':
            self._S.optRedTail = val
        elif name is 'optLazy':
            self._S.optLazy = val
        elif name is 'optLL':
            self._S.optLL = val
        elif name is 'optDelayNonMinimals':
            self._S.optDelayNonMinimals = val
        elif name is 'optBrutalReductions':
            self._S.optBrutalReductions = val
        elif name is 'optExchange':
            self._S.optExchange = val
        elif name is 'optAllowRecursion':
            self._S.optAllowRecursion = val
        elif name is 'optRedTailDegGrowth':
            self._S.optRedTailDegGrowth = val
        elif name is 'optStepBounded':
            self._S.optStepBounded = val
        elif name is 'optLinearAlgebraInLastBlock':
            self._S.optLinearAlgebraInLastBlock = val
        elif name is 'optRedTailInLastBlock':
            self._S.optRedTailInLastBlock = val
        elif name is 'redByReduced':
            self._S.reduceByTailReduced = val
        else:
            raise AttributeError, name

cdef inline CCuddNavigator new_CN_from_PBNavigator(PBNavigator juice):
    """
    Construct a new CCuddNavigator
    """
    cdef CCuddNavigator n
    n = <CCuddNavigator>PY_NEW(CCuddNavigator)
    n._N = juice
    return n

cdef class BooleVariable:
    def index(self):
        return self._V.index()

    def is_equal(self, BooleVariable other):
        return self._V.is_equal(other._V)

cdef inline BooleVariable new_BV_from_PBVar(PBVar juice):
    """
    Construct a new BooleVariable
    """
    cdef BooleVariable n
    n = <BooleVariable>PY_NEW(BooleVariable)
    n._V = juice
    return n

cdef inline BooleVariable new_BV_from_int(int juice):
    """
    Construct a new BooleVariable
    """
    cdef BooleVariable n
    n = <BooleVariable>PY_NEW(BooleVariable)
    PBVar_construct_int(&n._V, juice)
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

def init_M4RI():
    buildAllCodes()
    setupPackingMasks()

def recursively_insert(CCuddNavigator n, int ind, CCuddNavigator m):
    cdef PBSet b
    b = pb_recursively_insert((<CCuddNavigator>n)._N, ind,
                                                (<CCuddNavigator>m)._N)
    return new_BS_from_PBSet(b)

def ll_red_nf(BooleanPolynomial p, BooleSet reductors):
    cdef PBPoly t
    t = pb_ll_red_nf(p._P, reductors._S)
    return new_BP_from_PBPoly(p._parent, t)

def ll_red_nf_noredsb(BooleanPolynomial p, BooleSet reductors):
    cdef PBPoly t
    t = pb_ll_red_nf_noredsb(p._P, reductors._S)
    return new_BP_from_PBPoly(p._parent, t)

def mod_mon_set(BooleSet as, BooleSet vs):
    cdef PBSet b
    b = pb_mod_mon_set((<BooleSet>as)._S, (<BooleSet>vs)._S)
    return new_BS_from_PBSet(b)

def get_order_code():
    R = get_cring()
    return (<BooleanPolynomialRing>R)._R.getOrderCode()

def change_ordering(order):
    global cur_ring
    pb_change_ordering(order)
    cur_ring._R = get_current_ring()

def parallel_reduce(BooleanPolynomialVector inp, GroebnerStrategy strat, \
                                    int average_steps, double delay_f):
    return new_BPV_from_PBPolyVector(inp._parent, \
        pb_parallel_reduce(inp._vec, strat._S, average_steps, delay_f))

def have_degree_order():
    global cur_ring
    return cur_ring._R.isDegreeOrder()

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

init_M4RI()

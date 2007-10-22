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

from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.finite_field import GF
from sage.monoids.monoid import Monoid_class

order_dict= {"lp":      lp,
             "dlex":    dlex,
             "dp_asc":  dp_asc,
             "bdlex":   block_dlex,
             "bdp_asc": block_dp_asc,
             }

cdef class BooleanPolynomialRing(MPolynomialRing_generic):
    """
    Boolean Polynomial Ring.
    """
    def __init__(self, n, names, order='lp'):
        cdef char *_n

        PBRing_construct(&self._R, n, order_dict.get(order, lp))
        MPolynomialRing_generic.__init__(self, GF(2), n, names, order)

        for i in range(self.ngens()):
            _n = self._names[i]
            self._R.setRingVariableName(i,_n)

        self._zero_element = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._zero_element)._P, 0)
        self._one_element  = new_BP(self)
        PBPoly_construct_int(&(<BooleanPolynomial>self._one_element)._P, 1)

        self._monom_monoid = BooleanMonomialMonoid(self)

    def __dealloc__(self):
        PBRing_destruct(&self._R)

    def ngens(self):
        return self._R.nVariables()

    def gen(self, int n=0):
        """
        """
        return new_BP_from_DD(self, self._R.variable(n))

    def gens(self):
        """
        """
        return [new_BP_from_DD(self, self._R.variable(i)) \
                for i in xrange(self.ngens())]

    def _repr_(self):
        self._R.activate()
        gens = ", ".join(map(str,self.gens()))
        return "Boolean PolynomialRing in %s"%(gens)

    cdef _coerce_c_impl(self, other):
        """
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

    def __hash__(self):
        """
        """
        return hash(str(self))

    def ideal(self, gens, coerce=True):
        if coerce:
            gens = [self(p) for p in gens]
        return BooleanPolynomialIdeal(self, gens)

class BooleanMonomialMonoid(Monoid_class):
    def __init__(self, polring):
        self._ring = polring

    def __repr__(self):
        return "MonomialMonoid of %s" % (str(self._ring))

    def __hash__(self):
        """
        """
        return hash(str(self))

    def ngens(self):
        return self._ring.ngens()

    def gen(self, int n=0):
        """
        """
        return new_BM_from_DD(self,
                (<BooleanPolynomialRing>self._ring)._R.variable(n))

    def gens(self):
        """
        """
        return [new_BM_from_DD(self,
            (<BooleanPolynomialRing>self._ring)._R.variable(i)) \
                for i in xrange(self.ngens())]

    def __call__(self, other):
        """
        """
        cdef BooleanMonomial m
        if PY_TYPE_CHECK(other, BooleanPolynomial) and \
                        (<BooleanPolynomial>other)._parent is self._ring and \
                        (<BooleanPolynomial>other)._P.isSingleton():
            m = new_BM_from_PBMonom(self,
                    (<BooleanPolynomial>other)._P.lead())
            return m
        else:
            raise TypeError, "cannot coerce to BooleanMonoid"

cdef class BooleanMonomial(MonoidElement):
    def __init__(self, parent):
        PBMonom_construct(&self._M)
        _parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBMonom_destruct(&self._M)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        cdef comparecodes res

        res = left._M.compare((<BooleanMonomial>right)._M)
        if op == 0:
            return res == less_than
        elif op == 1:
            return res == less_or_equal_max
        elif op == 2:
            return res == equality
        elif op == 3:
            return res != equality
        elif op == 4:
            return res == greater_than
        elif op == 5:
            return res == greater_or_equal_min

    def __repr__(self):
        (<BooleanPolynomialRing>self._parent._ring)._R.activate()
        return PBMonom_to_str(&self._M)

    def __hash__(self):
        return self._M.hash()

    def __iter__(self):
        return new_BMI_from_PBMonomIter(self._M, self._M.begin())

    def _mul_c_impl(left, MonoidElement right):
        cdef BooleanMonomial m = new_BM_from_PBMonom(\
                (<BooleanMonomial>left)._parent, (<BooleanMonomial>left)._M)
        m._M.imul( (<BooleanMonomial>right)._M )
        return m

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
    def __init__(self, parent):
        PBPoly_construct(&self._P)
        self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBPoly_destruct(&self._P)

    def __repr__(self):
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

    def __pow__(BooleanPolynomial self, int exp, ignored):
        """
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
            raise NotImplementedError, "Negative exponents for non constant polynomials are not implemented."

    def __neg__(BooleanPolynomial self):
        """
        """
        return self

    def total_degree(BooleanPolynomial self):
        """
        """
        return self._P.deg()

    def lm(BooleanPolynomial self):
        """
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._P.lead())

    def lt(BooleanPolynomial self):
        """
        """
        return self.lm()

    def is_zero(BooleanPolynomial self):
        """
        """
        return self._P.isZero()

    def is_one(BooleanPolynomial self):
        """
        """
        return self._P.isOne()

    def is_unit(BooleanPolynomial self):
        """
        """
        return self._P.isOne()

    def is_constant(BooleanPolynomial self):
        """
        """
        return self._P.isConstant()

    def lm_degree(BooleanPolynomial self):
        return self._P.lmDeg()

    def vars(self):
        return new_BM_from_PBMonom(self._parent._monom_monoid, self._P.usedVariables())

    def __getattr__(self, name):
        # this provides compatibility with the boost-python wrappers
        # of PolyBoRi and prevents us from polluting the sage objects namespace
        # i.e., these don't show up in tab completion lists
        if name == 'diagram':
            return new_DD_from_PBDD(self._P.diagram())
        elif name == 'lead':
            return self.lm
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

class BooleanPolynomialIdeal(MPolynomialIdeal):
    def __init__(self, ring, gens=[], coerce=True):
        MPolynomialIdeal.__init__(self, ring, gens, coerce)

    def groebner_basis(self):
        return groebner_basis_c_impl(self.ring(), self.gens())

cdef groebner_basis_c_impl(BooleanPolynomialRing R, g):
    cdef int i
    cdef PBPoly t
    cdef BooleanPolynomial p, r
    cdef PBPoly_vector vec
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
        else:
            PBSet_construct(&self._S)

    def __dealloc__(self):
        PBSet_destruct(&self._S)

    def empty(self):
        return self._S.emptiness()

    def navigation(self):
        return new_CN_from_PBNavigator(self._S.navigation())

cdef inline BooleSet new_BS_from_PBSet(PBSet juice):
    """
    Construct a new BooleSet
    """
    cdef BooleSet s
    s = <BooleSet>PY_NEW(BooleSet)
    s._S = juice
    return s

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

cdef inline CCuddNavigator new_CN_from_PBNavigator(PBNavigator juice):
    """
    Construct a new CCuddNavigator
    """
    cdef CCuddNavigator n
    n = <CCuddNavigator>PY_NEW(CCuddNavigator)
    n._N = juice
    return n

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

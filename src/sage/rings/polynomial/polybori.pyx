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
             "block_dlex":   block_dlex,
             "block_dp_asc": block_dp_asc,
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

        set_cring(self)

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

    def __call__(self, other = None):
        """
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
            raise TypeError, "cannot coerce to BooleanMonoid"

cdef class BooleanMonomial(MonoidElement):
    def __init__(self, BooleanMonomial parent):
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

    def deg(BooleanMonomial self):
        """
        """
        return self._M.deg()

    def __len__(BooleanMonomial self):
        return self._M.deg()

    def __iter__(self):
        return new_BMI_from_PBMonomIter(self._M, self._M.begin())

    cdef MonoidElement _mul_c_impl(left, MonoidElement right):
        cdef BooleanMonomial m = new_BM_from_PBMonom(\
                (<BooleanMonomial>left)._parent, (<BooleanMonomial>left)._M)
        m._M.imul( (<BooleanMonomial>right)._M )
        return m

    cdef BooleanPolynomial _add_c_impl(BooleanMonomial left, BooleanMonomial right):
        cdef BooleanPolynomial res
        res = new_BP_from_PBMonom( (<BooleanMonomial>left)._parent._ring, \
                                        (<BooleanMonomial>left)._M)
        res._P.iadd_PBMonom((<BooleanMonomial>right)._M)
        return res

    def __add__(self, right):
        if PY_TYPE_CHECK(self, BooleanMonomial) and \
            PY_TYPE_CHECK(right, BooleanMonomial) and \
            (<BooleanMonomial>self)._parent._ring is \
                (<BooleanMonomial>right)._parent._ring:
            return (<BooleanMonomial>self)._add_c_impl(right)

        if PY_TYPE_CHECK(self, BooleanMonomial):
            if right == 0:
                return self
            elif right == 1:
                return (<BooleanMonomial>self)._add_c_impl(\
                    (<BooleanMonomial>self)._parent())

        elif PY_TYPE_CHECK(right, BooleanMonomial):
            if self == 0:
                return right
            elif self == 1:
                return (<BooleanMonomial>right)._add_c_impl(\
                    (<BooleanMonomial>right)._parent())

        raise NotImplementedError, "add operation not implemented for types %s and %s"%(type(self),type(right))


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

    def is_equal(self, BooleanPolynomial right):
        #FIXME: change this to replace == operator
        return self._P.is_equal(right._P)

    def __iter__(self):
        return new_BPI_from_PBPolyIter(self, self._P.orderedBegin())

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

    def lm_lex(BooleanPolynomial self):
        """
        """
        return new_BM_from_PBMonom(self._parent._monom_monoid,
                                                self._P.lexLead())

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

    def elimination_length(self):
        return self._P.eliminationLength()

    def __len__(self):
        return self._P.length()

    def __getattr__(self, name):
        # this provides compatibility with the boost-python wrappers
        # of PolyBoRi and prevents us from polluting the sage objects namespace
        # i.e., these don't show up in tab completion lists
        if name == 'diagram':
            return new_DD_from_PBDD(self._P.diagram())
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

    def empty(self):
        return self._S.emptiness()

    def navigation(self):
        return new_CN_from_PBNavigator(self._S.navigation())

    def unateProduct(self, rhs):
        return new_BS_from_PBSet(self._S.unateProduct((<BooleSet>rhs)._S))

    def diff(self, rhs):
        return new_BS_from_PBSet(self._S.diff((<BooleSet>rhs)._S))

    def change(self, ind):
        return new_BS_from_PBSet(self._S.change(ind))

    def usedVariables(self):
        return new_BM_from_PBMonom(get_cring()._monom_monoid,
                                            self._S.usedVariables())

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
        self._S.addGeneratorDelayed(p._P)

    def addGenerator(self, BooleanPolynomial p, bint is_impl=False):
        return self._S.addGenerator(p._P, is_impl)

    def addAsYouWish(self, BooleanPolynomial p):
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
        elif name is 'chainCriterions':
            self._S.chainCriterions = val
        elif name is 'variableChainCriterions':
            self._S.variableChainCriterions = val
        elif name is 'easyProductCriterions':
            self._S.easyProductCriterions = val
        elif name is 'extendedProductCriterions':
            self._S.extendedProductCriterions = val
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
        cdef PBVar v
        PBVar_construct_int(&v, self.offset+self.start_index+self.size-1-i)
        return new_BM_from_PBVar(get_cring()._monom_monoid, v)

cdef class VariableBlockFalse(VariableBlock_base):
    def __init__(self, size, start_index, offset):
        VariableBlock_base.__init__(self, size, start_index, offset)

    def __call__(self, int i):
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

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
        cdef int i
        cdef BooleanPolynomial p
        i = int(other)
        i = i % 2
        p = new_BP(self)
        PBPoly_construct_int(&p._P,i)
        return p

    def __call__(self, other):
        return self._coerce_c(other)

    def __hash__(self):
        """
        """
        return hash(str(self))

    def ideal(self, gens, coerce=True):
        if coerce:
            gens = [self(p) for p in gens]
        return BooleanPolynomialIdeal(self, gens)

cdef class BooleanMonomial(MonoidElement):
    def __init__(self, parent):
        PBMonom_construct(&self._M)
        _parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBMonom_destruct(&self._M)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        #FIXME type check?
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
        (<BooleanPolynomialRing>self._parent)._R.activate()
        return str(PBMonom_to_str(&self._M))

    def __hash__(self):
        return self._M.hash()

cdef inline BooleanMonomial new_BM_from_PBMonom(BooleanPolynomialRing parent,
                                                            PBMonom juice):
    """
    Construct a new BooleanMonomial
    """
    cdef BooleanMonomial m
    m = <BooleanMonomial>PY_NEW(BooleanMonomial)
    m._M = juice
    m._parent = parent
    return m

cdef class BooleanPolynomial(MPolynomial):
    def __init__(self, val, parent = None):
        if PY_TYPE_CHECK(val, BooleanPolynomial):
            PBPoly_construct_pbpoly(&self._P, (<BooleanPolynomial>val)._P)
            if parent is None:
                self._parent = (<BooleanPolynomial>val)._parent
        else:
            PBPoly_construct(&self._P)
            self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBPoly_destruct(&self._P)

    def __repr__(self):
        (<BooleanPolynomialRing>self._parent)._R.activate()
        return str(PBPoly_to_str(&self._P))

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        cdef BooleanPolynomial p = new_BP((<BooleanPolynomial>left)._parent)
        p._P = PBPoly_add((<BooleanPolynomial>left)._P, (<BooleanPolynomial>right)._P)
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
        cdef BooleanPolynomial p = new_BP((<BooleanPolynomial>left)._parent)
        p._P = PBPoly_mul((<BooleanPolynomial>left)._P, (<BooleanPolynomial>right)._P)
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
        return new_BM_from_PBMonom(self._parent, self._P.lead())

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

    def __getattr__(self, name):
        if name == 'diagram':
            return new_DD_from_PBDD(self._P.diagram())
        elif name == 'lead':
            return self.lm
        elif name == 'constant':
            return self.is_constant
        elif name == 'lmDeg':
            return self.lm_degree
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
    cdef BooleanPolynomial p
    p = <BooleanPolynomial>PY_NEW(BooleanPolynomial)
    p._parent = parent
    return p

cdef inline BooleanPolynomial new_BP_from_DD(BooleanPolynomialRing parent, PBDD juice):
    """
    Construct a new BooleanPolynomial element
    """
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_dd(&p._P,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBPoly(BooleanPolynomialRing parent, PBPoly juice):
    """
    Construct a new BooleanPolynomial element
    """
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbpoly(&p._P,juice)
    return p

cdef inline BooleanPolynomial new_BP_from_PBMonom(BooleanPolynomialRing parent, PBMonom juice):
    """
    Construct a new BooleanPolynomial element
    """
    cdef BooleanPolynomial p = new_BP(parent)
    PBPoly_construct_pbmonom(&p._P,juice)
    return p

cdef class DD:
    def __call__(self):
        return self

    def __dealloc__(self):
        PBDD_destruct(&self._D)

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
        else:
            PBSet_construct(&self._S)

    def __dealloc__(self):
        PBSet_destruct(&self._S)

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

def rec_insert(CCuddNavigator n, int ind, CCuddNavigator m):
    cdef PBSet b
    b = recursively_insert((<CCuddNavigator>n)._N, ind, (<CCuddNavigator>m)._N)
    return new_BS_from_PBSet(b)

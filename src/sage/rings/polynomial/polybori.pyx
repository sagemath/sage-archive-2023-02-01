include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
include '../../libs/polybori/decl.pxi'

from sage.structure.element cimport RingElement
from sage.structure.element cimport ModuleElement

from sage.rings.finite_field import GF

order_dict= {"lp":      lp,
             "dlex":    dlex,
             "dp_asc":  dp_asc,
             "bdlex":   block_dlex,
             "bdp_asc": block_dp_asc,
             }

cdef class BPRing(MPolynomialRing_generic):
    def __init__(self, base_ring, n, names, order):
        cdef char *_n

        if base_ring is not GF(2):
            raise TypeError, "base ring must be GF(2)"

        PBRing_construct(&self._R, n, order_dict.get(order, lp))
        MPolynomialRing_generic.__init__(self, base_ring, n, names, order)

        for i in range(self.ngens()):
            _n = self._names[i]
            self._R.setRingVariableName(i,_n)

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
        gens = ", ".join(map(str,self.gens()))
        return "Boolean PolynomialRing in %s"%(gens)

cdef class BPolynomial(MPolynomial):
    def __init__(self, parent):
        self._P = PBPoly_new()
        self._parent = <ParentWithBase>parent

    def __dealloc__(self):
        PBPoly_delete(self._P)

    def __repr__(self):
        return str(PBPoly_to_str(self._P))

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        cdef PBPoly *l, *r, *t
        l = (<BPolynomial>left)._P
        r = (<BPolynomial>right)._P
        t = PBPoly_add(l, r)
        return new_BP_from_pPBPoly(left._parent, t)

    cdef ModuleElement _sub_c_impl(left, ModuleElement right):
        return left._add_c_impl(right)

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        if not left:
            return new_BP_from_pPBPoly(left._parent, self._P)
        else:
            return 0

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        return self._rmul_c_impl(right)

    cdef RingElement _mul_c_impl(left, RingElement right):
        cdef PBPoly *l, *r, *t
        l = (<BPolynomial>left)._P
        r = (<BPolynomial>right)._P
        t = PBPoly_mul(l, r)
        return new_BP_from_pPBPoly(left._parent, t)


cdef inline BPolynomial new_BP_from_DD(BPRing parent, PBDD juice):
    """
    Construct a new BPolynomial element
    """
    cdef BPolynomial p
    p = PY_NEW(BPolynomial)
    p._parent = <ParentWithBase>parent
    p._P = PBPoly_new_dd(juice)
    return p

cdef inline BPolynomial new_BP_from_pPBPoly(BPRing parent, PBPoly *juice):
    """
    Construct a new BPolynomial element
    """
    cdef BPolynomial p
    p = PY_NEW(BPolynomial)
    p._parent = <ParentWithBase>parent
    p._P = juice
    return p


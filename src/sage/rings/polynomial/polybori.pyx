include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
include '../../libs/polybori/decl.pxi'

from sage.structure.element cimport RingElement
from sage.structure.element cimport ModuleElement

order_dict= {"lp":      lp,
             "dlex":    dlex,
             "dp_asc":  dp_asc,
             "bdlex":   block_dlex,
             "bdp_asc": block_dp_asc,
             }

cdef class BPRing(MPolynomialRing_generic):
    def __init__(self, n, order):
        pass

    def __new__(self, n, order):
        PBRing_construct(&self._R, n, order_dict.get(order, lp))

    def __dealloc__(self):
        PBRing_destruct(&self._R)

    def ngens(self):
        return self._R.nVariables()

    def gens(self):
        return [new_BP_from_DD(self._R.variable(i)) \
                                    for i in xrange(self.ngens())]

    def __repr__(self):
        return "PolyBoRi Ring with %s generators"%self.ngens()

#cdef inline BPolynomial new_BP(BRing parent, PBDD juice):
cdef inline BPolynomial new_BP_from_DD(PBDD juice):
    """
    Construct a new BPolynomial element
    """
    cdef BPolynomial p
    p = PY_NEW(BPolynomial)
    #p._parent = <ParentWithBase>parent
    p._P = PBPoly_new_dd(juice)
    return p

cdef inline BPolynomial new_BP_from_pPBPoly(PBPoly *juice):
    """
    Construct a new BPolynomial element
    """
    cdef BPolynomial p
    p = PY_NEW(BPolynomial)
    #p._parent = <ParentWithBase>parent
    p._P = juice
    return p

cdef class BPolynomial(MPolynomial):
    #def __init__(self, BPolynomial parent):
    def __init__(self):
        #self._parent = parent
        pass

    #def __new__(self, BPolynomial parent):
    def __new__(self):
        self._P = PBPoly_new()

    def __dealloc__(self):
        PBPoly_delete(self._P)

    def __repr__(self):
        return str(PBPoly_to_str(self._P))

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        cdef PBPoly *l, *r, *t
        l = (<BPolynomial>left)._P
        r = (<BPolynomial>right)._P
        t = PBPoly_add(l, r)
        return new_BP_from_pPBPoly(t)

    cdef ModuleElement _sub_c_impl(left, ModuleElement right):
        return left._add_c_impl(right)

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        if not left:
            return new_BP_from_pPBPoly(self._P)
        else:
            return 0

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        return self._rmul_c_impl(right)

    cdef RingElement _mul_c_impl(left, RingElement right):
        cdef PBPoly *l, *r, *t
        l = (<BPolynomial>left)._P
        r = (<BPolynomial>right)._P
        t = PBPoly_mul(l, r)
        return new_BP_from_pPBPoly(t)

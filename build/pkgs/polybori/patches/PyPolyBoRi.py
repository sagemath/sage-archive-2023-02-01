from sage import all
from sage.rings.polynomial.pbori import *

class OrderCode:
    pass

OrderCode.__dict__ = order_dict

def Ring(n, order='lp'):
    return BooleanPolynomialRing(n, 'x', order=order)

BoolePolynomialVector = BooleanPolynomialVector

_gauss_on_polys=gauss_on_polys
def gauss_on_polys(l):
    vec=BoolePolynomialVector(l)
    return list(_gauss_on_polys(vec))

_add_up_polynomials=add_up_polynomials
def add_up_polynomials(polys, init):
    """
    Adds up the polynomials in polys (which should be a BoolePolynomialVector or a sequence of ???
    """
    if not isinstance(polys, BoolePolynomialVector):
        vec = BoolePolynomialVector
        for p in polys:
            vec.append(p)
        polys=vec

    return _add_up_polynomials(polys, init)

old_ring_var=Ring.var
def ring_var(self, i):
    warnings.warn('Ring.var is deprectated')
    return old_ring_var(self, i)

Ring.var=ring_var

def weakringref_call(self):
    if self.is_valid():
        return self.deref()
    return None

WeakRingRef.__call__ = weakringref_call

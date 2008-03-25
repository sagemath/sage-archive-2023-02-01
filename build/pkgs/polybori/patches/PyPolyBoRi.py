from sage import all
from sage.rings.polynomial.pbori import *

class RingProxy:
    def __init__(self, ring, order):
        self.ring = ring
        self.ordercode = order

    def set(self):
        self.ring._change_ordering(self.ordercode)

def global_ring():
    return RingProxy(get_cring(), get_order_code())

Polynomial = lambda x = None: get_cring()._coerce_(x)
Monomial = lambda x = None: get_cring()._monom_monoid(x)

class OrderCode:
    pass

OrderCode.__dict__ = order_dict

Variable = lambda x: get_cring().gen(x)

def Ring(n, order='lp'):
    return BooleanPolynomialRing(n, 'x', order=order)

BoolePolynomialVector = BooleanPolynomialVector

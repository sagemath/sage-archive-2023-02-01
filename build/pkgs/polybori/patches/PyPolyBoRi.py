from sage import all
from sage.rings.polynomial.pbori import *

def global_ring():
    return get_cring()

Polynomial = PolynomialFactory()

Monomial = MonomialFactory()

class OrderCode:
    pass

OrderCode.__dict__ = order_dict

Variable = VariableFactory()

def Ring(n, order='lp'):
    return BooleanPolynomialRing(n, 'x', order=order)

BoolePolynomialVector = BooleanPolynomialVector


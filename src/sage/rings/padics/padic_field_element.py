class pAdicFieldElement(pAdicFieldGenericElement):
    def __cmp__(self, other): raise NotImplementedError
    def __init__(self, parent, x, big_oh, ordp=None, construct=False): raise NotImplementedError
    def __invert__(self, prec=infinity): raise NotImplementedError
    def __mod__(self, right): raise NotImplementedError
    def __neg__(self): raise NotImplementedError
    def __pos__(self): raise NotImplementedError
    def __pow__(self, right): raise NotImplementedError
    def _add_(self, right): raise NotImplementedError
    def _div_(self, right): raise NotImplementedError
    def _floordiv_(self, right): raise NotImplementedError #return an error?
    def _integer_(self): raise NotImplementedError
    def _latex_(self): raise NotImplementedError
    def _mul_(self, right): raise NotImplementedError
    def _pari_init_(self): raise NotImplementedError
    def _repr_(self, do_latex=False): raise NotImplementedError
    def _sub_(self, right): raise NotImplementedError
    def additive_order(self): raise NotImplementedError
    def algdep(self, n): raise NotImplementedError
    def get_big_oh(self): raise NotImplementedError
    def set_big_oh(self): raise NotImplementedError
    def copy(self): raise NotImplementedError
    def element_sequence(self): raise NotImplementedError
    def exp(self): raise NotImplementedError
    def exp_artin_hasse(self): raise NotImplementedError
    def gamma(self): raise NotImplementedError
    def is_integral(self): raise NotImplementedError
    def is_square(self): raise NotImplementedError
    def is_unit(self): raise NotImplementedError #different from integer
    def is_weakly_zero(self, modulus): raise NotImplementedError
    def are_weakly_equal(self, right): raise NotImplementedError
    def log(self): raise NotImplementedError
    def log_artin_hasse(self): raise NotImplementedError
    def min_poly(self): raise NotImplementedError
    def minimal_polynomial(self): raise NotImplementedError
    def multiplicative_order(self): raise NotImplementedError
    def norm(self, ground=None): raise NotImplementedError
    def ordp(self): raise NotImplementedError
    def prec(self): raise NotImplementedError
    def get_precision(self): raise NotImplementedError
    def set_precision(self, n): raise NotImplementedError
    def rational_reconstruction(self): raise NotImplementedError
    def sqrt(self): raise NotImplementedError
    def square_root(self): raise NotImplementedError
    def trace(self, ground=None): raise NotImplementedError
    def unit_part(self): raise NotImplementedError
    def v(self): raise NotImplementedError
    def valuation(self): raise NotImplementedError

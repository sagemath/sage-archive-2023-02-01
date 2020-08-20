if __name__ == '__main__':
    from sys import path as search_path
    from os import path as file_path
    search_path.append(file_path.join(file_path.dirname(__file__), '..'))

from .PyPolyBoRi import Monomial, Polynomial, BooleSet, BooleConstant
from .PyPolyBoRi import Variable as VariableType


class IntegerPolynomial(object):
    """Polynomial with positive integer coefficients"""

    def __init__(self, boolean_polys):
        super(IntegerPolynomial, self).__init__()
        if not isinstance(boolean_polys, dict):
            boolean_polys = dict([(0, Polynomial(boolean_polys))])
        self.boolean_polys = boolean_polys

    def __coerce__(self, other):
        #TODO handle long
        if isinstance(other, int) and other >= 0:
            i = 0
            res = []
            while 2 ** i <= other:
                and_ = 2 ** i & other
                if and_:
                    res.append(i)
            return (self, IntegerPolynomial(dict([(i, BooleConstant(1)) for i
                in res])))
        if not isinstance(other, IntegerPolynomial):
            other = Polynomial(other)
        return (self, IntegerPolynomial(dict([(0, other)])))

    def __add__(self, other):
        r"""
        
        TESTS::
        
            sage: from sage.rings.polynomial.pbori import *
            sage: r= declare_ring([Block("x",1000)], globals()) # doctest: +ELLIPSIS
            sage: p=IntegerPolynomial(x(1))
            sage: p
            {0: x(1)}
            sage: p=p+p;p
            {1: x(1)}
            sage: p=p+x(2); p
            {0: x(2), 1: x(1)}
            sage: p+p
            {1: x(2), 2: x(1)}
	"""
        if not isinstance(other, IntegerPolynomial):
            (self, other) = coerce(self, other)
            return self + other
        assert isinstance(other, IntegerPolynomial)

        def add_simple_poly(p, i):
            p = Polynomial(p)
            if p.is_zero():
                return
            res_p = Polynomial(res.get(i, Polynomial(0, p.ring())))
            res[i] = res_p + p
            if res[i].is_zero():
                del res[i]
            inter = BooleSet(res_p).intersect(BooleSet(p))
            add_simple_poly(inter, i + 1)
            return
        res = dict(self.boolean_polys.items())
        for (k, p) in other.boolean_polys.items():
            add_simple_poly(p=p, i=k)
        return IntegerPolynomial(res)

    def __unicode__(self):
        return unicode(self.boolean_polys)

    def __str__(self):
        return str(self.boolean_polys)

    def __repr__(self):
        try:
            return unicode(self)
        except NameError:
            return str(self)


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

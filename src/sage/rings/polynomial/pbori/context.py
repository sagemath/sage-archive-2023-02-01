if __name__ == '__main__':
    from sys import path as search_path
    from os import path as file_path
    search_path.append(file_path.join(file_path.dirname(__file__), '..'))


def _exists():
    r"""
    PolyBoRi convention: checking optional components for prerequisites here
    
    TESTS::

        sage: from sage.rings.polynomial.pbori.context import *
        sage: _exists()
        True
    """
    from distutils.sysconfig import get_python_version
    return float(get_python_version()) > 2.4

from .PyPolyBoRi import Ring, VariableFactory, MonomialFactory
from .PyPolyBoRi import PolynomialFactory, SetFactory
from .PyPolyBoRi import Variable, Monomial, Polynomial, BooleSet
import polybori


class FactoryContext(object):
    r"""
    Temporarily exchange the constructor of a given type with a compatible
    callable object. It is useful together with the with statement.

    EXAMPLES::
    
        sage: r = Ring(1000)
        sage: from sage.rings.polynomial.pbori import Variable
        sage: def var(idx): return Variable(idx, r)
        sage: with FactoryContext(Variable, var):
        ...     print Variable(17)
        x(17)
        sage: try:
        ...     print Variable(17)
        ...   except:
        ...     print "caught expected exception"
        caught expected exception
    """

    def __init__(self, original, factory):
        self.original = original
        self.factory = factory

    def __enter__(self):
        self.fallback = self.original.__init__

        def func(orig, *args):
            try:
                self.fallback(orig, self.factory(*args))
            except:
                self.fallback(orig, *args)

        self.original.__init__ = func
        return self

    def __exit__(self, type, value, traceback):
        self.original.__init__ = self.fallback
        return False


class RingContext(object):
    r"""
    Temporarily fix the ring for constructors of some ring-dependent types
    like Variable and Monomial to a given ring. It is useful together with
    the with statement.

    EXAMPLES::
    
        sage: r = Ring(1000)
        sage: from sage.rings.polynomial.pbori import Variable
        sage: print Variable(17, r)
        x(17)
        sage: with RingContext(r):
        ...     print Variable(17), Monomial(), Polynomial(0), BooleSet()
        x(17) 1 0 {}
        sage: try:
        ...     print Variable(17)
        ...   except:
        ...     print "caught expected exception"
        caught expected exception
    """

    def __init__(self, ring):
        self.contexts = (FactoryContext(Variable, VariableFactory(ring)),
                         FactoryContext(Monomial, MonomialFactory(ring)),
                         FactoryContext(Polynomial, PolynomialFactory(ring)),
                         FactoryContext(BooleSet, SetFactory(ring)))

    def __enter__(self):
        for elt in self.contexts:
            elt.__enter__()
        return self

    def __exit__(self, type, value, traceback):
        result = False
        for elt in reversed(self.contexts):
            result = result or elt.__exit__(type, value, traceback)
        return result


if __name__ == '__main__':
    import doctest
    doctest.testmod()

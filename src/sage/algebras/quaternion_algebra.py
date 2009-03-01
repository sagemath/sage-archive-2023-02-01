"""

"""

from sage.rings.ring import Algebra

def QuaternionAlgebra(base_ring, a, b):
    return QuaternionAlgebra_optimized(base_ring, a, b)

class QuaternionAlgebra_abstract(Algebra):
    def __init__(self, base_ring, names='i,j,k'):
        self.__base_ring = base_ring
        self.__names = names

    def base_ring(self):
        return self.__base_ring

    def _repr_(self):
        return "Quaternion Algebra with base ring %s"%self.__base_ring

from quaternion_algebra_element import QuaternionAlgebraElement_generic

class QuaternionAlgebra_optimized(QuaternionAlgebra_abstract):
    def __init__(self, base_ring, a, b, names='i,j,k'):
        QuaternionAlgebra_abstract.__init__(self, base_ring, names=names)
        self._a = a
        self._b = b
        self._populate_coercion_lists_()
        self.__gens = [self([0,1,0,0]), self([0,0,1,0]), self([0,0,0,1])]

    def _repr_(self):
        return "Optimized Quaternion Algebra (%s, %s) with base ring %s"%(self._a, self._b, self.base_ring())

    def _element_constructor_(self, x):
        """
        Make sure x defines a valid member of self, then return the
        constructed element.
        """
        if isinstance(x, QuaternionAlgebraElement_generic):
            return x  # ok, since they are immutable
        return QuaternionAlgebraElement_generic(self, x[0], x[1], x[2], x[3])

    def gen(self, i=0):
        """

        """
        return self.__gens[i]


def QuaternionAlgebraWithInnerProduct(K, norms, traces, names='i,j,k'):
    raise NotImplementedError

def QuaternionAlgebraWithGramMatrix(K, gram, names='i,j,k'):
    raise NotImplementedError

def QuaternionAlgebraWithDiscriminants(D1, D2, T, names=['i','j','k'], M=2):
    raise NotImplementedError


########################################################
# Utility functions
########################################################
from sage.rings.integer import Integer
def ramified_primes(a,b):
    """
    Return sorted list of the finite primes that ramify in Q(a,b).

    EXAMPLES::

        sage: sage.algebras.quaternion_algebra.ramified_primes(-1,-1)
        [2]
        sage: sage.algebras.quaternion_algebra.ramified_primes(-1,-7)
        [7]
        sage: sage.algebras.quaternion_algebra.ramified_primes(-1,7)
        [2, 7]
    """
    a = Integer(a); b = Integer(b)
    if a.is_square() or b.is_square() or (a+b).is_square():
        return [ ]
    a = a.squarefree_part()
    b = b.squarefree_part()
    c = Integer(a.gcd(b))
    if c != 1:
        p = c.factor()[0][0]
        ram_prms = ramified_primes(p,-(b//p))
        for p in ramified_primes(a//p,b):
            if p in ram_prms:
                ram_prms.remove(p)
            else:
                ram_prms.append(p)
        ram_prms.sort()
        return ram_prms
    ram_prms = [ ]
    S1 = [ p[0] for p in abs(a).factor() ]
    for p in S1:
        if p == 2 and b%4 == 3:
            if (a+b).kronecker(p) == -1:
                ram_prms.append(p)
        elif b.kronecker(p) == -1:
            ram_prms.append(p)
    S2 = [ p[0] for p in abs(b).factor() ]
    for q in S2:
        if q == 2 and a%4 == 3:
            if (a+b).kronecker(q) == -1:
                ram_prms.append(q)
        elif a.kronecker(q) == -1:
            ram_prms.append(q)
    if not 2 in ram_prms and a%4 == 3 and b%4 == 3:
        ram_prms.append(2)
    ram_prms.sort()
    return ram_prms

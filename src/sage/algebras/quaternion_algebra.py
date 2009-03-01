"""
Quaternion Algebras

AUTHORS:
    * Jon Bobber
    * David Kohel
    * William Stein
"""

from sage.rings.ring import Algebra
from sage.structure.parent_gens import ParentWithGens

def QuaternionAlgebra(base_ring, a, b, names='i,j,k'):
    return QuaternionAlgebra_optimized(base_ring, a, b, names)

class QuaternionAlgebra_abstract(Algebra):
    def _repr_(self):
        return "Quaternion Algebra with base ring %s"%self.base_ring()



from quaternion_algebra_element import QuaternionAlgebraElement_generic

class QuaternionAlgebra_optimized(QuaternionAlgebra_abstract):
    def __init__(self, base_ring, a, b, names='i,j,k'):
        ParentWithGens.__init__(self, base_ring, names)
        self._a = a
        self._b = b
        self._populate_coercion_lists_(coerce_list=[base_ring])
        self.__gens = [self([0,1,0,0]), self([0,0,1,0]), self([0,0,0,1])]

    def ngens(self):
        return 3

    def _repr_(self):
        return "Optimized Quaternion Algebra (%s, %s) with base ring %s"%(self._a, self._b, self.base_ring())

    def _element_constructor_(self, x):
        """
        Make sure x defines a valid member of self, then return the
        constructed element.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-1,-2)
            sage: a = Q._element_constructor_(2/3); a
            2/3
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
            sage: Q._element_constructor_(a)
            2/3
            sage: Q._element_constructor_([1,2,3,4])
            1 + 2*i + 3*j + 4*k
            sage: Q._element_constructor_((1,2,3,4))
            1 + 2*i + 3*j + 4*k
            sage: Q._element_constructor_(-3/5)
            -3/5
        """
        if isinstance(x, QuaternionAlgebraElement_generic):
            return x  # ok, since they are immutable
        R = self.base_ring()
        if isinstance(x, (list, tuple)):
            return QuaternionAlgebraElement_generic(self, R(x[0]), R(x[1]), R(x[2]), R(x[3]))
        else:
            return QuaternionAlgebraElement_generic(self, R(x), R(0), R(0), R(0))

    def gen(self, i=0):
        """
        EXAMPLES::

            sage: Q.<ii,jj,kk> = QuaternionAlgebra(QQ,-1,-2); Q
            Optimized Quaternion Algebra (-1, -2) with base ring Rational Field
            sage: Q.gen(0)
            ii
            sage: Q.gen(1)
            jj
            sage: Q.gen(2)
            kk
            sage:
            sage: Q.gens()
            (ii, jj, kk)
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
from sage.rings.arith import GCD, fundamental_discriminant, hilbert_symbol
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


def ramified_primes_from_discs(D1, D2, T):
    M = Integer(GCD([D1,D2,T]))
    D3 = (T**2 - D1*D2)//4
    facs = D3.factor()
    D1 = fundamental_discriminant(D1)
    D2 = fundamental_discriminant(D2)
    D3 = fundamental_discriminant(D3)
    ram_prms = []
    for pow in facs:
        p = pow[0]
        if pow[1] % 2 == 1:
            chi = (D.kronecker(p) for D in (D1,D2,D3))
            if -1 in chi:
                ram_prms.append(p)
            elif not 1 in chi and hilbert_symbol(D1,D3,p) == -1:
                ram_prms.append(p)
        elif D1%p == 0 and D2%p == 0:
            chi = (D.kronecker(p) for D in (D1,D2,D3))
            if hilbert_symbol(D1,D3,p) == -1:
                ram_prms.append(p)
    return ram_prms

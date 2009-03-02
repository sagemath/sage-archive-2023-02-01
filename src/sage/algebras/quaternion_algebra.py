"""
Quaternion Algebras

AUTHORS:

- Jon Bobber -- 2009 rewrite

- William Stein -- 2009 rewrite

This code is partly based on Sage code by David Kohel from 2005.

TESTS::

We test pickles::

    sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
    sage: Q == loads(dumps(Q))
    True
"""

#*****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#       Copyright (C) 2009 Jonathon Bober <jwbober@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.arith import (GCD, fundamental_discriminant, hilbert_symbol,
                              hilbert_conductor_inverse, hilbert_conductor)
from sage.rings.integer import Integer
from sage.rings.rational import Rational

from sage.rings.ring import Algebra, is_Field
from sage.rings.rational_field import is_RationalField
from sage.rings.number_field.number_field import is_NumberField
from sage.structure.parent_gens import ParentWithGens
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import diagonal_matrix
from sage.structure.sequence import Sequence
from sage.structure.element import is_Element
from sage.modules.free_module import VectorSpace

########################################################
# Constructor
########################################################

_cache = {}

def QuaternionAlgebra(arg0, arg1=None, arg2=None, names='i,j,k'):
    """
    There are three input formats::

        - QuaternionAlgebra(a, b): quaternion algebra with i^2 = a, j^2 = b, i*j=k

        - QuaternionAlgebra(K, a, b): same as above but over a field K

        - QuaternionAlgebra(D): a rational quaternion algebra with discriminant D

    EXAMPLES::

    QuaternionAlgebra(a, b) - return quaternion algebra over the
    smallest field containing a and b with generators i, j, k with
    i^2=a, j^2=b and i*j=-j*i.::

        sage: QuaternionAlgebra(-2,-3)
        Quaternion Algebra (-2, -3) with base ring Rational Field
        sage: QuaternionAlgebra(GF(5)(2), GF(5)(3))
        Quaternion Algebra (2, 3) with base ring Finite Field of size 5
        sage: QuaternionAlgebra(QQ[sqrt(2)](-1), -5)
        Quaternion Algebra (-1, -5) with base ring Number Field in sqrt2 with defining polynomial x^2 - 2
        sage: QuaternionAlgebra(sqrt(-1), sqrt(-3))
        Quaternion Algebra (I, sqrt(3)*I) with base ring Symbolic Ring

    QuaternionAlgebra(K, a, b) - return quaternion algebra over the
    field K with generators i, j, k with i^2=a, j^2=b and i*j=-j*i.::

        sage: QuaternionAlgebra(QQ, -7, -21)
        Quaternion Algebra (-7, -21) with base ring Rational Field
        sage: QuaternionAlgebra(QQ[sqrt(2)], -2,-3)
        Quaternion Algebra (-2, -3) with base ring Number Field in sqrt2 with defining polynomial x^2 - 2

    QuaternionAlgebra(D) - D a squarefree integer; returns a rational
    quaternion algebra of discriminant D.::

        sage: QuaternionAlgebra(1)
        Quaternion Algebra (-1, 1) with base ring Rational Field
        sage: QuaternionAlgebra(2)
        Quaternion Algebra (-1, -1) with base ring Rational Field
        sage: QuaternionAlgebra(7)
        Quaternion Algebra (-1, -7) with base ring Rational Field
        sage: QuaternionAlgebra(2*3*5*7)
        Quaternion Algebra (-22, 210) with base ring Rational Field

    If the coefficients `a` and `b` in the definition of the quaternion
    algebra are not integral, then a slower generic type is used for
    arithmetic::

        sage: type(QuaternionAlgebra(-1,-3).0)
        <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
        sage: type(QuaternionAlgebra(-1,-3/2).0)
        <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_generic'>
    """
    global _cache
    key = (arg0, arg1, arg2, names)
    if _cache.has_key(key):
        return _cache[key]

    # QuaternionAlgebra(D)
    if arg1 is None and arg2 is None:
        A = QuaternionAlgebra_disc(arg0, names=names)

    elif arg2 is None:
        if is_Element(arg0):
            # QuaternionAlgebra(a, b)
            v = Sequence([arg0, arg1])
            K = v.universe().fraction_field()
            A = QuaternionAlgebra_ab(K, v[0], v[1], names=names)
        else:
            raise ValueError, "unknown input"

    # QuaternionAlgebra(K, a, b)
    else:
        A = QuaternionAlgebra_ab(arg0, arg1, arg2, names=names)

    A._key = key
    _cache[key] = A
    return A

def QuaternionAlgebra_disc(D, names):
    """
    Return a rational quaternion algebra of discriminant D.

    INPUT:

         D -- square-free positive integer
         names -- variable names

    EXAMPLES::

        sage: sage.algebras.quaternion_algebra.QuaternionAlgebra_disc(2, 'i,j,k')
        Quaternion Algebra (-1, -1) with base ring Rational Field
    """
    a, b = hilbert_conductor_inverse(D)
    a = Rational(a); b = Rational(b)
    return QuaternionAlgebra_ab(a.parent(), a, b,  names=names)



########################################################
# Classes
########################################################

def is_QuaternionAlgebra(A):
    """
    Return True if A is of the QuaternionAlgebra data type.

    EXAMPLES::

        sage: sage.algebras.quaternion_algebra.is_QuaternionAlgebra(QuaternionAlgebra(QQ,-1,-1))
        True
        sage: sage.algebras.quaternion_algebra.is_QuaternionAlgebra(ZZ)
        False
    """
    return isinstance(A, QuaternionAlgebra_abstract)

class QuaternionAlgebra_abstract(Algebra):
    def _repr_(self):
        """
        EXAMPLES::

            sage: sage.algebras.quaternion_algebra.QuaternionAlgebra_abstract(QQ)._repr_()
            'Quaternion Algebra with base ring Rational Field'
        """
        return "Quaternion Algebra with base ring %s"%self.base_ring()

    def ngens(self):
        """
        Return number of generators of this quaternion algebra.

        Though quaternion algebras could be generated using 2 elements,
        we always chose three algebra generators.  E.g., in many cases
        the generators are i, j, and k. [[explain better.]]

        EXAMPLES::
            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: Q.ngens()
            3
            sage: Q.gens()
            [i, j, k]
        """
        return 3

    def basis(self):
        """
        Return the fixed basis of self, which is 1, i, j, k, where
        i,j,k are the gens of self.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: Q.basis()
            (1, i, j, k)

            sage: Q.<xyz,abc,theta> = QuaternionAlgebra(GF(9,'a'),-5,-2)
            sage: Q.basis()
            (1, xyz, abc, theta)

        The basis is cached:
            sage: Q.basis() is Q.basis()
            True
        """
        try:
            return self.__basis
        except AttributeError:
            self.__basis = tuple([self(1)] + list(self.gens()))
            return self.__basis

    def inner_product_matrix(self):
        """
        Return the inner product matrix associated to self.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(-5,-19)
            sage: sage.algebras.quaternion_algebra.QuaternionAlgebra_abstract.inner_product_matrix(Q)
            [  2   0   0   0]
            [  0  10   0   0]
            [  0   0  38   0]
            [  0   0   0 190]

        """
        try: return self.__inner_product_matrix
        except AttributeError: pass

        K = self.base_ring()
        M = MatrixSpace(K,4)(0)
        B = self.basis()
        for i in range(4):
            x = B[i]
            M[i,i] = 2*(x.reduced_norm())
            for j in range(i+1,4):
                y = B[j]
                c = (x * y.conjugate()).reduced_trace()
                M[i,j] = c
                M[j,i] = c
        M.set_immutable()
        self.__inner_product_matrix = M
        return M

    def is_commutative(self):
        """
        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -3,-7)
            sage: Q.is_commutative()
            False
        """
        return False

    def is_division_algebra(self):
        """
        Return True if the quaternion algebra is a division algebra (i.e. a
        ring, not necessarily commutative, in which every nonzero element
        is invertible). So if this returns False, the quaternion algebra is
        isomorphic to the 2x2 matrix algebra.

        EXAMPLES::

            sage: QuaternionAlgebra(QQ,-5,-2).is_division_algebra()
            True
            sage: QuaternionAlgebra(1).is_division_algebra()
            False
        """
        return self.discriminant() != 1

    def is_exact(self):
        """
        Return True if elements of this quaternion algebra are represented
        exactly, i.e. there is no precision loss when doing arithmetic. A
        quaternion algebra is exact if and only if its base field is
        exact.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -3, -7)
            sage: Q.is_exact()
            True
            sage: Q.<i,j,k> = QuaternionAlgebra(Qp(7), -3, -7)
            sage: Q.is_exact()
            False
        """
        return self.base_ring().is_exact()

    def is_field(self):
        """
        Return False always, since all quaternion algebras are
        noncommutative and all fields are commutative.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -3, -7)
            sage: Q.is_field()
            False
        """
        return False

    def is_finite(self):
        """
        Return True if the quaternion algebra is finite as a set.

        Algorithm: A quaternion algebra is finite if and only if the
        base field is finite.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -3, -7)
            sage: Q.is_finite()
            False
            sage: Q.<i,j,k> = QuaternionAlgebra(GF(5), -3, -7)
            sage: Q.is_finite()
            True
        """
        return self.base_ring().is_finite()

    def is_integral_domain(self):
        """
        Return False always, since all quaternion algebras are
        noncommutative and integral domains are commutative (in Sage).

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -3, -7)
            sage: Q.is_integral_domain()
            False
        """
        return False

    def is_noetherian(self):
        """
        Return True always, since any quaternion algebra is a noetherian
        ring (because it's a finitely-generated module over a field, which
        is noetherian).

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -3, -7)
            sage: Q.is_noetherian()
            True
        """
        return True

    def order(self):
        """
        Return the number of elements of the quaternion algebra, or
        +Infinity if the algebra is not finite.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -3, -7)
            sage: Q.order()
            +Infinity
            sage: Q.<i,j,k> = QuaternionAlgebra(GF(5), -3, -7)
            sage: Q.order()
            20
        """
        return 4*self.base_ring().order()

    def random_element(self, *args, **kwds):
        """
        Return a random element of this quaternion algebra.

        The args and kwds are passed to the random_element method of
        the base ring.

        EXAMPLES::

            sage: QuaternionAlgebra(QQ[sqrt(2)],-3,7).random_element()
            -4 + (-1/95*sqrt2 - 1/2)*i + (-12*sqrt2 + 1/2)*j + (1/2*sqrt2 - 1)*k

            sage: QuaternionAlgebra(-3,19).random_element()
            -1/4 + 2/3*i - 5/2*j

        Specify the numerator and denominator bounds::

            sage: QuaternionAlgebra(-3,19).random_element(10^6,10^6)
            -979933/553629 + 255525/657688*i - 3511/6929*j - 700105/258683*k



        """
        K = self.base_ring()
        return self([ K.random_element(*args, **kwds) for _ in range(4) ])

    def vector_space(self):
        """
        Return vector space with inner product associated to self.

        EXAMPLES::

            sage: QuaternionAlgebra(-3,19).vector_space()
            Ambient quadratic space of dimension 4 over Rational Field
            Inner product matrix:
            [   2    0    0    0]
            [   0    6    0    0]
            [   0    0  -38    0]
            [   0    0    0 -114]
        """
        try:
            return self.__vector_space
        except AttributeError:
            V = VectorSpace(self.base_ring(), 4, inner_product_matrix = self.inner_product_matrix())
            self.__vector_space = V
            return V

import quaternion_algebra_element

class QuaternionAlgebra_ab(QuaternionAlgebra_abstract):
    """
    The quaternion algebra of the form (a,b/R), where i^2=a, j^2 = b,
    and i*j=-j*i=k.

    EXAMPLES::

    """
    def __init__(self, base_ring, a, b, names='i,j,k'):
        """
        Create the quaternion algebra with i^2 = a, j^2 = b, and i*j =
        -j*i = k.

        INPUT:
            - ``base_ring`` -
            - ``a, b`` -
            - ``names`` -

        TESTS::

        Test making quaternion elements (using the element constructor)::

            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-1,-2)
            sage: a = Q(2/3); a
            2/3
            sage: type(a)
            <type 'sage.algebras.quaternion_algebra_element.QuaternionAlgebraElement_rational_field'>
            sage: Q(a)
            2/3
            sage: Q([1,2,3,4])
            1 + 2*i + 3*j + 4*k
            sage: Q((1,2,3,4))
            1 + 2*i + 3*j + 4*k
            sage: Q(-3/5)
            -3/5

        The base ring must be a field:
            sage: Q.<ii,jj,kk> = QuaternionAlgebra(ZZ,-5,-19)
            Traceback (most recent call last):
            ...
            TypeError: base ring of quaternion algebra must be a field
        """
        ParentWithGens.__init__(self, base_ring, names=names)
        self._a = a
        self._b = b
        if not is_Field(base_ring):
            raise TypeError, "base ring of quaternion algebra must be a field"
        if is_RationalField(base_ring) and a.denominator() == 1 and b.denominator() == 1:
            element_constructor = quaternion_algebra_element.QuaternionAlgebraElement_rational_field
        elif is_NumberField(base_ring) and base_ring.degree() > 2 and base_ring.is_absolute() and \
                 a.denominator() == 1 and b.denominator() == 1 and base_ring.defining_polynomial().is_monic():
            # This QuaternionAlgebraElement_number_field class is not
            # designed to work with elements of a quadratic field.  To
            # do that, the main thing would be to implement
            # __getitem__, etc.  This would maybe give a factor of 2
            # (or more?) speedup.  Much care must be taken because the
            # underlying representation of quadratic fields is a bit
            # tricky.
            element_constructor = quaternion_algebra_element.QuaternionAlgebraElement_number_field
        else:
            element_constructor = quaternion_algebra_element.QuaternionAlgebraElement_generic
        self._populate_coercion_lists_(coerce_list=[base_ring], element_constructor=element_constructor)
        self._gens = [self([0,1,0,0]), self([0,0,1,0]), self([0,0,0,1])]

    def invariants(self):
        """
        Return the structural invariants a, b of this quaternion
        algebra.  Thus are the numbers that are squares of the first
        two generators.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(15)
            sage: Q.invariants()
            (-3, 5)
            sage: i^2
            -3
            sage: j^2
            5
        """
        return self._a, self._b

    def __reduce__(self):
        """
        Internal method used for pickling.

        TESTS::

            sage: QuaternionAlgebra(QQ,-1,-2).__reduce__()
            (<function unpickle_QuaternionAlgebra_v0 at ...>, (Rational Field, -1, -2, 'i,j,k'))

        Test uniqueness of parent::

            sage: Q = QuaternionAlgebra(QQ,-1,-2)
            sage: loads(dumps(Q)) is Q
            True
        """
        return unpickle_QuaternionAlgebra_v0, self._key

    def gen(self, i=0):
        """
        EXAMPLES::

            sage: Q.<ii,jj,kk> = QuaternionAlgebra(QQ,-1,-2); Q
            Quaternion Algebra (-1, -2) with base ring Rational Field
            sage: Q.gen(0)
            ii
            sage: Q.gen(1)
            jj
            sage: Q.gen(2)
            kk
            sage: Q.gens()
            [ii, jj, kk]
        """
        return self._gens[i]

    def _repr_(self):
        """
        TESTS::
            sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
            sage: type(Q)
            <class 'sage.algebras.quaternion_algebra.QuaternionAlgebra_ab'>
            sage: Q._repr_()
            'Quaternion Algebra (-5, -2) with base ring Rational Field'
            sage: Q
            Quaternion Algebra (-5, -2) with base ring Rational Field
            sage: print Q
            Quaternion Algebra (-5, -2) with base ring Rational Field
            sage: str(Q)
            'Quaternion Algebra (-5, -2) with base ring Rational Field'
        """
        return "Quaternion Algebra (%r, %r) with base ring %s"%(self._a, self._b, self.base_ring())

    def inner_product_matrix(self):
        """
        Return the inner product (or Gram) matrix associated to self,
        which is diagonal since this quaternion algebra is of the form
        (a,b).

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(-5,-19)
            sage: Q.inner_product_matrix()
            [  2   0   0   0]
            [  0  10   0   0]
            [  0   0  38   0]
            [  0   0   0 190]

            sage: R.<a,b> = QQ[]; Q.<i,j,k> = QuaternionAlgebra(Frac(R),a,b)
            sage: Q.inner_product_matrix()
            [    2     0     0     0]
            [    0  -2*a     0     0]
            [    0     0  -2*b     0]
            [    0     0     0 2*a*b]
        """
        a, b = self._a, self._b
        return diagonal_matrix(self.base_ring(), [2, -2*a, -2*b, 2*a*b])

    def discriminant(self):
        """
        Given a quaternion algebra A defined over the field of
        rational numbers, return the discriminant of A, i.e. the
        product of the ramified primes of A.

        EXAMPLES::

            sage: QuaternionAlgebra(210,-22).discriminant()
            210
            sage: QuaternionAlgebra(19).discriminant()
            19

        This raises a NotImplementedError if the base field is not
        the rational numbers.::

            sage: QuaternionAlgebra(QQ[sqrt(2)],3,19).discriminant()
            Traceback (most recent call last):
            ...
            NotImplementedError: base field must be rational numbers
        """
        try: return self.__discriminant
        except AttributeError: pass
        if not is_RationalField(self.base_ring()):
            raise NotImplementedError, "base field must be rational numbers"
        self.__discriminant = hilbert_conductor(self._a, self._b)
        return self.__discriminant


############################################################
# Unpickling
############################################################
def unpickle_QuaternionAlgebra_v0(*key):
    """
    The 0th version of pickling for quaternion algebras.

    EXAMPLES:
        sage: Q = QuaternionAlgebra(-5,-19)
        sage: f, t = Q.__reduce__()
        sage: sage.algebras.quaternion_algebra.unpickle_QuaternionAlgebra_v0(*t)
        Quaternion Algebra (-5, -19) with base ring Rational Field
        sage: loads(dumps(Q)) == Q
        True
        sage: loads(dumps(Q)) is Q
        True
    """
    return QuaternionAlgebra(*key)

############################################################
# Leftover code that we may want to use for something later.
## def ramified_primes_from_discs(D1, D2, T):
##     M = Integer(GCD([D1,D2,T]))
##     D3 = (T**2 - D1*D2)//4
##     facs = D3.factor()
##     D1 = fundamental_discriminant(D1)
##     D2 = fundamental_discriminant(D2)
##     D3 = fundamental_discriminant(D3)
##     ram_prms = []
##     for pow in facs:
##         p = pow[0]
##         if pow[1] % 2 == 1:
##             chi = (D.kronecker(p) for D in (D1,D2,D3))
##             if -1 in chi:
##                 ram_prms.append(p)
##             elif not 1 in chi and hilbert_symbol(D1,D3,p) == -1:
##                 ram_prms.append(p)
##         elif D1%p == 0 and D2%p == 0:
##             chi = (D.kronecker(p) for D in (D1,D2,D3))
##             if hilbert_symbol(D1,D3,p) == -1:
##                 ram_prms.append(p)
##     return ram_prms



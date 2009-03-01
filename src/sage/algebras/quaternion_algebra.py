"""
Quaternion Algebras

AUTHORS:

- David Kohel: original version in 2005 that some of this code is based on

- Jon Bobber -- 2009 rewrite

- William Stein -- 2009 rewrite

EXAMPLES::

This currently fails miserably
    sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
    sage: Q == loads(dumps(Q))
    True
"""
from sage.rings.arith import GCD, fundamental_discriminant, hilbert_symbol
from sage.rings.integer import Integer

from sage.rings.ring import Algebra
from sage.structure.parent_gens import ParentWithGens
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.sequence import Sequence
from sage.structure.element import is_Element

########################################################
# Constructor
########################################################

def QuaternionAlgebra(arg0, arg1=None, arg2=None, names='i,j,k', **kwds):
    """
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

          sage: QuaternionAlgebra(2)
          ??

    QuaternionAlgebra(D1, D2, T) - D1, D2, T integers; returns
    quaternion algebra Q<x, y> where Z[x] and Z[y] are quadratic
    suborders of discriminants D1, D2, respectively, and Z[x*y - y*x]
    is a quadratic suborder of discriminant D1*D2-T^2.::

        - QuaternionAlgebra(K, norms, traces) - quaternion algebra with
          inner product having given norms and traces.

        - QuaternionAlgebra(K, gram) - quaternion algebra with given
          Gram matrix.
    """
    # QuaternionAlgebra(D)
    if arg1 is None and arg2 is None:
        return QuaternionAlgebra_disc(arg0, names=names, **kwds)

    if arg2 is None:
        if is_Element(arg0):
            # QuaternionAlgebra(a, b)
            v = Sequence([arg0, arg1])
            K = v.universe().fraction_field()
            return QuaternionAlgebra(K, v[0], v[1])
        else:
            # QuaternionAlgebra(K, gram)
            return QuaternionAlgebra_gram(arg0, arg1, names=names, **kwds)

    # QuaternionAlgebra(D1, D2, T)
    if isinstance(arg0, (int, long, Integer)):
        return QuaternionAlgebra_D1D2T(arg0, arg1, arg2, names=names, **kwds)

    # QuaternionAlgebra(K, norms, traces)
    if isinstance(arg1, (list, tuple)):
        return QuaternionAlgebra_norms_traces(arg0, arg1, arg2, names=names, **kwds)

    # QuaternionAlgebra(K, a, b)
    return QuaternionAlgebra_ab(arg0, arg1, arg2, names=names, **kwds)


def QuaternionAlgebra_disc(D, names):
    """
    Return a rational quaternion algebra of discriminant D.
    """
    raise NotImplementedError

def QuaternionAlgebra_gram(K, gram, names):
    """
    Returns a quaternion algebra with given Gram matrix.
    """
    raise NotImplementedError

def QuaternionAlgebra_D1D2T(D1, D2, T, names):
    """
    Returns a quaternion algebra Q<x, y> where Z[x] and Z[y] are
    quadratic suborders of discriminants D1, D2, respectively, and
    Z[x*y - y*x] is a quadratic suborder of discriminant D1*D2-T^2.
    """
    raise NotImplementedError

def QuaternionAlgebra_norms_traces(D1, D2, T, names):
    """
    Returns a quaternion algebra with inner product having given norms
    and traces.
    """
    raise NotImplementedError


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

            sage: sage.algebras.quaternion_algebra.QuaternionAlgebra_abstract(QQ)
            Quaternion Algebra with base ring Rational Field
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

    def discriminant(self):
        """
        Given a quaternion algebra A defined over the field of rational
        numbers, return the discriminant of A, i.e. the product of the
        ramified primes of A.
        """
        raise NotImplementedError

    def gram_matrix(self):
        """
        The Gram matrix of the inner product determined by the norm.

        This is an alias for self.inner_product_matrix().

        EXAMPLES::

        """
        return self.inner_product_matrix()

    def inner_product_matrix(self):
        """
        Return the inner product matrix associated to self.

        EXAMPLES::

        """
        return self.vector_space().inner_product_matrix()

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

        At the moment, this is implemented only for finite fields.

        EXAMPLES::

            sage: Q.<i,j,k> = QuaternionAlgebra(GF(7),-5,-2)
            sage: Q.is_division_algebra()
            False
            sage: 1/(j+2*k)
            Traceback (most recent call last):
            ...
            AttributeError: The second operand must be a unit
            sage: QuaternionAlgebra(QQ,-5,-2).is_division_algebra()
            True
        """
        R = self.base_ring()
        if R.characteristic() != 0:
            return False
        # In characteristic 0 over field anything is invertible.  If
        # not over field, then stuff in base ring not invertible
        # already.  (NB -- we actually don't allow quaternion algebras
        # not over a field, so this is a little silly.  Anyway, this
        # is a good "just in case" test.)
        return R.is_field()

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

    def ramified_primes(self):
        """
        Return the set of ramified primes of this quaternion algebra.

        EXAMPLES::

        """
        try:
            return self._ramified_primes
        except:
            raise AttributeError, "Ramified primes have not been computed."

    def random_element(self, *args, **kwds):
        """
        Return a random element of this quaternion algebra.

        The args and kwds are passed to the random_element method of
        the base ring.

        EXAMPLES::

        """
        K = self.base_ring()
        return self([ K.random_element(*args, **kwds) for _ in range(4) ])

    def vector_space(self):
        """
        Return vector space with inner product associated to self.

        EXAMPLES::
        """
        try:
            return self.__vector_space
        except AttributeError:
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
            V = VectorSpace(K,4,inner_product_matrix = M)
            self.__vector_space = V
            return V
        raise NotImplementedError


from quaternion_algebra_element import QuaternionAlgebraElement_generic

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

        EXAMPLES::
        """
        ParentWithGens.__init__(self, base_ring, names)
        self._a = a
        self._b = b
        self._populate_coercion_lists_(coerce_list=[base_ring])
        self._gens = [self([0,1,0,0]), self([0,0,1,0]), self([0,0,0,1])]

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



########################################################
# Utility functions
########################################################

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

r"""
Symmetric Group Algebra
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from combinatorial_algebra import CombinatorialAlgebra
from free_module import CombinatorialFreeModule
from sage.categories.all import FiniteDimensionalAlgebrasWithBasis
from sage.combinat.permutation import Permutation, Permutations, Permutations_nk, PermutationOptions
import partition
from tableau import Tableau, StandardTableaux_size, StandardTableaux_shape, StandardTableaux
from sage.interfaces.all import gap
from sage.rings.all import QQ, PolynomialRing
from sage.rings.arith import factorial
from sage.matrix.all import matrix
from sage.modules.all import vector
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.categories.all import GroupAlgebras

permutation_options = PermutationOptions

def SymmetricGroupAlgebra(R, n):
    """
    Return the symmetric group algebra of order ``n`` over the ring ``R``.

    EXAMPLES::

        sage: QS3 = SymmetricGroupAlgebra(QQ, 3); QS3
        Symmetric group algebra of order 3 over Rational Field
        sage: QS3(1)
        [1, 2, 3]
        sage: QS3(2)
        2*[1, 2, 3]
        sage: basis = [QS3(p) for p in Permutations(3)]
        sage: a = sum(basis); a
        [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
        sage: a^2
        6*[1, 2, 3] + 6*[1, 3, 2] + 6*[2, 1, 3] + 6*[2, 3, 1] + 6*[3, 1, 2] + 6*[3, 2, 1]
        sage: a^2 == 6*a
        True
        sage: b = QS3([3, 1, 2])
        sage: b
        [3, 1, 2]
        sage: b*a
        [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
        sage: b*a == a
        True

    The canonical embedding from the symmetric group algebra of order
    `n` to the symmetric group algebra of order `p > n` is available as
    a coercion::

        sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
        sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
        sage: QS4.coerce_map_from(QS3)
        Generic morphism:
          From: Symmetric group algebra of order 3 over Rational Field
          To:   Symmetric group algebra of order 4 over Rational Field

        sage: x3  = QS3([3,1,2]) + 2 * QS3([2,3,1]); x3
        2*[2, 3, 1] + [3, 1, 2]
        sage: QS4(x3)
        2*[2, 3, 1, 4] + [3, 1, 2, 4]

    This allows for mixed expressions::

        sage: x4  = 3*QS4([3, 1, 4, 2])
        sage: x3 + x4
        2*[2, 3, 1, 4] + [3, 1, 2, 4] + 3*[3, 1, 4, 2]

        sage: QS0 = SymmetricGroupAlgebra(QQ, 0)
        sage: QS1 = SymmetricGroupAlgebra(QQ, 1)
        sage: x0 = QS0([])
        sage: x1 = QS1([1])
        sage: x0 * x1
        [1]
        sage: x3 - (2*x0 + x1) - x4
        -3*[1, 2, 3, 4] + 2*[2, 3, 1, 4] + [3, 1, 2, 4] - 3*[3, 1, 4, 2]

    Caveat: to achieve this, constructing ``SymmetricGroupAlgebra(QQ,
    10)`` currently triggers the construction of all symmetric group
    algebras of smaller order. Is this a feature we really want to have?

    .. WARNING::

        The semantics of multiplication in symmetric group algebras is
        determined by the order in which permutations are multiplied,
        which currently defaults to "in such a way that multiplication
        is associative with permutations acting on integers from the
        right", but can be changed to the opposite order at runtime
        by setting a global variable (see
        :meth:`sage.combinat.permutation.Permutations.global_options` ).
        In view of this, it is recommended that code not rely on the
        usual multiplication function, but rather use the methods
        :meth:`left_action_product` and :meth:`right_action_product`
        for multiplying permutations (these methods don't depend on the
        setting). See :trac:`14885` for more information.

    TESTS::

        sage: TestSuite(QS3).run()
    """
    return SymmetricGroupAlgebra_n(R,n)

class SymmetricGroupAlgebra_n(CombinatorialFreeModule):

    def __init__(self, R, n):
        """
        TESTS::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: TestSuite(QS3).run()
        """
        self.n = n
        self._name = "Symmetric group algebra of order %s"%self.n
        CombinatorialFreeModule.__init__(self, R, Permutations(n), prefix='', latex_prefix='', category = (GroupAlgebras(R),FiniteDimensionalAlgebrasWithBasis(R)))
        # This is questionable, and won't be inherited properly
        if n > 0:
            S = SymmetricGroupAlgebra(R, n-1)
            self.register_coercion(S.canonical_embedding(self))

    # _repr_ customization: output the basis element indexed by [1,2,3] as [1,2,3]
    _repr_option_bracket = False

    def group(self):
        """
        Return the underlying group.

        EXAMPLES::

            sage: SymmetricGroupAlgebra(QQ,4).group()
            Symmetric group of order 4! as a permutation group
        """
        return SymmetricGroup(self.n)

    @cached_method
    def one_basis(self):
        """
        Return the identity of the symmetric group, as per
        ``AlgebrasWithBasis.ParentMethods.one_basis``.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.one_basis()
            [1, 2, 3]
        """
        P = self.basis().keys()
        return P(range(1,self.n+1))

    def product_on_basis(self, left, right):
        """
        Return the product of the basis elements indexed by ``left`` and
        ``right``.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: p1 = Permutation([1,2,3])
            sage: p2 = Permutation([2,1,3])
            sage: QS3.product_on_basis(p1,p2)
            [2, 1, 3]
        """
        return self.monomial(left * right)

    def left_action_product(self, left, right):
        """
        Return the product of two elements ``left`` and ``right`` of
        ``self``, where multiplication is defined in such a way that
        for two permutations `p` and `q`, the product `pq` is the
        permutation obtained by first applying `q` and then applying
        `p`. This definition of multiplication is tailored to make
        multiplication of permutations associative with their action on
        numbers if permutations are to act on numbers from the left.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: p1 = Permutation([2, 1, 3])
            sage: p2 = Permutation([3, 1, 2])
            sage: QS3.left_action_product(QS3(p1), QS3(p2))
            [3, 2, 1]
            sage: x = QS3([1, 2, 3]) - 2*QS3([1, 3, 2])
            sage: y = 1/2 * QS3([3, 1, 2]) + 3*QS3([1, 2, 3])
            sage: QS3.left_action_product(x, y)
            3*[1, 2, 3] - 6*[1, 3, 2] - [2, 1, 3] + 1/2*[3, 1, 2]
            sage: QS3.left_action_product(0, x)
            0

        The method coerces its input into the algebra ``self``::

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: QS4.left_action_product(QS3([1, 2, 3]), QS3([2, 1, 3]))
            [2, 1, 3, 4]
            sage: QS4.left_action_product(1, Permutation([4, 1, 2, 3]))
            [4, 1, 2, 3]

        .. WARNING::

            Note that coercion presently works from permutations of ``n``
            into the ``n``-th symmetric group algebra, and also from all
            smaller symmetric group algebras into the ``n``-th symmetric
            group algebra, but not from permutations of integers smaller
            than ``n`` into the ``n``-th symmetric group algebra.
        """
        a = self(left)
        b = self(right)
        P = self.basis().keys()
        return self.sum_of_terms([(P([p[i-1] for i in q]), x * y)
                                  for (p, x) in a for (q, y) in b])
        # Why did we use P([p[i-1] for i in q])
        # instead of p.left_action_product(q) ?
        # Because having cast a and b into self, we already know that
        # p and q are permutations of the same number of elements,
        # and thus we don't need to waste our time on the input
        # sanitizing of left_action_product.

    def right_action_product(self, left, right):
        """
        Return the product of two elements ``left`` and ``right`` of
        ``self``, where multiplication is defined in such a way that
        for two permutations `p` and `q`, the product `pq` is the
        permutation obtained by first applying `p` and then applying
        `q`. This definition of multiplication is tailored to make
        multiplication of permutations associative with their action on
        numbers if permutations are to act on numbers from the right.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: p1 = Permutation([2, 1, 3])
            sage: p2 = Permutation([3, 1, 2])
            sage: QS3.right_action_product(QS3(p1), QS3(p2))
            [1, 3, 2]
            sage: x = QS3([1, 2, 3]) - 2*QS3([1, 3, 2])
            sage: y = 1/2 * QS3([3, 1, 2]) + 3*QS3([1, 2, 3])
            sage: QS3.right_action_product(x, y)
            3*[1, 2, 3] - 6*[1, 3, 2] + 1/2*[3, 1, 2] - [3, 2, 1]
            sage: QS3.right_action_product(0, x)
            0

        The method coerces its input into the algebra ``self``::

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: QS4.right_action_product(QS3([1, 2, 3]), QS3([2, 1, 3]))
            [2, 1, 3, 4]
            sage: QS4.right_action_product(1, Permutation([4, 1, 2, 3]))
            [4, 1, 2, 3]

        .. WARNING::

            Note that coercion presently works from permutations of ``n``
            into the ``n``-th symmetric group algebra, and also from all
            smaller symmetric group algebras into the ``n``-th symmetric
            group algebra, but not from permutations of integers smaller
            than ``n`` into the ``n``-th symmetric group algebra.
        """
        a = self(left)
        b = self(right)
        P = self.basis().keys()
        return self.sum_of_terms([(P([q[i-1] for i in p]), x * y)
                                  for (p, x) in a for (q, y) in b])
        # Why did we use P([q[i-1] for i in p])
        # instead of p.right_action_product(q) ?
        # Because having cast a and b into self, we already know that
        # p and q are permutations of the same number of elements,
        # and thus we don't need to waste our time on the input
        # sanitizing of right_action_product.

    def canonical_embedding(self, other):
        """
        Return the canonical embedding of ``self`` into ``other``.

        INPUT:

        - ``other`` -- a symmetric group algebra with order `p`
          satisfying `p \leq n` where `n` is the order of ``self``.

        EXAMPLES::

            sage: QS2 = SymmetricGroupAlgebra(QQ, 2)
            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: phi = QS2.canonical_embedding(QS4); phi
            Generic morphism:
              From: Symmetric group algebra of order 2 over Rational Field
              To:   Symmetric group algebra of order 4 over Rational Field

            sage: x = QS2([2,1]) + 2 * QS2([1,2])
            sage: phi(x)
            2*[1, 2, 3, 4] + [2, 1, 3, 4]

            sage: loads(dumps(phi))
            Generic morphism:
              From: Symmetric group algebra of order 2 over Rational Field
              To:   Symmetric group algebra of order 4 over Rational Field
        """
        if not isinstance(other, SymmetricGroupAlgebra_n) or self.n > other.n:
            raise ValueError("There is no canonical embedding from {0} to {1}".format(other, self))
        return self.module_morphism(other.monomial_from_smaller_permutation, codomain = other) # category = self.category() (currently broken)

    def monomial_from_smaller_permutation(self, permutation):
        """
        Convert ``permutation`` into a permutation, possibly extending it
        to the appropriate size, and return the corresponding basis
        element of ``self``.

        EXAMPLES::

            sage: QS5 = SymmetricGroupAlgebra(QQ, 5)
            sage: QS5.monomial_from_smaller_permutation([])
            [1, 2, 3, 4, 5]
            sage: QS5.monomial_from_smaller_permutation(Permutation([3,1,2]))
            [3, 1, 2, 4, 5]
            sage: QS5.monomial_from_smaller_permutation([5,3,4,1,2])
            [5, 3, 4, 1, 2]

        TESTS::

            sage: QS5.monomial_from_smaller_permutation([5,3,4,1,2]).parent()
            Symmetric group algebra of order 5 over Rational Field
        """
        P = self.basis().keys()
        return self.monomial( P(permutation) )

    def antipode(self, x):
        r"""
        Return the image of the element ``x`` of ``self`` under the
        antipode of the Hopf algebra ``self`` (where the
        comultiplication is the usual one on a group algebra).

        Explicitly, this is obtained by replacing each permutation
        `\sigma` by `\sigma^{-1}` in ``x`` while keeping all
        coefficients as they are.

        EXAMPLES::

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: QS4.antipode(2 * QS4([1, 3, 4, 2]) - 1/2 * QS4([1, 4, 2, 3]))
            -1/2*[1, 3, 4, 2] + 2*[1, 4, 2, 3]
            sage: all( QS4.antipode(QS4(p)) == QS4(p.inverse())
            ....:      for p in Permutations(4) )
            True

            sage: ZS3 = SymmetricGroupAlgebra(ZZ, 3)
            sage: ZS3.antipode(ZS3.zero())
            0
            sage: ZS3.antipode(-ZS3(Permutation([2, 3, 1])))
            -[3, 1, 2]
        """
        return self.sum_of_terms([(p.inverse(), coeff) for
                                  (p, coeff) in self(x)],
                                 distinct=True)

    def retract_plain(self, f, m):
        r"""
        Return the plain retract of the element `f \in R S_n`
        to `R S_m`, where `m \leq n` (and where `R S_n` is ``self``).

        If `m` is a nonnegative integer less or equal to `n`, then the
        plain retract from `S_n` to `S_m` is defined as an `R`-linear
        map `S_n \to S_m` which sends every permutation `p \in S_n`
        to

        .. MATH::

            \begin{cases} \mbox{pret}(p) &\mbox{if } \mbox{pret}(p)\mbox{ is defined;} \\
            0 & \mbox{otherwise} \end{cases}.

        Here `\mbox{pret}(p)` denotes the plain retract of the
        permutation `p` to `S_m`, which is defined in
        :meth:`~sage.combinat.permutation.Permutation.retract_plain`.

        EXAMPLES::

            sage: SGA3 = SymmetricGroupAlgebra(QQ, 3)
            sage: SGA3.retract_plain(2*SGA3([1,2,3]) - 4*SGA3([2,1,3]) + 7*SGA3([1,3,2]), 2)
            2*[1, 2] - 4*[2, 1]
            sage: SGA3.retract_plain(2*SGA3([1,3,2]) - 5*SGA3([2,3,1]), 2)
            0

            sage: SGA5 = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 4)
            11*[3, 2, 1, 4]
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 3)
            11*[3, 2, 1]
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 2)
            0
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 1)
            0

            sage: SGA5.retract_plain(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 3)
            8*[1, 2, 3] - 6*[1, 3, 2]
            sage: SGA5.retract_plain(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 1)
            8*[1]
            sage: SGA5.retract_plain(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 0)
            8*[]

        .. SEEALSO::

            :meth:`retract_direct_product`, :meth:`retract_okounkov_vershik`
        """
        RSm = SymmetricGroupAlgebra(self.base_ring(), m)
        pairs = []
        for (p, coeff) in f.monomial_coefficients().iteritems():
            p_ret = p.retract_plain(m)
            if p_ret is not None:
                pairs.append((p_ret, coeff))
        return RSm.sum_of_terms(pairs, distinct=True)

    def retract_direct_product(self, f, m):
        r"""
        Return the direct-product retract of the element `f \in R S_n`
        to `R S_m`, where `m \leq n` (and where `R S_n` is ``self``).

        If `m` is a nonnegative integer less or equal to `n`, then the
        direct-product retract from `S_n` to `S_m` is defined as an
        `R`-linear map `S_n \to S_m` which sends every permutation
        `p \in S_n` to

        .. MATH::

            \begin{cases} \mbox{dret}(p) &\mbox{if } \mbox{dret}(p)\mbox{ is defined;} \\
            0 & \mbox{otherwise} \end{cases}.

        Here `\mbox{dret}(p)` denotes the direct-product retract of the
        permutation `p` to `S_m`, which is defined in
        :meth:`~sage.combinat.permutation.Permutation.retract_direct_product`.

        EXAMPLES::

            sage: SGA3 = SymmetricGroupAlgebra(QQ, 3)
            sage: SGA3.retract_direct_product(2*SGA3([1,2,3]) - 4*SGA3([2,1,3]) + 7*SGA3([1,3,2]), 2)
            2*[1, 2] - 4*[2, 1]
            sage: SGA3.retract_direct_product(2*SGA3([1,3,2]) - 5*SGA3([2,3,1]), 2)
            0

            sage: SGA5 = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 4)
            11*[3, 2, 1, 4]
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 3)
            -6*[1, 3, 2] + 11*[3, 2, 1]
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 2)
            0
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 1)
            2*[1]

            sage: SGA5.retract_direct_product(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 3)
            8*[1, 2, 3] - 6*[1, 3, 2]
            sage: SGA5.retract_direct_product(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 1)
            2*[1]
            sage: SGA5.retract_direct_product(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 0)
            2*[]

        .. SEEALSO::

            :meth:`retract_plain`, :meth:`retract_okounkov_vershik`
        """
        RSm = SymmetricGroupAlgebra(self.base_ring(), m)
        dct = {}
        for (p, coeff) in f.monomial_coefficients().iteritems():
            p_ret = p.retract_direct_product(m)
            if not (p_ret is None):
                if not p_ret in dct.keys():
                    dct[p_ret] = coeff
                else:
                    dct[p_ret] += coeff
        return RSm._from_dict(dct)

    def retract_okounkov_vershik(self, f, m):
        r"""
        Return the Okounkov-Vershik retract of the element `f \in R S_n`
        to `R S_m`, where `m \leq n` (and where `R S_n` is ``self``).

        If `m` is a nonnegative integer less or equal to `n`, then the
        Okounkov-Vershik retract from `S_n` to `S_m` is defined as an
        `R`-linear map `S_n \to S_m` which sends every permutation
        `p \in S_n` to the Okounkov-Vershik retract of the permutation
        `p` to `S_m`, which is defined in
        :meth:`~sage.combinat.permutation.Permutation.retract_okounkov_vershik`.

        EXAMPLES::

            sage: SGA3 = SymmetricGroupAlgebra(QQ, 3)
            sage: SGA3.retract_okounkov_vershik(2*SGA3([1,2,3]) - 4*SGA3([2,1,3]) + 7*SGA3([1,3,2]), 2)
            9*[1, 2] - 4*[2, 1]
            sage: SGA3.retract_okounkov_vershik(2*SGA3([1,3,2]) - 5*SGA3([2,3,1]), 2)
            2*[1, 2] - 5*[2, 1]

            sage: SGA5 = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 4)
            -6*[1, 3, 2, 4] + 8*[1, 4, 2, 3] + 11*[3, 2, 1, 4]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 3)
            2*[1, 3, 2] + 11*[3, 2, 1]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 2)
            13*[1, 2]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 1)
            13*[1]

            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 3)
            8*[1, 2, 3] - 6*[1, 3, 2]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 1)
            2*[1]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 0)
            2*[]

        .. SEEALSO::

            :meth:`retract_plain`, :meth:`retract_direct_product`
        """
        RSm = SymmetricGroupAlgebra(self.base_ring(), m)
        dct = {}
        for (p, coeff) in f.monomial_coefficients().iteritems():
            p_ret = p.retract_okounkov_vershik(m)
            if not p_ret in dct.keys():
                dct[p_ret] = coeff
            else:
                dct[p_ret] += coeff
        return RSm._from_dict(dct)

#     def _coerce_start(self, x):
#         """
#         Coerce things into the symmetric group algebra.

#         EXAMPLES::

#             sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
#             sage: QS3._coerce_start([])
#             [1, 2, 3]
#             sage: QS3._coerce_start([2,1])
#             [2, 1, 3]
#             sage: _.parent()
#             Symmetric group algebra of order 3 over Rational Field
#         """
#         if x == []:
#             return self( self._one )
#         if len(x) < self.n and x in permutation.Permutations():
#             return self( list(x) + range(len(x)+1, self.n+1) )
#         raise TypeError

    def cpis(self):
        """
        Return a list of the centrally primitive idempotents of
        ``self``.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: a = QS3.cpis()
            sage: a[0]  # [3]
            1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
            sage: a[1]  # [2, 1]
            2/3*[1, 2, 3] - 1/3*[2, 3, 1] - 1/3*[3, 1, 2]
        """
        return [self.cpi(p) for p in partition.Partitions_n(self.n)]

    def cpi(self, p):
        """
        Return the centrally primitive idempotent for the symmetric group
        of order `n` corresponding to the irreducible representation
        indexed by the partition ``p``.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3.cpi([2,1])
            2/3*[1, 2, 3] - 1/3*[2, 3, 1] - 1/3*[3, 1, 2]
            sage: QS3.cpi([3])
            1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
            sage: QS3.cpi([1,1,1])
            1/6*[1, 2, 3] - 1/6*[1, 3, 2] - 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] - 1/6*[3, 2, 1]

            sage: QS0 = SymmetricGroupAlgebra(QQ, 0)
            sage: QS0.cpi(Partition([]))
            []

        TESTS::

            sage: QS3.cpi([2,2])
            Traceback (most recent call last):
            ...
            TypeError: p (= [2, 2]) must be a partition of n (= 3)
        """
        if p not in partition.Partitions_n(self.n):
            raise TypeError("p (= {p}) must be a partition of n (= {n})".format(p=p, n=self.n))

        character_table = eval(gap.eval("Display(Irr(SymmetricGroup(%d)));"%self.n))

        np = partition.Partitions_n(self.n).list()
        np.reverse()
        p_index = np.index(p)

        big_coeff = character_table[p_index][0] / factorial(self.n)

        character_row = character_table[p_index]
        P = self.basis().keys()
        dct = { g : big_coeff * character_row[np.index(g.cycle_type())] for g in P }

        return self._from_dict(dct)

    def algebra_generators(self):
        r"""
        Return generators of this group algebra (as algebra) as a
        list of permutations.

        The generators used for the group algebra of `S_n` are the
        transposition `(2, 1)` and the `n`-cycle `(1, 2, \ldots, n)`,
        unless `n \leq 1` (in which case no generators are needed).

        EXAMPLES::

            sage: SymmetricGroupAlgebra(ZZ,5).algebra_generators()
            [[2, 1, 3, 4, 5], [2, 3, 4, 5, 1]]

            sage: SymmetricGroupAlgebra(QQ,0).algebra_generators()
            []

            sage: SymmetricGroupAlgebra(QQ,1).algebra_generators()
            []

        TESTS:

        Check that :trac:`15309` is fixed::

            sage: S3 = SymmetricGroupAlgebra(QQ, 3)
            sage: S3.algebra_generators()
            [[2, 1, 3], [2, 3, 1]]
            sage: C = CombinatorialFreeModule(ZZ, ZZ)
            sage: M = C.module_morphism(lambda x: S3.zero(), codomain=S3)
            sage: M.register_as_coercion()
        """
        if self.n <= 1:
            return []
        a = range(1, self.n+1)
        a[0] = 2
        a[1] = 1
        b = range(2, self.n+2)
        b[self.n-1] = 1
        return [self.monomial(self._basis_keys(a)), self.monomial(self._basis_keys(b))]

    def _conjugacy_classes_representatives_underlying_group(self):
        r"""
        Return a complete list of representatives of conjugacy
        classes of the underlying symmetric group.

        EXAMPLES::

            sage: SG=SymmetricGroupAlgebra(ZZ,3)
            sage: SG._conjugacy_classes_representatives_underlying_group()
            [[2, 3, 1], [2, 1, 3], [1, 2, 3]]
        """
        P = self.basis().keys()
        return [P.element_in_conjugacy_classes(nu) for nu in partition.Partitions(self.n)]

    def rsw_shuffling_element(self, k):
        r"""
        Return the `k`-th Reiner-Saliola-Welker shuffling element in
        the group algebra ``self``.

        The `k`-th Reiner-Saliola-Welker shuffling element in the
        symmetric group algebra `R S_n` over a ring `R` is defined as the
        sum `\sum_{\sigma \in S_n} \mathrm{noninv}_k(\sigma) \cdot \sigma`,
        where for every permutation `\sigma`, the number
        `\mathrm{noninv}_k(\sigma)` is the number of all
        `k`-noninversions of `\sigma` (that is, the number of all
        `k`-element subsets of `\{ 1, 2, \ldots, n \}` on which
        `\sigma` restricts to a strictly increasing map). See
        :meth:`sage.combinat.permutation.number_of_noninversions` for
        the `\mathrm{noninv}` map.

        This element is more or less the operator `\nu_{k, 1^{n-k}}`
        introduced in [RSW2011]_; more precisely, `\nu_{k, 1^{n-k}}`
        is the left multiplication by this element.

        It is a nontrivial theorem (Theorem 1.1 in [RSW2011]_) that
        the operators `\nu_{k, 1^{n-k}}` (for fixed `n` and varying
        `k`) pairwise commute. It is a conjecture (Conjecture 1.2 in
        [RSW2011]_) that all their eigenvalues are integers (which, in
        light of their commutativity and easily established symmetry,
        yields that they can be simultaneously diagonalized over `\QQ`
        with only integer eigenvalues).

        EXAMPLES:

        The Reiner-Saliola-Welker shuffling elements on `\QQ S_3`::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.rsw_shuffling_element(0)
            [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            sage: QS3.rsw_shuffling_element(1)
            3*[1, 2, 3] + 3*[1, 3, 2] + 3*[2, 1, 3] + 3*[2, 3, 1] + 3*[3, 1, 2] + 3*[3, 2, 1]
            sage: QS3.rsw_shuffling_element(2)
            3*[1, 2, 3] + 2*[1, 3, 2] + 2*[2, 1, 3] + [2, 3, 1] + [3, 1, 2]
            sage: QS3.rsw_shuffling_element(3)
            [1, 2, 3]
            sage: QS3.rsw_shuffling_element(4)
            0

        Checking the commutativity of Reiner-Saliola-Welker shuffling
        elements (we leave out the ones for which it is trivial)::

            sage: def test_rsw_comm(n):
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     rsws = [QSn.rsw_shuffling_element(k) for k in range(2, n)]
            ....:     return all( all( rsws[i] * rsws[j] == rsws[j] * rsws[i]
            ....:                      for j in range(i) )
            ....:                 for i in range(len(rsws)) )
            sage: test_rsw_comm(3)
            True
            sage: test_rsw_comm(4)
            True
            sage: test_rsw_comm(5)   # long time
            True

        .. NOTE::

            For large ``k`` (relative to ``n``), it might be faster to call
            ``QSn.left_action_product(QSn.semi_rsw_element(k), QSn.antipode(binary_unshuffle_sum(k)))``
            than ``QSn.rsw_shuffling_element(n)``.

        .. SEEALSO::

            :meth:`semi_rsw_element`, :meth:`binary_unshuffle_sum`
        """
        P = self.basis().keys()
        return self.sum_of_terms([(p, p.number_of_noninversions(k)) for p in P],
                                 distinct=True)

    def semi_rsw_element(self, k):
        r"""
        Return the `k`-th semi-RSW element in the group algebra ``self``.

        The `k`-th semi-RSW element in the symmetric group algebra
        `R S_n` over a ring `R` is defined as the sum of all permutations
        `\sigma \in S_n` satisfying
        `\sigma(1) < \sigma(2) < \cdots < \sigma(k)`.

        This element has the property that, if it is denoted by `s_k`,
        then `s_k S(s_k)` is `(n-k)!` times the `k`-th
        Reiner-Saliola-Welker shuffling element of `R S_n` (see
        :meth:`rsw_shuffling_element`). Here, `S` denotes the antipode
        of the group algebra `R S_n`.

        The `k`-th semi-RSW element is the image of the complete
        non-commutative symmetric function `S^{(k, 1^{n-k})}` in the
        ring of non-commutative symmetric functions under the canonical
        projection on the symmetric group algebra (through the descent
        algebra).

        EXAMPLES:

        The semi-RSW elements on `\QQ S_3`::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.semi_rsw_element(0)
            [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            sage: QS3.semi_rsw_element(1)
            [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            sage: QS3.semi_rsw_element(2)
            [1, 2, 3] + [1, 3, 2] + [2, 3, 1]
            sage: QS3.semi_rsw_element(3)
            [1, 2, 3]
            sage: QS3.semi_rsw_element(4)
            0

        Let us check the relation with the `k`-th Reiner-Saliola-Welker
        shuffling element stated in the docstring::

            sage: def test_rsw(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     for k in range(1, n):
            ....:         a = ZSn.semi_rsw_element(k)
            ....:         b = ZSn.left_action_product(a, ZSn.antipode(a))
            ....:         if factorial(n-k) * ZSn.rsw_shuffling_element(k) != b:
            ....:             return False
            ....:     return True
            sage: test_rsw(3)
            True
            sage: test_rsw(4)
            True
            sage: test_rsw(5)  # long time
            True

        Let us also check the statement about the complete
        non-commutative symmetric function::

            sage: def test_rsw_ncsf(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     NSym = NonCommutativeSymmetricFunctions(ZZ)
            ....:     S = NSym.S()
            ....:     for k in range(1, n):
            ....:         a = S(Composition([k] + [1]*(n-k))).to_symmetric_group_algebra()
            ....:         if a != ZSn.semi_rsw_element(k):
            ....:             return False
            ....:     return True
            sage: test_rsw_ncsf(3)
            True
            sage: test_rsw_ncsf(4)
            True
            sage: test_rsw_ncsf(5)  # long time
            True
        """
        n = self.n
        if n < k:
            return self.zero()
        def complement(xs):
            res = range(1, n+1)
            for x in xs:
                res.remove(x)
            return res
        P = Permutations()
        return self.sum_of_monomials([P(complement(q) + list(q))
                                      for q in Permutations_nk(n, n-k)])

    def binary_unshuffle_sum(self, k):
        r"""
        Return the `k`-th binary unshuffle sum in the group algebra
        ``self``.

        The `k`-th binary unshuffle sum in the symmetric group algebra
        `R S_n` over a ring `R` is defined as the sum of all permutations
        `\sigma \in S_n` satisfying
        `\sigma(1) < \sigma(2) < \cdots < \sigma(k)` and
        `\sigma(k+1) < \sigma(k+2) < \cdots < \sigma(n)`.

        This element has the property that, if it is denoted by `t_k`,
        and if the `k`-th semi-RSW element (see :meth:`semi_rsw_element`)
        is denoted by `s_k`, then `s_k S(t_k)` and `t_k S(s_k)` both
        equal the `k`-th Reiner-Saliola-Welker shuffling element of
        `R S_n` (see :meth:`rsw_shuffling_element`).

        The `k`-th binary unshuffle sum is the image of the complete
        non-commutative symmetric function `S^{(k, n-k)}` in the
        ring of non-commutative symmetric functions under the canonical
        projection on the symmetric group algebra (through the descent
        algebra).

        EXAMPLES:

        The binary unshuffle sums on `\QQ S_3`::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.binary_unshuffle_sum(0)
            [1, 2, 3]
            sage: QS3.binary_unshuffle_sum(1)
            [1, 2, 3] + [2, 1, 3] + [3, 1, 2]
            sage: QS3.binary_unshuffle_sum(2)
            [1, 2, 3] + [1, 3, 2] + [2, 3, 1]
            sage: QS3.binary_unshuffle_sum(3)
            [1, 2, 3]
            sage: QS3.binary_unshuffle_sum(4)
            0

        Let us check the relation with the `k`-th Reiner-Saliola-Welker
        shuffling element stated in the docstring::

            sage: def test_rsw(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     for k in range(1, n):
            ....:         a = ZSn.semi_rsw_element(k)
            ....:         b = ZSn.binary_unshuffle_sum(k)
            ....:         c = ZSn.left_action_product(a, ZSn.antipode(b))
            ....:         d = ZSn.left_action_product(b, ZSn.antipode(a))
            ....:         e = ZSn.rsw_shuffling_element(k)
            ....:         if c != e or d != e:
            ....:             return False
            ....:     return True
            sage: test_rsw(3)
            True
            sage: test_rsw(4)  # long time
            True
            sage: test_rsw(5)  # long time
            True

        Let us also check the statement about the complete
        non-commutative symmetric function::

            sage: def test_rsw_ncsf(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     NSym = NonCommutativeSymmetricFunctions(ZZ)
            ....:     S = NSym.S()
            ....:     for k in range(1, n):
            ....:         a = S(Composition([k, n-k])).to_symmetric_group_algebra()
            ....:         if a != ZSn.binary_unshuffle_sum(k):
            ....:             return False
            ....:     return True
            sage: test_rsw_ncsf(3)
            True
            sage: test_rsw_ncsf(4)
            True
            sage: test_rsw_ncsf(5)  # long time
            True
        """
        n = self.n
        if n < k:
            return self.zero()
        def complement(xs):
            res = range(1, n+1)
            for x in xs:
                res.remove(x)
            return res
        from sage.combinat.subset import Subsets
        P = Permutations()
        return self.sum_of_monomials([P(sorted(q) + complement(q)) for q in Subsets(n, k)])

    def jucys_murphy(self, k):
        r"""
        Return the Jucys-Murphy element `J_k` (also known as a
        Young-Jucys-Murphy element) for the symmetric group
        algebra ``self``.

        The Jucys-Murphy element `J_k` in the symmetric group algebra
        `R S_n` is defined for every `k \in \{ 1, 2, \ldots, n \}` by

        .. MATH::

            J_k = (1, k) + (2, k) + \cdots + (k-1, k) \in R S_n,

        where the addends are transpositions in `S_n` (regarded as
        elements of `R S_n`). We note that there is not a dependence on `n`,
        so it is often surpressed in the notation.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.jucys_murphy(1)
            0
            sage: QS3.jucys_murphy(2)
            [2, 1, 3]
            sage: QS3.jucys_murphy(3)
            [1, 3, 2] + [3, 2, 1]

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: j3 = QS4.jucys_murphy(3); j3
            [1, 3, 2, 4] + [3, 2, 1, 4]
            sage: j4 = QS4.jucys_murphy(4); j4
            [1, 2, 4, 3] + [1, 4, 3, 2] + [4, 2, 3, 1]
            sage: j3*j4 == j4*j3
            True

            sage: QS5 = SymmetricGroupAlgebra(QQ, 5)
            sage: QS5.jucys_murphy(4)
            [1, 2, 4, 3, 5] + [1, 4, 3, 2, 5] + [4, 2, 3, 1, 5]

        TESTS::

            sage: QS3.jucys_murphy(4)
            Traceback (most recent call last):
            ...
            ValueError: k (= 4) must be between 1 and n (= 3) (inclusive)
        """
        if k < 1 or k > self.n:
            raise ValueError("k (= {k}) must be between 1 and n (= {n}) (inclusive)".format(k=k, n=self.n))

        res = self.zero()

        for i in range(1, k):
            p = range(1, self.n+1)
            p[i-1] = k
            p[k-1] = i
            res += self.monomial(self._basis_keys(p))
        return res



    def seminormal_basis(self, mult='l2r'):
        r"""
        Return a list of the seminormal basis elements of ``self``.

        The seminormal basis of a symmetric group algebra is defined as
        follows:

        Let `n` be a nonnegative integer. Let `R` be a `\QQ`-algebra.
        In the following, we will use the "left action" convention for
        multiplying permutations. This means that for all permutations
        `p` and `q` in `S_n`, the product `pq` is defined in such a way
        that `(pq)(i) = p(q(i))` for each `i \in \{ 1, 2, \ldots, n \}`
        (this is the same convention as in :meth:`left_action_product`,
        but not the default semantics of the `*` operator on
        permutations in Sage). Thus, for instance, `s_2 s_1` is the
        permutation obtained by first transposing `1` with `2` and
        then transposing `2` with `3` (where `s_i = (i, i+1)`).

        For every partition `\lambda` of `n`, let

        .. MATH::

            \kappa_{\lambda} = \frac{n!}{f^{\lambda}}

        where `f^{\lambda}` is the number of standard Young tableaux
        of shape `\lambda`. Note that `\kappa_{\lambda}` is an integer,
        namely the product of all hook lengths of `\lambda` (by the
        hook length formula). In Sage, this integer can be computed by
        using :func:`sage.combinat.symmetric_group_algebra.kappa()`.

        Let `T` be a standard tableau of size `n`.

        Let `a(T)` denote the formal sum (in `R S_n`) of all
        permutations in `S_n` which stabilize the rows of `T` (as
        sets), i. e., which map each entry `i` of `T` to an entry in
        the same row as `i`. (See
        :func:`sage.combinat.symmetric_group_algebra.a()` for
        an implementation of this.)

        Let `b(T)` denote the signed formal sum (in `R S_n`) of all
        permutations in `S_n` which stabilize the columns of `T` (as
        sets). Here, "signed" means that each permutation is
        multiplied with its sign. (This is implemented in
        :func:`sage.combinat.symmetric_group_algebra.b()`.)

        Define an element `e(T)` of `R S_n` to be `a(T) b(T)`. (This
        is implemented in :func:`sage.combinat.symmetric_group_algebra.e()`
        for `R = \QQ`.)

        Let `\mathrm{sh}(T)` denote the shape of `T`.
        (See :meth:`~sage.combinat.tableau.Tableau.shape`.)

        Let `\overline{T}` denote the standard tableau of size `n-1`
        obtained by removing the letter `n` (along with its cell) from
        `T` (if `n \geq 1`).

        Now, we define an element `\epsilon(T)` of `R S_n`. We define
        it by induction on the size `n` of `T`, so we set
        `\epsilon(\emptyset) = 1` and only need to define `\epsilon(T)`
        for `n \geq 1`, assuming that `\epsilon(\overline{T})` is
        already defined. We do this by setting

        .. MATH::

            \epsilon(T) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                          \epsilon(\overline{T})
                          e(T) \epsilon(\overline{T}).

        This element `\epsilon(T)` is implemented as
        :func:`sage.combinat.symmetric_group_algebra.epsilon` for
        `R = \QQ`, but it is also a particular case of the elements
        `\epsilon(T, S)` defined below.

        Now let `S` be a further tableau of the same shape as `T`
        (possibly equal to `T`). Let `\pi_{T, S}` denote the
        permutation in `S_n` such that applying this permutation to
        the entries of `T` yields the tableau `S`. Define an element
        `\epsilon(T, S)` of `R S_n` by

        .. MATH::

            \epsilon(T, S) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                             \epsilon(\overline S) \pi_{T, S}
                             e(T) \epsilon(\overline T)
                           = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                             \epsilon(\overline S) a(S) \pi_{T, S}
                             b(T) \epsilon(\overline T).

        This element `\epsilon(T, S)` is called *Young's seminormal
        unit corresponding to the bitableau `(T, S)`*, and is the
        return value of :meth:`epsilon_ik` applied to ``T`` and
        ``S``. Note that `\epsilon(T, T) = \epsilon(T)`.

        If we let `\lambda` run through all partitions of `n`, and
        `(T, S)` run through all pairs of tableaux of shape
        `\lambda`, then the elements `\epsilon(T, S)` form a basis
        of `R S_n`. This basis is called *Young's seminormal basis*
        and has the properties that

        .. MATH::

            \epsilon(T, S) \epsilon(U, V) = \delta_{T, V} \epsilon(U, S)

        (where `\delta` stands for the Kronecker delta).

        .. WARNING::

            Because of our convention, we are multiplying our elements in
            reverse of those given in some papers, for example [Ram1997]_.
            Using the other convention of multiplying permutations, we would
            instead have
            `\epsilon(U, V) \epsilon(T, S) = \delta_{T, V} \epsilon(U, S)`.

        In other words, Young's seminormal basis consists of the matrix
        units in a (particular) Artin-Wedderburn decomposition of `R S_n`
        into a direct product of matrix algebras over `\QQ`.

        The output of :meth:`seminormal_basis` is a list of all
        elements of the seminormal basis of ``self``.

        INPUT:

        - ``mult`` -- string (default: ``'l2r'``). If set to ``'r2l'``,
          this causes the method to return the list of the
          antipodes (:meth:`antipode`) of all `\epsilon(T, S)`
          instead of the `\epsilon(T, S)` themselves.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3.seminormal_basis()
            [1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1],
            1/3*[1, 2, 3] + 1/6*[1, 3, 2] - 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] + 1/6*[3, 2, 1],
            1/3*[1, 3, 2] + 1/3*[2, 3, 1] - 1/3*[3, 1, 2] - 1/3*[3, 2, 1],
            1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1],
            1/3*[1, 2, 3] - 1/6*[1, 3, 2] + 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] - 1/6*[3, 2, 1],
            1/6*[1, 2, 3] - 1/6*[1, 3, 2] - 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] - 1/6*[3, 2, 1]]

        REFERENCES:

        .. [Ram1997] Arun Ram. *Seminormal representations of Weyl groups
           and Iwahori-Hecke algebras*. Proc. London Math. Soc. (3)
           **75** (1997). 99-133. :arxiv:`math/9511223v1`.
           http://www.ms.unimelb.edu.au/~ram/Publications/1997PLMSv75p99.pdf
        """
        basis = []
        for part in partition.Partitions_n(self.n):
            stp = StandardTableaux_shape(part)
            for t1 in stp:
                for t2 in stp:
                    basis.append(self.epsilon_ik(t1, t2, mult=mult))
        return basis


    def dft(self, form="seminormal", mult='l2r'):
        """
        Return the discrete Fourier transform for ``self``.

        INPUT:

        - ``mult`` -- string (default: `l2r`). If set to `r2l`,
          this causes the method to use the antipodes
          (:meth:`antipode`) of the seminormal basis instead of
          the seminormal basis.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.dft()
            [   1    1    1    1    1    1]
            [   1  1/2   -1 -1/2 -1/2  1/2]
            [   0  3/4    0  3/4 -3/4 -3/4]
            [   0    1    0   -1    1   -1]
            [   1 -1/2    1 -1/2 -1/2 -1/2]
            [   1   -1   -1    1    1   -1]
        """
        if form == "seminormal":
            return self._dft_seminormal(mult=mult)
        else:
            raise ValueError("invalid form (= %s)"%form)

    def _dft_seminormal(self, mult='l2r'):
        """
        Return the seminormal form of the discrete Fourier for ``self``.

        INPUT:

        - ``mult`` -- string (default: `l2r`). If set to `r2l`,
          this causes the method to use the antipodes
          (:meth:`antipode`) of the seminormal basis instead of
          the seminormal basis.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3._dft_seminormal()
            [   1    1    1    1    1    1]
            [   1  1/2   -1 -1/2 -1/2  1/2]
            [   0  3/4    0  3/4 -3/4 -3/4]
            [   0    1    0   -1    1   -1]
            [   1 -1/2    1 -1/2 -1/2 -1/2]
            [   1   -1   -1    1    1   -1]

        .. SEEALSO::

            :meth:`seminormal_basis`
        """
        snb = self.seminormal_basis(mult=mult)
        return matrix( [vector(b) for b in snb] ).inverse().transpose()

    def epsilon_ik(self, itab, ktab, star=0, mult='l2r'):
        r"""
        Return the seminormal basis element of ``self`` corresponding to the
        pair of tableaux ``itab`` and ``ktab`` (or restrictions of these
        tableaux, if the optional variable ``star`` is set).

        INPUT:

        - ``itab``, ``ktab`` -- two standard tableaux of size `n`.

        - ``star`` -- integer (default: `0`).

        - ``mult`` -- string (default: `l2r`). If set to `r2l`,
          this causes the method to return the antipode
          (:meth:`antipode`) of `\epsilon(I, K)` instead of
          `\epsilon(I, K)` itself.

        OUTPUT:

        The element `\epsilon(I, K)`, where `I` and `K` are the tableaux
        obtained by removing all entries higher than `n - \mathrm{star}`
        from ``itab`` and ``ktab``, respectively. Here, we are using the
        notations from :meth:`seminormal_basis`.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = QS3.epsilon_ik([[1,2,3]], [[1,2,3]]); a
            1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
            sage: QS3.dft()*vector(a)
            (1, 0, 0, 0, 0, 0)
            sage: a = QS3.epsilon_ik([[1,2],[3]], [[1,2],[3]]); a
            1/3*[1, 2, 3] - 1/6*[1, 3, 2] + 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] - 1/6*[3, 2, 1]
            sage: QS3.dft()*vector(a)
            (0, 0, 0, 0, 1, 0)

        Let us take some properties of the seminormal basis listed in
        the docstring of :meth:`seminormal_basis`, and verify them on
        the situation of `S_3`.

        First, check the formula

        .. MATH::

            \epsilon(T) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                          \epsilon(\overline{T})
                          e(T) \epsilon(\overline{T}).

        In fact::

            sage: from sage.combinat.symmetric_group_algebra import e
            sage: def test_sn1(n):
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     QSn1 = SymmetricGroupAlgebra(QQ, n - 1)
            ....:     for T in StandardTableaux(n):
            ....:         TT = T.restrict(n-1)
            ....:         eTT = QSn1.epsilon_ik(TT, TT)
            ....:         eT = QSn.epsilon_ik(T, T)
            ....:         kT = prod(T.shape().hooks())
            ....:         if kT * eT != eTT * e(T) * eTT:
            ....:             return False
            ....:     return True
            sage: test_sn1(3)
            True
            sage: test_sn1(4)   # long time
            True

        Next, we check the identity

        .. MATH::

            \epsilon(T, S) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                             \epsilon(\overline S) \pi_{T, S}
                             e(T) \epsilon(\overline T)

        which we used to define `\epsilon(T, S)`. In fact::

            sage: from sage.combinat.symmetric_group_algebra import e
            sage: def test_sn2(n):
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     mul = QSn.left_action_product
            ....:     QSn1 = SymmetricGroupAlgebra(QQ, n - 1)
            ....:     for lam in Partitions(n):
            ....:         k = prod(lam.hooks())
            ....:         for T in StandardTableaux(lam):
            ....:             for S in StandardTableaux(lam):
            ....:                 TT = T.restrict(n-1)
            ....:                 SS = S.restrict(n-1)
            ....:                 eTT = QSn1.epsilon_ik(TT, TT)
            ....:                 eSS = QSn1.epsilon_ik(SS, SS)
            ....:                 eTS = QSn.epsilon_ik(T, S)
            ....:                 piTS = [0] * n
            ....:                 for (i, j) in T.cells():
            ....:                     piTS[T[i][j] - 1] = S[i][j]
            ....:                 piTS = QSn(Permutation(piTS))
            ....:                 if k * eTS != mul(mul(eSS, piTS), mul(e(T), eTT)):
            ....:                     return False
            ....:     return True
            sage: test_sn2(3)
            True
            sage: test_sn2(4)   # long time
            True

        Let us finally check the identity

        .. MATH::

            \epsilon(T, S) \epsilon(U, V) = \delta_{T, V} \epsilon(U, S)

        In fact::

            sage: def test_sn3(lam):
            ....:     n = lam.size()
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     mul = QSn.left_action_product
            ....:     for T in StandardTableaux(lam):
            ....:         for S in StandardTableaux(lam):
            ....:             for U in StandardTableaux(lam):
            ....:                 for V in StandardTableaux(lam):
            ....:                     lhs = mul(QSn.epsilon_ik(T, S), QSn.epsilon_ik(U, V))
            ....:                     if T == V:
            ....:                         rhs = QSn.epsilon_ik(U, S)
            ....:                     else:
            ....:                         rhs = QSn.zero()
            ....:                     if rhs != lhs:
            ....:                         return False
            ....:     return True
            sage: all( test_sn3(lam) for lam in Partitions(3) )
            True
            sage: all( test_sn3(lam) for lam in Partitions(4) )   # long time
            True
        """
        it = Tableau(itab)
        kt = Tableau(ktab)

        stn = StandardTableaux_size(self.n)

        if it not in stn:
            raise TypeError("it must be a standard tableau of size %s"%self.n)

        if kt not in stn:
            raise TypeError("kt must be a standard tableau of size %s"%self.n)

        if it.shape() != kt.shape():
            raise ValueError("it and kt must be of the same shape")

        BR = self.base_ring()
        z_elts = {}
        epik = epsilon_ik(it, kt, star=star)
        for m,c in epik._monomial_coefficients.iteritems():
            z_elts[m] = BR(c)
        z = self._from_dict(z_elts)

        if mult == 'l2r':
            return z
        else:
            return z.map_support(lambda x: x.inverse())


epsilon_ik_cache = {}
def epsilon_ik(itab, ktab, star=0):
    """
    Return the seminormal basis element of the symmetric group
    algebra `\QQ S_n` corresponding to the pair of tableaux
    ``itab`` and ``ktab`` (or restrictions of these tableaux,
    if the optional variable ``star`` is set).

    INPUT:

    - ``itab``, ``ktab`` -- two standard tableaux of same size.

    - ``star`` -- integer (default: `0`).

    OUTPUT:

    The element `\epsilon(I, K) \in \QQ S_n`, where `I` and `K`
    are the tableaux obtained by removing all entries higher
    than `n - \mathrm{star}` from ``itab`` and ``ktab``,
    respectively (where `n` is the size of ``itab`` and
    ``ktab``). Here, we are using the notations from
    :meth:`~sage.combinat.symmetric_group_algebra.SymmetricGroupAlgebra_n.seminormal_basis`.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import epsilon_ik
        sage: epsilon_ik([[1,2],[3]], [[1,3],[2]])
        1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1]
        sage: epsilon_ik([[1,2],[3]], [[1,3],[2]], star=1)
        Traceback (most recent call last):
        ...
        ValueError: the two tableaux must be of the same shape
    """
    it = Tableau(itab)
    kt = Tableau(ktab)
    if star:
        it = it.restrict(it.size() - star)
        kt = kt.restrict(kt.size() - star)

    if it.shape() != kt.shape():
        raise ValueError("the two tableaux must be of the same shape")

    if kt == it:
        res = epsilon(itab)
    elif (it, kt) in epsilon_ik_cache:
        res =  epsilon_ik_cache[(it,kt)]
    else:
        eik = e_ik(it, kt, star)
        QSn = eik.parent()
        mul = QSn.right_action_product
        epsilon_ik_cache[(it,kt)] = mul(mul(epsilon(it, star+1), eik),
                                        epsilon(kt, star+1)) * (1/kappa(it.shape()))
        res =  epsilon_ik_cache[(it,kt)]

    return res


epsilon_cache = {}
def epsilon(tab, star=0):
    r"""
    The `(T, T)`-th element of the seminormal basis of the group
    algebra `\QQ[S_n]`, where `T` is the tableau ``tab`` (with its
    ``star`` highest entries removed if the optional variable
    ``star`` is set).

    See the docstring of
    :meth:`~sage.combinat.symmetric_group_algebra.SymmetricGroupAlgebra_n.seminormal_basis`
    for the notation used herein.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import epsilon
        sage: epsilon([[1,2]])
        1/2*[1, 2] + 1/2*[2, 1]
        sage: epsilon([[1],[2]])
        1/2*[1, 2] - 1/2*[2, 1]
    """
    t = Tableau(tab)

    if star:
        t = t.restrict(t.size() - star)

    if t in epsilon_cache:
        res = epsilon_cache[t]
    else:
        if t.size() == 2:
            epsilon_cache[t] = e(t)*(1 / kappa(t.shape()))
            res =  epsilon_cache[t]
        elif t == Tableau([[1]]):
            epsilon_cache[t] = e(t)
            res =  epsilon_cache[t]
        else:
            et = e(t)
            QSn = et.parent()
            mul = QSn.right_action_product
            epsilon_cache[t] = mul(mul(epsilon(t, 1), e(t)), epsilon(t, 1)) * (1 / kappa(t.shape()))
            res = epsilon_cache[t]

    return res


def pi_ik(itab, ktab):
    r"""
    Return the permutation `p` which sends every entry of the
    tableau ``itab`` to the respective entry of the tableau
    ``ktab``, as an element of the corresponding symmetric group
    algebra.

    This assumes that ``itab`` and ``ktab`` are tableaux (possibly
    given just as lists of lists) of the same shape.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import pi_ik
        sage: pi_ik([[1,3],[2]], [[1,2],[3]])
        [1, 3, 2]
    """
    it = Tableau(itab)
    kt = Tableau(ktab)

    p = [None]*kt.size()
    for i in range(len(kt)):
        for j in range(len(kt[i])):
            p[ it[i][j] -1 ] = kt[i][j]

    QSn = SymmetricGroupAlgebra(QQ, it.size())
    p = Permutation(p)
    return QSn(p)


def kappa(alpha):
    r"""
    Return `\kappa_\alpha`, which is `n!` divided by the number
    of standard tableaux of shape `\alpha` (where `\alpha` is a
    partition of `n`).

    INPUT:

    - ``alpha`` -- integer partition (can be encoded as a list).

    OUTPUT:

    The factorial of the size of ``alpha``, divided by the number of
    standard tableaux of shape ``alpha``. Equivalently, the product
    of all hook lengths of ``alpha``.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import kappa
        sage: kappa(Partition([2,1]))
        3
        sage: kappa([2,1])
        3
    """
    try:
        n = alpha.size()
    except AttributeError:
        n = sum(alpha)
    return factorial(n) / StandardTableaux(alpha).cardinality()

def a(tableau, star=0, base_ring=QQ):
    r"""
    The row projection operator corresponding to the Young tableau
    ``tableau`` (which is supposed to contain every integer from
    `1` to its size precisely once, but may and may not be standard).

    This is the sum (in the group algebra of the relevant symmetric
    group over `\QQ`) of all the permutations which preserve
    the rows of ``tableau``. It is called `a_{\text{tableau}}` in
    [EtRT]_, Section 4.2.

    REFERENCES:

    .. [EtRT] Pavel Etingof, Oleg Golberg, Sebastian Hensel, Tiankai
       Liu, Alex Schwendner, Dmitry Vaintrob, Elena Yudovina,
       "Introduction to representation theory",
       :arXiv:`0901.0827v5`.

    INPUT:

    - ``tableau`` -- Young tableau which contains every integer
      from `1` to its size precisely once.

    - ``star`` -- nonnegative integer (default: `0`). When this
      optional variable is set, the method computes not the row
      projection operator of ``tableau``, but the row projection
      operator of the restriction of ``tableau`` to the entries
      ``1, 2, ..., tableau.size() - star`` instead.

    - ``base_ring`` -- commutative ring (default: ``QQ``). When this
      optional variable is set, the row projection operator is
      computed over a user-determined base ring instead of `\QQ`.
      (Note that symmetric group algebras currently don't preserve
      coercion, so e. g. a symmetric group algebra over `\ZZ`
      does not coerce into the corresponding one over `\QQ`; so
      convert manually or choose your base rings wisely!)

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import a
        sage: a([[1,2]])
        [1, 2] + [2, 1]
        sage: a([[1],[2]])
        [1, 2]
        sage: a([])
        []
        sage: a([[1, 5], [2, 3], [4]])
        [1, 2, 3, 4, 5] + [1, 3, 2, 4, 5] + [5, 2, 3, 4, 1] + [5, 3, 2, 4, 1]
        sage: a([[1,4], [2,3]], base_ring=ZZ)
        [1, 2, 3, 4] + [1, 3, 2, 4] + [4, 2, 3, 1] + [4, 3, 2, 1]
    """
    t = Tableau(tableau)
    if star:
        t = t.restrict(t.size()-star)

    rs = t.row_stabilizer().list()
    n = t.size()

    sgalg = SymmetricGroupAlgebra(base_ring, n)
    one = base_ring.one()
    P = Permutation

    # Ugly hack for the case of an empty tableau, due to the
    # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
    # being [1] rather than [] (which seems to have its origins in
    # permutation group code).
    # TODO: Fix this.
    if len(tableau) == 0:
        return sgalg.one()

    rd = dict((P(h), one) for h in rs)
    return sgalg._from_dict(rd)

def b(tableau, star=0, base_ring=QQ):
    r"""
    The column projection operator corresponding to the Young tableau
    ``tableau`` (which is supposed to contain every integer from
    `1` to its size precisely once, but may and may not be standard).

    This is the signed sum (in the group algebra of the relevant
    symmetric group over `\QQ`) of all the permutations which
    preserve the column of ``tableau`` (where the signs are the usual
    signs of the permutations). It is called `b_{\text{tableau}}` in
    [EtRT]_, Section 4.2.

    INPUT:

    - ``tableau`` -- Young tableau which contains every integer
      from `1` to its size precisely once.

    - ``star`` -- nonnegative integer (default: `0`). When this
      optional variable is set, the method computes not the column
      projection operator of ``tableau``, but the column projection
      operator of the restriction of ``tableau`` to the entries
      ``1, 2, ..., tableau.size() - star`` instead.

    - ``base_ring`` -- commutative ring (default: ``QQ``). When this
      optional variable is set, the column projection operator is
      computed over a user-determined base ring instead of `\QQ`.
      (Note that symmetric group algebras currently don't preserve
      coercion, so e. g. a symmetric group algebra over `\ZZ`
      does not coerce into the corresponding one over `\QQ`; so
      convert manually or choose your base rings wisely!)

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import b
        sage: b([[1,2]])
        [1, 2]
        sage: b([[1],[2]])
        [1, 2] - [2, 1]
        sage: b([])
        []
        sage: b([[1, 2, 4], [5, 3]])
        [1, 2, 3, 4, 5] - [1, 3, 2, 4, 5] - [5, 2, 3, 4, 1] + [5, 3, 2, 4, 1]
        sage: b([[1, 4], [2, 3]], base_ring=ZZ)
        [1, 2, 3, 4] - [1, 2, 4, 3] - [2, 1, 3, 4] + [2, 1, 4, 3]
        sage: b([[1, 4], [2, 3]], base_ring=Integers(5))
        [1, 2, 3, 4] + 4*[1, 2, 4, 3] + 4*[2, 1, 3, 4] + [2, 1, 4, 3]

    With the ``l2r`` setting for multiplication, the unnormalized
    Young symmetrizer ``e(tableau)`` should be the product
    ``b(tableau) * a(tableau)`` for every ``tableau``. Let us check
    this on the standard tableaux of size 5::

        sage: from sage.combinat.symmetric_group_algebra import a, b, e
        sage: all( e(t) == b(t) * a(t) for t in StandardTableaux(5) )
        True
    """
    t = Tableau(tableau)
    if star:
        t = t.restrict(t.size()-star)

    cs = t.column_stabilizer().list()
    n = t.size()

    sgalg = SymmetricGroupAlgebra(base_ring, n)
    one = base_ring.one()
    P = Permutation

    # Ugly hack for the case of an empty tableau, due to the
    # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
    # being [1] rather than [] (which seems to have its origins in
    # permutation group code).
    # TODO: Fix this.
    if len(tableau) == 0:
        return sgalg.one()

    cd = dict((P(v), v.sign()*one) for v in cs)
    return sgalg._from_dict(cd)

e_cache = {}
def e(tableau, star=0):
    r"""
    The unnormalized Young projection operator corresponding to
    the Young tableau ``tableau`` (which is supposed to contain
    every integer from `1` to its size precisely once, but may
    and may not be standard).

    If `n` is a nonnegative integer, and `T` is a Young tableau
    containing every integer from `1` to `n` exactly once, then
    the unnormalized Young projection operator `e(T)` is defined by

    .. MATH::

        e(T) = a(T) b(T) \in \QQ S_n,

    where `a(T) \in \QQ S_n` is the sum of all permutations in `S_n`
    which fix the rows of `T` (as sets), and `b(T) \in \QQ S_n` is the
    signed sum of all permutations in `S_n` which fix the columns of
    `T` (as sets). Here, "signed" means that each permutation is
    multiplied with its sign; and the product on the group `S_n` is
    defined in such a way that `(pq)(i) = p(q(i))` for any
    permutations `p` and `q` and any `1 \leq i \leq n`.

    Note that the definition of `e(T)` is not uniform across
    literature. Others define it as `b(T) a(T)` instead, or include
    certain scalar factors (we do not, whence "unnormalized").

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import e
        sage: e([[1,2]])
        [1, 2] + [2, 1]
        sage: e([[1],[2]])
        [1, 2] - [2, 1]
        sage: e([])
        []

    There are differing conventions for the order of the symmetrizers
    and antisymmetrizers.  This example illustrates our conventions::

        sage: e([[1,2],[3]])
        [1, 2, 3] + [2, 1, 3] - [3, 1, 2] - [3, 2, 1]

    To obtain the product `b(T) a(T)`, one has to take the antipode
    of this::

        sage: QS3 = parent(e([[1,2],[3]]))
        sage: QS3.antipode(e([[1,2],[3]]))
        [1, 2, 3] + [2, 1, 3] - [2, 3, 1] - [3, 2, 1]

    .. SEEALSO::

        :func:`e_hat`
    """
    # TODO:
    # The current method only computes the e's over QQ. There should be
    # a way to compute them over other base rings as well. Be careful
    # with the cache.

    t = Tableau(tableau)
    if star:
        t = t.restrict(t.size()-star)

    if t in e_cache:
        res = e_cache[t]
    else:
        rs = t.row_stabilizer().list()
        cs = t.column_stabilizer().list()
        n = t.size()

        QSn = SymmetricGroupAlgebra(QQ, n)
        one = QQ.one()
        P = Permutation

        rd = dict((P(h), one) for h in rs)
        sym = QSn._from_dict(rd)

        cd = dict((P(v), v.sign()*one) for v in cs)
        antisym = QSn._from_dict(cd)

        res = QSn.right_action_product(antisym, sym)

        # Ugly hack for the case of an empty tableau, due to the
        # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
        # being [1] rather than [] (which seems to have its origins in
        # permutation group code).
        # TODO: Fix this.
        if len(tableau) == 0:
            res = QSn.one()

        e_cache[t] = res

    return res

ehat_cache = {}
def e_hat(tab, star=0):
    r"""
    The Young projection operator corresponding to the Young tableau
    ``tab`` (which is supposed to contain every integer from `1` to
    its size precisely once, but may and may not be standard). This
    is an idempotent in the rational group algebra.

    If `n` is a nonnegative integer, and `T` is a Young tableau
    containing every integer from `1` to `n` exactly once, then
    the Young projection operator `\widehat{e}(T)` is defined by

    .. MATH::

        \widehat{e}(T) = \frac{1}{\kappa_\lambda} a(T) b(T) \in \QQ S_n,

    where `\lambda` is the shape of `T`, where `\kappa_\lambda` is
    `n!` divided by the number of standard tableaux of shape
    `\lambda`, where `a(T) \in \QQ S_n` is the sum of all
    permutations in `S_n` which fix the rows of `T` (as sets), and
    where `b(T) \in \QQ S_n` is the signed sum of all permutations
    in `S_n` which fix the columns of `T` (as sets). Here, "signed"
    means that each permutation is multiplied with its sign; and
    the product on the group `S_n` is defined in such a way that
    `(pq)(i) = p(q(i))` for any permutations `p` and `q` and any
    `1 \leq i \leq n`.

    Note that the definition of `\widehat{e}(T)` is not uniform
    across literature. Others define it as
    `\frac{1}{\kappa_\lambda} b(T) a(T)` instead.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import e_hat
        sage: e_hat([[1,2,3]])
        1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
        sage: e_hat([[1],[2]])
        1/2*[1, 2] - 1/2*[2, 1]

    There are differing conventions for the order of the symmetrizers
    and antisymmetrizers.  This example illustrates our conventions::

        sage: e_hat([[1,2],[3]])
        1/3*[1, 2, 3] + 1/3*[2, 1, 3] - 1/3*[3, 1, 2] - 1/3*[3, 2, 1]

    .. SEEALSO::

        :func:`e`
    """
    t = Tableau(tab)
    if star:
        t = t.restrict(t.size()-star)
    if t in ehat_cache:
        res = ehat_cache[t]
    else:
        res = (1/kappa(t.shape()))*e(t)
    return res

e_ik_cache = {}
def e_ik(itab, ktab, star=0):
    """
    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import e_ik
        sage: e_ik([[1,2,3]], [[1,2,3]])
        [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
        sage: e_ik([[1,2,3]], [[1,2,3]], star=1)
        [1, 2] + [2, 1]
    """
    it = Tableau(itab)
    kt = Tableau(ktab)
    if star:
        it = it.restrict(it.size() - star)
        kt = kt.restrict(kt.size() - star)

    if it.shape() != kt.shape():
        raise ValueError("the two tableaux must be of the same shape")

    if kt == it:
        return e(it)
    if (it, kt) in e_ik_cache:
        return e_ik_cache[(it,kt)]

    pi = pi_ik(it,kt)
    QSn = pi.parent()
    res = QSn.right_action_product(e(it), pi)
    e_ik_cache[(it,kt)] = res
    return res

def seminormal_test(n):
    """
    Run a variety of tests to verify that the construction of the
    seminormal basis works as desired. The numbers appearing are
    results in James and Kerber's 'Representation Theory of the
    Symmetric Group' [JamesKerber]_.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import seminormal_test
        sage: seminormal_test(3)
        True
    """
    for part in partition.Partitions_n(n):
        for tab in StandardTableaux(part):
            #Theorem 3.1.10
            if not e(tab)*(1/kappa(part)) - e_hat(tab) == 0:
                raise ValueError("3.1.10 - %s"%tab)

            #Lemma 3.2.12 (ii)
            value = e(tab)*epsilon(tab,1)*e(tab) - e(tab)*(kappa(part))
            if value != 0:
                print value
                raise ValueError("3.2.12.2 - %s"%tab)

            for tab2 in StandardTableaux(part):
                #3.2.8 (i)
                if e_ik(tab, tab2) - e(tab)*pi_ik(tab, tab2)*e(tab2)*(1/kappa(part)) != 0:
                    raise ValueError("3.2.8.1 - %s, %s"%(tab, tab2))

                #3.2.8 (ii)
                if e(tab)*e_ik(tab, tab2) - e_ik(tab, tab2)*(kappa(part)) != 0:
                    raise ValueError("3.2.8.2 - %s, %s"%(tab, tab2))

                if tab == tab2:
                    continue

                if tab.last_letter_lequal(tab2):
                    #Lemma 3.1.20
                    if e(tab2)*e(tab) != 0:
                        raise ValueError("3.1.20 - %s, %s"%(tab, tab2))
                    if e_hat(tab2)*e_hat(tab) != 0:
                        raise ValueError("3.1.20 - %s, %s"%(tab, tab2))
    return True

#######################


def HeckeAlgebraSymmetricGroupT(R, n, q=None):
    r"""
    Return the Hecke algebra of the symmetric group `S_n` on the T-basis
    with quantum parameter ``q`` over the ring `R`.

    If `R` is a commutative ring and `q` is an invertible element of `R`,
    and if `n` is a nonnegative integer, then the Hecke algebra of the
    symmetric group `S_n` over `R` with quantum parameter `q` is defined
    as the algebra generated by the generators `T_1, T_2, \ldots, T_{n-1}`
    with relations

    .. MATH::

        T_i T_{i+1} T_i = T_{i+1} T_i T_{i+1}

    for all `i < n-1` ("braid relations"),

    .. MATH::

        T_i T_j = T_j T_i

    for all `i` and `j` such that `| i-j | > 1` ("locality relations"),
    and

    .. MATH::

        T_i^2 = q + (q-1) T_i

    for all `i` (the "quadratic relations", also known in the form
    `(T_i + 1) (T_i - q) = 0`). (This is only one of several existing
    definitions in literature, not all of which are fully equivalent.
    We are following the conventions of [GS93]_.) For any permutation
    `w \in S_n`, we can define an element `T_w` of this Hecke algebra by
    setting `T_w = T_{i_1} T_{i_2} \cdots T_{i_k}`, where
    `w = s_{i_1} s_{i_2} \cdots s_{i_k}` is a reduced word for `w`
    (with `s_i` meaning the transposition `(i, i+1)`, and the product of
    permutations being evaluated by first applying `s_{i_k}`, then
    `s_{i_{k-1}}`, etc.). This element is independent of the choice of
    the reduced decomposition, and can be computed in Sage by calling
    ``H[w]`` where ``H`` is the Hecke algebra and ``w`` is the
    permutation.

    The Hecke algebra of the symmetric group `S_n` with quantum parameter
    `q` over `R` can be seen as a deformation of the group algebra
    `R S_n`; indeed, it becomes `R S_n` when `q = 1`.

    .. WARNING::

        The multiplication on the Hecke algebra of the symmetric group
        does *not* follow the global option ``mult`` of the
        :class:`Permutations` class (see
        :meth:`~sage.combinat.permutation.Permutations.global_options`).
        It is always as defined above. It does not match the default
        option (``mult=l2r``) of the symmetric group algebra!

    REFERENCES:

    .. [GS93] David M. Goldschmidt.
       *Group characters, symmetric functions, and the Hecke algebras*.
       AMS 1993.

    EXAMPLES::

        sage: HeckeAlgebraSymmetricGroupT(QQ, 3)
        Hecke algebra of the symmetric group of order 3 on the T basis over Univariate Polynomial Ring in q over Rational Field

    ::

        sage: HeckeAlgebraSymmetricGroupT(QQ, 3, 2)
        Hecke algebra of the symmetric group of order 3 with q=2 on the T basis over Rational Field

    The multiplication on the Hecke algebra follows a different convention
    than the one on the symmetric group algebra does by default::

        sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
        sage: H3([1,3,2]) * H3([2,1,3])
        T[3, 1, 2]
        sage: S3 = SymmetricGroupAlgebra(QQ, 3)
        sage: S3([1,3,2]) * S3([2,1,3])
        [2, 3, 1]

        sage: TestSuite(H3).run()
    """

    return HeckeAlgebraSymmetricGroup_t(R, n, q)

class HeckeAlgebraSymmetricGroup_generic(CombinatorialAlgebra):
    def __init__(self, R, n, q=None):
        """
        TESTS::

            sage: HeckeAlgebraSymmetricGroupT(QQ, 3)
            Hecke algebra of the symmetric group of order 3 on the T basis over Univariate Polynomial Ring in q over Rational Field

        ::

            sage: HeckeAlgebraSymmetricGroupT(QQ, 3, q=1)
            Hecke algebra of the symmetric group of order 3 with q=1 on the T basis over Rational Field
        """
        self.n = n
        self._basis_keys = Permutations(n)
        self._name = "Hecke algebra of the symmetric group of order {}".format(n)
        self._one = self._basis_keys(range(1,n+1))

        if q is None:
            q = PolynomialRing(R, 'q').gen()
            R = q.parent()
        else:
            if q not in R:
                raise ValueError("q must be in R (= {})".format(R))
            self._name += " with q={}".format(q)

        self._q = q

        CombinatorialAlgebra.__init__(self, R)
        # _repr_ customization: output the basis element indexed by [1,2,3] as [1,2,3]
        self.print_options(prefix="")

    _repr_option_bracket = False

    def q(self):
        """
        EXAMPLES::

            sage: HeckeAlgebraSymmetricGroupT(QQ, 3).q()
            q
            sage: HeckeAlgebraSymmetricGroupT(QQ, 3, 2).q()
            2
        """
        return self._q


    def _coerce_start(self, x):
        """
        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: H3._coerce_start([2,1])
            T[2, 1, 3]
        """
        ###################################################
        # Coerce permutations of size smaller that self.n #
        ###################################################
        if x == []:
            return self.one()
        if len(x) < self.n and x in Permutations():
            return self.monomial(self._basis_keys( list(x) + range(len(x)+1, self.n+1) ))
        raise TypeError

class HeckeAlgebraSymmetricGroup_t(HeckeAlgebraSymmetricGroup_generic):

    def __init__(self, R, n, q=None):
        """
        TESTS::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: H3 == loads(dumps(H3))
            True
        """
        HeckeAlgebraSymmetricGroup_generic.__init__(self, R, n, q)
        self._name += " on the T basis"
        self.print_options(prefix="T")

    def t_action_on_basis(self, perm, i):
        r"""
        Return the product `T_i \cdot T_{perm}`, where ``perm`` is a
        permutation in the symmetric group `S_n`.

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: H3.t_action_on_basis(Permutation([2,1,3]), 1)
            q*T[1, 2, 3] + (q-1)*T[2, 1, 3]
            sage: H3.t_action_on_basis(Permutation([1,2,3]), 1)
            T[2, 1, 3]
            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3, 1)
            sage: H3.t_action_on_basis(Permutation([2,1,3]), 1)
            T[1, 2, 3]
            sage: H3.t_action_on_basis(Permutation([1,3,2]), 2)
            T[1, 2, 3]
        """
        if i not in range(1, self.n):
            raise ValueError("i (= %(i)d) must be between 1 and n (= %(n)d)" % {'i': i, 'n': self.n})

        t_i = Permutation( (i, i+1) )
        perm_i = t_i.right_action_product(perm)
        # This used to be perm_i = t_i * perm. I have changed it to
        # perm_i = t_i.right_action_product(perm) because it would
        # otherwise cause TestSuite(H3) to fail when
        # Permutations.global_options(mult) would be set to "r2l".
        # -- Darij, 19 Nov 2013

        if perm[i-1] < perm[i]:
            return self.monomial(self._basis_keys(perm_i))
        else:
            #Ti^2 = (q - q^(-1))*Ti - q1*q2
            q = self.q()
            z_elt = {perm_i:q, perm:q-1}
            return self._from_dict(z_elt)


    def t_action(self, a, i):
        r"""
        Return the product `T_i \cdot a`.

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: a = H3([2,1,3])+2*H3([1,2,3])
            sage: H3.t_action(a, 1)
            q*T[1, 2, 3] + (q+1)*T[2, 1, 3]
            sage: H3.t(1)*a
            q*T[1, 2, 3] + (q+1)*T[2, 1, 3]
        """
        t_i = lambda x: self.t_action_on_basis(x, i)
        return self._apply_module_endomorphism(a, t_i)


    def _multiply_basis(self, perm1, perm2):
        """
        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3, 1)
            sage: a = H3([2,1,3])+2*H3([1,2,3])-H3([3,2,1])
            sage: a^2 #indirect doctest
            6*T[1, 2, 3] + 4*T[2, 1, 3] - T[2, 3, 1] - T[3, 1, 2] - 4*T[3, 2, 1]

        ::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = QS3([2,1,3])+2*QS3([1,2,3])-QS3([3,2,1])
            sage: a^2
            6*[1, 2, 3] + 4*[2, 1, 3] - [2, 3, 1] - [3, 1, 2] - 4*[3, 2, 1]
        """
        res = self(perm1)
        for i in perm2.reduced_word():
            res = self.t_action(res, i)
        return res

    def t(self, i):
        """
        Return the element `T_i` of the Hecke algebra ``self``.

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ,3)
            sage: H3.t(1)
            T[2, 1, 3]
            sage: H3.t(2)
            T[1, 3, 2]
            sage: H3.t(0)
            Traceback (most recent call last):
            ...
            ValueError: i (= 0) must be between 1 and n-1 (= 2)
        """
        if i not in range(1, self.n):
            raise ValueError("i (= %(i)d) must be between 1 and n-1 (= %(nm)d)" % {'i': i, 'nm': self.n - 1})

        P = self.basis().keys()
        return self.monomial(P( range(1, i) + [i+1, i] + range(i+2, self.n+1) ))
        # The permutation here is simply the transposition (i, i+1).

    def algebra_generators(self):
        """
        Return the generators of the algebra.

        EXAMPLES::

            sage: HeckeAlgebraSymmetricGroupT(QQ,3).algebra_generators()
            [T[2, 1, 3], T[1, 3, 2]]
        """
        return map(self.t, range(1, self.n))

    def jucys_murphy(self, k):
        """
        Return the Jucys-Murphy element `J_k` of the Hecke algebra.

        These Jucys-Murphy elements are defined by

        .. MATH::

            J_k = (T_{k-1} T_{k-2} \cdots T_1) (T_1 T_2 \cdots T_{k-1}).

        More explicitly,

        .. MATH::

            J_k = q^{k-1} + \sum_{l=1}^{k-1} (q^l - q^{l-1}) T_{(l, k)}.

        For generic `q`, the `J_k` generate a maximal commutative
        sub-algebra of the Hecke algebra.

        .. WARNING::

            The specialization `q = 1` does *not* map these elements
            `J_k` to the Young-Jucys-Murphy elements of the group
            algebra `R S_n`. (Instead, it maps the "reduced"
            Jucys-Murphy elements `(J_k - q^{k-1}) / (q - 1)` to the
            Young-Jucys-Murphy elements of `R S_n`.)

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ,3)
            sage: j2 = H3.jucys_murphy(2); j2
            q*T[1, 2, 3] + (q-1)*T[2, 1, 3]
            sage: j3 = H3.jucys_murphy(3); j3
            q^2*T[1, 2, 3] + (q^2-q)*T[1, 3, 2] + (q-1)*T[3, 2, 1]
            sage: j2*j3 == j3*j2
            True
            sage: j0 = H3.jucys_murphy(1); j0 == H3.one()
            True
            sage: H3.jucys_murphy(0)
            Traceback (most recent call last):
            ...
            ValueError: k (= 0) must be between 1 and n (= 3)
        """
        if k not in range(2, self.n+1):
            if k == 1:
                return self.one()
            raise ValueError("k (= %(k)d) must be between 1 and n (= %(n)d)" % {'k': k, 'n': self.n})

        q = self.q()
        P = self._basis_keys
        v = self.sum_of_terms( ( ( P(range(1, l) + [k] + range(l+1, k) + [l]),
                                   q ** l - q ** (l-1) )
                                 for l in range(1, k) ),
                               distinct=True )
        v += q ** (k-1) * self.one()
        return v
        
        #old algorithm:
        # left = 1
        # right = 1
        # for j in range(1, k):
        #    left *= self.t(k-j)
        #    right *= self.t(j)
        # return left*right


# For unpickling backward compatibility (Sage <= 4.1)
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.symmetric_group_algebra', 'HeckeAlgebraSymmetricGroupElement_t',  CombinatorialFreeModule.Element)
register_unpickle_override('sage.combinat.symmetric_group_algebra', 'SymmetricGroupAlgebraElement_n',  CombinatorialFreeModule.Element)

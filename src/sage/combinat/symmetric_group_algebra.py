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
import permutation
import partition
from tableau import Tableau, StandardTableaux_size, StandardTableaux_shape, StandardTableaux
from sage.interfaces.all import gap
from sage.rings.all import factorial, QQ, PolynomialRing
from sage.matrix.all import matrix
from sage.modules.all import vector
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.categories.all import GroupAlgebras

permutation_options = permutation.PermutationOptions

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
        CombinatorialFreeModule.__init__(self, R, permutation.Permutations(n), prefix='', latex_prefix='', category = (GroupAlgebras(R),FiniteDimensionalAlgebrasWithBasis(R)))
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
        return permutation.Permutation(range(1,self.n+1))

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
        return self.monomial( self.basis().keys()(permutation) )


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

        cpi = self.zero()

        np = partition.Partitions_n(self.n).list()
        np.reverse()
        p_index = np.index(p)

        big_coeff = character_table[p_index][0] / factorial(self.n)

        character_row = character_table[p_index]
        dct = { g : big_coeff * character_row[np.index(g.cycle_type())]
                for g in permutation.StandardPermutations_n(self.n) }

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
        return [permutation.Permutations(self.n).element_in_conjugacy_classes(nu) for nu in partition.Partitions(self.n)]

    def jucys_murphy(self, k):
        """
        Return the Jucys-Murphy element `J_k` (also known as a
        Young-Jucys-Murphy element) for the symmetric group
        algebra ``self``.

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
        res = self.zero()

        if k < 1 or k > self.n:
            raise ValueError("k (= {k}) must be between 1 and n (= {n}) (inclusive)".format(k=k, n=self.n))

        for i in range(1, k):
            p = range(1, self.n+1)
            p[i-1] = k
            p[k-1] = i
            res += self.monomial(self._basis_keys(p))
        return res



    def seminormal_basis(self):
        """
        Return a list of the seminormal basis elements of ``self``.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3.seminormal_basis()
            [1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1],
            1/3*[1, 2, 3] + 1/6*[1, 3, 2] - 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] + 1/6*[3, 2, 1],
            1/3*[1, 3, 2] + 1/3*[2, 3, 1] - 1/3*[3, 1, 2] - 1/3*[3, 2, 1],
            1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1],
            1/3*[1, 2, 3] - 1/6*[1, 3, 2] + 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] - 1/6*[3, 2, 1],
            1/6*[1, 2, 3] - 1/6*[1, 3, 2] - 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] - 1/6*[3, 2, 1]]
        """
        basis = []
        for part in partition.Partitions_n(self.n):
            stp = StandardTableaux_shape(part)
            for t1 in stp:
                for t2 in stp:
                    basis.append(self.epsilon_ik(t1,t2))
        return basis


    def dft(self, form="seminormal"):
        """
        Return the discrete Fourier transform for ``self``.

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
            return self._dft_seminormal()
        else:
            raise ValueError("invalid form (= %s)"%form)

    def _dft_seminormal(self):
        """
        Return the seminormal form of the discrete Fourier for ``self``.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3._dft_seminormal()
            [   1    1    1    1    1    1]
            [   1  1/2   -1 -1/2 -1/2  1/2]
            [   0  3/4    0  3/4 -3/4 -3/4]
            [   0    1    0   -1    1   -1]
            [   1 -1/2    1 -1/2 -1/2 -1/2]
            [   1   -1   -1    1    1   -1]
        """
        snb = self.seminormal_basis()
        return matrix( [vector(b) for b in snb] ).inverse().transpose()






    def epsilon_ik(self, itab, ktab, star=0):
        """
        Return the seminormal basis element of ``self`` corresponding to the
        pair of tableaux ``itab`` and ``ktab``.

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

        if permutation_options['mult'] == 'l2r':
            return z
        else:
            return z.map_support(lambda x: x.inverse())


epsilon_ik_cache = {}
def epsilon_ik(itab, ktab, star=0):
    """
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

    mult = permutation_options['mult']
    permutation_options['mult'] = 'l2r'
    if kt == it:
        res = epsilon(itab)
    elif (it, kt) in epsilon_ik_cache:
        res =  epsilon_ik_cache[(it,kt)]
    else:
        epsilon_ik_cache[(it,kt)] = epsilon(it, star+1)*e_ik(it,kt,star)*epsilon(kt, star+1) * (1/kappa(it.shape()))
        res =  epsilon_ik_cache[(it,kt)]

    permutation_options['mult'] = mult
    return res


epsilon_cache = {}
def epsilon(tab, star=0):
    """
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

    mult = permutation_options['mult']
    permutation_options['mult'] = 'l2r'

    if t in epsilon_cache:
        res = epsilon_cache[t]
    else:
        if t.size() == 2:
            epsilon_cache[t] = e(t)*(1/kappa(t.shape()))
            res =  epsilon_cache[t]
        elif t == Tableau([[1]]):
            epsilon_cache[t] = e(t)
            res =  epsilon_cache[t]
        else:
            epsilon_cache[t] =  epsilon(t, 1)*e(t)*epsilon(t,1)*( 1 / kappa(t.shape()))
            res = epsilon_cache[t]

    permutation_options['mult'] = mult
    return res


def pi_ik(itab, ktab):
    """
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
    p = permutation.Permutation(p)
    return QSn(p)


def kappa(alpha):
    r"""
    Return `\kappa_\alpha`, which is `n!` divided by the number
    of standard tableaux of shape `\alpha`.

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


def a(tableau, star=0):
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
    """
    t = Tableau(tableau)
    if star:
        t = t.restrict(t.size()-star)

    rs = t.row_stabilizer().list()
    n = t.size()

    # This all should be over ZZ, not over QQ, but symmetric group
    # algebras don't seem to preserve coercion (the one over ZZ
    # doesn't coerce into the one over QQ even for the same n),
    # and the QQ version of this method is more important, so let
    # me stay with QQ.
    # TODO: Fix this.
    sgalg = SymmetricGroupAlgebra(QQ, n)
    one = QQ.one()
    P = permutation.Permutation

    # Ugly hack for the case of an empty tableau, due to the
    # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
    # being [1] rather than [] (which seems to have its origins in
    # permutation group code).
    # TODO: Fix this.
    if len(tableau) == 0:
        return sgalg.one()

    rd = dict((P(h), one) for h in rs)
    return sgalg._from_dict(rd)

def b(tableau, star=0):
    r"""
    The column projection operator corresponding to the Young tableau
    ``tableau`` (which is supposed to contain every integer from
    `1` to its size precisely once, but may and may not be standard).

    This is the signed sum (in the group algebra of the relevant
    symmetric group over `\QQ`) of all the permutations which
    preserve the column of ``tableau`` (where the signs are the usual
    signs of the permutations). It is called `b_{\text{tableau}}` in
    [EtRT]_, Section 4.2.

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

    With the `l2r` setting for multiplication, the unnormalized
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

    # This all should be over ZZ, not over QQ, but symmetric group
    # algebras don't seem to preserve coercion (the one over ZZ
    # doesn't coerce into the one over QQ even for the same n),
    # and the QQ version of this method is more important, so let
    # me stay with QQ.
    # TODO: Fix this.
    sgalg = SymmetricGroupAlgebra(QQ, n)
    one = QQ.one()
    P = permutation.Permutation

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
    """
    The unnormalized Young projection operator corresponding to
    the Young tableau ``tableau`` (which is supposed to contain
    every integer from `1` to its size precisely once, but may
    and may not be standard).

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
    """
    t = Tableau(tableau)
    if star:
        t = t.restrict(t.size()-star)

    mult = permutation_options['mult']
    permutation_options['mult'] = 'l2r'

    if t in e_cache:
        res = e_cache[t]
    else:
        rs = t.row_stabilizer().list()
        cs = t.column_stabilizer().list()
        n = t.size()

        QSn = SymmetricGroupAlgebra(QQ, n)
        one = QQ.one()
        P = permutation.Permutation

        rd = dict((P(h), one) for h in rs)
        sym = QSn._from_dict(rd)

        cd = dict((P(v), v.sign()*one) for v in cs)
        antisym = QSn._from_dict(cd)

        res = antisym*sym

        # Ugly hack for the case of an empty tableau, due to the
        # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
        # being [1] rather than [] (which seems to have its origins in
        # permutation group code).
        # TODO: Fix this.
        if len(tableau) == 0:
            res = QSn.one()

        e_cache[t] = res

    permutation_options['mult'] = mult

    return res

ehat_cache = {}
def e_hat(tab, star=0):
    """
    The Young projection operator corresponding to the Young tableau
    ``tab`` (which is supposed to contain every integer from `1` to
    its size precisely once, but may and may not be standard). This
    is an idempotent in the rational group algebra.

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

    mult = permutation_options['mult']
    permutation_options['mult'] = 'l2r'

    if kt == it:
        res =  e(it)
    elif (it, kt) in e_ik_cache:
        res = e_ik_cache[(it,kt)]
    else:
        pi = pi_ik(it,kt)
        e_ik_cache[(it,kt)] = e(it)*pi
        res = e_ik_cache[(it,kt)]

    permutation_options['mult'] = mult
    return res


def seminormal_test(n):
    """
    Runs a variety of tests to verify that the construction of the
    seminormal basis works as desired. The numbers appearing are
    Theorems in James and Kerber's 'Representation Theory of the
    Symmetric Group'.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import seminormal_test
        sage: seminormal_test(3)
        True
    """
    for part in partition.Partitions_n(n):
        for tab in StandardTableaux(part):
            #3.1.10
            if not e(tab)*(1/kappa(part)) - e_hat(tab) == 0:
                raise ValueError("3.1.10 - %s"%tab)

            #3.2.12.2
            value = e(tab)*epsilon(tab,1)*e(tab) - e(tab)*(kappa(part))
            if value != 0:
                print value
                raise ValueError("3.2.12.2 - %s"%tab)

            for tab2 in StandardTableaux(part):
                #3.2.8 1
                if e_ik(tab, tab2) - e(tab)*pi_ik(tab, tab2)*e(tab2)*(1/kappa(part)) != 0:
                    raise ValueError("3.2.8.1 - %s, %s"%(tab, tab2))

                #3.2.8.1
                if e(tab)*e_ik(tab, tab2) - e_ik(tab, tab2)*(kappa(part)) != 0:
                    raise ValueError("3.2.8.2 - %s, %s"%(tab, tab2))

                if tab == tab2:
                    continue

                if tab.last_letter_lequal(tab2):
                    #3.1.20
                    if e(tab2)*e(tab) != 0:
                        raise ValueError("3.1.20 - %s, %s"%(tab, tab2))
                    if e_hat(tab2)*e_hat(tab) != 0:
                        raise ValueError("3.1.20 - %s, %s"%(tab, tab2))
    return True

#######################


def HeckeAlgebraSymmetricGroupT(R, n, q=None):
    """
    Return the Hecke algebra of the symmetric group on the T basis.

    EXAMPLES::

        sage: HeckeAlgebraSymmetricGroupT(QQ, 3)
        Hecke algebra of the symmetric group of order 3 on the T basis over Univariate Polynomial Ring in q over Rational Field

    ::

        sage: HeckeAlgebraSymmetricGroupT(QQ, 3, 2)
        Hecke algebra of the symmetric group of order 3 with q=2 on the T basis over Rational Field
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
        self._basis_keys = permutation.Permutations(n)
        self._name = "Hecke algebra of the symmetric group of order %s"%self.n
        self._one = permutation.Permutation(range(1,n+1))

        if q is None:
            q = PolynomialRing(R, 'q').gen()
            R = q.parent()
        else:
            if q not in R:
                raise ValueError("q must be in R (= %s)"%R)
            self._name += " with q=%s"%q

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
        if len(x) < self.n and x in permutation.Permutations():
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
        """
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

        t_i = permutation.Permutation( (i, i+1) )
        perm_i = t_i * perm

        if perm[i-1] < perm[i]:
            return self.monomial(self._basis_keys(perm_i))
        else:
            #Ti^2 = (q - q^(-1))*Ti - q1*q2
            q = self.q()
            z_elt = {perm_i:q, perm:q-1}
            return self._from_dict(z_elt)


    def t_action(self, a, i):
        """
        Return the action of T_i on a.

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

        return self.monomial(self.basis().keys()(permutation.Permutation( (i, i+1) ) ))

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
        Return the Jucys-Murphy element `J_k` of the Hecke algebra. The
        Jucys-Murphy elements generate the maximal commutative sub-algebra
        of the Hecke algebra.

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

        left = 1
        right = 1
        for j in range(1, k):
            left *= self.t(k-j)
            right *= self.t(j)

        return left*right


# For unpickling backward compatibility (Sage <= 4.1)
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.symmetric_group_algebra', 'HeckeAlgebraSymmetricGroupElement_t',  CombinatorialFreeModule.Element)
register_unpickle_override('sage.combinat.symmetric_group_algebra', 'SymmetricGroupAlgebraElement_n',  CombinatorialFreeModule.Element)

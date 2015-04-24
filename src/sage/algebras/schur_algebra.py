r"""
Schur algebras for `GL_n`

This file implements:

- Schur algebras for `GL_n` over an arbitrary field,

- The canonical action of the Schur algebra on a tensor power of the standard
  representation,

- Using the above to calculate the characters of irreducible `GL_n` modules.

AUTHORS:

- Eric Webster (2010-07-01): implement Schur algebra

- Hugh Thomas (2011-05-08): implement action of Schur algebra and characters
  of irreducible modules

REFERENCES:

.. [GreenPoly] J. Green, Polynomial representations of `GL_n`, Springer Verlag.
"""
#*****************************************************************************
#  Copyright (C) 2010 Eric Webster
#  Copyright (C) 2011 Hugh Thomas (hugh.ross.thomas@gmail.com)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.categories.all import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations
from copy import copy
from sage.categories.rings import Rings
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method
from sage.combinat.sf.sf import SymmetricFunctions
from sage.rings.rational_field import QQ
from sage.combinat.words.word import Word
from sage.combinat.words.words import Words
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.combinat.tableau import SemistandardTableaux
from sage.misc.flatten import flatten
from sage.combinat.partition import Partitions, Partition
from sage.matrix.constructor import Matrix


def _schur_I_nr_representatives(n, r, element, index):
    """
    Internal function called by :func:`schur_representation_indices`.
    """
    if r == 0:
        return index

    if len(element) == r:
        index.append(copy(element))
        return

    if len(element) == 0:
        for i in range(1, n + 1):
            element.append(i)
            _schur_I_nr_representatives(n, r, element, index)
            element.pop()
    else:
        for i in range(element[-1], n + 1):
            element.append(i)
            _schur_I_nr_representatives(n, r, element, index)
            element.pop()

    return index


def schur_representative_indices(n, r):
    r"""
    Return a set which functions as a basis of `S_K(n,r)`.

    More specifically, the basis for `S_K(n,r)` consists of
    equivalence classes of pairs words of length ``r`` on the alphabet
    `1\dots n`, where the equivalence relation is simultaneous
    permutation of the two words.  We can therefore fix a
    representative for each equivalence class in which the entries of
    the first word weakly increase, and the entries of the second word
    whose corresponding values in the first word are equal, also
    weakly increase.

    EXAMPLES::

        sage: from sage.algebras.schur_algebra import schur_representative_indices
        sage: schur_representative_indices(2,2)
        [(word: 11, word: 11), (word: 11, word: 12), (word: 11, word: 22),
        (word: 12, word: 11), (word: 12, word: 12), (word: 12, word: 21),
        (word: 12, word: 22), (word: 22, word: 11), (word: 22, word: 12),
        (word: 22, word: 22)]
    """
    basis = []
    I_nr_repr = _schur_I_nr_representatives(n, r, [], [])
    for e in I_nr_repr:
        j = 0
        k = 0
        I1 = []
        l = len(e)
        while k < l:
            if e[k] != e[j]:
                I2 = []
                if j == 0:
                    I1 = _schur_I_nr_representatives(n, k - j, [], [])
                else:
                    I2 = _schur_I_nr_representatives(n, k - j, [], [])
                    I = []
                    for m1 in range(len(I1)):
                        for m2 in range(len(I2)):
                            I.append(I1[m1] + I2[m2])
                    I1 = I
                j = k
            elif k == l - 1:
                I2 = []
                k += 1
                if j == 0:
                    I1 = _schur_I_nr_representatives(n, k - j, [], [])
                else:
                    I2 = _schur_I_nr_representatives(n, k - j, [], [])
                    I = []
                    for m1 in range(len(I1)):
                        for m2 in range(len(I2)):
                            I.append(I1[m1] + I2[m2])
                    I1 = I
            else:
                k += 1

        for v in I1:
            basis.append((Word(e), Word(v)))

    return basis


def schur_representative_from_index(index):
    """
    Simultaneously reorder a pair of words to obtain the equivalent
    element of the distinguished basis of the Schur algebra.

    .. SEEALSO::

        :func:`schur_representative_indices`

    INPUT:

    - A pair of words from `Words (range(1,n+1),r)`

    OUTPUT:

    - The corresponding pair of words ordered correctly.

    EXAMPLES::

        sage: from sage.algebras.schur_algebra import schur_representative_from_index
        sage: w1 = Word([2,1,2,2])
        sage: w2 = Word([1,3,0,0])
        sage: schur_representative_from_index((w1,w2))
        (word: 1222, word: 3001)
    """
    w = []
    for i in range(len(index[0])):
        w.append((index[0][i], index[1][i]))
    w.sort()
    index = [[], []]
    for i in range(len(w)):
        index[0].append(w[i][0])
        index[1].append(w[i][1])
    return tuple(map(Word, index))


class SchurAlgebra(CombinatorialFreeModule):
    """
    This is the class that implements Schur algebras.

    EXAMPLES::

        sage: from sage.algebras.all import SchurAlgebra
        sage: S = SchurAlgebra(2, 2, ZZ); S
        Schur Algebra (2,2) over Integer Ring

    TESTS::

        sage: SchurAlgebra(-2, 2, ZZ)
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive integer (n=-2)
        sage: SchurAlgebra(2, -2, ZZ)
        Traceback (most recent call last):
        ...
        ValueError: r must be a non-negative integer (r=-2)
        sage: SchurAlgebra(2, 2, 'niet')
        Traceback (most recent call last):
        ...
        ValueError: R must be a commutative Ring (R=niet)
    """
    def __init__(self, n, r, R):
        if not isinstance(n, (int, Integer)) or n <= 0:
            raise ValueError("n must be a positive integer (n=%s)" % (n))
        if not isinstance(r, (int, Integer)) or r < 0:
            raise ValueError("r must be a non-negative integer (r=%s)" % (r))
        if not R in Rings.Commutative():
            raise ValueError("R must be a commutative Ring (R=%s)" % (R))

        self._n = n
        self._r = r

        CombinatorialFreeModule.__init__(self, R,
                                         schur_representative_indices(n, r),
                                         category=AlgebrasWithBasis(R))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.all import SchurAlgebra
            sage: S = SchurAlgebra(4, 4, ZZ)   # indirect doctest
            sage: repr(S)
            'Schur Algebra (4,4) over Integer Ring'
        """
        return "Schur Algebra (%s,%s) over %s" % (self._n, self._r,
                                                  self.base_ring())

    @cached_method
    def one_basis(self):
        """
        Return the index of the one of the algebra.

        THIS IS WRONG !

        EXAMPLES::

            sage: from sage.algebras.all import SchurAlgebra
            sage: SchurAlgebra(4, 4, ZZ).one()   # indirect doctest
            B[None]
        """
        return None

    def product_on_basis(self, e_ij, e_kl):
        r"""
        Return the product of basis elements.

        EXAMPLES::

            sage: from sage.algebras.all import SchurAlgebra
            sage: S = SchurAlgebra(2, 3, QQ)
            sage: B = S.basis()

        If we multiply two basis elements `x` and `y`, such that
        `x[1]` and `y[0]` are not permutations of each other, the
        result is zero ::

            sage: x = B[(Word((1, 1, 1)), Word((1, 1, 2)))]
            sage: y = B[(Word((1, 2, 2)),Word((1, 1, 2)))]
            sage: x * y
            0

        If we multiply a basis element `x` by a basis element which
        consists of the same tuple repeated twice (on either side),
        the result is either zero (if the previous case applies) or `x` ::

            sage: ww = B[(Word((1, 2, 2)), Word((1, 2, 2)))]
            sage: x = B[(Word((1, 2, 2)), Word((1, 1, 2)))]
            sage: ww * x
            B[(word: 122, word: 112)]

        An arbitrary product, on the other hand, may have multiplicities::

            sage: x = B[(Word((1, 1, 1)), Word((1, 1, 2)))]
            sage: y = B[(Word((1, 1, 2)), Word((1, 2, 2)))]
            sage: x * y
            2*B[(word: 111, word: 122)]
        """
        j = e_ij[1]

        i = e_ij[0]
        l = e_kl[1]

        l = sorted(l)

        # Find basis elements (p,q) such that p ~ i and q ~ l
        e_pq = []
        for v in self.basis().keys():
            if v[0] == i and sorted(v[1]) == l:
                e_pq.append(v)

        b = self.basis()
        product = self.zero()

        # Find s in I(n,r) such that (p,s) ~ (i,j) and (s,q) ~ (k,l)
        for e in e_pq:
            Z_ijklpq = self.base_ring()(0)
            for s in Permutations([xx for xx in j]):
                if (schur_representative_from_index((e[0], s)) == e_ij
                        and schur_representative_from_index((s, e[1])) == e_kl):
                    Z_ijklpq += self.base_ring().one()
            product += Z_ijklpq * b[e]

        return product


class TensorSpace(CombinatorialFreeModule):
    """
    This is the ``r``-fold tensor product of an ``n``-dimensional free
    module over ``R``, equipped with an action of the Schur algebra
    `S(n,r)` and the symmetric group `S_r`.

    EXAMPLES::

        sage: from sage.algebras.all import TensorSpace
        sage: TensorSpace(2, 3, QQ)
        The 3-fold tensor product of a free module of dimension 2
        over Rational Field
    """
    def __init__(self, n, r, R):

        self._n = n
        self._r = r
        self._R = R
        CombinatorialFreeModule.__init__(self, R, Words(range(1, n + 1), r))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.all import TensorSpace
            sage: TensorSpace(2, 3, QQ)
            The 3-fold tensor product of a free module of dimension 2
            over Rational Field
        """
        msg = "The {}-fold tensor product of a free module of dimension {}"
        msg += " over {}"
        return msg.format(self._r, self._n, self.base_ring())

    def _basis_elt_from_permuted_word(self, v, perm):
        """
        Return the basis element of ``self`` corresponding to applying
        the permutation perm to a word v.
        """
        return self.basis()[Word(v).apply_permutation_to_positions(perm)]

    def action_by_perm(self, t, perm):
        """
        Apply a permutation to an element `t` of ``self``

        INPUT:

        - ``perm`` -- an element of Permutations(self._r)
        - ``t`` -- an element of ``self``

        OUTPUT:

        - the output is the result of acting by ``perm`` on ``t``
        """
        h = self.module_morphism(self._basis_elt_from_permuted_word,
                                 codomain=self)
        return h(t, perm)

    def action_by_symmetric_group_algebra(self, t, z):
        """
        Return the action by an element of the symmetric group algebra.

        INPUT:

        - ``t`` -- an element of ``self``
        - ``z`` -- an element of ``SymmetricGroupAlgebra(self._R,self._r)``

        OUTPUT:

        result of action of ``z`` on ``t``.
        """
        S = SymmetricGroupAlgebra(self._R, self._r)
        assert z in S
        sym_action = S.module_morphism(self.action_by_perm, codomain=self,
                                       position=1)
        return sym_action(t, z)

    def _monomial_product(self, xi, v):
        """
        Result of acting by the basis element ``xi`` of ``S`` on the
        basis element ``v`` of ``self``.
        """
        x = self.zero()
        for i in Words(range(1, self._n + 1), self._r):
            if schur_representative_from_index((i, v)) == xi:
                x += self.basis()[i]
        return x

    def action_by_Schur_alg(self, nu, v):
        r"""
        Return the action of ``nu`` in Schur algebra on ``v`` in ``self``.
        """
        A = SchurAlgebra(self._n, self._r, self._R)
        assert nu in A
        g = A.module_morphism(self._monomial_product, codomain=self)
        action = self.module_morphism(g, codomain=self, position=1)
        return action(nu, v)


def bracket(r, X, S):
    r"""
    Given ``X`` a set of permutations of ``r`` in cycle notation,
    return the sum in the symmetric group algebra
    of those permutations, times their sign.

    This implements the notation `\{X\}` from just before (5.3a) of Green.

    EXAMPLES::

        sage: P = PermutationGroupElement
        sage: S2 = SymmetricGroupAlgebra(QQ,2)
        sage: sage.algebras.schur_algebra.bracket(2,[P(()),P((1,2))], S2)
        () - (1,2)
    """
    SG = SymmetricGroup(r)
    return sum([x.sign() * S.basis()[SG(x)] for x in X])


def GL_n_irred_character(n, mu, KK):
    r"""
    Return the character of the irreducible module indexed by ``mu``
    of `GL(n)` over the field ``KK``.

    INPUT:

    - ``n`` -- a positive integer.
    - ``mu`` -- a partition of at most ``n`` parts.
    - ``KK`` -- a field.

    OUTPUT:

    a symmetric function which should be interpreted in ``n``
    variables to be meaningful as a character

    EXAMPLES:

    Over `\QQ`, the irreducible character for ``mu`` is the Schur
    function associated to ``mu``, plus garbage terms (Schur
    functions associated to partitions with more than `n` parts) ::

        sage: from sage.algebras.schur_algebra import GL_n_irred_character
        sage: z = GL_n_irred_character(2, [2], QQ)
        sage: sbasis = SymmetricFunctions(QQ).s()
        sage: sbasis(z)
        s[2]

        sage: from sage.algebras.schur_algebra import GL_n_irred_character
        sage: z = GL_n_irred_character(4, [3, 2], QQ)  # long time
        sage: sbasis = SymmetricFunctions(QQ).s()  # long time
        sage: sbasis(z)  # long time
        -5*s[1, 1, 1, 1, 1] + s[3, 2]

    Over a Galois field, the irreducible character for `\mu` will
    in general be smaller.

    In characteristic `p`, for a one-part partition `(r)`, where
    `r= a_0 + p a_1 + p^2 a_2 + \dots`, the result is [Green,
    after 5.5d] the product of `h[a_0], h[a_1]( pbasis[p]), h[a_2]
    ( pbasis[p^2]),\dots,` which is consistent with the following ::

        sage: from sage.algebras.schur_algebra import GL_n_irred_character
        sage: GL_n_irred_character(2, [7], GF(3)) # long time
        m[4, 3] + m[6, 1] + m[7]
    """
    mbasis = SymmetricFunctions(QQ).m()
    r = sum(mu)
    A = SchurAlgebra(n, r, KK)
    M = TensorSpace(n, r, KK)
    S = SymmetricGroupAlgebra(KK, r)

    #make ST the superstandard tableau of shape mu
    from sage.combinat.tableau import from_shape_and_word
    ST = from_shape_and_word(mu, range(1, r + 1), convention='English')

    #make ell the reading word of the highest weight tableau of shape mu
    ell = []
    for i in range(len(mu)):
        for j in range(mu[i]):
            ell.append(i + 1)

    e = M.basis()[Word(ell)]  # the element e_l
    BracC = bracket(r, ST.column_stabilizer(), S)
    f = M.action_by_symmetric_group_algebra(e, BracC)

    # [Green, Theorem 5.3b] says that a basis of the Carter-Lusztig
    # module V_\mu is given by taking this f, and multiplying by all
    # xi_{i,ell} with ell as above and i semistandard.

    carter_lusztig = []
    for i in [Word(flatten(T)) for T in SemistandardTableaux(mu, max_entry=n)]:
        schur_rep = schur_representative_from_index((i, Word(ell)))
        y = M.action_by_Schur_alg(A.basis()[schur_rep], e)
        carter_lusztig.append(y.to_vector())

    #Therefore, we now have carter_lusztig as a list giving the basis
    #of `V_\mu`

    #We want to think of expressing this character as a sum of monomial
    #symmetric functions.

    #We will determine a basis element for each m_\lambda in the
    #character, and we want to keep track of them by \lambda.

    #That means that we only want to pick out the basis elements above for
    #those semistandard words whose content is a partition.

    contents = Partitions(r, max_length=n).list()
    # all partitions of r, length at most n

    # JJ will consist of a list for each element of `contents`,
    # recording the list
    # of semistandard tableaux words with that content

    # graded_basis will consist of the a corresponding basis element

    graded_basis = []
    JJ = []
    for i in range(len(contents)):
        graded_basis.append([])
        JJ.append([])
    for i in [Word(flatten(T), range(1, n + 1))
              for T in SemistandardTableaux(mu, max_entry=n)]:
        con = i.evaluation()
        if all([con[j + 1] <= con[j] for j in range(len(con) - 1)]):
            #don't test whether con is in Partitions, because con could
            #have trailing zeros
            JJ[contents.index(Partition(con))].append(i)
            schur_rep = schur_representative_from_index((i, Word(ell)))
            x = M.action_by_Schur_alg(A.basis()[schur_rep], f)
            graded_basis[contents.index(Partition(con))].append(x.to_vector())

    #There is an inner product on the Carter-Lusztig module V_\mu; its
    #maximal submodule is exactly the kernel of the inner product.

    #Now, for each possible partition content, we look at the graded piece of
    #that degree, and we record how these elements pair with each of the
    #elements of carter_lusztig.

    #The kernel of this pairing is the part of this graded piece which is
    #not in the irreducible module for \mu.

    length = len(carter_lusztig)

    Phi = mbasis.zero()
    for aa in range(len(contents)):
        Mat = []
        for kk in range(len(JJ[aa])):
            temp = []
            for j in range(length):
                temp.append(graded_basis[aa][kk].inner_product(carter_lusztig[j]))
            Mat.append(temp)
        Angle = Matrix(Mat)
        Phi += (len(JJ[aa]) - Angle.nullity()) * mbasis(contents[aa])
    return Phi

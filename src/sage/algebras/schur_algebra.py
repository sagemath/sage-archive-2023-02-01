r"""
Schur algebras for `GL_n`

This file implements:

- Schur algebras for `GL_n` over an arbitrary field.

- The canonical action of the Schur algebra on a tensor power of the standard
  representation.

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
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModule_Tensor
from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.integer_list import IntegerListsLex
from sage.combinat.partition import Partitions, Partition
from sage.combinat.permutation import Permutations
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
from sage.combinat.tableau import SemistandardTableaux
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.matrix.constructor import Matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.rings.all import ZZ, QQ

from copy import copy


def _schur_I_nr_representatives(n, r, element, index):
    """
    Internal function called by :func:`schur_representation_indices`.
    """
    if r == 0:
        return index

    if len(element) == r:
        index.append(copy(element))
        return

    if not element:
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
    `1 \dots n`, where the equivalence relation is simultaneous
    permutation of the two words.  We can therefore fix a
    representative for each equivalence class in which the entries of
    the first word weakly increase, and the entries of the second word
    whose corresponding values in the first word are equal, also
    weakly increase.

    EXAMPLES::

        sage: from sage.algebras.schur_algebra import schur_representative_indices
        sage: schur_representative_indices(2, 2)
        [((1, 1), (1, 1)), ((1, 1), (1, 2)),
         ((1, 1), (2, 2)), ((1, 2), (1, 1)),
         ((1, 2), (1, 2)), ((1, 2), (2, 1)),
         ((1, 2), (2, 2)), ((2, 2), (1, 1)),
         ((2, 2), (1, 2)), ((2, 2), (2, 2))]
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
                    I1 = _schur_I_nr_representatives(n, k, [], [])
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
                    I1 = _schur_I_nr_representatives(n, k, [], [])
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
            basis.append((tuple(e), tuple(v)))

    return basis


def schur_representative_from_index(i0, i1):
    """
    Simultaneously reorder a pair of words to obtain the equivalent
    element of the distinguished basis of the Schur algebra.

    .. SEEALSO::

        :func:`schur_representative_indices`

    INPUT:

    - A pair of words from ``Words(range(1,n+1), r)``

    OUTPUT:

    - The corresponding pair of words ordered correctly.

    EXAMPLES::

        sage: from sage.algebras.schur_algebra import schur_representative_from_index
        sage: schur_representative_from_index([2,1,2,2], [1,3,0,0])
        ((1, 2, 2, 2), (3, 0, 0, 1))
    """
    w = []
    for i, val in enumerate(i0):
        w.append((val, i1[i]))
    w.sort()
    i0 = []
    i1 = []
    for pair in w:
        i0.append(pair[0])
        i1.append(pair[1])
    return (tuple(i0), tuple(i1))


class SchurAlgebra(CombinatorialFreeModule):
    """
    A Schur algebra.

    EXAMPLES::

        sage: S = SchurAlgebra(ZZ, 2, 2); S
        Schur algebra (2, 2) over Integer Ring
    """
    def __init__(self, R, n, r):
        """
        Initialize ``self``.

        TESTS::

            sage: S = SchurAlgebra(ZZ, 2, 2)
            sage: TestSuite(S).run()

        ::

            sage: SchurAlgebra(ZZ, -2, 2)
            Traceback (most recent call last):
            ...
            ValueError: n (=-2) must be a positive integer
            sage: SchurAlgebra(ZZ, 2, -2)
            Traceback (most recent call last):
            ...
            ValueError: r (=-2) must be a non-negative integer
            sage: SchurAlgebra('niet', 2, 2)
            Traceback (most recent call last):
            ...
            ValueError: R (=niet) must be a commutative ring
        """
        if n not in ZZ or n <= 0:
            raise ValueError("n (={}) must be a positive integer".format(n))
        if r not in ZZ or r < 0:
            raise ValueError("r (={}) must be a non-negative integer".format(r))
        if not R in Rings.Commutative():
            raise ValueError("R (={}) must be a commutative ring".format(R))

        self._n = n
        self._r = r

        CombinatorialFreeModule.__init__(self, R,
                                         schur_representative_indices(n, r),
                                         category=AlgebrasWithBasis(R))

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: SchurAlgebra(ZZ, 4, 4)
            Schur algebra (4, 4) over Integer Ring
        """
        msg = "Schur algebra ({}, {}) over {}"
        return msg.format(self._n, self._r, self.base_ring())

    @cached_method
    def one(self):
        """
        Return the one of the algebra.

        EXAMPLES::

            sage: S = SchurAlgebra(ZZ, 2, 2)
            sage: e = S.one(); e
            B[((1, 1), (1, 1))] + B[((1, 2), (1, 2))] + B[((2, 2), (2, 2))]

            sage: x = S.an_element()
            sage: x * e == x
            True
            sage: all(e * x == x for x in S.basis())
            True

            sage: S = SchurAlgebra(ZZ, 4, 4)
            sage: e = S.one()
            sage: x = S.an_element()
            sage: x * e == x
            True
        """
        tt = IntegerListsLex(length=self._r, min_part=1, max_part=self._n,
                             min_slope=0)
        words = [tuple(u) for u in tt]
        return self.sum(self._monomial((w, w)) for w in words)

    def product_on_basis(self, e_ij, e_kl):
        r"""
        Return the product of basis elements.

        EXAMPLES::

            sage: S = SchurAlgebra(QQ, 2, 3)
            sage: B = S.basis()

        If we multiply two basis elements `x` and `y`, such that
        `x[1]` and `y[0]` are not permutations of each other, the
        result is zero ::

            sage: x = B[((1, 1, 1), (1, 1, 2))]
            sage: y = B[((1, 2, 2), (1, 1, 2))]
            sage: x * y
            0

        If we multiply a basis element `x` by a basis element which
        consists of the same tuple repeated twice (on either side),
        the result is either zero (if the previous case applies) or `x` ::

            sage: ww = B[((1, 2, 2), (1, 2, 2))]
            sage: x = B[((1, 2, 2), (1, 1, 2))]
            sage: ww * x
            B[((1, 2, 2), (1, 1, 2))]

        An arbitrary product, on the other hand, may have multiplicities::

            sage: x = B[((1, 1, 1), (1, 1, 2))]
            sage: y = B[((1, 1, 2), (1, 2, 2))]
            sage: x * y
            2*B[((1, 1, 1), (1, 2, 2))]
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
            Z_ijklpq = self.base_ring().zero()
            for s in Permutations([xx for xx in j]):
                if (schur_representative_from_index(e[0], s) == e_ij
                        and schur_representative_from_index(s, e[1]) == e_kl):
                    Z_ijklpq += self.base_ring().one()
            product += Z_ijklpq * b[e]

        return product


class TensorSpace(CombinatorialFreeModule_Tensor):
    """
    This is the ``r``-fold tensor product of an ``n``-dimensional free
    module over ``R``, equipped with an action of the Schur algebra
    `S(n,r)` and the symmetric group `S_r`.

    EXAMPLES::

        sage: from sage.algebras.all import TensorSpace
        sage: TensorSpace(QQ, 2, 3)
        The 3-fold tensor product of a free module of dimension 2
        over Rational Field
    """
    def __init__(self, R, n, r):
        """
        Initialize ``self``.
        """
        C = CombinatorialFreeModule(R, range(1,n+1))
        self._n = n
        self._r = r
        self._sga = SymmetricGroupAlgebra(R, r)
        self._schur = SchurAlgebra(R, n, r)
        CombinatorialFreeModule_Tensor.__init__(self, [C]*r)
        g = self._schur.module_morphism(self._monomial_product, codomain=self)
        self._schur_action = self.module_morphism(g, codomain=self, position=1)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.all import TensorSpace
            sage: TensorSpace(QQ, 2, 3)
            The 3-fold tensor product of a free module of dimension 2
            over Rational Field
        """
        msg = "The {}-fold tensor product of a free module of dimension {}"
        msg += " over {}"
        return msg.format(self._r, self._n, self.base_ring())

    def _monomial_product(self, xi, v):
        """
        Result of acting by the basis element ``xi`` of ``S`` on the
        basis element ``v`` of ``self``.
        """
        x = self.zero()
        for i in CartesianProduct(*[range(1,self._n+1)]*self._r):
            if schur_representative_from_index(i, v) == xi:
                x += self.basis()[tuple(i)]
        return x

    class Element(CombinatorialFreeModule_Tensor.Element):
        def _acted_upon_(self, elt, self_on_left=False):
            """
            Return the action of ``elt`` on ``self``.

            We add the *left* action of the Schur algebra, and the *right*
            actions of the symmetric group algebra and the symemtric group.
            """
            P = self.parent()
            if self_on_left:
                if elt in P._sga:
                    return P.sum_of_terms((tuple([m[i-1] for i in me]), c * ce)
                                          for m,c in self for me,ce in elt)

                if elt in P._sga._indices:
                    return P.sum_of_terms((tuple([m[i-1] for i in elt]), c)
                                          for m,c in self)

            elif elt in P._schur: # self_on_left is False
                return P._schur_action(elt, self)
            return super(TensorSpace.Element, self)._acted_upon_(elt, self_on_left)


def GL_n_irred_character(n, mu, KK):
    r"""
    Return the character of the irreducible module indexed by ``mu``
    of `GL(n)` over the field ``KK``.

    INPUT:

    - ``n`` -- a positive integer
    - ``mu`` -- a partition of at most ``n`` parts
    - ``KK`` -- a field

    OUTPUT:

    a symmetric function which should be interpreted in ``n``
    variables to be meaningful as a character

    EXAMPLES:

    Over `\QQ`, the irreducible character for `\mu` is the Schur
    function associated to `\mu`, plus garbage terms (Schur
    functions associated to partitions with more than `n` parts)::

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
    `r = a_0 + p a_1 + p^2 a_2 + \dots`, the result is [Green,
    after 5.5d] the product of `h[a_0], h[a_1]( pbasis[p]), h[a_2]
    ( pbasis[p^2]), \dots,` which is consistent with the following ::

        sage: from sage.algebras.schur_algebra import GL_n_irred_character
        sage: GL_n_irred_character(2, [7], GF(3)) # long time
        m[4, 3] + m[6, 1] + m[7]
    """
    mbasis = SymmetricFunctions(QQ).m()
    r = sum(mu)
    M = TensorSpace(KK, n, r)
    A = M._schur
    SGA = M._sga

    #make ST the superstandard tableau of shape mu
    from sage.combinat.tableau import from_shape_and_word
    ST = from_shape_and_word(mu, range(1, r + 1), convention='English')

    #make ell the reading word of the highest weight tableau of shape mu
    ell = [i+1 for i,l in enumerate(mu) for dummy in range(l)]

    e = M.basis()[tuple(ell)]  # the element e_l

    # This is the notation `\{X\}` from just before (5.3a) of [GreenPoly]_.
    S = SGA._indices
    BracC = SGA._from_dict({S(x.tuple()): x.sign() for x in ST.column_stabilizer()},
                           remove_zeros=False)
    f = e * BracC #M.action_by_symmetric_group_algebra(e, BracC)

    # [Green, Theorem 5.3b] says that a basis of the Carter-Lusztig
    # module V_\mu is given by taking this f, and multiplying by all
    # xi_{i,ell} with ell as above and i semistandard.

    carter_lusztig = []
    for T in SemistandardTableaux(mu, max_entry=n):
        i = tuple(flatten(T))
        schur_rep = schur_representative_from_index(i, tuple(ell))
        y = A.basis()[schur_rep] * e #M.action_by_Schur_alg(A.basis()[schur_rep], e)
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
    for T in SemistandardTableaux(mu, max_entry=n):
        i = tuple(flatten(T))
        # Get the content of T
        con = [0]*n
        for a in i:
            con[a - 1] += 1
        try:
            P = Partition(con)
            P_index = contents.index(P)
            JJ[P_index].append(i)
            schur_rep = schur_representative_from_index(i, tuple(ell))
            x = A.basis()[schur_rep] * f #M.action_by_Schur_alg(A.basis()[schur_rep], f)
            graded_basis[P_index].append(x.to_vector())
        except ValueError:
            pass

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


r"""
`\mathcal{B}(\infty)` Crystals of Tableaux in Nonexceptional Types and `G_2`

A tableau model for `\mathcal{B}(\infty)`. For more information, see
:class:`~sage.combinat.crystals.infinity_crystals.InfinityCrystalOfTableaux`.

AUTHORS:

- Ben Salisbury: Initial version

- Travis Scrimshaw: Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013 Ben Salisbury <bsalisbury1 at gmail.com>
#                          Travis Scrimshaw <tscrim at ucdavis.edu>
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
#****************************************************************************

from sage.structure.parent import Parent
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten

from sage.combinat.partition import Partition
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.crystals.tensor_product import CrystalOfWords, CrystalOfTableauxElement

class InfinityCrystalOfTableaux(CrystalOfWords):
    r"""
    `\mathcal{B}(\infty)` crystal of tableaux.

    A tableaux model `\mathcal{T}(\infty)` for the crystal
    `\mathcal{B}(\infty)` is introduced by Hong and Lee in [HL08]_. This model
    is currently valid for types `A_n`, `B_n`, `C_n`, `D_n`, and `G_2`, and
    builds on the tableaux model given by Kashiwara and Nakashima [KN94]_ in
    types `A_n`, `B_n`, `C_n`, and `D_n`, and by Kang and Misra [KM94]_ in
    type `G_2`.

    .. NOTE::

        We are using the English convention for our tableaux.

    We say a tableau `T` is *marginally large* if:

    - for each `1 \leq i \leq n`, the leftmost box in the `i`-th row
      from the top in `T` is an `i`-box,

    - for each `1 \leq i \leq n`, the number of `i`-boxes in the `i`-th row
      from the top in `T` is greater than the total number of boxes in the
      `(i+1)`-th row by exactly one.

    We now will describe this tableaux model type-by-type.

    .. rubric:: Type `A_n`

    `\mathcal{T}(\infty)` is the set of marginally large semistandard
    tableaux with exactly `n` rows over the alphabet `\{1 \prec 2 \prec
    \cdots \prec n+1 \}`.

    .. rubric:: Type `B_n`

    `\mathcal{T}(\infty)` is the set of marginally large semistandard
    tableaux with exactly `n` rows over the alphabet `\{1 \prec \cdots
    \prec n \prec 0 \prec \overline{n} \prec \cdots \prec \overline{1} \}`
    and subject to the following constraints:

    - for each `1 \le i \le n`, the contents of the boxes in the
      `i`-th row are `\preceq \overline{i}`,

    - the entry `0` can appear at most once in a single row.

    .. rubric:: Type `C_n`

    `\mathcal{T}(\infty)` is the set of marginally large semistandard
    tableaux with exactly `n` rows over the alphabet `\{1 \prec \cdots
    \prec n \prec \overline{n} \prec \cdots \prec \overline{1} \}` and
    for each `1 \leq i \leq n`, the contents of the boxes in the `i`-th
    row are `\preceq \overline{i}`.

    .. rubric:: Type `D_n`

    `\mathcal{T}(\infty)` is the set of marginally large semistandard
    tableaux with exactly `n-1` rows over the alphabet `\{1 \prec \cdots
    \prec n, \overline{n} \prec \cdots \prec \overline{1} \}` and subject
    to the following constaints:

    - for each `1 \le i \le n`, the contents of the boxes in the `i`-th
      row are `\preceq \overline{i}`,

    - the entries `n` and `\overline{n}` may not appear simultaneously in
      a single row.

    .. rubric:: Type `G_2`

    `\mathcal{T}(\infty)` is the set of marginally large semistandard
    tableaux with exactly `2` rows over the ordered alphabet `\{1 \prec
    2 \prec 3 \prec 0 \prec \overline{3} \prec \overline{2} \prec
    \overline{1}\}` and subject to the following constraints:

    - the contents of the boxes in the first row are `\preceq \overline{i}`,

    - the contents of the boxes in the second row are `\preceq 3`,

    - the entry `0` can appear at most once in the first row and not at
      all in the second row.

    In particular, the shape of the tableaux is not fixed in any instance of
    `\mathcal{T}(\infty)`; the row lengths of a tableau can be arbitrarily long.

    REFERENCES:

    .. [BN10] D. Bump and M. Nakasuji.
       Integration on `p`-adic groups and crystal bases.
       Proc. Amer. Math. Soc. 138(5), pp. 1595--1605.

    .. [LS12] K.-H. Lee and B. Salisbury.
       Young tableaux, canonical bases, and the Gindikin-Karpelevich formula.
       :arXiv:`1205.6006`.

    .. [HL08] J. Hong and H. Lee.
       Young tableaux and crystal `B(\infty)` for finite simple Lie algebras.
       J. Algebra 320, pp. 3680--3693, 2008.

    .. [KM94] S.-J. Kang and K. C. Misra.
       Crystal bases and tensor product decompositions of `U_q(G_2)`-modules.
       J. Algebra 163, pp. 675--691, 1994.

    INPUT:

    - ``cartan_type`` -- One of ``['A',n]``, ``['B',n]``, ``['C',n]``,
      ``['D',n]``, or ``['G',2]``, where ``n`` is a positive integer

    EXAMPLES::

        sage: B = crystals.infinity.Tableaux(['A',2])
        sage: b = B.highest_weight_vector(); b.pp()
        1  1
        2
        sage: b.f_string([2,1,1,2,2,2]).pp()
        1  1  1  1  1  2  3
        2  3  3  3

        sage: B = crystals.infinity.Tableaux(['G',2])
        sage: b = B(rows=[[1,1,1,1,1,2,3,3,0,-3,-1,-1,-1],[2,3,3,3]])
        sage: b.e_string([2,1,1,1,1,1,1]).pp()
        1  1  1  1  2  3  3  3  3 -2 -2 -2
        2  3  3
        sage: b.e_string([2,1,1,1,1,1,1,1])

    We check that a few classical crystals embed into `\mathcal{T}(\infty)`::

        sage: def crystal_test(B, C):
        ....:     g = {C.module_generators[0] : B.module_generators[0]}
        ....:     f = C.crystal_morphism(g)
        ....:     G = B.digraph(subset=[f(x) for x in C])
        ....:     return G.is_isomorphic(C.digraph(), edge_labels=True)
        sage: B = crystals.infinity.Tableaux(['A',2])
        sage: C = crystals.Tableaux(['A',2], shape=[2,1])
        sage: crystal_test(B, C)
        True
        sage: C = crystals.Tableaux(['A',2], shape=[6,2])
        sage: crystal_test(B, C)
        True
        sage: B = crystals.infinity.Tableaux(['B',2])
        sage: C = crystals.Tableaux(['B',2], shape=[3])
        sage: crystal_test(B, C)
        True
        sage: C = crystals.Tableaux(['B',2], shape=[2,1])
        sage: crystal_test(B, C)
        True
        sage: B = crystals.infinity.Tableaux(['C',3])
        sage: C = crystals.Tableaux(['C',3], shape=[2,1])
        sage: crystal_test(B, C)
        True
        sage: B = crystals.infinity.Tableaux(['D',4])
        sage: C = crystals.Tableaux(['D',4], shape=[2])
        sage: crystal_test(B, C)
        True
        sage: C = crystals.Tableaux(['D',4], shape=[1,1,1,1])
        sage: crystal_test(B, C)
        True
        sage: B = crystals.infinity.Tableaux(['G',2])
        sage: C = crystals.Tableaux(['G',2], shape=[3])
        sage: crystal_test(B, C)
        True
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['A',4])
            sage: B2 = crystals.infinity.Tableaux(CartanType(['A',4]))
            sage: B is B2
            True
        """
        cartan_type = CartanType(cartan_type)
        if cartan_type.type() == 'D':
            return InfinityCrystalOfTableauxTypeD(cartan_type)
        return super(InfinityCrystalOfTableaux, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['A',2])
            sage: TestSuite(B).run() # long time
        """
        Parent.__init__( self, category=(HighestWeightCrystals(), InfiniteEnumeratedSets()) )
        self._cartan_type = cartan_type
        self.letters = CrystalOfLetters(cartan_type)
        self.module_generators = (self.module_generator(),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['A',4]); B
            The infinity crystal of tableaux of type ['A', 4]
        """
        return "The infinity crystal of tableaux of type %s"%self._cartan_type

    @cached_method
    def module_generator(self):
        """
        Return the module generator (or highest weight element) of ``self``.

        The module generator is the unique tableau of shape `(n, n-1, \ldots,
        2, 1)` with weight `0`.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: T.module_generator()
            [[1, 1, 1], [2, 2], [3]]
        """
        n = self._cartan_type.rank()
        p = Partition([x for x in reversed(range(1, n+1))])
        # The column canonical tableau, read by columns
        module_generator = flatten([[p[j]-i for i in range(p[j])] for j in range(n)])
        return self(list=[self.letters(x) for x in module_generator])

    def _element_constructor_(self, *args, **options):
        """
        Construct an element of ``self`` from the input data.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: T(rows=[[1,2],[3,4]])
            [[1, 2], [3, 4]]
            sage: T(columns=[[3,1],[4,2]])
            [[1, 2], [3, 4]]
        """
        return self.element_class(self, *args, **options)

    class Element(CrystalOfTableauxElement):
        r"""
        Elements in `\mathcal{B}(\infty)` crystal of tableaux.
        """

        def e(self,i):
            r"""
            Return the action of `\widetilde{e}_i` on ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux(['B',3])
                sage: b = B(rows=[[1,1,1,1,1,1,1,2,0,-3,-1,-1,-1,-1],[2,2,2,2,-2,-2],[3,-3,-3]])
                sage: b.e(3).pp()
                1  1  1  1  1  1  1  2  0 -3 -1 -1 -1 -1
                2  2  2  2 -2 -2
                3  0 -3
                sage: b.e(1).pp()
                1  1  1  1  1  1  1  0 -3 -1 -1 -1 -1
                2  2  2  2 -2 -2
                3 -3 -3
            """
            if i not in self.index_set():
                raise ValueError('i is not in the index set.')
            position = self.positions_of_unmatched_plus(i)
            if position == []:
                return None
            k = position[0]
            ret = self.set_index(k, self[k].e(i))
            if k+i > len(self):
                return ret
            for j in reversed(range(1, i+1)):
                if ret[k+i-j].value != j:
                    return ret
            # We've found a column, so we need to remove it
            for j in range(i):
                ret._list.pop(k)
            return ret

        def f(self, i):
            r"""
            Return the action of `\widetilde{f}_i` on ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux(['C',4])
                sage: b = B.highest_weight_vector()
                sage: b.f(1).pp()
                1  1  1  1  2
                2  2  2
                3  3
                4
                sage: b.f(3).pp()
                1  1  1  1  1
                2  2  2  2
                3  3  4
                4
                sage: b.f(3).f(4).pp()
                1  1  1  1  1
                2  2  2  2
                3  3 -4
                4
            """
            if i not in self.index_set():
                raise ValueError('i is not in the index set.')
            position = self.positions_of_unmatched_minus(i)
            if position == []:
                return None
            k = position[len(position)-1]
            ret = self.set_index(k, self[k].f(i))
            if k+i > len(self):
                return ret
            for j in reversed(range(1,i+1)):
                if self[k+i-j].value != j:
                    return ret
            # We've found a full column, so we'll need to add a new column
            for j in range(i):
                ret._list.insert(k,self.parent().letters(j+1))
            return ret

        def phi(self,i):
            r"""
            Return `\varphi_i` of ``self``.

            Let `T \in \mathcal{B}(\infty)` Define `\varphi_i(T) :=
            \varepsilon_i(T) + \langle h_i, \mathrm{wt}(T) \rangle`, where `h_i`
            is the `i`-th simple coroot and `\mathrm{wt}(T)` is the :meth:`weight`
            of `T`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux("A3")
                sage: [B.highest_weight_vector().f_string([1,3,2,3,1,3,2,1]).phi(i) for i in B.index_set()]
                [-3, 4, -3]

                sage: B = crystals.infinity.Tableaux("G2")
                sage: [B.highest_weight_vector().f_string([2,2,1,2,1,1,1,2]).phi(i) for i in B.index_set()]
                [5, -3]
            """
            P = self.parent().weight_lattice_realization()
            h = P.simple_coroots()
            return self.epsilon(i) + P(self.weight()).scalar(h[i])

        @cached_method
        def weight(self):
            r"""
            Return the weight of ``self``.

            From the definition of a crystal and that the highest weight
            element `b_{\infty}` of `\mathcal{B}(\infty)` is `0`, the weight of
            `T \in \mathcal{B}(\infty)` can be defined as `\mathrm{wt}(T)
            := -\sum_j \alpha_{i_j}` where `\widetilde{e}_{i_1} \cdots
            \widetilde{e}_{i_{\ell}} T = b_{\infty}` and `\{\alpha_i\}` is the
            set of simple roots. (Note that the weight is independent of the
            path chosen to get to the highest weight.)

            However we can also take advantage of the fact that
            `\rho \colon R_{\lambda} \otimes \mathcal{B}(\infty) \longrightarrow
            B(\lambda)`, where `\lambda` is the shape of `T`, preserves the
            tableau representation of `T`. Therefore

            .. MATH::

                \mathrm{wt}(T) = \mathrm{wt}\bigl( \rho(T) \bigr) - \lambda

            where `\mathrm{wt}\bigl( \rho(T) \bigr)` is just the usual weight of
            the tableau `T`.

            Let `\Lambda_i` be the `i`-th fundamental weight. In type `D`, the
            height `n-1` columns corresponds to `\Lambda_{n-1} + \Lambda_n` and
            the in type `B`, the height `n` columns corresponds to
            `2 \Lambda_n`.

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux("C7")
                sage: b = B.highest_weight_vector().f_string([1,6,4,7,4,2,4,6,2,4,6,7,1,2,4,7])
                sage: b.weight()
                (-2, -1, 3, -5, 5, -3, -3)

            Check that the definitions agree::

                sage: P = B.weight_lattice_realization()
                sage: alpha = P.simple_roots()
                sage: b.weight() == -2*alpha[1] - 3*alpha[2] - 5*alpha[4] - 3*alpha[6] - 3*alpha[7]
                True

            Check that it works for type `B`::

                sage: B = crystals.infinity.Tableaux("B2")
                sage: B.highest_weight_vector().weight()
                (0, 0)
                sage: b = B.highest_weight_vector().f_string([1,2,2,2,1,2])
                sage: P = B.weight_lattice_realization()
                sage: alpha = P.simple_roots()
                sage: b.weight() == -2*alpha[1] - 4*alpha[2]
                True

            Check that it works for type `D`::

                sage: B = crystals.infinity.Tableaux("D4")
                sage: B.highest_weight_vector().weight()
                (0, 0, 0, 0)
                sage: b = B.highest_weight_vector().f_string([1,4,4,2,4,3,2,4,1,3,2,4])
                sage: P = B.weight_lattice_realization()
                sage: alpha = P.simple_roots()
                sage: b.weight() == -2*alpha[1] - 3*alpha[2] - 2*alpha[3] - 5*alpha[4]
                True
            """
            P = self.parent().weight_lattice_realization()
            La = P.fundamental_weights()
            cur_col_len = 1
            shape_wt = P.zero()
            n = self.cartan_type().rank()
            ty = self.cartan_type().type()
            for i in range(1, len(self)):
                if self[i-1] < self[i] or (self[i-1].value != 0 and self[i-1] == self[i]):
                    if (cur_col_len == n - 1 and ty == 'D') or \
                            (cur_col_len == n and ty == 'B'):
                        shape_wt += La[n]
                    shape_wt += La[cur_col_len]
                    cur_col_len = 1
                else:
                    cur_col_len += 1
            shape_wt += La[1] # Since we miss the last column (which is always height 1)
            return CrystalOfTableauxElement.weight(self) - shape_wt

        def reduced_form(self):
            r"""
            Return the reduced form of ``self``.

            The reduced form of a tableaux `T \in \mathcal{T}(\infty)` is the
            (not necessarily semistandard) tableaux obtained from `T` by
            removing all `i`-boxes in the `i`-th row, subject to the condition
            that if the row is empty, a `\ast` is put as a placeholder.
            This is described in [BN10]_ and [LS12]_.

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux(['A',3])
                sage: b = B.highest_weight_vector().f_string([2,2,2,3,3,3,3,3])
                sage: b.pp()
                1  1  1  1  1  1  1  1
                2  2  2  2  4  4  4
                3  4  4
                sage: b.reduced_form().pp()
                *
                4  4  4
                4  4
            """
            tab = self.to_tableau()
            for i in range(len(tab)):
                j=0
                while j < len(tab[i]):
                    if tab[i][j] == i+1:
                        tab[i].pop(j)
                        if tab[i] == []:
                            tab[i].append('*')
                    else:
                        j += 1
            return tab

        def seg(self):
            r"""
            Returns the statistic `\mathrm{seg}` of ``self.``

            More precisely, following [LS12]_, define a `k`-segment of a
            tableau `T` in `\mathcal{B}(\infty)` to be a maximal string
            of `k`-boxes in a single row of `T`.  Set `\mathrm{seg}^{\prime}(T)`
            to be the number of `k`-segments in `T`, as `k` varies over
            all possible values.  Then `\mathrm{seg}(T)` is determined
            type-by-type.

            - In types `A_n` and `C_n`, define `\mathrm{seg}(T) :=
              \mathrm{seg}^{\prime}(T)`.

            - In types `B_n` and `G_2`, set `e(T)` to be the number of rows in
              `T` which contain both a `0`-box and an `\overline{\imath}`-box.
              Define `\mathrm{seg}(T) := \mathrm{seg}^{\prime}(T) - e(T)`.

            - In type `D_n`, set `d(T)` to be the number of rows in `T` which
              contain an `\overline{\imath}`-box, but no `n`-box nor
              `\overline{n}`-box. Define `\mathrm{seg}(T) :=
              \mathrm{seg}^{\prime}(T) + d(T)`.

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux(['A',3])
                sage: b = B.highest_weight_vector().f_string([1,3,2,2,3,1,1,3])
                sage: b.pp()
                1  1  1  1  1  1  2  2  4
                2  2  2  2  3
                3  4  4
                sage: b.seg()
                4

                sage: B = crystals.infinity.Tableaux(['D',4])
                sage: b = B(rows=[[1,1,1,1,1,1,3,-2,-1],[2,2,2,4,-2],[3,3],[4]])
                sage: b.pp()
                1  1  1  1  1  1  3 -2 -1
                2  2  2  4 -2
                3  3
                4
                sage: b.seg()
                6

                sage: B = crystals.infinity.Tableaux(['G',2])
                sage: b = B.highest_weight_vector().f_string([2,1,1,1,2,1,2,2,1,2,2,2,1,2,2,1])
                sage: b.pp()
                1  1  1  1  1  1  1  1  2  3  0 -3
                2  3  3  3  3  3  3
                sage: b.seg()
                5
            """
            tab = self.to_tableau()
            segments = []
            for r in range(len(tab)):
                for c in range(len(tab[r])):
                    if tab[r][c] != r+1:
                        if [r,tab[r][c]] not in segments:
                            segments.append([r,tab[r][c]])
            if self.parent().cartan_type().type() == 'B':
                for r in range(len(tab)):
                    for c in range(len(tab[r])):
                        if tab[r][c] == 0 and tab[r][-1] == -r-1:
                            segments.remove([r,tab[r][c]])
            if self.parent().cartan_type().type() == 'D':
                n = self.parent().cartan_type().rank()
                add = []
                for r in range(len(tab)):
                    if tab[r][-1] == -1*(r+1):
                        for c in range(len(tab[r])):
                            if tab[r][c] != n and tab[r][c] != -n:
                                if [r,n] not in add:
                                    add.append([r,n])
                if len(add) > 0:
                    segments.append([r,n])
            if self.parent().cartan_type().type() == 'G':
                for c in range(len(tab[0])):
                    if tab[0][c] == 0 and tab[0][-1] == -1:
                        segments.remove([0,tab[0][c]])
            return len(segments)

        def content(self):
            r"""
            Return the content of ``self``.

            The content `|T|` of `T \in \mathcal{B}(\infty)` is the number of
            blocks added to the highest weight to obtain `T` with any
            `\overline{\imath}`-boxes in the `i`-th row counted with
            multiplicity `2` provided the underlying Cartan type is of type
            `B`, `D`, or `G`.

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux("D5")
                sage: b = B.highest_weight_vector().f_string([5,4,3,1,1,3,4,5,3,4,5,1,4,5,2,3,5,3,2,4])
                sage: b.content()
                13

                sage: B = crystals.infinity.Tableaux("B2")
                sage: b = B(rows=[[1,1,1,1,1,1,2,2,2,-2,-2],[2,0,-2,-2,-2]])
                sage: b.content()
                12

                sage: B = crystals.infinity.Tableaux("C2")
                sage: b = B(rows=[[1,1,1,1,1,1,2,2,2,-2,-2],[2,-2,-2,-2]])
                sage: b.content()
                8
            """
            tab = self.to_tableau()
            count = 0
            ct = self.parent().cartan_type().type()
            for i,row in enumerate(tab):
                for entry in row:
                    if entry == -i-1 and ct in ('B','D','G'):
                        count += 2
                    elif entry != i+1:
                        count += 1
            return count

        def string_parameters(self,word):
            r"""
            Return the string parametrization of ``self`` with respect to ``word``.

            For `w` in the Weyl group, let `R(w)` denote the set of reduced
            expressions for `w`; that is, if `w = s_{i_1} s_{i_2} \cdots
            s_{i_{\ell}}` is a reduced expression, then
            `(i_1, i_2, \ldots, i_{\ell}) \in R(w)`.  For brevity, such words
            are denoted by `\mathbf{i}`.

            For `T \in \mathcal{T}(\infty)` and
            `\mathbf{i} = (i_1, \ldots, i_{\ell}) \in R(w)`,
            the string parametrization `(a_1, \ldots, a_{\ell})` of `T` in the
            direction `\mathbf{i}` is defined recursively by `a_1 =
            \varepsilon_{i_1}(T)` and `a_j = \varepsilon_{i_j}
            (\widetilde{e}_{i_{j-1}}^{\, a_{j-1}} \cdots \widetilde{e}_{i_1}^{\,
            a_1} T)` for `j = 2, \ldots, \ell`. If `w = w_0` is the longest
            element of the Weyl group, then the path determined by the string
            parametrization terminates at the highest weight vector.

            INPUT:

            - ``word`` -- A word in the alphabet of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux("A5")
                sage: b = B(rows=[[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,6,6,6,6,6,6],[2,2,2,2,2,2,2,2,2,4,5,5,5,6],[3,3,3,3,3,3,3,5],[4,4,4,6,6,6],[5,6]])
                sage: b.string_parameters([1,2,1,3,2,1,4,3,2,1,5,4,3,2,1])
                [0, 1, 1, 1, 1, 0, 4, 4, 3, 0, 11, 10, 7, 7, 6]

                sage: B = crystals.infinity.Tableaux("G2")
                sage: b = B(rows=[[1,1,1,1,1,3,3,0,-3,-3,-2,-2,-1,-1,-1,-1],[2,3,3,3]])
                sage: b.string_parameters([2,1,2,1,2,1])
                [5, 13, 11, 15, 4, 4]
                sage: b.string_parameters([1,2,1,2,1,2])
                [7, 12, 15, 8, 10, 0]
            """
            ret = []
            for i in word:
                a = 0
                while self.e(i) is not None:
                    self = self.e(i)
                    a += 1
                ret.append(a)
            return ret

class InfinityCrystalOfTableauxTypeD(InfinityCrystalOfTableaux):
    r"""
    `\mathcal{B}(\infty)` crystal of tableaux for type `D_n`.

    This is the set `\mathcal{T}(\infty)` of marginally large semistandard
    tableaux with exactly `n-1` rows over the alphabet `\{1 \prec \cdots
    \prec n, \overline{n} \prec \cdots \prec \overline{1} \}` and subject
    to the following constraints:

    - for each `1 \le i \le n`, the contents of the boxes in the `i`-th
      row are `\preceq \overline{i}`,

    - the entries `n` and `\overline{n}` may not appear simultaneously in
      a single row.

    For more information, see
    :class:`~sage.combinat.crystals.infinity_crystals.InfinityCrystalOfTableaux`.

    EXAMPLES::

        sage: B = crystals.infinity.Tableaux("D4")
        sage: b = B.highest_weight_vector().f_string([4,3,2,1,4])
        sage: b.pp()
        1  1  1  1  1  1  2
        2  2  2  2  3
        3 -4 -3
        sage: b.weight()
        (-1, 0, -2, -1)
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: B = crystals.infinity.Tableaux(['D',4])
            sage: B2 = crystals.infinity.Tableaux(CartanType(['D',4]))
            sage: B is B2
            True
        """
        return super(InfinityCrystalOfTableauxTypeD, cls).__classcall__(cls, CartanType(cartan_type))

    @cached_method
    def module_generator(self):
        """
        Return the module generator (or highest weight element) of ``self``.

        The module generator is the unique tableau of shape `(n-1, \ldots, 2,
        1)` with weight `0`.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['D',4])
            sage: T.module_generator()
            [[1, 1, 1], [2, 2], [3]]
        """
        n = self._cartan_type.rank()
        p = Partition([x for x in reversed(range(1, n))])
        # The column canonical tableau, read by columns
        module_generator = flatten([[p[j]-i for i in range(p[j])] for j in range(n-1)])
        return self(list=[self.letters(x) for x in module_generator])

    class Element(InfinityCrystalOfTableaux.Element):
        r"""
        Elements in `\mathcal{B}(\infty)` crystal of tableaux for type `D_n`.
        """
        def e(self, i):
            r"""
            Return the action of `\widetilde{e}_i` on ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux(['D',4])
                sage: b = B.highest_weight_vector().f_string([1,4,3,1,2]); b.pp()
                1  1  1  1  2  3
                2  2  2
                3 -3
                sage: b.e(2).pp()
                1  1  1  1  2  2
                2  2  2
                3 -3
            """
            if i not in self.index_set():
                raise ValueError('i is not in the index set.')
            position = self.positions_of_unmatched_plus(i)
            if position == []:
                return None
            k = position[0]
            ret = self.set_index(k, self[k].e(i))
            if i == self.cartan_type().rank():
                i -= 1
            if k+i > len(self):
                return ret
            for j in reversed(range(1, i+1)):
                if ret[k+i-j].value != j:
                    return ret
            # We've found a column, so we need to remove it
            for j in range(i):
                ret._list.pop(k)
            return ret

        def f(self, i):
            r"""
            Return the action of `\widetilde{f}_i` on ``self``.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux(['D',5])
                sage: b = B.highest_weight_vector().f_string([1,4,3,1,5]); b.pp()
                1  1  1  1  1  1  2  2
                2  2  2  2  2
                3  3  3 -5
                4  5
                sage: b.f(1).pp()
                1  1  1  1  1  1  2  2  2
                2  2  2  2  2
                3  3  3 -5
                4  5
                sage: b.f(5).pp()
                1  1  1  1  1  1  2  2
                2  2  2  2  2
                3  3  3 -5
                4 -4
            """
            ret = InfinityCrystalOfTableaux.Element.f(self, i)
            if ret._list[0].value == -self.cartan_type().rank():
                # Exceptional case for f_n where we need to add a new column
                for j in range(i-1):
                    ret._list.insert(0,self.parent().letters(j+1))
            return ret


r"""
Fock Space

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.realizations import Realizations, Category_realization_of_parent

from sage.rings.all import ZZ, QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.partition import (Partition, _Partitions,
        RegularPartitions_all, RegularPartitions_truncated)
from sage.combinat.integer_list import IntegerListsLex
from sage.combinat.partition_tuple import PartitionTuples
from sage.algebras.quantum_groups.q_numbers import q_factorial

class FockSpace(CombinatorialFreeModule):
    r"""
    The Fock space of `U_q(\widehat{\mathfrak{sl}}_n)` of
    type `(\gamma_1, \ldots, \gamma_m)`.

    INPUT:

    - ``n`` -- the value `n`
    - ``r`` -- (default: `[0]`) the type
    - ``q`` -- (optional) the parameter `q`
    - ``base_ring`` -- (optional) the base ring containing ``q``
    """
    @staticmethod
    def __classcall_private__(cls, n, r=[0], q=None, base_ring=None, truncated=None):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: F1 = FockSpace(3, [0])
            sage: F2 = FockSpace(3, 0, q)
            sage: F3 = FockSpace(3, (0,), q, R)
            sage: F1 is F2 and F2 is F3
            True
        """
        if q is None:
            base_ring = PolynomialRing(QQ, 'q')
            q = base_ring.gen(0)
        if base_ring is None:
            base_ring = q.parent()
        base_ring = FractionField(base_ring)
        q = base_ring(q)
        M = IntegerModRing(n)
        if r in ZZ:
            r = (r,)
        r = tuple(map(M, r))
        if truncated is not None:
            return FockSpaceTruncated(n, truncated, q, base_ring)
        return super(FockSpace, cls).__classcall__(cls, n, r, q, base_ring)

    def __init__(self, n, r, q, base_ring, indices=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F = FockSpace(3, [0])
            sage: TestSuite(F).run()
            sage: F = FockSpace(3, [1, 2])
            sage: TestSuite(F).run()
        """
        self._n = n
        self._q = q
        self._r = r
        # If the cell x is above the cell y
        if len(r) == 1: # For partitions
            self._above = lambda x,y: x[0] < y[0]
        else: # For partition tuples
            self._above = lambda x,y: x[0] < y[0] or (x[0] == y[0] and x[1] < y[1])
        self._addable = lambda la,i: filter(lambda x: la.content(*x, multicharge=self._r) == i,
                                            la.outside_corners())
        self._removable = lambda la,i: filter(lambda x: la.content(*x, multicharge=self._r) == i,
                                              la.corners())
        if indices is None:
            indices = PartitionTuples(len(r))
        CombinatorialFreeModule.__init__(self, base_ring, indices,
                                         prefix='', bracket=['|', '>'],
                                         latex_bracket=['\\lvert', '\\rangle'],
                                         monomial_cmp=lambda x,y: -cmp(x,y),
                                         category=ModulesWithBasis(base_ring))

    def __getitem__(self, i):
        """
        Return the basis element indexed by ``i``.

        INPUT:

        - ``i`` -- a partition

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: F[[]]
            |>
            sage: F[1]
            |1>
            sage: F[2,2,1]
            |2, 2, 1>

            sage: F = FockSpace(3, [1, 2])
            sage: F[[], []]
            |[], []>
            sage: F[[2,1], [3,1,1]]
            |[2, 1], [3, 1, 1]>
        """
        if i in ZZ:
            i = [i]
        return self.monomial(self._indices(i))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FockSpace(2)
            Fock space of rank 2 of type (0,) over Fraction Field of
             Univariate Polynomial Ring in q over Rational Field
            sage: FockSpace(4, [2, 0, 1])
            Fock space of rank 4 of type (2, 0, 1) over Fraction Field of
             Univariate Polynomial Ring in q over Rational Field
        """
        return "Fock space of rank {} of type {} over {}".format(self._n, self._r, self.base_ring())

    def _repr_term(self, m):
        """
        Return a representation of the monomial indexed by ``m``.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: F._repr_term(Partition([2,1,1]))
            '|2, 1, 1>'
            sage: F.highest_weight_vector()
            |>
            sage: F = FockSpace(2, [2, 1, 1])
            sage: mg = F.highest_weight_vector(); mg
            |[], [], []>
            sage: mg.f(1).f(1).f(2) / (F._q^-1 + F._q)
            |[1], [1], [1]> + q*|[], [2], [1]> + q^2*|[], [1, 1], [1]>
             + q^3*|[], [1], [2]> + q^4*|[], [1], [1, 1]>
        """
        return '|' + repr(m)[1:-1] + ">" # Strip the outer brackets of m

    def highest_weight_representation(self):
        """
        Return the highest weight representation `B(\Lambda)` in ``self``.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: F.highest_weight_representation()
            Highest weight representation of ['A', 1, 1] of weight Lambda[0]

            sage: F = FockSpace(4, [2,0,1])
            sage: F.highest_weight_representation()
            Highest weight representation of ['A', 3, 1] of weight Lambda[0] + Lambda[1] + Lambda[2]
        """
        return HighestWeightRepresentation(self)

    def highest_weight_vector(self):
        """
        Return the module generator of ``self``.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: F.highest_weight_vector()
            |>
            sage: F = FockSpace(4, [2, 0, 1])
            sage: F.highest_weight_vector()
            |[], [], []>
        """
        if len(self._r) == 1:
            return self.monomial(self._indices([]))
        return self.monomial( self._indices([[]]*len(self._r)) )

    class Element(CombinatorialFreeModuleElement):
        """
        An element in the Fock space.
        """
        def e(self, i):
            """
            Apply the action of `e_i` on ``self``.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: F[2,1,1].e(1)
                1/q*|1, 1, 1>
                sage: F[2,1,1].e(0)
                |2, 1>
                sage: F[2,1,1].e(0).e(1)
                |2> + q*|1, 1>
                sage: F[2,1,1].e(0).e(1).e(1)
                ((q^2+1)/q)*|1>
                sage: F[2,1,1].e(0).e(1).e(1).e(1)
                0
                sage: F[2,1,1].e(0).e(1).e(1).e(0)
                ((q^2+1)/q)*|>
                sage: F[2,1,1].e(1).e(0).e(1).e(0)
                1/q*|>

                sage: F = FockSpace(4, [2, 0, 1])
                sage: F[[2,1],[1],[2]]
                |[2, 1], [1], [2]>
                sage: F[[2,1],[1],[2]].e(2)
                |[2, 1], [1], [1]>
                sage: F[[2,1],[1],[2]].e(1)
                1/q*|[2], [1], [2]>
                sage: F[[2,1],[1],[2]].e(0)
                1/q*|[2, 1], [], [2]>
                sage: F[[2,1],[1],[2]].e(3)
                1/q^2*|[1, 1], [1], [2]>
                sage: F[[2,1],[1],[2]].e(3).e(2).e(1)
                1/q^2*|[1, 1], [1], []> + 1/q^2*|[1], [1], [1]>
                sage: F[[2,1],[1],[2]].e(3).e(2).e(1).e(0).e(1).e(2)
                2/q^3*|[], [], []>
            """
            P = self.parent()
            N_left = lambda la, x, i: \
                sum(1 for y in P._addable(la, i) if P._above(x, y)) \
                - sum(1 for y in P._removable(la, i) if P._above(x, y))
            q = P._q
            return P.sum_of_terms( [( la.remove_cell(*x), c * q**(-N_left(la, x, i)) )
                                    for la,c in self for x in P._removable(la, i)] )

        def f(self, i):
            """
            Apply the action of `f_i` on ``self``.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: mg = F.highest_weight_vector()
                sage: mg.f(0)
                |1>
                sage: mg.f(0).f(1)
                |2> + q*|1, 1>
                sage: mg.f(0).f(0)
                0
                sage: mg.f(0).f(1).f(1)
                ((q^2+1)/q)*|2, 1>
                sage: mg.f(0).f(1).f(0)
                |3> + q*|1, 1, 1>

                sage: F = FockSpace(4, [2, 0, 1])
                sage: mg = F.highest_weight_vector()
                sage: mg.f(0)
                |[], [1], []>
                sage: mg.f(2)
                |[1], [], []>
                sage: mg.f(1)
                |[], [], [1]>
                sage: mg.f(1).f(0)
                |[], [1], [1]> + q*|[], [], [1, 1]>
                sage: mg.f(0).f(1)
                |[], [2], []> + q*|[], [1], [1]>
                sage: mg.f(0).f(1).f(3)
                |[], [2, 1], []> + q*|[], [1, 1], [1]>
                sage: mg.f(3)
                0
            """
            P = self.parent()
            N_right = lambda la, x, i: \
                sum(1 for y in P._addable(la, i) if P._above(y, x)) \
                - sum(1 for y in P._removable(la, i) if P._above(y, x))
            q = P._q
            return P.sum_of_terms( [(la.add_cell(*x), c * q**N_right(la, x, i))
                                    for la,c in self for x in P._addable(la, i)] )

        def h(self, i):
            """
            Apply the action of `h_i` on ``self``.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: F[2,1,1].h(0)
                q*|2, 1, 1>
                sage: F[2,1,1].h(1)
                |2, 1, 1>

                sage: F = FockSpace(4, [2,0,1])
                sage: F[[2,1],[1],[2]].h(0)
                q^2*|[2, 1], [1], [2]>
                sage: F[[2,1],[1],[2]].h(1)
                |[2, 1], [1], [2]>
                sage: F[[2,1],[1],[2]].h(2)
                |[2, 1], [1], [2]>
                sage: F[[2,1],[1],[2]].h(3)
                q*|[2, 1], [1], [2]>
            """
            P = self.parent()
            N_i = lambda la, i: len(P._addable(la, i)) - len(P._removable(la, i))
            q = P._q
            return P.sum_of_terms([(la, c * q**N_i(la, i))
                                   for la,c in self])

        def d(self):
            """
            Apply the action of `d` on ``self``.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: F.highest_weight_vector().d()
                |>
                sage: F[2,1,1].d()
                q^2*|2, 1, 1>
                sage: F[5,3,3,1,1,1].d()
                q^7*|5, 3, 3, 1, 1, 1>

                sage: F = FockSpace(4, [2,0,1])
                sage: F.highest_weight_vector().d()
                |[], [], []>
                sage: F[[2,1],[1],[2]].d()
                q*|[2, 1], [1], [2]>
                sage: F[[4,2,2,1],[1],[5,2]].d()
                q^5*|[4, 2, 2, 1], [1], [5, 2]>
            """
            P = self.parent()
            q = P._q
            return P.sum_of_terms(
                [(la, c * q**sum(1 for x in la.cells() if la.content(*x, multicharge=P._r) == 0))
                 for la,c in self])

class HighestWeightRepresentation(Parent, UniqueRepresentation):
    r"""
    The highest weight representation `B(\Lambda)` of
    `U_q(\widehat{\mathfrak{sl}}_n)`.

    We realize this a subspace of the corresponding
    :class:`Fock space <FockSpace>`. We have two bases:

    - The approximation basis which comes from LLT(-type) algorithms.
    - The lower global crystal basis.

    .. TODO::

        The basis indexing set actually gives the whole tensor product space
        `\bigotimes_i B(\Lambda_{\gamma_i})`. Once :trac:`15584` is done,
        we should change this to a new class ``TensorProductRepresentation``
        and have the highest weight representation living as a subspace.

    REFERENCES:

    .. [LLT1996] Alain Lascoux, Bernard Leclerc, and Jean-Yves Thibon.
       *Hecke algebras at roots of unity and crystal bases of quantum affine
       algebras*. Commun. Math. Phys. **181** (1996). pp 205-263.

    .. [Fayers2010] Matthew Fayers. *An LLT-type algorithm for computing
       higher-level canonical bases*. (2010). :arxiv:`0908.1749v3`.
    """
    @staticmethod
    def __classcall_private__(cls, n, r=[0], q=None, base_ring=None):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: B1 = F.highest_weight_representation()
            sage: from sage.algebras.quantum_groups.fock_space import HighestWeightRepresentation
            sage: B2 = HighestWeightRepresentation(2)
            sage: R.<q> = QQ[]
            sage: B3 = HighestWeightRepresentation(2, [0], q, R)
            sage: B1 is B2 and B2 is B3
            True
        """
        if isinstance(n, FockSpace):
            return super(HighestWeightRepresentation, cls).__classcall__(cls, n)
        F = FockSpace(n, tuple(r), q, base_ring)
        return super(HighestWeightRepresentation, cls).__classcall__(cls, F)

    def __init__(self, fock_space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: B = F.highest_weight_representation()
            sage: x = B.a_realization() # This should be removed, but is needed so the test suite correctly obtains an element
            sage: TestSuite(B).run()
        """
        self._fock = fock_space
        R = fock_space.base_ring()

        Parent.__init__(self, base=R, category=ModulesWithBasis(R).WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: F.highest_weight_representation()
            Highest weight representation of ['A', 1, 1] of weight Lambda[0]
        """
        from sage.combinat.root_system.cartan_type import CartanType
        ct = CartanType(['A', self._fock._n - 1, 1])
        P = ct.root_system().weight_lattice()
        wt = P.sum_of_monomials(self._fock._r)
        return "Highest weight representation of {} of weight {}".format(ct, wt)

    def a_realization(self):
        """
        Return a realization of ``self``.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: B = F.highest_weight_representation()
            sage: B.a_realization()
            Highest weight representation of ['A', 1, 1] of weight Lambda[0] in the lower global crystal basis
        """
        return self.G()

    class A(CombinatorialFreeModule, BindableClass):
        """
        The `A` basis of the Fock space which is the approximation of the
        lower global crystal basis.

        EXAMPLES:

        We construct Example 6.5 and 6.7 in [LLT1996]_::

            sage: F = FockSpace(2)
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: A = B.A()
            sage: F(A[5])
            |5> + |3, 2> + 2*q*|3, 1, 1> + q^2*|2, 2, 1> + q^2*|1, 1, 1, 1, 1>
            sage: F(A[4,1])
            |4, 1> + q*|2, 1, 1, 1>
            sage: F(A[3,2])
            |3, 2> + q*|3, 1, 1> + q^2*|2, 2, 1>
            sage: F(G[5])
            |5> + q*|3, 1, 1> + q^2*|1, 1, 1, 1, 1>

        We construct the examples in Section 5.1 of [Fayers2010]_::

            sage: F = FockSpace(2, [0, 0])
            sage: B = F.highest_weight_representation()
            sage: A = B.A()
            sage: F(A[[2,1],[1]])
            |[2, 1], [1]> + q*|[2], [2]> + q^2*|[2], [1, 1]> + q^2*|[1, 1], [2]>
             + q^3*|[1, 1], [1, 1]> + q^4*|[1], [2, 1]>
            sage: F(A[[4],[]])
            |[4], []> + q*|[3, 1], []> + q*|[2, 1, 1], []> + (q^2+1)*|[2, 1], [1]>
             + 2*q*|[2], [2]> + 2*q^2*|[2], [1, 1]> + q^2*|[1, 1, 1, 1], []>
             + 2*q^2*|[1, 1], [2]> + 2*q^3*|[1, 1], [1, 1]> + (q^4+q^2)*|[1], [2, 1]>
             + q^2*|[], [4]> + q^3*|[], [3, 1]> + q^3*|[], [2, 1, 1]> + q^4*|[], [1, 1, 1, 1]>
        """
        def __init__(self, basic):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: TestSuite(A).run()
            """
            self._basis_name = "approximation"
            indices = PartitionTuples(level=len(basic._fock._r), regular=basic._fock._n)
            CombinatorialFreeModule.__init__(self, basic.base_ring(), indices,
                                             prefix='A', bracket=False,
                                             monomial_cmp=lambda x,y: -cmp(x,y),
                                             category=HighestWeightRepresentationBases(basic))
            self.module_morphism(self._A_to_fock_basis,
                                 triangular='upper', unitriangular=True,
                                 codomain=basic._fock).register_as_coercion()

        @cached_method
        def _A_to_fock_basis(self, la):
            """
            Return the `A` basis indexed by ``la`` in the natural basis.

            EXAMPLES::

                sage: F = FockSpace(3)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: A._A_to_fock_basis(Partition([3]))
                |3> + q*|2, 1>
                sage: A._A_to_fock_basis(Partition([2,1]))
                |2, 1> + q*|1, 1, 1>
            """
            fock = self.realization_of()._fock

            if la.size() == 0:
                return fock.highest_weight_vector()

            if len(fock._r) > 1:
                # Find the first to be non-empty partition
                # Note this is one more than the first non-empty partition
                k = 1
                for p in la:
                    if p.size() != 0:
                        break
                    k += 1

                # Reduce down to the lower level Fock space and do the computation
                #   and then lift back up to us by prepending empty partitions
                if k == len(fock._r): # This means we get the empty partition
                    cur = self.highest_weight_vector()
                else:
                    F = FockSpace(fock._n, fock._r[k:], fock._q, fock.base_ring())
                    Gp = F.highest_weight_representation().G()
                    if k + 1 == len(fock._r):
                        cur = Gp._G_to_fock_basis(Gp._indices(la[k]))
                        cur = fock.sum_of_terms((fock._indices([[]]*k + [p]), c) for p,c in cur)
                    else:
                        cur = Gp._G_to_fock_basis(Gp._indices(la[k:]))
                        cur = fock.sum_of_terms((fock._indices([[]]*k + list(pt)), c) for pt,c in cur)
                la = la[k-1]
            else:
                cur = fock.highest_weight_vector()

            # Get the ladders and apply it to the current element
            corners = la.corners()
            cells = set(la.cells())
            q = fock._q
            k = fock._n - 1 # This is sl_{k+1}
            r = fock._r[0]
            b = ZZ.zero()
            while any(c[1]*k + c[0] >= b for c in corners): # While there is some cell left to count
                power = 0
                i = -b + r # This will be converted to a mod n number
                for x in range(0, b // k + 1):
                    if (b-x*k, x) in cells:
                        power += 1
                        cur = cur.f(i)
                cur /= q_factorial(power, q)
                b += 1
            return cur

    approximation = A

    class G(CombinatorialFreeModule, BindableClass):
        """
        The lower global crystal basis living inside of Fock space.

        EXAMPLES:

        We construct some of the tables/entries given in Section 10
        of [LLT1996]_. For `\widehat{\mathfrak{sl}}_2`::

            sage: F = FockSpace(2)
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: F(G[2])
            |2> + q*|1, 1>
            sage: F(G[3])
            |3> + q*|1, 1, 1>
            sage: F(G[2,1])
            |2, 1>
            sage: F(G[4])
            |4> + q*|3, 1> + q*|2, 1, 1> + q^2*|1, 1, 1, 1>
            sage: F(G[3,1])
            |3, 1> + q*|2, 2> + q^2*|2, 1, 1>
            sage: F(G[5])
            |5> + q*|3, 1, 1> + q^2*|1, 1, 1, 1, 1>
            sage: F(G[4,2])
            |4, 2> + q*|4, 1, 1> + q*|3, 3> + q^2*|3, 1, 1, 1> + q^2*|2, 2, 2> + q^3*|2, 2, 1, 1>
            sage: F(G[4,2,1])
            |4, 2, 1> + q*|3, 3, 1> + q^2*|3, 2, 2> + q^3*|3, 2, 1, 1>
            sage: F(G[6,2])
            |6, 2> + q*|6, 1, 1> + q*|5, 3> + q^2*|5, 1, 1, 1> + q*|4, 3, 1>
             + q^2*|4, 2, 2> + (q^3+q)*|4, 2, 1, 1> + q^2*|4, 1, 1, 1, 1>
             + q^2*|3, 3, 1, 1> + q^3*|3, 2, 2, 1> + q^3*|3, 1, 1, 1, 1, 1>
             + q^3*|2, 2, 2, 1, 1> + q^4*|2, 2, 1, 1, 1, 1>
            sage: F(G[5,3,1])
            |5, 3, 1> + q*|5, 2, 2> + q^2*|5, 2, 1, 1> + q*|4, 4, 1>
             + q^2*|4, 2, 1, 1, 1> + q^2*|3, 3, 3> + q^3*|3, 3, 1, 1, 1>
             + q^3*|3, 2, 2, 2> + q^4*|3, 2, 2, 1, 1>
            sage: F(G[4,3,2,1])
            |4, 3, 2, 1>
            sage: F(G[7,2,1])
            |7, 2, 1> + q*|5, 2, 1, 1, 1> + q^2*|3, 2, 1, 1, 1, 1, 1>
            sage: F(G[10,1])
            |10, 1> + q*|8, 1, 1, 1> + q^2*|6, 1, 1, 1, 1, 1>
             + q^3*|4, 1, 1, 1, 1, 1, 1, 1> + q^4*|2, 1, 1, 1, 1, 1, 1, 1, 1, 1>
            sage: F(G[6,3,2])
            |6, 3, 2> + q*|6, 3, 1, 1> + q^2*|6, 2, 2, 1> + q^3*|5, 3, 2, 1>
             + q*|4, 3, 2, 1, 1> + q^2*|4, 3, 1, 1, 1, 1>
             + q^3*|4, 2, 2, 1, 1, 1> + q^4*|3, 3, 2, 1, 1, 1>
            sage: F(G[5,3,2,1])
            |5, 3, 2, 1> + q*|4, 4, 2, 1> + q^2*|4, 3, 3, 1> + q^3*|4, 3, 2, 2> + q^4*|4, 3, 2, 1, 1>

        For `\widehat{\mathfrak{sl}}_3`::

            sage: F = FockSpace(3)
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: F(G[2])
            |2>
            sage: F(G[1,1])
            |1, 1>
            sage: F(G[3])
            |3> + q*|2, 1>
            sage: F(G[2,1])
            |2, 1> + q*|1, 1, 1>
            sage: F(G[4])
            |4> + q*|2, 2>
            sage: F(G[3,1])
            |3, 1>
            sage: F(G[2,2])
            |2, 2> + q*|1, 1, 1, 1>
            sage: F(G[2,1,1])
            |2, 1, 1>
            sage: F(G[5])
            |5> + q*|2, 2, 1>
            sage: F(G[2,2,1])
            |2, 2, 1> + q*|2, 1, 1, 1>
            sage: F(G[4,1,1])
            |4, 1, 1> + q*|3, 2, 1> + q^2*|3, 1, 1, 1>
            sage: F(G[5,2])
            |5, 2> + q*|4, 3> + q^2*|4, 2, 1>
            sage: F(G[8])
            |8> + q*|5, 2, 1> + q*|3, 3, 1, 1> + q^2*|2, 2, 2, 2>
            sage: F(G[7,2])
            |7, 2> + q*|4, 2, 2, 1>
            sage: F(G[6,2,2])
            |6, 2, 2> + q*|6, 1, 1, 1, 1> + q*|4, 4, 2> + q^2*|3, 3, 2, 1, 1>

        For `\widehat{\mathfrak{sl}}_4`::

            sage: F = FockSpace(4)
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: F(G[4])
            |4> + q*|3, 1>
            sage: F(G[3,1])
            |3, 1> + q*|2, 1, 1>
            sage: F(G[2,2])
            |2, 2>
            sage: F(G[2,1,1])
            |2, 1, 1> + q*|1, 1, 1, 1>
            sage: F(G[3,2])
            |3, 2> + q*|2, 2, 1>
            sage: F(G[2,2,2])
            |2, 2, 2> + q*|1, 1, 1, 1, 1, 1>
            sage: F(G[6,1])
            |6, 1> + q*|4, 3>
            sage: F(G[3,2,2,1])
            |3, 2, 2, 1> + q*|3, 1, 1, 1, 1, 1> + q*|2, 2, 2, 2> + q^2*|2, 1, 1, 1, 1, 1, 1>
            sage: F(G[7,2])
            |7, 2> + q*|6, 2, 1> + q*|5, 4> + q^2*|5, 3, 1>
            sage: F(G[5,2,2,1])
            |5, 2, 2, 1> + q*|5, 1, 1, 1, 1, 1> + q*|4, 2, 2, 1, 1> + q^2*|4, 2, 1, 1, 1, 1>

        We construct the examples in Section 5.1 of [Fayers2010]_::

            sage: F = FockSpace(2, [0, 0])
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: F(G[[2,1],[1]])
            |[2, 1], [1]> + q*|[2], [2]> + q^2*|[2], [1, 1]> + q^2*|[1, 1], [2]>
             + q^3*|[1, 1], [1, 1]> + q^4*|[1], [2, 1]>
            sage: F(G[[4],[]])
            |[4], []> + q*|[3, 1], []> + q*|[2, 1, 1], []> + q^2*|[2, 1], [1]>
             + q*|[2], [2]> + q^2*|[2], [1, 1]> + q^2*|[1, 1, 1, 1], []> + q^2*|[1, 1], [2]>
             + q^3*|[1, 1], [1, 1]> + q^2*|[1], [2, 1]> + q^2*|[], [4]> + q^3*|[], [3, 1]>
             + q^3*|[], [2, 1, 1]> + q^4*|[], [1, 1, 1, 1]>
        """
        def __init__(self, basic):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: B = F.highest_weight_representation()
                sage: G = B.G()
                sage: TestSuite(G).run()
            """
            self._basis_name = "lower global crystal"
            indices = PartitionTuples(level=len(basic._fock._r), regular=basic._fock._n)
            CombinatorialFreeModule.__init__(self, basic.base_ring(), indices,
                                             prefix='G', bracket=False,
                                             monomial_cmp=lambda x,y: -cmp(x,y),
                                             category=HighestWeightRepresentationBases(basic))
            self.module_morphism(self._G_to_fock_basis,
                                 triangular='upper', unitriangular=True,
                                 codomain=basic._fock).register_as_coercion()

        @cached_method
        def _G_to_fock_basis(self, la):
            """
            Return the `G` basis indexed by ``la`` in the natural basis.

            EXAMPLES::

                sage: F = FockSpace(3)
                sage: B = F.highest_weight_representation()
                sage: G = B.G()
                sage: G._G_to_fock_basis(Partition([3]))
                |3> + q*|2, 1>
                sage: G._G_to_fock_basis(Partition([2,1]))
                |2, 1> + q*|1, 1, 1>
            """
            # Special case for the empty partition
            if la.size() == 0:
                return self.realization_of()._fock.highest_weight_vector()

            # Special case for empty leading partitions
            if len(self.realization_of()._fock._r) > 1 and la[0].size() == 0:
                fock = self.realization_of()._fock
                # Find the first non-empty partition
                k = 0
                for p in la:
                    if p.size() != 0:
                        break
                    k += 1

                # Reduce down to the lower level Fock space and do the computation
                #   and then lift back up by prepending empty partitions
                # Note that this will never be for the empty partition, which
                #   is already taken care of
                F = FockSpace(fock._n, fock._r[k:], fock._q, fock.base_ring())
                Gp = F.highest_weight_representation().G()
                if k + 1 == len(fock._r):
                    cur = Gp._G_to_fock_basis(Gp._indices(la[k]))
                    return fock.sum_of_terms((fock._indices([[]]*k + [p]), c) for p,c in cur)
                cur = Gp._G_to_fock_basis(Gp._indices(la[k:]))
                return fock.sum_of_terms((fock._indices([[]]*k + list(pt)), c) for pt,c in cur)

            cur = self.realization_of().A()._A_to_fock_basis(la)
            dom_cmp = lambda x,y: -1 if y.dominates(x) else 1
            s = cur.support()
            s.sort(cmp=dom_cmp) # Sort via dominance order
            s.pop() # Remove the largest

            q = self.realization_of()._fock._q
            while len(s) != 0:
                mu = s.pop()
                d = cur[mu].denominator()
                k = d.degree()
                n = cur[mu].numerator()
                if k != 0 or n.constant_coefficient() != 0:
                    gamma = sum(n[i] * (q**(i-k) + q**(k-i))
                                for i in range(min(n.degree(), k)))
                    gamma += n[k]
                    cur -= gamma * self._G_to_fock_basis(mu)

                    # Add any new support elements
                    for x in cur.support():
                        if x == mu or not mu.dominates(x): # Add only things (strictly) dominated by mu
                            continue
                        for i in reversed(range(len(s))):
                            if not s[i].dominates(x):
                                s.insert(i+1, x)
                                break
            return cur

        #def _repr_(self):
        #    """
        #    Return a string representation of ``self``.
        #    """
        #    R = self.realization_of()
        #    return "Lower global crystal basis of type {} and weight {}".format(R._cartan_type, R._wt)

    lower_global_crystal_basis = G

class HighestWeightRepresentationBases(Category_realization_of_parent):
    r"""
    The category of bases of a highest weight representation in a (truncated)
    Fock space.
    """
    def __init__(self, base):
        r"""
        Initialize the bases of a highest weight representation.

        INPUT:

        - ``base`` -- a highest weight representation

        TESTS::

            sage: from sage.algebras.quantum_groups.fock_space import HighestWeightRepresentationBases
            sage: F = FockSpace(2)
            sage: B = F.highest_weight_representation()
            sage: bases = HighestWeightRepresentationBases(B)
            sage: TestSuite(bases).run()
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Returns the representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.fock_space import HighestWeightRepresentationBases
            sage: F = FockSpace(2)
            sage: B = F.highest_weight_representation()
            sage: HighestWeightRepresentationBases(B)
            Category of bases of Highest weight representation of ['A', 1, 1] of weight Lambda[0]
        """
        return "Category of bases of {}".format(self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.fock_space import HighestWeightRepresentationBases
            sage: F = FockSpace(2)
            sage: B = F.highest_weight_representation()
            sage: bases = HighestWeightRepresentationBases(B)
            sage: bases.super_categories()
            [Category of vector spaces with basis over Fraction Field of Univariate Polynomial Ring in q over Rational Field,
             Category of realizations of Highest weight representation of ['A', 1, 1] of weight Lambda[0]]
        """
        return [ModulesWithBasis(self.base().base_ring()), Realizations(self.base())]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of Fock space.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: B = F.highest_weight_representation()
                sage: B.A()
                Highest weight representation of ['A', 1, 1] of weight Lambda[0] in the approximation basis
                sage: B.G()
                Highest weight representation of ['A', 1, 1] of weight Lambda[0] in the lower global crystal basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        @cached_method
        def highest_weight_vector(self):
            """
            Return the highest weight vector of ``self``.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: A.highest_weight_vector()
                A[]
                sage: G = B.G()
                sage: G.highest_weight_vector()
                G[]
            """
            level = len(self.realization_of()._fock._r)
            if level == 1:
                return self.monomial(self._indices([]))
            return self.monomial(self._indices([[]]*level))

        def __getitem__(self, i):
            """
            Return the basis element indexed by ``i``.

            INPUT:

            - ``i`` -- a partition

            EXAMPLES::

                sage: F = FockSpace(3)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: A[[]]
                A[]
                sage: A[4]
                A[4]
                sage: A[2,2,1]
                A[2, 2, 1]
                sage: G = B.G()
                sage: G[[]]
                G[]
                sage: G[4]
                G[4]
                sage: G[2,2,1]
                G[2, 2, 1]

            For higher levels::

                sage: F = FockSpace(2, [0, 0])
                sage: B = F.highest_weight_representation()
                sage: G = B.G()
                sage: G[[2,1],[1]]
                G([2, 1], [1])

            TESTS::

                sage: F = FockSpace(3)
                sage: A = F.highest_weight_representation().A()
                sage: A[2,2,2,1]
                Traceback (most recent call last):
                ...
                ValueError: (2, 2, 2, 1) is not an element of 3-Regular Partitions

                sage: F = FockSpace(3, [0, 0])
                sage: A = F.highest_weight_representation().A()
                sage: A[[], [2,2,2,1]]
                Traceback (most recent call last):
                ...
                ValueError: [[], [2, 2, 2, 1]] is not a 3-Regular partition tuples of level 2
            """
            if i in ZZ:
                i = [i]

            i = self._indices(i)
            if i.size() == 0:
                return self.highest_weight_vector()
            return self.monomial(i)

###############################################################################
## Truncated Fock space

class FockSpaceTruncated(FockSpace):
    """
    This is the Fock space given by partitions of length no more than `k`.

    This can be formed as the quotient `\mathcal{F} / \mathcak{F}_k` where
    `\mathcal{F}_k` is the submodule spanned by all diagrams of length
    (strictly) more than `k`.

    EXAMPLES::

        sage: F = FockSpace(2, truncated=2)
        sage: mg = F.highest_weight_vector()
        sage: mg.f(0)
        |1>
        sage: mg.f(0).f(1)
        |2> + q*|1, 1>
        sage: mg.f(0).f(1).f(0)
        |3>

    Compare this to the full Fock space::

        sage: F = FockSpace(2)
        sage: mg = F.highest_weight_vector()
        sage: mg.f(0).f(1).f(0)
        |3> + q*|1, 1, 1>
    """
    @staticmethod
    def __classcall_private__(cls, n, k, q=None, base_ring=None):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: F1 = FockSpace(3, truncated=2)
            sage: F2 = FockSpace(3, q=q, truncated=2)
            sage: F3 = FockSpace(3, q=q, base_ring=R, truncated=2)
            sage: F1 is F2 and F2 is F3
            True
            sage: from sage.algebras.quantum_groups.fock_space import FockSpaceTruncated
            sage: F4 = FockSpaceTruncated(3, 2, q, R)
            sage: F1 is F4
            True
        """
        if q is None:
            base_ring = PolynomialRing(QQ, 'q')
            q = base_ring.gen(0)
        if base_ring is None:
            base_ring = q.parent()
        base_ring = FractionField(base_ring)
        q = base_ring(q)
        return super(FockSpace, cls).__classcall__(cls, n, k, q, base_ring)

    def __init__(self, n, k, q, base_ring):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F = FockSpace(2, truncated=3)
            sage: TestSuite(F).run()
        """
        M = IntegerModRing(n)
        self._k = k
        indices = RegularPartitions_truncated(n, k)
        FockSpace.__init__(self, n, (M(0),), q, base_ring, indices)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FockSpace(2, truncated=3)
            Fock space of rank 2 truncated at 3 over Fraction Field of
             Univariate Polynomial Ring in q over Rational Field
        """
        return "Fock space of rank {} truncated at {} over {}".format(self._n, self._k, self.base_ring())

    def highest_weight_representation(self):
        """
        Return the highest weight representation `B(\Lambda)` in ``self``.

        EXAMPLES::

            sage: F = FockSpace(2, truncated=3)
            sage: F.highest_weight_representation()
            Highest weight representation of ['A', 1, 1] of weight Lambda[0] truncated at 3
        """
        return HighestWeightRepresentationTruncated(self)

    class Element(FockSpace.Element):
        """
        An element in the trucated Fock space.
        """
        def f(self, i):
            """
            Apply the action of `f_i` on ``self``.

            EXAMPLES::

                sage: F = FockSpace(2, truncated=3)
                sage: mg = F.highest_weight_vector()
                sage: mg.f(0)
                |1>
                sage: mg.f(0).f(1)
                |2> + q*|1, 1>
                sage: mg.f(0).f(1).f(0)
                |3> + q*|1, 1, 1>
                sage: mg.f(0).f(1).f(0).f(0)
                0
                sage: mg.f(0).f(1).f(0).f(1)
                |4> + q*|3, 1> + q*|2, 1, 1>
                sage: mg.f(0).f(1).f(0).f(1).f(0)
                |5> + |3, 2> + 2*q*|3, 1, 1> + q^2*|2, 2, 1>
            """
            P = self.parent()
            N_right = lambda la, x, i: \
                sum(1 for y in P._addable(la, i) if P._above(y, x)) \
                - sum(1 for y in P._removable(la, i) if P._above(y, x))
            q = P._q
            return P.sum_of_terms( [(la.add_cell(*x), c * q**N_right(la, x, i))
                                    for la,c in self for x in P._addable(la, i) if x[0] < P._k] )


class HighestWeightRepresentationTruncated(Parent, UniqueRepresentation):
    """
    The highest weight representation `B(\Lambda)` of
    `U_q(\widehat{\mathfrak{sl}}_n)` truncated at `k`.

    We realize this a subspace of the corresponding
    :class:`truncated Fock space <FockSpaceTruncated>`. We have two bases:

    - The approximation basis which comes from LLT(-type) algorithms.
    - The lower global crystal basis.

    REFERENCES:

    .. [GW1998] Frederick M. Goodman and Hans Wenzl. *Crystal bases of quantum
       affine algebras and affine Kazhdan-Lusztig polyonmials*. (1998).
       :arxiv:`math/9807014v1`.
    """
    @staticmethod
    def __classcall_private__(cls, n, k=None, q=None, base_ring=None):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: F = FockSpace(2, truncated=3)
            sage: B1 = F.highest_weight_representation()
            sage: from sage.algebras.quantum_groups.fock_space import HighestWeightRepresentationTruncated
            sage: B2 = HighestWeightRepresentationTruncated(2, 3)
            sage: R.<q> = QQ[]
            sage: B3 = HighestWeightRepresentationTruncated(2, 3, q, R)
            sage: B1 is B2 and B2 is B3
            True
        """
        if isinstance(n, FockSpaceTruncated):
            return super(HighestWeightRepresentationTruncated, cls).__classcall__(cls, n)
        F = FockSpaceTruncated(n, k, q, base_ring)
        return super(HighestWeightRepresentationTruncated, cls).__classcall__(cls, F)

    def __init__(self, fock_space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F = FockSpace(2, truncated=3)
            sage: B = F.highest_weight_representation()
            sage: x = B.a_realization() # This should be removed, but is needed so the test suite correctly obtains an element
            sage: TestSuite(B).run()
        """
        self._fock = fock_space
        R = fock_space.base_ring()
        Parent.__init__(self, base=R, category=ModulesWithBasis(R).WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: F.highest_weight_representation()
            Highest weight representation of ['A', 1, 1] of weight Lambda[0]
        """
        from sage.combinat.root_system.cartan_type import CartanType
        ct = CartanType(['A', self._fock._n - 1, 1])
        P = ct.root_system().weight_lattice()
        return "Highest weight representation of {} of weight {} truncated at {}".format(ct,
                    P.monomial(0), self._fock._k)

    def a_realization(self):
        """
        Return a realization of ``self``.

        EXAMPLES::

            sage: F = FockSpace(2, truncated=3)
            sage: B = F.highest_weight_representation()
            sage: B.a_realization()
            Highest weight representation of ['A', 1, 1] of weight Lambda[0]
             truncated at 3 in the lower global crystal basis
        """
        return self.G()

    class A(CombinatorialFreeModule, BindableClass):
        """
        The `A` basis of the Fock space which is the approximation of the
        lower global crystal basis.

        INPUT:

        - ``algorithm`` -- (default ``'GW'``) the algorithm to use when
          computing this basis in the Fock space; the possible values are:

          * ``'GW'`` -- use the algorithm given by Goodman and Wenzl
            in [GW1998]_
          * ``'LLT'`` -- use the LLT algorithm given in [LLT1996]_

        .. NOTE::

            The bases produced by the two algorithms are not the same
            in general.

        EXAMPLES::

            sage: F = FockSpace(5, truncated=4)
            sage: B = F.highest_weight_representation()
            sage: A = B.A()

        We demonstrate that they are different bases, but both algorithms
        still compute the basis `G`::

            sage: A2 = B.A('LLT')
            sage: G = B.G()
            sage: F(A[12,9])
            |12, 9> + q*|12, 4, 4, 1> + q*|8, 8, 5> + (q^2+1)*|8, 8, 4, 1>
            sage: F(A2[12,9])
            |12, 9> + q*|12, 4, 4, 1> + q*|8, 8, 5> + (q^2+2)*|8, 8, 4, 1>
            sage: G._G_to_fock_basis(Partition([12,9]), 'GW')
            |12, 9> + q*|12, 4, 4, 1> + q*|8, 8, 5> + q^2*|8, 8, 4, 1>
            sage: G._G_to_fock_basis(Partition([12,9]), 'LLT')
            |12, 9> + q*|12, 4, 4, 1> + q*|8, 8, 5> + q^2*|8, 8, 4, 1>
        """
        def __init__(self, basic, algorithm='GW'):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: F = FockSpace(2, truncated=3)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: TestSuite(A).run()
                sage: A2 = B.A('LLT')
                sage: TestSuite(A2).run()
            """
            self._basis_name = "approximation"
            if algorithm not in ['GW', 'LLT']:
                raise ValueError("invalid algorithm")
            self._alg = algorithm
            indices = RegularPartitions_truncated(basic._fock._n, basic._fock._k)
            CombinatorialFreeModule.__init__(self, basic.base_ring(), indices,
                                             prefix='A', bracket=False,
                                             monomial_cmp=lambda x,y: -cmp(x,y),
                                             category=HighestWeightRepresentationBases(basic))
            self.module_morphism(self._A_to_fock_basis,
                                 triangular='upper', unitriangular=True,
                                 codomain=basic._fock).register_as_coercion()

        @cached_method
        def LLT(self, la):
            """
            Return the result from the regular LLT algorithm on the partition
            ``la`` to compute the `A`-basis element in the corresponding
            Fock space.

            EXAMPLES::

                sage: F = FockSpace(5, truncated=4)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: F(A[12,9])
                |12, 9> + q*|12, 4, 4, 1> + q*|8, 8, 5> + (q^2+1)*|8, 8, 4, 1>
                sage: A.LLT(Partition([12,9]))
                |12, 9> + q*|12, 4, 4, 1> + q*|8, 8, 5> + (q^2+2)*|8, 8, 4, 1>
            """
            fock = self.realization_of()._fock
            k = fock._k
            I = fock._indices

            cur = fock.highest_weight_vector()
            # Get the ladders and apply it to the current element
            corners = la.corners()
            cells = set(la.cells())
            q = fock._q
            k = fock._n - 1 # This is sl_{k+1}
            r = fock._r[0]
            b = ZZ.zero()
            while any(c[1]*k + c[0] >= b for c in corners): # While there is some cell left to count
                power = 0
                i = -b + r # This will be converted to a mod n number
                for x in range(0, b // k + 1):
                    if (b-x*k, x) in cells:
                        power += 1
                        cur = cur.f(i)
                cur /= q_factorial(power, q)
                b += 1
            return cur

        def _skew_tableau(self, cur, nu, d):
            """
            Return the action of the skew tableaux formed from ``nu`` by
            applying the ``d`` fundamental weights.

            EXAMPLES::

                sage: F = FockSpace(3, truncated=3)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: G = B.G()
                sage: sk = A._skew_tableau(F(G[6,3]), Partition([6,3]), [1,1]); sk
                |8, 4> + q*|8, 2, 2> + q*|7, 4, 1> + q^2*|7, 3, 2> + q*|6, 6>
                 + q^2*|6, 5, 1> + q^2*|5, 4, 3> + q^3*|4, 4, 4>
                sage: sk == F(A[8,4])
                True
            """
            fock = self.realization_of()._fock
            last = None
            power = 0
            q = fock._q
            for i in reversed(range(len(d))):
                for dummy in range(d[i]):
                    for j in range(i+1):
                        col = nu[j] if j < len(nu) else 0
                        res = nu.content(j, col, multicharge=fock._r)
                        if res != last:
                            cur /= q_factorial(power, q)
                            power = 1
                            last = res
                        else:
                            power += 1
                        cur = cur.f(res)
                        nu = nu.add_cell(j, col)
            return cur / q_factorial(power, q)

        @cached_method
        def _A_to_fock_basis(self, la):
            """
            Return the `A` basis indexed by ``la`` in the natural basis.

            EXAMPLES::

                sage: F = FockSpace(2, truncated=3)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: A._A_to_fock_basis(Partition([2,1]))
                |2, 1>
                sage: F(A[3,1])
                |3, 1> + q*|2, 2> + q^2*|2, 1, 1>
                sage: A._A_to_fock_basis(Partition([3]))
                |3> + q*|1, 1, 1>

                sage: F = FockSpace(2, truncated=5)
                sage: B = F.highest_weight_representation()
                sage: A = B.A()
                sage: F(A[7,4,3])
                |7, 4, 3> + q*|7, 4, 1, 1, 1> + q^2*|7, 2, 2, 2, 1> + |5, 4, 3, 2>
                 + 2*q*|5, 4, 3, 1, 1> + 2*q^2*|5, 4, 2, 2, 1>
                 + 2*q^3*|5, 3, 3, 2, 1> + q^4*|4, 4, 3, 2, 1>
                sage: F(A[7,4,3,2])
                |7, 4, 3, 2> + q*|7, 4, 3, 1, 1> + q^2*|7, 4, 2, 2, 1>
                 + q^3*|7, 3, 3, 2, 1> + q^4*|6, 4, 3, 2, 1>
            """
            if self._alg == 'LLT':
                return self.LLT(la)
            fock = self.realization_of()._fock
            k = fock._k

            # We do some special cases first
            # For the empty partition
            if la.size() == 0:
                return fock.highest_weight_vector()

            # For length k partitions
            if len(la) == k:
                G = self.realization_of().G()
                return G._G_to_fock_basis(la)

            # For critical partitions
            n = fock._n
            if len(la) == k-1 and all((la[i] - la[i+1] + 1) % n == 0 for i in range(k-2)) \
                    and (la[-1] + 1) % n == 0:
                return fock.monomial(la)

            # For interior partitions
            shifted = [la[i] - (n - 1)*(k - 1 - i) for i in range(len(la))]
            if len(la) == k - 1 and shifted in _Partitions:
                # Construct the d's and the critical partition
                d = [(la[i] - la[i+1] + 1) % n for i in range(len(la)-1)]
                d.append((la[-1] + 1) % n)
                crit = list(la)
                for i,d_i in enumerate(d):
                    for j in range(i+1):
                        crit[j] -= d_i
                nu = fock._indices(crit)
                return self._skew_tableau(fock.monomial(nu), nu, d)

            # For non-interior partitions
            # Construct the d's and the partition ``a``
            a = list(la) + [0]*(k - 1 - len(la)) # Add 0's to get the correct length
            a = [a[i] + (k - 1 - i) for i in range(k-1)] # Shift the diagram
            #shifted = list(a) # Make a copy of the shifted partition in case we need it later
            d = [(a[i] - a[i+1]) % n for i in range(k-2)]
            d.append(a[-1] % n)
            for i,d_i in enumerate(d):
                for j in range(i+1):
                    a[j] -= d_i
            if sum(a) == 0: # a is contained in the fundamental box
                return self.LLT(la)

            p = list(a) # Make a copy that we can change
            for i in range(k-2):
                if a[i] - a[i+1] == 0:
                    d[i] -= 1
                    for j in range(i+1):
                        p[j] += 1
            if a[-1] == 0:
                d[-1] -= 1
                for j in range(k-1):
                    p[j] += 1
            p = [p[i] - (k - 1 - i) for i in range(k-1)]
            I = self._indices
            nu = I(p)

            is_regular = lambda mu: max(mu.to_exp()) < fock._n
            while not is_regular(nu):
                # It is not regular, so we need to find the first regular
                #   partition after adding columns
                while d[-1] == 0:
                    d.pop()
                for j in range(d[-1]):
                    p[j] += 1
                nu = I(p)

            if la == nu:
                j = -1
                for i in range(k-2):
                    if p[i] - p[i+1] == 0:
                        j = -2
                        break
                    if p[i] > n and p[i] - p[i+1] > n:
                        j = i
                if j != -2 and p[-1] > n:
                    j = k - 1
                if j < 0:
                    return self.LLT(la)

                G = self.realization_of().G()
                nu = I([p[i] - n if i <= j else p[i] for i in range(k-1)])
                d = [0]*j + [n]
                return self._skew_tableau(G._G_to_fock_basis(nu), nu, d)

            G = self.realization_of().G()
            return self._skew_tableau(G._G_to_fock_basis(nu), nu, d)

    approximation = A

    class G(CombinatorialFreeModule, BindableClass):
        """
        The lower global crystal basis living inside of Fock space.

        EXAMPLES::

            sage: F = FockSpace(4, truncated=2)
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: F(G[3,1])
            |3, 1>
            sage: F(G[6,2])
            |6, 2> + q*|5, 3>
            sage: F(G[14])
            |14> + q*|11, 3>

            sage: F = FockSpace(3, truncated=4)
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: F(G[4,1])
            |4, 1> + q*|3, 2>
            sage: F(G[4,2,2])
            |4, 2, 2> + q*|3, 2, 2, 1>

        We check agianst the tables in [LLT1996]_ (after truncating)::

            sage: F = FockSpace(3, truncated=3)
            sage: B = F.highest_weight_representation()
            sage: G = B.G()
            sage: F(G[10])
            |10> + q*|8, 2> + q*|7, 2, 1>
            sage: F(G[6,4])
            |6, 4> + q*|6, 2, 2> + q^2*|4, 4, 2>
            sage: F(G[5,5])
            |5, 5> + q*|4, 3, 3>

            sage: F = FockSpace(4, truncated=3)
            sage: G = F.highest_weight_representation().G()
            sage: F(G[3,3,1])
            |3, 3, 1>
            sage: F(G[3,2,2])
            |3, 2, 2>
            sage: F(G[7])
            |7> + q*|3, 3, 1>
        """
        def __init__(self, basic):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: F = FockSpace(2, truncated=3)
                sage: G = F.highest_weight_representation().G()
                sage: TestSuite(G).run()
                sage: F = FockSpace(4, truncated=3)
                sage: G = F.highest_weight_representation().G()
                sage: TestSuite(G).run()
            """
            self._basis_name = "lower global crystal"
            indices = RegularPartitions_truncated(basic._fock._n, basic._fock._k)
            CombinatorialFreeModule.__init__(self, basic.base_ring(), indices,
                                             prefix='G', bracket=False,
                                             monomial_cmp=lambda x,y: -cmp(x,y),
                                             category=HighestWeightRepresentationBases(basic))
            self.module_morphism(self._G_to_fock_basis,
                                 triangular='upper', unitriangular=True,
                                 codomain=basic._fock).register_as_coercion()

        @cached_method
        def _G_to_fock_basis(self, la, algorithm='GW'):
            """
            Return the `G` basis indexed by ``la`` in the natural basis.

            EXAMPLES::

                sage: F = FockSpace(3, truncated=3)
                sage: B = F.highest_weight_representation()
                sage: G = B.G()
                sage: G._G_to_fock_basis(Partition([3]))
                |3> + q*|2, 1>
                sage: G._G_to_fock_basis(Partition([2,1]))
                |2, 1> + q*|1, 1, 1>
                sage: G._G_to_fock_basis(Partition([2,1]), 'LLT')
                |2, 1> + q*|1, 1, 1>
            """
            fock = self.realization_of()._fock

            # Special cases:
            # For the empty partition
            if la.size() == 0:
                return fock.highest_weight_vector()

            # For length k partitions
            if algorithm == 'GW':
                n = fock._n
                k = fock._k
                if len(la) == k:
                    I = self._indices
                    x = la[-1]
                    mu = Partition([p - x for p in la])
                    add_cols = lambda nu: I([ v + x for v in list(nu) + [0]*(k - len(nu)) ])
                    return fock.sum_of_terms((add_cols(nu), c) for nu,c in self._G_to_fock_basis(mu))

                # For critical partitions
                n = fock._n
                if len(la) == k-1 and all((la[i] - la[i+1] + 1) % n == 0 for i in range(k-2)) \
                        and (la[-1] + 1) % n == 0:
                    return fock.monomial(la)

            # Perform the triangular reduction
            cur = self.realization_of().A(algorithm)._A_to_fock_basis(la)
            dom_cmp = lambda x,y: -1 if y.dominates(x) else 1
            s = cur.support()
            s.sort(cmp=dom_cmp) # Sort via dominance order
            s.pop() # Remove the largest

            q = fock._q
            while len(s) != 0:
                mu = s.pop()
                d = cur[mu].denominator()
                k = d.degree()
                n = cur[mu].numerator()
                if k != 0 or n.constant_coefficient() != 0:
                    gamma = sum(n[i] * (q**(i-k) + q**(k-i))
                                for i in range(min(n.degree(), k)))
                    gamma += n[k]
                    cur -= gamma * self._G_to_fock_basis(mu, algorithm)

                    # Add any new support elements
                    for x in cur.support():
                        if x == mu or not mu.dominates(x): # Add only things (strictly) dominated by mu
                            continue
                        for i in reversed(range(len(s))):
                            if not s[i].dominates(x):
                                s.insert(i+1, x)
                                break
            return cur

    lower_global_crystal_basis = G


"""
Generic code for bases

This is a collection of code that is shared by bases of noncommutative
symmetric functions and quasisymmetric functions.

AUTHORS:

- Jason Bandlow
- Franco Saliola
- Chris Berg
"""
#*****************************************************************************
#       Copyright (C) 2010 Jason Bandlow <jbandlow@gmail.com>,
#                     2012 Franco Saliola <saliola@gmail.com>,
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

from sage.misc.cachefunc import cached_method
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.modules_with_basis import ModulesWithBasis, ModuleMorphismByLinearity
from sage.combinat.composition import Composition
from sage.combinat.partition import Partition
from sage.combinat.permutation import Permutations
from sage.rings.integer import Integer
from sage.categories.all import AlgebrasWithBasis
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.abstract_method import abstract_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.realizations import RealizationsCategory

class BasesOfQSymOrNCSF(Category_realization_of_parent):

    def _repr_(self):
        r"""
        String representation of this category

        TESTS::

            sage: R = NonCommutativeSymmetricFunctions(ZZ).R()
            sage: C = R.category().super_categories()[0]
            sage: C._repr_()
            'Category of bases of Non-Commutative Symmetric Functions or Quasisymmetric functions over the Integer Ring'

        """
        return "Category of bases of Non-Commutative Symmetric Functions or Quasisymmetric functions over the %s" % self.base().base_ring()

    def super_categories(self):
        r"""
        TESTS::

            sage: from sage.combinat.ncsf_qsym.generic_basis_code import BasesOfQSymOrNCSF
            sage: QSym = QuasiSymmetricFunctions(QQ)
            sage: BasesOfQSymOrNCSF(QSym).super_categories()
            [Category of realizations of Quasisymmetric functions over the Rational Field,
             Category of graded hopf algebras with basis over Rational Field,
             Join of Category of graded hopf algebras over Rational Field and Category of realizations of hopf algebras over Rational Field]

        """
        R = self.base().base_ring()
        from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
        from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
        return [self.base().Realizations(),
                GradedHopfAlgebrasWithBasis(R),
                GradedHopfAlgebras(R).Realizations()]

    class ParentMethods:

        def _repr_(self):
            """
            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: S._repr_()
                'Non-Commutative Symmetric Functions over the Rational Field in the Complete basis'
                sage: F = QuasiSymmetricFunctions(ZZ).Fundamental()
                sage: F._repr_()
                'Quasisymmetric functions over the Integer Ring in the Fundamental basis'
            """
            return "%s in the %s basis" % (self.realization_of(), self._realization_name())

        def __getitem__(self, c, *rest):
            """
            This method implements the abuses of notations::

                sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
                sage: Psi[2,1]
                Psi[2, 1]
                sage: Psi[[2,1]]
                Psi[2, 1]
                sage: Psi[Composition([2,1])]
                Psi[2, 1]

            .. todo::

                This should call ``super.monomial`` if the input can't
                be made into a composition so as not to interfere with
                the standard notation ``Psi['x,y,z']``.

                This could possibly be shared with Sym, FQSym, and
                other algebras with bases indexed by list-like objects
            """
            if isinstance(c, Composition):
                assert len(rest) == 0
            else:
                if len(rest) > 0 or isinstance(c, (int, Integer)):
                    c = Composition([c] + list(rest))
                else:
                    c = Composition(list(c))
            return self.monomial(c)

        # could go to Algebras(...).Graded().Connected() or Modules(...).Graded().Connected()
        @cached_method
        def one_basis(self):
            r"""
            Return the empty composition.

            OUTPUT

            - The empty composition.

            EXAMPLES::

                sage: L=NonCommutativeSymmetricFunctions(QQ).L()
                sage: parent(L)
                <class 'sage.combinat.ncsf_qsym.ncsf.NonCommutativeSymmetricFunctions.Elementary_with_category'>
                sage: parent(L).one_basis()
                []
            """
            return Composition([])

        # Combinatorial rules

        def sum_of_finer_compositions(self, composition):
            r"""
            Return the sum of all finer compositions.

            INPUT:

            - ``composition`` -- a composition

            OUTPUT:

            - The sum of all basis ``self`` elements which are indexed by
              compositions finer than ``composition``.

            EXAMPLES::

                sage: L=NonCommutativeSymmetricFunctions(QQ).L()
                sage: L.sum_of_finer_compositions(Composition([2,1]))
                L[1, 1, 1] + L[2, 1]
                sage: R=NonCommutativeSymmetricFunctions(QQ).R()
                sage: R.sum_of_finer_compositions(Composition([1,3]))
                R[1, 1, 1, 1] + R[1, 1, 2] + R[1, 2, 1] + R[1, 3]
            """

            return self.sum_of_monomials( compo for compo in composition.finer() )

        def sum_of_fatter_compositions(self, composition):
            r"""
            Return the sum of all fatter compositions.

            INPUT:

            - ``composition`` -- a composition

            OUTPUT:

            - the sum of all basis elements which are indexed by
              compositions fatter (coarser?) than ``composition``.

            EXAMPLES::

                sage: L=NonCommutativeSymmetricFunctions(QQ).L()
                sage: L.sum_of_fatter_compositions(Composition([2,1]))
                L[2, 1] + L[3]
                sage: R=NonCommutativeSymmetricFunctions(QQ).R()
                sage: R.sum_of_fatter_compositions(Composition([1,3]))
                R[1, 3] + R[4]
            """
            return self.sum_of_monomials( compo for compo in composition.fatter() )

        def alternating_sum_of_finer_compositions(self, composition, conjugate = False):
            """
            Return the alternating sum of finer compositions in a basis of the
            non-commutative symmetric functions.

            INPUT:

            - ``composition`` -- a composition
            - ``conjugate`` -- (default: ``False``) a boolean

            OUTPUT:

            - The alternating sum of the compositions finer than ``composition``,
              in the basis ``self``. The alternation is upon the length of the
              compositions, and is normalized so that ``composition`` has
              coefficient `1`. If the variable ``conjugate`` is set to ``True``,
              then the conjugate of ``composition`` is used instead of
              ``composition``.

            EXAMPLES::

                sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
                sage: elementary = NCSF.elementary()
                sage: elementary.alternating_sum_of_finer_compositions(Composition([2,2,1]))
                L[1, 1, 1, 1, 1] - L[1, 1, 2, 1] - L[2, 1, 1, 1] + L[2, 2, 1]
                sage: elementary.alternating_sum_of_finer_compositions(Composition([1,2]))
                -L[1, 1, 1] + L[1, 2]

            TESTS::

                sage: complete = NonCommutativeSymmetricFunctions(ZZ).complete()
                sage: I = Composition([2])
                sage: x = complete.alternating_sum_of_finer_compositions(I)
                sage: [c.parent() for c in x.coefficients()]
                [Integer Ring, Integer Ring]
            """
            if conjugate:
                composition = composition.conjugate()
            l = len(composition)
            ring = self.base_ring()
            return self.sum_of_terms( (compo, ring((-1)**(len(compo)-l))) for compo in composition.finer() )

        def alternating_sum_of_fatter_compositions(self, composition):
            """
            Return the alternating sum of fatter compositions in a basis of the
            non-commutative symmetric functions.

            INPUT:

            - ``composition`` -- a composition

            OUTPUT:

            - The alternating sum of the compositions fatter than ``composition``,
              in the basis ``self``. The alternation is upon the length of the
              compositions, and is normalized so that ``composition`` has
              coefficient `1`.

            EXAMPLES::

                sage: NCSF=NonCommutativeSymmetricFunctions(QQ)
                sage: elementary = NCSF.elementary()
                sage: elementary.alternating_sum_of_fatter_compositions(Composition([2,2,1]))
                L[2, 2, 1] - L[2, 3] - L[4, 1] + L[5]
                sage: elementary.alternating_sum_of_fatter_compositions(Composition([1,2]))
                L[1, 2] - L[3]

            TESTS::

                sage: complete = NonCommutativeSymmetricFunctions(ZZ).complete()
                sage: I = Composition([1,1])
                sage: x = complete.alternating_sum_of_fatter_compositions(I)
                sage: [c.parent() for c in x.coefficients()]
                [Integer Ring, Integer Ring]
            """
            l = len(composition)
            ring = self.base_ring()
            return self.sum_of_terms( (compo, ring((-1)**(len(compo)-l))) for compo in composition.fatter() )

        def sum_of_partition_rearrangements(self, par):
            """
            Return the sum of all basis elements indexed by compositions which can be
            sorted to obtain a given partition.

            INPUT:

            - ``par`` -- a partition

            OUTPUT:

            - The sum of all ``self`` basis elements indexed by compositions
              which are permutations of ``par`` (without multiplicity).

            EXAMPLES::

                sage: NCSF=NonCommutativeSymmetricFunctions(QQ)
                sage: elementary = NCSF.elementary()
                sage: elementary.sum_of_partition_rearrangements(Partition([2,2,1]))
                L[1, 2, 2] + L[2, 1, 2] + L[2, 2, 1]
                sage: elementary.sum_of_partition_rearrangements(Partition([3,2,1]))
                L[1, 2, 3] + L[1, 3, 2] + L[2, 1, 3] + L[2, 3, 1] + L[3, 1, 2] + L[3, 2, 1]
                sage: elementary.sum_of_partition_rearrangements(Partition([]))
                L[]
            """
            return self.sum_of_monomials( Composition(comp) for comp in Permutations(par) )

        def _comp_to_par(self, comp):
            """
            Return the partition if the composition is actually a partition. Otherwise
            returns nothing.

            INPUT:

            - ``comp`` -- a composition

            OUTPUT:

            - ``comp`` as a partition, if it is sorted; otherwise returns
              ``None`` (nothing).

            EXAMPLES::

                sage: NCSF=NonCommutativeSymmetricFunctions(QQ)
                sage: L = NCSF.elementary()
                sage: L._comp_to_par(Composition([1,1,3,1,2]))
                sage: L.sum_of_partition_rearrangements(Composition([]))
                L[]
                sage: L._comp_to_par(Composition([3,2,1,1]))
                [3, 2, 1, 1]
            """
            try:
                return Partition(comp)
            except ValueError:
                return None

        def degree_on_basis(self, I):
            r"""
            Return the degree of the basis element indexed by `I`.

            INPUT:

            - ``I`` -- a composition

            OUTPUT:

            - The degree of the non-commutative symmetric function basis
              element of ``self`` indexed by ``I``. By definition, this is
              the size of the composition ``I``.

            EXAMPLES::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: R.degree_on_basis(Composition([2,3]))
                5
                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: M.degree_on_basis(Composition([3,2]))
                5
                sage: M.degree_on_basis(Composition([]))
                0
            """
            return I.size()

        def skew(self, x, y, side='left'):
            r"""
            Return a function ``x`` in ``self`` skewed by a function
            ``y`` in the Hopf dual of ``self``.

            INPUT:

            - ``x`` -- a non-commutative or quasi-symmetric function; it is
              an element of ``self``
            - ``y`` -- a quasi-symmetric or non-commutative symmetric
              function; it is an element of the dual algebra of ``self``
            - ``side`` -- (default: ``'left'``)
              either ``'left'`` or ``'right'``

            OUTPUT:

            - The result of skewing the element ``x`` by the Hopf algebra
              element ``y`` (either from the left or from the right, as
              determined by ``side``), written in the basis ``self``.

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: S.skew(S[2,2,2], F[1,1])
                S[1, 1, 2] + S[1, 2, 1] + S[2, 1, 1]
                sage: S.skew(S[2,2,2], F[2])
                S[1, 1, 2] + S[1, 2, 1] + S[2, 1, 1] + 3*S[2, 2]

            ::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: R.skew(R[2,2,2], F[1,1])
                R[1, 1, 2] + R[1, 2, 1] + R[1, 3] + R[2, 1, 1] + 2*R[2, 2] + R[3, 1] + R[4]
                sage: R.skew(R[2,2,2], F[2])
                R[1, 1, 2] + R[1, 2, 1] + R[1, 3] + R[2, 1, 1] + 3*R[2, 2] + R[3, 1] + R[4]

            ::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M.skew(M[3,2], S[2])
                0
                sage: M.skew(M[3,2], S[2], side='right')
                M[3]
                sage: M.skew(M[3,2], S[3])
                M[2]
                sage: M.skew(M[3,2], S[3], side='right')
                0

            TESTS::

                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: R.skew([2,1], [1])
                Traceback (most recent call last):
                ...
                AssertionError: x must be an element of Non-Commutative Symmetric Functions over the Rational Field
                sage: R([2,1]).skew_by([1])
                Traceback (most recent call last):
                ...
                AssertionError: y must be an element of Quasisymmetric functions over the Rational Field
                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F([2,1]).skew_by([1])
                Traceback (most recent call last):
                ...
                AssertionError: y must be an element of Non-Commutative Symmetric Functions over the Rational Field
            """
            alg = self.realization_of()
            assert x in alg, "x must be an element of %s" % alg
            assert y in alg.dual(), "y must be an element of %s" % alg.dual()
            if hasattr(self, 'dual'):
                x = self(x)
                y = self.dual()(y)
                v = 1 if side == 'left' else 0
                return self.sum(coeff * y[IJ[1-v]] * self[IJ[v]] \
                                for (IJ, coeff) in x.coproduct() if IJ[1-v] in y)
            else:
                return self._skew_by_coercion(x, y, side=side)

        def _skew_by_coercion(self, x, y, side='left'):
            r"""
            Return a function ``x`` in ``self`` skewed by a function
            ``y`` in the Hopf dual of ``self`` using coercion.

            INPUT:

            - ``x`` -- a non-commutative or quasi-symmetric function; it is
              an element of ``self``
            - ``y`` -- a quasi-symmetric or non-commutative symmetric
              function; it is an element of the dual algebra of ``self``
            - ``side`` -- (default: ``'left'``)
              either ``'left'`` or ``'right'``

            OUTPUT:

            - The result of skewing the element ``x`` by the Hopf algebra
              element ``y`` (either from the left or from the right, as
              determined by ``side``), written in the basis ``self``.
              This uses coercion to a concreate realization (either the
              complete basis of non-commutative symmetric functions or
              the monomial basis of the quasi-symmetric functions).

            EXAMPLES::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M._skew_by_coercion(M[1,2,1,3], R[1])
                M[2, 1, 3]
                sage: M._skew_by_coercion(M[1,2,1,3], R[1],side='right')
                0
            """
            a_realization = self.realization_of().a_realization()
            return self(a_realization.skew(a_realization(x), y, side=side))

        def duality_pairing(self, x, y):
            r"""
            The duality pairing between elements of `NSym` and elements
            of `QSym`.

            This is a default implementation that uses
            ``self.realizations_of().a_realization()`` and its dual basis.

            INPUT:

            - ``x`` -- an element of ``self``
            - ``y`` -- an element in the dual basis of ``self``

            OUTPUT:

            - The result of pairing the function ``x`` from ``self`` with the function
              ``y`` from the dual basis of ``self``

            EXAMPLES::

                sage: R = NonCommutativeSymmetricFunctions(QQ).Ribbon()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: R.duality_pairing(R[1,1,2], F[1,1,2])
                1
                sage: R.duality_pairing(R[1,2,1], F[1,1,2])
                0
                sage: F.duality_pairing(F[1,2,1], R[1,1,2])
                0

            ::

                sage: S = NonCommutativeSymmetricFunctions(QQ).Complete()
                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: S.duality_pairing(S[1,1,2], M[1,1,2])
                1
                sage: S.duality_pairing(S[1,2,1], M[1,1,2])
                0
                sage: M.duality_pairing(M[1,1,2], S[1,1,2])
                1
                sage: M.duality_pairing(M[1,2,1], S[1,1,2])
                0

            ::

                sage: S = NonCommutativeSymmetricFunctions(QQ).Complete()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: S.duality_pairing(S[1,2], F[1,1,1])
                0
                sage: S.duality_pairing(S[1,1,1,1], F[4])
                1
            """
            if hasattr(self, 'dual'):
                x = self(x)
                y = self.dual()(y)
                return sum(coeff * y[I] for (I, coeff) in x)
            else:
                return self.duality_pairing_by_coercion(x, y)

        def duality_pairing_by_coercion(self, x, y):
            r"""
            The duality pairing between elements of NSym and elements of QSym.

            This is a default implementation that uses
            ``self.realizations_of().a_realization()`` and its dual basis.

            INPUT:

            - ``x`` -- an element of ``self``
            - ``y`` -- an element in the dual basis of ``self``

            OUTPUT:

            - The result of pairing the function ``x`` from ``self`` with the function
              ``y`` from the dual basis of ``self``

            EXAMPLES::

                sage: L = NonCommutativeSymmetricFunctions(QQ).Elementary()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: L.duality_pairing_by_coercion(L[1,2], F[1,2])
                0
                sage: F.duality_pairing_by_coercion(F[1,2], L[1,2])
                0
                sage: L.duality_pairing_by_coercion(L[1,1,1], F[1,2])
                1
                sage: F.duality_pairing_by_coercion(F[1,2], L[1,1,1])
                1
            """
            a_realization = self.realization_of().a_realization()
            x = a_realization(x)
            y = a_realization.dual()(y)
            return sum(coeff * y[I] for (I, coeff) in x)

        def duality_pairing_matrix(self, basis, degree):
            r"""
            The matrix of scalar products between elements of NSym and
            elements of QSym.

            INPUT:

            - ``basis`` -- A basis of the dual Hopf algebra
            - ``degree`` -- a non-negative integer

            OUTPUT:

            - The matrix of scalar products between the basis ``self`` and the basis
              ``basis`` in the dual Hopf algebra of degree ``degree``.

            EXAMPLES:

            The ribbon basis of NCSF is dual to the fundamental basis of
            QSym::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: R.duality_pairing_matrix(F, 3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]
                sage: F.duality_pairing_matrix(R, 3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]

            The complete basis of NCSF is dual to the monomial basis of
            QSym::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: S.duality_pairing_matrix(M, 3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]
                sage: M.duality_pairing_matrix(S, 3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]

            The matrix between the ribbon basis of NCSF and the monomial
            basis of QSym::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: M = QuasiSymmetricFunctions(QQ).Monomial()
                sage: R.duality_pairing_matrix(M, 3)
                [ 1 -1 -1  1]
                [ 0  1  0 -1]
                [ 0  0  1 -1]
                [ 0  0  0  1]
                sage: M.duality_pairing_matrix(R, 3)
                [ 1  0  0  0]
                [-1  1  0  0]
                [-1  0  1  0]
                [ 1 -1 -1  1]

            The matrix between the complete basis of NCSF and the
            fundamental basis of QSym::

                sage: S = NonCommutativeSymmetricFunctions(QQ).complete()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: S.duality_pairing_matrix(F, 3)
                [1 1 1 1]
                [0 1 0 1]
                [0 0 1 1]
                [0 0 0 1]

            A base case test::

                sage: R.duality_pairing_matrix(M,0)
                [1]
            """
            from sage.matrix.constructor import matrix
            # TODO: generalize to keys indexing the basis of the graded component
            from sage.combinat.composition import Compositions
            return matrix(self.base_ring(),
                    [[self.duality_pairing(self[I], basis[J]) \
                            for J in Compositions(degree)] \
                            for I in Compositions(degree)])

        def counit_on_basis(self, I):
            r"""
            The counit is defined by sending all elements of positive degree to zero.

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.counit_on_basis([1,3])
                0
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M.counit_on_basis([1,3])
                0

            TESTS::

                sage: S.counit_on_basis([])
                1
                sage: S.counit_on_basis(Composition([]))
                1
                sage: M.counit_on_basis([])
                1
                sage: M.counit_on_basis(Composition([]))
                1
            """
            if I != []:
                return self.base_ring().zero()
            else:
                return self.base_ring().one()

    class ElementMethods:

        def duality_pairing(self, y):
            r"""
            The duality pairing between elements of `NSym` and elements
            of `QSym`.

            The complete basis is dual to the monomial basis with respect
            to this pairing.

            INPUT:

            - ``y`` -- an element of the dual Hopf algebra of ``self``

            OUTPUT:

            - The result of pairing ``self`` with ``y``.

            EXAMPLES::

                sage: R = NonCommutativeSymmetricFunctions(QQ).Ribbon()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: R[1,1,2].duality_pairing(F[1,1,2])
                1
                sage: R[1,2,1].duality_pairing(F[1,1,2])
                0

            ::

                sage: L = NonCommutativeSymmetricFunctions(QQ).Elementary()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: L[1,2].duality_pairing(F[1,2])
                0
                sage: L[1,1,1].duality_pairing(F[1,2])
                1

            """
            return self.parent().duality_pairing(self, y)

        def skew_by(self, y, side='left'):
            r"""
            The operation which is dual to multiplication by ``y``, where ``y``
            is an element of the dual space of ``self``.

            This is calculated through the coproduct of ``self`` and the
            expansion of ``y`` in the dual basis.

            INPUT:

            - ``y`` -- an element of the dual Hopf algebra of ``self``
            - ``side`` -- (Default='left') Either 'left' or 'right'

            OUTPUT:

            - The result of skewing ``self`` by ``y``, on the side ``side``

            EXAMPLES:

            Skewing an element of NCSF by an element of QSym::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: F = QuasiSymmetricFunctions(QQ).Fundamental()
                sage: R([2,2,2]).skew_by(F[1,1])
                R[1, 1, 2] + R[1, 2, 1] + R[1, 3] + R[2, 1, 1] + 2*R[2, 2] + R[3, 1] + R[4]
                sage: R([2,2,2]).skew_by(F[2])
                R[1, 1, 2] + R[1, 2, 1] + R[1, 3] + R[2, 1, 1] + 3*R[2, 2] + R[3, 1] + R[4]

            Skewing an element of QSym by an element of NCSF::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F[3,2].skew_by(R[1,1])
                0
                sage: F[3,2].skew_by(R[1,1], side='right')
                0
                sage: F[3,2].skew_by(S[1,1,1], side='right')
                F[2]
                sage: F[3,2].skew_by(S[1,2], side='right')
                F[2]
                sage: F[3,2].skew_by(S[2,1], side='right')
                0
                sage: F[3,2].skew_by(S[1,1,1])
                F[2]
                sage: F[3,2].skew_by(S[1,1])
                F[1, 2]
                sage: F[3,2].skew_by(S[1])
                F[2, 2]

            ::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: M = QuasiSymmetricFunctions(QQ).M()
                sage: M[3,2].skew_by(S[2])
                0
                sage: M[3,2].skew_by(S[2], side='right')
                M[3]
                sage: M[3,2].skew_by(S[3])
                M[2]
                sage: M[3,2].skew_by(S[3], side='right')
                0
            """
            return self.parent().skew(self, y, side=side)

        def degree(self):
            """
            The maximum of the degrees of the homogeneous summands.

            .. SEEALSO:: :meth:`~sage.categories.graded_algebras_with_basis.GradedAlgebrasWithBasis.ElementMethods.homogeneous_degree`

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.degree()
                2
                sage: (x^3 + 4*y^2).degree()
                6
                sage: ((1 + x)^3).degree()
                6

            ::

                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: (x, y) = (F[2], F[3])
                sage: x.degree()
                2
                sage: (x^3 + 4*y^2).degree()
                6
                sage: ((1 + x)^3).degree()
                6

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: The zero element does not have a well-defined degree.
                sage: F = QuasiSymmetricFunctions(QQ).F()
                sage: F.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: The zero element does not have a well-defined degree.
            """
            return self.maximal_degree()

class AlgebraMorphism(ModuleMorphismByLinearity): # Find a better name
    """
    A class for algebra morphism defined on a free algebra from the image of the generators
    """
    def __init__(self, domain, on_generators, position = 0, codomain = None, category = None, anti = False):
        """
        Given a map on the multiplicative basis of a free algebra, this method
        returns the algebra morphism that is the linear extension of its image
        on generators.

        INPUT:

        - ``domain`` -- an algebra with a multiplicative basis
        - ``on_generators`` -- a function defined on the index set of the generators
        - ``codomain`` -- the codomain
        - ``position`` -- integer; default is 0
        - ``category`` -- a category; defaults to None
        - ``anti`` -- a boolean; defaults to False

        OUTPUT:

        - module morphism

        EXAMPLES:

        We construct explicitly an algbera morphism::

            sage: from sage.combinat.ncsf_qsym.generic_basis_code import AlgebraMorphism
            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: Psi = NCSF.Psi()
            sage: f = AlgebraMorphism(Psi, attrcall('conjugate'), codomain=Psi)
            sage: f
            Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis

        Usually, however, one constructs algebra morphisms
        using the ``algebra_morphism`` method for an algebra::

            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: Psi = NCSF.Psi()
            sage: def double(i) : return Psi[i,i]
            sage: f = Psi.algebra_morphism(double, codomain = Psi)
            sage: f
            Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
            sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
            2*Psi[] + 3*Psi[1, 1, 3, 3, 2, 2] + Psi[2, 2, 4, 4]
            sage: f.category()
            Join of Category of hom sets in Category of modules with basis over Rational Field and Category of hom sets in Category of rings

        When extra properties about the morphism are known, one
        can specify the category of which it is a morphism::

            sage: def negate(i): return -Psi[i]
            sage: f = Psi.algebra_morphism(negate, codomain = Psi, category = GradedHopfAlgebrasWithBasis(QQ))
            sage: f
            Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
            sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
            2*Psi[] - 3*Psi[1, 3, 2] + Psi[2, 4]
            sage: f.category()
            Join of Category of hom sets in Category of modules with basis over Rational Field and Category of hom sets in Category of rings

        If ``anti`` is true, this returns an anti-algebra morphism::

            sage: f = Psi.algebra_morphism(double, codomain = Psi, anti=True)
            sage: f
            Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Psi basis
            sage: f(2*Psi[[]] + 3 * Psi[1,3,2] + Psi[2,4] )
            2*Psi[] + 3*Psi[2, 2, 3, 3, 1, 1] + Psi[4, 4, 2, 2]
            sage: f.category()
            Category of hom sets in Category of modules with basis over Rational Field


        TESTS::

            sage: Psi = NonCommutativeSymmetricFunctions(QQ).Psi()
            sage: Phi = NonCommutativeSymmetricFunctions(QQ).Phi()
            sage: f = Psi.algebra_morphism(Phi.antipode_on_generators, codomain=Phi)
            sage: f(Psi[1, 2, 2, 1])
            Phi[1, 2, 2, 1]
            sage: f(Psi[3, 1, 2])
            -Phi[3, 1, 2]
            sage: f.__class__
            <class 'sage.combinat.ncsf_qsym.generic_basis_code.AlgebraMorphism'>
            sage: TestSuite(f).run(skip=['_test_nonzero_equal']) # known issue; see ModuleMorphismByLinearity.__init__
            Failure in _test_category:
            ...
            The following tests failed: _test_category
        """
        assert position == 0
        assert codomain is not None
        if category is None:
            if anti:
                category = ModulesWithBasis (domain.base_ring())
            else:
                category = AlgebrasWithBasis(domain.base_ring())
        self._anti = anti
        self._on_generators = on_generators
        ModuleMorphismByLinearity.__init__(self, domain = domain, codomain = codomain, position = position, category = category)

    def _on_basis(self, c):
        r"""
        Computes the image of this morphism on the basis element indexed by
        ``c``.

        INPUT:

        - ``c`` -- an iterable that spits out generators

        OUTPUT:

        - element of the codomain

        EXAMPLES::

            sage: from sage.combinat.ncsf_qsym.generic_basis_code import AlgebraMorphism
            sage: NCSF = NonCommutativeSymmetricFunctions(QQ)
            sage: Psi = NCSF.Psi()
            sage: Phi = NCSF.Phi()
            sage: f = AlgebraMorphism(Psi, lambda i : Phi[i,i], codomain=Phi)
            sage: f._on_basis([ 3, 2 ])
            Phi[3, 3, 2, 2]

        """
        if self._anti:
            c = reversed(c)
        return self.codomain().prod(self._on_generators(i) for i in c)

class GradedModulesWithInternalProduct(Category_over_base_ring):
    r"""
    Constructs the class of modules with internal product. This is used to give an internal
    product structure to the non-commutative symmetric functions.

    EXAMPLES::

        sage: from sage.combinat.ncsf_qsym.generic_basis_code import GradedModulesWithInternalProduct
        sage: N = NonCommutativeSymmetricFunctions(QQ)
        sage: R = N.ribbon()
        sage: R in GradedModulesWithInternalProduct(QQ)
        True
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.combinat.ncsf_qsym.generic_basis_code import GradedModulesWithInternalProduct
            sage: GradedModulesWithInternalProduct(ZZ).super_categories()
            [Category of graded modules over Integer Ring]
        """
        from sage.categories.graded_modules import GradedModules
        R = self.base_ring()
        return [GradedModules(R)]

    class ParentMethods:
        @abstract_method(optional=True)
        def internal_product_on_basis(self, I, J):
            """
            The internal product of the two basis elements indexed by ``I`` and
            ``J`` (optional)

            INPUT:

             - ``I``, ``J`` -- compositions indexing two elements of the basis of self

            Returns the internal product of the corresponding basis elements.
            If this method is implemented, the internal product is defined from
            it by linearity.

            EXAMPLES::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: S = N.complete()
                sage: S.internal_product_on_basis([2,2], [1,2,1])
                2*S[1, 1, 1, 1] + S[1, 1, 2] + S[2, 1, 1]
                sage: S.internal_product_on_basis([2,2], [2,1])
                0
            """

        @lazy_attribute
        def internal_product(self):
            r"""
            Internal product as an endomorphism of ``self``.

            This is constructed by extending the method
            :meth:`internal_product_on_basis` bilinearly, if available,
            or using the method :meth:`internal_product_by_coercion`.

            OUTPUT:

            - The internal product map of the algebra the non-commutative
              symmetric functions.

            EXAMPLES::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: S = N.complete()
                sage: S.internal_product
                Generic endomorphism of Non-Commutative Symmetric Functions over the Rational Field in the Complete basis
                sage: S.internal_product(S[2,2], S[1,2,1])
                2*S[1, 1, 1, 1] + S[1, 1, 2] + S[2, 1, 1]
                sage: S.internal_product(S[2,2], S[1,2])
                0

            ::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: R = N.ribbon()
                sage: R.internal_product
                <bound method ....internal_product_by_coercion ...>
                sage: R.internal_product_by_coercion(R[1, 1], R[1,1])
                R[2]
                sage: R.internal_product(R[2,2], R[1,2])
                0

            .. TODO::

                Despite the ``__repr__``, this is NOT an endomorphism!
            """
            if self.internal_product_on_basis is not NotImplemented:
                return self.module_morphism(
                                self.module_morphism(self.internal_product_on_basis,
                                                     position=0,
                                                     codomain=self),
                                position=1)
            else:
                return self.internal_product_by_coercion

        itensor = internal_product
        kronecker_product = internal_product

    class ElementMethods:
        def internal_product(self, other):
            r"""
            Return the internal product of two non-commutative
            symmetric functions.

            The internal product on the algebra of non-commutative symmetric
            functions is adjoint to the internal coproduct on the algebra of
            quasisymmetric functions with respect to the duality pairing
            between these two algebras. This means, explicitly, that any
            two non-commutative symmetric functions `f` and `g` and any
            quasi-symmetric function `h` satisfy

            .. MATH::

                \langle f * g, h \rangle = \sum_i \langle f, h^{\prime}_i
                \rangle \langle g, h^{\prime\prime}_i \rangle,

            where we write `\Delta^{\times}(h)` as `\sum_i h^{\prime}_i
            \otimes h^{\prime\prime}_i`.

            Aliases for :meth:`internal_product()` are :meth:`itensor()` and
            :meth:`kronecker_product()`.

            INPUT:

            - ``other`` -- another non-commutative symmetric function

            OUTPUT:

            - The result of taking the internal product of ``self`` with
              ``other``.

            EXAMPLES::

                sage: N = NonCommutativeSymmetricFunctions(QQ)
                sage: S = N.complete()
                sage: x = S.an_element(); x
                2*S[] + 2*S[1] + 3*S[1, 1]
                sage: x.internal_product(S[2])
                3*S[1, 1]
                sage: x.internal_product(S[1])
                2*S[1]
                sage: S[1,2].internal_product(S[1,2])
                S[1, 1, 1] + S[1, 2]

            Let us check the duality between the inner product and the inner
            coproduct in degree `4`::

                sage: M = QuasiSymmetricFunctions(FiniteField(29)).M()
                sage: S = NonCommutativeSymmetricFunctions(FiniteField(29)).S()
                sage: def tensor_incopr(f, g, h):  # computes \sum_i \left< f, h'_i \right> \left< g, h''_i \right>
                ....:     result = h.base_ring().zero()
                ....:     h_parent = h.parent()
                ....:     for partition_pair, coeff in h.internal_coproduct().monomial_coefficients().items():
                ....:         result += coeff * f.duality_pairing(h_parent[partition_pair[0]]) * g.duality_pairing(h_parent[partition_pair[1]])
                ....:     return result
                sage: def testall(n):
                ....:     return all( all( all( tensor_incopr(S[u], S[v], M[w]) == (S[u].itensor(S[v])).duality_pairing(M[w])
                ....:                           for w in Compositions(n) )
                ....:                      for v in Compositions(n) )
                ....:                 for u in Compositions(n) )
                sage: testall(2)
                True
                sage: testall(3)  # long time
                True
                sage: testall(4)  # long time
                True

            The internal product on the algebra of non-commutative symmetric
            functions commutes with the canonical commutative projection on
            the symmetric functions::

                sage: S = NonCommutativeSymmetricFunctions(ZZ).S()
                sage: e = SymmetricFunctions(ZZ).e()
                sage: def int_pr_of_S_in_e(I, J):
                ....:     return (S[I].internal_product(S[J])).to_symmetric_function()
                sage: all( all( int_pr_of_S_in_e(I, J)
                ....:           == S[I].to_symmetric_function().internal_product(S[J].to_symmetric_function())
                ....:           for I in Compositions(3) )
                ....:      for J in Compositions(3) )
                True
            """
            return self.parent().internal_product(self, other)

        itensor = internal_product
        kronecker_product = internal_product

    class Realizations(RealizationsCategory):
        class ParentMethods:
            def internal_product_by_coercion(self, left, right):
                r"""
                Internal product of ``left`` and ``right``.

                This is a default implementation that computes
                the internal product in the realization specified
                by ``self.realization_of().a_realization()``.

                INPUT:

                - ``left`` -- an element of the non-commutative symmetric functions
                - ``right`` -- an element of the non-commutative symmetric functions

                OUTPUT:

                - The internal product of ``left`` and ``right``.

                EXAMPLES::

                    sage: S=NonCommutativeSymmetricFunctions(QQ).S()
                    sage: S.internal_product_by_coercion(S[2,1], S[3])
                    S[2, 1]
                    sage: S.internal_product_by_coercion(S[2,1], S[4])
                    0
                """
                R = self.realization_of().a_realization()
                return self(R.internal_product(R(left), R(right)))


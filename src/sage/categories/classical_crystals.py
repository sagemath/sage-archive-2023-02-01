r"""
Classical Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals

class ClassicalCrystals(Category_singleton):
    """
    The category of classical crystals, that is crystals of finite Cartan type.

    EXAMPLES::

        sage: C = ClassicalCrystals()
        sage: C
        Category of classical crystals
        sage: C.super_categories()
        [Category of regular crystals,
         Category of finite crystals,
         Category of highest weight crystals]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = ClassicalCrystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: ClassicalCrystals().super_categories()
            [Category of regular crystals,
             Category of finite crystals,
             Category of highest weight crystals]
        """
        return [RegularCrystals(), FiniteCrystals(), HighestWeightCrystals()]

    def example(self, n = 3):
        """
        Returns an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = ClassicalCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        return Crystals().example(n)


    class ParentMethods:

        @cached_method
        def opposition_automorphism(self):
            r"""
            Returns the opposition automorphism

            The *opposition automorphism* is the automorphism
            `i \mapsto i^*` of the vertices Dynkin diagram such that,
            for `w_0` the longest element of the Weyl group, and any
            simple root `\alpha_i`, one has `\alpha_{i^*} = -w_0(\alpha_i)`.

            The automorphism is returned as a dictionary.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',5],shape=[1])
                sage: T.opposition_automorphism()
                {1: 5, 2: 4, 3: 3, 4: 2, 5: 1}

                sage: T = CrystalOfTableaux(['D',4],shape=[1])
                sage: T.opposition_automorphism()
                {1: 1, 2: 2, 3: 3, 4: 4}

                sage: T = CrystalOfTableaux(['D',5],shape=[1])
                sage: T.opposition_automorphism()
                {1: 1, 2: 2, 3: 3, 4: 5, 5: 4}

                sage: T = CrystalOfTableaux(['C',4],shape=[1])
                sage: T.opposition_automorphism()
                {1: 1, 2: 2, 3: 3, 4: 4}
            """
            L = self.cartan_type().root_system().root_lattice()
            W = L.weyl_group()
            w0 = W.long_element()
            alpha = L.simple_roots()
            return dict( (i, (w0.action(alpha[i])).leading_support()) for i in self.index_set() )

        def demazure_character(self, w, f = None):
            r"""
            Returns the Demazure character associated to ``w``.

            INPUT:

                - ``w`` -- an element of the ambient weight lattice
                  realization of the crystal, or a reduced word, or an element
                  in the associated Weyl group

            OPTIONAL:

                - ``f`` -- a function from the crystal to a module.

            This is currently only supported for crystals whose underlying
            weight space is the ambient space.

            The Demazure character is obtained by applying the Demazure operator
            `D_w` (see :meth:`sage.categories.crystals.Crystals.ParentMethods.demazure_operator`)
            to the highest weight element of the classical crystal. The simple
            Demazure operators `D_i` (see :meth:`sage.categories.crystals.Crystals.ElementMethods.demazure_operator_simple`)
            do not braid on the level of crystals, but on the level of characters they do.
            That is why it makes sense to input ``w`` either as a weight, a reduced word,
            or as an element of the underlying Weyl group.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape = [2,1])
                sage: e = T.weight_lattice_realization().basis()
                sage: weight = e[0] + 2*e[2]
                sage: weight.reduced_word()
                [2, 1]
                sage: T.demazure_character(weight)
                x1^2*x2 + x1*x2^2 + x1^2*x3 + x1*x2*x3 + x1*x3^2

                sage: T = CrystalOfTableaux(['A',3],shape=[2,1])
                sage: T.demazure_character([1,2,3])
                x1^2*x2 + x1*x2^2 + x1^2*x3 + x1*x2*x3 + x2^2*x3
                sage: W = WeylGroup(['A',3])
                sage: w = W.from_reduced_word([1,2,3])
                sage: T.demazure_character(w)
                x1^2*x2 + x1*x2^2 + x1^2*x3 + x1*x2*x3 + x2^2*x3

                sage: T = CrystalOfTableaux(['B',2], shape = [2])
                sage: e = T.weight_lattice_realization().basis()
                sage: weight = -2*e[1]
                sage: T.demazure_character(weight)
                x1^2 + x1*x2 + x2^2 + x1 + x2 + x1/x2 + 1/x2 + 1/x2^2 + 1

                sage: T = CrystalOfTableaux("B2",shape=[1/2,1/2])
                sage: b2=WeylCharacterRing("B2",base_ring=QQ).ambient()
                sage: T.demazure_character([1,2],f=lambda x:b2(x.weight()))
                b2(-1/2,1/2) + b2(1/2,-1/2) + b2(1/2,1/2)

            REFERENCES::

                .. [D1974] M. Demazure, Desingularisation des varietes de Schubert,
                   Ann. E. N. S., Vol. 6, (1974), p. 163-172

                .. [M2009] Sarah Mason, An Explicit Construction of Type A Demazure Atoms,
                   Journal of Algebraic Combinatorics, Vol. 29, (2009), No. 3, p.295-313
                   (arXiv:0707.4267)

            """
            from sage.misc.misc_c import prod
            from sage.rings.integer_ring import ZZ
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            if hasattr(w, 'reduced_word'):
                word = w.reduced_word()
            else:
                word = w
            n = self.weight_lattice_realization().n
            u = self.algebra(ZZ).sum_of_monomials(self.module_generators)
            u = self.demazure_operator(u, word)
            if f is None:
                x = ['x%s'%i for i in range(1,n+1)]
                P = PolynomialRing(ZZ, x)
                # TODO: use P.linear_combination when PolynomialRing will be a ModulesWithBasis
                return sum((coeff*prod((x[i]**(c.weight()[i]) for i in range(n)), P.one()) for c, coeff in u), P.zero())
            else:
                return sum((coeff*f(c)) for c, coeff in u)

        def character(self, R=None):
            """
            Returns the character of this crystal.

            INPUT:

              - ``R`` -- a :class:`WeylCharacterRing`
                (default: the default :class:`WeylCharacterRing` for this Cartan type)

            Returns the character of ``self`` as an element of ``R``.

            EXAMPLES::

                sage: C = CrystalOfTableaux("A2", shape=[2,1])
                sage: chi = C.character(); chi
                A2(2,1,0)

                sage: T = TensorProductOfCrystals(C,C)
                sage: chiT = T.character(); chiT
                A2(2,2,2) + 2*A2(3,2,1) + A2(3,3,0) + A2(4,1,1) + A2(4,2,0)
                sage: chiT == chi^2
                True

            One may specify an alternate :class:`WeylCharacterRing`::

                sage: R = WeylCharacterRing("A2", style="coroots")
                sage: chiT = T.character(R); chiT
                A2(0,0) + 2*A2(1,1) + A2(0,3) + A2(3,0) + A2(2,2)
                sage: chiT in R
                True

            It should have the same Cartan type and use the same
            realization of the weight lattice as ``self``::

                sage: R = WeylCharacterRing("A3", style="coroots")
                sage: T.character(R)
                Traceback (most recent call last):
                ...
                ValueError: Weyl character ring does not have the right Cartan type

            """
            from sage.combinat.root_system.weyl_characters import WeylCharacterRing
            if R == None:
                R = WeylCharacterRing(self.cartan_type())
            if not R.cartan_type() == self.cartan_type():
                raise ValueError, "Weyl character ring does not have the right Cartan type"
            assert R.basis().keys() == self.weight_lattice_realization()

            return R.sum_of_monomials( x.weight() for x in self.highest_weight_vectors() )

        @cached_method
        def q_dimension(self, q=None, use_product=False):
            r"""
            Return the `q`-dimension of ``self``.

            Let `B(\lambda)` denote a highest weight crystal. Recall that
            the degree of the weight space `\alpha` of `B(\lambda)` (under
            the principal gradation) is equal to
            `\langle \rho^{\vee}, \lambda - \alpha \rangle` where
            `\langle \rho^{\vee}, \alpha_i \rangle = 1` for all `i \in I`
            (in particular, take `\rho^{\vee} = \sum_{i \in I} h_i`).

            The `q`-dimension of a highest weight crystal `B(\Lambda)` is
            defined as

            .. MATH::

                \dim_q B(\lambda) := \sum_{j \geq 0} \dim(B_j) q^j,

            where `B_j` denotes the degree `j` portion of `B(\lambda)`. This
            can be expressed as the product

            .. MATH::

                \dim_q B(\lambda) = \prod_{\alpha \in \Delta_+^{\vee}}
                \left( \frac{1 - q^{\langle \lambda + \rho, \alpha \rangle}}
                {1 - q^{\langle \rho, \alpha \rangle}} \right)^{\mathrm{mult}
                \, \alpha},

            where `\Delta_+^{\vee}` denotes the set of positive coroots.
            Taking the limit as `q \to 1` gives the dimension of `B(\lambda)`.

            INPUT:

            - ``q`` -- the (generic) parameter `q`
            - ``use_product`` -- (default: ``False``) if ``True``, use the
              product formula

            EXAMPLES::

                sage: C = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: qdim = C.q_dimension(); qdim
                q^4 + 2*q^3 + 2*q^2 + 2*q + 1
                sage: qdim(1)
                8
                sage: len(C) == qdim(1)
                True
                sage: C.q_dimension(use_product=True) == qdim
                True

                sage: C = CrystalOfTableaux(['A',2], shape=[5,2])
                sage: C.q_dimension()
                q^10 + 2*q^9 + 4*q^8 + 5*q^7 + 6*q^6 + 6*q^5
                 + 6*q^4 + 5*q^3 + 4*q^2 + 2*q + 1

                sage: C = CrystalOfTableaux(['B',2], shape=[2,1])
                sage: C.q_dimension()
                q^7 + 3*q^6 + 5*q^5 + 7*q^4 + 7*q^3 + 6*q^2 + 4*q + 2

                sage: C = CrystalOfTableaux(['D',4], shape=[2,1])
                sage: C.q_dimension()
                q^16 + 2*q^15 + 4*q^14 + 7*q^13 + 10*q^12 + 13*q^11
                 + 16*q^10 + 18*q^9 + 18*q^8 + 18*q^7 + 16*q^6 + 13*q^5
                 + 10*q^4 + 7*q^3 + 4*q^2 + 2*q + 1

                sage: TP = TensorProductOfCrystals(C, C)
                sage: TP.cardinality()
                25600
                sage: qdim = TP.q_dimension(use_product=True); qdim
                q^32 + 2*q^31 + 8*q^30 + 15*q^29 + 34*q^28 + 63*q^27 + 110*q^26
                 + 175*q^25 + 276*q^24 + 389*q^23 + 550*q^22 + 725*q^21
                 + 930*q^20 + 1131*q^19 + 1362*q^18 + 1548*q^17 + 1736*q^16
                 + 1858*q^15 + 1947*q^14 + 1944*q^13 + 1918*q^12 + 1777*q^11
                 + 1628*q^10 + 1407*q^9 + 1186*q^8 + 928*q^7 + 720*q^6
                 + 498*q^5 + 342*q^4 + 201*q^3 + 117*q^2 + 48*q + 26
                sage: qdim(1)
                25600
                sage: TP.q_dimension() == qdim # long time
                True
            """
            from sage.rings.all import ZZ
            if q is None:
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                q = PolynomialRing(ZZ, 'q').gen(0)

            P = q.parent()
            WLR = self.weight_lattice_realization()
            rho = WLR.rho()

            if use_product:
                # since we are in the classical case, all roots occur with multiplicity 1
                pos_roots = WLR.positive_roots()
                ret = P.zero()
                for v in self.highest_weight_vectors():
                    hw = v.weight()
                    ret += P.prod((1 - q**(rho+hw).scalar(alpha)) / (1 - q**rho.scalar(alpha))
                                  for alpha in pos_roots)
                    # We do a cast since the result would otherwise live in the fraction field
                return P(ret)

            ret = P.zero()
            for v in self.highest_weight_vectors():
                hw = v.weight()
                S = self.subcrystal(generators=[v], direction='lower')
                ret += P.sum(q**rho.scalar(hw - b.weight()) for b in S)
            return ret

        def __iter__(self):
            r"""
            Returns an iterator over the elements of this crystal.

            This iterator uses little memory, storing only one element
            of the crystal at a time. For details on the complexity, see
            :class:`sage.combinat.crystals.crystals.CrystalBacktracker`.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: [x for x in C]
                [1, 2, 3, 4, 5, 6]

            TESTS::

                sage: C = CrystalOfLetters(['D',4])
                sage: D = CrystalOfSpinsPlus(['D',4])
                sage: E = CrystalOfSpinsMinus(['D',4])
                sage: T=TensorProductOfCrystals(D,E,generators=[[D.list()[0],E.list()[0]]])
                sage: U=TensorProductOfCrystals(C,E,generators=[[C(1),E.list()[0]]])
                sage: T.cardinality()
                56

                sage: TestSuite(T).run(verbose = True)
                running ._test_an_element() . . . pass
                running ._test_category() . . . pass
                running ._test_elements() . . .
                  Running the test suite of self.an_element()
                  running ._test_category() . . . pass
                  running ._test_eq() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                  running ._test_stembridge_local_axioms() . . . pass
                  pass
                running ._test_elements_eq_reflexive() . . . pass
                running ._test_elements_eq_symmetric() . . . pass
                running ._test_elements_eq_transitive() . . . pass
                running ._test_elements_neq() . . . pass
                running ._test_enumerated_set_contains() . . . pass
                running ._test_enumerated_set_iter_cardinality() . . . pass
                running ._test_enumerated_set_iter_list() . . . pass
                running ._test_eq() . . . pass
                running ._test_fast_iter() . . . pass
                running ._test_not_implemented_methods() . . . pass
                running ._test_pickling() . . . pass
                running ._test_some_elements() . . . pass
                running ._test_stembridge_local_axioms() . . . pass

                sage: TestSuite(U).run(verbose = True)
                running ._test_an_element() . . . pass
                running ._test_category() . . . pass
                running ._test_elements() . . .
                  Running the test suite of self.an_element()
                  running ._test_category() . . . pass
                  running ._test_eq() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                  running ._test_stembridge_local_axioms() . . . pass
                  pass
                running ._test_elements_eq_reflexive() . . . pass
                running ._test_elements_eq_symmetric() . . . pass
                running ._test_elements_eq_transitive() . . . pass
                running ._test_elements_neq() . . . pass
                running ._test_enumerated_set_contains() . . . pass
                running ._test_enumerated_set_iter_cardinality() . . . pass
                running ._test_enumerated_set_iter_list() . . . pass
                running ._test_eq() . . . pass
                running ._test_fast_iter() . . . pass
                running ._test_not_implemented_methods() . . . pass
                running ._test_pickling() . . . pass
                running ._test_some_elements() . . . pass
                running ._test_stembridge_local_axioms() . . . pass

            Bump's systematic tests::

                sage: fa3 = lambda a,b,c: CrystalOfTableaux(['A',3],shape=[a+b+c,b+c,c])
                sage: fb3 = lambda a,b,c: CrystalOfTableaux(['B',3],shape=[a+b+c,b+c,c])
                sage: fc3 = lambda a,b,c: CrystalOfTableaux(['C',3],shape=[a+b+c,b+c,c])
                sage: fb4 = lambda a,b,c,d: CrystalOfTableaux(['B',4],shape=[a+b+c+d,b+c+d,c+d,d])
                sage: fd4 = lambda a,b,c,d: CrystalOfTableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d])
                sage: fd5 = lambda a,b,c,d,e: CrystalOfTableaux(['D',5],shape=[a+b+c+d+e,b+c+d+e,c+d+e,d+e,e])
                sage: def fd4spinplus(a,b,c,d):
                ...     C = CrystalOfTableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d])
                ...     D = CrystalOfSpinsPlus(['D',4])
                ...     return TensorProductOfCrystals(C,D,generators=[[C[0],D[0]]])
                sage: def fb3spin(a,b,c):
                ...     C = CrystalOfTableaux(['B',3],shape=[a+b+c,b+c,c])
                ...     D = CrystalOfSpins(['B',3])
                ...     return TensorProductOfCrystals(C,D,generators=[[C[0],D[0]]])

            TODO: choose a good panel of values for a,b,c ... both for
            basic systematic tests and for conditionally run,
            computationally involved tests.

            ::

                sage: TestSuite(fb4(1,0,1,0)).run(verbose = True)  # long time (8s on sage.math, 2011)
                running ._test_an_element() . . . pass
                running ._test_category() . . . pass
                running ._test_elements() . . .
                  Running the test suite of self.an_element()
                  running ._test_category() . . . pass
                  running ._test_eq() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                  running ._test_stembridge_local_axioms() . . . pass
                  pass
                running ._test_elements_eq_reflexive() . . . pass
                running ._test_elements_eq_symmetric() . . . pass
                running ._test_elements_eq_transitive() . . . pass
                running ._test_elements_neq() . . . pass
                running ._test_enumerated_set_contains() . . . pass
                running ._test_enumerated_set_iter_cardinality() . . . pass
                running ._test_enumerated_set_iter_list() . . . pass
                running ._test_eq() . . . pass
                running ._test_fast_iter() . . . pass
                running ._test_not_implemented_methods() . . . pass
                running ._test_pickling() . . . pass
                running ._test_some_elements() . . . pass
                running ._test_stembridge_local_axioms() . . . pass

            ::

                #sage: fb4(1,1,1,1).check() # expensive: the crystal is of size 297297
                #True
            """
            from sage.combinat.crystals.crystals import CrystalBacktracker
            return iter(CrystalBacktracker(self))

        def _test_fast_iter(self, **options):
            r"""
            Tests whether the elements returned by :meth:`.__iter__`
            and ``Crystal.list(self)`` are the same (the two
            algorithms are different).

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C._test_fast_iter()
            """
            tester = self._tester(**options)
            S = list(self)
            SS  = list(Crystals().parent_class.__iter__(self))
            tester.assert_( len(S) == len(SS) )
            tester.assert_( len(S) == len(set(S)))
            tester.assert_( set(S) == set(SS) )

        def cardinality(self):
            r"""
            Returns the number of elements of the crystal, using Weyl's
            dimension formula on each connected component.

            EXAMPLES::

                sage: C = ClassicalCrystals().example(5)
                sage: C.cardinality()
                6
            """
            return sum(self.weight_lattice_realization().weyl_dimension(x.weight())
                       for x in self.highest_weight_vectors())

    class ElementMethods:

        def lusztig_involution(self):
            r"""
            Returns the Lusztig involution on the classical highest weight crystal self.

            The Lusztig involution on a finite-dimensional highest weight crystal `B(\lambda)` of highest weight `\lambda`
            maps the highest weight vector to the lowest weight vector and the Kashiwara operator `f_i` to
            `e_{i^*}`, where `i^*` is defined as `\alpha_{i^*} = -w_0(\alpha_i)`. Here `w_0` is the longest element
            of the Weyl group acting on the `i`-th simple root `\alpha_i`.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',3],shape=[2,1])
                sage: b = B(rows=[[1,2],[4]])
                sage: b.lusztig_involution()
                [[1, 4], [3]]
                sage: b.to_tableau().schuetzenberger_involution(n=4)
                [[1, 4], [3]]

                sage: all(b.lusztig_involution().to_tableau() == b.to_tableau().schuetzenberger_involution(n=4) for b in B)
                True

                sage: B = CrystalOfTableaux(['D',4],shape=[1])
                sage: [[b,b.lusztig_involution()] for b in B]
                [[[[1]], [[-1]]], [[[2]], [[-2]]], [[[3]], [[-3]]], [[[4]], [[-4]]], [[[-4]],
                [[4]]], [[[-3]], [[3]]], [[[-2]], [[2]]], [[[-1]], [[1]]]]

                sage: B = CrystalOfTableaux(['D',3],shape=[1])
                sage: [[b,b.lusztig_involution()] for b in B]
                [[[[1]], [[-1]]], [[[2]], [[-2]]], [[[3]], [[3]]], [[[-3]], [[-3]]],
                [[[-2]], [[2]]], [[[-1]], [[1]]]]

                sage: C=CartanType(['E',6])
                sage: La=C.root_system().weight_lattice().fundamental_weights()
                sage: T = HighestWeightCrystal(La[1])
                sage: t = T[3]; t
                [(-4, 2, 5)]
                sage: t.lusztig_involution()
                [(-2, -3, 4)]
            """
            hw = self.to_highest_weight()[1]
            hw.reverse()
            hw = [self.parent().opposition_automorphism()[i] for i in hw]
            return self.to_lowest_weight()[0].e_string(hw)


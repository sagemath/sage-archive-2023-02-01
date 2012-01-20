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
from sage.categories.category import Category
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals

class ClassicalCrystals(Category):
    """
    The category of classical crystals, that is crystals of finite Cartan type.

    EXAMPLES::

        sage: C = ClassicalCrystals()
        sage: C
        Category of classical crystals
        sage: C.super_categories()
        [Category of finite crystals, Category of highest weight crystals]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = FiniteCrystals().example()
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
        running ._test_elements_eq() . . . pass
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

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: ClassicalCrystals().super_categories()
            [Category of finite crystals, Category of highest weight crystals]
        """
        return [FiniteCrystals(), HighestWeightCrystals()]

    def example(self, n = 3):
        """
        Returns an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = ClassicalCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.categories.crystals import Crystals
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

        def demazure_character(self, weight, reduced_word = False):
            r"""
            Returns the Demazure character associated to the specified
            weight in the ambient weight lattice.

            INPUT:

                - ``weight`` -- an element of the weight lattice
                  realization of the crystal, or a reduced word
                - ``reduced_word`` -- a boolean (default: ``False``)
                  whether ``weight`` is given as a reduced word

            This is currently only supported for crystals whose
            underlying weight space is the ambient space.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',2], shape = [2,1])
                sage: e = T.weight_lattice_realization().basis()
                sage: weight = e[0] + 2*e[2]
                sage: weight.reduced_word()
                [2, 1]
                sage: T.demazure_character(weight)
                x1^2*x2 + x1^2*x3 + x1*x2^2 + x1*x2*x3 + x1*x3^2

                sage: T = CrystalOfTableaux(['A',3],shape=[2,1])
                sage: T.demazure_character([1,2,3], reduced_word = True)
                x1^2*x2 + x1^2*x3 + x1*x2^2 + x1*x2*x3 + x2^2*x3

                sage: T = CrystalOfTableaux(['B',2], shape = [2])
                sage: e = T.weight_lattice_realization().basis()
                sage: weight = -2*e[1]
                sage: T.demazure_character(weight)
                x1^2 + x1*x2 + x2^2 + x1 + x2 + x1/x2 + 1/x2 + 1/x2^2 + 1

            TODO: detect automatically if weight is a reduced word,
            and remove the (untested!) ``reduced_word`` option.

            REFERENCES::

                .. [D1974] M. Demazure, Desingularisation des varietes de Schubert,
                   Ann. E. N. S., Vol. 6, (1974), p. 163-172

                .. [M2009] Sarah Mason, An Explicit Construction of Type A Demazure Atoms,
                   Journal of Algebraic Combinatorics, Vol. 29, (2009), No. 3, p.295-313
                   (arXiv:0707.4267)

            """
            from sage.misc.misc_c import prod
            from sage.rings.rational_field import QQ
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            if reduced_word:
                word = weight
            else:
                word = weight.reduced_word()
            n = self.weight_lattice_realization().n
            u = list( self.module_generators )
            for i in reversed(word):
                u = u + sum((x.demazure_operator(i, truncated = True) for x in u), [])
            x = ['x%s'%i for i in range(1,n+1)]
            P = PolynomialRing(QQ, x)
            u = [b.weight() for b in u]
            return sum((prod((x[i]**(la[i]) for i in range(n)), P.one()) for la in u), P.zero())

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

            It should have the same cartan type and use the same
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

        def list(self):
            r"""
            Returns the list of the elements of ``self``, as per
            :meth:`FiniteEnumeratedSets.ParentMethods.list`

            EXAMPLES::

                sage: C = CrystalOfLetters(['D',4])
                sage: C.list()
                [1, 2, 3, 4, -4, -3, -2, -1]

            FIXME: this is just there to reinstate the default
            implementation of :meth:`.list` from :meth:`.__iter__`
            which is overriden in :class:`Crystals`.
            """
            return self._list_from_iterator()

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
                running ._test_elements_eq() . . . pass
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
                running ._test_elements_eq() . . . pass
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
                sage: def fd4spinplus(a,b,c,d):\
                     C = CrystalOfTableaux(['D',4],shape=[a+b+c+d,b+c+d,c+d,d]);\
                     D = CrystalOfSpinsPlus(['D',4]);\
                     return TensorProductOfCrystals(C,D,generators=[[C[0],D[0]]])
                sage: def fb3spin(a,b,c):\
                     C = CrystalOfTableaux(['B',3],shape=[a+b+c,b+c,c]);\
                     D = CrystalOfSpins(['B',3]);\
                     return TensorProductOfCrystals(C,D,generators=[[C[0],D[0]]])

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
                running ._test_elements_eq() . . . pass
                running ._test_enumerated_set_contains() . . . pass
                running ._test_enumerated_set_iter_cardinality() . . . pass
                running ._test_enumerated_set_iter_list() . . .Enumerated set too big; skipping test; see ``self.max_test_enumerated_set_loop``
                pass
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
            S = self._list_brute_force()
            SS = set(S)
            tester.assert_( len(S) == len(SS) )
            tester.assert_( set(self) == set(SS) )

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
                [[-4, 2, 5]]
                sage: t.lusztig_involution()
                [[-2, -3, 4]]
            """
            hw = self.to_highest_weight()[1]
            hw.reverse()
            hw = [self.parent().opposition_automorphism()[i] for i in hw]
            return self.to_lowest_weight()[0].e_string(hw)

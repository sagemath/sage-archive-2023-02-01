"""
Weight lattice realizations
"""
#*****************************************************************************
#       Copyright (C) 2007-2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#       (with contributions of many others)
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

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.misc import prod
from sage.categories.category_types import Category_over_base_ring
from sage.combinat.family import Family
from root_lattice_realizations import RootLatticeRealizations

class WeightLatticeRealizations(Category_over_base_ring):
    r"""
    The category of weight lattice realizations over a given base ring

    A *weight lattice realization* `L` over a base ring `R` is a free
    module (or vector space if `R` is a field) endowed with an embedding
    of the root lattice of some root system. By restriction, this
    embedding defines an embedding of the root lattice of this root
    system, which makes `L` a root lattice realization.

    Typical weight lattice realizations over `\ZZ` include the weight
    lattice, and ambient lattice. Typical weight lattice realizations
    over `\QQ` include the weight space, and ambient space.

    To describe the embedding, a weight lattice realization must
    implement a method
    :meth:`~RootLatticeRealizations.ParentMethods.fundamental_weight`(i)
    returning for each `i` in the index set the image of the fundamental
    weight `\Lambda_i` under the embedding.

    In order to be a proper root lattice realization, a weight lattice
    realization should also implement the scalar product with the coroot
    lattice; on the other hand, the embedding of the simple roots is
    given for free.

    .. seealso::

        - :class:`~sage.combinat.root_system.root_system.RootSystem`
        - :class:`~sage.combinat.root_system.root_lattice_realizations.RootLatticeRealizations`
        - :class:`~sage.combinat.root_system.weight_space.WeightSpace`
        - :class:`~sage.combinat.root_system.ambient_space.AmbientSpace`

    EXAMPLES:

    Here, we consider the root system of type `A_7`, and embed the weight
    lattice element `x = \Lambda_1 + 2 \Lambda_3` in several root lattice
    realizations::

        sage: R = RootSystem(["A",7])
        sage: Lambda = R.weight_lattice().fundamental_weights()
        sage: x = Lambda[2] + 2 * Lambda[5]

        sage: L = R.weight_space()
        sage: L(x)
        Lambda[2] + 2*Lambda[5]

        sage: L = R.ambient_lattice()
        sage: L(x)
        (3, 3, 2, 2, 2, 0, 0, 0)

    We embed the weight space element `x = \Lambda_1 + 1/2 \Lambda_3` in
    the ambient space::

        sage: Lambda = R.weight_space().fundamental_weights()
        sage: x = Lambda[2] + 1/2 * Lambda[5]

        sage: L = R.ambient_space()
        sage: L(x)
        (3/2, 3/2, 1/2, 1/2, 1/2, 0, 0, 0)

    Of course, one can't embed the weight space in the ambient lattice::

        sage: L = R.ambient_lattice()
        sage: L(x)
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= Lambda[2] + 1/2*Lambda[5]) an element of self (=Ambient lattice of the Root system of type ['A', 7])

    If `K_1` is a subring of `K_2`, then one could in theory have an
    embedding from the weight space over `K_1` to any weight lattice
    realization over `K_2`; this is not implemented::

        sage: K1 = QQ
        sage: K2 = QQ['q']
        sage: L = R.ambient_space(K2)

        sage: Lambda = R.weight_space(K2).fundamental_weights()
        sage: L(Lambda[1])
        (1, 0, 0, 0, 0, 0, 0, 0)

        sage: Lambda = R.weight_space(K1).fundamental_weights()
        sage: L(Lambda[1])
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= Lambda[1]) an element of self (=Ambient space of the Root system of type ['A', 7])
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.weight_lattice_realizations import WeightLatticeRealizations
            sage: WeightLatticeRealizations(QQ).super_categories()
            [Category of root lattice realizations over Rational Field]
        """
        return [RootLatticeRealizations(self.base_ring())]

    class ParentMethods:

        @abstract_method
        def fundamental_weight(self, i):
            """
            Returns the `i^{th}` fundamental weight

            INPUT:

            - ``i`` -- an element of the index set

            By a slight notational abuse, for an affine type this method
            should also accept ``"delta"`` as input, and return the image
            of `\delta` of the extended weight lattice in this
            realization.

            This should be overridden by any subclass, and typically
            be implemented as a cached method for efficiency.

            EXAMPLES::

                sage: L = RootSystem(["A",3]).ambient_lattice()
                sage: L.fundamental_weight(1)
                (1, 0, 0, 0)

                sage: L = RootSystem(["A",3,1]).weight_lattice(extended=True)
                sage: L.fundamental_weight(1)
                Lambda[1]
                sage: L.fundamental_weight("delta")
                delta

            TESTS::

                sage: super(sage.combinat.root_system.weight_space.WeightSpace, L).fundamental_weight(1)
                Traceback (most recent call last):
                ...
                NotImplementedError: <abstract method fundamental_weight at ...>
            """

        def is_extended(self):
          """
          Returns whether this is a realization of the extended weight lattice

          .. seealso:: :class:`sage.combinat.root_system.weight_space.WeightSpace`

          EXAMPLES::

              sage: RootSystem(["A",3,1]).weight_lattice().is_extended()
              False
              sage: RootSystem(["A",3,1]).weight_lattice(extended=True).is_extended()
              True

          This method is irrelevant for finite root systems, since the
          weight lattice need not be extended to ensure that the root
          lattice embeds faithfully::

              sage: RootSystem(["A",3]).weight_lattice().is_extended()
              False

          """
          return False

        def __init_extra__(self):
            """
            Registers the embedding of the weight lattice into ``self``

            Also registers the embedding of the weight space over the same
            base field `K` into ``self`` if `K` is not `\ZZ`.

            If ``self`` is a realization of the extended weight lattice,
            then the embeddings from the extended weight space/lattices
            are registered instead.

            EXAMPLES:

            We embed the fundamental weight `\Lambda_1` of the weight
            lattice in the ambient lattice::

                sage: R = RootSystem(["A",3])
                sage: Lambda = R.root_lattice().simple_roots()
                sage: L = R.ambient_space()
                sage: L(Lambda[2])
                (0, 1, -1, 0)

            .. note::

                More examples are given in :class:`WeightLatticeRealizations`;
                The embeddings are systematically tested in
                :meth:`_test_weight_lattice_realization`.
            """
            from sage.rings.all import ZZ
            from weight_space import WeightSpace
            K = self.base_ring()
            # If self is the root lattice or the root space, we don't want
            # to register its trivial embedding into itself. This builds
            # the domains from which we want to register an embedding.
            domains = []
            if not isinstance(self, WeightSpace) or K is not ZZ:
                domains.append(self.root_system.weight_lattice(extended=self.is_extended()))
            if not isinstance(self, WeightSpace):
                domains.append(self.root_system.weight_space(K,extended=self.is_extended()))
            # Build and register the embeddings
            for domain in domains:
                domain.module_morphism(self.fundamental_weight,
                                       codomain = self
                                       ).register_as_coercion()


        def _test_weight_lattice_realization(self, **options):
            """
            Runs sanity checks on this weight lattice realization

            - scalar products between the fundamental weights and simple coroots
            - embeddings from the weight lattice and weight space
            - rho, highest_root, ...

            .. seealso:: :class:`TestSuite`

            EXAMPLES::

                sage: RootSystem(['A',3]).weight_lattice()._test_weight_lattice_realization()
            """
            from sage.rings.all import ZZ
            tester     = self._tester(**options)
            Lambda     = self.fundamental_weights()
            alphacheck = self.simple_coroots()
            tester.assertEqual(tuple(Lambda.keys()), self.index_set())

            # Check the consistency between simple_root and simple_roots
            for i in self.index_set():
                tester.assertEqual(self.fundamental_weight(i), Lambda[i])

            # Check the embeddings from:
            # - the weight lattice
            # - the weight space over the same base ring
            #
            # For an affine root system, this will check the embedding of
            # the extended ones, and also of the non extended ones if this
            # realization is not extended
            domains = [self.root_system.weight_space(base_ring, extended = extended)
                       for base_ring in set([ZZ, self.base_ring()])
                       for extended  in set([self.cartan_type().is_affine(), self.is_extended()])]
            for domain in domains:
                tester.assert_(self._internal_coerce_map_from(domain) is not None)
                for i in self.index_set():
                    # This embedding maps fundamental weights to fundamental weights
                    tester.assertEqual(self(domain.fundamental_weight(i)), Lambda[i])
                if self.cartan_type().is_affine():
                    tester.assertEqual(self(domain.null_root()), self.null_root())
                    if self.is_extended():
                        a = self.cartan_type().col_annihilator()
                        # This could be an over specification; we
                        # could imagine realizations of the extended
                        # weight lattice where the null root would not
                        # be a (multiple of) basis element.
                        tester.assertEqual(self.null_root(), self.term("delta", a[0]))
                    for i in self.index_set():
                        # The level of the fundamental weights is consistent
                        tester.assertEqual(domain.fundamental_weight(i).level(), Lambda[i].level())

            # Check that the fundamental weights form the dual basis of the simple coroots
            for i in self.index_set():
                assert(Lambda[i].is_dominant())
                for j in self.index_set():
                    tester.assertEqual(Lambda[j].scalar(alphacheck[i]), (1 if i==j else 0))

            tester.assert_(self.rho().is_dominant())
            if self.root_system.is_finite() and self.root_system.is_irreducible():
                tester.assert_(self.highest_root().is_dominant())

        @cached_method
        def fundamental_weights(self):
            r"""
            Returns the family `(\Lambda_i)_{i\in I}` of the fundamental weights.

            EXAMPLES::

                sage: e = RootSystem(['A',3]).ambient_lattice()
                sage: f = e.fundamental_weights()
                sage: [f[i] for i in [1,2,3]]
                [(1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0)]
            """
            return Family(self.index_set(), self.fundamental_weight)
            # It would be nice to give this family a nice name with
            # ``rename``, but this currently break some doctests.

        @cached_method
        def simple_root(self, i):
            r"""
            Returns the `i`-th simple root

            This default implementation takes the `i`-th simple root in
            the weight lattice and embeds it in ``self``.

            EXAMPLES:

            Since all the weight lattice realizations in Sage currently
            implement a simple_root method, we have to call this one by
            hand::

                sage: from sage.combinat.root_system.weight_lattice_realizations import WeightLatticeRealizations
                sage: simple_root = WeightLatticeRealizations(QQ).parent_class.simple_root.f
                sage: L = RootSystem("A3").ambient_space()
                sage: simple_root(L, 1)
                (1, -1, 0, 0)
                sage: simple_root(L, 2)
                (0, 1, -1, 0)
                sage: simple_root(L, 3)
                (1, 1, 2, 0)

            Note that this last root differs from the one implemented in
            ``L`` by a multiple of the vector ``(1,1,1,1)``::

                sage: L.simple_roots()
                Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1)}

            This is a harmless artefact of the `SL` versus `GL`
            interpretation of type `A`; see the thematic tutorial on Lie
            Methods and Related Combinatorics in Sage for details.
            """
            assert i in self.index_set()
            alphai = self.root_system.weight_lattice().simple_root(i)
            # Note: it would be nicer to just return ``self(alpha[i])``,
            # However the embedding from the weight lattice is defined
            # after the embedding from the root lattice, and the later
            # uses the simple roots. So we compute that embedding by hand.
            Lambda = self.fundamental_weights()
            return self.linear_combination( (Lambda[j], c) for j,c in alphai )

        @cached_method
        def rho(self):
            """
            EXAMPLES::

                sage: RootSystem(['A',3]).ambient_lattice().rho()
                (3, 2, 1, 0)
            """
            return sum(self.fundamental_weights())

        def reduced_word_of_alcove_morphism(self, f):
            r"""
            INPUT:

            - `f` -- a linear map from ``self`` to ``self`` which
              preserves alcoves.

            Let `A` be the fundamental alcove. This returns a reduced word
            `i_1,...,i_k` such that the affine Weyl group element `w =
            s_{i_1} \circ \dots \circ s_{i_k}` maps the alcove `f(A)` back
            to `A`. In other words, the alcove walk `i_1,...,i_k` brings
            the fundamental alcove to the corresponding translated alcove.

            Let us throw in a bit of context to explain the main use
            case.  It is customary to realize the alcove picture in
            the coroot or coweight lattice `R^\vee`. The extended
            affine Weyl group is then the group of linear maps on
            `R^\vee` which preserve the alcoves. By
            [Kac "Infinite-dimensional Lie algebra", Proposition 6.5]
            the affine Weyl group is the semidirect product of the
            associated finite Weyl group and the group of translations
            in the coroot lattice (the extended affine Weyl group uses
            the coweight lattice instead). In other words, an element
            of the extended affine Weyl group admits a unique
            decomposition of the form:

            .. math:: f = d w ,

            where `w` is in the Weyl group, and `d` is a function which
            maps the fundamental alcove to itself. As `d` permutes the
            walls of the fundamental alcove, it permutes accordingly the
            corresponding simple roots, which induces an automorphism of
            the Dynkin diagram.

            This method returns a reduced word for `w`, whereas the method
            :meth:`dynkin_diagram_automorphism_of_alcove_morphism` returns
            `d` as a permutation of the nodes of the Dynkin diagram.

            Nota bene: recall that the coroot (resp. coweight) lattice is
            implemented as the root (resp weight) lattice of the dual root
            system. Hence, this method is implemented for weight lattice
            realizations, but in practice is most of the time used on the
            dual side.

            EXAMPLES:

            We start with type `A` which is simply laced; hence we do not
            have to worry about the distinction between the weight and
            coweight lattice::

                sage: R = RootSystem(["A",2,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: Lambda = R.fundamental_weights()

            We consider first translations by elements of the root lattice::

                sage: R.reduced_word_of_alcove_morphism(alpha[0].translation)
                [1, 2, 1, 0]
                sage: R.reduced_word_of_alcove_morphism(alpha[1].translation)
                [0, 2, 0, 1]
                sage: R.reduced_word_of_alcove_morphism(alpha[2].translation)
                [0, 1, 0, 2]

            We continue with translations by elements of the classical
            weight lattice, embedded at level `0`:

                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - Lambda[0]

                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 2]
                sage: R.reduced_word_of_alcove_morphism(omega2.translation)
                [0, 1]

            The following tests ensure that the code agrees with the tables
            in Kashiwara's private notes on affine quantum algebras (2008).

            TESTS::

                sage: R = RootSystem(['A',5,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 5, 4, 3, 2]
                sage: R.reduced_word_of_alcove_morphism(alpha[0].translation)
                [1, 2, 3, 4, 5, 4, 3, 2, 1, 0]

                sage: R = RootSystem(['C',3,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = 2*(Lambda[1] - Lambda[0])
                sage: omega2 = 2*(Lambda[2] - Lambda[0])
                sage: omega3 = Lambda[3] - Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 1, 2, 3, 2, 1]
                sage: R.reduced_word_of_alcove_morphism(omega2.translation)
                [0, 1, 0, 2, 1, 3, 2, 1, 3, 2]
                sage: R.reduced_word_of_alcove_morphism(omega3.translation)
                [0, 1, 0, 2, 1, 0]
                sage: W = WeylGroup(['C',3,1])
                sage: s = W.simple_reflections()
                sage: w = s[0]*s[1]*s[2]*s[3]*s[2]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega2.translation)) == w*w
                True
                sage: w = s[0]*s[1]*s[2]*s[0]*s[1]*s[0]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega3.translation)) == w
                True

                sage: R = RootSystem(['D',4,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: omega3 = Lambda[3] - Lambda[0]
                sage: omega4 = Lambda[4] - Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 2, 3, 4, 2, 0]
                sage: R.reduced_word_of_alcove_morphism(omega2.translation)
                [0, 2, 1, 3, 2, 4, 2, 1, 3, 2]
                sage: R.reduced_word_of_alcove_morphism(omega3.translation)
                [0, 2, 1, 4, 2, 0]
                sage: R.reduced_word_of_alcove_morphism(omega4.translation)
                [0, 2, 1, 3, 2, 0]
                sage: W = WeylGroup(['D',4,1])
                sage: s = W.simple_reflections()
                sage: w = s[0]*s[2]*s[3]*s[4]*s[2]
                sage: w1= s[1]*s[2]*s[3]*s[4]*s[2]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega2.translation)) == w*w1
                True

                sage: R = RootSystem(['D',5,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: R.reduced_word_of_alcove_morphism(omega1.translation)
                [0, 2, 3, 4, 5, 3, 2, 0]
                sage: W = WeylGroup(['D',5,1])
                sage: s = W.simple_reflections()
                sage: w = s[0]*s[2]*s[3]*s[4]*s[5]*s[3]*s[2]
                sage: w1= s[1]*s[2]*s[3]*s[4]*s[5]*s[3]*s[2]
                sage: W.from_reduced_word(R.reduced_word_of_alcove_morphism(omega2.translation)) == w*w1
                True
            """
            return f(self.rho()).reduced_word()

        def dynkin_diagram_automorphism_of_alcove_morphism(self, f):
            """
            Returns the Dynkin diagram automorphism induced by an alcove morphism

            INPUT:

            - ``f`` - a linear map from ``self`` to ``self`` which preserves alcoves

            This method returns the Dynkin diagram automorphism for
            the decomposition `f = d w` (see
            :meth:`reduced_word_of_alcove_morphism`), as a dictionnary
            mapping elements of the index set to itself.

            EXAMPLES::

                sage: R = RootSystem(["A",2,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: Lambda = R.fundamental_weights()

            Translations by elements of the root lattice induce a
            trivial Dynkin diagram automorphism::

                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[0].translation)
                {0: 0, 1: 1, 2: 2}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[1].translation)
                {0: 0, 1: 1, 2: 2}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[2].translation)
                {0: 0, 1: 1, 2: 2}

            This is no more the case for translations by general
            elements of the (classical) weight lattice at level 0::

                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - Lambda[0]

                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega1.translation)
                {0: 1, 1: 2, 2: 0}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega2.translation)
                {0: 2, 1: 0, 2: 1}

                sage: R = RootSystem(['C',2,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(alpha[1].translation)
                {0: 2, 1: 1, 2: 0}

                sage: R = RootSystem(['D',5,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega1.translation)
                {0: 1, 1: 0, 2: 2, 3: 3, 4: 5, 5: 4}
                sage: R.dynkin_diagram_automorphism_of_alcove_morphism(omega2.translation)
                {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}

            Algorithm: computes `w` of the decomposition, and see how
            `f\circ w^{-1}` permutes the simple roots.
            """
            alpha = self.simple_roots()
            w = self.weyl_group().from_reduced_word(self.reduced_word_of_alcove_morphism(f))
            # Now, we have d = f w^-1
            winv = ~w
            assert all(alpha[i].level().is_zero() for i in self.index_set())
            rank_simple_roots = dict( (alpha[i],i) for i in self.index_set())
            permutation = dict()
            for i in self.index_set():
                root = f(winv.action(alpha[i])) # This is d(alpha_i)
                assert root in rank_simple_roots
                permutation[i] = rank_simple_roots[root]
                assert set(permutation.values()), set(self.index_set())
            return permutation

        def reduced_word_of_translation(self, t):
            """
            Given an element of the root lattice, this returns a reduced
            word `i_1,...,i_k` such that the Weyl group element `s_{i_1}
            \circ \dots \circ s_{i_k}` implements the "translation"
            where `x` maps to `x + level(x)*t`. In other words, the alcove walk
            `i_1,...,i_k` brings the fundamental alcove to the
            corresponding translated alcove.

            Note: there are some technical conditions for `t` to actually
            be a translation; those are not tested (TODO: detail).

            EXAMPLES::

                sage: R = RootSystem(["A",2,1]).weight_lattice()
                sage: alpha = R.simple_roots()
                sage: R.reduced_word_of_translation(alpha[1])
                [0, 2, 0, 1]
                sage: R.reduced_word_of_translation(alpha[2])
                [0, 1, 0, 2]
                sage: R.reduced_word_of_translation(alpha[0])
                [1, 2, 1, 0]

                sage: R = RootSystem(['D',5,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: omega1 = Lambda[1] - Lambda[0]
                sage: omega2 = Lambda[2] - 2*Lambda[0]
                sage: R.reduced_word_of_translation(omega1)
                [0, 2, 3, 4, 5, 3, 2, 0]
                sage: R.reduced_word_of_translation(omega2)
                [0, 2, 1, 3, 2, 4, 3, 5, 3, 2, 1, 4, 3, 2]

            A non simply laced case::

                sage: R = RootSystem(["C",2,1]).weight_lattice()
                sage: Lambda = R.fundamental_weights()
                sage: c = R.cartan_type().translation_factors()
                sage: c
                Finite family {0: 1, 1: 2, 2: 1}
                sage: R.reduced_word_of_translation((Lambda[1]-Lambda[0]) * c[1])
                [0, 1, 2, 1]
                sage: R.reduced_word_of_translation((Lambda[2]-Lambda[0]) * c[2])
                [0, 1, 0]

            See also :meth:`_test_reduced_word_of_translation`.

            TODO:

             - Add a picture in the doc
             - Add a method which, given an element of the classical
               weight lattice, constructs the appropriate value for t
            """
            return self.reduced_word_of_alcove_morphism(t.translation)

        #    # This should be in a method to_weight_lattice()
        #    alphac = self.simple_coroots()
        #    Lambda = self.fundamental_weights()
        #    assert( t == self.plus(t.scalar(alphac[i]) * Lambda[i] for i in self.index_set() ) )
        #    t = self.plus( t.scalar(alphac[i]) * c[i] * Lambda[i] for i in self.index_set() )

        def _test_reduced_word_of_translation(self, elements=None, **options):
            r"""
            Tests the method :meth:`reduced_word_of_translation`.

            INPUT::

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES::

                sage: R = RootSystem(['D',4,1]).weight_lattice()
                sage: R._test_reduced_word_of_translation()

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            if not self.cartan_type().is_affine(): # won't be necessary anymore once root systems are categorified
                return
            alpha = self.simple_roots()
            Lambda = self.fundamental_weights()
            rho = self.rho()
            G = self.dynkin_diagram()
            permutations = []

            # Note: this uses a special set of default elements instead of
            # the usual tester.some_elements(), namely the smallest
            # elements in the weight lattice giving rise to translations
            # preserving the alcoves.
            if elements is None:
                c = self.cartan_type().c()
                elements = [ c[i] * Lambda[i] for i in self.cartan_type().classical().index_set() ]

            # When the null root is zero in this root lattice realization,
            # the roots correspond to the classical roots. We use that to
            # check that w permute the simple roots according to a Dynkin
            # diagram automorphism. This test currently requires the index
            # set to be of the form 0..n
            test_automorphism = self.null_root().is_zero() and set(self.index_set()) == set(i for i in range(len(self.index_set())))
            # dictionary assigning a simple root to its index
            rank_simple_roots = dict( (alpha[i],i) for i in self.index_set() )

            for t in elements:
                t = t - self.base_ring()(t.level()/Lambda[0].level()) * Lambda[0]
                w = self.weyl_group().from_reduced_word(self.reduced_word_of_translation(t))
                if self.null_root().is_zero():
                    # The following formula is only valid when the null root is zero
                    tester.assertEquals(w.action(rho), rho + rho.level()*t)
                    # TODO: fix this formula to take delta into account,
                    # and remove the above condition
                if test_automorphism:
                    permutation = [None for i in self.index_set()]
                    for i in self.index_set():
                        root = w.action(alpha[i])
                        tester.assert_(root in rank_simple_roots)
                        permutation[i] = rank_simple_roots[root]
                    tester.assertEquals(set(permutation), set(self.index_set()))
                    #print permutation
                    # It could be nicer to test equality of G and its relabelling
                    for i in self.index_set():
                        for j in self.index_set():
                            tester.assertEquals(G[permutation[i],permutation[j]], G[i,j])
                    permutations.append(permutation)

            if test_automorphism and elements is None: # note: the test on elements is broken
                # Check that, if we start from all fundamental weights, we
                # get the full automorphism group
                # Disabled: this should actually check that one gets all special
                # automorphisms, which are in bijection with the special nodes
                #from sage.groups.perm_gps.permgroup import PermutationGroup
                #P = PermutationGroup([[i+1 for i in permutation] for permutation in permutations])
                #print P, len(P)
                #tester.assertEquals(P, G.automorphism_group())
                pass

        def signs_of_alcovewalk(self, walk):
            r"""
            Let walk = `[i_1,\dots,i_n]` denote an alcove walk starting
            from the fundamental alcove `y_0`, crossing at step 1 the
            wall `i_1`, and so on.

            For each `k`, set `w_k = s_{i_1} \circ s_{i_k}`, and denote
            by `y_k = w_k(y_0)` the alcove reached after `k` steps. Then,
            `y_k` is obtained recursively from `y_{k-1}` by applying the
            following reflection:

            .. math::

                  y_k = s_{w_{k-1} \alpha_{i_k}} y_{k-1}

            The step is said positive if `w_{k-1} \alpha_{i_k}` is a
            negative root (considering `w_{k-1}` as element of the classical
            Weyl group and `\alpha_{i_k}` as a classical root) and
            negative otherwise.

            This function returns a list of the form `[+1,+1,-1,...]`,
            where the `k^{th}` entry denotes whether the `k^{th}` step was
            positive or negative.

            See equation 3.4, of Ram: Alcove walks ..., arxiv:math/0601343v1 [math.RT]

            EXAMPLES::

                sage: L = RootSystem(['C',2,1]).weight_lattice()
                sage: L.signs_of_alcovewalk([1,2,0,1,2,1,2,0,1,2])
                [-1, -1, 1, -1, 1, 1, 1, 1, 1, 1]
                sage: L = RootSystem(['A',2,1]).weight_lattice()
                sage: L.signs_of_alcovewalk([0,1,2,1,2,0,1,2,0,1,2,0])
                [1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1]
            """
            lattice_classical = self.root_system.cartan_type().classical().root_system().ambient_space()
            W = lattice_classical.weyl_group()
            simple_reflections = W.simple_reflections()
            alphacheck = lattice_classical.alphacheck()
            rho = lattice_classical.rho()
            word = W.one()
            signs = []
            for s in walk:
                if ((alphacheck[s]).scalar((word).action(rho)) > 0):
                    signs.append(-1)
                else:
                    signs.append(1)
                word = simple_reflections[s]*word
            return signs

        def rho_classical(self):
            """
            For an affine type in a weight space, rho_classical is the analog of
            rho in the classical parabolic subgroup. it lives in the level 0.

            EXAMPLES::

                sage: RootSystem(['C',4,1]).weight_space().rho_classical()
                -4*Lambda[0] + Lambda[1] + Lambda[2] + Lambda[3] + Lambda[4]
                sage: WS = RootSystem(['D',4,1]).weight_space()
                sage: WS.rho_classical().scalar(WS.null_coroot())
                0
            """
            rho = self.rho()
            Lambda = self.fundamental_weights()
            return rho - (rho.level()/Lambda[0].level()) * Lambda[0]



        # Should it be a method of highest_weight?
        def weyl_dimension(self, highest_weight):
            """
            EXAMPLES::

                sage: RootSystem(['A',3]).ambient_lattice().weyl_dimension([2,1,0,0])
                20

                sage: type(RootSystem(['A',3]).ambient_lattice().weyl_dimension([2,1,0,0]))
                <type 'sage.rings.integer.Integer'>
            """
            highest_weight = self(highest_weight)
            assert(highest_weight.is_dominant())
            rho = self.rho()
            n = prod([(rho+highest_weight).dot_product(x) for x in self.positive_roots()])
            d = prod([ rho.dot_product(x) for x in self.positive_roots()])
            from sage.rings.integer import Integer
            return Integer(n/d)

# -*- coding: utf-8 -*-
"""
Root lattice realizations
"""
#*****************************************************************************
#       Copyright (C) 2007-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#                          2012 Nicolas Borie  <nicolas.borie at univ-mlv.fr>
#
#       (with contributions of many others)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method, AbstractMethod
from sage.misc.misc import attrcall
from sage.misc.cachefunc import cached_method, cached_in_parent_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.superseded import deprecated_function_alias
from sage.categories.coxeter_groups import CoxeterGroups
from sage.categories.category_types import Category_over_base_ring
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.structure.element import Element
from sage.sets.family import Family
from sage.rings.all import ZZ, QQ
from sage.modules.free_module_element import vector
from sage.combinat.backtrack import TransitiveIdeal, TransitiveIdealGraded
from sage.combinat.root_system.plot import PlotOptions, barycentric_projection_matrix

class RootLatticeRealizations(Category_over_base_ring):
    r"""
    The category of root lattice realizations over a given base ring

    A *root lattice realization* `L` over a base ring `R` is a free
    module (or vector space if `R` is a field) endowed with an embedding
    of the root lattice of some root system.

    Typical root lattice realizations over `\ZZ` include the root
    lattice, weight lattice, and ambient lattice. Typical root lattice
    realizations over `\QQ` include the root space, weight space, and
    ambient space.

    To describe the embedding, a root lattice realization must
    implement a method
    :meth:`~RootLatticeRealizations.ParentMethods.simple_root`
    returning for each `i` in the index set the image of the simple root
    `\alpha_i` under the embedding.

    A root lattice realization must further implement a method on elements
    :meth:`~RootLatticeRealizations.ElementMethods.scalar`, computing
    the scalar product with elements of the coroot lattice or coroot space.

    Using those, this category provides tools for reflections, roots,
    the Weyl group and its action, ...

    .. SEEALSO::

        - :class:`~sage.combinat.root_system.root_system.RootSystem`
        - :class:`~sage.combinat.root_system.weight_lattice_realizations.WeightLatticeRealizations`
        - :class:`~sage.combinat.root_system.root_space.RootSpace`
        - :class:`~sage.combinat.root_system.weight_space.WeightSpace`
        - :class:`~sage.combinat.root_system.ambient_space.AmbientSpace`

    EXAMPLES:

    Here, we consider the root system of type `A_7`, and embed the root
    lattice element `x = \alpha_2 + 2 \alpha_6` in several root lattice
    realizations::

        sage: R = RootSystem(["A",7])
        sage: alpha = R.root_lattice().simple_roots()
        sage: x = alpha[2] + 2 * alpha[5]

        sage: L = R.root_space()
        sage: L(x)
        alpha[2] + 2*alpha[5]

        sage: L = R.weight_lattice()
        sage: L(x)
        -Lambda[1] + 2*Lambda[2] - Lambda[3] - 2*Lambda[4] + 4*Lambda[5] - 2*Lambda[6]

        sage: L = R.ambient_space()
        sage: L(x)
        (0, 1, -1, 0, 2, -2, 0, 0)

    We embed the root space element `x = \alpha_2 + 1/2 \alpha_6` in
    several root lattice realizations::

        sage: alpha = R.root_space().simple_roots()
        sage: x = alpha[2] + 1/2 * alpha[5]

        sage: L = R.weight_space()
        sage: L(x)
        -Lambda[1] + 2*Lambda[2] - Lambda[3] - 1/2*Lambda[4] + Lambda[5] - 1/2*Lambda[6]

        sage: L = R.ambient_space()
        sage: L(x)
        (0, 1, -1, 0, 1/2, -1/2, 0, 0)

    Of course, one can't embed the root space in the weight lattice::

        sage: L = R.weight_lattice()
        sage: L(x)
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= alpha[2] + 1/2*alpha[5]) an element of self (=Weight lattice of the Root system of type ['A', 7])

    If `K_1` is a subring of `K_2`, then one could in theory have
    an embedding from the root space over `K_1` to any root
    lattice realization over `K_2`; this is not implemented::

        sage: K1 = QQ
        sage: K2 = QQ['q']
        sage: L = R.weight_space(K2)

        sage: alpha = R.root_space(K2).simple_roots()
        sage: L(alpha[1])
        2*Lambda[1] - Lambda[2]

        sage: alpha = R.root_space(K1).simple_roots()
        sage: L(alpha[1])
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= alpha[1]) an element of self (=Weight space over the Univariate Polynomial Ring in q over Rational Field of the Root system of type ['A', 7])

    By a slight abuse, the embedding of the root lattice is not actually
    required to be faithful. Typically for an affine root system, the
    null root of the root lattice is killed in the non extended weight
    lattice::

        sage: R = RootSystem(["A", 3, 1])
        sage: delta = R.root_lattice().null_root()
        sage: L = R.weight_lattice()
        sage: L(delta)
        0

    TESTS::

        sage: TestSuite(L).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.root_lattice_realizations import RootLatticeRealizations
            sage: RootLatticeRealizations(QQ).super_categories()
            [Category of modules with basis over Rational Field]
        """
        return [ModulesWithBasis(self.base_ring())]

    class ParentMethods:

        def __init_extra__(self):
            r"""
            Register the embedding of the root lattice into ``self``.

            Also registers the embedding of the root space over the same
            base field `K` into ``self`` if `K` is not `\ZZ`.

            EXAMPLES:

            We embed the simple root `\alpha_2` of the root lattice in
            the weight lattice::

                sage: R = RootSystem(["A",3])
                sage: alpha = R.root_lattice().simple_roots()
                sage: L = R.weight_lattice()
                sage: L(alpha[2])
                -Lambda[1] + 2*Lambda[2] - Lambda[3]

            .. NOTE::

                More examples are given in :class:`RootLatticeRealizations`;
                The embeddings are systematically tested in
                :meth:`_test_root_lattice_realization`.
            """
            from root_space import RootSpace
            K = self.base_ring()
            # If self is the root lattice or the root space, we don't want
            # to register its trivial embedding into itself. This builds
            # the domains from which we want to register an embedding.
            domains = []
            if not isinstance(self, RootSpace) or K is not ZZ:
                domains.append(self.root_system.root_lattice())
            if not isinstance(self, RootSpace):
                domains.append(self.root_system.root_space(K))
            # Build and register the embeddings
            for domain in domains:
                domain.module_morphism(self.simple_root,
                                       codomain = self
                                       ).register_as_coercion()
            if self.cartan_type().is_affine():
                self._to_classical.register_as_conversion()

        def cartan_type(self):
            """
            EXAMPLES::

                sage: r = RootSystem(['A',4]).root_space()
                sage: r.cartan_type()
                ['A', 4]
            """
            return self.root_system.cartan_type()

        def index_set(self):
            """
            EXAMPLES::

                sage: r = RootSystem(['A',4]).root_space()
                sage: r.index_set()
                (1, 2, 3, 4)
            """
            return self.root_system.index_set()

        def dynkin_diagram(self):
            """
            EXAMPLES::

                sage: r = RootSystem(['A',4]).root_space()
                sage: r.dynkin_diagram()
                O---O---O---O
                1   2   3   4
                A4
            """
            return self.root_system.dynkin_diagram()

        def _name_string_helper(self, name, capitalize=True, base_ring=True, type=True, prefix=""):
            """
            EXAMPLES::

                sage: r = RootSystem(['A',4]).root_space()
                sage: r._name_string_helper("root")
                "Root space over the Rational Field of the Root system of type ['A', 4]"
                sage: r._name_string_helper("root", base_ring=False)
                "Root space of the Root system of type ['A', 4]"
                sage: r._name_string_helper("root", base_ring=False, type=False)
                'Root space'
                sage: r._name_string_helper("root", capitalize=False, base_ring=False, type=False)
                'root space'

                sage: r = RootSystem(['A',4]).coroot_space()
                sage: r._name_string_helper("weight", prefix="extended ")
                "Extended coweight space over the Rational Field of the Root system of type ['A', 4]"
            """
            s = prefix
            if self.root_system.dual_side:
                s += "co"

            s += name + " "

            if self.base_ring() == ZZ:
                s += "lattice "
            else:
                s += "space "
                if base_ring:
                    s += "over the %s "%self.base_ring()

            if type:
                s += "of the "
                if self.root_system.dual_side:
                    s += repr(self.root_system.dual)
                else:
                    s += repr(self.root_system)

            if capitalize:
                s = s[:1].upper() + s[1:]


            return s.strip()

        ##########################################################################
        # checks
        ##########################################################################

        def _test_root_lattice_realization(self, **options):
            """
            Runs sanity checks on this root lattice realization

            - embedding of the root lattice
            - embedding of the root space over the same base ring
            - scalar products between simple roots and simple coroots
            - ...

            .. seealso:: :class:`TestSuite`

            EXAMPLES::

                sage: RootSystem(['A',3]).root_lattice()._test_root_lattice_realization()
            """
            tester = self._tester(**options)
            alpha = self.simple_roots()
            alphacheck = self.simple_coroots()
            R = self.base_ring()
            tester.assertEqual(alpha     .keys(), self.index_set())
            tester.assertEqual(alphacheck.keys(), self.index_set())

            # Check the consistency between simple_root and simple_roots
            for i in self.index_set():
                tester.assertEqual(self.simple_root(i), alpha[i])

            # Check the embeddings from the root lattice and the root space over the same base ring
            root_lattice = self.root_system.root_lattice()
            root_space   = self.root_system.root_space  (R)
            tester.assert_(self.coerce_map_from(root_lattice) is not None)
            tester.assert_(self.coerce_map_from(root_space  ) is not None)
            for i in self.index_set():
                # This embedding maps simple roots to simple roots
                tester.assertEqual(self(root_lattice.simple_root(i)), alpha[i])
                tester.assertEqual(self(root_space   .simple_root(i)), alpha[i])

            # Check that the scalar products match with the Dynkin diagram
            dynkin_diagram = self.dynkin_diagram()
            for i in self.index_set():
                for j in self.index_set():
                    tester.assertEqual(alpha[j].scalar(alphacheck[i]), R(dynkin_diagram[i,j]))

            # Check associated_coroot, if it is implemented
            if not isinstance(self.element_class.associated_coroot, AbstractMethod):
                for i in self.index_set():
                    tester.assertEqual(alpha[i].associated_coroot(), alphacheck[i])

            if self.cartan_type().is_affine():
                # Check that the null root is orthogonal to all coroots
                # and similarly for the null coroot
                nullroot = self.null_root()
                nullcoroot = self.null_coroot()
                special_node = self.cartan_type().special_node()
                for i in alpha.keys():
                    tester.assert_(nullroot.scalar(alphacheck[i]).is_zero())
                    tester.assert_(alpha[i].scalar(nullcoroot).is_zero())
                # Check the projection on the classical space
                classical = self.classical()
                alpha_classical = classical.alpha()
                for i in alpha.keys():
                    if i != special_node or self.cartan_type().is_untwisted_affine():
                        tester.assertEqual(classical(alpha[i]), alpha_classical[i])

            # Todo: add tests of highest root, roots, has_descent, ...

        ##########################################################################
        # highest root
        ##########################################################################

        @cached_method
        def highest_root(self):
            """
            Returns the highest root (for an irreducible finite root system)

            EXAMPLES::

                sage: RootSystem(['A',4]).ambient_space().highest_root()
                (1, 0, 0, 0, -1)

                sage: RootSystem(['E',6]).weight_space().highest_root()
                Lambda[2]

            """
            if not self.root_system.is_finite():
                raise ValueError, "The root system of %s is not of finite Cartan type"%self
            if not self.root_system.is_irreducible():
                raise ValueError, "The root system of %s is reducible"%self
            return self.a_long_simple_root().to_dominant_chamber()

        @cached_method
        def a_long_simple_root(self):
            """
            Returns a long simple root, corresponding to the highest outgoing edge
            in the Dynkin diagram.

            Caveat: this may be break in affine type `A_{2n}^{(2)}`

            Caveat: meaningful/broken for non irreducible?

            TODO: implement CartanType.nodes_by_length as in
            MuPAD-Combinat (using CartanType.symmetrizer), and use it
            here.

            TESTS::

                sage: X=RootSystem(['A',1]).weight_space()
                sage: X.a_long_simple_root()
                2*Lambda[1]
                sage: X=RootSystem(['A',5]).weight_space()
                sage: X.a_long_simple_root()
                2*Lambda[1] - Lambda[2]
            """
            if self.dynkin_diagram().rank() == 1:
                return self.simple_roots()[self.index_set()[0]]
            longest=self.dynkin_diagram().edge_iterator().next()
            for j in self.dynkin_diagram().edge_iterator():
                if j[2]>longest[2]:
                    longest=j
            return self.simple_roots()[longest[0]]


        ##########################################################################
        # simple roots
        ##########################################################################

        @abstract_method
        def simple_root(self, i):
            """
            Returns the `i^{th}` simple root.

            This should be overridden by any subclass, and typically
            implemented as a cached method for efficiency.

            EXAMPLES::

                sage: r = RootSystem(["A",3]).root_lattice()
                sage: r.simple_root(1)
                alpha[1]

            TESTS::

                sage: super(sage.combinat.root_system.root_space.RootSpace, r).simple_root(1)
                Traceback (most recent call last):
                ...
                NotImplementedError: <abstract method simple_root at ...>
            """

        @cached_method
        def simple_roots(self):
            r"""
            Returns the family `(\alpha_i)_{i\in I}` of the simple roots.

            EXAMPLES::

                sage: alpha = RootSystem(["A",3]).root_lattice().simple_roots()
                sage: [alpha[i] for i in [1,2,3]]
                [alpha[1], alpha[2], alpha[3]]
            """
            if not hasattr(self,"_simple_roots"):
                self._simple_roots = Family(self.index_set(), self.simple_root)
                # Should we use rename to set a nice name for this family?
                # self._simple_roots.rename("alpha")
                # This break some doctests
            return self._simple_roots

        @cached_method
        def alpha(self):
            r"""
            Returns the family `(\alpha_i)_{i\in I}` of the simple roots,
            with the extra feature that, for simple irreducible root
            systems, `\alpha_0` yields the opposite of the highest root.

            EXAMPLES::

                sage: alpha = RootSystem(["A",2]).root_lattice().alpha()
                sage: alpha[1]
                alpha[1]
                sage: alpha[0]
                -alpha[1] - alpha[2]

            """
            if self.root_system.is_finite() and self.root_system.is_irreducible():
                return Family(self.index_set(), self.simple_root, \
                              hidden_keys = [0], hidden_function = lambda i: - self.highest_root())
            else:
                return self.simple_roots()

        ##########################################################################
        # roots
        ##########################################################################

        def roots(self):
            """
            Returns the roots of self.

            EXAMPLES::

                sage: RootSystem(['A',2]).ambient_lattice().roots()
                [(1, -1, 0), (1, 0, -1), (0, 1, -1), (-1, 1, 0), (-1, 0, 1), (0, -1, 1)]


            This matches with http://en.wikipedia.org/wiki/Root_systems::

                sage: for T in CartanType.samples(finite = True, crystallographic = True):
                ...       print "%s %3s %3s"%(T, len(RootSystem(T).root_lattice().roots()), len(RootSystem(T).weight_lattice().roots()))
                ['A', 1]   2   2
                ['A', 5]  30  30
                ['B', 1]   2   2
                ['B', 5]  50  50
                ['C', 1]   2   2
                ['C', 5]  50  50
                ['D', 2]   4   4
                ['D', 3]  12  12
                ['D', 5]  40  40
                ['E', 6]  72  72
                ['E', 7] 126 126
                ['E', 8] 240 240
                ['F', 4]  48  48
                ['G', 2]  12  12

            .. todo:: the result should be an enumerated set, and handle infinite root systems
            """
            return list(self.positive_roots()) + list(self.negative_roots())

        def positive_roots(self):
            r"""
            Returns the positive roots of self.

            EXAMPLES::

                sage: L = RootSystem(['A',3]).root_lattice()
                sage: sorted(L.positive_roots())
                [alpha[1], alpha[1] + alpha[2], alpha[1] + alpha[2] + alpha[3], alpha[2], alpha[2] + alpha[3], alpha[3]]

            Algorithm: generate them from the simple roots by applying
            successive reflections toward the positive chamber.
            """
            if not self.cartan_type().is_finite():
                raise NotImplementedError, "Only implemented for finite Cartan type"
            return TransitiveIdealGraded(attrcall('pred'), self.simple_roots())

        @cached_method
        def positive_roots_by_height(self, increasing = True):
            r"""
            Returns a list of positive roots in increasing order by height.

            If ``increasing`` is False, returns them in decreasing order.

            .. warning::

                Returns an error unless the Cartan type is finite.

            EXAMPLES::

                sage: RootSystem(['C',2]).root_lattice().positive_roots_by_height()
                [alpha[1], alpha[2], alpha[1] + alpha[2], 2*alpha[1] + alpha[2]]
                sage: RootSystem(['C',2]).root_lattice().positive_roots_by_height(increasing = False)
                [2*alpha[1] + alpha[2], alpha[1] + alpha[2], alpha[1], alpha[2]]
                sage: RootSystem(['A',2,1]).root_lattice().positive_roots_by_height()
                Traceback (most recent call last):
                ...
                NotImplementedError: Only implemented for finite Cartan type

            """

            if not self.cartan_type().is_finite():
                raise NotImplementedError, "Only implemented for finite Cartan type"
            ranks = self.root_poset().level_sets()
            if not increasing:
                ranks.reverse()
            roots = []
            for x in ranks:
                roots += x
            return [x.element for x in roots]

        @cached_method
        def positive_roots_parabolic(self, index_set = None):
            r"""
            Returns the set of positive roots for the parabolic subsystem with Dynkin node set ``index_set``.

            INPUT:

            - ``index_set`` -- (default:None) the Dynkin node set of the parabolic subsystem. It should be a tuple. The default value implies the entire Dynkin node set

            EXAMPLES::

                sage: lattice =  RootSystem(['A',3]).root_lattice()
                sage: PhiP = lattice.positive_roots_parabolic((1,3))
                sage: [x for x in PhiP]
                [alpha[1], alpha[3]]
                sage: PhiP = lattice.positive_roots_parabolic((2,3))
                sage: [x for x in PhiP]
                [alpha[2], alpha[3], alpha[2] + alpha[3]]
                sage: PhiP = lattice.positive_roots_parabolic()
                sage: [x for x in PhiP]
                [alpha[1], alpha[2], alpha[3], alpha[1] + alpha[2], alpha[2] + alpha[3], alpha[1] + alpha[2] + alpha[3]]

            warning::

                This returns an error if the cartan type is not finite.

            """
            if not self.cartan_type().is_finite():
                raise NotImplementedError, "Only implemented for finite Cartan type"
            if index_set is None:
                index_set = tuple(self.cartan_type().index_set())

            def parabolic_covers(alpha):
                return [x for x in alpha.pred() if x.is_parabolic_root(index_set)]

            generators = [x for x in self.simple_roots() if x.is_parabolic_root(index_set)]
            return TransitiveIdealGraded(parabolic_covers, generators)

        @cached_method
        def positive_roots_nonparabolic(self, index_set = None):
            r"""
            Returns the set of positive roots outside the parabolic subsystem with Dynkin node set ``index_set``.

            INPUT:

            - ``index_set`` -- (default:None) the Dynkin node set of the parabolic subsystem. It should be a tuple. The default value implies the entire Dynkin node set

            EXAMPLES::

                sage: lattice =  RootSystem(['A',3]).root_lattice()
                sage: lattice.positive_roots_nonparabolic((1,3))
                [alpha[2], alpha[1] + alpha[2], alpha[2] + alpha[3], alpha[1] + alpha[2] + alpha[3]]
                sage: lattice.positive_roots_nonparabolic((2,3))
                [alpha[1], alpha[1] + alpha[2], alpha[1] + alpha[2] + alpha[3]]
                sage: lattice.positive_roots_nonparabolic()
                []
                sage: lattice.positive_roots_nonparabolic((1,2,3))
                []

            warning::

                This returns an error if the cartan type is not finite.

            """
            if not self.cartan_type().is_finite():
                raise NotImplementedError, "Only implemented for finite Cartan type"
            if index_set is None:
                index_set = tuple(self.cartan_type().index_set())
            return [x for x in self.positive_roots() if not x.is_parabolic_root(index_set)]

        @cached_method
        def positive_roots_nonparabolic_sum(self, index_set = None):
            r"""
            Returns the sum of positive roots outside the parabolic subsystem with Dynkin node set ``index_set``.

            INPUT:

            - ``index_set`` -- (default:None) the Dynkin node set of the parabolic subsystem. It should be a tuple. The default value implies the entire Dynkin node set

            EXAMPLES::

                sage: lattice =  RootSystem(['A',3]).root_lattice()
                sage: lattice.positive_roots_nonparabolic_sum((1,3))
                2*alpha[1] + 4*alpha[2] + 2*alpha[3]
                sage: lattice.positive_roots_nonparabolic_sum((2,3))
                3*alpha[1] + 2*alpha[2] + alpha[3]
                sage: lattice.positive_roots_nonparabolic_sum(())
                3*alpha[1] + 4*alpha[2] + 3*alpha[3]
                sage: lattice.positive_roots_nonparabolic_sum()
                0
                sage: lattice.positive_roots_nonparabolic_sum((1,2,3))
                0

            warning::

                This returns an error if the cartan type is not finite.

            """

            if not self.cartan_type().is_finite():
                raise ValueError, "Cartan type %s is not finite"%(self.cartan_type())
            if index_set is None or index_set == tuple(self.cartan_type().index_set()):
                return self.zero()
            return sum(self.positive_roots_nonparabolic(index_set))

        def root_poset(self, restricted=False, facade=False):
            r"""
            Returns the (restricted) root poset associated to ``self``.

            The elements are given by the positive roots (resp. non-simple, positive roots), and
            `\alpha \leq \beta` iff `\beta - \alpha` is a non-negative linear combination of simple roots.

            INPUT:

            - ``restricted`` -- (default:False) if True, only non-simple roots are considered.
            - ``facade`` -- (default:False) passes facade option to the poset generator.

            EXAMPLES::

                sage: Phi = RootSystem(['A',1]).root_poset(); Phi
                Finite poset containing 1 elements
                sage: Phi.cover_relations()
                []

                sage: Phi = RootSystem(['A',2]).root_poset(); Phi
                Finite poset containing 3 elements
                sage: Phi.cover_relations()
                [[alpha[1], alpha[1] + alpha[2]], [alpha[2], alpha[1] + alpha[2]]]

                sage: Phi = RootSystem(['A',3]).root_poset(restricted=True); Phi
                Finite poset containing 3 elements
                sage: Phi.cover_relations()
                [[alpha[1] + alpha[2], alpha[1] + alpha[2] + alpha[3]], [alpha[2] + alpha[3], alpha[1] + alpha[2] + alpha[3]]]

                sage: Phi = RootSystem(['B',2]).root_poset(); Phi
                Finite poset containing 4 elements
                sage: Phi.cover_relations()
                [[alpha[1], alpha[1] + alpha[2]], [alpha[2], alpha[1] + alpha[2]], [alpha[1] + alpha[2], alpha[1] + 2*alpha[2]]]
            """
            from sage.combinat.posets.posets import Poset
            rels = []
            dim = self.dimension()
            pos_roots = set(self.positive_roots())
            simple_roots = self.simple_roots()
            if restricted:
                pos_roots = [ beta for beta in pos_roots if beta not in simple_roots ]
            for root in pos_roots:
                for i in range(1,dim+1):
                    root_cover = root + simple_roots[i]
                    if root_cover in pos_roots:
                        rels.append((root,root_cover))
            return Poset((pos_roots,rels),cover_relations=True,facade=facade)

        def almost_positive_roots(self):
            r"""
            Returns the almost positive roots of ``self``

            These are the positive roots together with the simple negative roots.

            .. seealso:: :meth:`almost_positive_root_decomposition`, :meth:`tau_plus_minus`

            EXAMPLES::

                sage: L = RootSystem(['A',2]).root_lattice()
                sage: L.almost_positive_roots()
                [-alpha[1], alpha[1], alpha[1] + alpha[2], -alpha[2], alpha[2]]
            """
            if not self.cartan_type().is_finite():
                raise ValueError, "%s is not a finite Cartan type"%(self.cartan_type())
            return sorted([ -beta for beta in self.simple_roots() ] + list(self.positive_roots()))

        def negative_roots(self):
            r"""
            Returns the negative roots of self.

            EXAMPLES::

                sage: L = RootSystem(['A', 2]).weight_lattice()
                sage: sorted(L.negative_roots())
                [-2*Lambda[1] + Lambda[2], -Lambda[1] - Lambda[2], Lambda[1] - 2*Lambda[2]]

            Algorithm: negate the positive roots

            """
            if not self.cartan_type().is_finite():
                raise ValueError, "%s is not a finite Cartan type"%(self.cartan_type())
            from sage.combinat.combinat import MapCombinatorialClass
            return MapCombinatorialClass(self.positive_roots(), attrcall('__neg__'), "The negative roots of %s"%self)
            # Todo: use this instead once TransitiveIdeal will be a proper enumerated set
            #return self.positive_roots().map(attrcall('__negate__'))

        ##########################################################################
        # coroots
        ##########################################################################

        def coroot_lattice(self):
            """
            Returns the coroot lattice.

            EXAMPLES::

                sage: RootSystem(['A',2]).root_lattice().coroot_lattice()
                Coroot lattice of the Root system of type ['A', 2]

            """
            return self.root_system.coroot_lattice()

        def coroot_space(self, base_ring = QQ):
            """
            Returns the coroot space over ``base_ring``

            INPUT:

            - ``base_ring`` -- a ring (default: `\QQ`)

            EXAMPLES::

                sage: RootSystem(['A',2]).root_lattice().coroot_space()
                Coroot space over the Rational Field of the Root system of type ['A', 2]

                sage: RootSystem(['A',2]).root_lattice().coroot_space(QQ['q'])
                Coroot space over the Univariate Polynomial Ring in q over Rational Field of the Root system of type ['A', 2]

            """
            return self.root_system.coroot_space(base_ring = base_ring)


        def simple_coroot(self, i):
            """
            Returns the `i^{th}` simple coroot.

            EXAMPLES::

                sage: RootSystem(['A',2]).root_lattice().simple_coroot(1)
                alphacheck[1]
            """
            return self.coroot_lattice().simple_root(i)

        @cached_method
        def simple_coroots(self):
            r"""
            Returns the family `( \alpha^\vee_i)_{i\in I}` of the simple coroots.

            EXAMPLES::

                sage: alphacheck = RootSystem(['A',3]).root_lattice().simple_coroots()
                sage: [alphacheck[i] for i in [1, 2, 3]]
                [alphacheck[1], alphacheck[2], alphacheck[3]]

            """
            if not hasattr(self,"cache_simple_coroots"):
                self.cache_simple_coroots = Family(self.index_set(), self.simple_coroot)
                # Should we use rename to set a nice name for this family?
                # self.cache_simple_coroots.rename("alphacheck")
                # break some doctests
            return self.cache_simple_coroots

        def alphacheck(self):
            r"""
            Returns the family `( \alpha^\vee_i)_{i\in I}` of the simple
            coroots, with the extra feature that,  for simple irreducible
            root systems, `\alpha^\vee_0` yields the coroot associated to
            the opposite of the highest root (caveat: for non simply laced
            root systems, this is not the opposite of the highest coroot!)

            EXAMPLES::

                sage: alphacheck = RootSystem(["A",2]).ambient_space().alphacheck()
                sage: alphacheck
                Finite family {1: (1, -1, 0), 2: (0, 1, -1)}

            Here is now `\alpha^\vee_0`:

                (-1, 0, 1)

            .. todo:: add a non simply laced example

            Finaly, here is an affine example::

                sage: RootSystem(["A",2,1]).weight_space().alphacheck()
                Finite family {0: alphacheck[0], 1: alphacheck[1], 2: alphacheck[2]}

                sage: RootSystem(["A",3]).ambient_space().alphacheck()
                Finite family {1: (1, -1, 0, 0), 2: (0, 1, -1, 0), 3: (0, 0, 1, -1)}

            """
            if self.root_system.is_finite() and self.root_system.is_irreducible():
                return Family(self.index_set(), self.simple_coroot, \
                              hidden_keys = [0], hidden_function = lambda i: - self.cohighest_root())
            else:
                return self.simple_coroots()

        @cached_method
        def cohighest_root(self):
            """
            Returns the associated coroot of the highest root.

            .. note:: this is usually not the highest coroot.

            EXAMPLES::

                sage: RootSystem(['A', 3]).ambient_space().cohighest_root()
                (1, 0, 0, -1)
            """
            return self.highest_root().associated_coroot()

        ##########################################################################
        # null_root
        ##########################################################################

        @cached_method
        def null_root(self):
            """
            Returns the null root of self. The null root is the smallest
            non trivial positive root which is orthogonal to all simple
            coroots. It exists for any affine root system.

            EXAMPLES::

                sage: RootSystem(['C',2,1]).root_lattice().null_root()
                alpha[0] + 2*alpha[1] + alpha[2]
                sage: RootSystem(['D',4,1]).root_lattice().null_root()
                alpha[0] + alpha[1] + 2*alpha[2] + alpha[3] + alpha[4]
                sage: RootSystem(['F',4,1]).root_lattice().null_root()
                alpha[0] + 2*alpha[1] + 3*alpha[2] + 4*alpha[3] + 2*alpha[4]
            """
            if self.cartan_type().is_affine():
                coef = self.cartan_type().a()
                return sum(coef[k]*self.simple_roots()[k] for k in coef.keys())

        ##########################################################################
        # null_coroot (Also called CanonicalCentralElement)
        ##########################################################################

        @cached_method
        def null_coroot(self):
            """
            Returns the null coroot of self.

            The null coroot is the smallest non trivial positive
            coroot which is orthogonal to all simple roots. It exists
            for any affine root system.

            EXAMPLES::

                sage: RootSystem(['C',2,1]).root_lattice().null_coroot()
                alphacheck[0] + alphacheck[1] + alphacheck[2]
                sage: RootSystem(['D',4,1]).root_lattice().null_coroot()
                alphacheck[0] + alphacheck[1] + 2*alphacheck[2] + alphacheck[3] + alphacheck[4]
                sage: RootSystem(['F',4,1]).root_lattice().null_coroot()
                alphacheck[0] + 2*alphacheck[1] + 3*alphacheck[2] + 2*alphacheck[3] + alphacheck[4]
            """
            if not self.cartan_type().is_affine():
                raise ValueError, "%s is not an affine Cartan type"%(self.cartan_type())
            coef = self.cartan_type().acheck()
            return sum(coef[k]*self.simple_coroots()[k] for k in coef.keys())


        ##########################################################################
        # fundamental weights
        ##########################################################################

        def fundamental_weights_from_simple_roots(self):
            r"""
            Return the fundamental weights.

            This is computed from the simple roots by using the
            inverse of the Cartan matrix. This method is therefore
            only valid for finite types and if this realization of the
            root lattice is large enough to contain them.

            EXAMPLES:

            In the root space, we retrieve the inverse of the Cartan matrix::

                sage: L = RootSystem(["B",3]).root_space()
                sage: L.fundamental_weights_from_simple_roots()
                Finite family {1:     alpha[1] +   alpha[2] +     alpha[3],
                               2:     alpha[1] + 2*alpha[2] +   2*alpha[3],
                               3: 1/2*alpha[1] +   alpha[2] + 3/2*alpha[3]}
                sage: ~L.cartan_type().cartan_matrix()
                [  1   1 1/2]
                [  1   2   1]
                [  1   2 3/2]

            In the weight lattice and the ambient space, we retrieve
            the fundamental weights::

                sage: L = RootSystem(["B",3]).weight_lattice()
                sage: L.fundamental_weights_from_simple_roots()
                Finite family {1: Lambda[1], 2: Lambda[2], 3: Lambda[3]}

                sage: L = RootSystem(["B",3]).ambient_space()
                sage: L.fundamental_weights()
                Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1/2, 1/2, 1/2)}
                sage: L.fundamental_weights_from_simple_roots()
                Finite family {1: (1, 0, 0), 2: (1, 1, 0), 3: (1/2, 1/2, 1/2)}

            However the fundamental weights do not belong to the root
            lattice::

                sage: L = RootSystem(["B",3]).root_lattice()
                sage: L.fundamental_weights_from_simple_roots()
                Traceback (most recent call last):
                ...
                ValueError: The fundamental weights do not live in this realization of the root lattice

            Beware of the usual `GL_n` vs `SL_n` catch in type `A`::

                sage: L = RootSystem(["A",3]).ambient_space()
                sage: L.fundamental_weights()
                Finite family {1: (1, 0, 0, 0), 2: (1, 1, 0, 0), 3: (1, 1, 1, 0)}
                sage: L.fundamental_weights_from_simple_roots()
                Finite family {1: (3/4, -1/4, -1/4, -1/4), 2: (1/2, 1/2, -1/2, -1/2), 3: (1/4, 1/4, 1/4, -3/4)}

                sage: L = RootSystem(["A",3]).ambient_lattice()
                sage: L.fundamental_weights_from_simple_roots()
                Traceback (most recent call last):
                ...
                ValueError: The fundamental weights do not live in this realization of the root lattice
            """
            # We first scale the inverse of the Cartan matrix to be
            # with integer coefficients; then the linear combination
            # of the simple roots is guaranteed to live in this space,
            # and then we rely on division by d to fail gracefuly.
            M = self.cartan_type().cartan_matrix()
            d = M.det()
            if not d:
                raise TypeError("The Cartan matrix is not invertible")
            M = d*~M
            fundamental_weights = [self.linear_combination(zip(self.simple_roots(), column))
                                   for column in M.columns()]
            try:
                fundamental_weights = [x/d for x in fundamental_weights]
            except ValueError:
                raise ValueError("The fundamental weights do not live in this realization of the root lattice")
            return Family(dict(zip(self.index_set(),fundamental_weights)))


        ##########################################################################
        # reflections
        ##########################################################################

        def reflection(self, root, coroot=None):
            """
            Returns the reflection along the root, and across the
            hyperplane define by coroot, as a function from
            self to self.

            EXAMPLES::

                sage: space = RootSystem(['A',2]).weight_lattice()
                sage: x=space.simple_roots()[1]
                sage: y=space.simple_coroots()[1]
                sage: s = space.reflection(x,y)
                sage: x
                2*Lambda[1] - Lambda[2]
                sage: s(x)
                -2*Lambda[1] + Lambda[2]
                sage: s(-x)
                2*Lambda[1] - Lambda[2]
            """
            if coroot is None:
                coroot = root.associated_coroot()
            return lambda v: v - v.scalar(coroot) * root

        @cached_method
        def simple_reflection(self, i):
            """
            Returns the `i^{th}` simple reflection, as a function from
            self to self.

            INPUT:

            - ``i`` - i is in self's index set

            EXAMPLES::

                sage: space = RootSystem(['A',2]).ambient_lattice()
                sage: s = space.simple_reflection(1)
                sage: x = space.simple_roots()[1]
                sage: x
                (1, -1, 0)
                sage: s(x)
                (-1, 1, 0)
            """
            return self.reflection(self.simple_root(i), self.simple_coroot(i))

        @cached_method
        def simple_reflections(self):
            """
            Returns the family `(s_i)_{i\in I}` of the simple reflections
            of this root system.

            EXAMPLES::

                sage: r = RootSystem(["A", 2]).root_lattice()
                sage: s = r.simple_reflections()
                sage: s[1]( r.simple_root(1) )
                -alpha[1]

            TEST::

                sage: s
                simple reflections
            """
            res =  self.alpha().zip(self.reflection, self.alphacheck())
            # Should we use rename to set a nice name for this family?
            res.rename("simple reflections")
            return res

        s = simple_reflections

        ##########################################################################
        # projections
        ##########################################################################

        def projection(self, root, coroot=None, to_negative=True):
            r"""
            Returns the projection along the root, and across the
            hyperplane define by coroot, as a function `\pi` from self to
            self. `\pi` is a half-linear map which stabilizes the negative
            half space, and acts by reflection on the positive half space.

            If to_negative is False, then this project onto the positive
            half space instead.

            EXAMPLES::

                sage: space = RootSystem(['A',2]).weight_lattice()
                sage: x=space.simple_roots()[1]
                sage: y=space.simple_coroots()[1]
                sage: pi = space.projection(x,y)
                sage: x
                2*Lambda[1] - Lambda[2]
                sage: pi(x)
                -2*Lambda[1] + Lambda[2]
                sage: pi(-x)
                -2*Lambda[1] + Lambda[2]
                sage: pi = space.projection(x,y,False)
                sage: pi(-x)
                2*Lambda[1] - Lambda[2]
            """
            if coroot is None:
                coroot = root.associated_coroot()

            return lambda v: v - v.scalar(coroot) * root if ((v.scalar(coroot) > 0) == to_negative) else v

        @cached_method
        def simple_projection(self, i, to_negative=True):
            """
            Returns the projection along the `i^{th}` simple root, and across the
            hyperplane define by the `i^{th}` simple coroot, as a function from
            self to self.

            INPUT:

            - ``i`` - i is in self's index set

            EXAMPLES::

                sage: space = RootSystem(['A',2]).weight_lattice()
                sage: x = space.simple_roots()[1]
                sage: pi = space.simple_projection(1)
                sage: x
                2*Lambda[1] - Lambda[2]
                sage: pi(x)
                -2*Lambda[1] + Lambda[2]
                sage: pi(-x)
                -2*Lambda[1] + Lambda[2]
                sage: pi = space.simple_projection(1,False)
                sage: pi(-x)
                2*Lambda[1] - Lambda[2]
            """
            return self.projection(self.simple_root(i), self.simple_coroot(i), to_negative)

        @cached_method
        def simple_projections(self, to_negative=True):
            r"""
            Returns the family `(s_i)_{i\in I}` of the simple projections
            of this root system

            EXAMPLES::

                sage: space = RootSystem(['A',2]).weight_lattice()
                sage: pi = space.simple_projections()
                sage: x = space.simple_roots()
                sage: pi[1](x[2])
                -Lambda[1] + 2*Lambda[2]

            TESTS:
                sage: pi
                pi
            """
            if to_negative is not True:
                raise NotImplementedError, "only implemented when 'to_negative' is True"
            res = self.alpha().zip(self.projection, self.alphacheck())
            # Should this use rename to set a nice name for this family?
            res.rename("pi")
            return res

        @lazy_attribute
        def pi(self):
            r"""
            The simple projections of ``self``

            .. seealso:: :meth:`simple_projections`

            .. warning:: this shortcut is deprecated

            EXAMPLES::

                sage: space = RootSystem(['A',2]).weight_lattice()
                sage: pi = space.pi
                sage: x = space.simple_roots()
                sage: pi[1](x[2])
                -Lambda[1] + 2*Lambda[2]
            """
            # _test_not_implemented_methods apparently evaluates all lazy
            # attributes, which means that we can't use deprecation here!
            # from sage.misc.superseded import deprecation
            # deprecation(trac_number, "The lazy attribute pi is deprecated; please use the simple_projections method.")
            return self.simple_projections()

        ##########################################################################
        # Weyl group
        ##########################################################################

        def weyl_group(self, prefix=None):
            """
            Returns the Weyl group associated to self.

            EXAMPLES::

                sage: RootSystem(['F',4]).ambient_space().weyl_group()
                Weyl Group of type ['F', 4] (as a matrix group acting on the ambient space)
                sage: RootSystem(['F',4]).root_space().weyl_group()
                Weyl Group of type ['F', 4] (as a matrix group acting on the root space)

            """
            from sage.combinat.root_system.weyl_group import WeylGroup
            return WeylGroup(self, prefix=prefix)

        ##########################################################################
        # The piecewise linear involutive operators tau_plus and tau_minus on self,
        # and the orbit decomposition of the almost positive roots
        # by the associated dihedral group
        ##########################################################################

        # TODO: find a better name; at least, this temporary one won't
        # create conflicts
        def tau_epsilon_operator_on_almost_positive_roots(self, J):
            r"""
            The `\tau_\epsilon` operator on almost positive roots

            Given a subset `J` of non adjacent vertices of the Dynkin
            diagram, this constructs the operator on the almost positive
            roots which fixes the negative simple roots `\alpha_i` for `i`
            not in `J`, and acts otherwise by:

            .. math::

                \tau_+( \beta ) = (\prod_{i \in J} s_i) (\beta)

            See Equation (1.2) of [CFZ]_.

            EXAMPLES::

                sage: L = RootSystem(['A',4]).root_lattice()
                sage: tau = L.tau_epsilon_operator_on_almost_positive_roots([1,3])
                sage: alpha = L.simple_roots()

            The action on a negative simple root not in `J`::

                sage: tau(-alpha[2])
                -alpha[2]

            The action on a negative simple root in `J`::

                sage: tau(-alpha[1])
                alpha[1]

            The action on all almost positive roots::

                sage: for root in L.almost_positive_roots():
                ...      print 'tau({:<41}) ='.format(root), tau(root)
                tau(-alpha[1]                                ) = alpha[1]
                tau(alpha[1]                                 ) = -alpha[1]
                tau(alpha[1] + alpha[2]                      ) = alpha[2] + alpha[3]
                tau(alpha[1] + alpha[2] + alpha[3]           ) = alpha[2]
                tau(alpha[1] + alpha[2] + alpha[3] + alpha[4]) = alpha[2] + alpha[3] + alpha[4]
                tau(-alpha[2]                                ) = -alpha[2]
                tau(alpha[2]                                 ) = alpha[1] + alpha[2] + alpha[3]
                tau(alpha[2] + alpha[3]                      ) = alpha[1] + alpha[2]
                tau(alpha[2] + alpha[3] + alpha[4]           ) = alpha[1] + alpha[2] + alpha[3] + alpha[4]
                tau(-alpha[3]                                ) = alpha[3]
                tau(alpha[3]                                 ) = -alpha[3]
                tau(alpha[3] + alpha[4]                      ) = alpha[4]
                tau(-alpha[4]                                ) = -alpha[4]
                tau(alpha[4]                                 ) = alpha[3] + alpha[4]

            This method works on any root lattice realization::

                sage: L = RootSystem(['B',3]).ambient_space()
                sage: tau = L.tau_epsilon_operator_on_almost_positive_roots([1,3])
                sage: for root in L.almost_positive_roots():
                ...      print 'tau({:<41}) ='.format(root), tau(root)
                tau((-1, 1, 0)                               ) = (1, -1, 0)
                tau((1, 0, 0)                                ) = (0, 1, 0)
                tau((1, -1, 0)                               ) = (-1, 1, 0)
                tau((1, 1, 0)                                ) = (1, 1, 0)
                tau((1, 0, -1)                               ) = (0, 1, 1)
                tau((1, 0, 1)                                ) = (0, 1, -1)
                tau((0, -1, 1)                               ) = (0, -1, 1)
                tau((0, 1, 0)                                ) = (1, 0, 0)
                tau((0, 1, -1)                               ) = (1, 0, 1)
                tau((0, 1, 1)                                ) = (1, 0, -1)
                tau((0, 0, -1)                               ) = (0, 0, 1)
                tau((0, 0, 1)                                ) = (0, 0, -1)

            .. seealso:: :meth:`tau_plus_minus`

            REFERENCES:

                .. [CFZ] Chapoton, Fomin, Zelevinsky - Polytopal realizations of generalized associahedra
            """
            W = self.weyl_group()
            t = W.from_reduced_word(J)
            simple_roots = self.simple_roots()
            other_negative_simple_roots = set(-simple_roots[i] for i in self.index_set() if i not in J)
            def tau_epsilon(alpha):
                if alpha in other_negative_simple_roots:
                    return alpha
                else:
                    return t.action(alpha)
            return tau_epsilon

        def tau_plus_minus(self):
            r"""
            Returns the `\tau^+` and `\tau^-` piecewise linear operators on ``self``

            Those operators are induced by the bipartition `\{L,R\}` of
            the simple roots of ``self``, and stabilize the almost
            positive roots. Namely, `\tau_+` fixes the negative simple
            roots `\alpha_i` for `i` in `R`, and acts otherwise by:

            .. math::

                \tau_+( \beta ) = (\prod_{i \in L} s_i) (\beta)

            `\tau_-` acts analogously, with `L` and `R` interchanged.

            Those operators are used to construct the associahedron, a
            polytopal realization of the cluster complex (see
            :class:`Associahedron`).

            .. seealso:: :meth:`tau_epsilon_operator_on_almost_positive_roots`

            EXAMPLES:

            We explore the example of [CFZ1]_ Eq.(1.3)::

                sage: S = RootSystem(['A',2]).root_lattice()
                sage: taup, taum = S.tau_plus_minus()
                sage: for beta in S.almost_positive_roots(): print beta, ",", taup(beta), ",", taum(beta)
                -alpha[1] , alpha[1] , -alpha[1]
                alpha[1] , -alpha[1] , alpha[1] + alpha[2]
                alpha[1] + alpha[2] , alpha[2] , alpha[1]
                -alpha[2] , -alpha[2] , alpha[2]
                alpha[2] , alpha[1] + alpha[2] , -alpha[2]

            REFERENCES:

                .. [CFZ1] Chapoton, Fomin, Zelevinsky - Polytopal realizations of generalized associahedra
            """
            ct = self.cartan_type()
            L,R = ct.index_set_bipartition()
            return self.tau_epsilon_operator_on_almost_positive_roots(L), self.tau_epsilon_operator_on_almost_positive_roots(R)

        def almost_positive_roots_decomposition(self):
            r"""
            Returns the decomposition of the almost positive roots of ``self``

            This is the list of the orbits of the almost positive roots
            under the action of the dihedral group generated by the
            operators `\tau_+` and `\tau_-`.

            .. SEEALSO::

                - :meth:`almost_positive_roots`
                - :meth:`tau_plus_minus`

            EXAMPLES::

                sage: RootSystem(['A',2]).root_lattice().almost_positive_roots_decomposition()
                [[-alpha[1], alpha[1], alpha[1] + alpha[2], alpha[2], -alpha[2]]]

                sage: RootSystem(['B',2]).root_lattice().almost_positive_roots_decomposition()
                [[-alpha[1], alpha[1], alpha[1] + 2*alpha[2]], [-alpha[2], alpha[2], alpha[1] + alpha[2]]]

                sage: RootSystem(['D',4]).root_lattice().almost_positive_roots_decomposition()
                [[-alpha[1], alpha[1], alpha[1] + alpha[2], alpha[2] + alpha[3] + alpha[4]],
                 [-alpha[2], alpha[2], alpha[1] + alpha[2] + alpha[3] + alpha[4], alpha[1] + 2*alpha[2] + alpha[3] + alpha[4]],
                 [-alpha[3], alpha[3], alpha[2] + alpha[3], alpha[1] + alpha[2] + alpha[4]],
                 [-alpha[4], alpha[4], alpha[2] + alpha[4], alpha[1] + alpha[2] + alpha[3]]]

            REFERENCES:

            .. [CFZ2] Chapoton, Fomin, Zelevinsky - Polytopal realizations of
               generalized associahedra
            """
            # TODO: this should use a generic function for computing
            # orbits under the action of a group:
            # def orbits(seeds, operators)
            #     INPUT:
            #     - seeds: a list of elements
            #     - operators: a list of functions
            #
            #     Returns the orbits generated by seeds under the action of the operators
            tau_plus, tau_minus = self.tau_plus_minus()

            I = set(self.index_set())
            Delta = self.simple_roots()
            L, R = self.cartan_type().index_set_bipartition()

            orbits = []
            while I:
                i = I.pop()
                alpha = -self.simple_root(i)
                orbit = [alpha]
                if i in L:
                    plus = False
                    beta = tau_plus(alpha)
                else:
                    plus = True
                    beta = tau_minus(alpha)
                while -beta not in Delta and beta not in orbit:
                    orbit.append(beta)
                    if beta in Delta:
                        j = beta.leading_support()
                        I.discard(j)
                    if plus:
                        beta = tau_plus(beta)
                    else:
                        beta = tau_minus(beta)
                    plus = not plus
                if -beta in Delta:
                    orbit.append(beta)
                orbits.append(orbit)
            return orbits


        ##########################################################################
        # Methods for affine root lattice realizations
        # Should eventually go in an Affine nested class
        ##########################################################################

        @cached_method
        def classical(self):
            """
            Return the corresponding root/weight/ambient lattice/space.

            EXAMPLES::

                sage: RootSystem(["A",4,1]).root_lattice().classical()
                Root lattice of the Root system of type ['A', 4]
                sage: RootSystem(["A",4,1]).weight_lattice().classical()
                Weight lattice of the Root system of type ['A', 4]
                sage: RootSystem(["A",4,1]).ambient_space().classical()
                Ambient space of the Root system of type ['A', 4]
            """
            from root_space import RootSpace
            from weight_space import WeightSpace
            R = self.cartan_type().classical().root_system()
            if isinstance(self, RootSpace):
                return R.root_space(self.base_ring())
            elif isinstance(self, WeightSpace):
                return R.weight_space(self.base_ring())
            else:
                return R.ambient_space(self.base_ring())

        @lazy_attribute
        def _to_classical(self):
            r"""
            The projection onto the classical ambient space.

            EXAMPLES::

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: e = L.basis()
                sage: L._to_classical(e["delta"])
                (0, 0, 0)
                sage: L._to_classical(e["deltacheck"])
                (0, 0, 0)
                sage: L._to_classical(e[0])
                (1, 0, 0)
                sage: L._to_classical(e[1])
                (0, 1, 0)
                sage: L._to_classical(e[2])
                (0, 0, 1)
            """
            return self.module_morphism(self._to_classical_on_basis, codomain = self.classical())

        def _classical_alpha_0(self):
            """
            Return the projection of `\alpha_0` in the classical space.

            This is used e.g. to construct the projections onto the
            classical space.

            EXAMPLES:

            This is the opposite of the highest root in the untwisted case::

                sage: L = RootSystem(["B",3,1]).root_space()
                sage: L._classical_alpha_0()
                -alpha[1] - 2*alpha[2] - 2*alpha[3]
                sage: L._to_classical_on_basis(0)
                -alpha[1] - 2*alpha[2] - 2*alpha[3]
                sage: L.classical().highest_root()
                alpha[1] + 2*alpha[2] + 2*alpha[3]

            But not in the other cases::

                sage: L = RootSystem(CartanType(["B",3,1]).dual()).root_space()
                sage: L._to_classical_on_basis(0)
                -alpha[1] - 2*alpha[2] - alpha[3]
                sage: L.classical().highest_root()
                2*alpha[1] + 2*alpha[2] + alpha[3]
            """
            cartan_type  = self.cartan_type()
            special_node = cartan_type.special_node()
            a = self.cartan_type().col_annihilator()
            classical = self.classical()
            return -classical.sum(a[i] * self.simple_root(i)
                                  for i in self.index_set() if i != special_node) \
                                  / a[special_node]

        ######################################################################
        # Root system plots

        def plot(self,
                 roots="simple",
                 coroots=False,
                 reflection_hyperplanes="simple",
                 fundamental_weights=None,
                 fundamental_chamber=None,
                 alcoves=None,
                 alcove_labels=False,
                 alcove_walk=None,
                 **options):
            r"""
            Return a picture of this root lattice realization.

            INPUT:

            - ``roots`` -- which roots to display, if any.
              Can be one of the following:

              * ``"simple"`` -- The simple roots (the default)
              * ``"classical"`` -- Not yet implemented
              * ``"all"`` -- Only works in the finite case
              * A list or tuple of roots
              * ``False``

            - ``coroots`` -- which coroots to display, if any.
              Can be one of the following:

              * ``"simple"`` -- The simple coroots (the default)
              * ``"classical"`` -- Not yet implemented
              * ``"all"`` -- Only works in the finite case
              * A list or tuple of coroots
              * ``False``

            - ``fundamental_weights`` -- a boolean or ``None`` (default: ``None``)
              whether to display the fundamental weights.
              If ``None``, the fundamental weights are drawn if available.

            - ``reflection_hyperplanes`` -- which reflection
              hyperplanes to display, if any. Can be one of the
              following:

              * ``"simple"`` -- The simple roots
              * ``"classical"`` -- Not yet implemented
              * ``"all"`` -- Only works in the finite case
              * A list or tuple of roots
              * ``False`` (the default)

            - ``fundamental_chamber`` -- whether and how to draw the
              fundamental chamber. Can be one of the following:

              * A boolean -- Set to ``True`` to draw the fundamental
                chamber
              * ``"classical"`` -- Draw the classical fundamental chamber
              * ``None`` -- (the default) The fundamental chamber is
                drawn except in the root lattice where this is not yet
                implemented. For affine types the classical
                fundamental chamber is drawn instead.

            - ``alcoves`` -- one of the following (default: ``True``):

              * A boolean -- Whether to display the alcoves
              * A list of alcoves -- The alcoves to be drawn. Each alcove is
                specified by the coordinates of its center in the root lattice
                (affine type only). Otherwise the alcoves that intersect the
                bounding box are drawn.

            - ``alcove_labels`` -- one of the following (default: ``False``):

              * A boolean -- Whether to display the elements of the Weyl group
                indexing the alcoves. This currently requires to also
                set the ``alcoves`` option.
              * A number `l` -- The label is drawn at level `l` (affine type
                only), which only makes sense if ``affine`` is ``False``.

            - ``bounding_box`` -- a rational number or a list of pairs
              thereof (default: 3)

              Specifies a bounding box, in the coordinate system for
              this plot, in which to plot alcoves and other infinite
              objects. If the bounding box is a number `a`, then the
              bounding box is of the form `[-a,a]` in all directions.
              Beware that there can be some border effects and the
              returned graphic is not necessarily strictly contained
              in the bounding box.

            - ``alcove_walk`` -- an alcove walk or ``None`` (default: ``None``)

              The alcove walk is described by a list (or iterable) of
              vertices of the Dynkin diagram which specifies which
              wall is crossed at each step, starting from the
              fundamental alcove.

            - ``projection`` -- one of the following (default: ``True``):

              * ``True`` -- The default projection for the root
                lattice realization is used.
              * ``False`` -- No projection is used.
              * ``barycentric`` -- A barycentric projection is used.
              * A function -- If a function is specified, it should implement a
                linear (or affine) map taking as input an element of
                this root lattice realization and returning its
                desired coordinates in the plot, as a vector with
                rational coordinates.

            - ``color`` -- a function mapping vertices of the Dynkin
              diagram to colors (default: ``"black"`` for 0,
              ``"blue"`` for 1, ``"red"`` for 2, ``"green"`` for 3)

              This is used to set the color for the simple roots,
              fundamental weights, reflection hyperplanes, alcove
              facets, etc. If the color is ``None``, the object is not
              drawn.

            - ``labels`` -- a boolean (default: ``True``)
              whether to display labels on the simple roots,
              fundamental weights, etc.

            EXAMPLES::

                sage: L = RootSystem(["A",2,1]).ambient_space().plot()

            .. SEEALSO::

                - :meth:`plot_parse_options`
                - :meth:`plot_roots`, :meth:`plot_coroots`
                - :meth:`plot_fundamental_weights`
                - :meth:`plot_fundamental_chamber`
                - :meth:`plot_reflection_hyperplanes`
                - :meth:`plot_alcoves`
                - :meth:`plot_alcove_walk`
            """
            plot_options = self.plot_parse_options(**options)
            G = plot_options.empty()

            if roots:
                G += self.plot_roots(roots, plot_options=plot_options)

            # if coroots is None:
            #    coroot_lattice = self.root_system.coroot_lattice()
            #    if self.has_coerce_map_from(coroot_lattice):
            #        coroots="simple"
            #    else:
            #        coroots=False
            if coroots:
                G += self.plot_coroots(coroots, plot_options=plot_options)

            if fundamental_weights is None:
                fundamental_weights = hasattr(self, "fundamental_weights")
            if fundamental_weights:
                G += self.plot_fundamental_weights(plot_options=plot_options)

            if reflection_hyperplanes:
                G += self.plot_reflection_hyperplanes(reflection_hyperplanes, plot_options=plot_options)

            if alcoves is None:
                alcoves = self.cartan_type().is_affine() and hasattr(self, "fundamental_weights")
            if alcoves:
                G += self.plot_alcoves(alcoves, alcove_labels=alcove_labels, plot_options=plot_options)

            if fundamental_chamber is None:
                if not hasattr(self, "fundamental_weights"):
                    fundamental_chamber = False
                elif self.cartan_type().is_affine():
                    fundamental_chamber = "classical"
                else:
                    fundamental_chamber = True
            if fundamental_chamber:
                G += self.plot_fundamental_chamber(fundamental_chamber, plot_options=plot_options)

            if alcove_walk is not None:
                G += self.plot_alcove_walk(alcove_walk, plot_options=plot_options)

            return plot_options.finalize(G)

        def plot_parse_options(self, **args):
            r"""
            Return an option object to be used for root system plotting.

            EXAMPLES::

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: options = L.plot_parse_options()
                sage: options
                <sage.combinat.root_system.plot.PlotOptions instance at ...>

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting
            """
            if len(args) == 1 and "plot_options" in args:
                return args["plot_options"]
            else:
                return PlotOptions(self, **args)

        def _plot_projection(self, x):
            r"""
            Implement the default projection to be used for plots.

            EXAMPLES:

            By default, this is just the identity::

                sage: L = RootSystem(["B",3]).root_lattice()
                sage: l = L.an_element(); l
                2*alpha[1] + 2*alpha[2] + 3*alpha[3]
                sage: L._plot_projection(l)
                2*alpha[1] + 2*alpha[2] + 3*alpha[3]

            In the ambient space of type `A_2`, this is the
            barycentric projection. In the ambient space of affine
            type this goes through the classical ambient space.

            .. SEEALSO::

                - :meth:`sage.combinat.root_system.type_A.AmbientSpace._plot_projection`
                - :meth:`sage.combinat.root_system.type_affine.AmbientSpace._plot_projection`
                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting
            """
            return x

        @cached_method
        def _plot_projection_barycentric_matrix(self):
            """
            A rational approximation of the matrix for the barycentric projection

            OUTPUT: a matrix with rational coefficients whose column sum is zero

            .. SEE_ALSO::

                - :func:`sage.combinat.root_system.plot.barycentric_projection_matrix`
                - :meth:`_plot_projection_barycentric`

            EXAMPLES::

                sage: RootSystem(["A",0]).ambient_space()._plot_projection_barycentric_matrix()
                []
                sage: m = RootSystem(["A",1]).ambient_space()._plot_projection_barycentric_matrix(); m
                [ 1 -1]
                sage: sum(m.columns())
                (0)
                sage: m = RootSystem(["A",2]).ambient_space()._plot_projection_barycentric_matrix(); m
                [      1/2        -1       1/2]
                [ 989/1142         0 -989/1142]
                sage: sum(m.columns())
                (0, 0)
                sage: m = RootSystem(["A",3]).ambient_space()._plot_projection_barycentric_matrix(); m
                [      1277/1564      -1277/1564               0               0]
                [1009460/2141389        849/1801      -1121/1189               0]
                [            1/3             1/3             1/3              -1]
                sage: sum(m.columns())
                (0, 0, 0)

            """
            from sage.matrix.constructor import matrix
            from sage.symbolic.constants import pi
            m = matrix(QQ, barycentric_projection_matrix(self.dimension()-1, angle=2*pi/3).n(20))
            # We want to guarantee that the sum of the columns of the
            # result is zero. This is close to be the case for the
            # original matrix and for the current rational
            # approximation. We tidy up the work by replacing the
            # first colum by the opposite of the sum of the others.
            if self.dimension()>1: # not needed in the trivial cases
                m.set_column(0, -sum(m[:,1:].columns()))
            m.set_immutable()
            return m

        def _plot_projection_barycentric(self, x):
            r"""
            Implement the barycentric projection to be used for plots.

            It is in fact a rational approximation thereof, but the
            sum of the basis vectors is guaranteed to be mapped to
            zero.

            EXAMPLES::

                sage: L = RootSystem(["A",2]).ambient_space()
                sage: e = L.basis()
                sage: L._plot_projection_barycentric(e[0])
                (1/2, 989/1142)
                sage: L._plot_projection_barycentric(e[1])
                (-1, 0)
                sage: L._plot_projection_barycentric(e[2])
                (1/2, -989/1142)

            .. SEEALSO::

                - :meth:`_plot_projection`, :meth:`plot`
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting
            """
            return self._plot_projection_barycentric_matrix()*vector(x)

        def plot_roots(self, collection="simple", **options):
            r"""
            Plot the (simple/classical) roots of this root lattice.

            INPUT:

            - ``collection`` -- which roots to display
              can be one of the following:

              * ``"simple"`` (the default)
              * ``"classical"``
              * ``"all"``

            - ``**options`` -- Plotting options

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES::

                sage: RootSystem(["B",3]).ambient_space().plot_roots()
                sage: RootSystem(["B",3]).ambient_space().plot_roots("all")

            TESTS::

                sage: list(RootSystem(["A",2]).root_lattice().plot_roots())
                [Arrow from (0.0,0.0) to (1.0,0.0),
                 Text '$\alpha_{1}$' at the point (1.05,0.0),
                 Arrow from (0.0,0.0) to (0.0,1.0),
                 Text '$\alpha_{2}$' at the point (0.0,1.05)]

                sage: list(RootSystem(["A",2]).weight_lattice().plot_roots(labels=False))
                [Arrow from (0.0,0.0) to (2.0,-1.0),
                 Arrow from (0.0,0.0) to (-1.0,2.0)]

                 sage: list(RootSystem(["A",2]).ambient_lattice().plot_roots())
                 [Arrow from (0.0,0.0) to (1.5,0.86...),
                  Text '$\alpha_{1}$' at the point (1.575,0.90...),
                  Arrow from (0.0,0.0) to (-1.5,0.86...),
                  Text '$\alpha_{2}$' at the point (-1.575,0.90...)]

                 sage: list(RootSystem(["B",2]).ambient_space().plot_roots())
                 [Arrow from (0.0,0.0) to (1.0,-1.0),
                  Text '$\alpha_{1}$' at the point (1.05,-1.05),
                  Arrow from (0.0,0.0) to (0.0,1.0),
                  Text '$\alpha_{2}$' at the point (0.0,1.05)]

                sage: list(RootSystem(["A",2]).root_lattice().plot_roots("all"))
                [Arrow from (0.0,0.0) to (1.0,0.0),
                 Text '$\alpha_{1}$' at the point (1.05,0.0),
                 Arrow from (0.0,0.0) to (0.0,1.0),
                 Text '$\alpha_{2}$' at the point (0.0,1.05),
                 Arrow from (0.0,0.0) to (1.0,1.0),
                 Text '$\alpha_{1} + \alpha_{2}$' at the point (1.05,1.05),
                 Arrow from (0.0,0.0) to (-1.0,0.0),
                 Text '$-\alpha_{1}$' at the point (-1.05,0.0),
                 Arrow from (0.0,0.0) to (0.0,-1.0),
                 Text '$-\alpha_{2}$' at the point (0.0,-1.05),
                 Arrow from (0.0,0.0) to (-1.0,-1.0),
                 Text '$-\alpha_{1} - \alpha_{2}$' at the point (-1.05,-1.05)]
            """
            plot_options = self.plot_parse_options(**options)
            root_lattice = self.root_system.root_lattice()
            if collection == "simple":
                roots = root_lattice.simple_roots()
            elif collection == "classical":
                if not self.cartan_type().is_affine():
                    raise ValueError("plotting classical roots only available in affine type")
                raise NotImplementedError("classical roots")
            elif collection == "all":
                assert self.cartan_type().is_finite(), "plotting all roots only available in finite type"
                roots = root_lattice.roots()
            elif isinstance(collection, (list, tuple)):
                roots = collection
            else:
                raise ValueError("Unknown value: %s"%collection)
            roots = Family(roots, self)
            return plot_options.family_of_vectors(roots)

        def plot_coroots(self, collection="simple", **options):
            r"""
            Plot the (simple/classical) coroots of this root lattice.

            INPUT:

            - ``collection`` -- which coroots to display.
              Can be one of the following:

              * ``"simple"`` (the default)
              * ``"classical"``
              * ``"all"``

            - ``**options`` -- Plotting options

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES::

                sage: RootSystem(["B",3]).ambient_space().plot_coroots()

            TESTS::

                 sage: list(RootSystem(["B",2]).ambient_space().plot_coroots())
                 [Arrow from (0.0,0.0) to (1.0,-1.0),
                  Text '$\alpha^\vee_{1}$' at the point (1.05,-1.05),
                  Arrow from (0.0,0.0) to (0.0,2.0),
                  Text '$\alpha^\vee_{2}$' at the point (0.0,2.1)]
            """
            # Functionally speaking, this is duplicated from plot_roots ...
            # Can we avoid that, say by going to the dual space?
            plot_options = self.plot_parse_options(**options)
            coroot_lattice = self.root_system.coroot_lattice()
            if not self.has_coerce_map_from(coroot_lattice):
                raise ValueError("Can't plot the coroots: there is no embedding of the coroot lattice to this space")
            if collection == "simple":
                coroots = coroot_lattice.simple_roots()
            elif collection == "classical":
                if not self.cartan_type().is_affine():
                    raise ValueError("plotting classical coroots only available in affine type")
                raise NotImplementedError("classical coroots")
            elif collection == "all":
                assert self.cartan_type().is_finite(), "plotting all coroots only available in finite type"
                coroots = coroot_lattice.roots()
            elif isinstance(collection, (list, tuple)):
                coroots = collection
            else:
                raise ValueError("Unknown value: %s"%collection)
            coroots = Family(coroots, self)
            return plot_options.family_of_vectors(coroots)

        def plot_fundamental_weights(self, **options):
            r"""
            Plot the fundamental weights of this root lattice.

            INPUT:

            - ``**options`` -- Plotting options

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES::

                sage: RootSystem(["B",3]).ambient_space().plot_fundamental_weights()

            TESTS::

                sage: list(RootSystem(["A",2]).weight_lattice().plot_fundamental_weights())
                [Arrow from (0.0,0.0) to (0.0,1.0),
                 Text '$\Lambda_{2}$' at the point (0.0,1.05),
                 Arrow from (0.0,0.0) to (1.0,0.0),
                 Text '$\Lambda_{1}$' at the point (1.05,0.0)]

                 sage: list(RootSystem(["A",2]).ambient_lattice().plot_fundamental_weights())
                 [Arrow from (0.0,0.0) to (-0.5,0.86...),
                  Text '$\Lambda_{2}$' at the point (-0.525,0.90...),
                  Arrow from (0.0,0.0) to (0.5,0.86...),
                  Text '$\Lambda_{1}$' at the point (0.525,0.90...)]
            """
            plot_options = self.plot_parse_options(**options)
            # We build the family of fundamental weights in this space,
            # indexed by the fundamental weights in the weight lattice.
            #
            # To this end, we don't use the embdding of the weight
            # lattice into self as for the roots or coroots because
            # the ambient space can define the fundamental weights
            # slightly differently (the usual GL_n vs SL_n catch).
            weight_lattice = self.root_system.weight_lattice()
            fundamental_weights = Family(dict(zip(weight_lattice.fundamental_weights(),
                                                  self.fundamental_weights())))
            return plot_options.family_of_vectors(fundamental_weights)

        def plot_reflection_hyperplanes(self, collection="simple", **options):
            r"""
            Plot the simple reflection hyperplanes.

            INPUT:

            - ``collection`` -- which reflection hyperplanes to display.
              Can be one of the following:

              * ``"simple"`` (the default)
              * ``"classical"``
              * ``"all"``

            - ``**options`` -- Plotting options

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES::

                sage: RootSystem(["A",2,1]).ambient_space().plot_reflection_hyperplanes()
                sage: RootSystem(["G",2,1]).ambient_space().plot_reflection_hyperplanes()
                sage: RootSystem(["A",3]).weight_space().plot_reflection_hyperplanes()
                sage: RootSystem(["B",3]).ambient_space().plot_reflection_hyperplanes()
                sage: RootSystem(["A",3,1]).weight_space().plot_reflection_hyperplanes()
                sage: RootSystem(["B",3,1]).ambient_space().plot_reflection_hyperplanes()
                sage: RootSystem(["A",2,1]).weight_space().plot_reflection_hyperplanes(affine=False, level=1)
                sage: RootSystem(["A",2]).root_lattice().plot_reflection_hyperplanes()

            TESTS::

                sage: L = RootSystem(["A",2]).ambient_space()
                sage: print L.plot_reflection_hyperplanes().description()
                Text '$H_{\alpha^\vee_{1}}$' at the point (-1.81...,3.15)
                Text '$H_{\alpha^\vee_{2}}$' at the point (1.81...,3.15)
                Line defined by 2 points: [(-1.73..., 3.0), (1.73..., -3.0)]
                Line defined by 2 points: [(1.73..., 3.0), (-1.73..., -3.0)]

                sage: print L.plot_reflection_hyperplanes("all").description()
                Text '$H_{\alpha^\vee_{1} + \alpha^\vee_{2}}$' at the point (3.15,0.0)
                Text '$H_{\alpha^\vee_{1}}$' at the point (-1.81...,3.15)
                Text '$H_{\alpha^\vee_{2}}$' at the point (1.81...,3.15)
                Line defined by 2 points: [(-1.73..., 3.0), (1.73..., -3.0)]
                Line defined by 2 points: [(1.73..., 3.0), (-1.73..., -3.0)]
                Line defined by 2 points: [(3.0, 0.0), (-3.0, 0.0)]

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: print L.plot_reflection_hyperplanes().description()
                Text '$H_{\alpha^\vee_{0}}$' at the point (3.15,0.90...)
                Text '$H_{\alpha^\vee_{1}}$' at the point (-1.81...,3.15)
                Text '$H_{\alpha^\vee_{2}}$' at the point (1.81...,3.15)
                Line defined by 2 points: [(-1.73..., 3.0), (1.73..., -3.0)]
                Line defined by 2 points: [(1.73..., 3.0), (-1.73..., -3.0)]
                Line defined by 2 points: [(3.0, 0.86...), (-3.0, 0.86...)]

            .. TODO:: Provide an option for transparency?
            """
            plot_options = self.plot_parse_options(**options)

            coroot_lattice = self.root_system.coroot_lattice()
            # Recall that the coroots are given by the roots of the coroot lattice
            if collection == "simple":
                coroots = coroot_lattice.simple_roots()
            elif collection == "classical":
                if not self.cartan_type().is_affine():
                    raise ValueError("plotting classical reflection hyperplanes only available in affine type")
                raise NotImplementedError("classical roots")
            elif collection == "all":
                assert self.cartan_type().is_finite(), "plotting all reflection hyperplanes only available in finite type"
                coroots = coroot_lattice.positive_roots()
            elif isinstance(collection, (list, tuple)):
                coroots = collection
            else:
                raise ValueError("Unknown value: %s"%collection)

            G = plot_options.empty()
            for coroot in coroots:
                G += plot_options.reflection_hyperplane(coroot)
            return plot_options.finalize(G)


        def plot_hedron(self, **options):
            r"""
            Plot the polyhedron whose vertices are given by the orbit
            of `\rho`.

            In type `A`, this is the usual permutohedron.

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES::

                sage: RootSystem(["A",2]).ambient_space().plot_hedron()
                sage: RootSystem(["A",3]).ambient_space().plot_hedron()
                sage: RootSystem(["B",3]).ambient_space().plot_hedron()
                sage: RootSystem(["C",3]).ambient_space().plot_hedron()
                sage: RootSystem(["D",3]).ambient_space().plot_hedron()

            Surprise: polyhedrons of large dimension know how to
            project themselves nicely::

                sage: RootSystem(["F",4]).ambient_space().plot_hedron() # long time

            TESTS::

                sage: L = RootSystem(["B",2]).ambient_space()
                sage: print L.plot_hedron().description()
                Polygon defined by 8 points: [(1.5, 0.5), (0.5, 1.5), (-0.5, 1.5), (-1.5, 0.5), (-1.5, -0.5), (-0.5, -1.5), (0.5, -1.5), (1.5, -0.5)]
                Line defined by 2 points: [(-0.5, -1.5), (0.5, -1.5)]
                Line defined by 2 points: [(-0.5, 1.5), (0.5, 1.5)]
                Line defined by 2 points: [(-1.5, -0.5), (-0.5, -1.5)]
                Line defined by 2 points: [(-1.5, -0.5), (-1.5, 0.5)]
                Line defined by 2 points: [(-1.5, 0.5), (-0.5, 1.5)]
                Line defined by 2 points: [(0.5, -1.5), (1.5, -0.5)]
                Line defined by 2 points: [(0.5, 1.5), (1.5, 0.5)]
                Line defined by 2 points: [(1.5, -0.5), (1.5, 0.5)]
                Point set defined by 8 point(s): [(-1.5, -0.5), (-1.5, 0.5), (-0.5, -1.5), (-0.5, 1.5), (0.5, -1.5), (0.5, 1.5), (1.5, -0.5), (1.5, 0.5)]
            """
            from sage.geometry.polyhedron.all import Polyhedron
            plot_options = self.plot_parse_options(**options)
            assert self.cartan_type().is_finite()
            vertices = [plot_options.projection(vertex)
                        for vertex in self.rho().orbit()]
            return Polyhedron(vertices=vertices).plot()

        def plot_fundamental_chamber(self, style="normal", **options):
            r"""
            Plot the (classical) fundamental chamber.

            INPUT:

            - ``style`` -- ``"normal"`` or ``"classical"`` (default: ``"normal"``)

            - ``**options`` -- Plotting options

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES:

            2D plots::

                sage: RootSystem(["B",2]).ambient_space().plot_fundamental_chamber()
                sage: RootSystem(["B",2,1]).ambient_space().plot_fundamental_chamber()
                sage: RootSystem(["B",2,1]).ambient_space().plot_fundamental_chamber("classical")

            3D plots::

                sage: RootSystem(["A",3,1]).weight_space() .plot_fundamental_chamber()
                sage: RootSystem(["B",3,1]).ambient_space().plot_fundamental_chamber()

            This feature is currently not available in the root lattice/space::

                sage: list(RootSystem(["A",2]).root_lattice().plot_fundamental_chamber())
                Traceback (most recent call last):
                ...
                TypeError: classical fundamental chamber not yet available in the root lattice

            TESTS::

                sage: L = RootSystem(["B",2,1]).ambient_space()
                sage: print L.plot_fundamental_chamber().description()
                Polygon defined by 3 points:     [(0.5, 0.5), (1.0, 0.0), (0.0, 0.0)]

                sage: print L.plot_fundamental_chamber(style="classical").description()
                Polygon defined by 3 points:     [(0.0, 0.0), (3.0, 3.0), (3.0, 0.0)]
            """
            plot_options = self.plot_parse_options(**options)
            if not hasattr(self, "fundamental_weights"):
                raise TypeError("classical fundamental chamber not yet available in the root lattice")
            Lambda = self.fundamental_weights()
            cartan_type = self.cartan_type()
            if style=="classical":
                if not cartan_type.is_affine():
                    raise TypeError("classical fundamental chamber only available in affine type")
                I = cartan_type.classical().index_set()
                lines = [Lambda[cartan_type.special_node()]]
            else:
                I = cartan_type.index_set()
                lines = []
            return plot_options.cone(rays = [Lambda[i] for i in I],
                                     lines=lines,
                                     color="lightgrey",
                                     alpha=.3)

        def plot_alcoves(self, alcoves=True, alcove_labels=False, wireframe=False, **options):
            r"""
            Plot the alcoves and optionaly their labels.

            INPUT:

            - ``alcoves`` -- a list of alcoves or ``True`` (default: ``True``)

            - ``alcove_labels`` -- a boolean or a number specifying at
              which level to put the label (default: ``False``)

            - ``**options`` -- Plotting options

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a
                  tutorial on root system plotting, and in particular
                  how the alcoves can be specified.

            EXAMPLES:

            2D plots::

                sage: RootSystem(["B",2,1]).ambient_space().plot_alcoves()                      # long time (3s)

            3D plots::

                sage: RootSystem(["A",2,1]).weight_space() .plot_alcoves(affine=False)          # long time (3s)
                sage: RootSystem(["G",2,1]).ambient_space().plot_alcoves(affine=False, level=1) # long time (3s)

            Here we plot a single alcove::

                sage: L = RootSystem(["A",3,1]).ambient_space()
                sage: W = L.weyl_group()
                sage: L.plot(alcoves=[W.one()], reflection_hyperplanes=False, bounding_box=2)

            TESTS::

                sage: L = RootSystem(["A",2,1]).weight_space()
                sage: p = L.plot_alcoves(alcoves=[[0,0]])
                sage: print p.description()
                Line defined by 2 points: [(-1.0, 0.0), (0.0, -1.0)]
                Line defined by 2 points: [(-1.0, 1.0), (-1.0, 0.0)]
                Line defined by 2 points: [(-1.0, 1.0), (0.0, 0.0)]
                Line defined by 2 points: [(0.0, 0.0), (-1.0, 0.0)]
                Line defined by 2 points: [(0.0, 0.0), (0.0, -1.0)]
                Line defined by 2 points: [(0.0, 0.0), (1.0, -1.0)]
                Line defined by 2 points: [(0.0, 1.0), (-1.0, 1.0)]
                Line defined by 2 points: [(0.0, 1.0), (0.0, 0.0)]
                Line defined by 2 points: [(0.0, 1.0), (1.0, 0.0)]
                Line defined by 2 points: [(1.0, -1.0), (0.0, -1.0)]
                Line defined by 2 points: [(1.0, 0.0), (0.0, 0.0)]
                Line defined by 2 points: [(1.0, 0.0), (1.0, -1.0)]
                sage: [(line.options()['rgbcolor'], line.options()['thickness']) for line in p]
                [('black', 2), ('blue', 1), ('red', 1), ('black', 2), ('blue', 1),
                 ('black', 2), ('red', 1), ('black', 2), ('blue', 1),
                 ('black', 2), ('red', 1), ('black', 2)]
            """
            plot_options = self.plot_parse_options(**options)
            if not hasattr(self, "fundamental_weights"):
                raise TypeError("alcoves not yet available in the root lattice")
            Lambda = self.fundamental_weights()
            cartan_type = self.cartan_type()
            I = cartan_type.index_set()
            W = self.weyl_group()
            if alcove_labels is not False:
                rho = self.rho()
                if alcove_labels is not True:
                    # The input is the desired level
                    rho = rho * alcove_labels / rho.level()
                else:
                    rho = plot_options.intersection_at_level_1(rho)
            # The rays of the fundamental alcove
            fundamental_alcove_rays = Lambda.map(plot_options.intersection_at_level_1)

            def alcove_in_bounding_box(w):
                return any(plot_options.in_bounding_box(w.action(fundamental_alcove_rays[i]))
                           for i in I)
            def alcove_facet(w, i):
                # Alcove facets with degenerate intersection with the
                # bounding box bring no information; we might as well
                # not draw them. Besides this avoids ugly fat points
                # in dimension 2.
                return plot_options.cone(rays=[w.action(fundamental_alcove_rays[j]) for j in I if j != i],
                                         color=plot_options.color(i),
                                         thickness=plot_options.thickness(i),
                                         wireframe=wireframe,
                                         draw_degenerate=False)
            def alcove_label(w):
                label = "$1$" if w.is_one() else "$s_{"+"".join(str(j) for j in w.reduced_word())+"}$"
                position = plot_options.projection(w.action(rho))
                if position in plot_options.bounding_box:
                    return plot_options.text(label, position)
                else:
                    return plot_options.empty()

            G = plot_options.empty()
            if alcoves is not True:
                alcoves = list(alcoves)
            if alcoves is True or (len(alcoves)>0 and W.is_parent_of(alcoves[0])):
                if alcoves is True:
                    alcoves = W.weak_order_ideal(alcove_in_bounding_box, side="right")
                # We assume that the fundamental alcove lies within
                # the bounding box, and explore the alcoves
                # intersecting the bounding box by going up right
                # order (i.e. going away from the fundamental alcove)
                for w in alcoves:
                    for i in w.descents(side="right", positive=True):
                        G += alcove_facet(w, i)
                    if alcove_labels is not False:
                        G += alcove_label(w)
            else:
                if not cartan_type.is_affine():
                    raise TypeError("alcoves=list only available in affine type")
                translation_factors = cartan_type.translation_factors()
                simple_roots = self.simple_roots()
                translation_vectors = Family({i: translation_factors[i]*simple_roots[i]
                                          for i in cartan_type.classical().index_set()})
                # The elements of the classical weyl group, as elements of W
                W0 = [W.from_reduced_word(w.reduced_word()) for w in self.weyl_group().classical()]
                for alcove in alcoves:
                    # The translation mapping the center of the
                    # fundamental polygon to polygon indexed by alcove
                    shift = sum(x*v for x,v in zip(alcove, translation_vectors))
                    shift = W.from_morphism(shift.translation)
                    for w in W0:
                        for i in w.descents(side="right", positive=True):
                            G += alcove_facet(shift * w, i)
                        if alcove_labels:
                            G += alcove_label(w)
            return plot_options.finalize(G)

            # In this alternative commented-out implementation, the
            # alcove picture is constructed directly in the
            # projection. It only works for rank 2+1 with, but it is
            # faster; we keep for reference for now. With #12553
            # (Cythoned PPL polytopes), the difference is likely to
            # disappear. If this is confirmed, the code below should be discarded.
            #
            # from sage.plot.line import line
            # translation_vectors = Family({i: translation_factors[i]*plot_options.projection(simple_roots[i])
            #                               for i in cartan_type.classical().index_set()})
            #
            # # For each polygon P to be drawn, alcoves_shift contains the translation
            # # from fundamental polygon to P in the plot coordinate system
            # def immutable_vector(x):
            #     # Takes care of possible numerical instabilities
            #     x = x.numerical_approx(8)
            #     x.set_immutable()
            #     return x
            #
            # # Construct the fundamental polygon
            # # The classical group acting on ``self``
            # W0 = self.weyl_group().classical().list()
            # # The coordinates of the vertices of the fundamental alcove
            # fundamental_alcove_rays = Lambda.map(plot_options.intersection_at_level_1)
            # # The coordinates of the vertices of the fundamental polygon
            # fundamental_polygon_rays = {
            #     (i, w): plot_options.projection(w.action(fundamental_alcove_rays[i]))
            #     for w in W0
            #     for i in I
            #     }
            #
            # # Get the center of the polygons
            # if alcoves is True:
            #     def neighbors(x):
            #         return filter(lambda y: plot_options.bounding_box.contains(plot_options.origin_projected+y),
            #                       [immutable_vector(x+epsilon*t) for t in translation_vectors for epsilon in [-1,1]])
            #     alcoves_shift = list(TransitiveIdeal(neighbors, [immutable_vector(plot_options.origin_projected)]))
            # else:
            #     alcoves_shift = [sum(x*v for x,v in zip(alcove, translation_vectors))
            #                      for alcove in alcoves]
            #
            # G = plot_options.empty()
            # for shift in alcoves_shift:
            #     # for each center of polygon and each element of classical
            #     # parabolic subgroup, we have to draw an alcove.
            #     polygon_center = plot_options.origin_projected + shift
            #
            #     for w in W0:
            #         for i in I:
            #             facet_indices = [j for j in I if j != i]
            #             assert len(facet_indices) == 2
            #             facet = [fundamental_polygon_rays[j, w] + shift for j in facet_indices]
            #             # This takes a bit of time; do we really want that feature?
            #             #if not all(bounding_box_as_polytope.contains(v) for v in facet):
            #             #    continue
            #             G += line(facet,
            #                       rgbcolor = plot_options.color(i),
            #                       thickness = 2 if i == special_node else 1)


        def plot_bounding_box(self, **options):
            r"""
            Plot the bounding box.

            INPUT:

            - ``**options`` -- Plotting options

            This is mostly for testing purposes.

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES::

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: L.plot_bounding_box()

            TESTS::

                sage: list(L.plot_bounding_box())
                [Polygon defined by 4 points]
            """
            plot_options = self.plot_parse_options(**options)
            return plot_options.bounding_box.plot(color="gray", alpha=0.5, wireframe=False)

        def plot_alcove_walk(self, word, start=None, foldings=None, color ="orange", **options):
            r"""
            Plot an alcove walk.

            INPUT:

            - ``word`` -- a list of elements of the index set
            - ``foldings`` -- a list of booleans or ``None`` (default: ``None``)
            - ``start`` -- an element of this space (default: ``None`` for `\rho`)
            - ``**options`` -- plotting options

            .. SEEALSO::

                - :meth:`plot` for a description of the plotting options
                - :ref:`sage.combinat.root_system.plot` for a tutorial
                  on root system plotting

            EXAMPLES:

            An alcove walk of type `A_2^{(1)}`::

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: w1 = [0,2,1,2,0,2,1,0,2,1,2,1,2,0,2,0,1,2,0]
                sage: p = L.plot_alcoves(bounding_box=5)           # long time (5s)
                sage: p += L.plot_alcove_walk(w1)                  # long time
                sage: p                                            # long time

            The same plot with another alcove walk::

                sage: w2 = [2,1,2,0,2,0,2,1,2,0,1,2,1,2,1,0,1,2,0,2,0,1,2,0,2]
                sage: p += L.plot_alcove_walk(w2, color="orange")  # long time

            And another with some foldings::

                sage: L.plot_alcoves(bounding_box=3) + \
                ...   L.plot_alcove_walk([0,1,2,0,2,0,1,2,0,1],
                ...                      foldings = [False, False, True, False, False, False, True, False, True, False],
                ...                      color="green")            # long time (3s)

            TESTS::

                sage: L = RootSystem(["A",2,1]).weight_space()
                sage: p = L.plot_alcove_walk([0,1,2,0,2,0,1,2,0,1],
                ...                          foldings = [False, False, True, False, False, False, True, False, True, False],
                ...                          color="green",
                ...                          start=L.rho())
                sage: print p.description()
                Line defined by 2 points: [(-1.0, 8.0), (-1.5, 9.0)]
                Line defined by 2 points: [(1.0, 4.0), (1.5, 4.5)]
                Line defined by 2 points: [(1.0, 7.0), (1.5, 6.0)]
                Arrow from (-1.0,5.0) to (-2.0,7.0)
                Arrow from (-1.0,8.0) to (1.0,7.0)
                Arrow from (-1.5,9.0) to (-1.0,8.0)
                Arrow from (-2.0,7.0) to (-1.0,8.0)
                Arrow from (1.0,1.0) to (2.0,2.0)
                Arrow from (1.0,4.0) to (-1.0,5.0)
                Arrow from (1.0,7.0) to (2.0,8.0)
                Arrow from (1.5,4.5) to (1.0,4.0)
                Arrow from (1.5,6.0) to (1.0,7.0)
                Arrow from (2.0,2.0) to (1.0,4.0)
            """
            from sage.plot.line import line
            from sage.plot.arrow import arrow
            plot_options = self.plot_parse_options(**options)
            W = self.weyl_group()
            s = W.simple_reflections()
            if start is None:
                start = plot_options.intersection_at_level_1(self.rho())
            if foldings is None:
                foldings = [False] * len(word)
            w = W.one()
            source  = plot_options.projection(start)
            G = plot_options.empty()
            for (i, folding) in zip(word, foldings):
                w = w * s[i]
                target = plot_options.projection(w.action(start))
                if folding:
                    middle = (source+target)/2
                    G += line ([source, middle], rgbcolor=color)
                    G += arrow(middle, source, rgbcolor=color, arrowsize=plot_options._arrowsize)
                    # reset w
                    w = w * s[i]
                else:
                    G += arrow(source, target, rgbcolor=color, arrowsize=plot_options._arrowsize)
                    source=target
            return G

        @cached_method
        def _maximum_root_length(self):
            r"""
            Return the square of the maximum of the root lengths for irreducible finite type root systems.

            EXAMPLES::

                sage: Q = RootSystem(['C',2]).root_lattice()
                sage: Q._maximum_root_length()
                4
                sage: Q = RootSystem(['G',2]).root_lattice()
                sage: Q._maximum_root_length()
                6
                sage: Q = RootSystem(['A',3]).root_lattice()
                sage: Q._maximum_root_length()
                2
            """
            ct = self.cartan_type()
            if not ct.is_irreducible():
                raise NotImplementedError, "Implemented only for irreducible finite root systems"
            if not ct.is_finite():
                raise NotImplementedError, "Implemented only for irreducible finite root systems"
            L = self.root_system.ambient_space() # uses peculiarities of ambient embedding
            return max([root.scalar(root) for root in L.simple_roots()])

    ##########################################################################

    class ElementMethods:

        @abstract_method
        def scalar(self, lambdacheck):
            """
            Implement the natural pairing with the coroot lattice.

            INPUT:

            - ``self`` -- an element of a root lattice realization
            - ``lambdacheck`` -- an element of the coroot lattice or coroot space

            OUTPUT: the scalar product of ``self`` and ``lambdacheck``

            EXAMPLES::

                sage: L = RootSystem(['A',4]).root_lattice()
                sage: alpha      = L.simple_roots()
                sage: alphacheck = L.simple_coroots()
                sage: alpha[1].scalar(alphacheck[1])
                2
                sage: alpha[1].scalar(alphacheck[2])
                -1
                sage: matrix([ [ alpha[i].scalar(alphacheck[j])
                ...              for i in L.index_set() ]
                ...            for j in L.index_set() ])
                [ 2 -1  0  0]
                [-1  2 -1  0]
                [ 0 -1  2 -1]
                [ 0  0 -1  2]

            TESTS::

                sage: super(sage.combinat.root_system.root_space.RootSpaceElement,alpha[1]).scalar(alphacheck[1])
                Traceback (most recent call last):
                ...
                NotImplementedError: <abstract method scalar at ...>
            """

        ##########################################################################
        # Action and orbits w.r.t. the Weyl group
        ##########################################################################

        def simple_reflection(self, i):
            """
            Returns the image of ``self`` by the `i`-th simple reflection.

            EXAMPLES::

                sage: alpha = RootSystem(["A", 3]).root_lattice().alpha()
                sage: alpha[1].simple_reflection(2)
                alpha[1] + alpha[2]

                sage: Q = RootSystem(['A', 3, 1]).weight_lattice(extended = True)
                sage: Lambda = Q.fundamental_weights()
                sage: L = Lambda[0] + Q.null_root()
                sage: L.simple_reflection(0)
                -Lambda[0] + Lambda[1] + Lambda[3]
            """
            # Subclasses should optimize whenever possible!
            return self.parent().simple_reflection(i)(self)

        def simple_reflections(self):
            """
            The images of self by all the simple reflections

            EXAMPLES::

                sage: alpha = RootSystem(["A", 3]).root_lattice().alpha()
                sage: alpha[1].simple_reflections()
                [-alpha[1], alpha[1] + alpha[2], alpha[1]]
            """
            return [s(self) for s in self.parent().simple_reflections()]

        def orbit(self):
            r"""
            The orbit of self under the action of the Weyl group

            EXAMPLES:

            `\rho` is a regular element whose orbit is in bijection with the Weyl group.
            In particular, it as 6 elements for the symmetric group `S_3`::

                sage: L = RootSystem(["A", 2]).ambient_lattice()
                sage: sorted(L.rho().orbit())               # the output order is not specified
                [(1, 2, 0), (1, 0, 2), (2, 1, 0), (2, 0, 1), (0, 1, 2), (0, 2, 1)]

                sage: L = RootSystem(["A", 3]).weight_lattice()
                sage: len(L.rho().orbit())
                24
                sage: len(L.fundamental_weights()[1].orbit())
                4
                sage: len(L.fundamental_weights()[2].orbit())
                6
            """
            return [x for x in TransitiveIdeal(attrcall('simple_reflections'), [self])]

        ##########################################################################
        #
        ##########################################################################

        @abstract_method(optional=True)
        def associated_coroot(self):
            """
            Returns the coroot associated to this root

            EXAMPLES::

                sage: alpha = RootSystem(["A", 3]).root_space().simple_roots()
                sage: alpha[1].associated_coroot()
                alphacheck[1]
            """

        def reflection(self, root, use_coroot = False):
            r"""
            Reflects ``self`` across the hyperplane orthogonal to ``root``.

            If ``use_coroot`` is True, ``root`` is interpreted as a coroot.

            EXAMPLES::

                sage: R = RootSystem(['C',4])
                sage: weight_lattice = R.weight_lattice()
                sage: mu = weight_lattice.from_vector(vector([0,0,1,2]))
                sage: coroot_lattice = R.coroot_lattice()
                sage: alphavee = coroot_lattice.from_vector(vector([0,0,1,1]))
                sage: mu.reflection(alphavee, use_coroot=True)
                6*Lambda[2] - 5*Lambda[3] + 2*Lambda[4]
                sage: root_lattice = R.root_lattice()
                sage: beta = root_lattice.from_vector(vector([0,1,1,0]))
                sage: mu.reflection(beta)
                Lambda[1] - Lambda[2] + 3*Lambda[4]
            """
            if use_coroot:
                return self - self.scalar(root) * root.associated_coroot()
            else:
                return self - self.scalar(root.associated_coroot()) * root


        ##########################################################################
        # Descents
        ##########################################################################

        def has_descent(self, i, positive=False):
            """
            Test if self has a descent at position `i`, that is if self is
            on the strict negative side of the `i^{th}` simple reflection
            hyperplane.

            If positive if True, tests if it is on the strict positive
            side instead.

            EXAMPLES::

                sage: space=RootSystem(['A',5]).weight_space()
                sage: alpha=RootSystem(['A',5]).weight_space().simple_roots()
                sage: [alpha[i].has_descent(1) for i in space.index_set()]
                [False, True, False, False, False]
                sage: [(-alpha[i]).has_descent(1) for i in space.index_set()]
                [True, False, False, False, False]
                sage: [alpha[i].has_descent(1, True) for i in space.index_set()]
                [True, False, False, False, False]
                sage: [(-alpha[i]).has_descent(1, True) for i in space.index_set()]
                [False, True, False, False, False]
                sage: (alpha[1]+alpha[2]+alpha[4]).has_descent(3)
                True
                sage: (alpha[1]+alpha[2]+alpha[4]).has_descent(1)
                False
                sage: (alpha[1]+alpha[2]+alpha[4]).has_descent(1, True)
                True
            """
            s = self.scalar(self.parent().simple_coroots()[i])
            if positive:
                return s > 0
            else:
                return s < 0

        def first_descent(self, index_set=None, positive=False):
            """
            Returns the first descent of pt

            One can use the index_set option to restrict to the parabolic
            subgroup indexed by index_set.

            EXAMPLES::

                sage: space=RootSystem(['A',5]).weight_space()
                sage: alpha=space.simple_roots()
                sage: (alpha[1]+alpha[2]+alpha[4]).first_descent()
                3
                sage: (alpha[1]+alpha[2]+alpha[4]).first_descent([1,2,5])
                5
                sage: (alpha[1]+alpha[2]+alpha[4]).first_descent([1,2,5,3,4])
                5
            """
            if index_set == None:
                index_set = self.parent().index_set()
            for i in index_set:
                if self.has_descent(i, positive):
                    return i
            return None

        def descents(self, index_set=None, positive=False):
            """
            Returns the descents of pt

            EXAMPLES::

                sage: space=RootSystem(['A',5]).weight_space()
                sage: alpha=space.simple_roots()
                sage: (alpha[1]+alpha[2]+alpha[4]).descents()
                [3, 5]
            """
            if index_set==None:
                index_set=self.parent().index_set()
            return [ i for i in index_set if self.has_descent(i, positive) ]

        def to_dominant_chamber(self, index_set = None, positive = True, reduced_word = False):
            r"""
            Returns the unique dominant element in the Weyl group orbit of the vector ``self``.

            If ``positive`` is False, returns the antidominant orbit element.

            With the ``index_set`` optional parameter, this is done with
            respect to the corresponding parabolic subgroup.

            If ``reduced_word`` is True, returns the 2-tuple (``weight``, ``direction``)
            where ``weight`` is the (anti)dominant orbit element and ``direction`` is a reduced word
            for the Weyl group element sending ``weight`` to ``self``.

            .. warning::

                In infinite type, an orbit may not contain a dominant element.
                In this case the function may go into an infinite loop.

                For affine root systems, errors are generated if
                the orbit does not contain the requested kind of representative.
                If the input vector is of positive (resp. negative)
                level, then there is a dominant (resp. antidominant) element in its orbit
                but not an antidominant (resp. dominant) one. If the vector is of level zero,
                then there are neither dominant nor antidominant orbit representatives, except
                for multiples of the null root, which are themselves both dominant and antidominant
                orbit representatives.

            EXAMPLES::

                sage: space=RootSystem(['A',5]).weight_space()
                sage: alpha=RootSystem(['A',5]).weight_space().simple_roots()
                sage: alpha[1].to_dominant_chamber()
                Lambda[1] + Lambda[5]
                sage: alpha[1].to_dominant_chamber([1,2])
                Lambda[1] + Lambda[2] - Lambda[3]
                sage: wl=RootSystem(['A',2,1]).weight_lattice(extended=True)
                sage: mu=wl.from_vector(vector([1,-3,0]))
                sage: mu.to_dominant_chamber(positive=False, reduced_word = True)
                (-Lambda[1] - Lambda[2] - delta, [0, 2])

                sage: R = RootSystem(['A',1,1])
                sage: rl = R.root_lattice()
                sage: nu = rl.zero()
                sage: nu.to_dominant_chamber()
                0
                sage: nu.to_dominant_chamber(positive=False)
                0
                sage: mu = rl.from_vector(vector([0,1]))
                sage: mu.to_dominant_chamber()
                Traceback (most recent call last):
                ...
                ValueError: alpha[1] is not in the orbit of the fundamental chamber
                sage: mu.to_dominant_chamber(positive=False)
                Traceback (most recent call last):
                ...
                ValueError: alpha[1] is not in the orbit of the negative of the fundamental chamber
            """

            if index_set is None:
                # default index set is the entire Dynkin node set
                index_set = self.parent().index_set()
            cartan_type = self.parent().cartan_type()
            # generate errors for infinite loop cases in affine type
            if cartan_type.is_affine():
                if index_set == self.parent().index_set():
                    # If the full affine Weyl group is being used
                    level = self.level()
                    if level > 0:
                        if not positive:
                            raise ValueError, "%s is not in the orbit of the fundamental chamber"%(self)
                    elif level < 0:
                        if positive:
                            raise ValueError, "%s is not in the orbit of the negative of the fundamental chamber"%(self)
                    elif not (self == self.parent().zero()):
                        # nonzero level zero weight
                        if positive:
                            raise ValueError, "%s is not in the orbit of the fundamental chamber"%(self)
                        else:
                            raise ValueError, "%s is not in the orbit of the negative of the fundamental chamber"%(self)
            if reduced_word:
                direction = []
            while True:
                # The first index where it is *not* yet on the positive side
                i = self.first_descent(index_set, positive=(not positive))
                if i is None:
                    if reduced_word:
                        return self, direction
                    else:
                        return self
                else:
                    if reduced_word:
                        direction.append(i)
                    self = self.simple_reflection(i)

        to_positive_chamber = deprecated_function_alias(12667, to_dominant_chamber)

        def reduced_word(self, index_set = None, positive = True):
            r"""
            Returns a reduced word for the inverse of the shortest Weyl group element that sends the vector ``self`` into the dominant chamber.

            With the ``index_set`` optional parameter, this is done with
            respect to the corresponding parabolic subgroup.

            If ``positive`` is False, use the antidominant chamber instead.

            EXAMPLES::

                sage: space=RootSystem(['A',5]).weight_space()
                sage: alpha=RootSystem(['A',5]).weight_space().simple_roots()
                sage: alpha[1].reduced_word()
                [2, 3, 4, 5]
                sage: alpha[1].reduced_word([1,2])
                [2]

            """
            return self.to_dominant_chamber(index_set=index_set,positive=positive,reduced_word = True)[1]


        def is_dominant(self, index_set = None, positive = True):
            r"""
            Returns whether self is dominant.

            This is done with respect to the subrootsystem indicated by the subset of Dynkin nodes
            index_set. If index_set is None then the entire Dynkin node set is used.
            If positive is False then the dominance condition is replaced by antidominance.

            EXAMPLES::

                sage: L = RootSystem(['A',2]).ambient_lattice()
                sage: Lambda = L.fundamental_weights()
                sage: [x.is_dominant() for x in Lambda]
                [True, True]
                sage: [x.is_dominant(positive=False) for x in Lambda]
                [False, False]
                sage: (Lambda[1]-Lambda[2]).is_dominant()
                False
                sage: (-Lambda[1]+Lambda[2]).is_dominant()
                False
                sage: (Lambda[1]-Lambda[2]).is_dominant([1])
                True
                sage: (Lambda[1]-Lambda[2]).is_dominant([2])
                False
                sage: [x.is_dominant() for x in L.roots()]
                [False, True, False, False, False, False]
                sage: [x.is_dominant(positive=False) for x in L.roots()]
                [False, False, False, False, True, False]
            """
            return self.first_descent(index_set, not positive) is None

        def is_dominant_weight(self): # Or is_dominant_integral_weight?
            """
            Tests whether ``self`` is a dominant element of the weight lattice

            EXAMPLES::

                sage: L = RootSystem(['A',2]).ambient_lattice()
                sage: Lambda = L.fundamental_weights()
                sage: [x.is_dominant() for x in Lambda]
                [True, True]
                sage: (3*Lambda[1]+Lambda[2]).is_dominant()
                True
                sage: (Lambda[1]-Lambda[2]).is_dominant()
                False
                sage: (-Lambda[1]+Lambda[2]).is_dominant()
                False

            .. warning::

                The current implementation tests that the scalar products
                with the coroots are all non negative integers, which is not
                sufficient. For example, if `x` is the sum of a dominant
                element of the weight lattice plus some other element
                orthogonal to all coroots, then the current implementation
                erroneously reports `x` to be a dominant weight::

                    sage: x = Lambda[1] + L([-1,-1,-1])
                    sage: x.is_dominant_weight()
                    True
            """
            alphacheck = self.parent().simple_coroots()
            from sage.rings.semirings.non_negative_integer_semiring import NN
            return all(self.inner_product(alphacheck[i]) in NN
                       for i in self.parent().index_set())


        ##########################################################################
        # weak order
        ##########################################################################

        def succ(self):
            r"""
            Returns the immediate successors of self for the weak order

            EXAMPLES::

                sage: L = RootSystem(['A',3]).weight_lattice()
                sage: Lambda = L.fundamental_weights()
                sage: Lambda[1].succ()
                [-Lambda[1] + Lambda[2]]
                sage: L.rho().succ()
                [-Lambda[1] + 2*Lambda[2] + Lambda[3], 2*Lambda[1] - Lambda[2] + 2*Lambda[3], Lambda[1] + 2*Lambda[2] - Lambda[3]]
                sage: (-L.rho()).succ()
                []
           """
            return [ self.simple_reflection(i) for i in self.descents(positive=True) ]

        def pred(self):
            r"""
            Returns the immediate predecessors of self for the weak order

            EXAMPLES::

                sage: L = RootSystem(['A',3]).weight_lattice()
                sage: Lambda = L.fundamental_weights()
                sage: Lambda[1].pred()
                []
                sage: L.rho().pred()
                []
                sage: (-L.rho()).pred()
                [Lambda[1] - 2*Lambda[2] - Lambda[3], -2*Lambda[1] + Lambda[2] - 2*Lambda[3], -Lambda[1] - 2*Lambda[2] + Lambda[3]]
            """
            return [ self.simple_reflection(i) for i in self.descents() ]

        def greater(self):
            r"""
            Returns the elements in the orbit of self which are
            greater than self in the weak order.

            EXAMPLES::

                sage: L = RootSystem(['A',3]).ambient_lattice()
                sage: e = L.basis()
                sage: e[2].greater()
                [(0, 0, 1, 0), (0, 0, 0, 1)]
                sage: len(L.rho().greater())
                24
                sage: len((-L.rho()).greater())
                1
                sage: sorted([len(x.greater()) for x in L.rho().orbit()])
                [1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 8, 8, 8, 8, 12, 12, 12, 24]
            """
            return [x for x in TransitiveIdeal(attrcall('succ'), [self])]

        def smaller(self):
            r"""
            Returns the elements in the orbit of self which are
            smaller than self in the weak order.

            EXAMPLES::

                sage: L = RootSystem(['A',3]).ambient_lattice()
                sage: e = L.basis()
                sage: e[2].smaller()
                [(0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, 0)]
                sage: len(L.rho().smaller())
                1
                sage: len((-L.rho()).smaller())
                24
                sage: sorted([len(x.smaller()) for x in L.rho().orbit()])
                [1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 8, 8, 8, 8, 12, 12, 12, 24]
            """
            return [x for x in TransitiveIdeal(attrcall('pred'), [self])]

        ##########################################################################
        # Level
        ##########################################################################

        def level(self):
            """
            EXAMPLES::

                sage: L = RootSystem(['A',2,1]).weight_lattice()
                sage: L.rho().level()
                3
            """
            if not self.parent().cartan_type().is_affine():
                raise ValueError, "%s does not belong to a lattice of affine Cartan type"%self
            return self.scalar(self.parent().null_coroot())

        @cached_in_parent_method
        def to_simple_root(self, reduced_word=False):
            r"""
            Return (the index of) a simple root in the orbit of the positive root ``self``.

            INPUT:

            - ``self`` -- a positive root
            - ``reduced_word`` -- a boolean (default: ``False``)

            OUTPUT:

            - The index `i` of a simple root `\alpha_i`.
              If ``reduced_word`` is True, this returns instead a pair
              ``(i, word)``, where word is a sequence of reflections
              mapping `\alpha_i` up the root poset to ``self``.

            EXAMPLES::

                sage: L = RootSystem(["A",3]).root_lattice()
                sage: for alpha in L.positive_roots():
                ...       print alpha, alpha.to_simple_root()
                alpha[1] 1
                alpha[2] 2
                alpha[3] 3
                alpha[1] + alpha[2] 2
                alpha[2] + alpha[3] 3
                alpha[1] + alpha[2] + alpha[3] 3
                sage: for alpha in L.positive_roots():
                ...        print alpha, alpha.to_simple_root(reduced_word=True)
                alpha[1] (1, ())
                alpha[2] (2, ())
                alpha[3] (3, ())
                alpha[1] + alpha[2] (2, (1,))
                alpha[2] + alpha[3] (3, (2,))
                alpha[1] + alpha[2] + alpha[3] (3, (1, 2))

            ALGORITHM:

            This method walks from ``self`` down to the antidominant
            chamber by applying successively the simple reflection
            given by the first descent. Since ``self`` is a positive
            root, each step goes down the root poset, and one must
            eventually cross a simple root `\alpha_i`.

            .. SEEALSO::

                - :meth:`first_descent`
                - :meth:`to_dominant_chamber`

            .. WARNING::

                The behavior is not specified if the input is not a
                positive root. For a finite root system, this is
                currently caught (albeit with a not perfect message)::

                    sage: alpha = L.simple_roots()
                    sage: (2*alpha[1]).to_simple_root()
                    Traceback (most recent call last):
                    ...
                    ValueError: -2*alpha[1] - 2*alpha[2] - 2*alpha[3] is not a positive root

                For an infinite root systems, this method may run into
                an infinite reccursion if the input is not a positive
                root.
            """
            F = self.parent().simple_roots().inverse_family()
            try:
                j = F[self]
                if reduced_word:
                    return (j, ())
                else:
                    return j
            except KeyError:
                pass
            j = self.first_descent(positive=True)
            if j is None:
                raise ValueError, "%s is not a positive root"%self
            result = self.simple_reflection(j).to_simple_root(reduced_word=reduced_word)
            if reduced_word:
                return (result[0], (j,) + result[1])
            else:
                return result

        @cached_in_parent_method
        def associated_reflection(self):
            r"""
            Given a positive root ``self``, returns a reduced word for the reflection orthogonal to ``self``.

            Since the answer is cached, it is a tuple instead of a list.

            EXAMPLES::

                sage: RootSystem(['C',3]).root_lattice().simple_root(3).weyl_action([1,2]).associated_reflection()
                (1, 2, 3, 2, 1)
                sage: RootSystem(['C',3]).root_lattice().simple_root(2).associated_reflection()
                (2,)

            """
            i, reduced_word = self.to_simple_root(reduced_word=True)
            return reduced_word + (i,) + tuple(reversed(reduced_word))

        def translation(self, x):
            """
            INPUT:
             - ``self`` - an element `t` at level `0`
             - ``x`` - an element of the same space

            Returns `x` translated by `t`, that is `x+level(x) t`

            EXAMPLES::

                sage: L = RootSystem(['A',2,1]).weight_lattice()
                sage: alpha = L.simple_roots()
                sage: Lambda = L.fundamental_weights()
                sage: t = alpha[2]

            Let us look at the translation of an element of level `1`::

                sage: Lambda[1].level()
                1
                sage: t.translation(Lambda[1])
                -Lambda[0] + 2*Lambda[2]
                sage: Lambda[1] + t
                -Lambda[0] + 2*Lambda[2]

            and of an element of level `0`::

                sage: alpha [1].level()
                0
                sage: t.translation(alpha [1])
                -Lambda[0] + 2*Lambda[1] - Lambda[2]
                sage: alpha[1] + 0*t
                -Lambda[0] + 2*Lambda[1] - Lambda[2]

            The arguments are given in this seemingly unnatural order to
            make it easy to construct the translation function::

                sage: f = t.translation
                sage: f(Lambda[1])
                -Lambda[0] + 2*Lambda[2]
            """
            if not self.level().is_zero():
                raise ValueError, "%s is not of level zero"%(self)
            return x + x.level() * self

        def weyl_action(self, element, inverse = False):
            r"""
            Acts on ``self`` by an element of the Coxeter or Weyl group.

            INPUT:

            - ``element`` -- an element of a Coxeter or Weyl group
              of the same Cartan type, or a tuple or a list (such as a
              reduced word) of elements from the index set.

            - ``inverse`` -- a boolean (default: False); whether to act by the inverse element.

            EXAMPLES::

                sage: wl = RootSystem(['A',3]).weight_lattice()
                sage: mu = wl.from_vector(vector([1,0,-2]))
                sage: mu
                Lambda[1] - 2*Lambda[3]
                sage: mudom, rw = mu.to_dominant_chamber(positive=False, reduced_word = True)
                sage: mudom, rw
                (-Lambda[2] - Lambda[3], [1, 2])

            Acting by a (reduced) word::

                sage: mudom.weyl_action(rw)
                Lambda[1] - 2*Lambda[3]
                sage: mu.weyl_action(rw, inverse = True)
                -Lambda[2] - Lambda[3]

            Acting by an element of the Coxeter or Weyl group on a vector in its own
            lattice of definition (implemented by matrix multiplication on a vector)::

                sage: w = wl.weyl_group().from_reduced_word([1, 2])
                sage: mudom.weyl_action(w)
                Lambda[1] - 2*Lambda[3]

            Acting by an element of an isomorphic Coxeter or Weyl group (implemented by the
            action of a corresponding reduced word)::

                sage: W = WeylGroup(['A',3], prefix="s")
                sage: w = W.from_reduced_word([1, 2])
                sage: wl.weyl_group() == W
                False
                sage: mudom.weyl_action(w)
                Lambda[1] - 2*Lambda[3]

            """

            # TODO, some day: accept an iterator
            if isinstance(element, (tuple, list)):
                # Action by a (reduced) word
                the_word = [x for x in element]
                I = self.parent().index_set()
                if not all(i in I for i in the_word):
                    raise ValueError, "Not all members of %s are in the index set of the %s"%(element, self.parent())
            else:
                if not isinstance(element, Element):
                    raise TypeError, "%s should be an element of a Coxeter group"%(element)
                W = element.parent()
                if W is self.parent().weyl_group():
                    # Action by an element of the Coxeter or Weyl group of ``self``
                    if inverse is True:
                        element = element.inverse()
                    return element.action(self)
                else:
                    # Action by an element of an isomorphic Coxeter or Weyl group
                    if not (W in CoxeterGroups() and W.cartan_type() == self.parent().cartan_type()):
                        raise TypeError, "%s should be an element of a Coxeter group of type %s"%(element, self.parent().cartan_type())
                    the_word = element.reduced_word()
            if inverse is False:
                the_word.reverse()
            for i in the_word:
                self = self.simple_reflection(i)
            return self

        def weyl_stabilizer(self, index_set=None):
            r"""
            Returns the subset of Dynkin nodes whose reflections fix ``self``.

            If ``index_set`` is not None, only consider nodes in this set.
            Note that if ``self`` is dominant or antidominant, then its stabilizer is the
            parabolic subgroup defined by the returned node set.

            EXAMPLES::

                sage: wl = RootSystem(['A',2,1]).weight_lattice(extended = True)
                sage: al = wl.null_root()
                sage: al.weyl_stabilizer()
                [0, 1, 2]
                sage: wl = RootSystem(['A',4]).weight_lattice()
                sage: mu = wl.from_vector(vector([1,1,0,0]))
                sage: mu.weyl_stabilizer()
                [3, 4]
                sage: mu.weyl_stabilizer(index_set = [1,2,3])
                [3]

            """

            if index_set is None:
                index_set = self.parent().cartan_type().index_set()
            alphavee = self.parent().coroot_lattice().basis()
            return [i for i in index_set if self.scalar(alphavee[i]) == 0]

        def is_parabolic_root(self, index_set):
            r"""
            Supposing that ``self`` is a root, is it in the parabolic subsystem with Dynkin nodes ``index_set``?

            INPUT:

            - ``index_set`` -- the Dynkin node set of the parabolic subsystem.

            .. TODO:: This implementation is only valid in the root or weight lattice

            EXAMPLES::

                sage: alpha = RootSystem(['A',3]).root_lattice().from_vector(vector([1,1,0]))
                sage: alpha.is_parabolic_root([1,3])
                False
                sage: alpha.is_parabolic_root([1,2])
                True
                sage: alpha.is_parabolic_root([2])
                False

            """
            for i in self.support():
                if i not in index_set:
                    return False
            return True

        def is_short_root(self):
            r"""
            Is ``self`` a short root?

            Returns False unless the parent is an irreducible root system of finite type
            having two root lengths and ``self`` is of the shorter length.
            There is no check of whether ``self`` is actually a root.

            EXAMPLES::

                sage: Q = RootSystem(['C',2]).root_lattice()
                sage: al = Q.simple_root(1).weyl_action([1,2]); al
                alpha[1] + alpha[2]
                sage: al.is_short_root()
                True
                sage: bt = Q.simple_root(2).weyl_action([2,1,2]); bt
                -2*alpha[1] - alpha[2]
                sage: bt.is_short_root()
                False
                sage: RootSystem(['A',2]).root_lattice().simple_root(1).is_short_root()
                False

            """
            ct = self.parent().cartan_type()
            if not ct.is_irreducible():
                raise ValueError, "Cartan type needs to be irreducible!"
            if not ct.is_finite():
                raise NotImplementedError, "Implemented only for irreducible finite root systems"
            L = self.parent().root_system.ambient_space() # uses peculiarities of ambient embedding
            ls = L(self)
            return ls.scalar(ls) < L._maximum_root_length()
            #Alternative implementation
            #if ct.is_simply_laced():
            #    return False
            #L = self.parent().root_system.ambient_space() # uses peculiarities of ambient embedding
            #ls = L(self)
            #lensq = ls.scalar(ls)
            #if lensq > 2:
            #    return False
            #if lensq == 1:
            #    return True
            ## now only types BCFG remain and the square length is 2
            #if ct.type() == 'C' or ct.type() == 'G':
            #    return True
            #return False

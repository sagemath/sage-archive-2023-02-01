"""
Root lattice realizations
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

from sage.misc.abstract_method import abstract_method, AbstractMethod
from sage.misc.all import attrcall
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_types import Category_over_base_ring
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.sets.family import Family
from sage.rings.all import ZZ, QQ
from sage.combinat.backtrack import TransitiveIdeal
from sage.misc.misc import deprecated_function_alias
from copy import copy

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
    :meth:`~RootLatticeRealizations.ParentMethods.simple_root` (i)
    returning for each `i` in the index set the image of the simple root
    `\alpha_i` under the embedding.

    A root lattice realization must further implement a method on elements
    :meth:`~RootLatticeRealizations.ElementMethods.scalar`, computing
    the scalar product with elements of the coroot lattice or coroot space.

    Using those, this category provides tools for reflections, roots,
    the Weyl group and its action, ...

    .. seealso::

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
            """
            Registers the embedding of the root lattice into ``self``

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

            .. note::

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
                [1, 2, 3, 4]
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
            tester.assertEqual(alpha     .keys(), self.index_set())
            tester.assertEqual(alphacheck.keys(), self.index_set())

            # Check the consistency between simple_root and simple_roots
            for i in self.index_set():
                tester.assertEqual(self.simple_root(i), alpha[i])

            # Check the embeddings from the root lattice and the root space over the same base ring
            root_lattice = self.root_system.root_lattice()
            root_space   = self.root_system.root_space  (self.base_ring())
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
                    tester.assertEqual(alpha[j].scalar(alphacheck[i]), dynkin_diagram[i,j])

            # Check associated_coroot, if it is implemented
            if not isinstance(self.element_class.associated_coroot, AbstractMethod):
                for i in self.index_set():
                    tester.assertEqual(alpha[i].associated_coroot(), alphacheck[i])

            if self.cartan_type().is_affine():
                # Check that the null root is orthogonal to all coroots
                # and similarly for the null coroot
                nullroot = self.null_root()
                nullcoroot = self.null_coroot()
                for k in alpha.keys():
                    assert (nullroot.scalar(alphacheck[k])).is_zero()
                    assert (alpha[k].scalar(nullcoroot)).is_zero()

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
            assert(self.root_system.is_finite())
            assert(self.root_system.is_irreducible())
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
            longest=self.dynkin_diagram().outgoing_edges()[0]
            for j in self.dynkin_diagram().outgoing_edges():
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

                sage: for T in CartanType.samples(finite = True, crystalographic = True):
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
            assert self.cartan_type().is_finite()
            return TransitiveIdeal(attrcall('pred'), self.simple_roots())

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
            assert self.cartan_type().is_finite()
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
            assert self.cartan_type().is_finite()
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
            assert(self.cartan_type().is_affine())
            coef = self.cartan_type().acheck()
            return sum(coef[k]*self.simple_coroots()[k] for k in coef.keys())

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
            assert to_negative == True # Not implemented otherwise!!!!!
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
            # from sage.misc.misc import deprecation
            # deprecation("The lazy attribute pi is deprecated; please use the simple_projections method.", 'Sage Version 5.0')
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

            .. seealso:: :meth:`almost_positive_roots`, :meth:`tau_plus_minus`

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

                .. [CFZ2] Chapoton, Fomin, Zelevinsky - Polytopal realizations of generalized associahedra
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

    class ElementMethods:

        @abstract_method
        def scalar(self, lambdacheck):
            """
            The natural pairing with the coroot lattice

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

        def to_dominant_chamber(self, index_set = None, positive = True, get_direction = False):
            r"""
            Returns the unique dominant element in the Weyl group orbit of the vector ``self``.

            If ``positive`` is False, returns the antidominant orbit element.

            With the ``index_set`` optional parameter, this is done with
            respect to the corresponding parabolic subgroup.

            If ``get_direction`` is True, returns the 2-tuple (``weight``, ``direction``)
            where ``weight`` is the (anti)dominant orbit element and ``direction`` is a reduced word
            for the Weyl group element sending ``weight`` to ``self``.

            .. warning::

                In infinite type, an orbit may not contain a dominant element.
                In this case the function may go into an infinite loop.

                For affine root systems, assertion errors are generated if
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
                sage: mu.to_dominant_chamber(positive=False, get_direction = True)
                (-Lambda[1] - Lambda[2] - delta, [0, 2])

                sage: R = RootSystem(['A',1,1])
                sage: rl = R.root_lattice()
                sage: mu = rl.from_vector(vector([0,1]))
                sage: mu.to_dominant_chamber()
                Traceback (most recent call last):
                ...
                AssertionError: This element is not in the orbit of the fundamental chamber

            """

            if index_set is None:
                # default index set is the entire Dynkin node set
                index_set = self.parent().index_set()
            cartan_type = self.parent().cartan_type()
            # generate assertion errors for infinite loop cases in affine type
            if cartan_type.is_affine():
                if index_set == self.parent().index_set():
                    # If the full affine Weyl group is being used
                    level = self.level()
                    if level > 0:
                        assert positive, "This element is not in the orbit of the fundamental chamber"
                    elif level < 0:
                        assert not positive, "This element is not in the orbit of the negative of the fundamental chamber"
                    else:
                        # level zero
                        if positive:
                            assert self.is_dominant(), "This element is not in the orbit of the fundamental chamber"
                        else:
                            assert self.is_dominant(), "This element is not in the orbit of the negative of the fundamental chamber"
            if get_direction:
                direction = []
            while True:
                # The first index where it is *not* yet on the positive side
                i = self.first_descent(index_set, positive=(not positive))
                if i is None:
                    if get_direction:
                        return self, direction
                    else:
                        return self
                else:
                    if get_direction:
                        direction.append(i)
                    self = self.simple_reflection(i)

        to_positive_chamber = deprecated_function_alias(to_dominant_chamber, "Sage 4.8")

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
            return self.to_dominant_chamber(index_set=index_set,positive=positive,get_direction = True)[1]


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
            assert(self.parent().cartan_type().is_affine())
            return self.scalar(self.parent().null_coroot())


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
            assert self.level().is_zero()
            return x + x.level() * self

        def weyl_action(self, w = None, reduced_word = None, inverse = False):
            r"""
            Acts on ``self`` by a Weyl group element.

            INPUT:
            - If ``w`` is not None, use it to act.
            - If ``reduced_word`` is not None, use it to act.
            - Exactly one of ``w`` and ``reduced_word`` should not be None.
            - If ``inverse`` is True, act by the inverse element.

            EXAMPLES::

                sage: wl = RootSystem(['A',2,1]).weight_lattice(extended = True)
                sage: mu = wl.from_vector(vector([1,-3,0]))
                sage: mu
                Lambda[0] - 3*Lambda[1]
                sage: mudom, rw = mu.to_dominant_chamber(positive=False, get_direction = True)
                sage: mudom, rw
                (-Lambda[1] - Lambda[2] - delta, [0, 2])
                sage: mudom.weyl_action(reduced_word = rw)
                Lambda[0] - 3*Lambda[1]
                sage: mu.weyl_action(reduced_word = rw, inverse = True)
                -Lambda[1] - Lambda[2] - delta

            """

            if w is None:
                assert reduced_word is not None
                rw = copy(reduced_word)
            else:
                rw = w.reduced_word()
            if not inverse:
                rw.reverse()
            for i in rw:
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

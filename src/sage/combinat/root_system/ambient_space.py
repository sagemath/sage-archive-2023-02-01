r"""
Ambient spaces
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from root_lattice_realization import RootLatticeRealizationElement
from weight_lattice_realization import WeightLatticeRealization
from sage.rings.all import QQ
from sage.misc.cachefunc import ClearCacheOnPickle

class AmbientSpace(ClearCacheOnPickle, CombinatorialFreeModule, WeightLatticeRealization):
    r"""
    Abstract class for ambient spaces

    Any implementation of this class should implement a class method
    smallest_base_ring as described below, and a method dimension
    working on a partially initialized instance with just root_system
    as attribute. There is no safe default implementation for the later,
    so none is provided.

    EXAMPLES::

        sage: AL = RootSystem(['A',2]).ambient_lattice()

    Caveat: Most of the ambient spaces currently have a basis indexed
    by `0,\dots, n`, unlike the usual mathematical convention::

        sage: e = AL.basis()
        sage: e[0], e[1], e[2]
        ((1, 0, 0), (0, 1, 0), (0, 0, 1))

    This will be cleaned up!

    TESTS::
        sage: types = CartanType.samples(finite=True, crystalographic = True)
        sage: for e in [ct.root_system().ambient_space() for ct in types]:
        ...       if e is not None:
        ...            TestSuite(e).run()

    """
    def __init__(self, root_system, base_ring):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: s = e.simple_reflections()

        """
        self.root_system = root_system
        CombinatorialFreeModule.__init__(self, base_ring,
                                         range(0,self.dimension()),
                                         element_class = AmbientSpaceElement,
                                         prefix='e')

        # FIXME: here for backward compatibility;
        # Should we use dimension everywhere?
        self.n = self.dimension()


    def _test_norm_of_simple_roots(self, **options):
        """
        Tests that the norm of the roots is, up to an overal constant factor,
        the norm of the roots is given by the symmetrizer of the Cartan matrix.

        Not yet implemented for reducible Cartan types.
        """
        tester = self._tester(**options)
        T = self.cartan_type()
        if T.is_reducible():
            return
        D = T.symmetrizer()
        alpha = self.simple_roots()
        tester.assertEquals(len( set( alpha[i].scalar(alpha[i]) / D[i] for i in self.index_set() ) ), 1)

    # FIXME: attribute or method?
    def dimension(self):
        """
        Returns the dimension of this ambient space.

        EXAMPLES::

            sage: from sage.combinat.root_system.ambient_space import AmbientSpace
            sage: e = RootSystem(['F',4]).ambient_space()
            sage: AmbientSpace.dimension(e)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    @classmethod
    def smallest_base_ring(cls):
        """
        Returns the smallest ground ring over which the ambient space can be realized.

        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()
            sage: e.smallest_base_ring()
            Rational Field
        """
        return QQ

    def _repr_(self):
        """
        EXAMPLES::

            sage: RootSystem(['A',4]).ambient_lattice()
            Ambient lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).ambient_space()
            Ambient space of the Root system of type ['B', 4]

        """
        return self._name_string()

    def _name_string(self, capitalize=True, base_ring=False, type=True):
        """
        EXAMPLES::

            sage: RootSystem(['A',4]).ambient_lattice()._name_string()
            "Ambient lattice of the Root system of type ['A', 4]"

        """
        return self._name_string_helper("ambient", capitalize=capitalize, base_ring=base_ring, type=type)

    def __call__(self, v):
        """
        TESTS::

            sage: R = RootSystem(['A',4]).ambient_lattice()
            sage: R([1,2,3,4,5])
            (1, 2, 3, 4, 5)
            sage: len(R([1,0,0,0,0]).monomial_coefficients())
            1
        """
        # This adds coercion from a list
        if isinstance(v, (list, tuple)):
            K = self.base_ring()
            return self._from_dict(dict((i,K(c)) for i,c in enumerate(v) if c))
        else:
            return CombinatorialFreeModule.__call__(self, v)

    def __getitem__(self,i):
        """
        Note that indexing starts at 1.

        EXAMPLES::

            sage: e = RootSystem(['A',2]).ambient_lattice()
            sage: e[1]
            (1, 0, 0)
        """
        return self.monomial(i-1)

    def coroot_lattice(self):
        """
        EXAMPLES::

            sage: e = RootSystem(["A", 3]).ambient_lattice()
            sage: e.coroot_lattice()
            Ambient lattice of the Root system of type ['A', 3]
        """
        return self

    def simple_coroot(self, i):
        r"""
        Returns the i-th simple coroot, as an element of this space

        EXAMPLES::

            sage: R = RootSystem(["A",3])
            sage: L = R.ambient_lattice()
            sage: L.simple_coroot(1)
            (1, -1, 0, 0)
            sage: L.simple_coroot(2)
            (0, 1, -1, 0)
            sage: L.simple_coroot(3)
            (0, 0, 1, -1)
        """
        return self.simple_root(i).associated_coroot()

    def reflection(self, root, coroot=None):
        """
        EXAMPLES::

            sage: e = RootSystem(["A", 3]).ambient_lattice()
            sage: a = e.simple_root(0); a
            (-1, 0, 0, 0)
            sage: b = e.simple_root(1); b
            (1, -1, 0, 0)
            sage: s_a = e.reflection(a)
            sage: s_a(b)
            (0, -1, 0, 0)

        """
        # TODO: get rid of this as one can use the generic implementation
        # (i.e. scalar and associated coroot are implemented)
        return lambda v: v - root.base_ring()(2*root.inner_product(v)/root.inner_product(root))*root

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: e1 = RootSystem(['A',3]).ambient_lattice()
            sage: e2 = RootSystem(['B',3]).ambient_lattice()
            sage: e1 == e1
            True
            sage: e1 == e2
            False
        """
        if self.__class__ != other.__class__:
            return cmp(self.__class__, other.__class__)
        if self.root_system != other.root_system:
            return cmp(self.root_system, other.root_system)
        return 0

class AmbientSpaceElement(CombinatorialFreeModuleElement, RootLatticeRealizationElement):
    def __hash__(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',2]).ambient_space()
            sage: hash(e.simple_root(0))
            -4601450286177489034          # 64-bit
            549810038                     # 32-bit
        """
        return hash(tuple(sorted([(m,c) for m,c in self._monomial_coefficients.iteritems()])))

    # For backward compatibility
    def _repr_(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['A',2]).ambient_space()
            sage: e.simple_root(0)
            (-1, 0, 0)
        """
        return str(self.to_vector())

    def inner_product(self, lambdacheck):
        """
        The scalar product with elements of the coroot lattice
        embedded in the ambient space.

        EXAMPLES::

            sage: e = RootSystem(['A',2]).ambient_space()
            sage: a = e.simple_root(0); a
            (-1, 0, 0)
            sage: a.inner_product(a)
            2
        """
        self_mc = self._monomial_coefficients
        lambdacheck_mc = lambdacheck._monomial_coefficients

        result = self.parent().base_ring().zero()
        for t,c in lambdacheck_mc.iteritems():
            if t not in self_mc:
                continue
            result += c*self_mc[t]
        return result

    scalar = inner_product
    dot_product = inner_product

    def associated_coroot(self):
        """
        EXAMPLES::

            sage: e = RootSystem(['F',4]).ambient_space()
            sage: a = e.simple_root(0); a
            (1/2, -1/2, -1/2, -1/2)
            sage: a.associated_coroot()
            (1, -1, -1, -1)

        """
        # FIXME: make it work over ZZ!
        return self * self.base_ring()(2/self.inner_product(self))

    def is_positive_root(self):
        """
        EXAMPLES::

            sage: R = RootSystem(['A',3]).ambient_space()
            sage: r=R.simple_root(1)+R.simple_root(2)
            sage: r.is_positive_root()
            True
            sage: r=R.simple_root(1)-R.simple_root(2)
            sage: r.is_positive_root()
            False
        """
        return self.parent().rho().scalar(self) > 0

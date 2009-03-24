from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from root_lattice_realization import RootLatticeRealizationElement
from weight_lattice_realization import WeightLatticeRealization
from sage.rings.all import ZZ, QQ

class AmbientSpace(CombinatorialFreeModule, WeightLatticeRealization):
    r"""
    Abstract class for ambient spaces

    Any implementation of this class should implement a class method
    smallest_base_ring as described below, and a method dimension
    working on a partially initialized instance with just root_system
    as attribute. There is no safe default implementation for the later,
    so none is provided.
    """
    def __init__(self, root_system, base_ring):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e == loads(dumps(e))
            True
        """
        self.root_system = root_system
        CombinatorialFreeModule.__init__(self, base_ring,
                                         range(0,self.dimension()),
                                         element_class = AmbientSpaceElement,
                                         prefix='e')

        # FIXME: here for backward compatibility;
        # Should we use dimension everywhere?
        self.n = self.dimension()



    # FIXME: attribute or method?
    def dimension(self):
        """
        Returns the dimension of this ambient space.

        EXAMPLES:
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

        EXAMPLES:
            sage: e = RootSystem(['F',4]).ambient_space()
            sage: e.smallest_base_ring()
            Rational Field
        """
        return QQ

    def __repr__(self):
        """
        EXAMPLES:
            sage: RootSystem(['A',4]).ambient_lattice()
            Ambient lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).ambient_space()
            Ambient space of the Root system of type ['B', 4]

        """
        return self._name_string()

    def _name_string(self, capitalize=True, base_ring=False, type=True):
        """
        EXAMPLES:
            sage: RootSystem(['A',4]).ambient_lattice()._name_string()
            "Ambient lattice of the Root system of type ['A', 4]"

        """
        return self._name_string_helper("ambient", capitalize=capitalize, base_ring=base_ring, type=type)

    def __call__(self, v):
        """
        TESTS:
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

        EXAMPLES:
            sage: e = RootSystem(['A',2]).ambient_lattice()
            sage: e[1]
            (1, 0, 0)
        """
        return self.term(i-1)

    def coroot_lattice(self):
        """
        EXAMPLES:
            sage: e = RootSystem(["A", 3]).ambient_lattice()
            sage: e.coroot_lattice()
            Ambient lattice of the Root system of type ['A', 3]
        """
        return self

    def simple_coroot(self, i):
        r"""
        Returns the i-th simple coroot, as an element of this space

        EXAMPLES:
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
        EXAMPLES:
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

    def _term(self, i):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['A',2]).ambient_lattice()
            sage: e._term(0)
            (1, 0, 0)
        """
        return self.term(i)

    def __cmp__(self, other):
        """
        EXAMPLES:
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
    def __repr__(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',2]).ambient_space()
            sage: e.simple_root(0)
            (-1, 0, 0)
        """
        return str(self.to_vector())

    def inner_product(self, lambdacheck):
        """
        The scalar product with elements of the coroot lattice
        embedded in the ambient space.

        EXAMPLES:
            sage: e = RootSystem(['A',2]).ambient_space()
            sage: a = e.simple_root(0); a
            (-1, 0, 0)
            sage: a.inner_product(a)
            2
        """
        self_mc = self._monomial_coefficients
        lambdacheck_mc = lambdacheck._monomial_coefficients

        result = self.parent().base_ring().zero_element()
        for t,c in lambdacheck_mc.iteritems():
            if t not in self_mc:
                continue
            result += c*self_mc[t]
        return result

    scalar = inner_product
    dot_product = inner_product

    def associated_coroot(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['F',4]).ambient_space()
            sage: a = e.simple_root(0); a
            (1/2, -1/2, -1/2, -1/2)
            sage: a.associated_coroot()
            (1, -1, -1, -1)

        """
        # FIXME: make it work over ZZ!
        return self * self.base_ring()(2/self.inner_product(self))


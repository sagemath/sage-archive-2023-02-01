from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from root_lattice_realization import RootLatticeRealizationElement
from weight_lattice_realization import WeightLatticeRealization
from sage.modules.free_module import FreeModule
from sage.rings.all import ZZ, QQ
from sage.modules.free_module_element import vector

class AmbientSpace(CombinatorialFreeModule, WeightLatticeRealization):
    r"""
    Abstract class for ambient spaces

    Any implementation of this class should implement a class method
    smallest_base_ring as described below, and a method dimension
    working on a partially initialized instance with just root_system
    as attribute. There is no safe default implementation for the later,
    so none is provided.
    """

    # FIXME: attribute or method?
    def dimension(self):
        """
        Returns the dimension of this ambient space
        """
        raise NotImplementedError

    @classmethod
    def smallest_base_ring(cls):
        """
        Returns the smallest ground ring over which the ambient space can be realized
        """
        return QQ;

    def __init__(self, root_system, base_ring):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e == loads(dumps(e))
            True
        """
        self.root_system = root_system
        basis_name = "alphacheck" if root_system.dualSide else "alpha"
        CombinatorialFreeModule.__init__(self, base_ring,\
                                         #range(1,self.dimension()+1),\
                                         range(0,self.dimension()),\
                                         element_class = AmbientSpaceElement,\
                                         prefix='e')

        # FIXME: here for backward compatibility;
        # Should we use dimension everywhere?
        self.n = self.dimension()

    def __repr__(self):
        """
        TEST:
            sage: RootSystem(['A',4]).ambient_lattice()
            Ambient lattice for the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).ambient_space()
            Ambient space for the Root system of type ['B', 4]

        """
        if self.base_ring() == ZZ:
            space = "lattice"
        else:
            space = "space"
        return "Ambient "+space+" for the %s"%self.root_system

    def __call__(self, v):
        """
        TESTS:
            sage: R = RootSystem(['A',4]).ambient_lattice()
            sage: R([1,2,3,4,5])
            (1, 2, 3, 4, 5)
        """
        # This adds coercion from a list
        if isinstance(v, list) or isinstance(v, tuple):
            K = self.base_ring()
            return self._from_dict(dict([(i,K(v[i])) for i in range(len(v))]))
        else:
            return CombinatorialFreeModule.__call__(self, v)

    # For backward compatibility
    def _term(self, i):
        self.term(i)

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
        return self

    def simple_coroot(self, i):
        r"""
        Returns the i-th simple coroot, as an element of this space

        EXAMPLE:
        sage: R = RootSystem(["A",3])
        sage: L = R.ambient_lattice ()
        sage: L.simple_coroot(1)
        (1, -1, 0, 0)
        sage: L.simple_coroot(2)
        (0, 1, -1, 0)
        sage: L.simple_coroot(3)
        (0, 0, 1, -1)
        """
        return self.simple_root(i).associated_coroot()

    def reflection(self, root, coroot=None):
        # TODO: get rid of this as one can use the generic implementation
        # (i.e. scalar and associated coroot are implemented)
        return lambda v: v-2*root.inner_product(v)/root.inner_product(root)*root

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

    # For backward compatibility
    def __repr__(self):
        return str(self.to_vector())

    def inner_product(self, lambdacheck):
        """
        The scalar product with elements of the coroot lattice
        embedded in the ambient space
        """
        assert(lambdacheck.parent() == self.parent())
        return sum((c*self[t] for (t,c) in lambdacheck),
                   self.parent().base_ring().zero_element())

    scalar = inner_product
    dot_product = inner_product

    def associated_coroot(self):
        # FIXME: make it work over ZZ!
        return self * (2/self.inner_product(self))


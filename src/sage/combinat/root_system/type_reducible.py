from sage.combinat.root_system.cartan_type import CartanType_abstract, CartanType_simple
from sage.misc.flatten import flatten
from sage.matrix.constructor import block_diagonal_matrix
from ambient_space import AmbientSpace

class CartanType(CartanType_abstract):
    r"""
    A class for reducible Cartan types
    """

    def __init__(self, types):
        """
        INPUT:
           types -- a list of simple Cartan types
        Reducible root systems are ones that can be factored as
        direct products. Strictly speaking type D2 (corresponding
        to orthogonal groups of degree 4) are reducible since they
        are isomorphic to A1xA1. However type D2 is considered
        irreducible for our purposes.
        EXAMPLES:
           sage: [t1,t2]=[CartanType(x) for x in ['A',1],['B',2]]
           sage: CartanType([t1,t2])
           A1xB2
           sage: t = CartanType("A2xB2")
           sage: t == loads(dumps(t))
           True
        """
        self._types = types
        self.affine = False
        self._spaces = [t.root_system().ambient_space() for t in types]
        self._shifts = [sum(l.n for l in self._spaces[:k]) for k in range(len(types))]
        self._rshifts = [sum(l[1] for l in types[:k]) for k in range(len(types))]

    def __repr__(self):
        """
        EXAMPLES:
           sage: ct = CartanType("A2","B2")
           sage: repr(ct)
           'A2xB2'
        """
        names = [t[0]+str(t[1]) for t in self._types]
        return  names[0]+"".join(flatten([["x",t] for t in names[1:]]))

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: ct1 = CartanType(['A',1],['B',2])
            sage: ct2 = CartanType(['B',2],['A',1])
            sage: ct3 = CartanType(['A',4])
            sage: ct = [ct1, ct2, ct3]
            sage: [[x.__cmp__(y) for x in ct] for y in ct]
            [[0, 1, -1], [-1, 0, -1], [1, 1, 0]]
            sage: sorted(ct)
            [['A', 4], A1xB2, B2xA1]
        """
        if isinstance(other, CartanType_simple):
            return 1
        if len(self._types) < len(other._types):
            return -1
        if len(self._types) > len(other._types):
            return 1
        for i in range(len(self._types)):
            if self._types[i] != other._types[i]:
                return self._types[i].__cmp__(other._types[i])
        return 0

    def component_types(self):
        """
        A list of Cartan types making up the reducible type
        EXAMPLES:
          sage: CartanType(['A',2],['B',2]).component_types()
          [['A', 2], ['B', 2]]
        """
        return self._types

    def type(self):
        """
        Returns "reducible" since the type is reducible
        """
        return "reducible"

    def is_finite(self):
        return all(t.is_finite() for t in self.component_types())

    def rank(self):
        """
        Returns the rank of self.
        EXAMPLES:
           sage: CartanType("A2","A1").rank()
           3
        """
        return sum(t.rank() for t in self._types)

    def root_system(self):
        """
        Returns the root system associated to self.

        EXAMPLES:
            sage: CartanType(['A',4]).root_system()
            Root system of type ['A', 4]
        """
        return root_system.RootSystem(self)

    def cartan_matrix(self, subdivide=True):
        """
        Returns the Cartan matrix associated with self. By default
        the Cartan matrix is a subdivided block matrix showing the
        reducibility but the subdivision can be suppressed with
        the option subdivide=False.

        EXAMPLES:
           sage: ct = CartanType("A2","B2")
           sage: ct.cartan_matrix()
           [ 2 -1| 0  0]
           [-1  2| 0  0]
           [-----+-----]
           [ 0  0| 2 -1]
           [ 0  0|-2  2]
           sage: ct.cartan_matrix(subdivide=False)
           [ 2 -1  0  0]
           [-1  2  0  0]
           [ 0  0  2 -1]
           [ 0  0 -2  2]
        """
        return block_diagonal_matrix([t.cartan_matrix() for t in self._types], subdivide=subdivide)

    def is_irreducible(self):
        """
        Report that this Cartan type is not irreducible
        """
        return False

    def dual(self):
        """
        EXAMPLE:
          sage: CartanType("A2xB2").dual()
          A2xC2
        """
        return CartanType([t.dual() for t in self._types])

    def type_string(self):
        """
        EXAMPLE:
          sage: CartanType(['A',2],['B',2]).type_string()
          'type_reducible'
        """
        return "type_reducible"

    def is_affine(self):
        """
        Report that this reducible Cartan type is not affine
        """
        return False


class ambient_space(AmbientSpace):
    """
    EXAMPLES:
      sage: RootSystem("A2xB2").ambient_space()
      Ambient space for the Root system of type A2xB2
    """
    def cartan_type(self):
        """
        EXAMPLES:
           sage: RootSystem("A2xB2").ambient_space().cartan_type()
           A2xB2
        """
        return self.root_system.cartan_type()

    def component_types(self):
        """
        EXAMPLES:
          sage: RootSystem("A2xB2").ambient_space().component_types()
          [['A', 2], ['B', 2]]
        """
        return self.root_system.cartan_type().component_types()

    def dimension(self):
        """
        EXAMPLES:
          sage: RootSystem("A2xB2").ambient_space().dimension()
          5
        """
        return sum(v.dimension() for v in self.ambient_spaces())

    def ambient_spaces(self):
        """
        Returns a list of the irreducible Cartan types of which the
        given reducible Cartan type is a product.

        EXAMPLES:
          sage: RootSystem("A2xB2").ambient_space().ambient_spaces()
          [Ambient space for the Root system of type ['A', 2],
           Ambient space for the Root system of type ['B', 2]]
        """
        return [t.root_system().ambient_space() for t in self.component_types()]

    def inject_weights(self, i, v):
        """
        INPUT:
          i -- an integer in range(self.components)
          v -- a vector in the i-th component weight lattice
        Produces the corresponding element of the lattice.
        EXAMPLES:
          sage: V = RootSystem("A2xB2").ambient_space()
          sage: [V.inject_weights(i,V.ambient_spaces()[i].fundamental_weights()[1]) for i in range(2)]
          [(1, 0, 0, 0, 0), (0, 0, 0, 1, 0)]
          sage: [V.inject_weights(i,V.ambient_spaces()[i].fundamental_weights()[2]) for i in range(2)]
          [(1, 1, 0, 0, 0), (0, 0, 0, 1/2, 1/2)]
        """
        shift = self.root_system.cartan_type()._shifts[i]
        return self._from_dict( dict([(shift+k, c) for (k,c) in v ]))

    def simple_roots(self):
        """
        EXAMPLES:
          sage: RootSystem("A1xA2").ambient_space().simple_roots()
          [(1, -1, 0, 0, 0), (0, 0, 1, -1, 0), (0, 0, 0, 1, -1)]
        """
        ret = []
        for i in range(len(self.component_types())):
            ret.extend(self.inject_weights(i, v) for v in self.ambient_spaces()[i].simple_roots())
        return ret

    def positive_roots(self):
        """
        EXAMPLE:
           sage: RootSystem("A1xA2").ambient_space().positive_roots()
           [(1, -1, 0, 0, 0), (0, 0, 1, -1, 0), (0, 0, 1, 0, -1), (0, 0, 0, 1, -1)]
        """
        ret = []
        for i in range(len(self.component_types())):
            ret.extend(self.inject_weights(i, v) for v in self.ambient_spaces()[i].positive_roots())
        return ret

    def negative_roots(self):
        """
        EXAMPLE:
           sage: RootSystem("A1xA2").ambient_space().negative_roots()
           [(-1, 1, 0, 0, 0), (0, 0, -1, 1, 0), (0, 0, -1, 0, 1), (0, 0, 0, -1, 1)]
        """
        ret = []
        for i in range(len(self.component_types())):
            ret.extend(self.inject_weights(i, v) for v in self.ambient_spaces()[i].negative_roots())
        return ret

    def fundamental_weights(self):
        """
        EXAMPLE:
           sage: RootSystem("A2xB2").ambient_space().fundamental_weights()
           [(1, 0, 0, 0, 0), (1, 1, 0, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 1/2, 1/2)]
        """
        ret = []
        for i in range(len(self.component_types())):
            ret.extend(self.inject_weights(i, v) for v in self.ambient_spaces()[i].fundamental_weights())
        return ret

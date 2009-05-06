from sage.combinat.root_system.cartan_type import CartanType_abstract, CartanType_simple
from sage.matrix.constructor import block_diagonal_matrix
from sage.sets.family import Family
import ambient_space
import sage.combinat.root_system as root_system
from sage.structure.sage_object import SageObject

class CartanType(SageObject, CartanType_abstract):
    r"""
    A class for reducible Cartan types
    """

    def __init__(self, types):
        """
        Reducible root systems are ones that can be factored as
        direct products. Strictly speaking type D2 (corresponding
        to orthogonal groups of degree 4) are reducible since they
        are isomorphic to A1xA1. However type D2 is considered
        irreducible for our purposes.

        INPUT:

        - ``types`` - a list of simple Cartan types

        EXAMPLES::

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
        # fails for dual root systems
        try:
            self._shifts.append(sum(t.root_system().ambient_space().dimension() for t in types))
        except:
            pass
        self._rshifts = [sum(l[1] for l in types[:k]) for k in range(len(types))]
        self.tools = root_system.type_reducible

    def _repr_(self, compact = True): # We should make a consistent choice here
        """
        EXAMPLES::

           sage: CartanType("A2","B2")
           A2xB2

           sage: CartanType("A2",CartanType("F4").dual())
           A2xF4*
        """
        return  "x".join(t._repr_(compact = True) for t in self._types)

    def __cmp__(self, other):
        """
        EXAMPLES:::

            sage: ct1 = CartanType(['A',1],['B',2])
            sage: ct2 = CartanType(['B',2],['A',1])
            sage: ct3 = CartanType(['A',4])
            sage: ct1 == ct1
            True
            sage: ct1 == ct2
            False
            sage: ct1 == ct3
            False
        """
        if isinstance(other, CartanType_simple):
            return 1
        return cmp(self._types, other._types)

    def component_types(self):
        """
        A list of Cartan types making up the reducible type.

        EXAMPLES::

            sage: CartanType(['A',2],['B',2]).component_types()
            [['A', 2], ['B', 2]]
        """
        return self._types

    def type(self):
        """
        Returns "reducible" since the type is reducible.

        EXAMPLES::

            sage: CartanType(['A',2],['B',2]).type()
            'reducible'
        """
        return "reducible"

    def is_finite(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',2],['B',2])
            sage: ct.is_finite()
            True
        """
        return all(t.is_finite() for t in self.component_types())

    def rank(self):
        """
        Returns the rank of self.

        EXAMPLES::

            sage: CartanType("A2","A1").rank()
            3
        """
        return sum(t.rank() for t in self._types)

    def index_set(self):
        """
        Implements :meth:`CartanType_abstract.index_set`.

        For the moment, the index set is always of the form `{1,\dots,n}`.

        EXAMPLES::

            sage: CartanType("A2","A1").index_set()
            [1, 2, 3]
        """
        return range(1, self.rank()+1)

    def root_system(self):
        """
        Returns the root system associated to self.

        EXAMPLES::

            sage: CartanType(['A',4], ['B', 2]).root_system()
            Root system of type A4xB2

        """
        import root_system
        return root_system.RootSystem(self)

    def cartan_matrix(self, subdivide=True):
        """
        Returns the Cartan matrix associated with self. By default
        the Cartan matrix is a subdivided block matrix showing the
        reducibility but the subdivision can be suppressed with
        the option subdivide=False.

        EXAMPLES::

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

    def dynkin_diagram(t):
        """
        Returns a Dynkin diagram for type reducible.

        EXAMPLES::

            sage: dd = CartanType("A2xB2xF4").dynkin_diagram()
            sage: dd
            O---O
            1   2
            O=>=O
            3   4
            O---O=>=O---O
            5   6   7   8
            A2xB2xF4
            sage: dd.edges()
            [(1, 2, 1), (2, 1, 1), (3, 4, 2), (4, 3, 1), (5, 6, 1), (6, 5, 1), (6, 7, 2), (7, 6, 1), (7, 8, 1), (8, 7, 1)]

            sage: CartanType("F4xA2").dynkin_diagram()
            O---O=>=O---O
            1   2   3   4
            O---O
            5   6
            F4xA2

        """
        from dynkin_diagram import DynkinDiagram, DynkinDiagram_class
        g = DynkinDiagram_class(t)
        for i in range(len(t._types)):
            for [e1, e2, l] in DynkinDiagram(t._types[i]).edges():
                shift = t._rshifts[i]
                g.add_edge(e1+shift, e2+shift, label=l)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns an ascii art representation of this reducible Cartan type

        EXAMPLES::
            sage: print CartanType("F4xA2").ascii_art(label = lambda x: x+2)
            O---O=>=O---O
            3   4   5   6
            O---O
            7   8
        """
        types = self.component_types()
        return "\n".join(types[i].ascii_art(label = lambda x: label(x+self._rshifts[i]))
                         for i in range(len(types)))

    def is_irreducible(self):
        """
        Report that this Cartan type is not irreducible.

        EXAMPLES::

            sage: ct = CartanType(['A',2],['B',2])
            sage: ct.is_irreducible()
            False
        """
        return False

    def dual(self):
        """
        EXAMPLES::

            sage: CartanType("A2xB2").dual()
            A2xC2
        """
        return CartanType([t.dual() for t in self._types])

    def is_affine(self):
        """
        Report that this reducible Cartan type is not affine

        EXAMPLES::

            sage: CartanType(['A',2],['B',2]).is_affine()
            False
        """
        return False




class AmbientSpace(ambient_space.AmbientSpace):
    """
    EXAMPLES::

        sage: RootSystem("A2xB2").ambient_space()
        Ambient space of the Root system of type A2xB2

    """
    def cartan_type(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().cartan_type()
            A2xB2
        """
        return self.root_system.cartan_type()

    def component_types(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().component_types()
            [['A', 2], ['B', 2]]
        """
        return self.root_system.cartan_type().component_types()

    def dimension(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().dimension()
            5
        """
        return sum(v.dimension() for v in self.ambient_spaces())

    def ambient_spaces(self):
        """
        Returns a list of the irreducible Cartan types of which the
        given reducible Cartan type is a product.

        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().ambient_spaces()
            [Ambient space of the Root system of type ['A', 2],
             Ambient space of the Root system of type ['B', 2]]
        """
        return [t.root_system().ambient_space() for t in self.component_types()]

    def inject_weights(self, i, v):
        """
        Produces the corresponding element of the lattice.

        INPUT:

        - ``i`` - an integer in range(self.components)

        - ``v`` - a vector in the i-th component weight lattice

        EXAMPLES::

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
        EXAMPLES::

            sage: RootSystem("A1xB2").ambient_space().simple_roots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 1)}
        """
        res = []
        for i, ambient_space in enumerate(self.ambient_spaces()):
            res.extend(self.inject_weights(i, v) for v in ambient_space.simple_roots())
        return Family(dict([i,res[i-1]] for i in range(1,len(res)+1)))

    def simple_coroots(self):
        """
        EXAMPLES:

            sage: RootSystem("A1xB2").ambient_space().simple_coroots()
            Finite family {1: (1, -1, 0, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 2)}
        """
        cr = []
        for i, ambient_space in enumerate(self.ambient_spaces()):
            cr.extend(self.inject_weights(i, v) for v in ambient_space.simple_coroots())
        return Family(dict([i,cr[i-1]] for i in range(1,len(cr)+1)))

    def positive_roots(self):
        """
        EXAMPLES::

            sage: RootSystem("A1xA2").ambient_space().positive_roots()
            [(1, -1, 0, 0, 0), (0, 0, 1, -1, 0), (0, 0, 1, 0, -1), (0, 0, 0, 1, -1)]
        """
        res = []
        for i, ambient_space in enumerate(self.ambient_spaces()):
            res.extend(self.inject_weights(i, v) for v in ambient_space.positive_roots())
        return res

    def negative_roots(self):
        """
        EXAMPLES::

            sage: RootSystem("A1xA2").ambient_space().negative_roots()
            [(-1, 1, 0, 0, 0), (0, 0, -1, 1, 0), (0, 0, -1, 0, 1), (0, 0, 0, -1, 1)]
        """
        ret = []
        for i, ambient_space in enumerate(self.ambient_spaces()):
            ret.extend(self.inject_weights(i, v) for v in ambient_space.negative_roots())
        return ret

    def fundamental_weights(self):
        """
        EXAMPLES::

            sage: RootSystem("A2xB2").ambient_space().fundamental_weights()
            Finite family {1: (1, 0, 0, 0, 0), 2: (1, 1, 0, 0, 0), 3: (0, 0, 0, 1, 0), 4: (0, 0, 0, 1/2, 1/2)}
        """
        fw = []
        for i, ambient_space in enumerate(self.ambient_spaces()):
            fw.extend(self.inject_weights(i, v) for v in ambient_space.fundamental_weights())
        return Family(dict([i,fw[i-1]] for i in range(1,len(fw)+1)))


CartanType.AmbientSpace = AmbientSpace

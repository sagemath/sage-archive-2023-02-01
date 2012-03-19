"""
Root system data for relabelled Cartan types
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.root_system.cartan_type import CartanType_abstract
from sage.sets.family import FiniteFamily
import sage

class CartanType(UniqueRepresentation, SageObject, CartanType_abstract):
    r"""
    A class for relabelled Cartan types
    """

    @staticmethod
    def __classcall__(cls, type, relabelling):
        """
        This standardizes the input of the constructor to ensure
        unique representation.

        EXAMPLES::

            sage: ct1 = CartanType(['B',2]).relabel({1:2, 2:1})    # indirect doctest
            sage: ct2 = CartanType(['B',2]).relabel(lambda x: 3-x)
            sage: ct3 = CartanType(['B',2]).relabel({1:3, 2: 4})
            sage: ct4 = CartanType(['D',4]).relabel(lambda x: 3-x)
            sage: ct1 == ct2
            True
            sage: ct1 == ct3
            False
            sage: ct1 == ct4
            False

        """

        if isinstance(relabelling, (list, tuple, dict, FiniteFamily)):
            # allows for using relabellings with more entries than in the index_set
            # and by the way makes a copy of relabelling
            relabelling = dict( (i, relabelling[i]) for i in type.index_set() )
        else:
            relabelling = dict( (i, relabelling(i)) for i in type.index_set() )

        if isinstance(type, CartanType): # type is already a relabelled type
            relabelling = dict( (i, relabelling[type._relabelling[i]]) for i in type._type.index_set() )
            type = type._type

        if all( relabelling[i] == i for i in type.index_set() ):
            return type

        relabelling = FiniteFamily(relabelling) # Hack to emulate a frozendict which would be hashable!!!!
        return super(CartanType, cls).__classcall__(cls, type, relabelling)


    def __init__(self, type, relabelling):
        """
        INPUT:
         - ``type``: a Cartan type
         - ``relabelling``: a function (or a list, or a dictionary)

        Returns an isomorphic Cartan type obtained by relabelling the
        nodes of the dynkin diagram. Namely the node with label ``i``
        is relabelled ``f(i)`` (or, by ``f[i]`` if f is a list or
        dictionary).

        EXAMPLES:

        We take the Cartan type `B_4`::

            sage: T = CartanType(['B',4])
            sage: T.dynkin_diagram()
            O---O---O=>=O
            1   2   3   4
            B4

        And relabel its nodes::

            sage: cycle = {1:2, 2:3, 3:4, 4:1}

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()
            O---O---O=>=O
            2   3   4   1
            B4 relabelled by {1: 2, 2: 3, 3: 4, 4: 1}
            sage: sorted(T.dynkin_diagram().edges())
            [(1, 4, 1), (2, 3, 1), (3, 2, 1), (3, 4, 1), (4, 1, 2), (4, 3, 1)]

        Multiple relabelling are recomposed into a single one::

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()
            O---O---O=>=O
            3   4   1   2
            B4 relabelled by {1: 3, 2: 4, 3: 1, 4: 2}

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()
            O---O---O=>=O
            4   1   2   3
            B4 relabelled by {1: 4, 2: 1, 3: 2, 4: 3}

        And trivial relabelling are honoured nicely::

            sage: T = T.relabel(cycle)
            sage: T.dynkin_diagram()
            O---O---O=>=O
            1   2   3   4
            B4


        TESTS::

            sage: T = CartanType(['B',4]).relabel(cycle)
            sage: T == loads(dumps(T))
            True
        """
        assert isinstance(relabelling, FiniteFamily)
        self._type = type
        self._relabelling = relabelling._dictionary
        self._index_set = sorted(relabelling[i] for i in type.index_set())
        if type.is_affine():
            self._add_abstract_superclass(CartanType_affine)

    def _repr_(self, compact = False):
        """
        EXAMPLES::

           sage: CartanType(['F', 4]).relabel(lambda x: 5-x)
           ['F', 4] relabelled by {1: 4, 2: 3, 3: 2, 4: 1}

           sage: CartanType(['F', 4]).relabel(lambda x: 5-x)._repr_(compact = True)
           'F4 relabelled by {1: 4, 2: 3, 3: 2, 4: 1}'
        """
        return self._type._repr_(compact = compact)+" relabelled by %s"%self._relabelling

    def ascii_art(self, label = lambda i: i):
        """
        Returns an ascii art representation of this Cartan type

        EXAMPLES::

            sage: print CartanType(["G", 2]).relabel({1:2,2:1}).ascii_art()
              3
            O=<=O
            2   1
            sage: print CartanType(["B", 3, 1]).relabel([1,3,2,0]).ascii_art()
                O 1
                |
                |
            O---O=>=O
            3   2   0
            sage: print CartanType(["F", 4, 1]).relabel(lambda n: 4-n).ascii_art()
            O---O---O=>=O---O
            4   3   2   1   0
        """
        return self._type.ascii_art(label = lambda i: label(self._relabelling[i]))

    def dynkin_diagram(self):
        """
        Returns the dynkin diagram for this Cartan type.

        EXAMPLES::

            sage: CartanType(["G", 2]).relabel({1:2,2:1}).dynkin_diagram()
              3
            O=<=O
            2   1
            G2 relabelled by {1: 2, 2: 1}

        TESTS:

        To be compared with the examples in :meth:`ascii_art`::

            sage: sorted(CartanType(["G", 2]).relabel({1:2,2:1}).dynkin_diagram().edges())
            [(1, 2, 3), (2, 1, 1)]
            sage: sorted(CartanType(["B", 3, 1]).relabel([1,3,2,0]).dynkin_diagram().edges())
            [(0, 2, 1), (1, 2, 1), (2, 0, 2), (2, 1, 1), (2, 3, 1), (3, 2, 1)]
            sage: sorted(CartanType(["F", 4, 1]).relabel(lambda n: 4-n).dynkin_diagram().edges())
            [(0, 1, 1), (1, 0, 1), (1, 2, 1), (2, 1, 2), (2, 3, 1), (3, 2, 1), (3, 4, 1), (4, 3, 1)]

        """
        # Maybe we want to move this up as a relabel method for dynkin diagram
        # We will have to be careful setting the cartan type of the result though
        from copy import copy
        result = copy(self._type.dynkin_diagram())
        # relabelling in place allows to keep the extra dynkin diagram structure
        super(result.__class__, result).relabel(self._relabelling, inplace = True)
        result._cartan_type = self
        return result

    def is_irreducible(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_irreducible()
           True
        """
        return self._type.is_irreducible()

    def is_finite(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_finite()
           True
        """
        return self._type.is_finite()

    def is_crystalographic(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_crystalographic()
           True
        """
        return self._type.is_crystalographic()

    def is_affine(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_affine()
           False
        """
        return self._type.is_affine()

    def rank(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.rank()
           4
        """
        return self._type.rank()

    def index_set(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.index_set()
           [1, 2, 3, 4]
        """
        return self._index_set

    def dual(self):
        """
        Implements :meth:`sage.combinat.root_system.cartan_type.CartanType_abstract.dual`,
        using that taking the dual and relabelling are commuting operations.

        EXAMPLES::

            sage: T = CartanType(["BC",3, 2])
            sage: cycle = {1:2, 2:3, 3:0, 0:1}
            sage: T.relabel(cycle).dual().dynkin_diagram()
            O=>=O---O=>=O
            1   2   3   0
            BC3~* relabelled by {0: 1, 1: 2, 2: 3, 3: 0}
            sage: T.dual().relabel(cycle).dynkin_diagram()
            O=>=O---O=>=O
            1   2   3   0
            BC3~* relabelled by {0: 1, 1: 2, 2: 3, 3: 0}
        """
        return self._type.dual().relabel(self._relabelling)


###########################################################################
class CartanType_affine(sage.combinat.root_system.cartan_type.CartanType_affine):
    def classical(self):
        """
        Returns the classical Cartan type associated with self (which should be affine)

        EXAMPLES::

            sage: A41 = CartanType(['A',4,1])
            sage: A41.dynkin_diagram()
            0
            O-----------+
            |           |
            |           |
            O---O---O---O
            1   2   3   4
            A4~

            sage: T = A41.relabel({0:1, 1:2, 2:3, 3:4, 4:0})
            sage: T
            ['A', 4, 1] relabelled by {0: 1, 1: 2, 2: 3, 3: 4, 4: 0}
            sage: T.dynkin_diagram()
            1
            O-----------+
            |           |
            |           |
            O---O---O---O
            2   3   4   0
            A4~ relabelled by {0: 1, 1: 2, 2: 3, 3: 4, 4: 0}

            sage: T0 = T.classical()
            sage: T0
            ['A', 4] relabelled by {1: 2, 2: 3, 3: 4, 4: 0}
            sage: T0.dynkin_diagram()
            O---O---O---O
            2   3   4   0
            A4 relabelled by {1: 2, 2: 3, 3: 4, 4: 0}

        """
        return self._type.classical().relabel(self._relabelling)

    def special_node(self):
        r"""
        Returns a special node of the Dynkin diagram

        .. seealso:: :meth:`~sage.combinat.root_system.CartanType_affine.special_node`

        It is obtained by relabelling of the special node of the non
        relabelled Dynkin diagram.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).special_node()
            0
            sage: CartanType(['B', 3, 1]).relabel({1:2, 2:3, 3:0, 0:1}).special_node()
            1
        """
        return self._relabelling[self._type.special_node()]

    def is_untwisted_affine(self):
        """
        Implements :meth:'CartanType_affine.is_untwisted_affine`

        A relabelled Cartan type is untwisted affine if the original is.

        EXAMPLES::

            sage: CartanType(['B', 3, 1]).relabel({1:2, 2:3, 3:0, 0:1}).is_untwisted_affine()
            True

        """
        return self._type.is_untwisted_affine()


#class ambient_space(AmbientSpace):
# todo?

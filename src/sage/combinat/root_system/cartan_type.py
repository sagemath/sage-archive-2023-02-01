"""
Cartan types
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from sage.combinat import root_system
from sage.rings.all import ZZ
from sage.structure.sage_object import SageObject

# TODO:
# Get rid of almost all runtype type checking by extending the class hierarchy with:
#  - type_relabel.CartanType (in type_relabel module, together with relabelled AmbientSpace)
#  - CartanType_crystalographic
#  - CartanType_simply_laced
#  - type_A.CartanType (in the type_A module)
#  - ...
#  - CartanType_untwisted_affine
#  - type_A_affine.CartanType
#  - type_BC_affine.CartanType
# ...
# Implement the Kac conventions by relabeling/dual/... of the above
# Implement coxeter diagrams for non crystalographic
# Implement dual ambient space

# Intention: we want simultaneously CartanType to be a factory for
# the various subtypes of CartanType_abstract, as in:
#     CartanType(["A",4,1])
# and to behaves as a "module" for some extra utilities:
#     CartanType.samples()
#
# Implementation: CartanType is the unique instance of this class
# CartanTypeFactory. Is there a better/more standard way to do it?
class CartanTypeFactory(SageObject):
    def __call__(self, *args):
        """
        Returns an object corresponding to the Cartan type t.

        INPUT: [letter, rank] where letter is one of
        'A','B','C','D','E','F','G' and rank is the rank. An alternative
        string notation is allowed. A third optional parameter is permitted
        for affine types. Reducible types may be entered by giving a list
        of irreducible types or by a single string

        EXAMPLES::

            sage: CartanType(['A',4])
            ['A', 4]
            sage: CartanType("A4")
            ['A', 4]
            sage: CartanType(['A',2],['B',2])
            A2xB2
            sage: CartanType(['A',2],['B',2]).is_reducible()
            True
            sage: CartanType("A2xB2")
            A2xB2
            sage: CartanType("A2","B2") == CartanType("A2xB2")
            True
            sage: CartanType(['A',4,1])
            ['A', 4, 1]
            sage: CartanType(['A',4,1]).is_affine()
            True
        """
        if len(args) == 1:
            t = args[0]
        else:
            t = args
        if isinstance(t, CartanType_abstract):
            return t

        if type(t)==str:
            if "x" in t:
                return root_system.type_reducible.CartanType([CartanType(u) for u in t.split("x")])
            else:
                return CartanType([t[0], eval(t[1:])])

        t = list(t)

        if type(t[0]) == str and t[1] in ZZ:
            if len(t) == 2:
                return CartanType_simple_finite(t)
            elif len(t) == 3:
                return CartanType_simple_affine(t)

        return root_system.type_reducible.CartanType([ CartanType(subt) for subt in t ])

    def samples(self, finite=False, affine=False, crystalographic=False):
        """
        Returns a sample of the implemented cartan types

        With finite=True resp. affine=True, one can restrict to finite
        resp. affine only cartan types

        EXAMPLES::

            sage: CartanType.samples(finite=True)
            [['A', 1], ['A', 5], ['B', 5], ['C', 5], ['D', 5], ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['I', 5], ['H', 3], ['H', 4]]

        ::

            sage: CartanType.samples(affine=True)
            [['A', 1, 1], ['A', 5, 1], ['B', 5, 1], ['C', 5, 1], ['D', 5, 1], ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1], ['A', 2, 2], ['A', 10, 2], ['A', 9, 2], ['D', 5, 2], ['D', 4, 3], ['E', 6, 2]]

        ::

            sage: CartanType.samples()
            [['A', 1], ['A', 5], ['B', 5], ['C', 5], ['D', 5], ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['I', 5], ['H', 3], ['H', 4], ['A', 1, 1], ['A', 5, 1], ['B', 5, 1], ['C', 5, 1], ['D', 5, 1], ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1], ['A', 2, 2], ['A', 10, 2], ['A', 9, 2], ['D', 5, 2], ['D', 4, 3], ['E', 6, 2]]
            sage: CartanType.samples(crystalographic=True)
            [['A', 1], ['A', 5], ['B', 5], ['C', 5], ['D', 5], ['E', 6], ['E', 7], ['E', 8], ['F', 4], ['G', 2], ['A', 1, 1], ['A', 5, 1], ['B', 5, 1], ['C', 5, 1], ['D', 5, 1], ['E', 6, 1], ['E', 7, 1], ['E', 8, 1], ['F', 4, 1], ['G', 2, 1], ['A', 2, 2], ['A', 10, 2], ['A', 9, 2], ['D', 5, 2], ['D', 4, 3], ['E', 6, 2]]
        """
        if crystalographic:
            return [ t for t in CartanType.samples(finite=finite, affine=affine) if t.is_crystalographic() ]
        if finite:
            return([CartanType(t) for t in [["A", 1], ["A", 5], ["B", 5], ["C", 5], ["D", 5],
                                            ["E", 6], ["E", 7], ["E", 8],
                                            ["F", 4],
                                            ["G", 2],
                                            ["I", 5],
                                            ["H", 3], ["H", 4]]])
        elif affine:
            return([t.affine() for t in CartanType.samples(finite=True, crystalographic=True)] +
                   [CartanType(t) for t in [["A", 2, 2], ["A", 10, 2], ["A", 9, 2],
                                            ["D", 5, 2],
                                            ["D", 4, 3],
                                            ["E", 6, 2]]])
        else:
            return CartanType.samples(finite=True) + CartanType.samples(affine=True)

CartanType = CartanTypeFactory()

class CartanType_abstract(SageObject):
    r"""
    Abstract class for cartan types

    Subclasses should implement:

    - type()

    - type_string()

    - dynkin_diagram()

    - cartan_matrix()

    - is_finite()

    - is_affine()

    - is_irreducible()
    """

    def type(self):
        r"""
        Returns the type of self, or None if unknown. This method should be
        overridden in any subclass.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.type() is None
            True
        """
        return None

    def rank(self):
        """
        Returns the rank of self.

        EXAMPLES::

            sage: CartanType(['A', 4]).rank()
            4
            sage: CartanType(['A', 7, 2]).rank()
            4
            sage: CartanType(['I', 8]).rank()
            2
        """
        raise NotImplementedError

    def dual(self):
        """
        Returns the dual cartan type, possibly just as a formal dual.

        EXAMPLES::

            sage: CartanType(['F',4]).dual()
            ['F', 4]^*
        """
        return root_system.type_dual.CartanType(self)

    def is_reducible(self):
        """
        Report whether the root system is reducible (i.e. not simple), that
        is whether it can be factored as a product of root systems.

        EXAMPLES::

            sage: CartanType("A2xB3").is_reducible()
            True
            sage: CartanType(['A',2]).is_reducible()
            False
        """
        return not self.is_irreducible()

    def is_irreducible(self):
        """
        Report whether this Cartan type is irreducible (i.e. simple). This
        should be overridden in any subclass.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.is_irreducible()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_finite(self):
        """
        Returns whether this Cartan type is finite. This should be
        overridden in any subclass.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.is_irreducible()
            Traceback (most recent call last):
            ...
            NotImplementedError

        ::

            sage: CartanType(['A',4]).is_finite()
            True
            sage: CartanType(['A',4, 1]).is_finite()
            False
        """
        raise NotImplementedError

    def is_affine(self):
        """
        Returns whether self is affine.

        EXAMPLES::

            sage: CartanType(['A', 3]).is_affine()
            False
            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """
        raise NotImplementedError

    def is_crystalographic(self):
        """
        Returns whether this Cartan type is simple laced

        EXAMPLES::

            sage: [ [t, t.is_crystalographic() ] for t in CartanType.samples(finite=True) ]
            [[['A', 1], True], [['A', 5], True],
            [['B', 5], True], [['C', 5], True], [['D', 5], True],
            [['E', 6], True], [['E', 7], True], [['E', 8], True],
            [['F', 4], True], [['G', 2], True],
            [['I', 5], False], [['H', 3], False], [['H', 4], False]]

        TESTS::

            sage: all(t.is_crystalographic() for t in CartanType.samples(affine=True))
            True
        """
        raise NotImplementedError

    def is_simple_laced(self):
        """
        Returns whether this Cartan type is simple laced

        EXAMPLES::

            sage: [ [t, t.is_simply_laced() ] for t in CartanType.samples() ]
            [[['A', 1], True], [['A', 5], True],
            [['B', 5], False], [['C', 5], False], [['D', 5], True],
            [['E', 6], True], [['E', 7], True], [['E', 8], True],
            [['F', 4], False], [['G', 2], False], [['I', 5], False], [['H', 3], False], [['H', 4], False],
            [['A', 1, 1], False], [['A', 5, 1], True],
            [['B', 5, 1], False], [['C', 5, 1], False], [['D', 5, 1], True],
            [['E', 6, 1], True], [['E', 7, 1], True], [['E', 8, 1], True],
            [['F', 4, 1], False], [['G', 2, 1], False],
            [['A', 2, 2], False], [['A', 10, 2], False], [['A', 9, 2], False], [['D', 5, 2], False], [['D', 4, 3], False], [['E', 6, 2], False]]
        """
        raise NotImplementedError

    def index_set(self):
        """
        Returns the index set for self.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).index_set()
            [0, 1, 2, 3]
            sage: CartanType(['D', 4]).index_set()
            [1, 2, 3, 4]
        """
        n = self.rank()
        if self.is_affine():
            return range(n+1)
        else:
            return range(1, n+1)

    def root_system(self):
        """
        Returns the root system associated to self.

        EXAMPLES::

            sage: CartanType(['A',4]).root_system()
            Root system of type ['A', 4]
        """
        return root_system.root_system.RootSystem(self)

# Maybe we want a separate class for affine

class CartanType_simple(CartanType_abstract):
    def __cmp__(self, other):
        """
        TESTS::

            sage: ct1 = CartanType(['A',4])
            sage: ct2 = CartanType(['A',4])
            sage: ct3 = CartanType(['A',5])
            sage: ct1 == ct2
            True
            sage: ct1 != ct3
            True
        """
        if other.__class__ != self.__class__:
            return cmp(self.__class__, other.__class__)
        if other.letter != self.letter:
            return cmp(self.letter, other.letter)
        return cmp(self.n, other.n)

    def __hash__(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',2])
            sage: hash(ct) #random
            -5684143898951441983
        """
        return hash(str(self))

    def __getitem__(self, x):
        """
        EXAMPLES::

            sage: t = CartanType(['A', 3, 1])
            sage: t[0]
            'A'
            sage: t[1]
            3
            sage: t[2]
            1
            sage: t[3]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.t[x]

    def is_irreducible(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3]).is_irreducible()
            True
        """
        return True

    def dynkin_diagram(self):
        """
        Returns the Dynkin diagram associated with self.

        EXAMPLES::

            sage: CartanType(['A',4]).dynkin_diagram()
            O---O---O---O
            1   2   3   4
            A4
        """
        return root_system.dynkin_diagram.DynkinDiagram(self)

    def cartan_matrix(self):
        """
        Returns the Cartan matrix associated with self.

        EXAMPLES::

            sage: CartanType(['A',4]).cartan_matrix()
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]
        """
        return root_system.cartan_matrix.cartan_matrix(self)

    def type(self):
        """
        Returns the type of self.

        EXAMPLES::

            sage: CartanType(['A', 4]).type()
            'A'
            sage: CartanType(['A', 4, 1]).type()
            'A'
        """
        return self.letter

    def dual(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).dual()
            ['A', 3]
            sage: CartanType(["B", 3]).dual()
            ['C', 3]
        """
        if self.type() in ["A", "D", "E"]:
            return self
        else:
            return CartanType_abstract.dual(self)


class CartanType_simple_finite(CartanType_simple):
    r"""
    A class for finite simple Cartan types
    """

    def __init__(self, t):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4])
            sage: ct == loads(dumps(ct))
            True
        """
        assert(len(t) == 2)
        assert(t[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'])
        assert(t[1] in ZZ and t[1] >= 0)
        if t[0] in ['B', 'C']:
            assert(t[1] >= 2)
        if t[0] == 'D':
            assert(t[1] >= 3)
        if t[0] == 'E':
            assert(t[1] <= 8)
        if t[0] == 'F':
            assert(t[1] <= 4)
        if t[0] == 'G':
            assert(t[1] <= 2)
        if t[0] == 'H':
            assert(t[1] <= 4)

        self.t = t
        self.letter = t[0]
        self.n = t[1]

        # Get the python module containing type-specific information,
        # if it exists
        self.tools = getattr(root_system,
                             "type_"+self.type(),
                             root_system.type_None)


    def __repr__(self):
        """
        TESTS::

            sage: ct = CartanType(['A',3])
            sage: repr(ct)
            "['A', 3]"
        """
        return "['%s', %s]"%(self.letter, self.n)

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CartanType(['A',4]))
            2
        """
        return 2

    def rank(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).rank()
            3
        """
        if self.letter == "I":
            return 2
        else:
            return self.n

    def is_finite(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).is_finite()
            True
        """
        return True

    def is_affine(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).is_affine()
            False
        """
        return False

    def is_crystalographic(self):
        """
        EXAMPLES::

            sage: CartanType(["A", 3]).is_crystalographic()
            True
            sage: CartanType(["I", 2]).is_crystalographic()
            False
        """
        return self.letter in ["A", "B", "C", "D", "E", "F", "G"]

    def is_simply_laced(self):
        """
        EXAMPLES::

            sage: CartanType(['A',3]).is_simply_laced()
            True
            sage: CartanType(['B',3]).is_simply_laced()
            False
        """
        return self.letter in  ["A", "D", "E"]

    def affine(self):
        """
        Returns the corresponding untwisted affine Cartan type

        EXAMPLES::

            sage: CartanType(['A',3]).affine()
            ['A', 3, 1]
        """
        return CartanType([self.letter, self.n, 1])

    def dual(self):
        """
        EXAMPLES::

            sage: CartanType(['A',3]).dual()
            ['A', 3]
            sage: CartanType(['D',4]).dual()
            ['D', 4]
            sage: CartanType(['E',8]).dual()
            ['E', 8]
            sage: CartanType(['B',3]).dual()
            ['C', 3]
            sage: CartanType(['C',2]).dual()
            ['B', 2]
        """
        if self.type() == "B":
            return CartanType(["C",self.n])
        elif self.type() == "C":
            return CartanType(["B",self.n])
        else:
            return CartanType_simple.dual(self)

##########################################################################
class CartanType_simple_affine(CartanType_simple):
    r"""
    A class for affine simple Cartan types
    """
    def __init__(self, t):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4])
            sage: ct == loads(dumps(ct))
            True
        """
        assert(len(t) == 3)
        assert(t[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        assert(t[1] in ZZ and t[1] >= 1)
        assert(t[2] in [1,2,3])
        if t[0] in ['B', 'C']:
            assert(t[1] >= 2)
        if t[0] == 'D':
            assert(t[1] >= 4)
        if t[0] == 'E':
            assert(t[1] <= 8)
        if t[0] == 'F':
            assert(t[1] <= 4)
        if t[0] == 'G':
            assert(t[1] <= 2)
        if t[2] == 3:
            assert(t[0] == "D")
            assert(t[1] == 4)
        if t[2] == 2:
            assert(t[0] in ['A', 'D', 'E'])
            if t[0] == 'E':
                assert(t[1] == 6)

        self.t = t
        self.letter = t[0]
        self.n = t[1]
        self.affine = t[2]

        # Get the python module containing type-specific information,
        # if it exists
        self.tools = getattr(root_system,
                             'type_'+self.type(),
                             root_system.type_None);

    def __repr__(self):
        """
        TESTS::

            sage: ct = CartanType(['A',3, 1])
            sage: repr(ct)
            "['A', 3, 1]"
        """
        return "['%s', %s, %s]"%(self.letter, self.n, self.affine)

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: ct1 = CartanType(['A',3, 1])
            sage: ct2 = CartanType(['B',3, 1])
            sage: ct3 = CartanType(['A',3])
            sage: ct1 == ct1
            True
            sage: ct1 == ct2
            False
            sage: ct1 == ct3
            False
        """
        c = CartanType_simple.__cmp__(self, other)
        if c != 0:
            return c
        else:
            return cmp(self.affine, other.affine)

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CartanType(['A',4,1]))
            3
        """
        return 3

    def rank(self):
        """
        EXAMPLES::

            sage: CartanType(['D', 4, 3]).rank()
            3
            sage: CartanType(['B', 4, 1]).rank()
            4
        """
        if self.affine == 3 and self.letter == 'D':
            return self.n-1
        elif self.affine == 2 and self.letter == 'A':
            ## FIXME: check in the literature what should be the
            ## appropriate definition for rank
            return int(self.n+1)/2
        else:
            return self.n

    def is_finite(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_finite()
            False
        """
        return False

    def is_affine(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """
        return True

    def is_crystalographic(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_crystalographic()
            True
        """
        return True

    def is_simply_laced(self):
        """
        EXAMPLES::

            sage: CartanType(['A', 3, 1]).is_simply_laced()
            True
            sage: CartanType(['D', 4, 3]).is_simply_laced()
            False
            sage: CartanType(['D', 4, 1]).is_simply_laced()
            True
            sage: CartanType(['B', 4, 1]).is_simply_laced()
            False
        """
        if self.affine != 1:
            return False
        if self.letter == "A":
            return self.n > 1
        return self.letter in  ["D", "E"]

    def classical(self):
        r"""
        Returns the classical Cartan type associated with self (which
        should be affine)

        Caveat: only implemented for untwisted

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).classical()
            ['A', 3]
            sage: CartanType(['B', 3, 1]).classical()
            ['B', 3]
        """

        if self.affine == 1:
            return CartanType([self.letter,self.n])
        else:
            raise notImplementedError

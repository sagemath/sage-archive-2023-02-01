r"""
Spin Crystals

These are the crystals associated with the three spin
representations: the spin representations of odd
orthogonal groups (or rather their double covers); and
the + and - spin representations of the even orthogonal
groups.

We follow Kashiwara and Nakashima (Journal of Algebra
165, 1994) in representing the elements of the spin Crystal
by sequences of signs +/-. Two other representations
are available as attributes internal_repn and signature
of the crystal element.

* A numerical internal representation, an integer N
such that if N-1 is written in binary and the 1's are
replaced by -, the 0's by +

* The signature, which is a list in which + is replaced
by +1 and - by -1.
"""

#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
#                          Daniel Bump    <bump at match.stanford.edu>
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
#****************************************************************************

## TODO: proper latex'ing of spins

from sage.rings.integer        import Integer
from sage.structure.element    import Element
from sage.combinat.cartan_type import CartanType
from crystals                  import Crystal, ClassicalCrystal, CrystalElement


#########################
# Type B spin
#########################

def CrystalOfSpins(type):
    r"""
    Return the spin crystal of the given type B.

    INPUT:
        ['B',n] -- A CartanType of type B.

    This is a combinatorial model for the crystal with
    highest weight $Lambda_n$ (the n-th fundamental
    weight). It has $2^n$ elements, here called Spins.
    See also CrystalOfLetters, CrystalOfSpinsPlus
    and CrystalOfSpinsMinus.

    TESTS:

        sage: C = CrystalOfSpins(['B',3])
        sage: C.list()
        [+++, ++-, +-+, -++, +--, -+-, --+, ---]
        sage: C.check()
        True

        sage: [x.signature() for x in C]
        [[1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]]

        sage: len(TensorProductOfCrystals(C,C,generators=[[C(1),C(1)]]).list())
        35
    """
    type = CartanType(type)
    if type[0] == 'B':
	return GenericCrystalOfSpins(type, Spin_crystal_type_B_element, 1)

    else:
	raise NotImplementedError

#########################
# Type D spins
#########################

def CrystalOfSpinsPlus(type):
    r"""
    Return the plus spin crystal of the given type D.

    INPUT:
        ['D',n] -- A CartanType of type D.

    This is the crystal with highest weight $Lambda_n$
    (the n-th fundamental weight).

    TESTS:

        sage: D = CrystalOfSpinsPlus(['D',4])
        sage: D.list()
        [++++, ++--, +-+-, -++-, +--+, -+-+, --++, ----]

#       TODO: update for D4
#        sage: [x.signature() for x in D]
#        [[1, 1, 1],
#         [1, 1, -1],
#         [1, -1, 1],
#         [-1, 1, 1],
#         [1, -1, -1],
#         [-1, 1, -1],
#         [-1, -1, 1],
#         [-1, -1, -1]]

        sage: D.check()
	True
    """
    type = CartanType(type)
    if type[0] == 'D':
	return GenericCrystalOfSpins(type, Spin_crystal_type_D_element, 1)
    else:
	raise NotImplementedError

def CrystalOfSpinsMinus(type):
    r"""
    Return the minus spin crystal of the given type D.

    INPUT:
        ['D',n] -- A CartanType of type D.

    This is the crystal with highest weight $Lambda_{n-1}$
    (the (n-1)-st fundamental weight).

    Caution: to get the highest weight vector use
    E = CrystalOfSpinsMinus(['D',4])
    E(2) NOT E(1).

    TESTS:

        sage: D = CrystalOfSpinsMinus(['D',4])
        sage: D.list()
        [+++-, ++-+, +-++, -+++, +---, -+--, --+-, ---+]

#       TODO: update for D4
#        sage: [x.signature() for x in C]
#        [[1, 1, 1],
#         [1, 1, -1],
#         [1, -1, 1],
#         [-1, 1, 1],
#         [1, -1, -1],
#         [-1, 1, -1],
#         [-1, -1, 1],
#         [-1, -1, -1]]

        sage: D.check()
	True
    """
    type = CartanType(type)
    if type[0] == 'D':
	return GenericCrystalOfSpins(type, Spin_crystal_type_D_element, 2)
    else:
	raise NotImplementedError

class GenericCrystalOfSpins(ClassicalCrystal):
    r"""
    Derived from ClassicalCrystalOfLetters
    """
    def __init__(self, type, element_class, generator):
        self.cartanType = CartanType(type)
        if type[0] == 'B':
            self._name = "The crystal of spins for type %s"%type
        elif generator == 1:
            self._name = "The plus crystal of spins for type %s"%type
        else:
            self._name = "The minus crystal of spins for type %s"%type

        self.index_set = self.cartanType.index_set()
        self.element_class = element_class
        self.module_generators = [self(generator)]
        self._list = ClassicalCrystal.list(self)
        self._digraph = ClassicalCrystal.digraph(self)
        self._digraph_closure = self.digraph().transitive_closure()

    def __call__(self, value):
        if value.__class__ == self.element_class and value.parent == self:
            return self
        else: # Should do sanity checks!
            return self.element_class(self, value)

    def list(self):
        return self._list

    def digraph(self):
        return self._digraph

    def cmp_elements(self, x,y):
        r"""
        The crystal being classical, it is an acyclic digraph
        which can be interpreted as a poset. This implements the
        comparison function of this poset.

        Returns whether there is a path from x to y in the crystal graph
        """
        assert x.parent() == self and y.parent() == self;
        if   self._digraph_closure.has_edge(x,y):
            return -1;
        elif self._digraph_closure.has_edge(y,x):
            return 1;
        else:
            return 0;

class Spin(Element):
    r"""
    A class for spins. Derived from Letter
    """

    def __init__(self, parent, value):
#        Element.__init__(self, parent);
        self._parent = parent
        self.value = value

    def parent(self):
        return self._parent  # Should be inherited from Element!

    def __repr__(self):
        sword = ""
        m = self.value-1
        for i in range(self._parent.cartanType.n,0,-1):
            if m >= 2**(i-1):
                sword = sword + "-"
                m = m - 2**(i-1)
            else:
                sword = sword + "+"
        return(sword)

    def signature(self):
        spinlist = []
        m = self.value-1
        for i in range(self._parent.cartanType.n,0,-1):
            if m >= 2**(i-1):
                spinlist.append(-1)
                m = m - 2**(i-1)
            else:
                spinlist.append(1)
        return(spinlist)


    def internal_repr(self):
        return "%s"%self.value

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
               self.parent()  == other.parent()   and \
               self.value     == other.value

    def __cmp__(self, other):
        assert self.parent() == other.parent()
        return self.parent().cmp_elements(self, other)

class Spin_crystal_type_B_element(Spin, CrystalElement):
    r"""
    Type B spin representation crystal element
    """
    def e(self, i):
        r"""
        TEST:
            sage: C = CrystalOfSpins(['B',3])
            sage: C(1).e(1) == None
            True
            sage: C(1).e(2) == None
            True
            sage: C(1).e(3) == None
            True
            sage: C(2).e(1) == None
            True
            sage: C(2).e(2) == None
            True
            sage: C(2).e(3) == C(1)
            True
            sage: C(3).e(1) == None
            True
            sage: C(3).e(2) == C(2)
            True
            sage: C(3).e(3) == None
            True
            sage: C(4).e(1) == None
            True
            sage: C(4).e(2) == None
            True
            sage: C(4).e(3) == C(3)
            True
            sage: C(5).e(1) == C(3)
            True
            sage: C(5).e(2) == None
            True
            sage: C(5).e(3) == None
            True
            sage: C(6).e(1) == C(4)
            True
            sage: C(6).e(2) == None
            True
            sage: C(6).e(3) == C(5)
            True
            sage: C(7).e(1) == None
            True
            sage: C(7).e(2) == C(6)
            True
            sage: C(7).e(3) == None
            True
            sage: C(8).e(1) == None
            True
            sage: C(8).e(2) == None
            True
            sage: C(8).e(3) == C(7)
            True
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if int((self.value-1)/2**(rank-i-1))%4 == 2:
                return self._parent(self.value-2**(rank-i-1))
        elif i == rank:
            if (self.value-1)%2 == 1:
                return self._parent(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        TESTS:
            sage: C = CrystalOfSpins(['B',3])
            sage: C(1).f(1) == None
            True
            sage: C(1).f(2) == None
            True
            sage: C(1).f(3) == C(2)
            True
            sage: C(2).f(1) == None
            True
            sage: C(2).f(2) == C(3)
            True
            sage: C(2).f(3) == None
            True
            sage: C(3).f(1) == C(5)
            True
            sage: C(3).f(2) == None
            True
            sage: C(3).f(3) == C(4)
            True
            sage: C(4).f(1) == C(6)
            True
            sage: C(4).f(2) == None
            True
            sage: C(4).f(3) == None
            True
            sage: C(5).f(1) == None
            True
            sage: C(5).f(2) == None
            True
            sage: C(5).f(3) == C(6)
            True
            sage: C(6).f(1) == None
            True
            sage: C(6).f(2) == C(7)
            True
            sage: C(6).f(3) == None
            True
            sage: C(7).f(1) == None
            True
            sage: C(7).f(2) == None
            True
            sage: C(7).f(3) == C(8)
            True
            sage: C(8).f(1) == None
            True
            sage: C(8).f(2) == None
            True
            sage: C(8).f(3) == None
            True
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if int((self.value-1)/2**(rank-i-1))%4 == 1:
                return self._parent(self.value+2**(rank-i-1))
        elif i == rank:
            if (self.value-1)%2 == 0:
                return self._parent(self.value+1)
        else:
            return None

class Spin_crystal_type_D_element(Spin, CrystalElement):
    r"""
    Type D spin representation crystal element
    """
    def e(self, i):
        r"""
        TEST:
            sage: D = CrystalOfSpinsPlus(['D',4])
            sage: D(1).e(1) == None
            True
            sage: D(1).e(2) == None
            True
            sage: D(1).e(3) == None
            True
            sage: D(1).e(4) == None
            True
            sage: D(4).e(1) == None
            True
            sage: D(4).e(2) == None
            True
            sage: D(4).e(3) == None
            True
            sage: D(4).e(4) == D(1)
            True
            sage: D(6).e(1) == None
            True
            sage: D(6).e(2) == D(4)
            True
            sage: D(6).e(3) == None
            True
            sage: D(6).e(4) == None
            True
            sage: D(7).e(1) == None
            True
            sage: D(7).e(2) == None
            True
            sage: D(7).e(3) == D(6)
            True
            sage: D(7).e(4) == None
            True
            sage: D(10).e(1) == D(6)
            True
            sage: D(10).e(2) == None
            True
            sage: D(10).e(3) == None
            True
            sage: D(10).e(4) == None
            True
            sage: D(11).e(1) == D(7)
            True
            sage: D(11).e(2) == None
            True
            sage: D(11).e(3) == D(10)
            True
            sage: D(11).e(4) == None
            True
            sage: D(13).e(1) == None
            True
            sage: D(13).e(2) == D(11)
            True
            sage: D(13).e(3) == None
            True
            sage: D(13).e(4) == None
            True
            sage: D(16).e(1) == None
            True
            sage: D(16).e(2) == None
            True
            sage: D(16).e(3) == None
            True
            sage: D(16).e(4) == D(13)
            True
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if int((self.value-1)/2**(rank-i-1))%4 == 2:
                return self._parent(self.value-2**(rank-i-1))
        elif i == rank:
            if (self.value-1)%4 == 3:
                return self._parent(self.value-3)
        else:
            return None

    def f(self, i):
        r"""
        TESTS:
            sage: D = CrystalOfSpinsPlus(['D',4])
            sage: D(1).f(1) == None
            True
            sage: D(1).f(2) == None
            True
            sage: D(1).f(3) == None
            True
            sage: D(1).f(4) == D(4)
            True
            sage: D(4).f(1) == None
            True
            sage: D(4).f(2) == D(6)
            True
            sage: D(4).f(3) == None
            True
            sage: D(4).f(4) == None
            True
            sage: D(6).f(1) == D(10)
            True
            sage: D(6).f(2) == None
            True
            sage: D(6).f(3) == D(7)
            True
            sage: D(6).f(4) == None
            True
            sage: D(7).f(1) == D(11)
            True
            sage: D(7).f(2) == None
            True
            sage: D(7).f(3) == None
            True
            sage: D(7).f(4) == None
            True
            sage: D(10).f(1) == None
            True
            sage: D(10).f(2) == None
            True
            sage: D(10).f(3) == D(11)
            True
            sage: D(10).f(4) == None
            True
            sage: D(11).f(1) == None
            True
            sage: D(11).f(2) == D(13)
            True
            sage: D(11).f(3) == None
            True
            sage: D(11).f(4) == None
            True
            sage: D(13).f(1) == None
            True
            sage: D(13).f(2) == None
            True
            sage: D(13).f(3) == None
            True
            sage: D(13).f(4) == D(16)
            True
            sage: D(16).f(1) == None
            True
            sage: D(16).f(2) == None
            True
            sage: D(16).f(3) == None
            True
            sage: D(16).f(4) == None
            True
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if int((self.value-1)/2**(rank-i-1))%4 == 1:
                return self._parent(self.value+2**(rank-i-1))
        elif i == rank:
            if (self.value-1)%4 == 0:
                return self._parent(self.value+3)
        else:
            return None

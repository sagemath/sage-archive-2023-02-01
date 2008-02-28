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
    [[1, 1, 1],
    [1, 1, -1],
    [1, -1, 1],
    [-1, 1, 1],
    [1, -1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
    [-1, -1, -1]]

    sage: [x.signature() for x in C]
    ['+++', '++-', '+-+', '-++', '+--', '-+-', '--+', '---']

    sage: len(TensorProductOfCrystals(C,C,generators=[[C.list()[0],C.list()[0]]]))
    35
    """
    type = CartanType(type)
    if type[0] == 'B':
        return GenericCrystalOfSpins(type, Spin_crystal_type_B_element, "spins")
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
      [[1, 1, 1, 1],
       [1, 1, -1, -1],
       [1, -1, 1, -1],
       [-1, 1, 1, -1],
       [1, -1, -1, 1],
       [-1, 1, -1, 1],
       [-1, -1, 1, 1],
       [-1, -1, -1, -1]]

      sage: [x.signature() for x in D]
      ['++++', '++--', '+-+-', '-++-', '+--+', '-+-+', '--++', '----']

      sage: D.check()
      True
    """
    type = CartanType(type)
    if type[0] == 'D':
        generator = [1 for i in range(type[1])]
        return GenericCrystalOfSpins(type, Spin_crystal_type_D_element, "plus")
    else:
        raise NotImplementedError

def CrystalOfSpinsMinus(type):
    r"""
    Return the minus spin crystal of the given type D.

    INPUT:
        ['D',n] -- A CartanType of type D.

    This is the crystal with highest weight $Lambda_{n-1}$
    (the (n-1)-st fundamental weight).

    TESTS:

        sage: E = CrystalOfSpinsMinus(['D',4])
        sage: E.list()
         [[1, 1, 1, -1],
         [1, 1, -1, 1],
         [1, -1, 1, 1],
         [-1, 1, 1, 1],
         [1, -1, -1, -1],
         [-1, 1, -1, -1],
         [-1, -1, 1, -1],
         [-1, -1, -1, 1]]
        sage: [x.signature() for x in E]
        ['+++-', '++-+', '+-++', '-+++', '+---', '-+--', '--+-', '---+']
        sage: len(TensorProductOfCrystals(E,E,generators=[[E.list()[0],E.list()[0]]]).list())
        35
        sage: D = CrystalOfSpinsPlus(['D',4])
        sage: len(TensorProductOfCrystals(D,E,generators=[[D.list()[0],E.list()[0]]]).list())
        56
    """
    type = CartanType(type)
    if type[0] == 'D':
        return GenericCrystalOfSpins(type, Spin_crystal_type_D_element, "minus")
    else:
        raise NotImplementedError

class GenericCrystalOfSpins(ClassicalCrystal):
    r"""
    Derived from ClassicalCrystalOfLetters
    """
    def __init__(self, type, element_class, case):
        self.cartanType = CartanType(type)
        if case == "spins":
            self._name = "The crystal of spins for type %s"%type
        elif case == "plus":
            self._name = "The plus crystal of spins for type %s"%type
        else:
            self._name = "The minus crystal of spins for type %s"%type

        self.index_set = self.cartanType.index_set()
        self.element_class = element_class
        if case == "minus":
            generator = [1 for i in range(type[1]-1)]
            generator.append(-1)
        else:
            generator = [1 for i in range(type[1])]
        self.module_generators = [self(generator)]
        self._list = ClassicalCrystal.list(self)
        self._digraph = ClassicalCrystal.digraph(self)
        self._digraph_closure = self.digraph().transitive_closure()

    def __call__(self, value):
        if value.__class__ == self.element_class and value.parent == self:
            return self
        else:
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
        self._parent = parent
        self.value = value

    def parent(self):
        return self._parent  # Should be inherited from Element!

    def __repr__(self):
        return str(self.value)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
               self.parent()  == other.parent()   and \
               self.value     == other.value

    def __cmp__(self, other):
        assert self.parent() == other.parent()
        return self.parent().cmp_elements(self, other)

    def signature(self):
        sword = ""
        for x in range(self._parent.cartanType.n):
            if self.value[x] == 1:
                sword = sword + "+"
            else:
                sword = sword + "-"
        return(sword)

class Spin_crystal_type_B_element(Spin, CrystalElement):
    r"""
    Type B spin representation crystal element

    TESTS:
        sage: C = CrystalOfSpins(['B',3])
        sage: [[C.list()[m].e(i) for i in range(1,4)] for m in range(8)]
        [[None, None, None],
        [None, None, [1, 1, 1]],
        [None, [1, 1, -1], None],
        [[1, -1, 1], None, None],
        [None, None, [1, -1, 1]],
        [[1, -1, -1], None, [-1, 1, 1]],
        [None, [-1, 1, -1], None],
        [None, None, [-1, -1, 1]]]

        sage: [[C.list()[m].f(i) for i in range(1,4)] for m in range(8)]
        [[None, None, [1, 1, -1]],
        [None, [1, -1, 1], None],
        [[-1, 1, 1], None, [1, -1, -1]],
        [None, None, [-1, 1, -1]],
        [[-1, 1, -1], None, None],
        [None, [-1, -1, 1], None],
        [None, None, [-1, -1, -1]],
        [None, None, None]]
    """
    def e(self, i):
        r"""
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if self.value[i-1] == -1 and self.value[i] == 1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-1] = 1
                ret[i] = -1
                return self._parent(ret)
        elif i == rank:
            if self.value[i-1] == -1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-1] = 1
                return self._parent(ret)
        else:
            return None

    def f(self, i):
        r"""
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if self.value[i-1] == 1 and self.value[i] == -1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-1] = -1
                ret[i] = 1
                return self._parent(ret)
        elif i == rank:
            if self.value[i-1] == 1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-1] = -1
                return self._parent(ret)
        else:
            return None

class Spin_crystal_type_D_element(Spin, CrystalElement):
    r"""
    Type D spin representation crystal element

    TESTS:

        sage: D = CrystalOfSpinsPlus(['D',4])
        sage: [[D.list()[m].e(i) for i in range(1,4)] for m in range(8)]
        [[None, None, None],
         [None, None, None],
         [None, [1, 1, -1, -1], None],
         [[1, -1, 1, -1], None, None],
         [None, None, [1, -1, 1, -1]],
         [[1, -1, -1, 1], None, [-1, 1, 1, -1]],
         [None, [-1, 1, -1, 1], None],
         [None, None, None]]

        sage: [[D.list()[m].f(i) for i in range(1,4)] for m in range(8)]
        [[None, None, None],
         [None, [1, -1, 1, -1], None],
         [[-1, 1, 1, -1], None, [1, -1, -1, 1]],
         [None, None, [-1, 1, -1, 1]],
         [[-1, 1, -1, 1], None, None],
         [None, [-1, -1, 1, 1], None],
         [None, None, None],
         [None, None, None]]

        sage: E = CrystalOfSpinsMinus(['D',4])
        sage: [[E.list()[m].e(i) for i in range(1,4)] for m in range(8)]
        [[None, None, None],
         [None, None, [1, 1, 1, -1]],
         [None, [1, 1, -1, 1], None],
         [[1, -1, 1, 1], None, None],
         [None, None, None],
         [[1, -1, -1, -1], None, None],
         [None, [-1, 1, -1, -1], None],
         [None, None, [-1, -1, 1, -1]]]

        sage: [[E.list()[m].f(i) for i in range(1,4)] for m in range(8)]
        [[None, None, [1, 1, -1, 1]],
         [None, [1, -1, 1, 1], None],
         [[-1, 1, 1, 1], None, None],
         [None, None, None],
         [[-1, 1, -1, -1], None, None],
         [None, [-1, -1, 1, -1], None],
         [None, None, [-1, -1, -1, 1]],
         [None, None, None]]
    """
    def e(self, i):
        r"""
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if self.value[i-1] == -1 and self.value[i] == 1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-1] = 1
                ret[i] = -1
                return self._parent(ret)
        elif i == rank:
            if self.value[i-2] == -1 and self.value[i-1] == -1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-2] = 1
                ret[i-1] = 1
                return self._parent(ret)
        else:
            return None

    def f(self, i):
        r"""
        TESTS:
        """
        assert i in self.index_set()
        rank = self._parent.cartanType.n
        if i < rank:
            if self.value[i-1] == 1 and self.value[i] == -1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-1] = -1
                ret[i] = 1
                return self._parent(ret)
        elif i == rank:
            if self.value[i-2] == 1 and self.value[i-1] == 1:
                ret = [self.value[x] for x in range(rank)]
                ret[i-2] = -1
                ret[i-1] = -1
                return self._parent(ret)
        else:
            return None

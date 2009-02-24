r"""
Crystals of letters
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

from sage.structure.element import Element
from sage.combinat.root_system.cartan_type import CartanType
from crystals import ClassicalCrystal, CrystalElement

def CrystalOfLetters(cartan_type):
    r"""
    Returns the crystal of letters of the given type.

    For classical types, this is a combinatorial model for the crystal
    with highest weight Lambda_1 (the first fundamental weight).

    Any irreducible classical crystal appears as the irreducible
    component of the tensor product of several copies of this crystal
    (plus possibly one copy of the spin crystal, see CrystalOfSpins).
    See M. Kashiwara, T. Nakashima,
    *Crystal graphs for representations of   the `q`-analogue of classical Lie algebras*,
    J. Algebra **165** (1994), no. 2, 295-345. Elements of this
    irreducible component have a fixed shape, and can be fit inside a
    tableau shape. Otherwise said, any irreducible classical crystal is
    isomorphic to a crystal of tableaux with cells filled by elements
    of the crystal of letters (possibly tensored with the crystal of
    spins).

    INPUT:


    -  ``T`` - A CartanType


    EXAMPLES::

        sage: C = CrystalOfLetters(['A',5])
        sage: C.list()
        [1, 2, 3, 4, 5, 6]
    """
    ct = CartanType(cartan_type)
    if ct[0] == 'A':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_A_element)
    elif ct[0] == 'B':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_B_element)
    elif ct[0] == 'C':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_C_element)
    elif ct[0] == 'D':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_D_element)
    elif ct[0] == 'G':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_G_element)
    else:
        raise NotImplementedError


class ClassicalCrystalOfLetters(ClassicalCrystal):
    r"""
    A generic class for classical crystals of letters.

    All classical crystals of letters should be instances of this class
    or of subclasses. To define a new crystal of letters, one only
    needs to implement a class for the elements (which subclasses
    Letter and CrystalElement), with appropriate e and f operations. If
    the module generator is not 1, one also needs to define the
    subclass ClassicalCrystalOfLetters for the crystal itself.

    The basic assumption is that crystals of letters are small, but
    used intensively as building blocks. Therefore, we explicitly build
    in memory the list of all elements, the crystal graph and its
    transitive closure, so as to make the following operations constant
    time: list, cmp, (todo: phi, epsilon, e, f with caching)
    """
    def __init__(self, cartan_type, element_class):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: C == loads(dumps(C))
            True
        """
        self.cartan_type = CartanType(cartan_type)
        self._name = "The crystal of letters for type %s"%cartan_type
        self.index_set = self.cartan_type.index_set()
        self.element_class = element_class
        self.module_generators = [self(1)]
        self._list = ClassicalCrystal.list(self)
        self._digraph = ClassicalCrystal.digraph(self)
        self._digraph_closure = self.digraph().transitive_closure()

    def __call__(self, value):
        """
        Coerces value into self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: c = C(1); c
            1
            sage: c.parent()
            The crystal of letters for type ['A', 5]
            sage: c is C(c)
            True
        """
        if value.__class__ == self.element_class and value.parent() == self:
            return value
        else: # Should do sanity checks!
            return self.element_class(self, value)


    def list(self):
        """
        Returns a list of the elements of self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: C.list()
            [1, 2, 3, 4, 5, 6]
        """
        return self._list


    def digraph(self):
        """
        Returns the directed graph associated to self.

        EXAMPLES::

            sage: CrystalOfLetters(['A',5]).digraph()
            Digraph on 6 vertices
        """
        return self._digraph

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: 1 in C
            False
            sage: C(1) in C
            True
        """
        return x in self._list

    def cmp_elements(self, x, y):
        r"""
        Returns True if and only if there is a path from x to y in the
        crystal graph.

        Because the crystal graph is classical, it is a directed acyclic
        graph which can be interpreted as a poset. This function implements
        the comparison function of this poset.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A', 5])
            sage: x = C(1)
            sage: y = C(2)
            sage: C.cmp_elements(x,y)
            -1
            sage: C.cmp_elements(y,x)
            1
            sage: C.cmp_elements(x,x)
            0
        """
        assert x.parent() == self and y.parent() == self
        if self._digraph_closure.has_edge(x,y):
            return -1
        elif self._digraph_closure.has_edge(y,x):
            return 1
        else:
            return 0

    # TODO: cmp, in, ...

# Utility. Note: much of this class should be factored out at some point!
class Letter(Element):
    r"""
    A class for letters
    """

    def __init__(self, parent, value):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.letters import Letter
            sage: a = Letter(ZZ, 1)
            sage: a == loads(dumps(a))
            True
        """
        self._parent = parent
        self.value = value

    def parent(self):
        """
        Returns the parent of self.

        EXAMPLES::

            sage: from sage.combinat.crystals.letters import Letter
            sage: Letter(ZZ, 1).parent()
            Integer Ring
        """
        return self._parent  # Should be inherited from Element!

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.letters import Letter
            sage: Letter(ZZ, 1).__repr__()
            '1'
        """
        return "%s"%self.value

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.letters import Letter
            sage: parent1 = 1  # Any fake value ...
            sage: parent2 = 2  # Any fake value ...
            sage: l11 = Letter(parent1, 1)
            sage: l12 = Letter(parent1, 2)
            sage: l21 = Letter(parent2, 1)
            sage: l22 = Letter(parent2, 2)
            sage: l11 == l11
            True
            sage: l11 == l12
            False
            sage: l11 == l21
            False
        """
        return self.__class__ is other.__class__ and \
               self.parent() == other.parent() and \
               self.value == other.value

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A', 5])
            sage: C(1) < C(2)
            True
            sage: C(2) < C(1)
            False
            sage: C(2) > C(1)
            True
            sage: C(1) <= C(1)
            True
        """
        if type(self) is not type(other):
            return cmp(type(self), type(other))
        if self.parent() != other.parent():
            return cmp(self.parent(), other.parent())
        return self.parent().cmp_elements(self, other)

#########################
# Type A
#########################

class Crystal_of_letters_type_A_element(Letter, CrystalElement):
    r"""
    Type A crystal of letters elements

    TESTS::

        sage: C = CrystalOfLetters (['A',3])
        sage: C.list()
        [1, 2, 3, 4]
        sage: [ [x < y for y in C] for x in C ]
        [[False, True, True, True],
         [False, False, True, True],
         [False, False, False, True],
         [False, False, False, False]]

    ::

        sage: C = CrystalOfLetters(['A',5])
        sage: C(1) < C(1), C(1) < C(2), C(1) < C(3), C(2) < C(1)
        (False, True, True, False)

    ::

        sage: C.check()
        True
    """

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: [v.weight() for v in CrystalOfLetters(['A',3])]
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        """
        return self._parent.weight_lattice_realization()._term(self.value-1)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',4])
            sage: [[c,i,c.e(i)] for i in C.index_set for c in C if c.e(i) is not None]
            [[2, 1, 1], [3, 2, 2], [4, 3, 3], [5, 4, 4]]
        """
        assert i in self.index_set()
        if self.value == i+1:
            return self._parent(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',4])
            sage: [[c,i,c.f(i)] for i in C.index_set for c in C if c.f(i) is not None]
            [[1, 1, 2], [2, 2, 3], [3, 3, 4], [4, 4, 5]]
        """
        assert i in self.index_set()
        if self.value == i:
            return self._parent(self.value+1)
        else:
            return None

#########################
# Type B
#########################

class Crystal_of_letters_type_B_element(Letter, CrystalElement):
    r"""
    Type B crystal of letters elements

    TESTS::

        sage: C = CrystalOfLetters (['B',3])
        sage: C.check()
        True
    """

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: [v.weight() for v in CrystalOfLetters(['B',3])]
            [(1, 0, 0),
             (0, 1, 0),
             (0, 0, 1),
             (0, 0, 0),
             (0, 0, -1),
             (0, -1, 0),
             (-1, 0, 0)]
        """
        if self.value > 0:
            return self._parent.weight_lattice_realization()._term(self.value-1)
        elif self.value < 0:
            return -self._parent.weight_lattice_realization()._term(-self.value-1)
        else:
            return self._parent.weight_lattice_realization()(0)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['B',4])
            sage: [[c,i,c.e(i)] for i in C.index_set for c in C if c.e(i) is not None]
            [[2, 1, 1],
             [-1, 1, -2],
             [3, 2, 2],
             [-2, 2, -3],
             [4, 3, 3],
             [-3, 3, -4],
             [0, 4, 4],
             [-4, 4, 0]]
        """
        assert i in self.index_set()
        if self.value == i+1:
            return self._parent(i)
        elif self.value == 0 and i == self._parent.cartan_type.n:
            return self._parent(self._parent.cartan_type.n)
        elif self.value == -i:
            if i == self._parent.cartan_type.n:
                return self._parent(0)
            else:
                return self._parent(-i-1)
        else:
            return None

    def f(self, i):
        r"""
        Returns the actions of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['B',4])
            sage: [[c,i,c.f(i)] for i in C.index_set for c in C if c.f(i) is not None]
            [[1, 1, 2],
             [-2, 1, -1],
             [2, 2, 3],
             [-3, 2, -2],
             [3, 3, 4],
             [-4, 3, -3],
             [4, 4, 0],
             [0, 4, -4]]
        """
        assert i in self.index_set()
        if self.value == i:
            if i < self._parent.cartan_type.n:
                return self._parent(i+1)
            else:
                return self._parent(0)
        elif self.value == 0 and i == self._parent.cartan_type.n:
            return self._parent(-self._parent.cartan_type.n)
        elif self.value == -i-1:
            return(self._parent(-i))
        else:
            return None

#########################
# Type C
#########################

class Crystal_of_letters_type_C_element(Letter, CrystalElement):
    r"""
    Type C crystal of letters elements

    TESTS::

        sage: C = CrystalOfLetters (['C',3])
        sage: C.list()
        [1, 2, 3, -3, -2, -1]
        sage: [ [x < y for y in C] for x in C ]
        [[False, True, True, True, True, True],
         [False, False, True, True, True, True],
         [False, False, False, True, True, True],
         [False, False, False, False, True, True],
         [False, False, False, False, False, True],
         [False, False, False, False, False, False]]
        sage: C.check()
        True
    """

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: [v.weight() for v in CrystalOfLetters(['C',3])]
            [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1), (0, -1, 0), (-1, 0, 0)]
        """
        if self.value > 0:
            return self._parent.weight_lattice_realization()._term(self.value-1)
        elif self.value < 0:
            return -self._parent.weight_lattice_realization()._term(-self.value-1)
        else:
            return self._parent.weight_lattice_realization()(0)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['C',4])
            sage: [[c,i,c.e(i)] for i in C.index_set for c in C if c.e(i) is not None]
            [[2, 1, 1],
             [-1, 1, -2],
             [3, 2, 2],
             [-2, 2, -3],
             [4, 3, 3],
             [-3, 3, -4],
             [-4, 4, 4]]
        """
        assert i in self.index_set()
        if self.value == -self._parent.cartan_type.n and self.value == -i:
            return self._parent(-self.value)
        elif self.value == i+1 or self.value == -i:
            return self._parent(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        Retursn the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['C',4])
            sage: [[c,i,c.f(i)] for i in C.index_set for c in C if c.f(i) is not None]
            [[1, 1, 2],
             [-2, 1, -1],
             [2, 2, 3],
             [-3, 2, -2],
             [3, 3, 4],
             [-4, 3, -3],
             [4, 4, -4]]
        """
        assert i in self.index_set()
        if self.value == self._parent.cartan_type.n and self.value == i:
            return  self._parent(-self.value)
        elif self.value == i or self.value == -i-1:
            return self._parent(self.value+1)
        else:
            return None

#########################
# Type D
#########################

class Crystal_of_letters_type_D_element(Letter, CrystalElement):
    r"""
    Type D crystal of letters elements

    TESTS::

        sage: C = CrystalOfLetters(['D',4])
        sage: C.list()
        [1, 2, 3, 4, -4, -3, -2, -1]
        sage: C.check()
        True
    """

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: [v.weight() for v in CrystalOfLetters(['D',4])]
            [(1, 0, 0, 0),
             (0, 1, 0, 0),
             (0, 0, 1, 0),
             (0, 0, 0, 1),
             (0, 0, 0, -1),
             (0, 0, -1, 0),
             (0, -1, 0, 0),
             (-1, 0, 0, 0)]
        """
        if self.value > 0:
            return self._parent.weight_lattice_realization()._term(self.value-1)
        elif self.value < 0:
            return -self._parent.weight_lattice_realization()._term(-self.value-1)
        else:
            return self._parent.weight_lattice_realization()(0)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['D',5])
            sage: [[c,i,c.e(i)] for i in C.index_set for c in C if c.e(i) is not None]
            [[2, 1, 1],
             [-1, 1, -2],
             [3, 2, 2],
             [-2, 2, -3],
             [4, 3, 3],
             [-3, 3, -4],
             [5, 4, 4],
             [-4, 4, -5],
             [-5, 5, 4],
             [-4, 5, 5]]
        """
        assert i in self.index_set()
        if i == self._parent.cartan_type.n:
            if self.value == -i:
                return self._parent(i-1)
            elif self.value == -(i-1):
                return self._parent(i)
            else:
                return None
        elif self.value == i+1:
            return self._parent(i)
        elif self.value == -i:
            return self._parent(-(i+1))
        else:
            return None

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['D',5])
            sage: [[c,i,c.f(i)] for i in C.index_set for c in C if c.f(i) is not None]
            [[1, 1, 2],
             [-2, 1, -1],
             [2, 2, 3],
             [-3, 2, -2],
             [3, 3, 4],
             [-4, 3, -3],
             [4, 4, 5],
             [-5, 4, -4],
             [4, 5, -5],
             [5, 5, -4]]
        """
        assert i in self.index_set()
        if i == self.value:
            if i == self._parent.cartan_type.n:
                return self._parent(-(i-1))
            else:
                return self._parent(i+1)
        elif self.value == -(i+1):
            return self._parent(-i)
        elif self.value == self._parent.cartan_type.n-1 and i == self.value+1:
            return self._parent(-i)
        else:
            return None

#########################
# Type G2
#########################

class Crystal_of_letters_type_G_element(Letter, CrystalElement):
    r"""
    Type G2 crystal of letters elements

    TESTS::

        sage: C = CrystalOfLetters(['G',2])
        sage: C.list()
        [1, 2, 3, 0, -3, -2, -1]
        sage: C.check()
        True
    """

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: [v.weight() for v in CrystalOfLetters(['G',2])]
            [(1, 0, -1), (1, -1, 0), (0, 1, -1), (0, 0, 0), (0, -1, 1), (-1, 1, 0), (-1, 0, 1)]
        """
        if self.value == 1:
            return self._parent.weight_lattice_realization()((1, 0, -1))
        elif self.value == 2:
            return self._parent.weight_lattice_realization()((1, -1, 0))
        elif self.value == 3:
            return self._parent.weight_lattice_realization()((0, 1, -1))
        elif self.value == 0:
            return self._parent.weight_lattice_realization()((0, 0, 0))
        elif self.value == -3:
            return self._parent.weight_lattice_realization()((0, -1, 1))
        elif self.value == -2:
            return self._parent.weight_lattice_realization()((-1, 1, 0))
        elif self.value == -1:
            return self._parent.weight_lattice_realization()((-1, 0, 1))
        else:
            raise RuntimeError, "G2 crystal of letters element %d not valid"%self.value

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['G',2])
            sage: [[c,i,c.e(i)] for i in C.index_set for c in C if c.e(i) is not None]
            [[2, 1, 1],
             [0, 1, 3],
             [-3, 1, 0],
             [-1, 1, -2],
             [3, 2, 2],
             [-2, 2, -3]]
        """
        assert i in self.index_set()
        if i == 1:
            if self.value == 2:
                return self._parent(1)
            elif self.value == 0:
                return self._parent(3)
            elif self.value == -3:
                return self._parent(0)
            elif self.value == -1:
                return self._parent(-2)
            else:
                return None
        else:
            if self.value == 3:
                return self._parent(2)
            elif self.value == -2:
                return self._parent(-3)
            else:
                return None

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['G',2])
            sage: [[c,i,c.f(i)] for i in C.index_set for c in C if c.f(i) is not None]
            [[1, 1, 2],
             [3, 1, 0],
             [0, 1, -3],
             [-2, 1, -1],
             [2, 2, 3],
             [-3, 2, -2]]
        """
        assert i in self.index_set()
        if i == 1:
            if self.value == 1:
                return self._parent(2)
            elif self.value == 3:
                return self._parent(0)
            elif self.value == 0:
                return self._parent(-3)
            elif self.value == -2:
                return self._parent(-1)
            else:
                return None
        else:
            if self.value == 2:
                return self._parent(3)
            elif self.value == -3:
                return self._parent(-2)
            else:
                return None



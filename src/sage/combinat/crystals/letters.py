r"""
Crystals of letters
"""

#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
#                          Daniel Bump    <bump at match.stanford.edu>
#                          Brant Jones    <brant at math.ucdavis.edu>
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

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.combinat.root_system.cartan_type import CartanType
from crystals import ClassicalCrystal, CrystalElement

def CrystalOfLetters(cartan_type, element_print_style = None, dual = None):
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
	sage: C.cartan_type()
	['A', 5]

    For type E6, one can also specify how elements are printed.
    This option is usually set to None and the default representation is used.
    If one chooses the option 'compact', the elements are printed in the more
    compact convention with 27 letters +abcdefghijklmnopqrstuvwxyz and
    the 27 letters -ABCDEFGHIJKLMNOPQRSTUVWXYZ for the dual crystal.

    EXAMPLES::

        sage: C = CrystalOfLetters(['E',6], element_print_style = 'compact')
	sage: C.list()
	[+, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z]
        sage: C = CrystalOfLetters(['E',6], element_print_style = 'compact', dual = True)
	sage: C.list()
	[-, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z]
    """
    ct = CartanType(cartan_type)
    if ct.letter == 'A':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_A_element)
    elif ct.letter == 'B':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_B_element)
    elif ct.letter == 'C':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_C_element)
    elif ct.letter == 'D':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_D_element)
    elif ct.letter == 'E' and ct[1] == 6:
	if dual is None:
	    return ClassicalCrystalOfLetters(ct,
					     Crystal_of_letters_type_E6_element, element_print_style)
	else:
	    return ClassicalCrystalOfLetters(ct,
					     Crystal_of_letters_type_E6_element_dual, element_print_style, dual = True)
    elif ct.letter == 'G':
        return ClassicalCrystalOfLetters(ct,
                                         Crystal_of_letters_type_G_element)
    else:
        raise NotImplementedError


class ClassicalCrystalOfLetters(UniqueRepresentation, ClassicalCrystal):
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
    def __init__(self, cartan_type, element_class, element_print_style = None, dual = None):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: C == loads(dumps(C))
            True
        """
        self._cartan_type = CartanType(cartan_type)
        self._name = "The crystal of letters for type %s"%cartan_type
        self.element_class = element_class
	if cartan_type == CartanType(['E',6]):
            if dual:
                self.module_generators = [self([6])]
                self._ambient = CrystalOfLetters(CartanType(['E',6]))
            else:
                self.module_generators = [self([1])]
       	else:
	    self.module_generators = [self(1)]
        self._list = ClassicalCrystal.list(self)
	self._element_print_style = element_print_style
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

    def lt_elements(self, x, y):
        r"""
        Returns True if and only if there is a path from x to y in the
        crystal graph, when x is not equal to y.

        Because the crystal graph is classical, it is a directed acyclic
        graph which can be interpreted as a poset. This function implements
        the comparison function of this poset.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A', 5])
            sage: x = C(1)
            sage: y = C(2)
            sage: C.lt_elements(x,y)
            True
            sage: C.lt_elements(y,x)
            False
            sage: C.lt_elements(x,x)
            False
	    sage: C = CrystalOfLetters(['D', 4])
	    sage: C.lt_elements(C(4),C(-4))
	    False
	    sage: C.lt_elements(C(-4),C(4))
	    False
        """
        assert x.parent() == self and y.parent() == self
        if self._digraph_closure.has_edge(x,y):
            return True
	return False

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

    def __ne__(self, other):
	"""
	EXAMPLES::

	    sage: C = CrystalOfLetters(['B', 3])
	    sage: C(0) <> C(0)
	    False
	    sage: C(1) <> C(-1)
	    True
        """
	return not self == other

    def __lt__(self, other):
	"""
	EXAMPLES::

	    sage: C = CrystalOfLetters(['D', 4])
	    sage: C(-4) < C(4)
	    False
	    sage: C(4) < C(-3)
	    True
	    sage: C(4) < C(4)
	    False
        """
	return self.parent().lt_elements(self, other)

    def __gt__(self, other):
	"""
	EXAMPLES::

	    sage: C = CrystalOfLetters(['D', 4])
	    sage: C(-4) > C(4)
	    False
	    sage: C(4) > C(-3)
	    False
	    sage: C(4) < C(4)
	    False
	    sage: C(-1) > C(1)
	    True
        """
	return other.__lt__(self)

    def __le__(self, other):
	"""
	EXAMPLES::

	    sage: C = CrystalOfLetters(['D', 4])
	    sage: C(-4) <= C(4)
	    False
	    sage: C(4) <= C(-3)
	    True
	    sage: C(4) <= C(4)
	    True
        """
	return self.__lt__(other) or self == other

    def __ge__(self, other):
	"""
	EXAMPLES::

	    sage: C = CrystalOfLetters(['D', 4])
	    sage: C(-4) >= C(4)
	    False
	    sage: C(4) >= C(-3)
	    False
	    sage: C(4) >= C(4)
	    True
        """
	return other.__le__(self)

#    def __cmp__(self, other):
#        """
#        EXAMPLES::
#
#            sage: C = CrystalOfLetters(['A', 5])
#            sage: C(1) < C(2)
#            True
#            sage: C(2) < C(1)
#            False
#            sage: C(2) > C(1)
#            True
#            sage: C(1) <= C(1)
#            True
#        """
#        if type(self) is not type(other):
#            return cmp(type(self), type(other))
#        if self.parent() != other.parent():
#            return cmp(self.parent(), other.parent())
#        return self.parent().cmp_elements(self, other)

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
        return self.parent().weight_lattice_realization()._term(self.value-1)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',4])
            sage: [[c,i,c.e(i)] for i in C.index_set() for c in C if c.e(i) is not None]
            [[2, 1, 1], [3, 2, 2], [4, 3, 3], [5, 4, 4]]
        """
        assert i in self.index_set()
        if self.value == i+1:
            return self.parent()(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',4])
            sage: [[c,i,c.f(i)] for i in C.index_set() for c in C if c.f(i) is not None]
            [[1, 1, 2], [2, 2, 3], [3, 3, 4], [4, 4, 5]]
        """
        assert i in self.index_set()
        if self.value == i:
            return self.parent()(self.value+1)
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
            return self.parent().weight_lattice_realization()._term(self.value-1)
        elif self.value < 0:
            return -self.parent().weight_lattice_realization()._term(-self.value-1)
        else:
            return self.parent().weight_lattice_realization()(0)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['B',4])
            sage: [[c,i,c.e(i)] for i in C.index_set() for c in C if c.e(i) is not None]
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
            return self.parent()(i)
        elif self.value == 0 and i == self.parent().cartan_type().n:
            return self.parent()(self.parent().cartan_type().n)
        elif self.value == -i:
            if i == self.parent().cartan_type().n:
                return self.parent()(0)
            else:
                return self.parent()(-i-1)
        else:
            return None

    def f(self, i):
        r"""
        Returns the actions of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['B',4])
            sage: [[c,i,c.f(i)] for i in C.index_set() for c in C if c.f(i) is not None]
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
            if i < self.parent().cartan_type().n:
                return self.parent()(i+1)
            else:
                return self.parent()(0)
        elif self.value == 0 and i == self.parent().cartan_type().n:
            return self.parent()(-self.parent().cartan_type().n)
        elif self.value == -i-1:
            return(self.parent()(-i))
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
            return self.parent().weight_lattice_realization()._term(self.value-1)
        elif self.value < 0:
            return -self.parent().weight_lattice_realization()._term(-self.value-1)
        else:
            return self.parent().weight_lattice_realization()(0)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['C',4])
            sage: [[c,i,c.e(i)] for i in C.index_set() for c in C if c.e(i) is not None]
            [[2, 1, 1],
             [-1, 1, -2],
             [3, 2, 2],
             [-2, 2, -3],
             [4, 3, 3],
             [-3, 3, -4],
             [-4, 4, 4]]
        """
        assert i in self.index_set()
        if self.value == -self.parent().cartan_type().n and self.value == -i:
            return self.parent()(-self.value)
        elif self.value == i+1 or self.value == -i:
            return self.parent()(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        Retursn the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['C',4])
            sage: [[c,i,c.f(i)] for i in C.index_set() for c in C if c.f(i) is not None]
            [[1, 1, 2],
             [-2, 1, -1],
             [2, 2, 3],
             [-3, 2, -2],
             [3, 3, 4],
             [-4, 3, -3],
             [4, 4, -4]]
        """
        assert i in self.index_set()
        if self.value == self.parent().cartan_type().n and self.value == i:
            return  self.parent()(-self.value)
        elif self.value == i or self.value == -i-1:
            return self.parent()(self.value+1)
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
            return self.parent().weight_lattice_realization()._term(self.value-1)
        elif self.value < 0:
            return -self.parent().weight_lattice_realization()._term(-self.value-1)
        else:
            return self.parent().weight_lattice_realization()(0)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['D',5])
            sage: [[c,i,c.e(i)] for i in C.index_set() for c in C if c.e(i) is not None]
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
        if i == self.parent().cartan_type().n:
            if self.value == -i:
                return self.parent()(i-1)
            elif self.value == -(i-1):
                return self.parent()(i)
            else:
                return None
        elif self.value == i+1:
            return self.parent()(i)
        elif self.value == -i:
            return self.parent()(-(i+1))
        else:
            return None

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['D',5])
            sage: [[c,i,c.f(i)] for i in C.index_set() for c in C if c.f(i) is not None]
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
            if i == self.parent().cartan_type().n:
                return self.parent()(-(i-1))
            else:
                return self.parent()(i+1)
        elif self.value == -(i+1):
            return self.parent()(-i)
        elif self.value == self.parent().cartan_type().n-1 and i == self.value+1:
            return self.parent()(-i)
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
            return self.parent().weight_lattice_realization()((1, 0, -1))
        elif self.value == 2:
            return self.parent().weight_lattice_realization()((1, -1, 0))
        elif self.value == 3:
            return self.parent().weight_lattice_realization()((0, 1, -1))
        elif self.value == 0:
            return self.parent().weight_lattice_realization()((0, 0, 0))
        elif self.value == -3:
            return self.parent().weight_lattice_realization()((0, -1, 1))
        elif self.value == -2:
            return self.parent().weight_lattice_realization()((-1, 1, 0))
        elif self.value == -1:
            return self.parent().weight_lattice_realization()((-1, 0, 1))
        else:
            raise RuntimeError, "G2 crystal of letters element %d not valid"%self.value

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['G',2])
            sage: [[c,i,c.e(i)] for i in C.index_set() for c in C if c.e(i) is not None]
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
                return self.parent()(1)
            elif self.value == 0:
                return self.parent()(3)
            elif self.value == -3:
                return self.parent()(0)
            elif self.value == -1:
                return self.parent()(-2)
            else:
                return None
        else:
            if self.value == 3:
                return self.parent()(2)
            elif self.value == -2:
                return self.parent()(-3)
            else:
                return None

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['G',2])
            sage: [[c,i,c.f(i)] for i in C.index_set() for c in C if c.f(i) is not None]
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
                return self.parent()(2)
            elif self.value == 3:
                return self.parent()(0)
            elif self.value == 0:
                return self.parent()(-3)
            elif self.value == -2:
                return self.parent()(-1)
            else:
                return None
        else:
            if self.value == 2:
                return self.parent()(3)
            elif self.value == -3:
                return self.parent()(-2)
            else:
                return None


#########################
# Type E6
#########################

class Crystal_of_letters_type_E6_element(Letter, CrystalElement):
    r"""
    Type `E_6` crystal of letters elements. This crystal corresponds to the highest weight
    crystal `B(\Lambda_1)`.

    TESTS::

        sage: C = CrystalOfLetters(['E',6])
	sage: C.module_generators
	[[1]]
	sage: C.list()
	[[1], [-1, 3], [-3, 4], [-4, 2, 5], [-2, 5], [-5, 2, 6], [-2, -5, 4, 6],
	[-4, 3, 6], [-3, 1, 6], [-1, 6], [-6, 2], [-2, -6, 4], [-4, -6, 3, 5],
	[-3, -6, 1, 5], [-1, -6, 5], [-5, 3], [-3, -5, 1, 4], [-1, -5, 4], [-4, 1, 2],
	[-1, -4, 2, 3], [-3, 2], [-2, -3, 4], [-4, 5], [-5, 6], [-6], [-2, 1], [-1, -2, 3]]
        sage: C.check()
        True
	sage: all(b.f(i).e(i) == b for i in C.index_set() for b in C if b.f(i) is not None)
	True
	sage: all(b.e(i).f(i) == b for i in C.index_set() for b in C if b.e(i) is not None)
	True
	sage: G = C.digraph()
	sage: G.show(edge_labels=true, figsize=12, vertex_size=1)
    """

    def __repr__(self):
	"""
	In their full representation, the vertices of this crystal are labeled
	by their weight. For example vertex [-5,2,6] indicates that a 5-arrow
	is coming into this vertex, and a 2-arrow and 6-arrow is leaving the vertex.
	Specifying element_print_style = 'compact' for a given crystal C, labels the
	vertices of this crystal by the 27 letters +abcdefghijklmnopqrstuvwxyz.

	EXAMPLES::

            sage: C = CrystalOfLetters(['E',6], element_print_style = 'compact')
	    sage: C.list()
	    [+, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z]
        """
	if self.parent()._element_print_style == 'compact':
	    l=['+','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
	    return l[self.parent().list().index(self)]
	return "%s"%self.value

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: [v.weight() for v in CrystalOfLetters(['E',6])]
	    [(0, 0, 0, 0, 0, -2/3, -2/3, 2/3),
	    (-1/2, 1/2, 1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
	    (1/2, -1/2, 1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
	    (1/2, 1/2, -1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
	    (-1/2, -1/2, -1/2, 1/2, 1/2, -1/6, -1/6, 1/6),
	    (1/2, 1/2, 1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
	    (-1/2, -1/2, 1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
	    (-1/2, 1/2, -1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
	    (1/2, -1/2, -1/2, -1/2, 1/2, -1/6, -1/6, 1/6),
	    (0, 0, 0, 0, 1, 1/3, 1/3, -1/3),
	    (1/2, 1/2, 1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
	    (-1/2, -1/2, 1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
	    (-1/2, 1/2, -1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
	    (1/2, -1/2, -1/2, 1/2, -1/2, -1/6, -1/6, 1/6),
	    (0, 0, 0, 1, 0, 1/3, 1/3, -1/3),
	    (-1/2, 1/2, 1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
	    (1/2, -1/2, 1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
	    (0, 0, 1, 0, 0, 1/3, 1/3, -1/3),
	    (1/2, 1/2, -1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
	    (0, 1, 0, 0, 0, 1/3, 1/3, -1/3),
	    (1, 0, 0, 0, 0, 1/3, 1/3, -1/3),
	    (0, -1, 0, 0, 0, 1/3, 1/3, -1/3),
	    (0, 0, -1, 0, 0, 1/3, 1/3, -1/3),
	    (0, 0, 0, -1, 0, 1/3, 1/3, -1/3),
	    (0, 0, 0, 0, -1, 1/3, 1/3, -1/3),
	    (-1/2, -1/2, -1/2, -1/2, -1/2, -1/6, -1/6, 1/6),
	    (-1, 0, 0, 0, 0, 1/3, 1/3, -1/3)]
        """
	R=self.parent().weight_lattice_realization().fundamental_weights()
	return sum(cmp(i,0)*R[abs(i)] for i in self.value)

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

	    sage: C = CrystalOfLetters(['E',6])
	    sage: C([-1,3]).e(1)
	    [1]
	    sage: C([-2,-3,4]).e(2)
	    [-3, 2]
	    sage: C([1]).e(1)
        """
        assert i in self.index_set()

        if self.value == [-1, 3] and i == 1:
            return self.parent()([1])
        if self.value == [-3, 4] and i == 3:
            return self.parent()([-1, 3])
        if self.value == [-4, 2, 5] and i == 4:
            return self.parent()([-3, 4])
        if self.value == [-5, 2, 6] and i == 5:
            return self.parent()([-4, 2, 5])
        if self.value == [-2, 5] and i == 2:
            return self.parent()([-4, 2, 5])
        if self.value == [-6, 2] and i == 6:
            return self.parent()([-5, 2, 6])
        if self.value == [-2, -5, 4, 6] and i == 2:
            return self.parent()([-5, 2, 6])
        if self.value == [-2, -6, 4] and i == 2:
            return self.parent()([-6, 2])
        if self.value == [-2, -5, 4, 6] and i == 5:
            return self.parent()([-2, 5])
        if self.value == [-2, -6, 4] and i == 6:
            return self.parent()([-2, -5, 4, 6])
        if self.value == [-4, 3, 6] and i == 4:
            return self.parent()([-2, -5, 4, 6])
        if self.value == [-4, -6, 3, 5] and i == 4:
            return self.parent()([-2, -6, 4])
        if self.value == [-4, -6, 3, 5] and i == 6:
            return self.parent()([-4, 3, 6])
        if self.value == [-3, 1, 6] and i == 3:
            return self.parent()([-4, 3, 6])
        if self.value == [-5, 3] and i == 5:
            return self.parent()([-4, -6, 3, 5])
        if self.value == [-3, -6, 1, 5] and i == 3:
            return self.parent()([-4, -6, 3, 5])
        if self.value == [-3, -5, 1, 4] and i == 3:
            return self.parent()([-5, 3])
        if self.value == [-3, -6, 1, 5] and i == 6:
            return self.parent()([-3, 1, 6])
        if self.value == [-1, 6] and i == 1:
            return self.parent()([-3, 1, 6])
        if self.value == [-3, -5, 1, 4] and i == 5:
            return self.parent()([-3, -6, 1, 5])
        if self.value == [-1, -6, 5] and i == 1:
            return self.parent()([-3, -6, 1, 5])
        if self.value == [-4, 1, 2] and i == 4:
            return self.parent()([-3, -5, 1, 4])
        if self.value == [-1, -5, 4] and i == 1:
            return self.parent()([-3, -5, 1, 4])
        if self.value == [-2, 1] and i == 2:
            return self.parent()([-4, 1, 2])
        if self.value == [-1, -4, 2, 3] and i == 1:
            return self.parent()([-4, 1, 2])
        if self.value == [-1, -2, 3] and i == 1:
            return self.parent()([-2, 1])
        if self.value == [-1, -6, 5] and i == 6:
            return self.parent()([-1, 6])
        if self.value == [-1, -5, 4] and i == 5:
            return self.parent()([-1, -6, 5])
        if self.value == [-1, -4, 2, 3] and i == 4:
            return self.parent()([-1, -5, 4])
        if self.value == [-1, -2, 3] and i == 2:
            return self.parent()([-1, -4, 2, 3])
        if self.value == [-3, 2] and i == 3:
            return self.parent()([-1, -4, 2, 3])
        if self.value == [-2, -3, 4] and i == 3:
            return self.parent()([-1, -2, 3])
        if self.value == [-2, -3, 4] and i == 2:
            return self.parent()([-3, 2])
        if self.value == [-4, 5] and i == 4:
            return self.parent()([-2, -3, 4])
        if self.value == [-5, 6] and i == 5:
            return self.parent()([-4, 5])
        if self.value == [-6] and i == 6:
            return self.parent()([-5, 6])

        else:
            return None

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

	    sage: C = CrystalOfLetters(['E',6])
	    sage: C([1]).f(1)
	    [-1, 3]
	    sage: C([-6]).f(1)
        """
        assert i in self.index_set()

        if self.value == [1] and i == 1:
            return self.parent()([-1, 3])
        if self.value == [-1, 3] and i == 3:
            return self.parent()([-3, 4])
        if self.value == [-3, 4] and i == 4:
            return self.parent()([-4, 2, 5])
        if self.value == [-4, 2, 5] and i == 5:
            return self.parent()([-5, 2, 6])
        if self.value == [-4, 2, 5] and i == 2:
            return self.parent()([-2, 5])
        if self.value == [-5, 2, 6] and i == 6:
            return self.parent()([-6, 2])
        if self.value == [-5, 2, 6] and i == 2:
            return self.parent()([-2, -5, 4, 6])
        if self.value == [-6, 2] and i == 2:
            return self.parent()([-2, -6, 4])
        if self.value == [-2, 5] and i == 5:
            return self.parent()([-2, -5, 4, 6])
        if self.value == [-2, -5, 4, 6] and i == 6:
            return self.parent()([-2, -6, 4])
        if self.value == [-2, -5, 4, 6] and i == 4:
            return self.parent()([-4, 3, 6])
        if self.value == [-2, -6, 4] and i == 4:
            return self.parent()([-4, -6, 3, 5])
        if self.value == [-4, 3, 6] and i == 6:
            return self.parent()([-4, -6, 3, 5])
        if self.value == [-4, 3, 6] and i == 3:
            return self.parent()([-3, 1, 6])
        if self.value == [-4, -6, 3, 5] and i == 5:
            return self.parent()([-5, 3])
        if self.value == [-4, -6, 3, 5] and i == 3:
            return self.parent()([-3, -6, 1, 5])
        if self.value == [-5, 3] and i == 3:
            return self.parent()([-3, -5, 1, 4])
        if self.value == [-3, 1, 6] and i == 6:
            return self.parent()([-3, -6, 1, 5])
        if self.value == [-3, 1, 6] and i == 1:
            return self.parent()([-1, 6])
        if self.value == [-3, -6, 1, 5] and i == 5:
            return self.parent()([-3, -5, 1, 4])
        if self.value == [-3, -6, 1, 5] and i == 1:
            return self.parent()([-1, -6, 5])
        if self.value == [-3, -5, 1, 4] and i == 4:
            return self.parent()([-4, 1, 2])
        if self.value == [-3, -5, 1, 4] and i == 1:
            return self.parent()([-1, -5, 4])
        if self.value == [-4, 1, 2] and i == 2:
            return self.parent()([-2, 1])
        if self.value == [-4, 1, 2] and i == 1:
            return self.parent()([-1, -4, 2, 3])
        if self.value == [-2, 1] and i == 1:
            return self.parent()([-1, -2, 3])
        if self.value == [-1, 6] and i == 6:
            return self.parent()([-1, -6, 5])
        if self.value == [-1, -6, 5] and i == 5:
            return self.parent()([-1, -5, 4])
        if self.value == [-1, -5, 4] and i == 4:
            return self.parent()([-1, -4, 2, 3])
        if self.value == [-1, -4, 2, 3] and i == 2:
            return self.parent()([-1, -2, 3])
        if self.value == [-1, -4, 2, 3] and i == 3:
            return self.parent()([-3, 2])
        if self.value == [-1, -2, 3] and i == 3:
            return self.parent()([-2, -3, 4])
        if self.value == [-3, 2] and i == 2:
            return self.parent()([-2, -3, 4])
        if self.value == [-2, -3, 4] and i == 4:
            return self.parent()([-4, 5])
        if self.value == [-4, 5] and i == 5:
            return self.parent()([-5, 6])
        if self.value == [-5, 6] and i == 6:
            return self.parent()([-6])

        else:
            return None

class Crystal_of_letters_type_E6_element_dual(Letter, CrystalElement):
    r"""
    Type `E_6` crystal of letters elements. This crystal corresponds to the highest weight
    crystal `B(\Lambda_6)`. This crystal is dual to `B(\Lambda_1)` of type `E_6`.

    TESTS::

        sage: C = CrystalOfLetters(['E',6], dual = True)
	sage: C.module_generators
	[[6]]
	sage: all(b==b.retract(b.lift()) for b in C)
	True
	sage: C.list()
	[[6], [5, -6], [4, -5], [2, 3, -4], [3, -2], [1, 2, -3], [2, -1], [1, 4, -2, -3],
	[4, -1, -2], [1, 5, -4], [3, 5, -1, -4], [5, -3], [1, 6, -5], [3, 6, -1, -5], [4, 6, -3, -5],
	[2, 6, -4], [6, -2], [1, -6], [3, -1, -6], [4, -3, -6], [2, 5, -4, -6], [5, -2, -6], [2, -5],
	[4, -2, -5], [3, -4], [1, -3], [-1]]
        sage: C.check()
        True
	sage: all(b.f(i).e(i) == b for i in C.index_set() for b in C if b.f(i) is not None)
	True
	sage: all(b.e(i).f(i) == b for i in C.index_set() for b in C if b.e(i) is not None)
	True
	sage: G = C.digraph()
	sage: G.show(edge_labels=true, figsize=12, vertex_size=1)
    """

    def __repr__(self):
	"""
	In their full representation, the vertices of this crystal are labeled
	by their weight. For example vertex [-2,1] indicates that a 2-arrow
	is coming into this vertex, and a 1-arrow is leaving the vertex.
	Specifying the option element_print_style = 'compact' for a given crystal C,
	labels the vertices of this crystal by the 27 letters -ABCDEFGHIJKLMNOPQRSTUVWXYZ

	EXAMPLES::

            sage: C = CrystalOfLetters(['E',6], element_print_style = 'compact', dual = True)
	    sage: C.list()
	    [-, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z]
	    """
	if self.parent()._element_print_style == 'compact':
	    l=['-','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
	    return l[self.parent().list().index(self)]
	return "%s"%self.value

    def lift(self):
	"""
	Lifts an element of self to the crystal of letters CrystalOfLetters(['E',6])
	by taking its inverse weight.

	EXAMPLES::

	    sage: C = CrystalOfLetters(['E',6], dual = True)
	    sage: b=C.module_generators[0]
	    sage: b.lift()
	    [-6]
	"""
	return self.parent()._ambient([-i for i in self.value])

    def retract(self, p):
	"""
	Retracts element p, which is an element in CrystalOfLetters(['E',6]) to
	an element in CrystalOfLetters(['E',6], dual = True) by taking its inverse weight.

	EXAMPLES::

	    sage: C = CrystalOfLetters(['E',6])
	    sage: Cd = CrystalOfLetters(['E',6], dual = True)
	    sage: b = Cd.module_generators[0]
	    sage: p = C([-1,3])
	    sage: b.retract(p)
	    [1, -3]
	    sage: b.retract(None)
	"""
	if p is None:
	    return None
	return self.parent()([-i for i in p.value])

    def e(self, i):
        r"""
        Returns the action of `e_i` on self.

        EXAMPLES::

	    sage: C = CrystalOfLetters(['E',6], dual = True)
	    sage: C([-1]).e(1)
	    [1, -3]
        """
	return self.retract(self.lift().f(i))

    def f(self, i):
        r"""
        Returns the action of `f_i` on self.

        EXAMPLES::

	    sage: C = CrystalOfLetters(['E',6], dual = True)
	    sage: C([6]).f(6)
	    [5, -6]
	    sage: C([6]).f(1)
        """
	return self.retract(self.lift().e(i))

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['E',6], dual = True)
	    sage: b=C.module_generators[0]
	    sage: b.weight()
	    (0, 0, 0, 0, 1, -1/3, -1/3, 1/3)
            sage: [v.weight() for v in C]
	    [(0, 0, 0, 0, 1, -1/3, -1/3, 1/3),
	    (0, 0, 0, 1, 0, -1/3, -1/3, 1/3),
	    (0, 0, 1, 0, 0, -1/3, -1/3, 1/3),
	    (0, 1, 0, 0, 0, -1/3, -1/3, 1/3),
	    (-1, 0, 0, 0, 0, -1/3, -1/3, 1/3),
	    (1, 0, 0, 0, 0, -1/3, -1/3, 1/3),
	    (1/2, 1/2, 1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
	    (0, -1, 0, 0, 0, -1/3, -1/3, 1/3),
	    (-1/2, -1/2, 1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
	    (0, 0, -1, 0, 0, -1/3, -1/3, 1/3),
	    (-1/2, 1/2, -1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
	    (1/2, -1/2, -1/2, 1/2, 1/2, 1/6, 1/6, -1/6),
	    (0, 0, 0, -1, 0, -1/3, -1/3, 1/3),
	    (-1/2, 1/2, 1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
	    (1/2, -1/2, 1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
	    (1/2, 1/2, -1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
	    (-1/2, -1/2, -1/2, -1/2, 1/2, 1/6, 1/6, -1/6),
	    (0, 0, 0, 0, -1, -1/3, -1/3, 1/3),
	    (-1/2, 1/2, 1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
	    (1/2, -1/2, 1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
	    (1/2, 1/2, -1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
	    (-1/2, -1/2, -1/2, 1/2, -1/2, 1/6, 1/6, -1/6),
	    (1/2, 1/2, 1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
	    (-1/2, -1/2, 1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
	    (-1/2, 1/2, -1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
	    (1/2, -1/2, -1/2, -1/2, -1/2, 1/6, 1/6, -1/6),
	    (0, 0, 0, 0, 0, 2/3, 2/3, -2/3)]
        """
	return -self.lift().weight()

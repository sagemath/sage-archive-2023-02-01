r"""
Crystals of letters
"""

#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
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

from sage.structure.element    import Element
from sage.combinat.cartan_type import CartanType
from crystals                  import Crystal, ClassicalCrystal, CrystalElement
from sage.misc.flatten         import flatten

def CrystalOfLetters(type):
    r"""
    Return the crystal of letters of the given type

    INPUT:
        T -- A CartanType

    EXAMPLES:

        sage: C = CrystalOfLetters(['A',5])
        sage: C.list()
        [1, 2, 3, 4, 5, 6]

    TEST:

        sage: C.unrank(0) == C(1)  # todo: fix this test
        True

    """
    if type[0] == 'A':
	return ClassicalCrystalOfLetters(type,
                                         Crystal_of_letters_type_A_element)
    elif type[0] == 'C':
	return ClassicalCrystalOfLetters(type,
                                         Crystal_of_letters_type_C_element)
    else:
	print('not yet implemented')


class ClassicalCrystalOfLetters(ClassicalCrystal):
    r"""
    A generic class for classical crystals of letters.

    All classical crystals of letters should be instance of this class
    or of subclasses. To define a new crystal of letters, one only
    need to implement a class for the elements (which subclasses
    Letter and CrystalElement), with appropriate e and f
    operations. If the module generator is not 1, one also need to
    subclass ClassicalCrystalOfLetters for the crystal itself.

    The basic assumption is that crystals of letters are small, but
    used intensivelly as building blocks. All other operations (list,
    comparison, ...) can be derived automatically in the brute force
    way, with appropriate caching).

    """
    def __init__(self, type, element_class):
        self.cartanType = CartanType(type)
        self._name = "The crystal of letters for type %s"%type
        self.index_set = self.cartanType.index_set()
        self.element_class = element_class
        self.module_generators = [self(1)]
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
        assert x.parent() == self and y.parent() == self;
        if   self._digraph_closure.has_edge(x,y):
            return -1;
        elif self._digraph_closure.has_edge(y,x):
            return 1;
        else:
            return 0;

    # TODO: cmp, in, ...

# Utility. Note: much of this class should be factored out at some point!
class Letter(Element):
    r"""
    A class for letters

    TEST:
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

    def __init__(self, parent, value):
#        Element.__init__(self, parent);
        self._parent = parent
        self.value = value

    def parent(self):
        return self._parent  # Should be inherited from Element!

    def __repr__(self):
        return "%s"%self.value

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
               self.parent()  == other.parent()   and \
               self.value     == other.value

    def __cmp__(self, other):
        assert self.parent() == other.parent()
        return self.parent().cmp_elements(self, other)

#########################
# Type A
#########################

class Crystal_of_letters_type_A_element(Letter, CrystalElement):
    r"""
    Type A crystal of letters elements

    TEST:
        sage: C = CrystalOfLetters (['A',3])
        sage: C.list()
        [1, 2, 3, 4]
        sage: [ [x < y for y in C] for x in C ]
        [[False, True, True, True],
         [False, False, True, True],
         [False, False, False, True],
         [False, False, False, False]]

        sage: C = CrystalOfLetters(['A',5])
        sage: C(1) < C(1), C(1) < C(2), C(1) < C(3), C(2) < C(1)
        (False, True, True, False)

    """
    def e(self, i):
        r"""
        TEST:
            sage: C = CrystalOfLetters(['A',5])
            sage: C(1).e(1) == None
            True
            sage: C(2).e(1) == C(1)
            True
            sage: C(3).e(1) == None
            True
            sage: C(1).e(2) == None
            True
            sage: C(2).e(2) == None
            True
            sage: C(3).e(2) == C(2)
            True
        """
        assert i in self.index_set()
        if self.value == i+1:
            return self._parent(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        TESTS:
            sage: C = CrystalOfLetters(['A',5])
            sage: C(1).f(1) == C(2)
            True
            sage: C(2).f(1) == None
            True
            sage: C(3).f(1) == None
            True
            sage: C(1).f(2) == None
            True
            sage: C(2).f(2) == C(3)
            True
            sage: C(3).f(2) == None
            True
        """
        assert i in self.index_set()
        if self.value == i:
            return self._parent(self.value+1)
        else:
            return None


#########################
# Type C
#########################

class Crystal_of_letters_type_C_element(Letter, CrystalElement):
    r"""
    Type C crystal of letters elements

    TEST:
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

    """
    def e(self, i):
        r"""
        TEST:
            sage: C = CrystalOfLetters(['C',2])
            sage: C(1).e(1) == None
            True
            sage: C(2).e(1) == C(1)
            True
            sage: C(-2).e(1) == None
            True
            sage: C(-1).e(1) == C(-2)
            True
            sage: C(-2).e(2) == C(2)
            True
        """
        assert i in self.index_set()
    	if self.value == -self._parent.cartanType.n and self.value == -i:
	    return self._parent(-self.value)
	elif self.value == i+1 or self.value == -i:
	    return self._parent(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        TESTS:
            sage: C = CrystalOfLetters(['C',2])
            sage: C(1).f(1) == C(2)
            True
            sage: C(2).f(1) == None
            True
            sage: C(-2).f(1) == C(-1)
            True
            sage: C(-1).f(2) == None
            True
        """
        assert i in self.index_set()
	if self.value == self._parent.cartanType.n and self.value == i:
	    return  self._parent(-self.value)
        elif self.value == i or self.value == -i-1:
            return self._parent(self.value+1)
        else:
            return None

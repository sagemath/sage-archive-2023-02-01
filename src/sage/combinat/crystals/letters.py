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
from crystals                  import Crystal, CrystalElement
from sage.misc.flatten         import flatten

def CrystalOfLetters(type):
    r"""
    Return the crystal of letters of the given type

    INPUT:
        T -- A CartanType

    EXAMPLES:

        sage: C = CrystalOfLetters(['A',5])
        sage: C.list()
        [1, 2, 3, 4, 5]

    TEST:

        sage: C.unrank(0) == C(1)  # todo: fix this test
        True

    """
    if type[0] == 'A':
	return Crystal_of_letters_type_A(type)
    elif type[0] == 'C':
	return Crystal_of_letters_type_C(type)
    else:
	print('not yet implemented')

# This class should be factored out at some point!
#
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

#########################
# Type A
#########################

class Crystal_of_letters_type_A(Crystal):
    r"""
    Type A crystal of letters
    """
    def __init__(self, type):
        self.cartanType = CartanType(type)
        self._name = "The crystal of letters for type %s"%type
        self.index_set = self.cartanType.index_set()
        self.module_generators = [self(1)]

    def list(self):
        return [self(i) for i in range(1,self.cartanType.n+1)]

    def __call__(self, value):
        return Crystal_of_letters_type_A_element(self, value)

class Crystal_of_letters_type_A_element(Letter, CrystalElement):
    r"""
    Type A crystal of letters elements
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

class Crystal_of_letters_type_C(Crystal):
    r"""
    Type C crystal of letters
    """
    def __init__(self, type):
        self.cartanType = CartanType(type)
        self._name = "The crystal of letters for type %s"%type
        self.index_set = self.cartanType.index_set()
        self.module_generators = [self(1)]

    def list(self):
        return flatten([[self(i) for i in range(1,self.cartanType.n+1)], [self(-i) for i in range(1,self.cartanType.n+1)]])

    def __call__(self, value):
        return Crystal_of_letters_type_C_element(self, value)

class Crystal_of_letters_type_C_element(Letter, CrystalElement):
    r"""
    Type C crystal of letters elements
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

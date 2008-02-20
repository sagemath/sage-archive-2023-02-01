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
    if type[0] == 'B':
	return Crystal_of_letters_type_B(type)
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
        return [self(i) for i in range(1,self.cartanType.n+2)]

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


#########################
# Type B
#########################

class Crystal_of_letters_type_B(Crystal):
    r"""
    Type A crystal of letters
    """
    def __init__(self, type):
        self.cartanType = CartanType(type)
        self._name = "The crystal of letters for type %s"%type
        self.index_set = self.cartanType.index_set()
        self.module_generators = [self(1)]

    def list(self):
        return [self(i) for i in range(-self.cartanType.n,self.cartanType.n+1)]

    def __call__(self, value):
        return Crystal_of_letters_type_B_element(self, value)

class Crystal_of_letters_type_B_element(Letter, CrystalElement):
    r"""
    Type B crystal of letters elements
    """
    def e(self, i):
        r"""
        TEST:
            sage: C = CrystalOfLetters(['B',2])
            sage: C(1).e(1) == None
            True
            sage: C(1).e(2) == None
            True
            sage: C(2).e(1) == C(1)
            True
            sage: C(2).e(2) == None
            True
            sage: C(0).e(1) == None
            True
            sage: C(0).e(2) == C(2)
            True
            sage: C(-2).e(1) == None
            True
            sage: C(-2).e(2) == C(0)
            True
            sage: C(-1).e(1) == C(-2)
            True
            sage: C(-1).e(2) == None
            True
        """
        assert i in self.index_set()
        if self.value == i+1:
            return self._parent(i)
        elif self.value == 0 and i == self._parent.cartanType.n:
            return self._parent(self._parent.cartanType.n)
        elif self.value == -i:
           if i == self._parent.cartanType.n:
               return self._parent(0)
           else:
               return self._parent(-i-1)
        else:
            return None

    def f(self, i):
        r"""
        TESTS:
            sage: C = CrystalOfLetters(['B',2])
            sage: C(1).f(1) == C(2)
            True
            sage: C(1).f(2) == None
            True
            sage: C(2).f(1) == None
            True
            sage: C(2).f(2) == C(0)
            True
            sage: C(0).f(1) == None
            True
            sage: C(0).f(2) == C(-2)
            True
            sage: C(-2).f(1) == C(-1)
            True
            sage: C(-2).f(2) == None
            True
            sage: C(-1).f(1) == None
            True
            sage: C(-1).f(2) == None
            True
        """
        assert i in self.index_set()
        if self.value == i:
            if i < self._parent.cartanType.n:
                return self._parent(i+1)
            else:
                return self._parent(0)
        elif self.value == 0 and i == self._parent.cartanType.n:
            return self._parent(-self._parent.cartanType.n)
        elif self.value == -i-1:
            return(self._parent(-i))
        else:
            return None


#########################
# Type B spin
#########################

def SpinCrystal(type):
    r"""
    Return the spin crystal of the given type

    INPUT:
        T -- A CartanType (either B or D)

    EXAMPLES:

        sage: C = SpinCrystal(['B',3])
        sage: C.list()
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    """
    if type[0] == 'B':
	return Spin_crystal_type_B(type)
    else:
	print('not yet implemented')

    return("not yet implemented")


class Spin_crystal_type_B(Crystal):
    r"""
    Type B spin representation crystal
    """
    def __init__(self, type):
        self.cartanType = CartanType(type)
        self._name = "The spin crystal for type %s"%type
        self.index_set = self.cartanType.index_set()
        self.module_generators = [self(1)]

    def list(self):
        return [self(i) for i in range(1, 1+2**(self.cartanType.n))]

    def __call__(self, value):
        return Spin_crystal_type_B_element(self, value)

class Spin_crystal_type_B_element(Letter, CrystalElement):
    r"""
    Type B spin representation crystal element
    """
    def e(self, i):
        r"""
        TEST:
            sage: C = SpinCrystal(['B',3])
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
            if ((self.value-1)/2**(rank-i-1)).floor()%4 == 2:
                return self._parent(self.value-2**(rank-i-1))
        elif i == rank:
            if (self.value-1)%2 == 1:
                return self._parent(self.value-1)
        else:
            return None

    def f(self, i):
        r"""
        TESTS:
            sage: C = Crystal_type_B_spin(3)
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
            if ((self.value-1)/2**(rank-i-1)).floor()%4 == 1:
                return self._parent(self.value+2**(rank-i-1))
        elif i == rank:
            if (self.value-1)%2 == 0:
                return self._parent(self.value+1)
        else:
            return None

    def signature(self):
        sword = ""
        m = self.value-1
        for i in range(self._parent.cartanType.n,0,-1):
            if m >= 2**(i-1):
                sword = sword + "-"
                m = m - 2**(i-1)
            else:
                sword = sword + "+"
        return(sword)

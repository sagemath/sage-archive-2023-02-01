r"""
Spin Crystals
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

from sage.structure.element    import Element
from sage.combinat.cartan_type import CartanType
from crystals                  import Crystal, ClassicalCrystal, CrystalElement
from letters                   import Letter

#########################
# Type B spin
#########################

def CrystalOfSpins(type):
    r"""
    Return the spin crystal of the given type

    INPUT:
        T -- A CartanType (either B or D)

    Kashiwara and Nakashima represent the elements by sequences
    of signs +/-. In this implementation, the elements are represented
    by positive integers, but their notation can be recovered using
    signature()

    EXAMPLES:

        sage: C = CrystalOfSpins(['B',3])
        sage: C.list()
        [1, 2, 3, 4, 5, 6, 7, 8]

        sage: [x.signature() for x in C]
        ['+++', '++-', '+-+', '+--', '-++', '-+-', '--+', '---']

    """
    type = CartanType(type)
    if type[0] == 'B':
	return Spin_crystal_type_B(type)
    else:
	raise NotImplementedError

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
            sage: C = SpinCrystal(['B',3])
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

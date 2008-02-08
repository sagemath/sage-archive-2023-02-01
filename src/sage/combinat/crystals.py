"""
Crystals
"""
#*****************************************************************************
#       Copyright (C) 2007 Nicolas Thiery <nthiery at users.sf.net>,
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
from sage.combinat.combinat import CombinatorialClass
from sage.combinat.combinat import CombinatorialObject
from sage.combinat.cartan_type import CartanType

## Our Cat::Crystal
class Crystal(CombinatorialClass, Parent):
    # should implement module generators
    def bla(self):
        return 1;
    # list / ...

class CrystalElement(CombinatorialObject):
    # should implement e, f, ...
    def epsilon(self, i):
        x = self;
        eps = 0;
        while x is not None:
            x = x.e(i)
            eps = eps+1
        return eps


# crystals::letters dispatcher
def CrystalOfLetters(type):
    return Crystals_of_letters_type_A(type)


#
class Crystals_of_letters_type_A(Crystal):
    def __init__(self, type):
        self.cartanType = CartanType(type)
        self._name = "Crystal of letters for type %s"%type
        self.indexSet = type.index_set()

    def list(self):
        return range(1,self.cartanType.n+1)

#    def new(self, value):
#        return Crystals_of_letters_type_A_element(value,self)

    def moduleGenerators(self):
        return [1];

    def __call__(self, value):
        return Crystals_of_letters_type_A_element(self, value);


class Crystals_of_letters_type_A_element(Element):
    def __init__(self, parent, value):
#        Element.__init__(self, parent);
        self._parent = parent
        self.value = value

    def _repr_(self):
        return "%s"%self.value

    def e(self, i):
        assert i in self._parent.indexSet
        if self.value == i:
            return self._parent(self.value-1)
        else:
            return None

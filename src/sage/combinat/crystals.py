r"""
Crystals

Let $T$ be a CartanType with index set I, and W be a realization of the
type $T$ weight lattice.

A type $T$ crystal $C$ is an oriented graph equiped with a weight
function the nodes to some realization of the type $T$ weight lattice
such that:
\begin{itemize}
\item each edge has a label in I
\item for each i in I, each node x has:
    - at most one i-successor f_i(x)
    - at most one i-predecessor e_i(x)
   Furthermore, when the exists,
    - f_i(x).weight() = x.weight() - \alpha_i
    - e_i(x).weight() = x.weight() + \alpha_i

This crystal actually models a representation of a Lie algebra if it
satisfies some further local conditions due to Stembridge.

EXAMPLES:

We construct the type $A_5$ crystal on letters

    sage: C = CrystalOfLetters(['A',5]); C
    The crystal of letters for type ['A', 5]

It has a single highest weight element:
    sage: C.module_generators
    [1]

A crystal is a CombinatorialClass; and we can count and list its elements
in the usual way:
    sage: C.count()    # todo: not implemented
    5
    sage: C.list()
    [1, 2, 3, 4, 5]
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

from sage.structure.parent     import Parent
from sage.structure.element    import Element
from sage.combinat.combinat    import CombinatorialClass
from sage.combinat.combinat    import CombinatorialObject
from sage.combinat.cartan_type import CartanType

## MuPAD-Combinat's Cat::Crystal
class Crystal(CombinatorialClass, Parent):
    r"""
    The abstract class of crystals

    elements of a class implementing Crystal should have the following
    attributes
    \begin{itemize}
    \item cartan_type
    \item module_generators
    \item weight_lattice_realization
    \end{itemize}
    """

    def bla(self):
        return 1;
    # list / ...

class CrystalElement(Element):
    r"""
    The abstract class of crystal elements

    Sub classes should implement:
    \begin{itemize}
    \item x.e(i)        (returning $e_i(x)$)
    \item x.f(i)        (returning $f_i(x)$)
    \item x.weight()
    \end{itemize}
    """

    def e(self, i):
        r"""
        Returns $e_i(x)$ if it exists or None otherwise
        """
        raise NotImplementedError

    # should implement e, f, ...
    def epsilon(self, i):
        x = self;
        eps = 0;
        while x is not None:
            x = x.e(i)
            eps = eps+1
        return eps


def CrystalOfLetters(type):
    r"""
    Return the crystal of letters of the given type

    INPUT:
        T -- A CartanType

    EXAMPLES:

        sage: C = CrystalOfLetters(['A',5])
        sage: C.list()
        [1, 2, 3, 4, 5]
    """
    return Crystals_of_letters_type_A(type)

class Crystals_of_letters_type_A(Crystal):
    r"""
    Type A crystal of letters
    """
    def __init__(self, type):
        self.cartanType = CartanType(type)
        self._name = "The crystal of letters for type %s"%type
        self.indexSet = self.cartanType.index_set()
        self.module_generators = [self(1)]

    def list(self):
        return range(1,self.cartanType.n+1)

    def __call__(self, value):
        return Crystals_of_letters_type_A_element(self, value);

class Letters(Element):
    r"""
    A generic class for letters
    """

    def __init__(self, parent, value):
#        Element.__init__(self, parent);
        self._parent = parent
        self.value = value

    def __repr__(self):
        return "%s"%self.value


class Crystals_of_letters_type_A_element(Letters, CrystalElement):
    r"""
    Type A crystal of letters elements
    """
    def e(self, i):
        assert i in self.indexSet()
        if self.value == i:
            return self._parent(self.value-1)
        else:
            return None

#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
#*****************************************************************************
from sage.rings.all import Integer
import root_system
from cartan_matrix import cartan_matrix

def CartanType(t):
    """
    Returns an object corresponding to the Cartan type t.

    EXAMPLES:
        sage: CartanType(['A',4])
        ['A', 4]
    """
    if isinstance(t, CartanType_simple):
        return t;
    else:
        return CartanType_simple(t)

class CartanType_simple:
    def __init__(self, t):
        """
        EXAMPLES:
            sage: ct = CartanType(['A',4])
            sage: ct == loads(dumps(ct))
            True
        """
        self.t = t

        self.letter = t[0]
        self.n = t[1]

        if len(t) > 2:
            self.affine = t[2]
        else:
            self.affine = None

    def __repr__(self):
        """
        TESTS:
            sage: ct = CartanType(['A',3])
            sage: repr(ct)
            "['A', 3]"
        """
        if self.affine is None:
            return "['%s', %s]"%(self.letter, self.n)
        else:
            return "['%s', %s, %s]"%(self.letter, self.n, self.affine)

    def __getitem__(self, x):
        """
        EXAMPLES:
            sage: t = CartanType(['A', 3, 1])
            sage: t[0]
            'A'
            sage: t[1]
            3
            sage: t[2]
            1
            sage: t[3]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self.t[x]

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: ct1 = CartanType(['A',4])
            sage: ct2 = CartanType(['A',4])
            sage: ct3 = CartanType(['A',5])
            sage: ct1 == ct2
            True
            sage: ct1 != ct3
            True
        """
        if other.__class__ != self.__class__:
            return cmp(self.__class__, other.__class__)

        if other.letter != self.letter:
            return cmp(self.letter, other.letter)
        if other.affine != self.affine:
            return cmp(self.affine, other.affine)
        if other.n != self.n:
            return cmp(self.n, other.n)

        return 0


    def __len__(self):
        """
        EXAMPLES:
            sage: len(CartanType(['A',4]))
            2
            sage: len(CartanType(['A',4,1]))
            3
        """
        if self.affine:
            return 3
        else:
            return 2

    def is_finite(self):
        """
        Returns True if self is finite.

        EXAMPLES:
            sage: CartanType(['A',4]).is_finite()
            True
            sage: CartanType(['A',4,1]).is_finite()
            False
        """
        if self.affine is not None:
            return False

        if self.letter in ['A', 'B', 'C', 'D', 'I']:
            return True

        if self.letter == 'E':
            return self.n <= 8

        if self.letter == 'F':
            return self.n <= 4

        if self.letter == 'G':
            return self.n <= 2

    def is_affine(self):
        """
        Returns True if self is affine.

        EXAMPLES:
            sage: CartanType(['A', 3]).is_affine()
            False
            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """
        return self.affine is not None


    def rank(self):
        """
        Returns the rank of self.
        EXAMPLES:
            sage: CartanType(['A', 4]).rank()
            4
            sage: CartanType(['A', 7, 2]).rank()
            4
            sage: CartanType(['I', 8]).rank()
            2
        """
        if self.is_affine():
            if self.affine == 3 and self.letter == 'D':
                return self.n-1
            elif self.affine == 2 and self.letter == 'A':
                ## FIXME: check in the literature what should be the
		## appropriate definition for rank
                return int(self.n+1)/2
            else:
                return self.n
        else:
            if self.letter == "I":
                return 2
            else:
                return self.n

    def index_set(self):
        """
        Returns the index set for self.

        EXAMPLES:
            sage: CartanType(['A', 3, 1]).index_set()
            [0, 1, 2, 3]
            sage: CartanType(['D', 4]).index_set()
            [1, 2, 3, 4]
        """
        n = self.rank()
        if self.is_affine():
            return range(n+1)
        else:
            return range(1, n+1)

    def root_system(self):
        """
        Returns the root system associated to self.

        EXAMPLES:
            sage: CartanType(['A',4]).root_system()
            Root system of type ['A', 4]
        """
        return root_system.RootSystem(self)

    def cartan_matrix(self):
        """
        Returns the Cartan matrix associated with self.

        EXAMPLES:
            sage: CartanType(['A',4]).cartan_matrix()
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]

        """
        return cartan_matrix(self)

    def type(self):
        """
        Returns the type of self.

        EXAMPLES:
            sage: CartanType(['A', 4]).type()
            'A'
        """
        return self.letter

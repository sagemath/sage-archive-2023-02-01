"""
Root system data for dual Cartan types
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Anne Schilling <anne at math.ucdavis.edu>
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.sage_object import SageObject
from sage.misc.misc import attrcall
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.root_system.cartan_type import CartanType_crystalographic
import sage

class CartanType(UniqueRepresentation, SageObject, CartanType_crystalographic):
    r"""
    A class for dual Cartan types

    The dual of a (crystalographic) Cartan type is a Cartan type with
    the same index set, but all arrows reversed in the Dynkin diagram
    (otherwise said, the Cartan matrix is transposed). It shares a lot
    of properties in common with its dual. In particular, the Weyl
    group is isomorphic to that of the dual as a Coxeter group.

    Note: ``CartanType(['BC',4,2]).dual()`` is implemented by
    relabelling, not duality!

    """
    def __init__(self, type):
        """
        INPUT:

        - ``type`` -- a Cartan type

        EXAMPLES::

           sage: ct = CartanType(['F',4]).dual()
           sage: ct == loads(dumps(ct))
           True

        TESTS::

            sage: ct1 = CartanType(['B',2]).dual()
            sage: ct2 = CartanType(['B',2]).dual()
            sage: ct3 = CartanType(['D',4]).dual()
            sage: ct1 == ct2
            True
            sage: ct1 == ct3
            False

        """
        assert type.is_crystalographic()
        self._dual = type
        if type.is_affine():
            self._add_abstract_superclass(CartanType_affine)

    def _repr_(self, compact = False):
        """
        EXAMPLES::

           sage: CartanType(['F', 4]).dual()
           ['F', 4]^*

           sage: CartanType(['F', 4]).dual()._repr_(compact = True)
           'F4*'
        """
        return self.dual()._repr_(compact)+("*" if compact else "^*")

    def __reduce__(self):
        """
        TESTS::

            sage: CartanType(['F', 4]).dual().__reduce__()
            (*.dual(), (['F', 4],))

        """
        return (attrcall("dual"), (self._dual,))

    def ascii_art(self, label = lambda x: x):
        """
        Returns an ascii art representation of this Cartan type

        (by hacking the ascii art representation of the dual cartan type)

        EXAMPLES::

            sage: print CartanType(["G", 2]).dual().ascii_art()
              3
            O=>=O
            1   2
            sage: print CartanType(["F", 4]).dual().ascii_art()
            O---O=<=O---O
            1   2   3   4
            sage: print CartanType(["B", 3, 1]).dual().ascii_art()
                O 0
                |
                |
            O---O=<=O
            1   2   3
            sage: print CartanType(["C", 4, 1]).dual().ascii_art()
            O=<=O---O---O=>=O
            0   1   2   3   4
            sage: print CartanType(["G", 2, 1]).dual().ascii_art()
              3
            O=>=O---O
            1   2   0
            sage: print CartanType(["F", 4, 1]).dual().ascii_art()
            O---O---O=<=O---O
            0   1   2   3   4

            sage: print CartanType(["BC", 4, 2]).dual().ascii_art()
            O=>=O---O---O=>=O
            0   1   2   3   4
        """
        res = self.dual().ascii_art(label)
        # swap, like a computer science freshman!
        # This assumes that the oriented multiple arrows are always ascii arted as =<= or =>=
        res = res.replace("=<=", "=?=")
        res = res.replace("=>=", "=<=")
        res = res.replace("=?=", "=>=")
        return res

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: B4     = CartanType(['B', 4])
            sage: B4dual = CartanType(['B', 4]).dual()
            sage: F4dual = CartanType(['F', 4]).dual()
            sage: cmp(F4dual, F4dual)
            0

        Whether ``cmp()`` returns 1 or -1 doesn't matter, just check
        that the following are non-zero::

            sage: cmp(F4dual, B4dual) != 0
            True
            sage: cmp(B4dual, F4dual) * cmp(F4dual, B4dual) < 0
            True
            sage: cmp(B4dual, B4) != 0
            True
        """
        if other.__class__ != self.__class__:
            return cmp(self.__class__, other.__class__)
        return cmp(self._dual, other._dual)

    def dual(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.dual()
           ['F', 4]
        """
        return self._dual

    def is_irreducible(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_irreducible()
           True
        """
        return self._dual.is_irreducible()

    def is_finite(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_finite()
           True
        """
        return self._dual.is_finite()

    def is_crystalographic(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_crystalographic()
           True
        """
        return self._dual.is_crystalographic()

    def is_affine(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.is_affine()
           False
        """
        return self._dual.is_affine()

    def rank(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.rank()
           4
        """
        return self._dual.rank()

    def index_set(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.index_set()
           [1, 2, 3, 4]
        """
        return self._dual.index_set()

    def dynkin_diagram(self):
        """
        EXAMPLES::

           sage: ct = CartanType(['F', 4]).dual()
           sage: ct.dynkin_diagram()
           O---O=<=O---O
           1   2   3   4
           F4*
        """
        return self._dual.dynkin_diagram().dual()


###########################################################################
from cartan_type import CartanType_affine
class CartanType_affine(CartanType_affine):
    def classical(self):
        """
        Returns the classical Cartan type associated with self (which should be affine)

        EXAMPLES::

            sage: CartanType(['A',3,1]).dual().classical()
            ['A', 3]
            sage: CartanType(['B',3,1]).dual().classical()
            ['C', 3]
            sage: CartanType(['F',4,1]).dual().classical()
            ['F', 4]^*
            sage: CartanType(['BC',4,2]).dual().classical()
            ['B', 4]
        """
        return self.dual().classical().dual()

    def special_node(self):
        """
        Implements :meth:`CartanType_affine.special_node`

        The special node of the dual of an affine type `T` is the
        special node of `T`.

        EXAMPLES::

            sage: CartanType(['A',3,1]).dual().special_node()
            0
            sage: CartanType(['B',3,1]).dual().special_node()
            0
            sage: CartanType(['F',4,1]).dual().special_node()
            0
            sage: CartanType(['BC',4,2]).dual().special_node()
            0
        """
        return self.dual().special_node()

#class ambient_space(AmbientSpace):
# todo?

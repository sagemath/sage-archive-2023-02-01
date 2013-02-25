"""
Root system data for dual Cartan types
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Anne Schilling <anne at math.ucdavis.edu>
#       Copyright (C) 2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.misc import attrcall
from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject
from sage.combinat.root_system.cartan_type import CartanType_crystalographic, CartanType_finite, CartanType_affine, CartanType_simple
import sage
import ambient_space

class CartanType(UniqueRepresentation, SageObject, CartanType_crystalographic):
    r"""
    A class for dual Cartan types

    The dual of a (crystalographic) Cartan type is a Cartan type with
    the same index set, but all arrows reversed in the Dynkin diagram
    (otherwise said, the Cartan matrix is transposed). It shares a lot
    of properties in common with its dual. In particular, the Weyl
    group is isomorphic to that of the dual as a Coxeter group.

    EXAMPLES:

    For most finite Cartan types, and in particular the simply laced
    ones, the dual Cartan type is given by another preexisting Cartan
    type::

        sage: CartanType(['A',4]).dual()
        ['A', 4]
        sage: CartanType(['B',4]).dual()
        ['C', 4]
        sage: CartanType(['C',4]).dual()
        ['B', 4]

    So to exercise this class we need to consider the non simply laced
    exceptionnal or affine Cartan types::

        sage: F4d = CartanType(['F', 4]).dual(); F4d
        ['F', 4]^*
        sage: G21d = CartanType(['G',2,1]).dual(); G21d
        ['G', 2, 1]^*
        sage: B41d = CartanType(['B',4,1]).dual(); B41d
        ['B', 4, 1]^*

    They share many properties with their original Cartan types::

        sage: F4d.is_irreducible()
        True
        sage: F4d.is_crystalographic()
        True
        sage: F4d.is_simply_laced()
        False
        sage: F4d.is_finite()
        True
        sage: G21d.is_finite()
        False
        sage: F4d.is_affine()
        False
        sage: G21d.is_affine()
        True

    TESTS::

        sage: TestSuite(F4d).run()
    """
    def __init__(self, type):
        """
        INPUT:

        - ``type`` -- a Cartan type

        EXAMPLES::

           sage: ct = CartanType(['F',4]).dual()
           sage: TestSuite(ct).run()

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
        # TODO: design an appropriate infrastructure to handle this
        # automatically? Maybe using categories and axioms?
        if type.is_finite():
            self._add_abstract_superclass(CartanType_finite)
        elif type.is_affine():
            self._add_abstract_superclass(CartanType_affine)
        if type.is_irreducible():
            self._add_abstract_superclass(CartanType_simple)
        # No need to check non crystalographic types (they have no dual)
        # or simply laced types (they are self-dual)

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
        Return an ascii art representation of this Cartan type

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

class AmbientSpace(ambient_space.AmbientSpace):
    """
    Ambient space for a dual finite Cartan type.

    It is constructed in the canonical way from the ambient space of
    the original Cartan type by switching the roles of simple roots,
    fundamental weights, etc.

    EXAMPLES::

        sage: F4 = CartanType(["F",4])
        sage: L = F4.dual().root_system().ambient_space(); L
        Ambient space of the Root system of type ['F', 4]^*
        sage: TestSuite(L).run()
    """

    def dimension(self):
        """
        Return the dimension of this ambient space.

        .. SEEALSO:: :meth:`sage.combinat.root_system.ambient_space.AmbientSpace.dimension`

        EXAMPLES::

            sage: F4 = CartanType(["F",4])
            sage: F4.dual().root_system().ambient_space().simple_root(1)
            (0, 1, -1, 0)
        """
        return self.root_system.dual.ambient_space().dimension()

    @cached_method
    def simple_root(self, i):
        """
        Return the ``i``-th simple root.

        It is constructed by looking up the corresponding simple
        coroot in the ambient space for the dual Cartan type.

        EXAMPLES::

            sage: F4 = CartanType(["F",4])
            sage: F4.dual().root_system().ambient_space().simple_root(1)
            (0, 1, -1, 0)

            sage: F4.dual().root_system().ambient_space().simple_roots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 2), 4: (1, -1, -1, -1)}

            sage: F4.root_system().ambient_space().simple_coroots()
            Finite family {1: (0, 1, -1, 0), 2: (0, 0, 1, -1), 3: (0, 0, 0, 2), 4: (1, -1, -1, -1)}
        """
        K = self.base_ring()
        dual_coroot = self.cartan_type().dual().root_system().ambient_space(K).simple_coroot(i)
        return self.sum_of_terms(dual_coroot)

    @cached_method
    def fundamental_weights(self):
        """
        Return the fundamental weights.

        They are computed from the simple roots by inverting the
        Cartan matrix. This is acceptable since this is only about
        ambient spaces for finite Cartan types. Also, we do not have
        to worry about the usual `GL_n` vs `SL_n` catch because type
        `A` is self dual.

        An alternative would have been to start from the fundamental
        coweights in the dual ambient space, but those are not yet
        implemented.

        EXAMPLES::

            sage: L = CartanType(["F",4]).dual().root_system().ambient_space()
            sage: L.fundamental_weights()
            Finite family {1: (1, 1, 0, 0), 2: (2, 1, 1, 0), 3: (3, 1, 1, 1), 4: (2, 0, 0, 0)}
        """
        return self.fundamental_weights_from_simple_roots()

class CartanType_finite(sage.combinat.root_system.cartan_type.CartanType_finite):
    AmbientSpace = AmbientSpace

###########################################################################
class CartanType_affine(sage.combinat.root_system.cartan_type.CartanType_affine):
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
        Implement :meth:`CartanType_affine.special_node`

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

"""
Root system data for type A infinity
"""
#*****************************************************************************
# Copyright (C) 2016 Andrew Mathas <Andrew dot Mathas at Sydney dot edu dot au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from cartan_type import CartanType_simple
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

class CartanType(UniqueRepresentation, SageObject, CartanType_simple):
    r"""
    Define the Cartan type `A_oo`.

    We do not inherit from CartanType_crystallographic because it provides
    methods that are not yet implemented.
    """
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct=CartanType(['A',oo])
            sage: CartanType(['A',oo]) is CartanType(['A', Infinity])
            True
            sage: ct
            ['A', oo]
            sage: ct._repr_(compact = True)
            'Aoo'
            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            False
            sage: ct.is_affine()
            False
            sage: ct.is_untwisted_affine()
            False
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.dual()
            ['A', oo]

        TESTS::

            sage: TestSuite(ct).run()
        """
        super(CartanType, self).__init__()
        assert n == Infinity
        self.letter = 'A'
        self.n = n

    def _repr_(self, compact = False):
        """
        TESTS::

            sage: CartanType(['A',oo])
            ['A', oo]
            sage: CartanType(['A',oo])._repr_(compact=True)
            'Aoo'
        """
        format = '%s%s' if compact else "['%s', %s]"
        return format%(self.letter, 'oo')

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: ct = CartanType(['A',oo])
            sage: latex(ct)
            A_{\infty}
        """
        return 'A_{\infty}'

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: t = CartanType(['A', oo])
            sage: t[0]
            'A'
            sage: t[1]
            +Infinity
            sage: t[2]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i == 0:
            return self.letter
        elif i==1:
            return self.n
        else:
            raise IndexError("index out of range")

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CartanType(['A', oo]))
            2
        """
        return 2

    def ascii_art(self, label=lambda i: i, node=None):
        """
        Return an ascii art representation of the extended Dynkin diagram.

        EXAMPLES::

            sage: print(CartanType(['A', oo]).ascii_art())
            ...---O---O---O---O---O---O---O---...
                 -3  -2  -1   0   1   2   3

        """
        if node is None:
            node = self._ascii_art_node
        ret  = '...---'+'---'.join(node(label(i)) for i in range(7))+'---...\n'
        ret += '   '+''.join("{:4}".format(label(i)) for i in range(-3,4))
        return ret

    def dynkin_diagram(self):
        """
        Not implemented!

        The Dynkin diagram associated with ``self`` is the graph with vertex set
        `ZZ` and edges `i --> i+1`, for all `i\in ZZ`. This is not implemented
        because graphs with infinite vertex sets are not yet supported.

        EXAMPLES::

            sage: CartanType(['A',oo]).dynkin_diagram()
            Traceback (most recent call last):
            ...
            NotImplementedError: The Dynkin diagram of type A_oo is not implemented

        """
        raise NotImplementedError('The Dynkin diagram of type A_oo is not implemented')

    def cartan_matrix(self):
        """
        Not implemented!

        The Cartan matrix of ``self`` is the matrix with rows and columns indexed
        by `ZZ` and entries `(a_{ij})` given by:

        .. MATH::

               a_{ij} = \begin{cases}
                     2,& \text{if }i=j,\\
                    -1,& \text{if }|i-j|=1,\\
                     0,& \text{otherwise}.
                \end{cases}

        This is not implemented because matrices with an infinite number of rows
        or columns are not yet supported.

        EXAMPLES::

            sage: CartanType(['A', oo]).cartan_matrix()
            Traceback (most recent call last):
            ...
            NotImplementedError: The Cartan matrix of type A_oo is not implemented
        """
        raise NotImplementedError('The Cartan matrix of type A_oo is not implemented')

    def is_simply_laced(self):
        """
        Return whether ``self`` is simply laced, which is ``True``.

        EXAMPLES::

            sage: CartanType(['A', oo]).is_simply_laced()
            True
        """
        return True

    def dual(self):
        """
        Simply laced Cartan types are self-dual, so return ``self``.

        EXAMPLES::

            sage: CartanType(["A", oo]).dual()
            ['A', oo]
        """
        return self

    def index_set(self):
        """
        Implements :meth:`CartanType_abstract.index_set`.

        The index set for all standard finite Cartan types is of the form
        `\{1, \ldots, n\}`. (See :mod:`~sage.combinat.root_system.type_I`
        for a slight abuse of this).

        EXAMPLES::

            sage: CartanType(['A', oo]).index_set()
            Integer Ring
        """
        return ZZ

    def is_crystallographic(self):
        """
        Implements :meth:`CartanType_abstract.is_crystallographic`
        by returning ``True``.

        EXAMPLES::

            sage: CartanType(['A', oo]).is_crystallographic()
            True
        """
        return True

    def is_finite(self):
        """
        EXAMPLES::

            sage: CartanType(['A', oo]).is_finite()
            False
        """
        return False

    def is_affine(self):
        """
        EXAMPLES::

            sage: CartanType(['A', oo]).is_affine()
            False
        """
        return False

    def is_untwisted_affine(self):
        """
        Return whether ``self`` is untwisted affine

        A Cartan type is untwisted affine if it is the canonical
        affine extension of some finite type. Every affine type is
        either untwisted affine, dual thereof, or of type ``BC``.

        EXAMPLES::

            sage: CartanType(['A', oo]).is_untwisted_affine()
            False
        """
        return False

    def rank(self):
        """
        Return the rank of ``self`` which for type `X_n` is `n`.

        EXAMPLES::

            sage: CartanType(['A', oo]).rank()
            +Infinity
        """
        return self.n

    def type(self):
        """
        Returns the type of ``self``.

        EXAMPLES::

            sage: CartanType(['A', oo]).type()
            'A'
        """
        return self.letter

"""
Special Linear Groups

AUTHORS:

- William Stein: initial version

- David Joyner (2006-05): added examples, _latex_, __str__, gens,
  as_matrix_group

- William Stein (2006-12-09): rewrite

EXAMPLES::

    sage: SL(2, ZZ)
    Special Linear Group of degree 2 over Integer Ring
    sage: G = SL(2,GF(3)); G
    Special Linear Group of degree 2 over Finite Field of size 3
    sage: G.is_finite()
    True
    sage: G.conjugacy_class_representatives()
    [
    [1 0]
    [0 1],
    [0 2]
    [1 1],
    [0 1]
    [2 1],
    [2 0]
    [0 2],
    [0 2]
    [1 2],
    [0 1]
    [2 2],
    [0 2]
    [1 0]
    ]
    sage: G = SL(6,GF(5))
    sage: G.gens()
    [
    [2 0 0 0 0 0]
    [0 3 0 0 0 0]
    [0 0 1 0 0 0]
    [0 0 0 1 0 0]
    [0 0 0 0 1 0]
    [0 0 0 0 0 1],
    [4 0 0 0 0 1]
    [4 0 0 0 0 0]
    [0 4 0 0 0 0]
    [0 0 4 0 0 0]
    [0 0 0 4 0 0]
    [0 0 0 0 4 0]
    ]
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import is_FiniteField, Integer, FiniteField
from matrix_group import MatrixGroup_gap, MatrixGroup_gap_finite_field
from matrix_group_element import MatrixGroupElement

def SL(n, R, var='a'):
    r"""
    Return the special linear group of degree `n` over the ring
    `R`.

    EXAMPLES::

        sage: SL(3,GF(2))
        Special Linear Group of degree 3 over Finite Field of size 2
        sage: G = SL(15,GF(7)); G
        Special Linear Group of degree 15 over Finite Field of size 7
        sage: G.order()
        1956712595698146962015219062429586341124018007182049478916067369638713066737882363393519966343657677430907011270206265834819092046250232049187967718149558134226774650845658791865745408000000
        sage: len(G.gens())
        2
        sage: G = SL(2,ZZ); G
        Special Linear Group of degree 2 over Integer Ring
        sage: G.gens()
        [
        [ 0  1]
        [-1  0],
        [1 1]
        [0 1]
        ]

    Next we compute generators for `\mathrm{SL}_3(\mathbb{Z})`.

    ::

        sage: G = SL(3,ZZ); G
        Special Linear Group of degree 3 over Integer Ring
        sage: G.gens()
        [
        [0 1 0]
        [0 0 1]
        [1 0 0],
        [ 0  1  0]
        [-1  0  0]
        [ 0  0  1],
        [1 1 0]
        [0 1 0]
        [0 0 1]
        ]
    """
    if isinstance(R, (int, long, Integer)):
        R = FiniteField(R, var)
    if is_FiniteField(R):
        return SpecialLinearGroup_finite_field(n, R)
    else:
        return SpecialLinearGroup_generic(n, R)

class SpecialLinearGroup_generic(MatrixGroup_gap):
    def _gap_init_(self):
        """
        String to create this grop in GAP.

        EXAMPLES::

            sage: G = SL(6,GF(5)); G
            Special Linear Group of degree 6 over Finite Field of size 5
            sage: G._gap_init_()
            'SL(6, GF(5))'
        """
        return "SL(%s, %s)"%(self.degree(), self.base_ring()._gap_init_())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: G = SL(6,GF(5))
            sage: latex(G)
            \text{SL}_{6}(\mathbf{F}_{5})
        """
        return "\\text{SL}_{%s}(%s)"%(self.degree(), self.field_of_definition()._latex_())

    def _repr_(self):
        """
        Text representation of self.

        EXAMPLES::

            sage: SL(6,GF(5))
            Special Linear Group of degree 6 over Finite Field of size 5
        """
        return "Special Linear Group of degree %s over %s"%(self.degree(), self.base_ring())

    def __call__(self, x):
        """
        Construct a new element in this group, i.e. try to coerce x into
        self if at all possible.

        EXAMPLES::

            sage: G = SL(3, ZZ)
            sage: x = [[1,0,1], [0,1,0], [0,0,1]]
            sage: G(x)
            [1 0 1]
            [0 1 0]
            [0 0 1]
        """
        if isinstance(x, MatrixGroupElement) and x.parent() is self:
            return x
        try:
            m = self.matrix_space()(x)
        except TypeError:
            raise TypeError, "Cannot coerce %s to a %s-by-%s matrix over %s"%(x,self.degree(),self.degree(),self.base_ring())
        if m.determinant() == self.base_ring()(1):
            return MatrixGroupElement(m, self)
        else:
            raise TypeError, "%s does not have determinant 1"%(x)

    def __contains__(self, x):
        """
        Return True if x is an element of self, False otherwise.

        EXAMPLES::

            sage: G = SL(2, GF(101))
            sage: x = [[1,1], [0,1]]
            sage: x in G
            True

        ::

            sage: G = SL(3, ZZ)
            sage: x = [[1,0,1], [0,-1,0], [0,0,1]]
            sage: x in G
            False
        """
        try:
            x = self(x)
        except TypeError:
            return False
        return True



class SpecialLinearGroup_finite_field(SpecialLinearGroup_generic, MatrixGroup_gap_finite_field):
    pass


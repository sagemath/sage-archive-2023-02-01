r"""
General Linear Groups

EXAMPLES:
    sage: GL(4,QQ)
    General Linear Group of degree 4 over Rational Field
    sage: GL(1,ZZ)
    General Linear Group of degree 1 over Integer Ring
    sage: GL(100,RR)
    General Linear Group of degree 100 over Real Field with 53 bits of precision
    sage: GL(3,GF(49))
    General Linear Group of degree 3 over Finite Field in a of size 7^2

AUTHORS:
    -- David Joyner (2006-01)
    -- William Stein (2006-01)
    -- David Joyner (2006-05) - added _latex_, __str__, examples

TODO:
  Write a method to coerce GL into a MatrixGroup...

"""

##TODO: Rework this and \code{special_linear} into MatrixGroup class for any
##field, wrapping all of GAP's matrix group commands in chapter 41
##Matrix Groups of the GAP reference manual.


#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

from sage.groups.group import Group
from sage.rings.all import Integer, is_Ring
from sage.rings.finite_field import is_FiniteField
from sage.matrix.matrix_space import MatrixSpace
from linear import LinearGroup_generic, LinearGroup_finite_field

def GL(n, R):
    """
    Return the general linear group of degree $n$ over the ring $R$.

    EXAMPLES:
        sage: G = GL(6,GF(5))
        sage: G.order()
        11064475422000000000000000L
        sage: G.base_ring()
        Finite Field of size 5

        sage: F = GF(3); MS = MatrixSpace(F,2,2)
        sage: gens = [MS([[0,1],[1,0]]),MS([[1,1],[0,1]])]
        sage: G = MatrixGroup(gens)
        sage: G.order()
        48
        sage: H = GL(2,F)
        sage: H.order()
        48
        sage: H == G
        False
        sage: H.as_matrix_group() == G
        False
        sage: H.gens()
        [[2 0]
         [0 1], [2 1]
                [2 0]]

    """
    if not is_Ring(R):
        raise TypeError, "R (=%) must be a ring"%R
    if is_FiniteField(R):
        return GeneralLinearGroup_finite_field(n, R)
    else:
        return GeneralLinearGroup_generic(n, R)

class GeneralLinearGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        """
        EXAMPLES:
            sage: G = GL(6,GF(5))
            sage: G._gap_init_()
            'GL(6, 5)'

        """
        return "GL(%s, %s)"%(self.degree(), self.base_ring().order())

    def _latex_(self):
        """
        EXAMPLES:
            sage: G = GL(6,GF(5))
            sage: G._latex_()
            'GL$(6, GF(5))$'

        """
        return "GL$(%s, GF(%s))$"%(self.degree(), self.base_ring().order())

    def __str__(self):
        """
        EXAMPLES:
            sage: G = GL(6,GF(5))
            sage: print G
            GL(6, GF(5))

        """
        return "GL(%s, GF(%s))"%(self.degree(), self.base_ring().order())

    def __repr__(self):
        return "General Linear Group of degree %s over %s"%(self.degree(), self.base_ring())

    def gens(self):
        """
        EXAMPLES:
            sage: G = GL(6,GF(5))
            sage: G.gens()
            [[2 0 0 0 0 0]
             [0 1 0 0 0 0]
             [0 0 1 0 0 0]
             [0 0 0 1 0 0]
             [0 0 0 0 1 0]
             [0 0 0 0 0 1],
             [4 0 0 0 0 1]
             [4 0 0 0 0 0]
             [0 4 0 0 0 0]
             [0 0 4 0 0 0]
             [0 0 0 4 0 0]
             [0 0 0 0 4 0]]
        """
        from sage.interfaces.all import gap
        F = self.base_ring()
        G = self._gap_init_()
        n = eval(gap.eval("Length(GeneratorsOfGroup(%s))"%G))
        gens = [gap("GeneratorsOfGroup(%s)[%s]"%(G,i))._matrix_(F) for i in range(1,n+1)]
        return gens

    def as_matrix_group(self):
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        gns = self.gens()
        G = MatrixGroup(gns)
        return G

class GeneralLinearGroup_finite_field(LinearGroup_finite_field, GeneralLinearGroup_generic):
    pass

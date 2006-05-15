r"""
 Unitary Groups $GU(n,q)$ and $SU(n,q)$

These are $n \times n$ unitary matrices with entries in $GF(q^2)$.

AUTHOR:
    -- David Joyner: initial version (2006-3), modified from
                     special_linear (by W. Stein)
    -- David Joyner (2006-05): minor additions (examples, _latex_,
                                            __str__, gens, as_matrix_group)

EXAMPLES:
    sage: G = SU(3,GF(5))
    sage: G.order()
    378000
    sage: G
    Special Unitary Group of degree 3 over Finite Field of size 5
    sage: G._gap_init_()
    'SU(3, 5)'
    sage: G.random()
    [      3 2*a + 1   a + 1]
    [2*a + 2       4       a]
    [4*a + 2   a + 4       1]
    sage: G.base_ring()
    Finite Field of size 5
    sage: G.field_of_definition()
    Finite Field in a of size 5^2

"""

#*********************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*********************************************************************************

from sage.rings.all import IntegerRing, is_FiniteField, GF
from linear import LinearGroup_generic, LinearGroup_finite_field

def GU(n, R):
    if is_FiniteField(R):
        return GeneralUnitaryGroup_finite_field(n, R)
    else:
        return GeneralUnitaryGroup_generic(n, R)

class GeneralUnitaryGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        return "GU(%s, %s)"%(self.degree(), self.base_ring().order())

    def _latex_(self):
        """
        EXAMPLES:
            sage: G = GU(3,GF(5))
            sage: G._latex_()
            'GU$(3, 5)$'

        """
        return "GU$(%s, %s)$"%(self.degree(), self.base_ring().order())

    def __str__(self):
        """
        EXAMPLES:
            sage: G = GU(3,GF(5))
            sage: print G
            GU(3, GF(5))

        """
        return "GU(%s, GF(%s))"%(self.degree(), self.base_ring().order())

    def __repr__(self):
        """
        EXAMPLES:
            sage: G = GU(3,GF(5))
            sage: G
            General Unitary Group of degree 3 over Finite Field of size 5

        """
        return "General Unitary Group of degree %s over %s"%(self.degree(), self.base_ring())

    def gens(self):
        """
        EXAMPLES:
            sage: G = GU(4,GF(5))
            sage: G.gens()
            [[  a   0   0   0]
            [  0   1   0   0]
            [  0   0   1   0]
            [  0   0   0 3*a], [      1       0 4*a + 3       0]
                               [      1       0       0       0]
                               [      0 2*a + 4       0       1]
                               [      0 3*a + 1       0       0]]
        """
        from sage.interfaces.all import gap
        q = self.base_ring().order()
        F = GF(q**2)
        G = self._gap_init_()
        n = eval(gap.eval("Length(GeneratorsOfGroup(%s))"%G))
        gens = [gap("GeneratorsOfGroup(%s)[%s]"%(G,i))._matrix_(F) for i in range(1,n+1)]
        return gens

    def as_matrix_group(self):
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        gns = self.gens()
        G = MatrixGroup(gns)
        return G

class GeneralUnitaryGroup_finite_field(GeneralUnitaryGroup_generic, LinearGroup_finite_field):
    pass

def SU(n, R):
    if is_FiniteField(R):
        return SpecialUnitaryGroup_finite_field(n, R)
    else:
        return SpecialUnitaryGroup_generic(n, R)

class SpecialUnitaryGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        return "SU(%s, %s)"%(self.degree(), self.base_ring().order())

    def _latex_(self):
        """
        EXAMPLES:
            sage: G = SU(3,GF(5))
            sage: G._latex_()
            'SU$(3, 5)$'

        """
        return "SU$(%s, %s)$"%(self.degree(), self.base_ring().order())

    def __str__(self):
        """
        EXAMPLES:
            sage: G = SU(3,GF(5))
            sage: print G
            SU(3, GF(5))
        """
        return "SU(%s, GF(%s))"%(self.degree(), self.base_ring().order())

    def __repr__(self):
        """
        EXAMPLES:
            sage: G = SU(3,GF(5))
            sage: G
            Special Unitary Group of degree 3 over Finite Field of size 5

        """
        return "Special Unitary Group of degree %s over %s"%(self.degree(), self.base_ring())

    def gens(self):
        """
        EXAMPLES:
            sage: G = SU(4,GF(5))
            sage: G.gens()
            [[      a       0       0       0]
             [      0 2*a + 3       0       0]
             [      0       0 4*a + 1       0]
             [      0       0       0     3*a],
             [      1       0 4*a + 3       0]
             [      1       0       0       0]
             [      0 2*a + 4       0       1]
             [      0 3*a + 1       0       0]]
        """
        from sage.interfaces.all import gap
        q = self.base_ring().order()
        F = GF(q**2)
        G = self._gap_init_()
        n = eval(gap.eval("Length(GeneratorsOfGroup(%s))"%G))
        gens = [gap("GeneratorsOfGroup(%s)[%s]"%(G,i))._matrix_(F) for i in range(1,n+1)]
        return gens

    def as_matrix_group(self):
        """
        EXAMPLES:
        sage: G = SU(4,GF(5))
        sage: G.as_matrix_group()
        Matrix group over Finite Field in a of size 5^2 with 2 generators:
        [[[a, 0, 0, 0], [0, 2*a + 3, 0, 0], [0, 0, 4*a + 1, 0], [0, 0, 0, 3*a]], [[1, 0, 4*a + 3, 0], [1, 0, 0, 0], [0, 2*a + 4, 0, 1], [0, 3*a + 1, 0, 0]]]
        """
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        gns = self.gens()
        G = MatrixGroup(gns)
        return G

class SpecialUnitaryGroup_finite_field(SpecialUnitaryGroup_generic, LinearGroup_finite_field):
    pass


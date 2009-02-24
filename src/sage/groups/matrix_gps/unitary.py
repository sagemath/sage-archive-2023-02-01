r"""
Unitary Groups `GU(n,q)` and `SU(n,q)`

These are `n \times n` unitary matrices with entries in
`GF(q^2)`.

AUTHORS:

- David Joyner (2006-03): initial version, modified from
  special_linear (by W. Stein)

- David Joyner (2006-05): minor additions (examples, _latex_, __str__,
  gens)

- William Stein (2006-12): rewrite

EXAMPLES::

    sage: G = SU(3,GF(5))
    sage: G.order()
    378000
    sage: G
    Special Unitary Group of degree 3 over Finite Field of size 5
    sage: G._gap_init_()
    'SU(3, 5)'
    sage: G.random_element()
    [      1 4*a + 4 4*a + 1]
    [2*a + 4 2*a + 1       0]
    [      4     3*a 4*a + 2]
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

from sage.rings.all import IntegerRing, is_FiniteField, GF, Integer
from sage.interfaces.all import gap
from matrix_group import MatrixGroup_gap, MatrixGroup_gap_finite_field

###############################################################################
# General Unitary Group
###############################################################################

def GU(n, F, var='a'):
    """
    Return the general unitary group of degree n over the finite field
    F.

    INPUT:


    -  ``n`` - a positive integer

    -  ``F`` - finite field

    -  ``var`` - variable used to represent generator of
       quadratic extension of F, if needed.


    EXAMPLES::

        sage: G = GU(3,GF(7)); G
        General Unitary Group of degree 3 over Finite Field of size 7
        sage: G.gens()
        [
        [  a   0   0]
        [  0   1   0]
        [  0   0 5*a],
        [6*a   6   1]
        [  6   6   0]
        [  1   0   0]
        ]
        sage: G = GU(2,QQ)
        Traceback (most recent call last):
        ...
        NotImplementedError: general unitary group only implemented over finite fields

    ::

        sage: G = GU(3,GF(5), var='beta')
        sage: G.gens()
        [
        [  beta      0      0]
        [     0      1      0]
        [     0      0 3*beta],
        [4*beta      4      1]
        [     4      4      0]
        [     1      0      0]
        ]
    """
    if isinstance(F, (int, long, Integer)):
        F = GF(F,var)
    if is_FiniteField(F):
        return GeneralUnitaryGroup_finite_field(n, F, var)
    else:
        raise NotImplementedError, "general unitary group only implemented over finite fields"

class UnitaryGroup_finite_field(MatrixGroup_gap_finite_field):
    def field_of_definition(self):
        """
        Return the field of definition of this general unity group.

        EXAMPLES::

            sage: G = GU(3,GF(5))
            sage: G.field_of_definition()
            Finite Field in a of size 5^2
            sage: G.base_field()
            Finite Field of size 5
        """
        try:
            return self._field_of_definition
        except AttributeError:
            if self.base_ring().degree() % 2 == 0:
                k = self.base_ring()
            else:
                k = GF(self.base_ring().order()**2, names=self._var)
            self._field_of_definition = k
            return k

class GeneralUnitaryGroup_finite_field(UnitaryGroup_finite_field):
    def _gap_init_(self):
        """
        Return string that evaluates to creates this group as an element of
        GAP.

        EXAMPLES::

            sage: G = GU(3,GF(7)); G
            General Unitary Group of degree 3 over Finite Field of size 7
            sage: G._gap_init_()
            'GU(3, 7)'
            sage: gap(G._gap_init_())
            GU(3,7)
        """
        return "GU(%s, %s)"%(self.degree(), self.base_ring().order())

    def _latex_(self):
        r"""
        Return LaTeX string representation of this group.

        EXAMPLES::

            sage: G = GU(3,GF(7)); G
            General Unitary Group of degree 3 over Finite Field of size 7
            sage: latex(G)
            \text{GU}_{3}(\mathbf{F}_{7^{2}})
        """
        return "\\text{GU}_{%s}(%s)"%(self.degree(), self.field_of_definition()._latex_())

    def _repr_(self):
        """
        Return text representatin of self.

        EXAMPLES::

            sage: G = GU(3,GF(5))
            sage: G
            General Unitary Group of degree 3 over Finite Field of size 5
        """
        return "General Unitary Group of degree %s over %s"%(self.degree(), self.base_ring())


###############################################################################
# Special Unitary Group
###############################################################################

def SU(n, F, var='a'):
    """
    Return the special unitary group of degree `n` over
    `F`.

    EXAMPLES::

        sage: SU(3,5)
        Special Unitary Group of degree 3 over Finite Field of size 5
        sage: SU(3,QQ)
        Traceback (most recent call last):
        ...
        NotImplementedError: special unitary group only implemented over finite fields
    """
    if isinstance(F, (int, long, Integer)):
        F = GF(F,var)
    if is_FiniteField(F):
        return SpecialUnitaryGroup_finite_field(n, F, var=var)
    else:
        raise NotImplementedError, "special unitary group only implemented over finite fields"

class SpecialUnitaryGroup_finite_field(UnitaryGroup_finite_field):

    def _gap_init_(self):
        """
        Return string that creates this group in GAP.

        EXAMPLES::

            sage: SU(3,5)._gap_init_()
            'SU(3, 5)'
        """
        return "SU(%s, %s)"%(self.degree(), self.base_ring().order())

    def _latex_(self):
        """
        Return latex representatin of this group.

        EXAMPLES::

            sage: G = SU(3,GF(5))
            sage: latex(G)
            \text{SU}_{3}(\mathbf{F}_{5^{2}})
        """
        return "\\text{SU}_{%s}(%s)"%(self.degree(), self.field_of_definition()._latex_())


    def _repr_(self):
        """
        Return text representation of this special unitary group.

        EXAMPLES::

            sage: G = SU(3,GF(5))
            sage: G
            Special Unitary Group of degree 3 over Finite Field of size 5
        """
        return "Special Unitary Group of degree %s over %s"%(self.degree(), self.base_ring())






"""
Orthogonal Linear Groups


Paraphrased from the GAP manual:
The general orthogonal group GO(e,d,q) consists of those dxd matrices
over the field GF(q) that respect a non-singular quadratic
form specified by e. (Use the GAP command InvariantQuadraticForm to
determine this form explicitly.) The value of e must be 0 for odd d
(and can optionally be omitted in this case), respectively one of
1 or -1 for even d.

SpecialOrthogonalGroup returns a group isomorphic to the special
orthogonal group SO(e,d,q), which is the subgroup of all those
matrices in the general orthogonal group that have determinant one.
(The index of SO(e,d,q) in GO(e,d,q) is 2 if q is odd,
but SO(e,d,q) = GO(e,d,q) if q is even.)

AUTHOR:
    -- David Joyner: initial version (2006-3)


"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.modular.all import SL2Z
from sage.rings.all import IntegerRing, is_FiniteField
from linear import LinearGroup_generic, LinearGroup_finite_field

def SO(d, R, e=0):
    if d%2!=0 and e!=0:
        raise ValueError, "\n Must have e = 0 for d even.\n"
    if d%2 == 0 and e**2!=1:
            raise ValueError, "\n Must have e=-1 or e=1 if d is even.\n"
    if is_FiniteField(R):
        return SpecialOrthogonalGroup_finite_field(d, R, e)
    else:
        return SpecialOrthogonalGroup_generic(d, R, e)

class SpecialOrthogonalGroup_generic(LinearGroup_generic):
    """

    EXAMPLES:
        sage: G = SO(4,GF(7),1)
        sage: G
        Special Orthogonal Group of degree 4, form parameter 1, over the Finite Field of size 7
        sage: G._gap_init_()
        'SO(1, 4, 7)'
        sage: G.random()
        [4 2 5 6]
	[0 3 2 4]
	[5 3 5 2]
	[1 1 6 2]

    """
    def _gap_init_(self):
        return "SO(%s, %s, %s)"%( self.invariant_form(), self.degree(), self.base_ring().order())

    def _repr_(self):
        return "Special Orthogonal Group of degree %s, form parameter %s, over the %s"%( self.degree(), self.invariant_form(), self.base_ring())

    def invariant_quadratic_form(self):
        """
        This wraps GAP's command "InvariantQuadraticForm". From the GAP documentation:

        INPUT:
           self -- a matrix group G
        OUTPUT:
           Q -- the matrix satisfying the property: The quadratic form q on the natural vector space V on
                which G acts is given by $q(v) = v Q v^t$, and the invariance under G is given by the equation
                $q(v) = q(v M)$ for all $v \in V$ and $M \in G$.

        EXAMPLES:
            sage: G = SO(4,GF(7),1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 3 0]
            [0 0 0 1]
            """
        from sage.interfaces.gap import gap
        F = self.base_ring()
        G = self._gap_init_()
        cmd = "r := InvariantQuadraticForm("+G+")"
        gap.eval(cmd)
        cmd = "r.matrix"
        Q = gap(cmd)
        return Q._matrix_(F)

class SpecialOrthogonalGroup_finite_field(SpecialOrthogonalGroup_generic, LinearGroup_finite_field):
    pass

def GO(d, R, e=0):
    if d%2!=0 and e!=0:
        raise ValueError, "\n if e = 0 then d must be even.\n"
    if d%2 == 0 and e**2!=1:
            raise ValueError, "\n Must have e=-1 or e=1 if d is even.\n"
    if is_FiniteField(R):
        return GeneralOrthogonalGroup_finite_field(d, R, e)
    else:
        return GeneralOrthogonalGroup_generic(d, R, e)

class GeneralOrthogonalGroup_generic(LinearGroup_generic):
    """
    EXAMPLES:
        sage: GO(3,GF(7),0)
        General Orthogonal Group of degree 3, form parameter 0, over the Finite Field of size 7
        sage: GO(3,GF(7),0).order()
        672
        sage: GO(3,GF(7),0).random()
        [1 6 6]
        [3 2 6]
        [3 6 5]

    """
    def _gap_init_(self):
        return "GO(%s, %s, %s)"%(  self.invariant_form(), self.degree(), (self.base_ring()).order() )

    def _repr_(self):
        return "General Orthogonal Group of degree %s, form parameter %s, over the %s"%( self.degree(), self.invariant_form(), self.base_ring())

    def invariant_quadratic_form(self):
        """
        This wraps GAP's command "InvariantQuadraticForm". From the GAP documentation:

        INPUT:
           self -- a matrix group G
        OUTPUT:
           Q -- the matrix satisfying the property: The quadratic form q on the natural vector space V on
                which G acts is given by $q(v) = v Q v^t$, and the invariance under G is given by the equation
                $q(v) = q(v M)$ for all $v \in V$ and $M \in G$.

        EXAMPLES:
            sage: G = GO(4,GF(7),1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 3 0]
            [0 0 0 1]
        """
        from sage.interfaces.gap import gap
        F = self.base_ring()
        G = self._gap_init_()
        cmd = "r := InvariantQuadraticForm("+G+")"
        gap.eval(cmd)
        cmd = "r.matrix"
        Q = gap(cmd)
        return Q._matrix_(F)

class GeneralOrthogonalGroup_finite_field(GeneralOrthogonalGroup_generic, LinearGroup_finite_field):
    pass


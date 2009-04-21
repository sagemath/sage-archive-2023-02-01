r"""
Orthogonal Linear Groups

Paraphrased from the GAP manual: The general orthogonal group
`GO(e,d,q)` consists of those `d\times d` matrices
over the field `GF(q)` that respect a non-singular
quadratic form specified by `e`. (Use the GAP command
InvariantQuadraticForm to determine this form explicitly.) The
value of `e` must be `0` for odd `d` (and
can optionally be omitted in this case), respectively one of
`1` or `-1` for even `d`.

SpecialOrthogonalGroup returns a group isomorphic to the special
orthogonal group `SO(e,d,q)`, which is the subgroup of all
those matrices in the general orthogonal group that have
determinant one. (The index of `SO(e,d,q)` in
`GO(e,d,q)` is `2` if `q` is odd, but
`SO(e,d,q) = GO(e,d,q)` if `q` is even.)

.. warning::

   GAP notation: GO([e,] d, q), SO([e,] d, q) ([...] denotes and
   optional value)

Sage notation: GO(d, GF(q), e=0), SO( d, GF(q), e=0)

There is no Python trick I know of to allow the first argument to
have the default value e=0 and leave the other two arguments as
non-default. This forces us into non-standard notation.

AUTHORS:

- David Joyner (2006-03): initial version

- David Joyner (2006-05): added examples, _latex_, __str__, gens,
  as_matrix_group

- William Stein (2006-12-09): rewrite
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import IntegerRing, is_FiniteField, GF, Integer, FiniteField
from matrix_group import MatrixGroup_gap, MatrixGroup_gap_finite_field

from sage.interfaces.gap import gap

def SO(n, R, e=0, var='a'):
    """
    Return the special orthogonal group of degree `n` over the
    ring `R`.

    INPUT:


    -  ``n`` - the degree

    -  ``R`` - ring

    -  ``e`` - a parameter for orthogonal groups only
       depending on the invariant form


    EXAMPLES::

        sage: G = SO(3,GF(5))
        sage: G.gens()
        [
        [2 0 0]
        [0 3 0]
        [0 0 1],
        [3 2 3]
        [0 2 0]
        [0 3 1],
        [1 4 4]
        [4 0 0]
        [2 0 4]
        ]
        sage: G = SO(3,GF(5))
        sage: G.as_matrix_group()
        Matrix group over Finite Field of size 5 with 3 generators:
        [[[2, 0, 0], [0, 3, 0], [0, 0, 1]], [[3, 2, 3], [0, 2, 0], [0, 3, 1]], [[1, 4, 4], [4, 0, 0], [2, 0, 4]]]
    """
    if isinstance(R, (int, long, Integer)):
        R = FiniteField(R, var)
    if n%2!=0 and e != 0:
        raise ValueError, "must have e = 0 for n even"
    if n%2 == 0 and e**2 != 1:
            raise ValueError, "must have e=-1 or e=1 if n is even"
    if is_FiniteField(R):
        return SpecialOrthogonalGroup_finite_field(n, R, e)
    else:
        return SpecialOrthogonalGroup_generic(n, R, e)

class OrthogonalGroup(MatrixGroup_gap):
    def __init__(self, n, R, e=0, var='a'):
        """
        INPUT:


        -  ``n`` - the degree

        -  ``R`` - the base ring

        -  ``e`` - a parameter for orthogonal groups only
           depending on the invariant form

        -  ``var`` - variable used to define field of
           definition of actual matrices in this group.
        """
        MatrixGroup_gap.__init__(self, n, R, var)
        self.__form = e

    def invariant_form(self):
        """
        Return the invariant form of this orthogonal group.

        TODO: What is the point of this? What does it do? How does it
        work?

        EXAMPLES::

            sage: G = SO( 4, GF(7), 1)
            sage: G.invariant_form()
            1
        """
        return self.__form

class SpecialOrthogonalGroup_generic(OrthogonalGroup):
    """
    EXAMPLES::

        sage: G = SO( 4, GF(7), 1); G
        Special Orthogonal Group of degree 4, form parameter 1, over the Finite Field of size 7
        sage: G._gap_init_()
        'SO(1, 4, 7)'
        sage: G.random_element()
        [1 2 5 0]
        [2 2 1 0]
        [1 3 1 5]
        [1 3 1 3]
    """
    def _gap_init_(self):
        """
        EXAMPLES::

            sage: G = SO(3,GF(5))
            sage: G._gap_init_()
            'SO(0, 3, 5)'
        """
        return "SO(%s, %s, %s)"%(  self.invariant_form(), self.degree(), self.base_ring().order())

    def _repr_(self):
        """
        EXAMPLES::

            sage: G = SO(3,GF(5))
            sage: G
            Special Orthogonal Group of degree 3, form parameter 0, over the Finite Field of size 5
        """
        return "Special Orthogonal Group of degree %s, form parameter %s, over the %s"%( self.degree(), self.invariant_form(), self.base_ring())

    def _latex_(self):
        """
        EXAMPLES::

            sage: G = SO(3,GF(5))
            sage: latex(G)
            \text{SO}_{3}(\Bold{F}_{5}, 0)
        """
        return "\\text{SO}_{%s}(%s, %s)"%(self.degree(), self.base_ring()._latex_(), self.invariant_form())

    def invariant_quadratic_form(self):
        r"""
        Return the quadratic form `q(v) = v Q v^t` on the space on
        which this group `G` that satisfies the equation
        `q(v) = q(v M)` for all `v \in V` and
        `M \in G`.

        .. note::

           Uses GAP's command InvariantQuadraticForm.

        OUTPUT:


        -  ``Q`` - matrix that defines the invariant quadratic
           form.


        EXAMPLES::

            sage: G = SO( 4, GF(7), 1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 3 0]
            [0 0 0 1]
        """
        F = self.base_ring()
        I = gap(self).InvariantQuadraticForm()
        Q = I.attribute('matrix')
        return Q._matrix_(F)


class SpecialOrthogonalGroup_finite_field(SpecialOrthogonalGroup_generic, MatrixGroup_gap_finite_field):
    pass


########################################################################
# General Orthogonal Group
########################################################################

def GO( n , R , e=0 ):
    """
    Return the general orthogonal group.

    EXAMPLES:
    """
    if n%2!=0 and e!=0:
        raise ValueError, "if e = 0, then n must be even."
    if n%2 == 0 and e**2!=1:
        raise ValueError, "must have e=-1 or e=1, if d is even."
    if isinstance(R, (int, long, Integer)):
        R = FiniteField(R)
    if is_FiniteField(R):
        return GeneralOrthogonalGroup_finite_field(n, R, e)
    else:
        return GeneralOrthogonalGroup_generic(n, R, e)

class GeneralOrthogonalGroup_generic(OrthogonalGroup):
    """
    EXAMPLES::

        sage: GO( 3, GF(7), 0)
        General Orthogonal Group of degree 3, form parameter 0, over the Finite Field of size 7
        sage: GO( 3, GF(7), 0).order()
        672
        sage: GO( 3, GF(7), 0).random_element()
        [5 1 4]
        [1 0 0]
        [6 0 1]
    """
    def _gap_init_(self):
        """
        EXAMPLES::

            sage: GO( 3, GF(7), 0)._gap_init_()
            'GO(0, 3, 7)'
        """
        return "GO(%s, %s, %s)"%( self.invariant_form(), self.degree(), (self.base_ring()).order() )

    def _repr_(self):
        """
        String representation of self.

        EXAMPLES::

            sage: GO(3,7)
            General Orthogonal Group of degree 3, form parameter 0, over the Finite Field of size 7
        """
        return "General Orthogonal Group of degree %s, form parameter %s, over the %s"%( self.degree(), self.invariant_form(), self.base_ring())

    def _latex_(self):
        """
        EXAMPLES::

            sage: G = GO(3,GF(5))
            sage: latex(G)
            \text{GO}_{3}(5, 0)
        """
        return "\\text{GO}_{%s}(%s, %s)"%( self.degree(), self.base_ring().order(),self.invariant_form() )

    def invariant_quadratic_form(self):
        """
        This wraps GAP's command "InvariantQuadraticForm". From the GAP
        documentation:

        INPUT:


        -  ``self`` - a matrix group G


        OUTPUT:


        -  ``Q`` - the matrix satisfying the property: The
           quadratic form q on the natural vector space V on which G acts is
           given by `q(v) = v Q v^t`, and the invariance under G is
           given by the equation `q(v) = q(v M)` for all
           `v \in V` and `M \in G`.


        EXAMPLES::

            sage: G = GO( 4, GF(7), 1)
            sage: G.invariant_quadratic_form()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 3 0]
            [0 0 0 1]
        """
        F = self.base_ring()
        G = self._gap_init_()
        cmd = "r := InvariantQuadraticForm("+G+")"
        gap.eval(cmd)
        cmd = "r.matrix"
        Q = gap(cmd)
        return Q._matrix_(F)

class GeneralOrthogonalGroup_finite_field(GeneralOrthogonalGroup_generic, MatrixGroup_gap_finite_field):
    pass


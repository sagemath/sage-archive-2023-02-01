r"""
Hecke triangle group elements

AUTHORS:

- Jonas Jermann (2014): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.matrix_gps.group_element import MatrixGroupElement_generic


class HeckeTriangleGroupElement(MatrixGroupElement_generic):
    r"""
    Elements of HeckeTriangleGroup.
    """

    def __init__(self, parent, M, check=True, **kwargs):
        r"""
        An element of HeckeTriangle group given by a matrix ``M``.

        INPUT:

        - ``parent``  -- A ``HeckeTriangleGroup``.

        - ``M``       -- A matrix which coerces into the matrix space
                         of ``parent``. For example with entries in a
                         polynomial ring over ``ZZ`` with parameter ``lam``.

        - ``check``   -- ``True`` (default) or ``False``. If ``True``
                         then a (possibly long) check is performed
                         whether ``M`` really corresponds to a group
                         element of ``parent``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup, HeckeTriangleGroupElement
            sage: lam = PolynomialRing(ZZ, 'lam').gen()
            sage: M = matrix([[-1, 0], [1-lam+lam^5, -1]])
            sage: G = HeckeTriangleGroup(10)
            sage: HeckeTriangleGroupElement(G, M)
            [                 -1                   0]
            [5*lam^3 - 6*lam + 1                  -1]
            sage: G = HeckeTriangleGroup(4)
            sage: G(M)
            [       -1         0]
            [3*lam + 1        -1]
            sage: G(M).matrix().parent()
            Full MatrixSpace of 2 by 2 dense matrices over Order in Number Field in lam with defining polynomial x^2 - 2

            sage: M = matrix([[-1, lam], [0, 1]])
            sage: HeckeTriangleGroupElement(G, M)
            Traceback (most recent call last):
            ...
            ValueError: The matrix is not an element of Hecke triangle group for n = 4, it has determinant -1 != 1.

            sage: G.T().inverse()
            [   1 -lam]
            [   0    1]
            sage: G.U() == G.T()*G.S()
            True
            sage: G.U()^(-4) == -G.I()
            True
        """

        MatrixGroupElement_generic.__init__(self, parent, M, check, convert=True)

        if (check):
            if self._matrix.determinant() != 1:
                raise ValueError("The matrix is not an element of {}, it has determinant {} != 1.".format(parent, self._matrix.determinant()))

            #TODO: check if ``self`` really is a group element...

    def __neg__(self):
        r"""
        Return the group element corresponding to the negative of the underlying matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: -U
            [-lam    1]
            [  -1    0]
        """

        return self.parent().element_class(self.parent(), -self._matrix, check=False)

    def __getitem__(self, key):
        r"""
        Return the corresponding rows/entries of the underlying matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U[0]
            (lam, -1)
            sage: U[0].parent()
            Ambient free module of rank 2 over the principal ideal domain Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
            sage: U[1][0]
            1
            sage: U[1][0].parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        return self._matrix.__getitem__(key)

    def a(self):
        r"""
        Return the upper left entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.a()
            lam
            sage: U.a().parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        return self._matrix[0][0]

    def b(self):
        r"""
        Return the upper right entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.b()
            -1
            sage: U.b().parent()
            Maximal Order in Number Field in lam with defining polynomial x^3 - x^2 - 2*x + 1
        """

        return self._matrix[0][1]

    def c(self):
        r"""
        Return the lower left entry of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.c()
            1
        """

        return self._matrix[1][0]

    def d(self):
        r"""
        Return the lower right of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: U = HeckeTriangleGroup(n=7).U()
            sage: U.d()
            0
        """

        return self._matrix[1][1]

    def trace(self):
        r"""
        Return the trace of ``self``, which is the sum of the diagonal entries.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.U().trace()
            lam
            sage: G.S().trace()
            0
        """

        return self._matrix.trace()

    def discriminant(self):
        r"""
        Return the discriminant of ``self`` which corresponds to
        the discriminant of the corresponding quadratic form of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: G.V(3).discriminant()
            4*lam^2 + 4*lam - 4
            sage: AA(G.V(3).discriminant())
            16.19566935808922?
        """

        return self.trace()**2 - 4

    def is_hyperbolic(self):
        r"""
        Return whether ``self`` is a hyperbolic matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_hyperbolic() for k in range(1,8) ]
            [False, True, True, True, True, False, False]
            sage: G.U().is_hyperbolic()
            False
        """

        return self.trace().abs() > 2

    def is_parabolic(self, exclude_one=False):
        r"""
        Return whether ``self`` is a parabolic matrix.

        If ``exclude_one`` is set, then +- the identity
        element is not considered parabolic.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_parabolic() for k in range(1,8) ]
            [True, False, False, False, False, True, False]
            sage: G.U().is_parabolic()
            False
            sage: G.V(6).is_parabolic(exclude_one=True)
            True
            sage: G.V(7).is_parabolic(exclude_one=True)
            False
        """

        if (exclude_one and self.is_identity()):
            return False

        return self.trace().abs() == 2

    def is_identity(self):
        r"""
        Return whether ``self`` is the identity or minus the identity.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_identity() for k in range(1,8) ]
            [False, False, False, False, False, False, False]
            sage: G.U().is_identity()
            False
        """

        if self == self.parent().I() or self == -self.parent().I():
            return True
        else:
            return False

    def is_elliptic(self):
        r"""
        Return whether ``self`` is an elliptic matrix.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(n=7)
            sage: [ G.V(k).is_elliptic() for k in range(1,8) ]
            [False, False, False, False, False, False, True]
            sage: G.U().is_elliptic()
            True
        """

        return self.trace().abs() < 2

    def acton(self, z):
        r"""
        Return the image of ``z`` under the action of ``self``
        by linear fractional transformations.

        INPUT:

        - ``z``     -- A complex number or an element of AlgebraicField().

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
            sage: G = HeckeTriangleGroup(5)
            sage: G.S().acton(1 + i/2)
            2/5*I - 4/5
            sage: G.S().acton(1 + i/2).parent()
            Symbolic Ring
            sage: G.S().acton(i + exp(-2))
            -1/(e^(-2) + I)
            sage: G.S().acton(i + exp(-2)).parent()
            Symbolic Ring
        """

        return (self.a()*z + self.b()) / (self.c()*z + self.d())

    # def _act_on_(self, other, self_on_left):
    #     TODO: implement default actions for "suitable" x
    #     if (self_on_left):
    #         return self.acton(other)

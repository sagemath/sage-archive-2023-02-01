# -*- coding: utf-8 -*-
"""
Braid groups

Braid groups are implemented as a particular case of finitely presented groups,
but with a lot of specific methods for braids.

A braid group can be created by giving the number of strands, and the name of the generators::

    sage: BraidGroup(3)
    Braid group on 3 strands
    sage: BraidGroup(3,'a')
    Braid group on 3 strands
    sage: BraidGroup(3,'a').gens()
    (a0, a1)
    sage: BraidGroup(3,'a,b').gens()
    (a, b)

The elements can be created by operating with the generators, or by passing a list
with the indices of the letters to the group::

    sage: B.<s0,s1,s2> = BraidGroup(4)
    sage: s0*s1*s0
    s0*s1*s0
    sage: B([1,2,1])
    s0*s1*s0

The mapping class action of the braid group over the free group is
also implemented, see :class:`MappingClassGroupAction` for an
explanation. This action is left multiplication of a free group
element by a braid::

    sage: B.<b0,b1,b2> = BraidGroup()
    sage: F.<f0,f1,f2,f3> = FreeGroup()
    sage: B.strands() == F.rank()   # necessary for the action to be defined
    True
    sage: f1 * b1
    f1*f2*f1^-1
    sage: f0 * b1
    f0
    sage: f1 * b1
    f1*f2*f1^-1
    sage: f1^-1 * b1
    f1*f2^-1*f1^-1

AUTHORS:

- Miguel Angel Marco Buzunariz
- Volker Braun
- Søren Fuglede Jørgensen
- Robert Lipshitz
- Thierry Monteil: add a ``__hash__`` method consistent with the word
  problem to ensure correct Cayley graph computations.
"""

##############################################################################
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

import six
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.groups.free_group import FreeGroup, is_FreeGroup
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.matrix.constructor import identity_matrix, matrix
from sage.combinat.permutation import Permutation
from sage.categories.action import Action
from sage.sets.set import Set
from sage.groups.finitely_presented import FinitelyPresentedGroup, FinitelyPresentedGroupElement


class Braid(FinitelyPresentedGroupElement):
    """
    Class that models elements of the braid group.

    It is a particular case of element of a finitely presented group.

    EXAMPLES::

        sage: B.<s0,s1,s2> = BraidGroup(4)
        sage: B
        Braid group on 4 strands
        sage: s0*s1/s2/s1
        s0*s1*s2^-1*s1^-1
        sage: B((1, 2, -3, -2))
        s0*s1*s2^-1*s1^-1
    """
    def _cmp_(self, other):
        """
        Compare ``self`` and ``other``

        TESTS::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 1])
            sage: c = B([2, 1, 2])
            sage: b == c #indirect doctest
            True
            sage: b._cmp_(c^(-1))
            -1
            sage: B([])._cmp_(B.one())
            0
        """
        if self.Tietze()==other.Tietze():
            return 0
        nfself = [i.Tietze() for i in self.left_normal_form()]
        nfother = [i.Tietze() for i in other.left_normal_form()]
        return cmp(nfself, nfother)

    __cmp__ = _cmp_

    def __hash__(self):
        r"""
        Return a hash value for ``self``.

        EXAMPLES::

            sage: B.<s0,s1,s2> = BraidGroup(4)
            sage: hash(s0*s2) == hash(s2*s0)
            True
            sage: hash(s0*s1) == hash(s1*s0)
            False
        """
        return hash(tuple(i.Tietze() for i in self.left_normal_form()))

    def _latex_(self):
        """
        Return a LaTeX representation

        OUTPUT:

        String. A valid LaTeX math command sequence.

        TESTS::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 3, -1, 2, -3])
            sage: b._latex_()
            '\\sigma_{1}\\sigma_{2}\\sigma_{3}\\sigma_{1}^{-1}\\sigma_{2}\\sigma_{3}^{-1}'
        """
        latexrepr = ''
        for i in self.Tietze():
            if i>0:
                latexrepr = latexrepr+"\sigma_{%s}"%i
            if i<0:
                latexrepr = latexrepr+"\sigma_{%s}^{-1}"%(-i)
        return latexrepr

    def strands(self):
        """
        Return the number of strands in the braid.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, -1, 3, -2])
            sage: b.strands()
            4
        """
        return self.parent().strands()

    def exponent_sum(self):
        """
        Return the exponent sum of the braid.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: b = B([1, 4, -3, 2])
            sage: b.exponent_sum()
            2
            sage: b = B([])
            sage: b.exponent_sum()
            0
        """
        tietze = self.Tietze()
        return sum([s.sign() for s in tietze])

    def components_in_closure(self):
        """
        Return the number of components of the trace closure of the braid.

        OUTPUT:

        Positive integer.

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: b = B([1, -3])  # Three disjoint unknots
            sage: b.components_in_closure()
            3
            sage: b = B([1, 2, 3, 4])  # The unknot
            sage: b.components_in_closure()
            1
            sage: B = BraidGroup(4)
            sage: K11n42 = B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2])
            sage: K11n42.components_in_closure()
            1
        """
        n = self.strands()
        perm = self.permutation()
        cycles = perm.to_cycles(singletons=False)
        return n-sum([len(c)-1 for c in cycles])

    def burau_matrix(self, var='t', reduced=False):
        """
        Return the Burau matrix of the braid.

        INPUT:

        - ``var`` -- string (default: ``'t'``); the name of the
          variable in the entries of the matrix
        - ``reduced`` -- boolean (default: ``False``); whether to
          return the reduced or unreduced Burau representation

        OUTPUT:

        The Burau matrix of the braid. It is a matrix whose entries
        are Laurent polynomials in the variable ``var``. If ``reduced``
        is ``True``, return the matrix for the reduced Burau representation
        instead.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.inject_variables()
            Defining s0, s1, s2
            sage: b = s0*s1/s2/s1
            sage: b.burau_matrix()
            [       1 - t            0      t - t^2          t^2]
            [           1            0            0            0]
            [           0            0            1            0]
            [           0         t^-2 -t^-2 + t^-1    -t^-1 + 1]
            sage: s2.burau_matrix('x')
            [    1     0     0     0]
            [    0     1     0     0]
            [    0     0 1 - x     x]
            [    0     0     1     0]
            sage: s0.burau_matrix(reduced=True)
            [-t  0  0]
            [-t  1  0]
            [-t  0  1]

        REFERENCES:

        - :wikipedia:`Burau_representation`
        """
        R = LaurentPolynomialRing(IntegerRing(), var)
        t = R.gen()
        n = self.strands()
        if not reduced:
            M = identity_matrix(R, n)
            for i in self.Tietze():
                A = identity_matrix(R, n)
                if i > 0:
                    A[i-1, i-1] = 1-t
                    A[i, i] = 0
                    A[i, i-1] = 1
                    A[i-1, i] = t
                if i < 0:
                    A[-1-i, -1-i] = 0
                    A[-i, -i] = 1-t**(-1)
                    A[-1-i, -i] = 1
                    A[-i, -1-i] = t**(-1)
                M = M * A
        else:
            M = identity_matrix(R, n - 1)
            for j in self.Tietze():
                A = identity_matrix(R, n - 1)
                if j > 1:
                    i = j-1
                    A[i-1, i-1] = 1-t
                    A[i, i] = 0
                    A[i, i-1] = 1
                    A[i-1, i] = t
                if j < -1:
                    i = j+1
                    A[-1-i, -1-i] = 0
                    A[-i, -i] = 1-t**(-1)
                    A[-1-i, -i] = 1
                    A[-i, -1-i] = t**(-1)
                if j == 1:
                    for k in range(n - 1):
                        A[k,0] = -t
                if j == -1:
                    A[0,0] = -t**(-1)
                    for k in range(1, n - 1):
                        A[k,0] = -1
                M = M * A
        return M

    def alexander_polynomial(self, var='t', normalized=True):
        r"""
        Return the Alexander polynomial of the closure of the braid.

        INPUT:

        - ``var`` -- string (default: ``'t'``); the name of the
          variable in the entries of the matrix
        - ``normalized`` -- boolean (default: ``True``); whether to
          return the normalized Alexander polynomial

        OUTPUT:

        The Alexander polynomial of the braid closure of the braid.

        This is computed using the reduced Burau representation. The
        unnormalized Alexander polynomial is a Laurent polynomial,
        which is only well-defined up to multiplication by plus or
        minus times a power of `t`.

        We normalize the polynomial by dividing by the largest power
        of `t` and then if the resulting constant coefficient
        is negative, we multiply by `-1`.

        EXAMPLES:

        We first construct the trefoil::

            sage: B = BraidGroup(3)
            sage: b = B([1,2,1,2])
            sage: b.alexander_polynomial(normalized=False)
            1 - t + t^2
            sage: b.alexander_polynomial()
            t^-2 - t^-1 + 1

        Next we construct the figure 8 knot::

            sage: b = B([-1,2,-1,2])
            sage: b.alexander_polynomial(normalized=False)
            -t^-2 + 3*t^-1 - 1
            sage: b.alexander_polynomial()
            t^-2 - 3*t^-1 + 1

        Our last example is the Kinoshita-Terasaka knot::

            sage: B = BraidGroup(4)
            sage: b = B([1,1,1,3,3,2,-3,-1,-1,2,-1,-3,-2])
            sage: b.alexander_polynomial(normalized=False)
            -t^-1
            sage: b.alexander_polynomial()
            1

        REFERENCES:

        - :wikipedia:`Alexander_polynomial`
        """
        n = self.strands()
        p = (self.burau_matrix(reduced=True) - identity_matrix(n - 1)).det()
        K, t = LaurentPolynomialRing(IntegerRing(), var).objgen()
        if p == 0:
            return K.zero()
        qn = sum(t ** i for i in range(n))
        p //= qn
        if normalized:
            p *= t ** (-p.degree())
            if p.constant_coefficient() < 0:
                p = -p
        return p

    def permutation(self):
        """
        Return the permutation induced by the braid in its strands.

        OUTPUT:

        A permutation.

        EXAMPLES::

            sage: B.<s0,s1,s2> = BraidGroup()
            sage: b = s0*s1/s2/s1
            sage: b.permutation()
            [4, 1, 3, 2]
            sage: b.permutation().cycle_string()
            '(1,4,2)'
        """
        per = Permutation((()))
        for i in self.Tietze():
            j = abs(i)
            per = per*Permutation(((j, j+1)))
        return per

    def plot(self, color='rainbow', orientation='bottom-top', gap=0.05, aspect_ratio=1, axes=False, **kwds):
        """
        Plot the braid

        The following options are available:

        - ``color`` -- (default: ``'rainbow'``) the color of the
          strands. Possible values are:

            * ``'rainbow'``, uses :meth:`~sage.plot.colors.rainbow`
              according to the number of strands.

            * a valid color name for :meth:`~sage.plot.bezier_path`
              and :meth:`~sage.plot.line`. Used for all strands.

            * a list or a tuple of colors for each individual strand.

        - ``orientation`` -- (default: ``'bottom-top'``) determines how
          the braid is printed. The possible values are:

            * ``'bottom-top'``, the braid is printed from bottom to top

            * ``'top-bottom'``, the braid is printed from top to bottom

            * ``'left-right'``, the braid is printed from left to right

        - ``gap`` -- floating point number (default: 0.05). determines
          the size of the gap left when a strand goes under another.

        - ``aspect_ratio`` -- floating point number (default:
          ``1``). The aspect ratio.

        - ``**kwds`` -- other keyword options that are passed to
          :meth:`~sage.plot.bezier_path` and :meth:`~sage.plot.line`.

        EXAMPLES::

            sage: B = BraidGroup(4, 's')
            sage: b = B([1, 2, 3, 1, 2, 1])
            sage: b.plot()
            Graphics object consisting of 30 graphics primitives
            sage: b.plot(color=["red", "blue", "red", "blue"])
            Graphics object consisting of 30 graphics primitives

            sage: B.<s,t> = BraidGroup(3)
            sage: b = t^-1*s^2
            sage: b.plot(orientation="left-right", color="red")
            Graphics object consisting of 12 graphics primitives
        """
        from sage.plot.bezier_path import bezier_path
        from sage.plot.plot import Graphics, line
        from sage.plot.colors import rainbow
        if orientation=='top-bottom':
            orx = 0
            ory = -1
            nx = 1
            ny = 0
        elif orientation=='left-right':
            orx = 1
            ory = 0
            nx = 0
            ny = -1
        elif orientation=='bottom-top':
            orx = 0
            ory = 1
            nx = 1
            ny = 0
        else:
            raise ValueError('unknown value for "orientation"')
        n = self.strands()
        if isinstance(color, (list, tuple)):
            if len(color) != n:
                raise TypeError("color (=%s) must contain exactly %d colors" % (color, n))
            col = list(color)
        elif color == "rainbow":
            col = rainbow(n)
        else:
            col = [color]*n
        braid = self.Tietze()
        a = Graphics()
        op = gap
        for i, m in enumerate(braid):
            for j in range(n):
                if m==j+1:
                    a += bezier_path([[(j*nx+i*orx, i*ory+j*ny), (j*nx+orx*(i+0.25), j*ny+ory*(i+0.25)),
                                       (nx*(j+0.5)+orx*(i+0.5), ny*(j+0.5)+ory*(i+0.5))],
                                      [(nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)),
                                       (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col[j], **kwds)
                elif m==j:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j-0.5+4*op)+orx*(i+0.5-2*op), ny*(j-0.5+4*op)+ory*(i+0.5-2*op)),
                                       (nx*(j-0.5+2*op)+orx*(i+0.5-op), ny*(j-0.5+2*op)+ory*(i+0.5-op))]],
                                     color=col[j], **kwds)
                    a += bezier_path([[(nx*(j-0.5-2*op)+orx*(i+0.5+op), ny*(j-0.5-2*op)+ory*(i+0.5+op)),
                                       (nx*(j-0.5-4*op)+orx*(i+0.5+2*op), ny*(j-0.5-4*op)+ory*(i+0.5+2*op)),
                                       (nx*(j-1)+orx*(i+0.75), ny*(j-1)+ory*(i+0.75)),
                                       (nx*(j-1)+orx*(i+1), ny*(j-1)+ory*(i+1))]], color=col[j], **kwds)
                    col[j], col[j-1] = col[j-1], col[j]
                elif -m==j+1:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j+0.5-4*op)+orx*(i+0.5-2*op), ny*(j+0.5-4*op)+ory*(i+0.5-2*op)),
                                       (nx*(j+0.5-2*op)+orx*(i+0.5-op), ny*(j+0.5-2*op)+ory*(i+0.5-op))]],
                                     color=col[j], **kwds)
                    a += bezier_path([[(nx*(j+0.5+2*op)+orx*(i+0.5+op), ny*(j+0.5+2*op)+ory*(i+0.5+op)),
                                       (nx*(j+0.5+4*op)+orx*(i+0.5+2*op), ny*(j+0.5+4*op)+ory*(i+0.5+2*op)),
                                       (nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)),
                                       (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col[j], **kwds)
                elif -m==j:
                    a += bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)),
                                       (nx*(j-0.5)+orx*(i+0.5), ny*(j-0.5)+ory*(i+0.5))],
                                      [(nx*(j-1)+orx*(i+0.75), ny*(j-1)+ory*(i+0.75)),
                                       (nx*(j-1)+orx*(i+1), ny*(j-1)+ory*(i+1))]], color=col[j], **kwds)
                    col[j], col[j-1] = col[j-1], col[j]
                else:
                    a += line([(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+1), ny*j+ory*(i+1))], color=col[j], **kwds)
        a.set_aspect_ratio(aspect_ratio)
        a.axes(axes)
        return a

    def plot3d(self, color='rainbow'):
        """
        Plots the braid in 3d.

        The following option is available:

        - ``color`` -- (default: ``'rainbow'``) the color of the
          strands. Possible values are:

            * ``'rainbow'``, uses :meth:`~sage.plot.colors.rainbow`
              according to the number of strands.

            * a valid color name for :meth:`~sage.plot.plot3d.bezier3d`.
              Used for all strands.

            * a list or a tuple of colors for each individual strand.

        EXAMPLES::

            sage: B = BraidGroup(4, 's')
            sage: b = B([1, 2, 3, 1, 2, 1])
            sage: b.plot3d()
            Graphics3d Object
            sage: b.plot3d(color="red")
            Graphics3d Object
            sage: b.plot3d(color=["red", "blue", "red", "blue"])
            Graphics3d Object
        """
        from sage.plot.plot3d.shapes2 import bezier3d
        from sage.plot.colors import rainbow
        b = []
        n = self.strands()
        if isinstance(color, (list, tuple)):
            if len(color) != n:
                raise TypeError("color (=%s) must contain exactly %d colors" % (color, n))
            col = list(color)
        elif color == "rainbow":
            col = rainbow(n)
        else:
            col = [color]*n
        braid = self.Tietze()

        for i, m in enumerate(braid):
            for j in range(n):
                if m==j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j+0.5, i+0.5)],
                                       [(0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]], color=col[j]))
                elif -m==j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j+0.5, i+0.5)],
                                       [(-0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]], color=col[j]))
                elif m==j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j-0.5, i+0.5)],
                                       [(-0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]], color=col[j]))
                    col[j], col[j-1] = col[j-1], col[j]
                elif -m==j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j-0.5, i+0.5)],
                                       [(0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]], color=col[j]))
                    col[j], col[j-1] = col[j-1], col[j]
                else:
                    b.append(bezier3d([[(0, j, i), (0, j, i+1)]], color=col[j]))
        return sum(b)

    def LKB_matrix(self, variables='x,y'):
        """
        Return the Lawence-Krammer-Bigelow representation matrix.

        The matrix is expressed in the basis $\{e_{i, j} \mid 1\\leq i
        < j \leq n\}$, where the indices are ordered
        lexicographically.  It is a matrix whose entries are in the
        ring of Laurent polynomials on the given variables.  By
        default, the variables are ``'x'`` and ``'y'``.

        INPUT:

        - ``variables`` -- string (default: ``'x,y'``). A string
          containing the names of the variables, separated by a comma.

        OUTPUT:

        The matrix corresponding to the Lawence-Krammer-Bigelow representation of the braid.

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1, 2, 1])
            sage: b.LKB_matrix()
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: c = B([2, 1, 2])
            sage: c.LKB_matrix()
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]

        REFERENCES:

        .. [Bigelow] Bigelow, Stephen J. The Lawrence-Krammer representation. arXiv:math/0204057v1
        """
        return self.parent()._LKB_matrix_(self.Tietze(), variab=variables)

    def TL_matrix(self, drain_size, variab='A', ring=IntegerRing(), sparse=True):
        """
        Return the matrix representation of the Temperley--Lieb--Jones
        representation of the braid in a certain basis.
        
        The basis is given by non-intersecting pairings of $(n+d)$ points,
        where `n` is the number of strands, and `d` is given by
        ``drain_size``, and the pairings satisfy certain rules. See
        ``BraidGroup_class.TL_basis_with_drain()`` for details.

        We use the convention that the eigenvalues of the standard generators
        are `1` and `-A^4`, where `A` is the generator of the Laurent
        polynomial ring.

        When `d = n - 2` and the variables are picked appropriately, the
        resulting representation is equivalent to the reduced Burau
        representation.
        
        The parameter ``sparse`` can be set to False if it is expected that the
        result will not be sparse. We make no attempt at guessing this.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands (both
          included)
        - ``variab`` -- string (default: ``'A'``); the name of the
          variable in the entries of the matrices
        - ``ring`` -- ring (default: ``IntegerRing()``); the ring to which the
          coefficients of the polynomial entries belong
        - ``sparse`` -- boolean (default: ``True``); whether or not the result
          should be given as a sparse matrix

        OUTPUT:

        The matrix of the TL representation of the braid.


        EXAMPLES:

        Let us calculate a few examples for `B_4` with `d = 0`::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, -3])
            sage: b.TL_matrix(0)
            [-A^4 + 1   1/-A^2]
            [    -A^6        0]
            sage: 2*b.TL_matrix(0, ring=GF(2))
            [0 0]
            [0 0]
            sage: b = B([])
            sage: b.TL_matrix(0)
            [1 0]
            [0 1]

        Test of one of the relations in $B_8$::

            sage: B = BraidGroup(8)
            sage: d = 0
            sage: B([4,5,4]).TL_matrix(d) == B([5,4,5]).TL_matrix(d)
            True

        An element of the kernel of the Burau representation, following
        [Big99]_::

            sage: B = BraidGroup(6)
            sage: psi1 = B([4, -5, -2, 1])
            sage: psi2 = B([-4, 5, 5, 2, -1, -1])
            sage: w1 = psi1^(-1) * B([3]) * psi1
            sage: w2 = psi2^(-1) * B([3]) * psi2
            sage: (w1 * w2 * w1^(-1) * w2^(-1)).TL_matrix(4)
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]

        REFERENCES:

        .. [Big99] Stephen J. Bigelow. The Burau representation is not faithful
                   for n = 5. Geom. Topol., 3:397–404, 1999.
        .. [Jon] Vaughan Jones. The Jones Polynomial.
                 https://math.berkeley.edu/~vfr/jones.pdf
        """
        R = LaurentPolynomialRing(ring, variab)
        A = R.gens()[0]
        n = self.strands()
        d = drain_size
        B = BraidGroup(n)
        rep = B.create_TL_rep(d, variab, ring)
        M = identity_matrix(R, B.dim_of_TL_space(d), sparse=sparse)
        for i in self.Tietze():
            if i > 0:
                M = M*rep[i-1]
            if i < 0:
                # Note: calculating the inverse might seem a bit silly at first
                # glance since it's obtained simply by replacing A by A^(-1).
                # However, in my testing, doing this calculation actually turns
                # out to be faster than obtaining the inverse through lookups.
                M = M*rep[-i-1]**(-1)
        return M

    def tropical_coordinates(self):
        r"""
        Return the tropical coordinates of ``self`` in the braid group `B_n`.

        OUTPUT:

        - a list of `2n` tropical integers

        EXAMPLES::

            sage: B = BraidGroup(3)
            sage: b = B([1])
            sage: tc = b.tropical_coordinates(); tc
            [1, 0, 0, 2, 0, 1]
            sage: tc[0].parent()
            Tropical semiring over Integer Ring

            sage: b = B([-2, -2, -1, -1, 2, 2, 1, 1])
            sage: b.tropical_coordinates()
            [1, -19, -12, 9, 0, 13]

        REFERENCES:

        .. [Dynnikov07] I. Dynnikov and B. Wiest, On the complexity of braids,
           J. Europ. Math. Soc. 9 (2007)
        .. [Dehornoy] P. Dehornoy, Le probleme d'isotopie des tresses, in
           lecons de mathematiques d'aujourd'hui vol. 4
        """
        coord = [0, 1] * self.strands()
        for s in self.Tietze():
            k = 2*(abs(s)-1)
            x1, y1, x2, y2 = coord[k:k+4]
            if s > 0:
                sign = 1
                z = x1 - min(y1, 0) - x2 + max(y2, 0)
                coord[k+1] = y2 - max(z, 0)
                coord[k+3] = y1 + max(z, 0)
            else:
                sign = -1
                z = x1 + min(y1, 0) - x2 - max(y2, 0)
                coord[k+1] = y2 + min(z, 0)
                coord[k+3] = y1 - min(z, 0)

            coord[k] = x1 + sign*(max(y1, 0) + max(max(y2, 0) - sign*z, 0))
            coord[k+2] = x2 + sign*(min(y2, 0) + min(min(y1, 0) + sign*z, 0))

        from sage.rings.semirings.tropical_semiring import TropicalSemiring
        T = TropicalSemiring(IntegerRing())
        return [T(_) for _ in coord]

    def markov_trace(self, variab='A', ring=IntegerRing()):
        """
        Return the Markov trace of the braid.
        
        The normalisation is so that in the underlying braid group
        representation, the eigenvalues of the standard generators of the braid
        group are `1` and `-A^4`.

        INPUT:

        - ``variab`` -- string (default: ``'A'``); the name of the variable in
          the resulting Laurent polynomial
        - ``ring`` -- ring (default: ``IntegerRing()``); the ring to which the
          coefficients of the polynomial entries belong

        OUTPUT:

        Quotient of Laurent polynomials over ``ring`` in the variable
        ``variab``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, -3])
            sage: b.markov_trace().factor()
            A^4 * (A^4 + 1)^-3
        """
        R = LaurentPolynomialRing(ring, variab)
        A = R.gens()[0]
        delta = -A**2 - A**(-2)
        n = self.strands()
        drains = [d for d in range(n+1) if (n+d) % 2 == 0]

        def weighted_trace(d):
            quantum_integer = (A**(2*(d+1)) - A**(-2*(d+1)))/(A**2 - A**(-2))
            return quantum_integer * self.TL_matrix(d, variab, ring).trace()

        traces = [weighted_trace(d) for d in drains]
        return sum(traces) / (-delta)**n

    def jones_polynomial(self, skein_variable=False):
        """
        Return the Jones polynomial of the trace closure of the braid.
        
        The normalization is so that the unknot has Jones polynomial `1`. If
        ``skein_variable`` is True, give the result in terms of a variable
        ``'A'`` so that the result agrees with the conventions of [Lic]_ (which
        in particular differs slightly from the conventions used otherwise in
        this class). If ``skein_variable`` is False, return the result in terms
        of the variable ``'t'``, also used in [Lic]_.

        The computation uses the representation of the braid group on the
        Temperley--Lieb algebra.

        INPUT:

        - ``skein_variable`` -- boolean (default: ``False``); determines the
          variable of the resulting polynomial.

        OUTPUT:

        Jones polynomial of the closure as an element in the symbolic ring,
        in the variable ``'A'`` or ``'t'`` depending on the value
        of ``skein_variable``. Might have fractional powers if
        ``skein_variable`` is False and the closure of the braid is not a knot.

        EXAMPLES:

        The unknot::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: b.jones_polynomial(skein_variable=True)
            1

        Two unlinked unknots::

            sage: B = BraidGroup(2)
            sage: b = B([])
            sage: b.jones_polynomial(skein_variable=True)
            -A^2 - 1/A^2

        The Hopf link::

            sage: B = BraidGroup(2)
            sage: b = B([-1,-1])
            sage: b.jones_polynomial()
            -1/sqrt(t) - 1/t^(5/2)

        Two different representations of the trefoil and one of its mirror::

            sage: B = BraidGroup(2)
            sage: b = B([-1, -1, -1])
            sage: b.jones_polynomial(skein_variable=True)
            1/A^4 + 1/A^12 - 1/A^16
            sage: b.jones_polynomial()
            1/t + 1/t^3 - 1/t^4
            sage: B = BraidGroup(3)
            sage: b = B([-1, -2, -1, -2])
            sage: b.jones_polynomial(skein_variable=True)
            1/A^4 + 1/A^12 - 1/A^16
            sage: B = BraidGroup(3)
            sage: b = B([1, 2, 1, 2])
            sage: b.jones_polynomial(skein_variable=True)
            -A^16 + A^12 + A^4

        K11n42 (the mirror of the "Kinoshita-Terasaka" knot) and K11n34 (the
        mirror of the "Conway" knot)::

            sage: B = BraidGroup(4)
            sage: b11n42 = B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2])
            sage: b11n34 = B([1, 1, 2, -3, 2, -3, 1, -2, -2, -3, -3])
            sage: cmp(b11n42.jones_polynomial(), b11n34.jones_polynomial())
            0

        REFERENCES:

        .. [Lic] William B. Raymond Lickorish. An Introduction to Knot Theory,
                 volume 175 of Graduate Texts in Mathematics. Springer-Verlag,
                 New York, 1997. ISBN 0-387-98254-X
        """
        from sage.symbolic.ring import SR
        from sage.rings.integer_ring import ZZ
        variab = 'A'
        ring = IntegerRing()
        R = LaurentPolynomialRing(ring, variab)
        A = R.gens()[0]
        D = -A**2 - A**(-2)
        n = self.strands()
        exp_sum = self.exponent_sum()
        num_comp = self.components_in_closure()
        trace = self.markov_trace(variab, ring)
        jones_pol = (-1)**(num_comp-1) * (-D)**(n-1) * A**(2*exp_sum) * trace
        if skein_variable:
            output_var = SR('A')
            return jones_pol.subs(A=output_var).expand()
        else:
            output_var = SR('t')
            return jones_pol.subs(A=output_var**(ZZ(1)/ZZ(4))).expand()

    @cached_method
    def left_normal_form(self):
        """
        Return the left normal form of the braid.

        OUTPUT:

        A tuple of braid generators in the left normal form. The first
        element is a power of $\Delta$, and the rest are permutation
        braids.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: b = B([1, 2, 3, -1, 2, -3])
            sage: b.left_normal_form()
            (s0^-1*s1^-1*s2^-1*s0^-1*s1^-1*s0^-1, s0*s1*s2*s1*s0, s0*s2*s1)
            sage: c = B([1])
            sage: c.left_normal_form()
            (1, s0)
        """
        lnfp = self._left_normal_form_perm_()
        a = lnfp[0]
        l = lnfp[1:]
        n = self.strands()
        delta = Permutation([n-i for i in range(n)])
        P = self.parent()
        return tuple( [P._permutation_braid(delta).__pow__(a)] +
                      [P._permutation_braid(i) for i in l] )

    def _left_normal_form_perm_(self):
        """
        Return the left normal form of the braid, in permutation form.

        OUTPUT:

        A tuple whose first element is the power of $\Delta$, and the
        rest are the permutations corresponding to the simple factors.

        EXAMPLES::

            sage: B = BraidGroup(12)
            sage: B([2, 2, 2, 3, 1, 2, 3, 2, 1, -2])._left_normal_form_perm_()
            (-1,
             [12, 11, 10, 9, 8, 7, 6, 5, 2, 4, 3, 1],
             [4, 1, 3, 2, 5, 6, 7, 8, 9, 10, 11, 12],
             [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12],
             [3, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12],
             [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12])
            sage: C=BraidGroup(6)
            sage: C([2, 3, -4, 2, 3, -5, 1, -2, 3, 4, 1, -2])._left_normal_form_perm_()
            (-2, [3, 5, 4, 2, 6, 1], [1, 6, 3, 5, 2, 4], [5, 6, 2, 4, 1, 3],
             [3, 2, 4, 1, 5, 6], [1, 5, 2, 3, 4, 6])
        """
        n = self.parent().strands()
        delta = 0
        Delta = Permutation([n-i for i in range(n)])
        l = self.Tietze()
        if l==():
            return (0,)
        form = []
        for i in l:
            if i>0:
                form.append(Permutation((i, i+1)))
            else:
                delta = delta+1
                form = [Delta*a*Delta for a in form]
                form.append(Delta*Permutation((-i, -i+1)))
        i = j = 0
        while j<len(form):
            while i<len(form)-j-1:
                e = form[i].inverse().descents()
                s = form[i+1].descents()
                S = set(s).difference(set(e))
                while S!=set([]):
                    a = list(S)[0]
                    form[i] = form[i]*Permutation((a+1, a+2))
                    form[i+1] = Permutation((a+1, a+2))*form[i+1]
                    e = form[i].inverse().descents()
                    s = form[i+1].descents()
                    S = set(s).difference(set(e))
                if form[i+1].length()==0:
                    form.pop(i+1)
                    i = 0
                else:
                    i += 1
            j += 1
            i = 0
        form = [a for a in form if a.length()>0]
        while form!=[] and form[0]==Delta:
            form.pop(0)
            delta = delta-1
        return tuple([-delta]+form)


class BraidGroup_class(FinitelyPresentedGroup):
    """
    The braid group on `n` strands.

    EXAMPLES::

        sage: B1 = BraidGroup(5)
        sage: B1
        Braid group on 5 strands
        sage: B2 = BraidGroup(3)
        sage: B1==B2
        False
        sage: B2 is BraidGroup(3)
        True
    """
    Element = Braid

    def __init__(self, names):
        """
        Python constructor.

        INPUT:

        - ``names`` -- a tuple of strings. The names of the
          generators.

        TESTS::

            sage: B1 = BraidGroup(5) # indirect doctest
            sage: B1
            Braid group on 5 strands
            sage: TestSuite(B1).run()


        Check that :trac:`14081` is fixed::

            sage: BraidGroup(2)
            Braid group on 2 strands
            sage: BraidGroup(('a',))
            Braid group on 2 strands

        Check that :trac:`15505` is fixed::

            sage: B=BraidGroup(4)
            sage: B.relations()
            (s0*s1*s0*s1^-1*s0^-1*s1^-1, s0*s2*s0^-1*s2^-1, s1*s2*s1*s2^-1*s1^-1*s2^-1)
            sage: B=BraidGroup('a,b,c,d,e,f')
            sage: B.relations()
            (a*b*a*b^-1*a^-1*b^-1,
             a*c*a^-1*c^-1,
             a*d*a^-1*d^-1,
             a*e*a^-1*e^-1,
             a*f*a^-1*f^-1,
             b*c*b*c^-1*b^-1*c^-1,
             b*d*b^-1*d^-1,
             b*e*b^-1*e^-1,
             b*f*b^-1*f^-1,
             c*d*c*d^-1*c^-1*d^-1,
             c*e*c^-1*e^-1,
             c*f*c^-1*f^-1,
             d*e*d*e^-1*d^-1*e^-1,
             d*f*d^-1*f^-1,
             e*f*e*f^-1*e^-1*f^-1)
        """
        n = len(names)
        if n<1: #n is the number of generators, not the number of strands (see ticket 14081)
            raise ValueError("the number of strands must be an integer bigger than one")
        free_group = FreeGroup(names)
        rels = []
        for i in range(1, n):
            rels.append(free_group([i, i+1, i, -i-1, -i, -i-1]))
            for j in range(i+2, n+1):
                rels.append(free_group([i, j, -i, -j]))
        FinitelyPresentedGroup.__init__(self, free_group, tuple(rels))
        self._nstrands_ = n+1

    def __reduce__(self):
        """
        TESTS::

            sage: B = BraidGroup(3)
            sage: B.__reduce__()
            (<class 'sage.groups.braid.BraidGroup_class'>, (('s0', 's1'),))
            sage: B = BraidGroup(3, 'sigma')
            sage: B.__reduce__()
            (<class 'sage.groups.braid.BraidGroup_class'>, (('sigma0', 'sigma1'),))
        """
        return (BraidGroup_class, (self.variable_names(), ))

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        TESTS::

            sage: B1 = BraidGroup(5)
            sage: B1 # indirect doctest
            Braid group on 5 strands
        """
        return "Braid group on %s strands"%self._nstrands_

    def cardinality(self):
        """
        Return the number of group elements.

        OUTPUT:

        Infinity.

        TESTS::

            sage: B1 = BraidGroup(5)
            sage: B1.cardinality()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    order = cardinality

    def as_permutation_group(self):
        """
        Return an isomorphic permutation group.

        OUTPUT:

        Raises a ``ValueError`` error since braid groups are infinite.

        TESTS::

            sage: B = BraidGroup(4, 'g')
            sage: B.as_permutation_group()
            Traceback (most recent call last):
            ...
            ValueError: the group is infinite
        """
        raise ValueError("the group is infinite")

    def strands(self):
        """
        Return the number of strands.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: B.strands()
            4
        """
        return self._nstrands_

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: B = BraidGroup(4)
            sage: B([1, 2, 3]) # indirect doctest
            s0*s1*s2
        """
        return self.element_class(self, x)

    def an_element(self):
        """
        Return an element of the braid group.

        This is used both for illustration and testing purposes.

        EXAMPLES::

            sage: B=BraidGroup(2)
            sage: B.an_element()
            s
        """
        return self.gen(0)

    def some_elements(self):
        """
        Return a list of some elements of the braid group.

        This is used both for illustration and testing purposes.

        EXAMPLES::

            sage: B=BraidGroup(3)
            sage: B.some_elements()
            [s0, s0*s1, (s0*s1)^3]
        """
        elements_list = [self.gen(0)]
        elements_list.append(self(range(1,self.strands())))
        elements_list.append(elements_list[-1]**self.strands())
        return elements_list

    def _permutation_braid_Tietze(self, p):
        """
        Helper for :meth:`_permutation_braid`.

        INPUT:

        - ``p`` -- a permutation.

        OUTPUT:

        The lexicographically smallest word that represents the braid,
        in Tietze list form.

        EXAMPLES::

            sage: B=BraidGroup(5)
            sage: P=Permutation([5, 3, 1, 2, 4])
            sage: B._permutation_braid_Tietze(P)
            (1, 2, 1, 3, 2, 4)
        """
        if p.length() == 0:
            return ()
        pl = p
        l = []
        while pl.length()>0:
            i = 1
            while i<max(pl):
                if pl(i)>pl(i+1):
                    l.append(i)
                    pl = Permutation([(i, i+1)])*pl
                    i = 1
                else:
                    i = i+1
        return tuple(l)

    @cached_method
    def _permutation_braid(self, p):
        """
        Return the braid that corresponds to the given permutation.

        It is the only braid with the following properties:

        - The braid induces the given permutation.

        - The braid is positive (that is, it can be writen without using the inverses of the generators).

        - Every two strands cross each other at most once.

        INPUT:

        - ``p`` -- a permutation.

        OUTPUT:

        The braid that corresponds to the permutation.

        EXAMPLES::

            sage: B = BraidGroup(5)
            sage: P = Permutation([5, 3, 1, 2, 4])
            sage: B._permutation_braid(P)
            s0*s1*s0*s2*s1*s3
        """
        return self(self._permutation_braid_Tietze(p))

    @cached_method
    def _LKB_matrix_(self, braid, variab):
        """
        Compute the Lawrence-Krammer-Bigelow representation matrix.

        The variables of the matrix must be given. This actual
        computation is done in this helper method for caching
        purposes.

        INPUT:

        - ``braid`` -- tuple of integers. The Tietze list of the
          braid.

        - ``variab`` -- string. the names of the variables that will
          appear in the matrix. They must be given as a string,
          separated by a comma

        OUTPUT:

        The LKB matrix of the braid, with respect to the variables.

        TESTS::

            sage: B=BraidGroup(3)
            sage: B._LKB_matrix_((2, 1, 2), 'x, y')
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: B._LKB_matrix_((1, 2, 1), 'x, y')
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: B._LKB_matrix_((-1, -2, -1, 2, 1, 2), 'x, y')
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        n = self.strands()
        if len(braid)>1:
            A = self._LKB_matrix_(braid[:1], variab)
            for i in braid[1:]:
                A = A*self._LKB_matrix_((i,), variab)
            return A
        l = list(Set(range(n)).subsets(2))
        R = LaurentPolynomialRing(IntegerRing(), variab)
        q = R.gens()[0]
        t = R.gens()[1]
        if len(braid)==0:
            return identity_matrix(R, len(l), sparse=True)
        A = matrix(R, len(l), sparse=True)
        if braid[0]>0:
            i = braid[0]-1
            for m in range(len(l)):
                j = min(l[m])
                k = max(l[m])
                if i==j-1:
                    A[l.index(Set([i, k])), m] = q
                    A[l.index(Set([i, j])), m] = q*q-q
                    A[l.index(Set([j, k])), m] = 1-q
                elif i==j and not j==k-1:
                    A[l.index(Set([j, k])), m] = 0
                    A[l.index(Set([j+1, k])), m] = 1
                elif k-1==i and not k-1==j:
                    A[l.index(Set([j, i])), m] = q
                    A[l.index(Set([j, k])), m] = 1-q
                    A[l.index(Set([i, k])), m] = (1-q)*q*t
                elif i==k:
                    A[l.index(Set([j, k])), m] = 0
                    A[l.index(Set([j, k+1])), m] = 1
                elif i==j and j==k-1:
                    A[l.index(Set([j, k])), m] = -t*q*q
                else:
                    A[l.index(Set([j, k])), m] = 1
            return A
        else:
            i = -braid[0]-1
            for m in range(len(l)):
                j = min(l[m])
                k = max(l[m])
                if i==j-1:
                    A[l.index(Set([j-1, k])), m] = 1
                elif i==j and not j==k-1:
                    A[l.index(Set([j+1, k])), m] = q**(-1)
                    A[l.index(Set([j, k])), m] = 1-q**(-1)
                    A[l.index(Set([j, j+1])), m] = t**(-1)*q**(-1)-t**(-1)*q**(-2)
                elif k-1==i and not k-1==j:
                    A[l.index(Set([j, k-1])), m] = 1
                elif i==k:
                    A[l.index(Set([j, k+1])), m] = q**(-1)
                    A[l.index(Set([j, k])), m] = 1-q**(-1)
                    A[l.index(Set([k, k+1])), m] = -q**(-1)+q**(-2)
                elif i==j and j==k-1:
                    A[l.index(Set([j, k])), m] = -t**(-1)*q**(-2)
                else:
                    A[l.index(Set([j, k])), m] = 1
            return A

    def dim_of_TL_space(self, drain_size):
        """
        Return the dimension of a particular Templerley--Lieb representation
        summand.
        
        Following the notation of ``TL_basis_with_drain``, the summand is the
        one corresponding to the number of drains being fixed to be
        ``drain_size``.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands (both
          included)

        OUTPUT:

        The dimension of the representation corresponding to the number of
        drains given by ``drain_size``.

        EXAMPLES:

        Calculation of the dimension of the representation of `B_8`
        corresponding to having `2` drains::

            sage: B = BraidGroup(8)
            sage: B.dim_of_TL_space(2)
            28

        The direct sum of endomorphism spaces of these vector spaces make up
        the entire Temperley--Lieb algebra::

            sage: import sage.combinat.diagram_algebras as da
            sage: B = BraidGroup(6)
            sage: dimensions = [B.dim_of_TL_space(d)**2 for d in [0, 2, 4, 6]]
            sage: total_dim = sum(dimensions)
            sage: total_dim == len(da.temperley_lieb_diagrams(6))
            True
        """
        n = self.strands()
        d = drain_size
        if d > n:
            raise ValueError("Number of drains must not exceed number of strands")
        if (n+d) % 2 == 1:
            raise ValueError("Parity of strands and drains must agree")
        if n % 2 == 0:
            k = n / 2
            g = k - 1
            m = (2*k - d) / 2
            return Integer(2*g + 1).binomial(m) - Integer(2*g + 1).binomial(m - 2)
        else:
            k = (n - 1) / 2
            g = k
            m = (2*k + 1 - d) / 2
            return Integer(2*g).binomial(m) - Integer(2*g).binomial(m - 2)

    def TL_basis_with_drain(self, drain_size):
        """
        Return a basis of a summand of the Temperley--Lieb--Jones
        representation.
        
        The basis elements are given by non-intersecting pairings of `n+d`
        points in a square with `n` points marked 'on the top' and `d` points
        'on the bottom' so that every bottom point is paired with a top point.
        Here, `n` is the number of strands of the braid group, and `d` is
        specified by ``drain_size``.

        A basis element is specified as a list of integers obtained by
        considering the pairings as obtained as the 'highest term' of trivalent
        trees marked by Jones--Wenzl projectors (see e.g. [Wan]_). In practice,
        this is a list of non-negative integers whose first element is
        ``drain_size``, whose last element is `0`, and satisfying that
        consecutive integers have difference `1`. Moreover, the length of each
        basis element is `n+1`.

        Given these rules, the list of lists is constructed recursively in the
        natural way.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands (both
          included)

        OUTPUT:

        A list of basis elements, each of which is a list of integers.

        EXAMPLES:

        We calculate the basis for the appropriate vector space for `B_5` when
        `d = 3`::

            sage: B = BraidGroup(5)
            sage: B.TL_basis_with_drain(3)
            [[3, 4, 3, 2, 1, 0],
             [3, 2, 3, 2, 1, 0],
             [3, 2, 1, 2, 1, 0],
             [3, 2, 1, 0, 1, 0]]

        The number of basis elements hopefully correponds to the general
        formula for the dimension of the representation spaces::

            sage: B = BraidGroup(10)
            sage: d = 2
            sage: B.dim_of_TL_space(d) == len(B.TL_basis_with_drain(d))
            True

        REFERENCES:

        .. [Wan] Zhenghan Wang. Tolological quantum computation. Providence,
                 RI: American Mathematical Society (AMS), 2010.
                 ISBN 978-0-8218-4930-9
        """
        def fill_out_forest(forest, treesize):
            if len(forest) == 0:
                raise ValueError("Forest has to start with a tree")
            if forest[0][0]+treesize % 2 == 0:
                raise ValueError("Parity mismatch in forest creation")
            # Loop over all trees
            newforest = list(forest)
            for tree in forest:
                if len(tree) < treesize:
                    newtreeup = list(tree)
                    newtreedown = list(tree)
                    newforest.remove(tree)  # Cut down the original tree
                    # Add two greater trees, admissibly
                    if tree[-1] < treesize - len(tree) + 1:
                        newtreeup.append(tree[-1] + 1)
                        newforest.append(newtreeup)
                    if tree[-1] > 0:
                        newtreedown.append(tree[-1] - 1)
                        newforest.append(newtreedown)
            # Are we there yet?
            if len(newforest[0]) == treesize:
                return newforest
            else:
                return fill_out_forest(newforest, treesize)

        n = self.strands()
        d = drain_size
        if d > n:
            raise ValueError("Number of drains must not exceed number of strands")
        if (n+d) % 2 == 1:
            raise ValueError("Parity of strands and drains must agree")
        basis = [[d]]  # Let's start out with no elements and recursively fill out
        forest = fill_out_forest(basis, n-1)
        for tree in forest:
            tree.extend([1, 0])
        return forest

    @cached_method
    def create_TL_rep(self, drain_size, variab='A', ring=IntegerRing()):
        """
        Return representation matrices of the Temperley--Lieb--Jones
        representation of standard braid group generators.
        
        The basis is given by non-intersecting pairings of `(n+d)` points,
        where `n` is the number of strands, and `d` is given by ``drain_size``,
        and the pairings satisfy certain rules. See
        ``BraidGroup_class.TL_basis_with_drain()`` for details. This basis has
        the useful property that all resulting entries can be regarded as
        Laurent polynomials.

        We use the convention that the eigenvalues of the standard generators
        are `1` and `-A^4`, where `A` is the generator of the Laurent
        polynomial ring.

        When `d = n - 2` and the variables are picked appropriately, the
        resulting representation is equivalent to the reduced Burau
        representation. When `d = n`, the resulting representation is trivial
        and 1-dimensional.

        Store the result of the calculation as part of the braid group.

        INPUT:

        - ``drain_size`` -- integer between 0 and the number of strands (both
          included)
        - ``variab`` -- string (default: ``'A'``); the name of the
          variable in the entries of the matrices
        - ``ring`` -- ring (default: ``IntegerRing()``); the ring to which the
          coefficients of the polynomial entries belong

        OUTPUT:

        A list of matrices corresponding to the representations of each of the
        standard generators.

        EXAMPLES::
            sage: B = BraidGroup(4)
            sage: B.create_TL_rep(0)
            [
            [   1    0]  [-A^4  A^2]  [   1    0]
            [ A^2 -A^4], [   0    1], [ A^2 -A^4]
            ]
            sage: B.create_TL_rep(0, ring=GF(2))
            [
            [  1   0]  [A^4 A^2]  [  1   0]
            [A^2 A^4], [  0   1], [A^2 A^4]
            ]
            sage: B = BraidGroup(8)
            sage: B.create_TL_rep(8)
            [[1], [1], [1], [1], [1], [1], [1]]
        """
        n = self.strands()
        d = drain_size
        basis = self.TL_basis_with_drain(d)
        auxmat = matrix(n-1, len(basis))
        for i in range(1, n):
            for v in range(len(basis)):
                tree = basis[v]
                if tree[i-1] < tree[i] and tree[i+1] < tree[i]:
                    auxmat[i-1, v] = v
                if tree[i-1] > tree[i] and tree[i+1] > tree[i]:
                    newtree = list(tree)
                    newtree[i] += 2
                    auxmat[i-1, v] = basis.index(newtree)
                if tree[i-1] > tree[i] and tree[i+1] < tree[i]:
                    newtree = list(tree)
                    newtree[i-1] -= 2
                    j = 2
                    while newtree[i-j] != newtree[i] and i-j >= 0:
                        newtree[i-j] -= 2
                        j += 1
                    if newtree in basis:
                        auxmat[i-1, v] = basis.index(newtree)
                    else:
                        auxmat[i-1, v] = -1
                if tree[i-1] < tree[i] and tree[i+1] > tree[i]:
                    newtree = list(tree)
                    newtree[i+1] -= 2
                    j = 2
                    while newtree[i+j] != newtree[i] and i+j <= n:
                        newtree[i+j] -= 2
                        j += 1
                    if newtree in basis:
                        auxmat[i-1, v] = basis.index(newtree)
                    else:
                        auxmat[i-1, v] = -1
        repmat = []
        R = LaurentPolynomialRing(ring, variab)
        A = R.gens()[0]
        for i in range(1, n):
            repmatnew = identity_matrix(R, len(basis), sparse=True)
            for v in range(len(basis)):
                newmatentry = auxmat[i-1, v]
                if newmatentry == v:
                    repmatnew[v, v] = -A**4
                elif newmatentry >= 0:
                    repmatnew[newmatentry, v] = A**2
            repmat.append(repmatnew)
        return repmat

    def mapping_class_action(self, F):
        """
        Return the action of self in the free group F as mapping class group.

        This action corresponds to the action of the braid over the
        punctured disk, whose fundamental group is the free group on
        as many generators as strands.

        In Sage, this action is the result of multiplying a free group
        element with a braid. So you generally do not have to
        construct this action yourself.

        OUTPUT:

        A :class:`MappingClassGroupAction`.

        EXAMPLES ::

            sage: B = BraidGroup(3)
            sage: B.inject_variables()
            Defining s0, s1
            sage: F.<a,b,c> = FreeGroup(3)
            sage: A = B.mapping_class_action(F)
            sage: A(a,s0)
            a*b*a^-1
            sage: a * s0    # simpler notation
            a*b*a^-1
        """
        return MappingClassGroupAction(self, F)

    def _get_action_(self, S, op, self_on_left):
        """
        Let the coercion system discover actions of the braid group on free groups.

            sage: B.<b0,b1,b2> = BraidGroup()
            sage: F.<f0,f1,f2,f3> = FreeGroup()
            sage: f1 * b1
            f1*f2*f1^-1

            sage: from sage.structure.all import get_coercion_model
            sage: cm = get_coercion_model()
            sage: cm.explain(f1, b1, operator.mul)
            Action discovered.
                Right action by Braid group on 4 strands on Free Group on generators {f0, f1, f2, f3}
            Result lives in Free Group on generators {f0, f1, f2, f3}
            Free Group on generators {f0, f1, f2, f3}
            sage: cm.explain(b1, f1, operator.mul)
            Will try _r_action and _l_action
            Unknown result parent.
        """
        import operator
        if is_FreeGroup(S) and op==operator.mul and not self_on_left:
            return self.mapping_class_action(S)
        return None



def BraidGroup(n=None, names='s'):
    """
    Construct a Braid Group

    INPUT:

    - ``n`` -- integer or ``None`` (default). The number of
      strands. If not specified the ``names`` are counted and the
      group is assumed to have one more strand than generators.

    - ``names`` -- string or list/tuple/iterable of strings (default:
      ``'x'``). The generator names or name prefix.

    EXAMPLES::

        sage: B.<a,b> = BraidGroup();  B
        Braid group on 3 strands
        sage: H = BraidGroup('a, b')
        sage: B is H
        True
        sage: BraidGroup(3)
        Braid group on 3 strands

    The entry can be either a string with the names of the generators,
    or the number of generators and the prefix of the names to be
    given. The default prefix is ``'s'`` ::

        sage: B=BraidGroup(3); B.generators()
        (s0, s1)
        sage: BraidGroup(3, 'g').generators()
        (g0, g1)

    Since the word problem for the braid groups is solvable, their Cayley graph
    can be localy obtained as follows (see :trac:`16059`)::

        sage: def ball(group, radius):
        ....:     ret = set()
        ....:     ret.add(group.one())
        ....:     for length in range(1, radius):
        ....:         for w in Words(alphabet=group.gens(), length=length):
        ....:              ret.add(prod(w))
        ....:     return ret
        sage: B = BraidGroup(4)
        sage: GB = B.cayley_graph(elements=ball(B, 4), generators=B.gens()); GB
        Digraph on 31 vertices

    Since the braid group has nontrivial relations, this graph contains less
    vertices than the one associated to the free group (which is a tree)::

        sage: F = FreeGroup(3)
        sage: GF = F.cayley_graph(elements=ball(F, 4), generators=F.gens()); GF
        Digraph on 40 vertices

    TESTS::

        sage: G1 = BraidGroup(3, 'a,b')
        sage: G2 = BraidGroup('a,b')
        sage: G3.<a,b> = BraidGroup()
        sage: G1 is G2, G2 is G3
        (True, True)
    """
    # Support Freegroup('a,b') syntax
    if n is not None:
        try:
            n = Integer(n)-1
        except TypeError:
            names = n
            n = None
    # derive n from counting names
    if n is None:
        if isinstance(names, six.string_types):
            n = len(names.split(','))
        else:
            names = list(names)
            n = len(names)
    from sage.structure.parent import normalize_names
    names = tuple(normalize_names(n, names))
    return BraidGroup_class(names)



class MappingClassGroupAction(Action):
    r"""
    The action of the braid group the free group as the mapping class
    group of the punctured disk.

    That is, this action is the action of the braid over the punctured
    disk, whose fundamental group is the free group on as many
    generators as strands.

    This action is defined as follows:

    .. MATH::

        x_j \cdot \sigma_i=\begin{cases}
        x_{j}\cdot x_{j+1}\cdot {x_j}^{-1} & \text{if $i=j$} \\
        x_{j-1} & \text{if $i=j-1$} \\
        x_{j} & \text{otherwise}
        \end{cases},

    where $\sigma_i$ are the generators of the braid group on $n$
    strands, and $x_j$ the generators of the free group of rank $n$.

    You should left multiplication of the free group element by the
    braid to compute the action. Alternatively, use the
    :meth:`~sage.groups.braid.BraidGroup_class.mapping_class_action`
    method of the braid group to constuct this action.

    EXAMPLES::

        sage: B.<s0,s1,s2> = BraidGroup(4)
        sage: F.<x0,x1,x2,x3> = FreeGroup(4)
        sage: x0 * s1
        x0
        sage: x1 * s1
        x1*x2*x1^-1
        sage: x1^-1 * s1
        x1*x2^-1*x1^-1

        sage: A = B.mapping_class_action(F)
        sage: A
        Right action by Braid group on 4 strands on Free Group on generators {x0, x1, x2, x3}
        sage: A(x0, s1)
        x0
        sage: A(x1, s1)
        x1*x2*x1^-1
        sage: A(x1^-1, s1)
        x1*x2^-1*x1^-1
    """
    def __init__(self, G, M, is_left=0):
        """
        TESTS::

            sage: B = BraidGroup(3)
            sage: G = FreeGroup('a, b, c')
            sage: B.mapping_class_action(G) # indirect doctest
            Right action by Braid group on 3 strands on Free Group on generators {a, b, c}
        """
        import operator
        Action.__init__(self, G, M, is_left, operator.mul)

    def _call_(self, x, b):
        """
        Return the action of ``b`` on ``x``.

        INPUT:

        - ``x`` -- a free group element.

        - ``b`` -- a braid.

        OUTPUT:

        A new braid.

        TESTS::

            sage: B = BraidGroup(3)
            sage: G = FreeGroup('a, b, c')
            sage: A = B.mapping_class_action(G)
            sage: A(G.0, B.0) # indirect doctest
            a*b*a^-1
            sage: A(G.1, B.0) # indirect doctest
            a
        """
        t = x.Tietze()
        for j in b.Tietze():
            s=[]
            for i in t:
                if j==i and i>0:
                    s += [i, i+1, -i]
                elif j==-i and i<0:
                    s += [-i, i-1, i]
                elif j==-i and i>0:
                    s += [i+1]
                elif j==i and i<0:
                    s += [i-1]
                elif i>0 and j==i-1:
                    s += [i-1]
                elif i<0 and j==-i-1:
                    s += [i+1]
                elif i>0 and -j==i-1:
                    s += [-i, i-1, i]
                elif i<0 and j==i+1:
                    s += [i, i+1, -i]
                else:
                    s += [i]
            t = s
        return self.codomain()(t)

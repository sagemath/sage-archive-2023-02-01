r"""
Hyperplanes

In Sage, these are used as a notation for linear
constraints (or hyperplanes), that is a linear expression `x+2y+3z-5`
stands for the hyperplane with the equation `x+2y+3z=5`

For example, to create a single hyperplane with equation `3x + 2y - 5z
= 7` ::

    sage: h = Hyperplane([3,2,-5,7])
    sage: h.equation()
    Rational Field
    (3, 2, -5)*X = 7
    [3, 2, -5, 7]
    sage: h.equation(true)  # suppress_printing = true
    [3, 2, -5, 7]
    sage: h.normal()
    (3, 2, -5)
    sage: h.point()
    (7/3, 0, 0)
    sage: h.linear_part()
    Vector space of degree 3 and dimension 2 over Rational Field
    Basis matrix:
    [  1   0 3/5]
    [  0   1 2/5]
    sage: h.show()
    sage: h = Hyperplane([1,1/2,1/2,3])
    sage: h  # denominators are cleared
    Hyperplane over Rational Field
    (2, 1, 1)*X = 6
    sage: h = Hyperplane([2,3,4,0])
    sage: h.base_field()
    Rational Field
    sage: h.change_base_field(GF(3))
    Hyperplane over Finite Field of size 3
    (2, 0, 1)*X = 0
"""
#*****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from colorsys import hsv_to_rgb
from copy import copy, deepcopy
from sage.calculus.functional import expand
from sage.calculus.var import var
from sage.combinat.combinat import stirling_number2
from sage.combinat.posets.posets import Poset
from sage.functions.generalized import sign
from sage.functions.other import sqrt
from sage.geometry.polyhedron.all import Polyhedron
from sage.graphs.all import graphs
from sage.matrix.constructor import matrix, random_matrix, zero_matrix
from sage.misc.flatten import flatten
from sage.misc.prandom import random
from sage.misc.misc import powerset
from sage.misc.misc_c import prod
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.plot.line import line
from sage.plot.colors import Color
from sage.plot.graphics import Graphics
from sage.plot.plot import plot, parametric_plot
from sage.plot.point import point
from sage.plot.text import text
from sage.plot.plot3d.parametric_plot3d import parametric_plot3d
from sage.plot.plot3d.shapes2 import text3d
from sage.rings.arith import lcm, binomial
from sage.rings.finite_rings.constructor import GF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SR, var


class AffineSubspace(SageObject):
    r"""
    Class for an affine space (a translation of a linear subspace).

    INPUT:

    - ``W`` -- vector space
    - ``p`` -- list representing a point in W

    OUTPUT:

    - AffineSubspace

    EXAMPLES::

        sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
        sage: a.dim()
        4
        sage: a.point()
        (1, 0, 0, 0)
        sage: a.linear_part()
        Vector space of dimension 4 over Rational Field
        sage: a
        Affine space p + W where:
           p = (1, 0, 0, 0)
           W = Vector space of dimension 4 over Rational Field.
        sage: b = AffineSubspace((1,0,0,0),matrix(QQ, [[1,2,3,4]]).right_kernel())
        sage: c = AffineSubspace((0,2,0,0),matrix(QQ, [[0,0,1,2]]).right_kernel())
        sage: b.intersection(c)
        Affine space p + W where:
           p = (-3, 2, 0, 0)
           W = Vector space of degree 4 and dimension 2 over Rational Field
        Basis matrix:
        [  1   0  -1 1/2]
        [  0   1  -2   1].
        sage: b < a
        True
        sage: c < b
        False
        sage: A = AffineSubspace([8,38,21,250],VectorSpace(GF(19),4))
        sage: A
        Affine space p + W where:
            p = (8, 0, 2, 3)
            W = Vector space of dimension 4 over Finite Field of size 19.
        sage: A = AffineSubspace([2],VectorSpace(QQ,0))
        sage: A.point()
        (2)
        sage: A.linear_part()
        Vector space of dimension 0 over Rational Field
        sage: A.linear_part().basis_matrix()
        []
    """

    def __init__(self, p, V):
        r"""
        Initializes an AffineSubspace.  See `AffineSubspace` for more examples.

        INPUT:

        - ``W`` -- vector space
        - ``p`` -- list representing a point in W

        OUTPUT:

        - AffineSubspace

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
            sage: a
            Affine space p + W where:
               p = (1, 0, 0, 0)
               W = Vector space of dimension 4 over Rational Field.
        """
        if V.base_ring()==ZZ:
            V = V.change_ring(QQ)
        self._linear_part = V
        self._point = vector(V.base_field(), p)

    def __repr__(self):
        r"""
        String representation for an AffineSubspace.

        INPUT:

        - None

        OUTPUT:

        - string

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
            sage: a.__repr__()
            'Affine space p + W where:\n   p = (1, 0, 0, 0)\n   W = Vector space of dimension 4 over Rational Field.'
        """
        return "Affine space p + W where:\n   p = "+str(self._point)+"\n   W = "+str(self._linear_part)+"."

    def __eq__(self, other):
        r"""
        Tests whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- AffineSubspace

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0],matrix([[1,0,0]]).right_kernel())
            sage: b = AffineSubspace([2,0,0],matrix([[1,0,0]]).right_kernel())
            sage: c = AffineSubspace([1,1,0],matrix([[1,0,0]]).right_kernel())
            sage: a == b
            False
            sage: a == c
            True
        """
        try:
            V = self._linear_part
            W = other._linear_part
            if V == W and self._point - other._point in V:
                return True
        except:
            return False
        return False

    def __ne__(self, other):
        r"""
        Tests whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- AffineSubspace

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0],matrix([[1,0,0]]).right_kernel())
            sage: b = AffineSubspace([2,0,0],matrix([[1,0,0]]).right_kernel())
            sage: a == b
            False
            sage: a != b
            True
            sage: a != a
            False
        """
        return not self == other

    def __le__(self, other):
        r"""
        Tests whether ``self`` is an affine subspace of ``other``.

        INPUT:

        - ``other`` -- AffineSubspace

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: W1 = V.subspace([[1,0,0],[0,1,0]])
            sage: W2 = V.subspace([[1,0,0]])
            sage: a = AffineSubspace([1,2,3],W1)
            sage: b = AffineSubspace([1,2,3],W2)
            sage: a <= b
            False
            sage: a <= a
            True
            sage: b <= a
            True
        """
        V = self._linear_part
        W = other._linear_part
        if V.is_subspace(W) and self._point-other._point in W:
            return True
        else:
            return False

    def __lt__(self, other):
        r"""
        Tests whether ``self`` is a proper affine subspace of ``other``.

        INPUT:

        - ``other`` -- AffineSubspace

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: W1 = V.subspace([[1,0,0],[0,1,0]])
            sage: W2 = V.subspace([[1,0,0]])
            sage: a = AffineSubspace([1,2,3],W1)
            sage: b = AffineSubspace([1,2,3],W2)
            sage: a < b
            False
            sage: a < a
            False
            sage: b < a
            True
        """
        if self <= other and not self == other:
            return True
        else:
            return False

    def __ge__(self, other):
        r"""
        Tests whether ``other`` is an affine subspace of ``self``.

        INPUT:

        - ``other`` -- AffineSubspace

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: W1 = V.subspace([[1,0,0],[0,1,0]])
            sage: W2 = V.subspace([[1,0,0]])
            sage: a = AffineSubspace([1,2,3],W1)
            sage: b = AffineSubspace([1,2,3],W2)
            sage: a >= b
            True
            sage: a >= a
            True
            sage: b >= a
            False
        """
        V = self._linear_part
        W = other._linear_part
        if W.is_subspace(V) and self._point-other._point in V:
            return True
        else:
            return False

    def __gt__(self, other):
        r"""
        Tests whether ``other`` is a proper affine subspace of ``self``.

        INPUT:

        - ``other`` -- AffineSubspace

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: W1 = V.subspace([[1,0,0],[0,1,0]])
            sage: W2 = V.subspace([[1,0,0]])
            sage: a = AffineSubspace([1,2,3],W1)
            sage: b = AffineSubspace([1,2,3],W2)
            sage: a > b
            True
            sage: a > a
            False
            sage: b > a
            False
        """
        if self >= other and not self == other:
            return True
        else:
            return False

    def __contains__(self, q):
        r"""
        Tests whether the point ``q`` is in the affine space.

        INPUT:

        - ``q`` -- point (as a vector, list, or tuple)

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0],matrix([[1,0,0]]).right_kernel())
            sage: (1,1,0) in a
            True
            sage: (0,0,0) in a
            False
        """
        if type(q) in [list, tuple]:
            q = vector(self.base_field(),q)
        return self._point - q in self._linear_part

    def base_field(self):
        r"""
        Returns the base field of the affine space.

        INPUT:

        - None

        OUTPUT:

        - field

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
            sage: a.base_field()
            Rational Field
            sage: b = AffineSubspace([1,0,0,0],VectorSpace(GF(5),4))
            sage: b.base_field()
            Finite Field of size 5
        """
        return self.linear_part().base_field()

    def linear_part(self):
        r"""
        The linear part of the affine space.

        INPUT:

        - None

        OUTPUT:

        - vector space

        EXAMPLES::

            sage: A = AffineSubspace([2,3,1], matrix(QQ,[[1,2,3]]).right_kernel())
            sage: A.linear_part()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 -1/3]
            [   0    1 -2/3]
        """
        return copy(self._linear_part)

    def point(self):
        r"""
        A point ``p`` in the affine space.

        INPUT:

        - None

        OUTPUT:

        - vector

        EXAMPLES::

            sage: A = AffineSubspace([2,3,1], VectorSpace(QQ,3))
            sage: A.point()
            (2, 3, 1)
        """
        return copy(self._point)

    def dim(self):
        r"""
        The dimension of the affine space.

        INPUT:

        - None

        OUTPUT:

        - Integer

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
            sage: a.dim()
            4
        """
        return self.linear_part().dimension()

    def intersection(self, other):
        r"""
        The intersection of ``self`` with ``other``.

        INPUT:

        - ``other`` -- AffineSubspace

        OUTPUT:

        - AffineSubspace  (or -1 if the intersection is empty)

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: U = V.subspace([(1,0,0),(0,1,0)])
            sage: W = V.subspace([(0,1,0),(0,0,1)])
            sage: A = AffineSubspace((0,0,0),U)
            sage: B = AffineSubspace((1,1,1),W)
            sage: A.intersection(B)
            Affine space p + W where:
               p = (1, 1, 0)
               W = Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0].
            sage: C = AffineSubspace((0,0,1),U)
            sage: A.intersection(C)
            -1
            sage: D = AffineSubspace([1,2,3],VectorSpace(GF(5),3))
            sage: E = AffineSubspace([3,4,5],VectorSpace(GF(5),3))
            sage: D.intersection(E)
            Affine space p + W where:
               p = (3, 4, 0)
               W = Vector space of dimension 3 over Finite Field of size 5.
        """
        if self.linear_part().ambient_vector_space()!=other.linear_part().ambient_vector_space():
            raise UserWarning('incompatible ambient vector spaces')
        elif self.dim()==0:
            if self<=other:
                return copy(self)
            else:
                return -1 # empty intersection
        elif other.dim()==0:
            if other<=self:
                return copy(other)
            else:
                return -1 # empty intersection
        else:
            m = self.linear_part().matrix()
            n = other.linear_part().matrix()
            p = self.point()
            q = other.point()
            M = m.stack(n)
            v = q-p
            if M.rank() != (M.stack(v)).rank():
                return -1   # the intersection is empty
            else:
                t = M.solve_left(v)
                new_p = p + t[:m.nrows()]*m
                new_V = self.linear_part().intersection(other._linear_part)
                return AffineSubspace(new_p, new_V)

    def _isomorphism_with_Kn(self, v):
        r"""


        INPUT:

        - ``v`` -- vector in the ambient space

        OUTPUT:

        - vector

        EXAMPLES::

            sage: a = AffineSubspace([1,0,0,0],VectorSpace(QQ,4))
            sage: v = AffineSubspace([2,1,3,2],VectorSpace(QQ,4))._point
            sage: a._isomorphism_with_Kn(v)
            (1, 1, 3, 2)
        """
        v = vector(self.base_field(),v) # in case v was just a list
        # If the base field is QQ, approximate an orthogonal projection 
        # for the sake of visualizations.
        if self.base_field()==QQ:
            W = self.linear_part().basis_matrix()
            #begin to orthonormalize W
            G, M = W.gram_schmidt()
            g = G*G.transpose()
            q = g.apply_map(sqrt).inverse()
            # to keep things defined over QQ:
            q = q.apply_map(lambda x: QQ(round(RR(x),2)))
            u = q*G  # u is the approximately orthonormalized W
            return u.solve_left(v-self.point())
        else: # finite field
            p = self.point()
            W = self.linear_part()
            V = VectorSpace(W.base(), W.dimension())
            f = W.hom(V.gens())
            v = vector(v)
            return f(v-p)


class Hyperplane(AffineSubspace):
    def __init__(self, H, K=QQ):
        r"""
        The argument ``H`` is a hyperplane, given as a list `[a_1, ..., a_n, a]`
        representing the equation `a_1x_1 + ... + a_nx_n = a`. An optional
        field, ``K``, may also be provided.

        INPUT:

        - ``H`` -- list of integers representing a hyperplane
        - ``K`` -- field (default: ``QQ``)

        OUTPUT:

        - Hyperplane

        EXAMPLES::

            sage: h = Hyperplane([1,3,2,2,8])
            sage: h
            Hyperplane over Rational Field
            (1, 3, 2, 2)*X = 8
            sage: h.normal()
            (1, 3, 2, 2)
            sage: h.dim()
            3
            sage: h.point()
            (8, 0, 0, 0)
            sage: h.linear_part()
            Vector space of degree 4 and dimension 3 over Rational Field
            Basis matrix:
            [   1    0    0 -1/2]
            [   0    1    0 -3/2]
            [   0    0    1   -1]
            sage: H = Hyperplane([1,1/2,1/2,3])
            sage: H.equation()  # denominators are cleared
            Rational Field
            (2, 1, 1)*X = 6
            [2, 1, 1, 6]
            sage: Hyperplane([1,3,2,2,8],GF(5))
            Hyperplane over Finite Field of size 5
            (1, 3, 2, 2)*X = 3
        """
        m = matrix(K, H[:-1])
        if m.is_zero():
            raise UserWarning('not a hyperplane')
        else:
            p = m.solve_right(vector([H[-1]]))
            AffineSubspace.__init__(self, p, m.right_kernel())
            self._raw_data = H
            # if working over a finite field:
            self._equation = map(K,H)
            # Force the equations to have integer coefficients
            if K == QQ:
                v = vector(self._equation)
                self._equation = list(lcm([i.denom() for i in v])*v)
            self._normal = vector(K,self._equation[:-1])
            self._base_field = K

    def __repr__(self):
        r"""
        String representation of a hyperplane.

        INPUT:

        - None

        OUTPUT:

        - string

        EXAMPLES::

            sage: H = Hyperplane([1,0,1])
            sage: H.__repr__()
            'Hyperplane over Rational Field\n(1, 0)*X = 1'
        """
        return "Hyperplane over "+str(self.base_field())+"\n"+str(self._normal)+"*X = "+str(self._equation[-1])

    def base_field(self):
        r"""
        The base field of the hyperplane.

        INPUT:

        - None

        OUTPUT:

        - field

        EXAMPLES::

            sage: H = Hyperplane([1,2,3,-7])
            sage: H.base_field()
            Rational Field
            sage: h = Hyperplane([3,2,6,8],GF(7))
            sage: h.base_field()
            Finite Field of size 7

        .. SEEALSO::

            :meth:`change_base_field`
        """
        return self._base_field

    def equation(self, suppress_printing=False):
        r"""
        Returns the list `[a_1,\dots,a_n,a]` representing the hyperplane with
        equation `a_1x_1 + \dots + a_nx_n = a`.

        INPUT:

        - suppress_printing -- boolean (default: False)

        OUTPUT:

        - list

        EXAMPLES::

            sage: h = Hyperplane([1,2,3,4])
            sage: r = h.equation()
            Rational Field
            (1, 2, 3)*X = 4
            sage: r
            [1, 2, 3, 4]
            sage: h.equation(True)  # suppress printing
            [1, 2, 3, 4]
            sage: h = Hyperplane([10,29,115,24],GF(11))
            sage: h.equation()
            Finite Field of size 11
            (10, 7, 5)*X = 2
            [10, 7, 5, 2]
        """
        if not suppress_printing:
            F = self.base_field()
            print str(F)+"\n"+str(self._normal)+"*X = "+str(self._equation[-1])
        return deepcopy(self._equation)

    def change_base_field(self, F):
        r"""
        Returns a hyperplane defined over ``F``.

        INPUT:

        - ``F`` -- the rational field or a finite field

        OUTPUT:

        - hyperplane defined over ``F``

        EXAMPLES::

            sage: a = Hyperplane([8,3,7])
            sage: r = a.change_base_field(GF(7))
            sage: r
            Hyperplane over Finite Field of size 7
            (1, 3)*X = 0
        """
        eq = map(F,self.equation(True))
        return Hyperplane(eq,F)

    def normal(self):
        r"""
        A vector perpendicular to the hyperplane.

        INPUT:

        - None

        OUTPUT:

        - vector

        EXAMPLES::

            sage: h = Hyperplane([1,2,3,4])
            sage: h.normal()
            (1, 2, 3)
            sage: h = Hyperplane([12,62,45,69],GF(7))
            sage: h.normal()
            (5, 6, 3)
        """
        return self._normal

    def _pretty_print_equation(self, latex=True):
        r"""
        Returns a nice string for the equation of the hyperplane.  If ``latex``
        is ``True``, returns a latex-formatted string.  This function is only
        used in the ``plot`` and ``show`` functions, so it only applies to
        hyperplanes of dimension 3 or less.

        INPUT:

        - ``latex`` -- boolean (default: True)

        OUTPUT:

        - string

        EXAMPLES::

            sage: Hyperplane([-3,2])._pretty_print_equation()
            '$-3x = 2$'
            sage: Hyperplane([-3,2])._pretty_print_equation(false)
            '-3x = 2'
            sage: Hyperplane([1,3,-5])._pretty_print_equation()
            '$x + 3y = -5$'
            sage: Hyperplane([1,0,-1,4])._pretty_print_equation()
            '$x - z = 4$'
        """
        e = self.equation(True)
        if self.dim()==0:
            x = SR('x')  # this helps with +/- characters
            s = str(e[0]*x).replace('*','') + ' = ' + str(e[1])
        elif self.dim()==1:
            x, y = SR('x'), SR('y')
            s = str(e[0]*x+e[1]*y).replace('*','') + ' = ' + str(e[2])
        elif self.dim()==2:
            x, y, z = SR('x'), SR('y'), SR('z')
            s = str(e[0]*x+e[1]*y+e[2]*z).replace('*','') + ' = ' + str(e[3])
        else:
            s = '' # return blank string for dimensions > 2
        if latex:
            s = '$' + s + '$'
        return s

    def show(self, **kwds):
        r"""
        Displays the hyperplane.

        INPUT:

        - **kwds -- show options: see below

        OUTPUT:

        - None

        PLOT OPTIONS::

            Beside the usual options for show (enter show?), the show command for
            hyperplanes includes the following:

            - hyperplane_label -- Boolean value or string (default: ``True``).
              If ``True``, the hyperplane is labeled with its equation, if a
              string, it is labeled by that string, if ``False``, it is not
              labeled.

            - label_color -- Color for hyperplane_label (default: black).

            - label_fontsize -- Size for hyperplane_label font (default: 14).
              (Does not work in 3d, yet.)

            - label_offset -- Amount by which label is offset from self.point()
              (default: 0-dim: 0.1, 1-dim: (0,1), 2-dim: (0,0,0.2))

            - point_size -- Size of points in a zero-dimensional arrangement or
              of an arrangement over a finite field (default: 50).

            - ranges -- Range for the parameters for the parametric plot of the
              hyperplane. If a single positive number ``r`` is given for the
              value of ``ranges``, then the ranges for all parameters are set to
              [-r,r].  Otherwise, for a line in the plane, ``ranges`` has the
              form [a,b] (default: [-3,3]), and for a plane in 3-space, the
              ``ranges`` has the form [[a,b],[c,d]] (default: [[-3,3],[-3,3]]).
              (The ranges are centered around self.point().)

        EXAMPLES::

            sage: a = Hyperplane([3,4])
            sage: a.show()
            sage: a.show(point_size=100)
            sage: b = Hyperplane([3,4,5])
            sage: b.show()
            sage: b.show(ranges=(1,5),label_offset=(2,-1))
            sage: b.show(axes=false,hyperplane_label='blue line',label_offset=(0,1))
            sage: c = Hyperplane([2,3,4,5])
            sage: c.show()
            sage: c.show(label_offset=(1,0,1), color='green', label_color='red', frame=False)

        NOTES::

            For more examples, see :meth:`plot` (``self.plot?`` or ``Hyperplane.plot?``).
        """
        p = self.plot(**kwds)
        for k in ['hyperplane_label','label_color', 'label_fontsize','label_offset',
                'point_size', 'ranges']:
            if kwds.has_key(k):
                del kwds[k]
        if self.dim() == 0 and not kwds.has_key('ymin'):
            kwds['ymin'] = -0.3
        p.show(**kwds)

    def plot(self, **kwds):
        r"""
        Returns the plot of a given hyperplane.

        INPUT:

        - **kwds -- plot options: see below

        OUTPUT:

        - Graphics

        PLOT OPTIONS::

            Beside the usual plot options (enter plot?), the plot command for
            hyperplanes includes the following:

            - hyperplane_label -- Boolean value or string (default: ``True``).
              If ``True``, the hyperplane is labeled with its equation, if a
              string, it is labeled by that string, if ``False``, it is not
              labeled.

            - label_color -- Color for hyperplane_label (default: black).

            - label_fontsize -- Size for hyperplane_label font (default: 14).
              (Does not work in 3d, yet.)

            - label_offset -- Amount by which label is offset from self.point()
              (default: 0-dim: 0.1, 1-dim: (0,1), 2-dim: (0,0,0.2))

            - point_size -- Size of points in a zero-dimensional arrangement or
              of an arrangement over a finite field (default: 50).

            - ranges -- Range for the parameters for the parametric plot of the
              hyperplane. If a single positive number ``r`` is given for the
              value of ``ranges``, then the ranges for all parameters are set to
              [-r,r].  Otherwise, for a line in the plane, ``ranges`` has the
              form [a,b] (default: [-3,3]), and for a plane in 3-space, the
              ``ranges`` has the form [[a,b],[c,d]] (default: [[-3,3],[-3,3]]).
              (The ranges are centered around self.point().)

        EXAMPLES::

            sage: a = Hyperplane([3,4])
            sage: a.plot()
            sage: a.plot(point_size=100,hyperplane_label='hello')
            sage: b = Hyperplane([3,4,5])
            sage: b.plot()
            sage: b.plot(ranges=(1,5),label_offset=(2,-1))
            sage: c = Hyperplane([2,3,4,5])
            sage: c.plot()
            sage: c.plot(label_offset=(1,0,1), color='green', label_color='red', frame=False)
            sage: d = Hyperplane([-3,2,2,3])
            sage: d.plot(opacity=0.8)
            sage: e = Hyperplane([4,0,2,3])
            sage: e.plot(ranges=[[-1,1],[0,8]], label_offset=(2,2,1), aspect_ratio=1)
            sage: opts = {'hyperplane_label':True, 'label_color':'green',
            ....: 'label_fontsize':24, 'label_offset':(0,1.5)}
            sage: Hyperplane([3,4,5]).plot(**opts)

        NOTES::

            For more examples, see :meth:`show` (``self.show?`` or
            ``Hyperplane.show?``).
        """
        if self.base_field() != QQ:
            raise NotImplementedError('Field must be QQ')
        elif self.dim() not in [0,1,2]: # dimension of self, not ambient space
            return # silently
        # handle extra keywords
        if kwds.has_key('hyperplane_label'):
            hyp_label = kwds.pop('hyperplane_label')
            if hyp_label == False:
                has_hyp_label = False
            else:
                has_hyp_label = True
        else: # default
            hyp_label = True
            has_hyp_label = True
        if has_hyp_label:
            if hyp_label == True: # then label hyperplane with its equation
                if self.dim() == 2: # jmol does not like latex
                    label = self._pretty_print_equation(latex=False)
                else:
                    label = self._pretty_print_equation()
            else:
                label = hyp_label # a string
        if kwds.has_key('label_color'):
            label_color = kwds.pop('label_color')
        else:
            label_color = 'black'
        if kwds.has_key('label_fontsize'):
            label_fontsize = kwds.pop('label_fontsize')
        else:
            label_fontsize = 14
        if kwds.has_key('label_offset'):
            has_offset = True
            label_offset = kwds.pop('label_offset')
        else:
            has_offset = False # give default values below
        if kwds.has_key('point_size'):
            pt_size = kwds.pop('point_size')
        else:
            pt_size = 50
        if kwds.has_key('ranges'):
            ranges_set = True
            ranges = kwds.pop('ranges')
        else:
            ranges_set = False # give default values below
        # the extra keywords have now been handled
        # now create the plot
        if self.dim() == 0: # a point on a line
            x, d = self.equation(True)
            p = point((d/x,0), size = pt_size, **kwds)
            if has_hyp_label:
                if not has_offset:
                    label_offset = 0.1
                p += text(label, (d/x,label_offset),
                        color=label_color,fontsize=label_fontsize)
                p += text('',(d/x,label_offset+0.4)) # add space at top
            if not kwds.has_key('ymax'):
                kwds['ymax'] = 0.5
        elif self.dim() == 1: # a line in the plane
            pnt = self.point()
            w = self.linear_part().matrix()
            x, y, d = self.equation(True)
            t = SR.var('t')
            if ranges_set:
                if type(ranges) in [list,tuple]:
                    t0, t1 = ranges
                else:  # ranges should be a single positive number
                    t0, t1 = -ranges, ranges
            else: # default
                t0, t1 = -3, 3
            p = parametric_plot(pnt+t*w[0], (t,t0,t1), **kwds)
            if has_hyp_label:
                if has_offset:
                    b0, b1 = label_offset
                else:
                    b0, b1 = 0, 0.2
                label = text(label,(pnt[0]+b0,pnt[1]+b1),
                        color=label_color,fontsize=label_fontsize)
                p += label
        elif self.dim() == 2: # a plane in 3-space
            pnt = self.point()
            w = self.linear_part().matrix()
            a, b, c, d = self.equation(True)
            s,t = SR.var('s t')
            if ranges_set:
                if type(ranges) in [list,tuple]:
                    s0, s1 = ranges[0]
                    t0, t1 = ranges[1]
                else: # ranges should be a single positive integers
                    s0, s1 = -ranges, ranges
                    t0, t1 = -ranges, ranges
            else: # default
                s0, s1 = -3, 3
                t0, t1 = -3, 3
            p = parametric_plot3d(pnt+s*w[0]+t*w[1],(s,s0,s1),(t,t0,t1),**kwds)
            if has_hyp_label: 
                if has_offset:
                    b0, b1, b2 = label_offset
                else:
                    b0, b1, b2 = 0, 0, 0
                label = text3d(label,(pnt[0]+b0,pnt[1]+b1,pnt[2]+b2),
                        color=label_color,fontsize=label_fontsize)
                p += label
        return p


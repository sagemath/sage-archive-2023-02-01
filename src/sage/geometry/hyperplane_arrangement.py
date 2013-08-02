r"""
Hyperplane arrangements

AUTHORS:

- David Perkinson (2013-06): initial version

- Qiaoyu Yang (2013-07)

- Kuai Yu (2013-07)

VERSION: 0.5

This module implements hyperplane arrangements defined over the
rationals or over finite fields.  The original motivation was to make a
companion to Richard Stanley's notes on hyperplane arrangements.

REFERENCES::

.. [RS] Stanley, Richard
   "Hyperplane Arrangements"
    Geometric Combinatorics (E. Miller, V. Reiner, and B. Sturmfels, eds.),
    IAS/Park City Mathematics Series, vol. 13, American Mathematical Society,
    Providence, RI, 2007, pp. 389-496.

EXAMPLES:

HYPERPLANES
Create a single hyperplane with equation `3x + 2y - 5z = 7`::

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

ARRANGEMENTS
There are several ways to create hyperplane arrangements: 
(i) a list of equations::

    sage: a = HyperplaneArrangement([[1,0,1],[1,1,0],[0,1,2]])
    sage: a.show()  # an arrangement of lines in the plane
    sage: a.show(axes=false)
    sage: a = HyperplaneArrangement([[2,4,6],[1,0,1]])
    sage: b = HyperplaneArrangement([[1,2,3],[2,0,2]])
    sage: a == b
    True

The default base field is QQ, the rational numbers.  Finite fields are also
supported::

    sage: a = HyperplaneArrangement([[1,2,3,4],[5,6,7,8]],GF(5))
    sage: a.hyperplanes()
    [Hyperplane over Finite Field of size 5
    (1, 2, 3)*X = 4,
     Hyperplane over Finite Field of size 5
    (0, 1, 2)*X = 3]

(ii) a list of hyperplanes::

    sage: k = [Hyperplane([1,0,0,i]) for i in range(4)]
    sage: a = HyperplaneArrangement(k)
    sage: a.show()

(iii) using the library of arrangements::

    sage: a = hyperplane_arrangements.braid(4)
    sage: a.show()
    Displaying the essentialization.
    sage: hyperplane_arrangements.semiorder(3)
    Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 3, rank 2.
    sage: hyperplane_arrangements.graphical(graphs.PetersenGraph())
    Hyperplane arrangement of 15 hyperplanes over Rational Field of dimension 10, rank 9.
    sage: hyperplane_arrangements.Ish(5)
    Hyperplane arrangement of 20 hyperplanes over Rational Field of dimension 5, rank 4.

(iv) from the bounding hyperplanes of a polyhedron::

    sage: a = HyperplaneArrangement(polytopes.n_cube(3))
    sage: a
    Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 3, rank 3.
    sage: a.show(hyperplane_legend=false)
    sage: a.num_regions()
    27

New arrangements from old::

    sage: a = hyperplane_arrangements.braid(3)
    sage: b = a.add_hyperplane(Hyperplane([1,2,3,4]))
    sage: b
    Hyperplane arrangement of 4 hyperplanes over Rational Field of dimension 3, rank 3.
    sage: c = b.deletion(Hyperplane([1,2,3,4]))
    sage: a == c
    True
    sage: a = hyperplane_arrangements.braid(3)
    sage: b = a.union(hyperplane_arrangements.semiorder(3))
    sage: b == hyperplane_arrangements.Catalan(3)
    True
    sage: a
    Hyperplane arrangement of 3 hyperplanes over Rational Field of dimension 3, rank 2.
    sage: a = hyperplane_arrangements.coordinate(4)
    sage: h = a.hyperplanes()[0]
    sage: b = a.restriction(h)
    sage: b == hyperplane_arrangements.coordinate(3)
    True
    sage: a = HyperplaneArrangement([[1,1],[1,3],[1,6]])
    sage: a.show()
    sage: a.cone().show(axes=false)

A hyperplane arrangement is *essential* is the normals to its hyperplane span
the ambient space.  Otherwise, it is *inessential*.  The essentialization is
formed by intersecting the hyperplanes by this normal space (actually, it is a
bit more complicated over finite fields)::


    sage: a = hyperplane_arrangements.braid(4)
    sage: a
    Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 4, rank 3.
    sage: a.is_essential()
    False
    sage: a.rank() < a.dim()  # double-check
    True
    sage: a.essentialization()
    Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 3, rank 3.

The connected components of the complement of the hyperplanes of an arrangement
in `\mathbb{R}^n` are called the *regions* of the arrangement::

    sage: a = hyperplane_arrangements.semiorder(3)
    sage: b = a.essentialization()
    sage: b
    Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 2, rank 2.
    sage: b.num_regions()
    19
    sage: b.regions()
    [A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays]
    sage: b.bounded_regions()
    [A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
     A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]
    sage: b.num_bounded_regions()
    7
    sage: a.unbounded_regions()
    [A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
     A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line]
    sage: b.region_containing_point((0,0)).show()
    sage: b.region_containing_point((2,1)).show(xmax=4,ymax=4)
    sage: r1 = b.regions()[0]
    sage: r2 = b.regions()[-1]  
    sage: b.distance_between_regions(r1,r2)  # number of hyps. separating r1, r2
    6
    sage: b.distance_enumerator(r1)  # generating fnc. for distances from r1
    x^6 + 2*x^5 + 5*x^4 + 3*x^3 + 5*x^2 + 2*x + 1

Note: *bounded region* really mean *relatively bounded* here.  A region is
relatively bounded if its intersection with space spanned by the normals of the
hyperplanes in the arrangement is bounded.

The intersection poset of a hyperplane arrangement is the collection of all
nonempty intersections of hyperplanes in the arrangement, ordered by reverse
inclusion.  It includes the ambient space of the arrangement (as the
intersection over the empty set)::

    sage: a = hyperplane_arrangements.braid(3)
    sage: p = a.intersection_poset()
    sage: p.is_ranked()
    True
    sage: p.show()
    sage: p.order_polytope()
    A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 10 vertices

The characteristic polynomial is a basic invariant of a hyperplane arrangement.
It is defined as `\chi(x) := \sum_{w\in P \mu(w) x^{\dim(w)}` where the sum is
`P` is the intersection poset of the arrangement and `\mu` is the Moebius
function of `P`::

    sage: a = hyperplane_arrangements.semiorder(5)
    sage: a.characteristic_polynomial()
    x^5 - 20*x^4 + 180*x^3 - 790*x^2 + 1380*x
    sage: a.poincare_polynomial()
    1380*x^4 + 790*x^3 + 180*x^2 + 20*x + 1
    sage: a.num_regions()
    2371
    sage: a.characteristic_polynomial(-1)
    -2371
    sage: a.num_bounded_regions()
    751
    sage: a.characteristic_polynomial(1)
    751

For finer invariants derived from the intersection poset, see
:meth:`whitney_number` and :meth:`doubly_indexed_whitney_number`.

Miscellaneous methods (see documentation for an explanation)::

    sage: a = hyperplane_arrangements.semiorder(3)
    sage: a.has_good_reduction(5)
    True
    sage: b = a.change_base_field(GF(5))
    sage: pa = a.intersection_poset()
    sage: pb = b.intersection_poset()
    sage: pa.is_isomorphic(pb)
    True
    sage: a.face_vector()
    (0, 12, 30, 19)
    sage: a.face_vector()
    (0, 12, 30, 19)
    sage: a.is_central()
    False
    sage: a.is_linear()
    False
    sage: a.sign_vector((1,1,1))
    [1, -1, 1, -1, 1, -1]
    sage: a.varchenko_matrix()
    [          1          h0       h0*h1       h0*h2    h0*h1*h2 h0*h1*h2*h3]
    [         h0           1          h1          h2       h1*h2    h1*h2*h3]
    [      h0*h1          h1           1       h1*h2          h2       h2*h3]
    [      h0*h2          h2       h1*h2           1          h1       h1*h3]
    [   h0*h1*h2       h1*h2          h2          h1           1          h3]
    [h0*h1*h2*h3    h1*h2*h3       h2*h3       h1*h3          h3           1]

There are extensive methods for visualizing hyperplane arrangements in low
dimensions.  See :meth:`plot` and :meth:`show` for details.
"""
#*****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Possible extensions for hyperplane_arrangement.py:
# - the big face lattice
# - Orlik-Solomon algebras
# - create ties with the Sage matroid methods
# - hyperplane arrangements over other fields

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

class HyperplaneArrangement(SageObject):

    def __init__(self, A, K=None, check_for_duplicates=True):
        r"""
        The argument ``A`` is normally a list of hyperplanes or a list of lists
        representing hyperplanes. If the latter, each hyperplane is given as a
        list `[a_1, \dots, a_n, a]` representing the equation `a_1 x_1 + \dots +
        a_n x_n = a`. An optional field, ``K``, may also be provided.
        Arrangements with multiple copies of a hyperplane are not supported, and
        a check is performed by default to remove duplicate hyperplanes. It is
        also possible for ``A`` to be a polyhedron, i.e., the intersection of a
        finite number of half-planes.  In that case, the hyperplanes are those
        bounding ``A``.

        INPUT:

        - ``A`` -- list of hyperplanes, list of lists of integers, or polyhedron
        - ``K`` -- field (defaults to ``QQ`` if none given and not inherited from
                   hyperplanes in A)
        - check_for_duplicates -- boolean

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLE::

            sage: A = HyperplaneArrangement(([1,0,0],[2,1,5]))
            sage: A
            Hyperplane arrangement of 2 hyperplanes over Rational Field of dimension 2, rank 2.
            sage: G = graphs.CycleGraph(4)
            sage: B = hyperplane_arrangements.graphical(G)
            sage: B.is_essential()
            False
            sage: C = B.essentialization()
            sage: B.characteristic_polynomial()
            x^4 - 4*x^3 + 6*x^2 - 3*x
            sage: G.chromatic_polynomial()
            x^4 - 4*x^3 + 6*x^2 - 3*x
            sage: C.characteristic_polynomial()
            x^3 - 4*x^2 + 6*x - 3
            sage: B.num_regions()
            14
            sage: C.num_regions()
            14
            sage: C.num_bounded_regions()
            0
            sage: D = hyperplane_arrangements.semiorder(3)
            sage: D.is_essential()
            False
            sage: D.num_regions()
            19
            sage: D.num_bounded_regions()
            7
            sage: D.face_vector()
            (0, 12, 30, 19)
            sage: D.doubly_indexed_whitney_number(1,2)
            -24
            sage: p = polytopes.n_simplex(4)
            sage: a = HyperplaneArrangement(p)
            sage: a.num_regions()
            31
            sage: a = Hyperplane([1,3,2,4],GF(19))
            sage: b = Hyperplane([7,9,25,19],GF(19))
            sage: c = Hyperplane([11,19,13,29],GF(19))
            sage: A = HyperplaneArrangement([a,b,c])
            sage: A
            Hyperplane arrangement of 3 hyperplanes over Finite Field of size 19 of dimension 3, rank 3.
            sage: A.hyperplanes()
            [Hyperplane over Finite Field of size 19
            (1, 3, 2)*X = 4,
            Hyperplane over Finite Field of size 19
            (7, 9, 6)*X = 0,
            Hyperplane over Finite Field of size 19
            (11, 0, 13)*X = 10]
        """
        self._raw_data = A
        if hasattr(A,'Hrepresentation'):
            K = QQ
            self._from_polyhedron = True
            self._polyhedron = deepcopy(A)
            self._hrepresentation = A.Hrepresentation()
            eqns = [list(h[1:])+[h[0]] for h in A.Hrepresentation()]
        else:
            # determine the base field
            if K==None:
                fields = [h.base_field() for h in A if type(h)==Hyperplane]
                if len(set(fields))==0: # default field is QQ
                    K = QQ
                elif len(set(fields))==1:
                    K = fields[0]
                else:
                    raise TypeError('incompatible field types')
            self._from_polyhedron = False
            eqns = []
            for h in A:
                if type(h)==Hyperplane:
                    eqns.append(map(K,h.equation(True)))
                elif type(h) in [list,tuple]:
                    eqns.append(map(K,h))
                else:
                    raise TypeError('input must be a Hyperplane, list, or tuple')
        if check_for_duplicates:
            self._hyperplanes = []
            for h in eqns:
                H = Hyperplane(h,K)
                if not any (H==k for k in self._hyperplanes):
                    self._hyperplanes.append(H)
        else:
            self._hyperplanes = [Hyperplane(h,K) for h in eqns]
        self._base_field = K
        self._dim = len(eqns[0])-1 # we trust the user on this point
        self._ambient_space = VectorSpace(K,self._dim)

    def __getattr__(self, name):
        r"""
        Calculates the entered variable if it has not already been stored.

        INPUT:

        - string (representing a variable)

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)
            sage: a.__getattr__('_characteristic_polynomial')
        """
        if not name in self.__dict__:
            if name == '_bounded_regions':
                self._set_bounded_regions()
                return deepcopy(self.__dict__[name])
            elif name == '_characteristic_polynomial':
                self._set_characteristic_polynomial()
                return deepcopy(self.__dict__[name])
            elif name == '_essentialization':
                self._set_essentialization()
                return self.__dict__[name]
            elif name == '_intersection_poset':
                self._set_intersection_poset()
                return deepcopy(self.__dict__[name])
            elif name == '_is_central':
                self._set_is_central()
                return self.__dict__[name]
            elif name == '_poincare_polynomial':
                self._set_poincare_polynomial()
                return deepcopy(self.__dict__[name])
            elif name == '_rank':
                self._set_rank()
                return self.__dict__[name]
            elif name == '_regions':
                self._set_regions()
                return self.__dict__[name]
            elif name == '_whitney_numbers':
                self._set_whitney_numbers()
                return self.__dict__[name]
            else:
                raise AttributeError(name)

    def __repr__(self):
        r"""
        String representation for a HyperplaneArrangement.

        INPUT:

        - None

        OUTPUT:

        - string

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(3)
            sage: a.__repr__()
            'Hyperplane arrangement of 3 hyperplanes over Rational Field of dimension 3, rank 2.'
        """
        return "Hyperplane arrangement of "+str(len(self._hyperplanes))+ " hyperplanes over "+str(self._base_field)+" of dimension "+str(self._dim)+", rank "+str(self.rank())+"."

    def __eq__(self, other):
        r"""
        Tests whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- HyperplaneArrangement

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,2,3],[4,5,6]])
            sage: b = HyperplaneArrangement([[4,5,6],[2,4,6]])
            sage: c = HyperplaneArrangement([[1,2,3],[4,5,6]], GF(7))
            sage: a == b
            True
            sage: a == c
            False
            sage: d = a.change_base_field(GF(7))
            sage: c == d
            True
        """
        if self.base_field() != other.base_field():
            return False
        elif len(self.hyperplanes()) != len(other.hyperplanes()):
                return False
        else:
            return all([h in other.hyperplanes() for h in self.hyperplanes()])

    def __le__(self, other):
        r"""
        Tests whether each hyperplane in ``self`` is a hyperplane in ``other``.

        INPUT:

        - ``other`` -- HyperplaneArrangement

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,2,3],[4,5,6]])
            sage: b = a.add_hyperplane([7,8,9])
            sage: a <= b
            True
            sage: b <= a
            False
            sage: a <= a
            True
        """
        if self.base_field() != other.base_field():
            return False
        else:
            return all([h in other.hyperplanes() for h in self.hyperplanes()])

    def __lt__(self, other):
        r"""
        Tests whether the set of hyperplanes of ``self`` is a proper subset of the
        hyperplanes in ``other``.

        INPUT:

        - ``other`` -- HyperplaneArrangement

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,2,3],[4,5,6]])
            sage: b = a.add_hyperplane([7,8,9])
            sage: a < b
            True
            sage: b < a
            False
            sage: a < a
            False
        """
        if self <= other and not self == other:
            return True
        else:
            return False

    def __ge__(self, other):
        r"""
        Tests whether each hyperplane in ``other`` is a hyperplane in ``self``.

        INPUT:

        - ``other`` -- HyperplaneArrangement

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,2,3],[4,5,6]])
            sage: b = a.add_hyperplane([7,8,9])
            sage: b >= a
            True
            sage: a >= b
            False
            sage: a >= a
            True
        """
        if self.base_field() != other.base_field():
            return False
        else:
            return all([h in self.hyperplanes() for h in other.hyperplanes()])

    def __gt__(self, other):
        r"""
        Tests whether the set of hyperplanes of ``other`` is a proper subset of the
        hyperplanes in ``self``.

        INPUT:

        - ``other`` -- HyperplaneArrangement

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,2,3],[4,5,6]])
            sage: b = a.add_hyperplane([7,8,9])
            sage: b > a
            True
            sage: a > b
            False
            sage: a > a
            False
        """
        if self >= other and not self == other:
            return True
        else:
            return False

    def __contains__(self, H):
        r"""
        Tests whether the hyperplane is in ``self``.

        INPUT:

        - ``H`` -- Hyperplane or list representing a hyperplane

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,2,3],[4,5,6]])
            sage: h = Hyperplane([2,4,6])
            sage: h in a
            True
            sage: [8,10,12] in a
            True
            sage: [1,2,4] in a
            False
        """
        if type(H) == Hyperplane and H.base_field() != self.base_field():
            return False
        if type(H) in [list, tuple]:
            H = Hyperplane(H,self.base_field())
        return H in self.hyperplanes()

    def add_hyperplane(self, H):
        r"""
        Return a hyperplane arrangement that is the union of the hyperplanes in
        ``self`` and ``H``.

        INPUT:

        - ``H`` -- Hyperplane

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,2,3],[4,5,6]])
            sage: h = Hyperplane([7,8,9])
            sage: c = a.add_hyperplane(h)
            sage: c.hyperplanes()
            [Hyperplane over Rational Field
            (1, 2)*X = 3,
             Hyperplane over Rational Field
            (4, 5)*X = 6,
             Hyperplane over Rational Field
            (7, 8)*X = 9]
        """
        hyps = self.hyperplanes()
        hyps.append(H)
        return HyperplaneArrangement(hyps, self.base_field())

    def hyperplanes(self):
        r"""
        The hyperplanes of the arrangement.

        INPUT:

        - None

        OUTPUT:

        - list of Hyperplanes

        EXAMPLES::

            sage: H = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,-1],[1,-1,0],[1,1,0]])
            sage: H.hyperplanes()
            [Hyperplane over Rational Field
            (1, 0)*X = 0,
             Hyperplane over Rational Field
            (0, 1)*X = 1,
             Hyperplane over Rational Field
            (0, 1)*X = -1,
             Hyperplane over Rational Field
            (1, -1)*X = 0,
             Hyperplane over Rational Field
            (1, 1)*X = 0]
        """
        return copy(self._hyperplanes)

    def cone(self):
        r"""
        Returns the cone over the hyperplane arrangement.  Its equations consist
        of `[a_1,...,a_n,d,0]` for each `[a_1,...,a_n,d]` in the original
        arrangement and the equation `[0,...,0,1,0]`.

        INPUT:

        - None

        OUTPUT:

        - hyperplane arrangement

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: b = a.cone()
            sage: a.characteristic_polynomial()
            x^3 - 6*x^2 + 12*x
            sage: factor(b.characteristic_polynomial())
            (x - 1) * x * (x^2 - 6*x + 12)
            sage: a.hyperplanes()
            [Hyperplane over Rational Field
            (1, -1, 0)*X = -1,
             Hyperplane over Rational Field
            (1, -1, 0)*X = 1,
             Hyperplane over Rational Field
            (1, 0, -1)*X = -1,
             Hyperplane over Rational Field
            (1, 0, -1)*X = 1,
             Hyperplane over Rational Field
            (0, 1, -1)*X = -1,
             Hyperplane over Rational Field
            (0, 1, -1)*X = 1]
            sage: b.hyperplanes()
            [Hyperplane over Rational Field
            (1, -1, 0, 1)*X = 0,
             Hyperplane over Rational Field
            (1, -1, 0, -1)*X = 0,
             Hyperplane over Rational Field
            (1, 0, -1, 1)*X = 0,
             Hyperplane over Rational Field
            (1, 0, -1, -1)*X = 0,
             Hyperplane over Rational Field
            (0, 1, -1, 1)*X = 0,
             Hyperplane over Rational Field
            (0, 1, -1, -1)*X = 0,
             Hyperplane over Rational Field
            (0, 0, 0, 1)*X = 0]
        """
        hyps = []
        for h in self.hyperplanes():
            new = h.equation(suppress_printing=True)
            new[-1] = -new[-1]
            new.append(0)
            hyps.append(new)
        hyps.append([0]*self.dim()+[1, 0])
        return HyperplaneArrangement(hyps, self.base_field())

    def dim(self):
        r"""
        The dimension of the ambient space for the hyperplane arrangement.

        INPUT:

        - None

        OUTPUT:

        - integer

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,2,3,0],[4,5,6,-3]])
            sage: A.dim()
            3
            sage: a = hyperplane_arrangements.Shi(4)
            sage: a.dim()
            4
            sage: p = polytopes.cyclic_polytope(7,5)
            sage: H = HyperplaneArrangement(p)
            sage: H.dim()
            7
        """
        return self._dim

    def _set_rank(self):
        r"""
        Sets the rank of the arrangement.

        See the documentation for :meth:`rank`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a._set_rank()
        """
        self._rank = matrix([list(H.normal()) for H in self._hyperplanes]).rank()

    def rank(self):
        r"""
        The dimension of the span of the normals to the hyperplanes in the
        arrangement.

        INPUT:

        - None

        OUTPUT:

        - integer

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,2,3,0],[4,5,6,-3]])
            sage: A.dim()
            3
            sage: A.rank()
            2
            sage: B = hyperplane_arrangements.braid(3)
            sage: B.hyperplanes()
            [Hyperplane over Rational Field
            (1, -1, 0)*X = 0,
             Hyperplane over Rational Field
            (1, 0, -1)*X = 0,
             Hyperplane over Rational Field
            (0, 1, -1)*X = 0]
            sage: B.dim()
            3
            sage: B.rank()
            2
            sage: p = polytopes.n_simplex(5)
            sage: H = HyperplaneArrangement(p)
            sage: H.rank()
            5
        """
        return self._rank

    def has_good_reduction(self, p):
        r"""
        Let ``A`` be a hyperplane arrangement with equations defined over the
        integers, and let ``B`` be the hyperplane arrangement defined by
        reducing these equations modulo a prime ``p``.  Then ``A`` has good
        reduction modulo ``p`` if the intersection posets of ``A`` and ``B`` are
        isomorphic.

        INPUT:

        - ``p`` -- prime number

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a.has_good_reduction(5)
            True
            sage: a.has_good_reduction(3)
            False
            sage: b = a.change_base_field(GF(3))
            sage: a.characteristic_polynomial()
            x^3 - 6*x^2 + 12*x
            sage: b.characteristic_polynomial()  # not equal to that for a
            x^3 - 6*x^2 + 10*x
        """
        if self.base_field()!=QQ:
            raise TypeError('Arrangement must be defined over QQ')
        elif not p.is_prime:
            raise TypeError('Must reduce modulo a prime number.')
        else:
            a = self.change_base_field(GF(p))
            p = self.intersection_poset()
            q = a.intersection_poset()
            return p.is_isomorphic(q)

    def _set_is_central(self):
        r"""
        Is the intersection of all the hyperplanes in the arrangement nonempty?

        See the documentation for :meth:`is_central`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)
            sage: a._set_is_central()
        """
        m = matrix(self.base_field(),[i.normal() for i in self.hyperplanes()])
        b = vector(self.base_field(),[i.equation(True)[-1] for i in
            self.hyperplanes()])
        try:
            x = m.solve_right(b)
            self._is_central =  True
        except ValueError:
            self._is_central = False

    def is_central(self):
        r"""
        Is the intersection of all the hyperplanes in the arrangement nonempty?

        INPUT:

        - None

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: s = hyperplane_arrangements.semiorder(3)
            sage: s.is_central()
            False
            sage: b = hyperplane_arrangements.braid(3)
            sage: b.is_central()
            True
        """
        return self._is_central

    def is_essential(self):
        r"""
        Indicates whether the hyperplane arrangement is essential. A hyperplane
        arrangement is essential if the span of the normals of its hyperplanes
        spans the ambient space.

        INPUT:

        - None

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: HyperplaneArrangement([[1,0,0],[1,0,1]]).is_essential()
            False
            sage: HyperplaneArrangement([[1,0],[2,0]]).is_essential()
            True

        .. SEEALSO::

            :meth:`essentialization`
        """
        r = matrix([list(H.normal()) for H in self._hyperplanes]).rank()
        return r == self._dim

    def is_linear(self):
        r"""
        Do all the hyperplanes pass through the origin?

        INPUT:

        - None

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a.is_linear()
            False
            sage: b = hyperplane_arrangements.braid(3)
            sage: b.is_linear()
            True
            sage: c = HyperplaneArrangement([[1,0,1],[0,1,1]])
            sage: c.is_linear()
            False
            sage: c.is_central()
            True
        """
        return all(i.equation(True)[-1]==0 for i in self.hyperplanes())

    def _set_intersection_poset(self):
        r"""
        Calculates the intersection poset of the hyperplane arrangement.

        See the documentation for :meth:`intersection_poset`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a._set_intersection_poset()
        """
        # initiate with codim = 0, K^n
        K = self.base_field()
        Kn = AffineSubspace(vector(K, [0 for i in range(self._dim)]),
                      self._ambient_space)
        L = [[Kn]]
        active = True
        codim = 0
        while active:
            active = False
            new_level = []
            for T in L[codim]:
                for H in self._hyperplanes:
                    I = H.intersection(T)
                    if type(I) != type(-1) and I != T and I not in new_level:
                        new_level.append(I)
                        active = True
            if active:
                L.append(new_level)
            codim += 1
        L = flatten(L)
        t = {}
        for i in range(len(L)):
            t[i] = L[i]
        cmp_fn = lambda p, q: t[q] < t[p]
        self._intersection_poset = Poset((t, cmp_fn))

    def intersection_poset(self):
        """
        Return the intersection poset of the arrangement.

        OUTPUT:

        - a poset

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.intersection_poset()
            Finite poset containing 19 elements
        """
        return deepcopy(self._intersection_poset)

    def _set_characteristic_polynomial(self):
        r"""
        Calculates characteristic polynomial of the hyperplane arrangement.

        See the documentation for :meth:`characteristic_polynomial`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a._set_characteristic_polynomial()
        """
        #P = self._intersection_poset
        #n = self._dim
        #self._characteristic_polynomial = sum([P.mobius_function(0, p)*x^(n-P.rank(p)) for p in P])
        x = polygen(QQ, 'x')
        if self.rank() == 1:
            self._characteristic_polynomial = x**(self.dim()-1)*(x-self.num_hyperplanes())
        else:
            H = self.hyperplanes()[0]
            R = self.restriction(H)
            r = R.characteristic_polynomial()
            D = self.deletion(H)
            d = D.characteristic_polynomial()
            self._characteristic_polynomial = expand(d - r)
            return

    def _set_poincare_polynomial(self):
        r"""
        Calculates Poincare polynomial of the hyperplane arrangement.

        See the documentation for :meth:`poincare_polynomial`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a._set_poincare_polynomial()
        """
        x = polygen(QQ, 'x')
        p = (-x)**self.dim()*self.characteristic_polynomial(-ZZ(1)/x)
        self._poincare_polynomial = p

    def change_base_field(self, F):
        r"""
        Returns a hyperplane arrangement defined over F.

        INPUT:

        - ``F`` -- ``Q`` or a finite field

        OUTPUT:

        - hyperplane arrangement defined over ``F``

        EXAMPLES::

            sage: a = HyperplaneArrangement([[1,3,7],[6,2,14]])
            sage: b = a.change_base_field(GF(5))
            sage: c = b.change_base_field(GF(3))
            sage: d = c.change_base_field(QQ)
            sage: a.hyperplanes()
            [Hyperplane over Rational Field
            (1, 3)*X = 7,
             Hyperplane over Rational Field
            (6, 2)*X = 14]
            sage: b.hyperplanes()
            [Hyperplane over Finite Field of size 5
            (1, 3)*X = 2,
             Hyperplane over Finite Field of size 5
            (1, 2)*X = 4]
            sage: c.hyperplanes()
            [Hyperplane over Finite Field of size 3
            (1, 0)*X = 2,
             Hyperplane over Finite Field of size 3
            (1, 2)*X = 1]
            sage: d.hyperplanes()
            [Hyperplane over Rational Field
            (1, 0)*X = 2,
             Hyperplane over Rational Field
            (1, 2)*X = 1]
            sage: a == d
            False
        """
        eqs = [i.equation(True) for i in self.hyperplanes()]
        eqs = [map(F, i) for i in eqs]
        return HyperplaneArrangement(eqs, F)

    def characteristic_polynomial(self, a=None):
        r"""
        Returns the characteristic polynomial of the given hyperplane
        arrangement: `\chi(x) := \sum_{w\in P \mu(w) x^{\dim(w)}` where the sum
        is `P` is the intersection poset of the arrangement and `\mu` is the
        Moebius function of `P`.

        If the argument ``a`` is given, the characteristic polynomial
        is evaluated at ``a``.

        INPUT:

        - ``a`` -- element of base field (default: None)

        OUTPUT:

        - polynomial or base field element

        EXAMPLES::

            sage: A = HyperplaneArrangement(([1,0,0],[2,1,5],[3,1,1]))
            sage: A.characteristic_polynomial()
            x^2 - 3*x + 3
            sage: A.characteristic_polynomial(-1)
            7
            sage: A.num_regions()
            7
        """
        if a is None:
            return self._characteristic_polynomial
        else:
            return self._characteristic_polynomial.subs(x=a)

    def poincare_polynomial(self, a=None):
        r"""
        Returns the Poincare polynomial of the given hyperplane arrangement.

        If the argument ``a`` is given, the characteristic polynomial
        is evaluated at ``a``.

        INPUT:

        - ``a`` -- element of base field (default: None)

        OUTPUT:

        - polynomial or base field element

        EXAMPLES::

            sage: A = HyperplaneArrangement(([1,-1,3],[2,1,5],[3,2,0]))
            sage: A.poincare_polynomial()
            3*x^2 + 3*x + 1
            sage: A.characteristic_polynomial()
            x^2 - 3*x + 3
            sage: A.poincare_polynomial(1)
            7
        """
        if a is None:
            return self._poincare_polynomial
        else:
            return self._poincare_polynomial(a)

    def num_bounded_regions(self):
        r"""
        Returns the number of (relatively) bounded regions of the hyperplane
        arrangement as a subset of `\mathbb{RR}`.

        INPUT:

        - None

        OUTPUT:

        - integer

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.num_bounded_regions()
            7

        TESTS::

            sage: A = HyperplaneArrangement([[1,1,0],[2,3,-1],[4,5,3]])
            sage: B = A.change_base_field(FiniteField(7))
            sage: B.num_bounded_regions()
            Traceback (most recent call last):
            ...
            TypeError: base field must have characteristic zero
        """
        if self.base_field().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        return (-1)**self.rank()*self.characteristic_polynomial(1)

    def num_regions(self):
        r"""
        The number of regions of the hyperplane arrangement.

        INPUT:

        - None

        OUTPUT:

        - integer

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.num_regions()
            19

        TESTS::

            sage: A = HyperplaneArrangement([[1,1,0],[2,3,-1],[4,5,3]])
            sage: B = A.change_base_field(FiniteField(7))
            sage: B.num_regions()
            Traceback (most recent call last):
            ...
            TypeError: base field must have characteristic zero
        """
        if self.base_field().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        return (-1)**self._dim*self.characteristic_polynomial(-1)

    def num_hyperplanes(self):
        r"""
        The number of hyperplanes in the arrangement.

        INPUT:

        - None

        OUTPUT:

        - integer

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,1,0],[2,3,-1],[4,5,3]])
            sage: A.num_hyperplanes()
            3
        """
        return len(self.hyperplanes())

    def ambient_space(self):
        r"""
        The ambient space of the hyperplane arrangement.

        INPUT:

        - None

        OUTPUT:

        - vector space

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,-1],[1,-1,0],[1,1,0]])
            sage: A.ambient_space()
            Vector space of dimension 2 over Rational Field
            sage: A.change_base_field(GF(7)).ambient_space()
            Vector space of dimension 2 over Finite Field of size 7
        """
        return deepcopy(self._ambient_space)

    def base_field(self):
        r"""
        The base field of the hyperplane arrangement.

        INPUT:

        - None

        OUTPUT:

        - field

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,-1],[1,-1,0],[1,1,0]])
            sage: A.base_field()
            Rational Field
            sage: B = HyperplaneArrangement([[1,0,0],[0,1,1]],GF(5))
            sage: B.base_field()
            Finite Field of size 5

        .. SEEALSO::

            :meth:`change_base_field`
        """
        return self._base_field

    def deletion(self, h):
        r"""
        The hyperplane arrangement obtained by removing ``h`` from the
        arrangement.

        INPUT:

        - ``h`` -- Hyperplane

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,-1],[1,-1,0],[1,1,0]])
            sage: h = Hyperplane([1,0,0])
            sage: A
            Hyperplane arrangement of 5 hyperplanes over Rational Field of dimension 2, rank 2.
            sage: A.deletion(h)
            Hyperplane arrangement of 4 hyperplanes over Rational Field of dimension 2, rank 2.

        TESTS::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,-1],[1,-1,0],[1,1,0]])
            sage: h = Hyperplane([4,0,0])
            sage: A.deletion(h)
            Hyperplane arrangement of 4 hyperplanes over Rational Field of dimension 2, rank 2.
            sage: l = Hyperplane([1,2,3])
            sage: A.deletion(l)
            Traceback (most recent call last):
            ...
            TypeError: hyperplane is not in the arrangement

        .. SEEALSO::

            :meth:`restriction`
        """
        if h in self:
            b = self.hyperplanes()
            b.remove(h)
            return HyperplaneArrangement(b, self.base_field())
        else:
            raise TypeError('hyperplane is not in the arrangement')

    def _set_essentialization(self):
        r"""
        Calculates essentialization of the hyperplane arrangement.

        See the documentation for :meth:`essentialization`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(3)
            sage: a.is_essential()
            False
            sage: a._set_essentialization()

        TESTS::

            sage: b = hyperplane_arrangements.coordinate(2)
            sage: b.is_essential()
            True
            sage: b._set_essentialization()
            Traceback (most recent call last):
            ...
            UserWarning: already essential
        """
        if self.rank() == self.dim():
            raise UserWarning('already essential')
            self._essentialization = HyperplaneArrangement(self.hyperplanes())
        else:
            n = [h.normal() for h in self.hyperplanes()]
            N = self.ambient_space().subspace(n)
            # We need to be careful finding complementary spaces when the
            # characteristic is not 0.
            if N.base_field().characteristic() != 0:
                m = N.echelonized_basis_matrix()
                c = [m.nonzero_positions_in_row(i)[0] for i in range(m.nrows())]
                new_basis = [self.ambient_space().basis()[i] for i in c]
                N = self.ambient_space().subspace(new_basis)
            N = AffineSubspace(vector([0]*self.dim()), N)
            new_hyperplanes = []
            for h in self.hyperplanes():
                I = h.intersection(N)
                basis = I.linear_part().matrix()
                p = N._isomorphism_with_Kn(I.point())
                new_basis = matrix(self.base_field(),
                                   [N._isomorphism_with_Kn(b) for b in basis])
                if new_basis.nrows() == 0:
                    eqn = [1, p[0]]
                else:
                    eqn = new_basis.right_kernel_matrix()[0]
                    q = [eqn*p]
                    eqn = list(eqn)+q
                new_hyperplanes.append(eqn)
            self._essentialization = HyperplaneArrangement(new_hyperplanes,
                                                           self.base_field())

    def essentialization(self):
        r"""
        The essentialization of a hyperplane arrangement.

        If the characteristic of the base field is 0, this returns the
        hyperplane arrangement obtained by intersecting the
        hyperplanes by the space spanned by their normal vectors.

        INPUT:

        - None

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: A = hyperplane_arrangements.braid(3)
            sage: A.is_essential()
            False
            sage: A.essentialization()
            Hyperplane arrangement of 3 hyperplanes over Rational Field of dimension 2, rank 2.
            sage: B = HyperplaneArrangement([[1,0,1],[1,0,-1]])
            sage: B.is_essential()
            False
            sage: B.essentialization()
            Hyperplane arrangement of 2 hyperplanes over Rational Field of dimension 1, rank 1.
            sage: B.essentialization().hyperplanes()
            [Hyperplane over Rational Field
            (1)*X = 1,
            Hyperplane over Rational Field
            (1)*X = -1]
            sage: C = HyperplaneArrangement([[1,1,1],[1,1,0]],GF(2))
            sage: C.essentialization().hyperplanes()
            [Hyperplane over Finite Field of size 2
            (1)*X = 1,
            Hyperplane over Finite Field of size 2
            (1)*X = 0]
        """
        return self._essentialization

    def polyhedron(self):
        r"""
        A polyhedron is the intersection of a finite number of half-planes.

        INPUT:

        - None

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: a = HyperplaneArrangement(p)
            sage: a.num_regions()
            59
            sage: a.polyhedron()
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: b = HyperplaneArrangement([[2, 3, 1], [1, 1, 0]])
            sage: b.polyhedron()  # returns None since b was not defined using a polyhedron
        """
        if self._from_polyhedron:
            return deepcopy(self._polyhedron)

    def _set_regions(self):
        r"""
        Calculates the regions of the hyperplane arrangement.

        See the documentation for :meth:`regions`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)
            sage: a._set_regions()
        """
        if self.base_field().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        result = []
        num = self.num_hyperplanes()
        h = [i.equation(suppress_printing=True) for i in self.hyperplanes()]
        for pos in powerset(range(num)):
            q = []
            for i in range(num):
                if i in pos:
                    new = [-h[i][-1]]+h[i][:-1]
                    q.append(new)
                else:
                    new = [h[i][-1]]+[-t for t in h[i][:-1]]
                    q.append(new)
            P = Polyhedron(ieqs=q)
            if P.dim() == self.dim():
                result.append(P)
        self._regions = result

    def regions(self):
        r"""
        The regions of the hyperplane arrangement

        The regions are the connected components of the complement of
        the union of the hyperplanes as a subset of `\mathbb{R}^n`.

        INPUT:

        - None

        OUTPUT:

        - list of polyhedra

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1]])
            sage: A.regions()
            [A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays]
        """
        return self._regions

    def region_containing_point(self, p):
        r"""
        The region in the hyperplane arrangement containing a given point.  It
        is assumed that the arrangement is defined over the reals.

        INPUT:

        - ``p`` -- point

        OUTPUT:

        - polyhedron

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,-1],[1,-1,0],[1,1,0]])
            sage: A.region_containing_point((1,2))
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 2 rays

        TESTS::

            sage: A = HyperplaneArrangement([[1,1,0],[2,3,-1],[4,5,3]])
            sage: B = A.change_base_field(FiniteField(7))
            sage: B.region_containing_point((1,2))
            Traceback (most recent call last):
            ...
            TypeError: base field must have characteristic zero

            sage: A = HyperplaneArrangement([[1,1,0],[2,3,-1],[4,5,3]])
            sage: A.region_containing_point((1,-1))
            Traceback (most recent call last):
            ...
            UserWarning: point sits on a hyperplane
        """
        if self.base_field().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        h = [i.equation(suppress_printing=True) for i in self.hyperplanes()]
        q = []
        sign_vector = self.sign_vector(p)
        if 0 in sign_vector:
            raise UserWarning('point sits on a hyperplane')
            return
        for i in range(self.num_hyperplanes()):
            if sign_vector[i] == 1:
                new = [-h[i][-1]]+h[i][:-1]
                q.append(new)
            else:
                new = [h[i][-1]]+[-t for t in h[i][:-1]]
                q.append(new)
        return Polyhedron(ieqs=q)

    def repr_point(self, region, scale=1):
        r"""
        A representative point from the hyperplane ``region``.  It is assumed
        that the base field of the arrangment is ``QQ``.

        INPUT:

        - ``region`` -- polyhedron (region of the arrangement)
        - ``scale`` -- number (default: 1)

        OUTPUT:

        - vector

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: r = a.region_containing_point((1,1))
            sage: a.repr_point(r)
            (70711/100000, 70711/100000)

        NOTES::

        Some rounding occurs in order to stay in the rationals.  This might be a
        source of error.  Adjusting the scale might then help (see the source
        code.
        """
        c = region.center()
        # rounding off here, so this function might not be completely reliable
        c.apply_map(lambda x: QQ(round(RR(x),5)))
        if region.n_rays() == 0:
            return c
        else:
            r = sum([vector(QQ, i) for i in region.rays()])
            r = (scale*r/r.norm()).apply_map(lambda x: QQ(round(RR(x),5)))
            return c + r

    def _set_bounded_regions(self):
        r"""
        Find the relatively bounded regions of the arrangements.  This function
        returns the indices of these regions in self._regions.

        See the documentation for :meth:`bounded_regions`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a._set_bounded_regions()
        """
        if not self.is_essential():
            R = self.essentialization().regions()
        else:
            R = self.regions()
        b =[i for i in range(self.num_regions()) if R[i].is_compact()]
        self._bounded_regions = b

    def bounded_regions(self):
        r"""
        The relatively bounded regions of the arrangement.  A region is
        relatively bounded if its intersection with the space spanned by the
        normals to the hyperplanes is bounded.  This is the same as being
        bounded in the case that the hyperplane arrangement is essential.  It is
        assumed that the arrangement is defined over the rationals.

        INPUT:

        - None

        OUTPUT:

        - list of polyhedra

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.bounded_regions()
            [A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line,
             A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 1 line]
            sage: A.bounded_regions()[0].is_compact() # the regions are only *relatively* bounded
            False
            sage: A.is_essential()
            False

        .. SEEALSO::

            :meth:`unbounded_regions`
        """
        return [self.regions()[i] for i in self._bounded_regions]

    def unbounded_regions(self):
        r"""
        The regions of the arrangement that are not relatively bounded.  It is assumed that the
        arrangement is defined over the rationals.

        See the documentation for :meth:`bounded_regions`.

        INPUT:

        - None

        OUTPUT:

        - list of polyhedra

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: B = A.essentialization()
            sage: B.num_regions() - B.num_bounded_regions()
            12
            sage: B.unbounded_regions()
            [A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays]

        .. SEEALSO::

            :meth:`bounded_regions`
        """
        s = set(range(self.num_regions())).difference(set(self._bounded_regions))
        return [self.regions()[i] for i in s]

    def is_separating_hyperplane(self, region1, region2, hyp):
        r"""
        Does the hyperplane ``hyp`` separate the given regions?

        INPUT:

        - ``region1``, ``region2`` -- regions of the arrangement or representative points of regions
        - ``hyp`` -- hyperplane

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: A = hyperplane_arrangements.coordinate(2)
            sage: h = A.hyperplanes()[0]
            sage: h
            Hyperplane over Rational Field
            (1, 0)*X = 0
            sage: A.is_separating_hyperplane([1,1],[2,1],h)
            False
            sage: A.is_separating_hyperplane([1,1],[-1,1],h)
            True
            sage: r = A.region_containing_point([1,1])
            sage: s = A.region_containing_point([-1,1])
            sage: A.is_separating_hyperplane(r,s,h)
            True
        """
        if self.base_field() != QQ:
            raise TypeError('function requires base field QQ')
        if hasattr(region1, 'Hrepresentation'):
            p1 = list(self.repr_point(region1))
        else:
            p1 = list(region1)
        if hasattr(region2, 'Hrepresentation'):
            p2 = list(self.repr_point(region2))
        else:
            p2 = list(region2)
        p1.append(-1)
        p2.append(-1)
        p1 = vector(QQ, p1)
        p2 = vector(QQ, p2)
        e = vector(hyp.equation(True))
        s = sign(p1*e)*sign(p2*e)
        if s < 0:
            return True
        elif s > 0:
            return False
        else:
            raise UserWarning('point lies on hyperplane')

    def distance_between_regions(self, region1, region2):
        r"""
        Returns the number of hyperplanes separating the two regions.

        INPUT:

        - ``region1``, ``region2`` -- regions of the arrangement or representative points of regions

        OUTPUT:

        - integer

        EXAMPLES::

            sage: c = hyperplane_arrangements.coordinate(2)
            sage: r = c.region_containing_point([-1,-1])
            sage: s = c.region_containing_point([1,1])
            sage: c.distance_between_regions(r,s)
            2
            sage: c.distance_between_regions(s,s)
            0
        """
        r = [self.is_separating_hyperplane(region1,region2,h) for h in
                self.hyperplanes()]
        return r.count(True)

    def distance_enumerator(self, base_region):
        r"""
        Returns a polynomial `f(x)` for which the coefficient of `x^i` is the
        number of hyperplanes of distance `i` from ``base_region``, i.e., the 
        number of hyperplanes separated by `i` hyperplanes from ``base_region``.

        INPUT:

        - ``base_region`` -- region of arrangement or point in region

        OUTPUT:

        - polynomial

        EXAMPLES::

            sage: c = hyperplane_arrangements.coordinate(3)
            sage: c.distance_enumerator(c.region_containing_point([1,1,1]))
            x^3 + 3*x^2 + 3*x + 1
        """
        d = [self.distance_between_regions(r,base_region) for r in
                self.regions()]
        d = [d.count(i) for i in range(max(d)+1)]
        x = polygen(QQ, 'x')
        return sum([d[i]*x**i for i in range(len(d))])

    def sign_vector(self, p):
        r"""
        Sign_vector indicates on which side of each hyperplane the given point `p` lies.

        INPUT:

        - ``p`` -- point

        OUTPUT:

        - vector

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1]])
            sage: A.sign_vector((2,-2))
            [1, -1]
            sage: A.sign_vector((1,1))
            [1, 0]

        TESTS::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1]],GF(3))
            sage: A.sign_vector((1,2))
            Traceback (most recent call last):
            ...
            UserWarning: characteristic must be zero
        """
        if self.base_field().characteristic() != 0:
            raise UserWarning('characteristic must be zero')
        else:
            p = vector(p)
            eqns = [i.equation(suppress_printing=True)
                    for i in self.hyperplanes()]
            return [sign(p*vector(e[:-1])-e[-1]) for e in eqns]

    def face_vector(self):
        r"""
        The number of faces of each dimension: `face_vector(d)` is the number of
        faces of dimension `d`.  A *face* is is the intersection of a region
        with a hyperplane.  It is assumed that the arrangement is defined over
        the reals.

        INPUT:

        - None

        OUTPUT:

        - vector

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.face_vector()
            (0, 6, 21, 16)
        """
        m = self._whitney_numbers[0]
        v = list(sum(m.transpose().apply_map(abs)))
        v.reverse()
        v = [0]*(self.dim()-self.rank()) + v
        return vector(v)

    # Why not allow restricting to arbitrary affine subspaces?
    def restriction(self, H):
        r"""
        The restriction of the hyperplane arrangement to the given hyperplane, `H`.
        This is obtained by intersecting  `H` with all of the other distinct
        hyperplanes.  The hyperplane `H` must be in the arrangement.

        INPUT:

        - ``H`` -- Hyperplane

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: A = hyperplane_arrangements.braid(4)
            sage: H = A.hyperplanes()[0]
            sage: H
            Hyperplane over Rational Field
            (1, -1, 0, 0)*X = 0
            sage: R = A.restriction(H)
            sage: D = A.deletion(H)
            sage: ca = A.characteristic_polynomial()
            sage: cr = R.characteristic_polynomial()
            sage: cd = D.characteristic_polynomial()
            sage: ca
            x^4 - 6*x^3 + 11*x^2 - 6*x
            sage: cd - cr
            x^4 - 6*x^3 + 11*x^2 - 6*x

        .. SEEALSO::

            :meth:`deletion`
        """
        if H not in self.hyperplanes():
            raise UserWarning('hyperplane not in arrangement')
        else:
            new_hyperplanes = []
            for T in self.hyperplanes():
                if T != H:
                    I = H.intersection(T)
                    if I != -1 and I not in new_hyperplanes:
                        new_hyperplanes.append(I)
        final_hyperplanes = []
        for T in new_hyperplanes:
            basis = T.linear_part().matrix()
            p = H._isomorphism_with_Kn(T.point())
            if basis.nrows() == 0:
                eqn = [1] + list(p)
            else:
                new_basis = matrix(self.base_field(),
                                   [H._isomorphism_with_Kn(b+H.point())
                    for b in basis])
                eqn = new_basis.right_kernel_matrix()[0]
                q = eqn*p
                eqn = list(eqn)+[q]
            final_hyperplanes.append(eqn)
        return HyperplaneArrangement(final_hyperplanes, self.base_field())

    def show(self, **kwds):
        r"""
        Displays the hyperplane arrangement.

        INPUT:

        - **kwds -- show options: see below

        OUTPUT:

        - Graphics

        PLOT OPTIONS::

            Beside the usual show options (enter show?), the show command for
            hyperplanes includes the following:

            - hyperplane_colors -- Color or list of colors, one for each
              hyperplane (default: equally spread range of hues).

            - hyperplane_labels -- Boolean, 'short', 'long' (default: False).
              If False, no labels are shown; if 'short' or 'long', the
              hyperplanes are given short or long labels, respectively.  If
              ``True``, the hyperplanes are given long labels.

            - label_colors -- Color or list of colors, one for each hyperplane
              (default: black).

            - label_fontsize -- Size for hyperplane_label font (default: 14).
              This does not work for 3d plots.

            - label_offsets -- Amount be which labels are offset from h.point()
              for each hyperplane h.  The format is different for each
              dimension: if the hyperplanes have dimension 0, the offset can be
              a single number or a list of numbers, one for each hyperplane; if
              the hyperplanes have dimension 1, the offset can be a single
              2-tuple, or a list of 2-tuples, one for each hyperplane; if the
              hyperplanes have dimension 2, the offset can be a single 3-tuple
              or a list of 3-tuples, one for each hyperplane.  (Defaults:
              0-dim: 0.1, 1-dim: (0,1), 2-dim: (0,0,0.2)).

            - hyperplane_legend -- Boolean, 'short', 'long' (default: 'long').
              If ``False``, no legend is shown; if ``True``, 'short', or 'long',
              the legend is shown with the default, long, or short labeling,
              respectively. (For arrangements of lines or planes, only.)

            - hyperplane_opacities -- a number or list of numbers, one for each
              hyperplane, between 0 and 1.  Only applies to 3d plots.

            - point_sizes -- number or list of numbers, one for each hyperplane
              giving the sizes of points in a zero-dimensional arrangement
              (default: 50).

            - ranges -- Range for the parameters or a list of ranges of
              parameters, one for each hyperplane, for the parametric plots of
              the hyperplanes.  If a single positive number `r` is given for
              ``ranges``, then all parameters run from -r to r.  Otherwise, for
              a line in the plane, the range has the form [a,b] (default:
              [-3,3]), and for a plane in 3-space, the range has the form
              [[a,b],[c,d]] (default: [[-3,3],[-3,3]]). (The ranges are centered
              around self.point().)

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0,0],[0,0,1,5]])
            sage: A.show(hyperplane_opacities=0.5,hyperplane_labels=True,
            ....: hyperplane_legend=False, frame=False)
            sage: c = hyperplane_arrangements.Catalan(4)
            sage: c.show(frame=false)
            Displaying the essentialization.
            sage: pts = HyperplaneArrangement([[1,2],[5,2],[3,1]])
            sage: opts = {'point_sizes': 100}
            sage: opts['hyperplane_labels'] = True
            sage: opts['label_offsets'] = [0.05,-0.05,0.05]
            sage: pts.show(axes=false, xmin=-0.2, xmax=2.5, **opts)
            sage: opts['hyperplane_labels'] = 'short'
            sage: opts['hyperplane_legend'] = 'short'
            sage: pts.show(axes=false, xmin=-0.2, xmax=2.5, **opts)
            sage: br = hyperplane_arrangements.braid(3)
            sage: bre = br.essentialization()
            sage: bre.show(axes=false,hyperplane_legend='short')
            sage: br.hyperplanes()
            [Hyperplane over Rational Field
            (1, -1, 0)*X = 0,
             Hyperplane over Rational Field
            (1, 0, -1)*X = 0,
             Hyperplane over Rational Field
            (0, 1, -1)*X = 0]

        NOTES::

            For more examples, see :meth:`plot` (``self.plot?`` or
            ``HyperplaneArrangement.plot?).
        """
        result = self.plot(**kwds)
        for k in ['hyperplane_colors','hyperplane_labels','label_colors',
               'label_fontsize','label_offsets','hyperplane_legend',
               'hyperplane_opacities', 'point_sizes', 'ranges']:
            if kwds.has_key(k):
                del kwds[k]
        if self.dim() == 1 and not kwds.has_key('ymin'):
            kwds['ymin'] = -0.3
        if type(result) == tuple: # so dim=3 and hyperplane_legend!=False
            result[0].show(**kwds)
            result[1].show(**kwds)
        else:
            result.show(**kwds)

    def plot(self, **kwds):
        r"""
        Returns the plot of the hyperplane arrangement.  If the arrangement is
        in 4 dimensions but inessential, a plot of the essentialization is
        returned.

        INPUT:

        - **kwds -- plot options: see below

        OUTPUT:

        - Graphics

        PLOT OPTIONS::

            Beside the usual plot options (enter plot?), the plot command for
            hyperplanes includes the following:

            - hyperplane_colors -- Color or list of colors, one for each
              hyperplane (default: equally spread range of hues).

            - hyperplane_labels -- Boolean, 'short', 'long' (default: False).
              If False, no labels are shown; if 'short' or 'long', the
              hyperplanes are given short or long labels, respectively.  If
              ``True``, the hyperplanes are given long labels.

            - label_colors -- Color or list of colors, one for each hyperplane
              (default: black).

            - label_fontsize -- Size for hyperplane_label font (default: 14).
              This does not work for 3d plots.

            - label_offsets -- Amount be which labels are offset from h.point()
              for each hyperplane h.  The format is different for each
              dimension: if the hyperplanes have dimension 0, the offset can be
              a single number or a list of numbers, one for each hyperplane; if
              the hyperplanes have dimension 1, the offset can be a single
              2-tuple, or a list of 2-tuples, one for each hyperplane; if the
              hyperplanes have dimension 2, the offset can be a single 3-tuple
              or a list of 3-tuples, one for each hyperplane.  (Defaults:
              0-dim: 0.1, 1-dim: (0,1), 2-dim: (0,0,0.2)).

            - hyperplane_legend -- Boolean, 'short', 'long' (default: 'long').
              If ``False``, no legend is shown; if ``True``, 'short', or 'long',
              the legend is shown with the default, long, or short labeling,
              respectively. (For arrangements of lines or planes, only.)

            - hyperplane_opacities -- a number or list of numbers, one for each
              hyperplane, between 0 and 1.  Only applies to 3d plots.

            - point_sizes -- number or list of numbers, one for each hyperplane
              giving the sizes of points in a zero-dimensional arrangement
              (default: 50).

            - ranges -- Range for the parameters or a list of ranges of
              parameters, one for each hyperplane, for the parametric plots of
              the hyperplanes.  If a single positive number `r` is given for
              ``ranges``, then all parameters run from -r to r.  Otherwise, for
              a line in the plane, the range has the form [a,b] (default:
              [-3,3]), and for a plane in 3-space, the range has the form
              [[a,b],[c,d]] (default: [[-3,3],[-3,3]]). (The ranges are centered
              around self.point().)

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0,0],[0,0,1,5]])
            sage: A.plot(hyperplane_opacities=0.5,hyperplane_labels=True,hyperplane_legend=False)

            sage: B = hyperplane_arrangements.semiorder(4)
            sage: B.plot()
            Displaying the essentialization.
            (, )
            sage: a = hyperplane_arrangements.coordinate(3)
            sage: opts = {'hyperplane_colors':['yellow','green','blue']}
            sage: opts['hyperplane_labels'] = True
            sage: opts['label_offsets'] = [(0,2,2),(2,0,2),(2,2,0)]
            sage: opts['hyperplane_legend'] = False
            sage: opts['hyperplane_opacities'] = 0.7
            sage: a.plot(**opts)
            sage: opts['hyperplane_labels'] = 'short'
            sage: a.plot(**opts)
            sage: pts = HyperplaneArrangement([[3,4],[2,5],[7,1]])
            sage: pts.plot(hyperplane_colors=['yellow','black','blue'])
            sage: a = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,2]])
            sage: a.plot(hyperplane_labels=True,label_colors='blue',label_fontsize=18)
            sage: a.plot(hyperplane_labels=True,label_colors=['red','green','black'])
            sage: pts = HyperplaneArrangement([[3,4],[2,5],[7,1]])
            sage: pts.plot(point_sizes=[50,100,200],hyperplane_colors='blue')
            sage: c = HyperplaneArrangement([[1,0,0,0],[0,0,1,5]])
            sage: c.plot(ranges=10)
            (, )
            sage: c.plot(ranges=[[9.5,10],[-3,3]])
            (, )
            sage: c.plot(ranges=[[[9.5,10],[-3,3]],[[-6,6],[-5,5]]])
            (, )
            sage: s = HyperplaneArrangement([[1,1,0],[1,-1,0],[0,1,2]])
            sage: s.plot(ranges=20)

            sage: s.plot(ranges=[-1,10])

            sage: s.plot(ranges=[[-1,1],[-5,5],[-1,10]])

        NOTES::

            For more examples, see :meth:`show` (``self.show?`` or
            ``HyperplaneArrangement.show?).
        """
        if self.base_field() != QQ:
            raise NotImplementedError('Field must be QQ')
        elif self.dim() == 4:
            if not self.is_essential():
                print 'Displaying the essentialization.'
                self = self.essentialization() # I hope this is OK!
        elif self.dim() not in [1,2,3]: # revise to handle 4d
            return # silently
        # handle extra keywords
        if kwds.has_key('hyperplane_colors'):
            hyp_colors = kwds.pop('hyperplane_colors')
            if not type(hyp_colors) == list: # we assume its a single color then
                hyp_colors = [hyp_colors]*self.num_hyperplanes()
        else:
            N = self.num_hyperplanes()
            HSV_tuples = [(i*1.0/N, 0.8, 0.9) for i in range(N)]
            hyp_colors = map(lambda x: hsv_to_rgb(*x), HSV_tuples)
        if kwds.has_key('hyperplane_labels'):
            hyp_labels = kwds.pop('hyperplane_labels')
            has_hyp_label = True
            if not type(hyp_labels) == list: # we assume its a boolean then
                hyp_labels = [hyp_labels]*self.num_hyperplanes()
            relabeled = []
            for i in range(self.num_hyperplanes()):
                if hyp_labels[i] in [True,'long']:
                    relabeled.append(True)
                else:
                    relabeled.append(str(i))
            hyp_labels = relabeled
        else:
            has_hyp_label = False
        if kwds.has_key('label_colors'):
            label_colors = kwds.pop('label_colors')
            has_label_color = True
            if not type(label_colors) == list: # we assume its a single color then
                label_colors = [label_colors]*self.num_hyperplanes()
        else:
            has_label_color = False
        if kwds.has_key('label_fontsize'):
            label_fontsize = kwds.pop('label_fontsize')
            has_label_fontsize = True
            if not type(label_fontsize) == list: # we assume its a single size then
                label_fontsize = [label_fontsize]*self.num_hyperplanes()
        else:
            has_label_fontsize = False
        if kwds.has_key('label_offsets'):
            has_offsets = True
            offsets = kwds.pop('label_offsets')
        else:
            has_offsets = False # give default values below
        if kwds.has_key('hyperplane_legend'):
            hyperplane_legend = kwds.pop('hyperplane_legend')
        else:
            hyperplane_legend = 'long'
        if kwds.has_key('hyperplane_opacities'):
            hyperplane_opacities = kwds.pop('hyperplane_opacities')
            has_opacity = True
            if not type(hyperplane_opacities) == list: # we assume a single number then
                hyperplane_opacities = [hyperplane_opacities]*self.num_hyperplanes()
        else:
            has_opacity = False
        if kwds.has_key('point_sizes'):
            point_sizes = kwds.pop('point_sizes')
        else:
            point_sizes = 50
        if not type(point_sizes) == list:
            point_sizes = [point_sizes]*self.num_hyperplanes()
        if kwds.has_key('ranges'):
            ranges_set = True
            ranges = kwds.pop('ranges')
            if not type(ranges) in [list,tuple]: # ranges is a single number
                ranges = [ranges]*self.num_hyperplanes()
            # So ranges is some type of list.
            elif self.dim() == 2: # arrangement of lines in the plane
                if not type(ranges[0]) in [list,tuple]: # a single interval
                    ranges = [ranges]*self.num_hyperplanes()
            elif self.dim() == 3: # arrangement of planes in 3-space
                if not type(ranges[0][0]) in [list,tuple]:
                    ranges = [ranges]*self.num_hyperplanes()
            elif self.dim() not in [2,3]: # ranges is not an option unless dim is 2 or 3
                ranges_set = False
            else: # a list of intervals, one for each hyperplane is given
                pass # ranges does not need to be modified
        else:
            ranges_set = False # give default values below
        # the extra keywords have now been handled
        # now handle the legend
        if self.dim() in [1,2]: # points on a line or lines in the plane
            N = self.num_hyperplanes()
            if hyperplane_legend in [True,'long']:
                hyps = self.hyperplanes()
                legend_labels = [hyps[i]._pretty_print_equation(latex=True) for i
                        in range(N)]
            elif hyperplane_legend == 'short' :
                legend_labels = [str(i) for i in range(N)]
        else: # self.dim()==3,  arrangement of planes in 3-space
            if hyperplane_legend in [True, 'long']:
                legend3d = self._3d_legend(hyp_colors, 'long')
            elif hyperplane_legend == 'short':
                legend3d = self._3d_legend(hyp_colors, 'short')
        ## done handling the legend
        ## now create the plot
        p = Graphics()
        for i in range(self.num_hyperplanes()):
            newk = copy(kwds)
            if has_hyp_label:
                newk['hyperplane_label'] = hyp_labels[i]
                if has_offsets:
                    if type(offsets) != list:
                        newk['label_offset'] = offsets
                    else:
                        newk['label_offset'] = offsets[i]
            else:
                newk['hyperplane_label'] = False
            if has_label_color:
                newk['label_color'] = label_colors[i]
            if has_label_fontsize:
                newk['label_fontsize'] = label_fontsize[i]
            if has_opacity:
                newk['opacity'] = hyperplane_opacities[i]
            if self.dim() == 1:
                newk['point_size'] = point_sizes[i]
            if self.dim() in [1,2] and hyperplane_legend != False: # more options than T/F
                newk['legend_label'] = legend_labels[i]
            if ranges_set:
                newk['ranges'] = ranges[i]
            h = self.hyperplanes()[i]
            p += h.plot(rgbcolor=hyp_colors[i], **newk)
        if self.dim() == 1:
            if hyperplane_legend != False: # there are more options than T/F
                p.legend(True)
            return p
        elif self.dim() == 2:
            if hyperplane_legend != False: # there are more options than T/F
                p.legend(True)
            return p
        else: # self.dim()==3
            if hyperplane_legend != False: # there are more options than T/F
                return p, legend3d
            else:
                return p

    def _3d_legend(self, hyperplane_colors, length):
        r"""
        Create plot of a 3d legend for an arrangement of planes in 3-space.  The
        ``length`` parameter determines whether short or long labels are used in
        the legend.

        INPUT:

        - hyperplane_colors -- list of colors
        - length -- 'short' or 'long'

        OUTPUT:

        - Graphics

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a._3d_legend(colors.values()[:6],length='long')

            sage: b = hyperplane_arrangements.semiorder(4)
            sage: c = b.essentialization()
            sage: c._3d_legend(colors.values()[:12],length='long')

            sage: c._3d_legend(colors.values()[:12],length='short')

            sage: p = c._3d_legend(colors.values()[:12],length='short')
            sage: p.set_legend_options(ncol=4)
            sage: p
        """
        if self.dim() != 3:
            raise TypeError('Only implemented for arrangements of planes in 3-space')
        hyps = self.hyperplanes()
        N = self.num_hyperplanes()
        if length == 'short':
            labels = ['  ' + str(i) for i in range(N)]
        else:
            labels = ['  ' + hyps[i]._pretty_print_equation() for i in
                    range(N)]
        p = Graphics()
        for i in range(N):
            p += line([(0,0),(0,0)],color=hyperplane_colors[i], thickness=8,
                    legend_label=labels[i], axes=False)
        p.set_legend_options(title='Hyperplanes', loc='center', labelspacing=0.4, 
                fancybox=True, font_size='x-large', ncol=2)
        p.legend(True)
        return p

    def union(self, other):
        r"""
        The union of ``self`` with ``other``.  The base fields must agree.

        INPUT:

        - ``other`` -- HyperplaneArrangement

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: A = HyperplaneArrangement([[1,0,0],[0,1,1],[0,1,-1],[1,-1,0],[1,1,0]])
            sage: B = HyperplaneArrangement([[1,1,1],[1,-1,1],[1,0,-1]])
            sage: A.union(B)
            Hyperplane arrangement of 8 hyperplanes over Rational Field of dimension 2, rank 2.
        """
        hyps = [i.equation(True) for i in self.hyperplanes()]
        hyps += [i.equation(True) for i in other.hyperplanes()]
        return HyperplaneArrangement(hyps, self.base_field())

        # self._characteristic_polynomial while constructing the special hyperplane
        # arrangement.  Also set poset?

    def varchenko_matrix(self):
        r"""
        Returns the Varchenko matrix of the arrangement. Let `H_1, ..., H_s`
        and `R_1, ..., R_t` denote the hyperplanes and regions, respectively, of
        the arrangement.  Let `S = QQ[h_1,...,h_s]`, a polynomial ring with
        indeterminate `h_i` corresponding to hyperplane `H_i`.  The Varchenko
        matrix is the `t \times t` matrix with `i,j`th entry the product of those
        `h_k` such that `H_k` separates `R_i` and `R_j`.

        It is assumed that the base field for the arrangement is ``QQ``.
        INPUT:

        - None

        OUTPUT:

        - matrix

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(3)
            sage: v = a.varchenko_matrix()
            sage: v
            [    1    h0    h1]
            [   h0     1 h0*h1]
            [   h1 h0*h1     1]
            sage: factor(det(v))
            (h1 - 1) * (h1 + 1) * (h0 - 1) * (h0 + 1)
        """
        hyps = self.hyperplanes()
        k = len(hyps)
        r = self.regions()
        # a polynomial ring with one indeterminate for each hyperplane
        R = PolynomialRing(QQ,'h',k)
        h = R.gens()
        v = zero_matrix(R, k,k)
        for i in range(k):
            v[i,i] = 1
        for i in range(k):
            for j in range(i+1,k):
                t = [h[p] for p in range(k) if
                        self.is_separating_hyperplane(r[i],r[j],hyps[p])]
                t = prod(t)
                v[i,j] = t
                v[j,i] = t
        return v

    def _set_whitney_numbers(self):
        r"""
        Computes doubly-indexed Whitney numbers of the first and second kind.

        See the documentation for :meth:`whitney_numbers`.

        INPUT:

        - None

        OUTPUT:

        - None

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(3)
            sage: a._set_whitney_numbers()
        """
        p = self._intersection_poset
        r = p.rank_function()
        top = r(p.maximal_elements()[0])  # this will always be self.dim(), right?
        m1 = zero_matrix(top+1, top+1)
        m2 = zero_matrix(top+1, top+1)
        for i, j in p.relations_iterator():
            m1[r(i), r(j)] += p.mobius_function(i, j)
            m2[r(i), r(j)] += 1
        self._whitney_numbers = [m1, m2]

    def whitney_number(self, k, kind=1):
        r"""
        The ``k``-th Whitney number.

        If ``kind=1``, this number is obtained by summing the Moebius function
        values `mu(0, x)` over all `x` in the intersection poset with
        `\mathrm{rank}(x) = k`.

        If ``kind=2``, this number is the number of elements `x, y` in the
        intersection poset such that `x \leq y` with ranks `i` and `j`,
        respectively.

        INPUT:

        - ``k`` -- integer
        - ``kind`` -- 1 or 2 (default: 1)

        OUTPUT:

        - integer

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.whitney_number(0)
            1
            sage: A.whitney_number(1)
            -6
            sage: A.whitney_number(2)
            9
            sage: A.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x
            sage: A.whitney_number(1,kind=2)
            6
            sage: p = A.intersection_poset()
            sage: r = p.rank_function()
            sage: len([i for i in p if r(i) == 1])
            6

        REFERENCES::

        .. [GZ] Greene; Zaslavsky
           "On the Interpretation of Whitney Numbers Through Arrangements of
           Hyperplanes, Zonotopes, Non-Radon Partitions, and Orientations of
           Graphs"
           Transactions of the American Mathematical Society, Vol. 280, No. 1.
           (Nov., 1983), pp. 97-126.

        .. SEEALSO::

            :meth:`doubly_indexed_whitney_number`
            :meth:`whitney_data`
        """
        if k >= 0 and k <= self.dim():
            if kind == 1:
                return self._whitney_numbers[0][0, k]
            elif kind == 2:
                return self._whitney_numbers[1][0, k]
        else:
            raise UserWarning('argument out of range')

    def doubly_indexed_whitney_number(self, i, j, kind=1):
        r"""
        The `i,j`-th  doubly-indexed Whitney number.

        If ``kind=1``, this number is obtained by adding the Moebius function
        values `mu(x,y)` over all `x, y` in the intersection poset with
        `\mathrm{rank}(x) = i` and `\mathrm{rank}(y) = j`.

        If `kind=2`, this number is the number of elements `x,y` in the
        intersection poset such that `x \leq y` with ranks `i` and `j`,
        respectively.

        INPUT:

        - ``i``, ``j`` -- integers
        - ``kind`` -- 1 or 2 (default: 1)

        OUTPUT:

        - integer

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.doubly_indexed_whitney_number(0,2)
            9
            sage: A.whitney_number(2)
            9
            sage: A.doubly_indexed_whitney_number(1,2)
            -15

        REFERENCES::

        .. [GZ] Greene; Zaslavsky
           "On the Interpretation of Whitney Numbers Through Arrangements of
           Hyperplanes, Zonotopes, Non-Radon Partitions, and Orientations of
           Graphs"
           Transactions of the American Mathematical Society, Vol. 280, No. 1.
           (Nov., 1983), pp. 97-126.

        .. SEEALSO::

            :meth:`whitney_number`
            :meth:`whitney_data`
        """
        if 0 <= i and j <= self.dim():
            if kind == 1:
                return self._whitney_numbers[0][i, j]
            elif kind == 2:
                return self._whitney_numbers[1][i, j]
        else:
            raise UserWarning('argument out of range')

    def whitney_data(self, kind=1):
        r"""
        The matrix doubly-indexed Whitney numbers of the first or second kind.
        The `i,j`-th entry is the `i,j`-th doubly-indexed Whitney number.

        See the documentation for :meth:`doubly_indexed_whitney_number`.

        INPUT:

        - ``kind`` -- 1 or 2 (default: 1)

        OUTPUT:

        - matrix

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.whitney_data()
            [  1  -6   9]
            [  0   6 -15]
            [  0   0   6]
            sage: A.whitney_data(2)
            [ 1  6  6]
            [ 0  6 15]
            [ 0  0  6]

        .. SEEALSO::

            :meth:`whitney_number`
            :meth:`doubly_indexed_whitney_number`
        """
        if kind == 1:
            return deepcopy(self._whitney_numbers[0])
        elif kind == 2:
            return deepcopy(self._whitney_numbers[1])
        else:
            raise UserWarning('argument out of range')


class HyperplaneArrangementGenerators():

    def braid(self, n, K=QQ):
        r"""
        The braid arrangement

        This is the set of `n(n-1)/2` hyperplanes: `\{ x_i - x_j = 0 :
        1\leq i \leq j\leq n\\}.`

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default: ``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.braid(4)
            Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 4, rank 3.
        """
        x = polygen(QQ, 'x')
        A = self.graphical(graphs.CompleteGraph(n), K)
        A._characteristic_polynomial = prod(x-i for i in range(n))
        return A

    def bigraphical(self, G, A=None, K=QQ):
        r"""
        The hyperplane arrangement with hyperplanes `x_i - x_j = A[i,j]` and
        `x_j - x_i = A[j,i]` for each edge `v_i, v_j` of ``G``.  The indices
        `i,j` are the indices of elements of ``G.vertices()``.

        INPUT:

        - ``G`` -- Graph
        - ``A`` -- list, matrix, or dictionary (default: None gives semiorder), 'generic'
        - ``K`` -- field (default: `QQ`)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CycleGraph(4)
            sage: G.edges()
            [(0, 1, None), (0, 3, None), (1, 2, None), (2, 3, None)]
            sage: G.edges(labels=False)
            [(0, 1), (0, 3), (1, 2), (2, 3)]
            sage: A = {0:{1:1, 3:2}, 1:{0:3, 2:0}, 2:{1:2, 3:1}, 3:{2:0, 0:2}}
            sage: HA = hyperplane_arrangements.bigraphical(G,A)
            sage: HA.num_regions()
            63
            sage: hyperplane_arrangements.bigraphical(G,'generic').num_regions()
            65
            sage: hyperplane_arrangements.bigraphical(G).num_regions()
            59

        REFERENCES::

        .. [HP] S. Hopkins, D. Perkinson
           "Bigraphical Arrangements"
           :arxiv:`1212.4398`
        """
        n = G.num_verts()
        if A is None:  # default to G-semiorder arrangement
            A = matrix(K, n, lambda i, j: 1)
        elif A == 'generic':
            A = random_matrix(ZZ, n, x=10000)
            A = matrix(K, A)
        hyperplanes = []
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            new[-1] = A[i][j]
            hyperplanes.append(new)
            new = [0]*(n+1)
            new[j] = 1
            new[i] = -1
            new[-1] = A[j][i]
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def Catalan(self, n, K=QQ):
        r"""
        The Catalan arrangement

        This is the set of `3n(n-1)/2` hyperplanes: `\{ x_i - x_j =
        -1,0,1 : 1\leq i \leq j\leq n\\}.`

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default: ``QQ`)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.Catalan(5)
            Hyperplane arrangement of 30 hyperplanes over Rational Field of dimension 5, rank 4.
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*n
                new[i] = 1
                new[j] = -1
                for k in [-1, 0, 1]:
                    h = deepcopy(new)
                    h.append(k)
                    hyperplanes.append(h)
        Cn = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        Cn._characteristic_polynomial = x*prod([x-n-i for i in range(1, n)])
        return Cn

    def coordinate(self, n, K=QQ):
        r"""
        The coordinate hyperplane arrangement is the central hyperplane
        arrangement consisting of the coordinate hyperplanes `x_i=0`.

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.coordinate(5)
            Hyperplane arrangement of 5 hyperplanes over Rational Field of dimension 5, rank 5.
        """
        hyperplanes = []
        for i in range(n):
            new = [0]*(n+1)
            new[i] = 1
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def G_semiorder(self, G, K=QQ):
        r"""
        The semiorder hyperplane arrangement of a graph G is the arrangement `\{
        x_i - x_j = -1,1\}` where `ij` is an edge of ``G``.

        INPUT:

        - ``G`` -- graph
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.G_semiorder(G)
            Hyperplane arrangement of 20 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.G_semiorder(g)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 5, rank 4.
        """
        hyperplanes = []
        n = G.num_verts()
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            new[-1] = -1
            hyperplanes.append(new)
            new = deepcopy(new)
            new[-1]=1
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def G_Shi(self, G, K=QQ):
        r"""
        Return the Shi hyperplane arrangement of a graph ``G``.

        INPUT:

        - ``G`` -- graph

        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.G_Shi(G)
            Hyperplane arrangement of 20 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.G_Shi(g)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: a = hyperplane_arrangements.G_Shi(graphs.WheelGraph(4))
            sage: a.show(frame=false,hyperplane_legend=false,hyperplane_opacities=0.8)
            Displaying the essentialization.
        """
        hyperplanes = []
        n = G.num_verts()
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            hyperplanes.append(new)
            new = deepcopy(new)
            new[-1]=1
            hyperplanes.append(new)
        return HyperplaneArrangement(hyperplanes, K)

    def graphical(self, G, K=QQ):
        r"""
        The graphical hyperplane arrangement of a graph G is the arrangement `\{
        x_i - x_j = 0\}` where `ij` is an edge of ``G``.

        INPUT:

        - ``G`` -- graph
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: G = graphs.CompleteGraph(5)
            sage: hyperplane_arrangements.graphical(G)
            Hyperplane arrangement of 10 hyperplanes over Rational Field of dimension 5, rank 4.
            sage: g = graphs.HouseGraph()
            sage: hyperplane_arrangements.graphical(g)
            Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 5, rank 4.
        """
        hyperplanes = []
        n = G.num_verts()
        for e in G.edges():
            i = G.vertices().index(e[0])
            j = G.vertices().index(e[1])
            new = [0]*(n+1)
            new[i] = 1
            new[j] = -1
            hyperplanes.append(new)
        A = HyperplaneArrangement(hyperplanes, K)
        A._characteristic_polynomial = G.chromatic_polynomial()
        return A

    def Ish(self, n, K=QQ):
        r"""
        The Ish arrangement is the set of `n(n-1)` hyperplanes: ``\{ x_i - x_j =
        0 : 1\leq i \leq j\leq n\} \cup \{x_1 - x_j = i : 1\leq i \leq j\leq n\}:.``

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: a = hyperplane_arrangements.Ish(3)
            sage: a
            Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 3, rank 2.
            sage: a.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x
            sage: b = hyperplane_arrangements.Shi(3)
            sage: b.characteristic_polynomial()
            x^3 - 6*x^2 + 9*x

        REFERENCES::

        .. [AR] D. Armstrong, B. Rhoades
           "The Shi arrangement and the Ish arrangement"
           :arxiv:`1009.1655`
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*(n+1)
                new[i] = 1
                new[j] = -1
                # x_i - x_j = 0
                hyperplanes.append(new)
                # x_1 - x_j = i
                new = [0]*(n+1)
                new[0] = 1
                new[j] = -1
                new[-1] = i + 1
                hyperplanes.append(new)
                A = HyperplaneArrangement(hyperplanes, K)
                x = polygen(QQ, 'x')
                cp = sum([(-1)**k*stirling_number2(n,n-k)*prod([(x-1-j) for j in range(k,n-1)]) for k in range(0,n)])
                cp = x*cp
                cp = expand(cp)
                A._characteristic_polynomial = cp
        return A

    def linial(self, n, K=QQ):
        r"""
        The linial hyperplane arrangement is the set of hyperplanes
        ``\{x_i - x_j = 1 : 1\leq i < j \leq n\}`` 

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: a = hyperplane_arrangements.linial(4)
            sage: a.characteristic_polynomial()
            x^4 - 6*x^3 + 15*x^2 - 14*x
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1,n):
                new = [0]*(n+1)
                new[i] = 1
                new[j] = -1
                new[-1] = 1
                hyperplanes.append(new)
        A = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        cp = expand(x*sum(binomial(n,k)*(x-k)**(n-1) for k in range(n+1))/2**n)
        A._characteristic_polynomial = cp
        return A

    def semiorder(self, n, K=QQ):
        r"""
        The semiorder arrangement is the set of `n(n-1)` hyperplanes: `\{ x_i -
        x_j = -1,1 : 1\leq i \leq j\leq n\}.`

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.semiorder(4)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 4, rank 3.
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*n
                new[i] = 1
                new[j] = -1
                for k in [-1, 1]:
                    h = deepcopy(new)
                    h.append(k)
                    hyperplanes.append(h)
        A = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        cp = x*sum([stirling_number2(n,k)*prod([x-k-i for i in range(1, k)]) for k in range(1,n+1)])
        cp = expand(cp)
        A._characteristic_polynomial = cp
        return A

    def Shi(self, n, K=QQ):
        r"""
        The Shi arrangement is the set of `n(n-1)` hyperplanes: ``\{ x_i - x_j =
        0,1 : 1\leq i \leq j\leq n\}.``

        INPUT:

        - ``n`` -- integer
        - ``K`` -- field (default:``QQ``)

        OUTPUT:

        - HyperplaneArrangement

        EXAMPLES::

            sage: hyperplane_arrangements.Shi(4)
            Hyperplane arrangement of 12 hyperplanes over Rational Field of dimension 4, rank 3.
        """
        hyperplanes = []
        for i in range(n):
            for j in range(i+1, n):
                new = [0]*n
                new[i] = 1
                new[j] = -1
                for k in [0, 1]:
                    h = deepcopy(new)
                    h.append(k)
                    hyperplanes.append(h)
        A = HyperplaneArrangement(hyperplanes, K)
        x = polygen(QQ, 'x')
        cp = sum([(-1)**k*stirling_number2(n,n-k)*prod([(x-1-j) for j in range(k,n-1)]) for k in range(0,n)])
        cp = x*cp
        cp = expand(cp)
        A._characteristic_polynomial = cp
        return A

hyperplane_arrangements = HyperplaneArrangementGenerators()


r"""
Hyperplane Arrangements

Before talking about hyperplane arrangements, let us start with
individual hyperplanes. This package uses certain linear expressions
to represent hyperplanes, that is, a linear expression `3x+3y-5z-7`
stands for the hyperplane with the equation `x+3y-5z=7`. To create it
in Sage, you first have to create a :class:`HyperplaneArrangenments`
object to define the variables `x`, `y`, `z`::

    sage: H.<x,y,z> = HyperplaneArrangements(QQ)
    sage: h = 3*x + 2*y - 5*z - 7;  h
    Hyperplane 3*x + 2*y - 5*z - 7
    sage: h.normal()
    (3, 2, -5)
    sage: h.constant_term()
    -7

The individual hyperplanes behave like the linear expression with
regard to addition and scalar multiplication, which is why you can do
linear combinations of the coordinates::

    sage: -2*h
    Hyperplane -6*x - 4*y + 10*z + 14
    sage: x, y, z
    (Hyperplane x + 0*y + 0*z + 0, Hyperplane 0*x + y + 0*z + 0, Hyperplane 0*x + 0*y + z + 0)

ARRANGEMENTS

There are several ways to create hyperplane arrangements:

(i) by passing individual hyperplanes to the
:class:`HyperplaneArrangenments` object::

    sage: H.<x,y> = HyperplaneArrangements(QQ)
    sage: box = H(x, y, x-1, y-1);  box
    Arrangement <y - 1 | y | x - 1 | x>

(ii) by passing anything that defines a hyperplane, for example a
coefficient vector and constant term::

    sage: H = HyperplaneArrangements(QQ, ('x', 'y'))
    sage: triangle = H([(1, 0), 0], [(0, 1), 0], [(1,1), -1]);  triangle
    Arrangement <y | x | x + y - 1>

    sage: H.inject_variables()
    Defining x, y
    sage: triangle == H(x, y, x+y-1)
    True

The default base field is QQ, the rational numbers.  Finite fields are also
supported::

    sage: H.<x,y,z> = HyperplaneArrangements(GF(5))
    sage: a = H([(1,2,3), 4], [(5,6,7), 8]);  a
    Arrangement <y + 2*z + 3 | x + 2*y + 3*z + 4>

(iii) a list or tuple of hyperplanes::

    sage: H.<x,y,z> = HyperplaneArrangements(GF(5))
    sage: k = [x+i for i in range(4)];  k
    [Hyperplane x + 0*y + 0*z + 0, Hyperplane x + 0*y + 0*z + 1, 
     Hyperplane x + 0*y + 0*z + 2, Hyperplane x + 0*y + 0*z + 3]
    sage: H(k)
    Arrangement <x | x + 1 | x + 2 | x + 3>

(iv) using the library of arrangements::

    sage: hyperplane_arrangements.braid(4)
    Arrangement of 6 hyperplanes of dimension 4 and rank 3
    sage: hyperplane_arrangements.semiorder(3)
    Arrangement of 6 hyperplanes of dimension 3 and rank 2
    sage: hyperplane_arrangements.graphical(graphs.PetersenGraph())
    Arrangement of 15 hyperplanes of dimension 10 and rank 9
    sage: hyperplane_arrangements.Ish(5)
    Arrangement of 20 hyperplanes of dimension 5 and rank 4

(v) from the bounding hyperplanes of a polyhedron::

    
    sage: a = polytopes.n_cube(3).hyperplane_arrangement();  a
    Arrangement of 6 hyperplanes of dimension 3 and rank 3
    sage: a.n_regions()
    27

New arrangements from old::

    sage: a = hyperplane_arrangements.braid(3)
    sage: b = a.add_hyperplane([4, 1, 2, 3])     
    sage: b
    Arrangement <t1 - t2 | t0 - t1 | t0 - t2 | t0 + 2*t1 + 3*t2 + 4>
    sage: c = b.deletion([4, 1, 2, 3])
    sage: a == c
    True

    sage: a = hyperplane_arrangements.braid(3)
    sage: b = a.union(hyperplane_arrangements.semiorder(3))
    sage: b == hyperplane_arrangements.Catalan(3)
    True

    sage: a
    Arrangement <t1 - t2 | t0 - t1 | t0 - t2>
    sage: a = hyperplane_arrangements.coordinate(4)
    sage: h = a.hyperplanes()[0]
    sage: b = a.restriction(h)
    sage: b == hyperplane_arrangements.coordinate(3)
    True

A hyperplane arrangement is *essential* is the normals to its
hyperplane span the ambient space.  Otherwise, it is *inessential*.
The essentialization is formed by intersecting the hyperplanes by this
normal space (actually, it is a bit more complicated over finite
fields)::


#     sage: a = hyperplane_arrangements.braid(4)
#     sage: a
#     Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 4, rank 3.
#     sage: a.is_essential()
#     False
#     sage: a.rank() < a.dim()  # double-check
#     True
#     sage: a.essentialization()
#     Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 3, rank 3.

# The connected components of the complement of the hyperplanes of an arrangement
# in `\mathbb{R}^n` are called the *regions* of the arrangement::

#     sage: a = hyperplane_arrangements.semiorder(3)
#     sage: b = a.essentialization()
#     sage: b
#     Hyperplane arrangement of 6 hyperplanes over Rational Field of dimension 2, rank 2.
#     sage: b.num_regions()
#     19
#     sage: b.regions()
#     [A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices and 1 ray,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays]
#     sage: b.bounded_regions()
#     [A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
#      A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]
#     sage: b.num_bounded_regions()
#     7
#     sage: a.unbounded_regions()
#     [A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices, 1 ray, 1 line,
#      A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 2 rays, 1 line]
#     sage: b.region_containing_point((0,0)).show()
#     sage: b.region_containing_point((2,1)).show(xmax=4,ymax=4)
#     sage: r1 = b.regions()[0]
#     sage: r2 = b.regions()[-1]  
#     sage: b.distance_between_regions(r1,r2)  # number of hyps. separating r1, r2
#     6
#     sage: b.distance_enumerator(r1)  # generating fnc. for distances from r1
#     x^6 + 2*x^5 + 5*x^4 + 3*x^3 + 5*x^2 + 2*x + 1

# Note: *bounded region* really mean *relatively bounded* here.  A region is
# relatively bounded if its intersection with space spanned by the normals of the
# hyperplanes in the arrangement is bounded.

The intersection poset of a hyperplane arrangement is the collection
of all nonempty intersections of hyperplanes in the arrangement,
ordered by reverse inclusion.  It includes the ambient space of the
arrangement (as the intersection over the empty set)::

    sage: a = hyperplane_arrangements.braid(3)
    sage: p = a.intersection_poset()
    sage: p.is_ranked()
    True
    sage: p.order_polytope()
    A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 10 vertices

The characteristic polynomial is a basic invariant of a hyperplane arrangement.
It is defined as `\chi(x) := \sum_{w\in P \mu(w) x^{\dim(w)}` where the sum is
`P` is the intersection poset of the arrangement and `\mu` is the Moebius
function of `P`::

     sage: a = hyperplane_arrangements.semiorder(5)
     sage: a.characteristic_polynomial()               # long time (about a second on Core i7)
     x^5 - 20*x^4 + 180*x^3 - 790*x^2 + 1380*x
     sage: a.poincare_polynomial()                     # long time
     1380*x^4 + 790*x^3 + 180*x^2 + 20*x + 1
     sage: a.n_regions()                               # long time
     2371
     sage: charpoly = a.characteristic_polynomial()    # long time
     sage: charpoly(-1)                                # long time
     -2371
     sage: a.n_bounded_regions()                       # long time
     751
     sage: charpoly(1)                                 # long time
     751

# For finer invariants derived from the intersection poset, see
# :meth:`whitney_number` and :meth:`doubly_indexed_whitney_number`.

# Miscellaneous methods (see documentation for an explanation)::

#     sage: a = hyperplane_arrangements.semiorder(3)
#     sage: a.has_good_reduction(5)
#     True
#     sage: b = a.change_base_field(GF(5))
#     sage: pa = a.intersection_poset()
#     sage: pb = b.intersection_poset()
#     sage: pa.is_isomorphic(pb)
#     True
#     sage: a.face_vector()
#     (0, 12, 30, 19)
#     sage: a.face_vector()
#     (0, 12, 30, 19)
#     sage: a.is_central()
#     False
#     sage: a.is_linear()
#     False
#     sage: a.sign_vector((1,1,1))
#     [1, -1, 1, -1, 1, -1]
#     sage: a.varchenko_matrix()
#     [          1          h0       h0*h1       h0*h2    h0*h1*h2 h0*h1*h2*h3]
#     [         h0           1          h1          h2       h1*h2    h1*h2*h3]
#     [      h0*h1          h1           1       h1*h2          h2       h2*h3]
#     [      h0*h2          h2       h1*h2           1          h1       h1*h3]
#     [   h0*h1*h2       h1*h2          h2          h1           1          h3]
#     [h0*h1*h2*h3    h1*h2*h3       h2*h3       h1*h3          h3           1]

# There are extensive methods for visualizing hyperplane arrangements in low
# dimensions.  See :meth:`plot` and :meth:`show` for details.


# AUTHORS:

# - David Perkinson (2013-06): initial version

# - Qiaoyu Yang (2013-07)

# - Kuai Yu (2013-07)

# This module implements hyperplane arrangements defined over the
# rationals or over finite fields.  The original motivation was to make a
# companion to Richard Stanley's notes on hyperplane arrangements.

# REFERENCES::

# .. [RS] Stanley, Richard
#    "Hyperplane Arrangements"
#     Geometric Combinatorics (E. Miller, V. Reiner, and B. Sturmfels, eds.),
#     IAS/Park City Mathematics Series, vol. 13, American Mathematical Society,
#     Providence, RI, 2007, pp. 389-496.
"""

#*****************************************************************************
#       Copyright (C) 2013 David Perkinson <davidp@reed.edu>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Possible extensions for hyperplane_arrangement.py:
# - the big face lattice
# - Orlik-Solomon algebras
# - create ties with the Sage matroid methods
# - hyperplane arrangements over other fields

from sage.structure.parent import Parent
from sage.structure.element import Element, coerce_binop
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import QQ, ZZ
from sage.misc.cachefunc import cached_method
from sage.misc.misc import uniq
from sage.matrix.constructor import matrix, vector

from sage.geometry.hyperplane_arrangement.hyperplane import AmbientVectorSpace, Hyperplane



class HyperplaneArrangementElement(Element):
    
    def __init__(self, parent, hyperplanes, check=True):
        super(HyperplaneArrangementElement, self).__init__(parent)
        self._hyperplanes = hyperplanes
        if check:
            assert isinstance(hyperplanes, tuple)
            assert all(isinstance(h, Hyperplane) for h in hyperplanes)
            assert all(h.parent() is self.parent().ambient_space() for h in hyperplanes)

    def _first_ngens(self, n):
        """
        Workaround to support the construction with names.

        INPUT/OUTPUT:

        See :meth:`HyperplaneArrangements._first_ngens`

        EXAMPLES::

            sage: a.<x,y,z> = hyperplane_arrangements.braid(3)
        """
        return self.parent()._first_ngens(n)

    def __getitem__(self, i):
        return self._hyperplanes[i]
        
    def n_hyperplanes(self):
        r"""
        Return the number of hyperplanes in the arrangement.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,1,0], [2,3,-1], [4,5,3])
            sage: A.n_hyperplanes()
            3
            sage: len(A)    # equivalent
            3
        """
        return len(self._hyperplanes)

    __len__ = n_hyperplanes

    def hyperplanes(self):
        r"""
        Return the number of hyperplanes in the arrangement.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,1,0], [2,3,-1], [4,5,3])
            sage: A.hyperplanes()
            (Hyperplane x + 0*y + 1, Hyperplane 3*x - y + 2, Hyperplane 5*x + 3*y + 4)

        Note that the hyperplanes can be indexed as if they were a list::

            sage: A[0]
            Hyperplane x + 0*y + 1
        """
        return self._hyperplanes

    def _repr_(self):
        r"""
        String representation for a hyperplane arrangement.

        OUTPUT:

        String.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: H(x, y, x-1, y-1)
            Arrangement <y - 1 | y | x - 1 | x>
            sage: H(x, y, x - 1, y - 1, x + y, x - y)
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
        """
        if len(self) < 5:
            hyperplanes = ' | '.join(h._repr_linear(include_zero=False) for h in self._hyperplanes)
            return 'Arrangement <{0}>'.format(hyperplanes)
        return 'Arrangement of {0} hyperplanes of dimension {1} and rank {2}'.format(
            len(self), self.dimension(), self.rank())

    def dimension(self):
        return self.parent().ngens()

    def rank(self):
        """
        Return the rank.

        OUTPUT:

        Integer. The dimension of the span of the normals to the
        hyperplanes in the arrangement.


        EXAMPLES::

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H([[0, 1, 2, 3],[-3, 4, 5, 6]])
            sage: A.dimension()
            3
            sage: A.rank()
            2

            sage: B = hyperplane_arrangements.braid(3)
            sage: B.hyperplanes()
            (Hyperplane 0*t0 + t1 - t2 + 0, 
             Hyperplane t0 - t1 + 0*t2 + 0, 
             Hyperplane t0 + 0*t1 - t2 + 0)
            sage: B.dimension()
            3
            sage: B.rank()
            2

            sage: p = polytopes.n_simplex(5)
            sage: H = p.hyperplane_arrangement()
            sage: H.rank()
            5
        """
        R = self.parent().base_ring()
        normals = [h.normal() for h in self]
        return matrix(R, normals).rank()

    def __cmp__(self, other):
        assert (type(self) is type(other)) and (self.parent() is other.parent()) # guaranteed by framework
        return cmp(self._hyperplanes, other._hyperplanes)

    def union(self, other):
        r"""
        The union of ``self`` with ``other``.

        INPUT:

        - ``other`` -- A hyperplane arrangement or something that can
          be converted into a hyperplane arrangement.

        OUTPUT:

        A new hyperplane arrangement.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([1,2,3], [0,1,1], [0,1,-1], [1,-1,0], [1,1,0])
            sage: B = H([1,1,1], [1,-1,1], [1,0,-1])
            sage: A.union(B)
            Arrangement of 8 hyperplanes of dimension 2 and rank 2

        A single hyperplane is coerced into a hyperplane arrangement if necessary::

            sage: A.union(x+y-1)
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
            sage: A.add_hyperplane(x+y-1)    # alias
            Arrangement of 6 hyperplanes of dimension 2 and rank 2

            sage: P.<x,y> = HyperplaneArrangements(RR)
            sage: C = P(2*x + 4*y + 5)
            sage: C.union(A)
            Arrangement of 6 hyperplanes of dimension 2 and rank 2
        """
        P = self.parent()
        other = P(other)
        hyperplanes = self._hyperplanes + other._hyperplanes
        return P(*hyperplanes)

    add_hyperplane = union

    def plot(self, **kwds):
        """
        Plot the hyperplane arrangement.

        OUTPUT:

        Graphics object.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: L(x, y, x+y-2).plot()
        """
        from sage.geometry.hyperplane_arrangement.plot import plot
        return plot(self, **kwds)

    def cone(self, variable='t'):
        r"""
        Return the cone over the hyperplane arrangement.  

        INPUT:

        - ``variable`` -- string. The name of the additional variable.

        OUTPUT:

        A new yperplane arrangement. Its equations consist of
        `[0,-d,a_1,...,a_n]` for each `[d,a_1,...,a_n]` in the original
        arrangement and the equation `[0,1,0,...,0]`.

        EXAMPLES::

            sage: a.<x,y,z> = hyperplane_arrangements.semiorder(3)
            sage: b = a.cone()
            sage: a.characteristic_polynomial().factor()
            x * (x^2 - 6*x + 12)
            sage: b.characteristic_polynomial().factor()
            (x - 1) * x * (x^2 - 6*x + 12)
            sage: a.hyperplanes()
            (Hyperplane 0*x + y - z - 1, 
             Hyperplane 0*x + y - z + 1, 
             Hyperplane x - y + 0*z - 1, 
             Hyperplane x - y + 0*z + 1, 
             Hyperplane x + 0*y - z - 1, 
             Hyperplane x + 0*y - z + 1)
            sage: b.hyperplanes()
            (Hyperplane -t + 0*x + y - z + 0, 
             Hyperplane -t + x - y + 0*z + 0, 
             Hyperplane -t + x + 0*y - z + 0, 
             Hyperplane t + 0*x + 0*y + 0*z + 0, 
             Hyperplane t + 0*x + y - z + 0, 
             Hyperplane t + x - y + 0*z + 0, 
             Hyperplane t + x + 0*y - z + 0)
        """
        AA = self.parent().ambient_space()
        hyperplanes = []
        for h in self.hyperplanes():
            new_h = [0] + [h.b()] + list(h.A())
            hyperplanes.append(new_h)
        hyperplanes.append([0, 1] + [0] * self.dimension())
        P = self.parent()
        names = (variable,) + P._names
        H = HyperplaneArrangements(self.parent().base_ring(), names)
        return H(*hyperplanes)
        
    @cached_method
    def intersection_poset(self):
        r"""
        Return the intersection poset of the hyperplane arrangement.

        OUTPUT:

        The poset of non-empty intersections of hyperplanes.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a.intersection_poset()
            Finite poset containing 4 elements

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.intersection_poset()
            Finite poset containing 19 elements
        """
        K = self.base_ring()
        from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
        from sage.modules.all import VectorSpace
        whole_space = AffineSubspace(0, VectorSpace(K, self.dimension()))
        L = [[whole_space]]
        active = True
        codim = 0
        while active:
            active = False
            new_level = []
            for T in L[codim]:
                for H in self:
                    I = H._affine_subspace().intersection(T)
                    if I is not None and I != T and I not in new_level:
                        new_level.append(I)
                        active = True
            if active:
                L.append(new_level)
            codim += 1
        from sage.misc.flatten import flatten
        L = flatten(L)
        t = {}
        for i in range(len(L)):
            t[i] = L[i]
        cmp_fn = lambda p, q: t[q] < t[p]
        from sage.combinat.posets.posets import Poset
        return Poset((t, cmp_fn))

    def _slow_characteristic_polynomial(self):
        """
        Return the characteristic polynomial of the hyperplane arrangement.

        This is the slow computation directly from the definition. For
        educational use only.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a._slow_characteristic_polynomial()
            x^2 - 2*x + 1
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        x = polygen(QQ, 'x')
        P = self.intersection_poset()
        n = self.dimension()
        return sum([P.mobius_function(0, p) * x**(n - P.rank(p)) for p in P])
        
    @cached_method
    def characteristic_polynomial(self):
        r"""
        Return the characteristic polynomial of the hyperplane arrangement.

        OUTPUT:

        The characteristic polynomial in `\QQ[x]`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a.characteristic_polynomial()
            x^2 - 2*x + 1
        
        TESTS::
          
            sage: H.<s,t,u,v> = HyperplaneArrangements(QQ)
            sage: m = matrix([(0, -1, 0, 1, -1), (0, -1, 1, -1, 0), (0, -1, 1, 0, -1),
            ....:   (0, 1, 0, 0, 0), (0, 1, 0, 1, -1), (0, 1, 1, -1, 0), (0, 1, 1, 0, -1)])
            sage: R.<x> = QQ[]
            sage: expected_charpoly = (x - 1) * x * (x^2 - 6*x + 12)
            sage: for s in SymmetricGroup(4):   # long time (about a second on a Core i7)
            ....:     m_perm = [m.column(i) for i in [0, s(1), s(2), s(3), s(4)]]
            ....:     m_perm = matrix(m_perm).transpose()
            ....:     charpoly = H(m_perm.rows()).characteristic_polynomial()
            ....:     assert charpoly == expected_charpoly
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        x = polygen(QQ, 'x')
        if self.rank() == 1:
            return x**(self.dimension() - 1) * (x - len(self))
        else:
            H = self[0]
            R = self.restriction(H)
            charpoly_R = R.characteristic_polynomial()
            D = self.deletion(H)
            charpoly_D = D.characteristic_polynomial()
            return charpoly_D - charpoly_R

    @cached_method
    def poincare_polynomial(self):
        r"""
        Return the Poincare polynomial of the hyperplane arrangement.

        OUTPUT:

        The Poincare polynomial in `\QQ[x]`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.coordinate(2)
            sage: a.poincare_polynomial()
            x^2 + 2*x + 1
        """
        charpoly = self.characteristic_polynomial()
        R = charpoly.parent()
        x = R.gen(0)
        poincare = (-x)**self.dimension() * charpoly(-QQ(1)/x)
        return R(poincare)

    def deletion(self, hyperplanes):
        r"""
        Return the hyperplane arrangement obtained by removing ``h``.

        INPUT:

        - ``h`` -- a hyperplane or hyperplane arrangement.

        OUTPUT:

        A new hyperplane arrangement with the given hyperplane(s)
        ``h`` removed.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([0,1,0], [1,0,1], [-1,0,1], [0,1,-1], [0,1,1]);  A
            Arrangement of 5 hyperplanes of dimension 2 and rank 2
            sage: A.deletion(x)
            Arrangement <y - 1 | y + 1 | x - y | x + y>
            sage: h = H([0,1,0], [0,1,1])
            sage: A.deletion(h)
            Arrangement <y - 1 | y + 1 | x - y>

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([0,1,0], [1,0,1], [-1,0,1], [0,1,-1], [0,1,1])
            sage: h = H([0,4,0])
            sage: A.deletion(h)
            Arrangement <y - 1 | y + 1 | x - y | x + y>
            sage: l = H([1,2,3])
            sage: A.deletion(l)
            Traceback (most recent call last):
            ...
            ValueError: hyperplane is not in the arrangement

        .. SEEALSO::

            :meth:`restriction`
        """
        parent = self.parent()
        hyperplanes = parent(hyperplanes)
        planes = list(self)
        for hyperplane in hyperplanes:
            try:
                planes.remove(hyperplane)
            except ValueError:
                raise ValueError('hyperplane is not in the arrangement')
        return parent(*planes)

    def restriction(self, hyperplane):
        r"""
        Return the restriction to a hyperplane.

        INPUT:

        - ``hyperplane`` -- a hyperplane of the hyperplane arrangement.

        OUTPUT:

        The restriction of the hyperplane arrangement to the given
        ``hyperplane``.

        EXAMPLES::

            sage: A.<u,x,y,z> = hyperplane_arrangements.braid(4);  A
            Arrangement of 6 hyperplanes of dimension 4 and rank 3
            sage: H = A[0];  H
            Hyperplane 0*u + 0*x + y - z + 0
            sage: R = A.restriction(H);  R
            Arrangement <x - z | u - x | u - z>
            sage: D = A.deletion(H);  D
            Arrangement of 5 hyperplanes of dimension 4 and rank 3
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
        parent = self.parent()
        hyperplane = parent(hyperplane)[0]
        if hyperplane not in self.hyperplanes():
            raise ValueError('hyperplane not in arrangement')
        pivot = hyperplane._normal_pivot()
        n = hyperplane / hyperplane.A()[pivot]
        hyperplanes = []
        for h in self:
            rescale = h.A()[pivot] / hyperplane.A()[pivot]
            h = h - rescale * hyperplane
            A = list(h.A())
            A_pivot = A.pop(pivot)
            assert A_pivot == 0
            if all(a == 0 for a in A):
                continue
            b = h.b()
            hyperplanes.append([A, b])
        names = list(parent._names)
        names.pop(pivot)
        H = HyperplaneArrangements(parent.base_ring(), tuple(names))
        return H(*hyperplanes, signed=False)

    def change_ring(self, base_ring):
        """
        Return hyperplane arrangement over the new base ring.
        
        INPUT:

        - ``base_ring`` -- the new base ring. Must be a field for
          hyperplane arrangements.

        OUTPUT:

        The hyperplane arrangement obtained by changing the base
        field, as a new hyperplane arrangement.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,1), 0], [(2,3), -1])
            sage: A.change_ring(FiniteField(2))
            Arrangement <y + 1 | x + y>
        """
        parent = self.parent().change_ring(base_ring)
        return parent(self)

    @cached_method
    def n_regions(self):
        r"""
        The number of regions of the hyperplane arrangement.

        OUTPUT:

        Integer

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.n_regions()
            19

        TESTS::
         
            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,1), 0], [(2,3), -1], [(4,5), 3])
            sage: B = A.change_ring(FiniteField(7))
            sage: B.n_regions()
            Traceback (most recent call last):
            ...
            TypeError: base field must have characteristic zero
        """
        if self.base_ring().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        charpoly = self.characteristic_polynomial()
        return (-1)**self.dimension() * charpoly(-1)

    @cached_method
    def n_bounded_regions(self):
        r"""
        Return the number of (relatively) bounded regions.

        OUTPUT:

        Integer. The number of relatively bounded regions of the
        hyperplane arrangement.

        EXAMPLES::

            sage: A = hyperplane_arrangements.semiorder(3)
            sage: A.n_bounded_regions()
            7

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,1),0], [(2,3),-1], [(4,5),3])
            sage: B = A.change_ring(FiniteField(7))
            sage: B.n_bounded_regions()
            Traceback (most recent call last):
            ...
            TypeError: base field must have characteristic zero
        """
        if self.base_ring().characteristic() != 0:
            raise TypeError('base field must have characteristic zero')
        charpoly = self.characteristic_polynomial()
        return (-1)**self.rank() * charpoly(1)

    def has_good_reduction(self, p):
        r"""
        Return whether the hyperplane arrangement has good reduction mod `p`.

        INPUT:

        - ``p`` -- prime number

        OUTPUT:

        Boolean. 

        Let `A` be a hyperplane arrangement with equations defined
        over the integers, and let `B` be the hyperplane arrangement
        defined by reducing these equations modulo a prime `p`.  Then
        `A` has good reduction modulo `p` if the intersection posets
        of `A` and `B` are isomorphic.

        EXAMPLES::

            sage: a = hyperplane_arrangements.semiorder(3)
            sage: a.has_good_reduction(5)
            True
            sage: a.has_good_reduction(3)
            False
            sage: b = a.change_ring(GF(3))
            sage: a.characteristic_polynomial()
            x^3 - 6*x^2 + 12*x
            sage: b.characteristic_polynomial()  # not equal to that for a
            x^3 - 6*x^2 + 10*x
        """
        if self.base_ring()!=QQ:
            raise TypeError('Arrangement must be defined over QQ')
        if not p.is_prime():
            raise TypeError('Must reduce modulo a prime number.')
        from sage.rings.all import GF
        a = self.change_ring(GF(p))
        p = self.intersection_poset()
        q = a.intersection_poset()
        return p.is_isomorphic(q)

    def is_essential(self):
        r"""
        Test whether the hyperplane arrangement is essential.

        .. SEEALSO::

            :meth:`essentialization`

        OUTPUT:

        Boolean. Whethe the hyperplane arrangement is essential. A
        hyperplane arrangement is essential if the span of the normals
        of its hyperplanes spans the ambient space.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: H(x, x+1).is_essential()
            False
            sage: H(x, y).is_essential()
            True
        """
        return self.rank() == self.dimension()

    @cached_method
    def is_central(self):
        r"""
        Test whether the intersection of all the hyperplanes is nonempty.

        OUTPUT:

        Boolean. Whether the hyperplane arrangement is such that the
        intersection of all the hyperplanes in the arrangement is
        nonempty.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)
            sage: a.is_central()
            True
        """
        R = self.base_ring()
        m = matrix(R, [h.normal() for h in self])
        b = vector(R, [h.b() for h in self])
        try:
            m.solve_right(b)
            return True
        except ValueError:
            return False

    @cached_method
    def essentialization(self):
        r"""
        Return the essentialization of the hyperplane arrangement.

        See the documentation for :meth:`essentialization`.

        OUTPUT:

        THe essentialization as a new hyperplane arrangement. If the
        characteristic of the base field is 0, this returns the
        hyperplane arrangement obtained by intersecting the
        hyperplanes by the space spanned by their normal vectors.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(3)
            sage: a.is_essential()
            False
            sage: a.essentialization()
            Arrangement <t1 - t2 | t1 + 2*t2 | 2*t1 + t2>

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: B = H([(1,0),1], [(1,0),-1])
            sage: B.is_essential()
            False
            sage: B.essentialization()
            Arrangement <-x + 1 | x + 1>
            sage: B.essentialization().parent()
            Hyperplane arrangements in 1-dimensional linear space over 
            Rational Field with coordinate x
          
            sage: H.<x,y> = HyperplaneArrangements(GF(2))
            sage: C = H([(1,1),1], [(1,1),0])
            sage: C.essentialization()
            Arrangement <y | y + 1>

            sage: h = hyperplane_arrangements.semiorder(4)
            sage: h.essentialization()
            Arrangement of 12 hyperplanes of dimension 3 and rank 3

        TESTS::

            sage: b = hyperplane_arrangements.coordinate(2)
            sage: b.is_essential()
            True
            sage: b.essentialization() is b
            True
        """
        def echelon_col_iter(row_iter):
            """helper to iterat over the echelon pivot column indices"""
            for row in row_iter:
                if row == 0:
                    raise StopIteration
                for pivot in range(self.dimension()):
                    if row[pivot] != 0:
                        break
                assert row[pivot] == 1
                yield pivot, row

        if self.is_essential():
            return self
        parent = self.parent()
        H = parent.ambient_space()
        R = parent.base_ring()
        hyperplanes = self.hyperplanes()
        normals = matrix(R, [h.normal() for h in self]).transpose()
        # find a (any) complement to the normals
        if R.characteristic() == 0:
            complement_basis = normals.kernel().echelonized_basis()
        else:
            # we don't necessarily have an orthogonal complement, pick any complement
            complement_basis = []
            for pivot, row in echelon_col_iter(normals.echelon_form().rows()):
                v = [0] * self.dimension()
                v[pivot] = 1
                complement_basis.append(vector(R, v))
        # reduce the hyperplane equations
        echelon_pivots = []   # the column indices where N has 1's from the echelonization
        for pivot, row in echelon_col_iter(complement_basis):
            assert row[pivot] == 1
            echelon_pivots.append(pivot)
            hyperplanes = [h - h.A()[pivot] * H(row, 0) for h in hyperplanes]
        # eliminate the pivot'ed coodinates
        restricted = []
        for h in hyperplanes:
            A = h.A()
            if A == 0:
                continue
            A = [A[i] for i in range(self.dimension()) if i not in echelon_pivots]
            b = h.b()
            restricted.append([A, b])
        names = tuple(name for i, name in enumerate(parent._names) if i not in echelon_pivots)
        # Construct the result
        restricted_parent = HyperplaneArrangements(R, names)
        return restricted_parent(*restricted, signed=False)

    def sign_vector(self, p):
        r"""
        Sign_vector indicates on which side of each hyperplane the given point `p` lies.

        The base field must have characteristic zero.

        INPUT:

        - ``p`` -- point as a list/tuple/iterable.

        OUTPUT:

        A vector whose entries are in `[-1, 0, +1]`.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,0), 0], [(0,1), 1]);  A
            Arrangement <y + 1 | x>
            sage: A.sign_vector([2, -2])
            (-1, 1)
            sage: A.sign_vector((-1, -1))
            (0, -1)

        TESTS::

            sage: H.<x,y> = HyperplaneArrangements(GF(3))
            sage: A = H(x, y)
            sage: A.sign_vector([1, 2])
            Traceback (most recent call last):
            ...
            ValueError: characteristic must be zero
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('characteristic must be zero')
        from sage.functions.generalized import sign
        p = vector(p)
        values = [p*h.A() + h.b() for h in self]
        return vector(ZZ, map(sign, values))
    
    def face_vector(self):
        r"""
        Return the face vector.

        OUTPUT:

        A vector of integers.

        The `d`-th entry is the number of faces of dimension `d`.  A
        *face* is is the intersection of a region with a hyperplane.

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.face_vector()
            (0, 6, 21, 16)
        """
        m = self.whitney_data()[0]
        v = list(sum(m.transpose().apply_map(abs)))
        v.reverse()
        v = [0]*(self.dimension() - self.rank()) + v
        return vector(ZZ, v)

    def _make_region(self, hyperplanes):
        """
        Helper method to construct a region

        INPUT:

        - ``hyperplanes`` -- a list/tuple/iterable of hyperplanes.

        OUTPUT:

        The polyhedron constructed from taking the linear expressions as inequalities.
        
        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: h = H(x)
            sage: h._make_region([x, 1-x, y, 1-y])
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        """
        ieqs = [h.coefficients() for h in hyperplanes]
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(ieqs=ieqs, ambient_dim=self.dimension(), 
                          base_ring=self.parent().base_ring())

    @cached_method
    def regions(self):
        r"""
        Return the regions of the hyperplane arrangement.

        The base field must have characteristic zero.

        OUTPUT:

        A tuple containing the regions as polyhedra.

        The regions are the connected components of the complement of
        the union of the hyperplanes as a subset of `\mathbb{R}^n`.

        EXAMPLES::

            sage: a = hyperplane_arrangements.braid(2)
            sage: a.regions()
            (A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex, 1 ray, 1 line, 
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex, 1 ray, 1 line)

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H(x, y+1)
            sage: A.regions()
            (A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays)
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('base field must have characteristic zero')
        num = self.n_hyperplanes()
        result = []
        from sage.misc.misc import powerset
        for pos in powerset(range(num)):
            ieqs = []
            for i in range(num):
                if i in pos:
                    ieqs.append(self[i])
                else:
                    ieqs.append(-self[i])
            P = self._make_region(ieqs)
            if P.dim() == self.dimension():
                result.append(P)
        return tuple(result)

    def region_containing_point(self, p):
        r"""
        The region in the hyperplane arrangement containing a given point.

        The base field must have characteristic zero.

        INPUT:

        - ``p`` -- point

        OUTPUT:

        A polyhedron. A ``ValueError`` is raised if the point is not
        interior to a region, that is, sits on a hyperplane.

        EXAMPLES::

            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: A = H([(1,0), 0], [(0,1), 1], [(0,1), -1], [(1,-1), 0], [(1,1), 0])
            sage: A.region_containing_point([1,2])
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 2 rays

        TESTS::

            sage: A = H([(1,1),0], [(2,3),-1], [(4,5),3])
            sage: B = A.change_ring(FiniteField(7))
            sage: B.region_containing_point((1,2))
            Traceback (most recent call last):
            ...
            ValueError: base field must have characteristic zero

            sage: A = H([(1,1),0], [(2,3),-1], [(4,5),3])
            sage: A.region_containing_point((1,-1))
            Traceback (most recent call last):
            ...
            ValueError: point sits on a hyperplane
        """
        if self.base_ring().characteristic() != 0:
            raise ValueError('base field must have characteristic zero')
        sign_vector = self.sign_vector(p)
        ieqs = []
        for i, hyperplane in enumerate(self):
            sign = sign_vector[i]
            if sign == 1:
                ieqs.append(hyperplane)
            elif sign == -1:
                ieqs.append(-hyperplane)
            else:
                assert sign == 0
                raise ValueError('point sits on a hyperplane')
        return self._make_region(ieqs)

    @cached_method
    def whitney_data(self):
        r"""
        Return the Whitney numbers.

        .. SEEALSO::

            :meth:`whitney_number`
            :meth:`doubly_indexed_whitney_number`

        OUTPUT:

        A pair of integer matrices. The two matrices are the
        doubly-indexed Whitney numbers of the first or second kind,
        respectively. The `i,j`-th entry is the `i,j`-th
        doubly-indexed Whitney number.

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.whitney_data()
            (
            [  1  -6   9]  [ 1  6  6]
            [  0   6 -15]  [ 0  6 15]
            [  0   0   6], [ 0  0  6]
            )
        """
        p = self.intersection_poset()
        r = p.rank_function()
        top = r(p.maximal_elements()[0])
        from sage.matrix.constructor import zero_matrix
        m1 = zero_matrix(ZZ, top+1, top+1)
        m2 = zero_matrix(ZZ, top+1, top+1)
        for i, j in p.relations_iterator():
            m1[r(i), r(j)] += p.mobius_function(i, j)
            m2[r(i), r(j)] += 1
        m1.set_immutable()
        m2.set_immutable()
        return (m1, m2)

    def doubly_indexed_whitney_number(self, i, j, kind=1):
        r"""
        Return the `i,j`-th  doubly-indexed Whitney number.

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

        Integer. The `(i,j)`-th entry of the ``kind`` Whitney number.

        EXAMPLES::

            sage: A = hyperplane_arrangements.Shi(3)
            sage: A.doubly_indexed_whitney_number(0, 2)
            9
            sage: A.whitney_number(2)
            9
            sage: A.doubly_indexed_whitney_number(1, 2)
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
        if 0 <= i and j <= self.dimension():
            if kind == 1:
                return self.whitney_data()[0][i, j]
            elif kind == 2:
                return self.whitney_data()[1][i, j]
        raise ValueError('argument out of range')

    def whitney_number(self, k, kind=1):
        r"""
        Return the ``k``-th Whitney number.

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

        Integer. The ``k``-th Whitney number.

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
        if k >= 0 and k <= self.dimension():
            if kind == 1:
                return self.whitney_data()[0][0, k]
            elif kind == 2:
                return self.whitney_data()[1][0, k]
        raise ValueError('argument out of range')


class HyperplaneArrangements(Parent, UniqueRepresentation):
    
    Element = HyperplaneArrangementElement

    def __init__(self, base_ring, names=tuple()):
        """
        Hyperplane arrangements

        EXAMPLES::
        
            sage: H.<x,y> = HyperplaneArrangements(QQ)
            sage: x
            Hyperplane x + 0*y + 0
            sage: x + y
            Hyperplane x + y + 0
            sage: H(x, y, x-1, y-1)
            Arrangement <y - 1 | y | x - 1 | x>

        TESTS::
          
            sage: H = HyperplaneArrangements(QQ, ('x', 'y'))
            sage: type(H)
            <class 'sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangements_with_category'>
            sage: H.change_ring(RR).gen(0)
            Hyperplane 1.00000000000000*x + 0.000000000000000*y + 0.000000000000000
        """
        from sage.categories.all import Fields, Sets
        if not base_ring in Fields: 
            raise ValueError('base ring must be a field')
        super(HyperplaneArrangements, self).__init__(category=Sets())
        self._base_ring = base_ring
        self._names = names

    def base_ring(self):
        return self._base_ring

    def change_ring(self, base_ring):
        """
        Return hyperplane arrangements over a different base ring.

        INPUT:

        - ``base_ring`` -- a ring. The new base ring.

        OUTPUT:
        
        A new :class:`HyperplaneArrangements` instance over the new
        base ring.

        EXAMPLES::

            sage: L = HyperplaneArrangements(QQ, ('x', 'y'))
            sage: L.gen(0)
            Hyperplane x + 0*y + 0
            sage: L.change_ring(RR).gen(0)
            Hyperplane 1.00000000000000*x + 0.000000000000000*y + 0.000000000000000      

        TESTS::

            sage: L.change_ring(QQ) is L
            True
        """
        return HyperplaneArrangements(base_ring, self._names)

    @cached_method
    def ambient_space(self):
        return AmbientVectorSpace(self.base_ring(), self._names)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 2-dimensional linear space over Rational Field with coordinates x, y
        """
        return 'Hyperplane arrangements in {0}'.format(self.ambient_space())

    def _element_constructor_(self, *args, **kwds):
        """
        INPUT:

        - ``*args`` -- positional arguments, each defining a
          hyperplane. Alternatively, a single polytope or a single
          hyperplane arrangement.

        - ``signed`` -- boolean keyword argument (optional, default:
          ``True``). Whether to preserve signs of hyperplane
          equations.

        - ``warn_duplicates`` -- boolean keyword argument (optional,
          default: ``False``). Whether to issue a warning if duplicate
          hyperplanes were passed. Note that duplicate hyperplanes are
          always removed, whether or not there is a warning shown.

        - ``check`` -- boolean keyword argument (optional, default:
          ``True``). Whether to perform argument checking.

        EXAMPLES::

            sage: L.<x, y> = HyperplaneArrangements(QQ)
            sage: L._element_constructor_(x, y)
            Arrangement <y | x>
            sage: L._element_constructor_([x, y])
            Arrangement <y | x>
            sage: L._element_constructor_([0, 1, 0], [0, 0, 1])
            Arrangement <y | x>
            sage: L._element_constructor_([[0, 1, 0], [0, 0, 1]])
            Arrangement <y | x>

            sage: L._element_constructor(polytopes.n_cube(2))
            Arrangement <-x + 1 | -y + 1 | y + 1 | x + 1>

            sage: L(x, x, warn_duplicates=True)
            doctest:...: UserWarning: Input contained 2 hyperplanes, but only 1 are distinct.
            Arrangement <x>
            sage: L(-x, x + y - 1, signed=False)
            Arrangement <-x - y + 1 | x>
        """
        if len(args) == 1 and isinstance(args[0], HyperplaneArrangementElement) and args[0].parent() is self:
            # optimization if argument is already a hyperplane arrangement
            return args[0]
        # process keyword arguments
        not_char2 = (self.base_ring().characteristic() != 2)
        signed = kwds.pop('signed', not_char2)
        warn_duplicates = kwds.pop('warn_duplicates', False)
        check = kwds.pop('check', True)
        if len(kwds) > 0:
            raise ValueError('unknown keyword argument')
        # process positional arguments
        AA = self.ambient_space()
        try:
            hyperplanes = map(AA, args)
        except (TypeError, ValueError, AttributeError):
            if len(args) > 1:
                raise
            arg = args[0]
            if hasattr(arg, 'Hrepresentation'):
                hyperplanes = [AA(h) for h in arg.Hrepresentation()]
            else:
                hyperplanes = map(AA, arg)
        hyperplanes = [h.primitive(signed) for h in hyperplanes]
        n = len(hyperplanes)
        hyperplanes = tuple(uniq(hyperplanes))
        if warn_duplicates and n != len(hyperplanes):
            from warnings import warn
            warn('Input contained {0} hyperplanes, but only {1} are distinct.'.format(n, len(hyperplanes)))
        # argument checking (optional but recommended)
        if check:
            if signed and not not_char2:
                raise ValueError('cannot be signed in characteristic 2')
            hyperplane_set = set(hyperplanes)
            for h in hyperplanes:
                if h.A() == 0:
                    raise ValueError('linear expression must be non-constant to define a hyperplane')
                if not_char2 and -h in hyperplane_set:
                    raise ValueError('arrangement cannot simultaneouly have h and -h as hyperplane')
        return self.element_class(self, hyperplanes)

    @cached_method
    def ngens(self):
        """
        Return the number of linear variables

        OUTPUT:
        
        Integer.

        EXAMPLES::

            sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 3-dimensional linear space over Rational Field with coordinates x, y, z
            sage: L.ngens()
            3
        """
        return len(self._names)

    @cached_method
    def gens(self):
        """
        Return the coordinate hyperplanes
        
        OUTPUT:

        A tuple of linear expressions, one for each linear variable.

        EXAMPLES::

            sage: L = HyperplaneArrangements(QQ, ('x', 'y', 'z'))
            sage: L.gens()
            (Hyperplane x + 0*y + 0*z + 0, 
             Hyperplane 0*x + y + 0*z + 0, 
             Hyperplane 0*x + 0*y + z + 0)
        """
        return self.ambient_space().gens()

    def gen(self, i):
        """
        Return the `i`-th coordinate hyperplane

        INPUT:

        - ``i`` -- integer.

        OUTPUT:

        A linear expression.

        EXAMPLES::

            sage: L.<x, y, z> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 3-dimensional linear space over Rational Field with coordinates x, y, z
            sage: L.gen(0)
            Hyperplane x + 0*y + 0*z + 0
        """
        return self.gens()[i]

    def _coerce_map_from_(self, P):
        """
        Return whether there is a coercion.
   
        TESTS::
        
            sage: L.<x> = HyperplaneArrangements(QQ);  L
            Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
            sage: M.<y> = HyperplaneArrangements(RR);  M
            Hyperplane arrangements in 1-dimensional linear space over Real Field with 53 bits of precision with coordinate y

            sage: L.coerce_map_from(ZZ)
            Conversion map:
              From: Integer Ring
              To:   Hyperplane arrangements in 1-dimensional linear space over Rational Field with coordinate x
            sage: M.coerce_map_from(L)
            Conversion map:
              From: Hyperplane arrangements in 1-dimensional linear space over
                    Rational Field with coordinate x
              To:   Hyperplane arrangements in 1-dimensional linear space over
                    Real Field with 53 bits of precision with coordinate y
            sage: L.coerce_map_from(M)
        """
        if self.ambient_space().has_coerce_map_from(P):
            return True
        if isinstance(P, HyperplaneArrangements):
            return self.base_ring().has_coerce_map_from(P.base_ring())
        return super(HyperplaneArrangements, self)._coerce_map_from_(P)


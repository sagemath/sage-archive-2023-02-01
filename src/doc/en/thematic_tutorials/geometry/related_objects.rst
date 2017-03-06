.. -*- coding: utf-8 -*-

.. linkall

.. _related_objects:

==============================================================
How to obtain related ...
==============================================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>

Once one constructed the polyhedron object, one would like to know some
combinatorial and geometric information about this object.


... enumerative properties
==============================================================

Dimensions
--------------------------------------------------------------

The ambient dimension is the dimension of the space in which the object is
defined:

::

    sage: P1 = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]])
    sage: P1.ambient_dim()
    2

.. end of output

Whereas the dimension of the object is the dimension of the smallest affine
subspace containing it.

::

    sage: Polyhedron(rays = [[1, 1]])
    A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
    and 1 ray
    sage: Polyhedron(rays = [[1, 1]]).dim()
    1
    sage: Polyhedron(rays = [[1, 1]]).dimension()
    1

.. end of output

:math:`f`-vector
--------------------------------------------------------------

The :math:`f`-vector contains the number of faces of the object ordered by
increasing dimension:

The cube has 8 vertices, 12 edges and 6 polygons:

::

    sage: Cube = polytopes.cube()
    sage: Cube.f_vector()
    (1, 8, 12, 6, 1)

.. end of output

One can also ask the :math:`f`-vector of unbounded polyhedron. :code:`P1` has 2
vertices and 3 edges.

::

    sage: P1.f_vector()
    (1, 2, 3, 1)

.. end of output

Neighborliness
--------------------------------------------------------------

The *neighborliness* is the highest cardinality :math:`k` for which all
:math:`k`-subsets of the vertices are faces of the polyhedron.

::

    sage: Cube.neighborliness()
    1
    sage: CP = polytopes.cyclic_polytope(5,9)
    sage: CP.neighborliness()
    2

.. end of output

Number of representation objects
--------------------------------------------------------------

The number of objects used in each representations is stored in 
several methods:

::

    sage: P1.n_Hrepresentation()  # The number of elements in the H-representation
    3
    sage: P1.n_Vrepresentation()  # The number of elements in the V-representation
    3

    sage: P1.n_equations()
    0
    sage: P1.n_inequalities()
    3
    sage: P1.n_lines()
    0
    sage: P1.n_rays()
    1

    sage: P1.n_vertices()
    2
    sage: P1.n_facets()
    3

.. end of output

... geometric objects and properties
==============================================================

Base ring
--------------------------------------------------------------

The first important object related to a polyhedron is its base ring.

::

    sage: P1.base_ring()
    Integer Ring

.. end of output

Bounded edges
--------------------------------------------------------------

The bounded edges are edges between two vertices of the polyhedron. This
function returns an iterator.

::

    sage: list(P1.bounded_edges())
    [(A vertex at (0, 1), A vertex at (1, 0))]

.. end of output

Bounding box
--------------------------------------------------------------

The bounding box returns the two corners (the ones with minimal and maximal 
coordinates) of a box with integer coordinates that contains the polyhedron.

::

    sage: Cube.bounding_box()
    ((-1, -1, -1), (1, 1, 1))

.. end of output

Center and Representative point
--------------------------------------------------------------

The :code:`center` returns the average of the vertices while the
:code:`representative_point` returns a point in the interior as far as it is
possible; if the polyhedron is not full dimensional a point in the relative
interior is returned.

::

    sage: P1.center()
    (1/2, 1/2)
    sage: P1.representative_point()
    (3/2, 3/2)

    sage: P2 = Polyhedron(vertices = [[0, 0], [3/2, 0], [3/2, 3/2], [0, 3]])
    sage: P2.representative_point()
    (3/4, 9/8)

.. end of output

Containment
--------------------------------------------------------------

Testing if a polyhedron contains a point is done as follows.

::

    sage: P3 = Polyhedron(vertices=[(2, 3), (3, 2), (2, 1), (1, 2)])
    sage: P3.interior_contains([2, 2])
    True
    sage: P3.interior_contains([2, 1])
    False
    sage: P3.contains([2, 1])
    True

    sage: P4 = Polyhedron(vertices = [[1/2, 0], [0, 1/2]])
    sage: P4.relative_interior_contains([1/4, 1/4])
    True
    sage: P4.interior_contains([1/4, 1/4])
    False

.. end of output

Ehrhart polynomial
--------------------------------------------------------------

The Ehrhart polynomial can be computed using the :code:`latte_int` package.

::

    sage: Cube.ehrhart_polynomial()  # optional - latte_int
    8*t^3 + 12*t^2 + 6*t + 1

.. end of output

Face and Normal fans
--------------------------------------------------------------

The *face fan* and the *normal fan* are two structures encoding geometrical
data of the polyhedron.

::

    sage: FaceFan(Cube)
    Rational polyhedral fan in 3-d lattice M
    sage: NormalFan(P3)
    Rational polyhedral fan in 2-d lattice N

.. end of output

Gale transform
--------------------------------------------------------------

The Gale transform -- also called *Gale dual* -- is useful to study polytopes
with few vertices. It allows to visualize polytopes and linear relations
between the vertices in a relatively small dimensional space.

::

    sage: CP = polytopes.cyclic_polytope(5,8)  # A 5-dim. polytope with 8 vertices
    sage: CP.gale_transform()
    [(1, 0), (0, 1), (-21, -6), (70, 15), (-105, -20), (84, 15), (-35, -6), (6, 1)]

.. end of ouput

Hyperplane arrangement
--------------------------------------------------------------

You can obtain the hyperplane arrangement given by the
:math:`H`-representation as an hyperplane arrangement object.

::

    sage: CP.hyperplane_arrangement()
    Arrangement of 30 hyperplanes of dimension 5 and rank 5

.. end of output

Integral points
--------------------------------------------------------------

You can count integer points as follows. The package :code:`latte_int` is
a useful addition in this kind of computations. You can install it by typing

.. CODE::

    sage -i latte_int

.. end of output

in a console.

::
    
    sage: Square = Polyhedron(vertices = [[1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1]])
    sage: Square.integral_points()
    ((1, -1, -1),
     (1, -1, 0),
     (1, -1, 1),
     (1, 0, -1),
     (1, 0, 0),
     (1, 0, 1),
     (1, 1, -1),
     (1, 1, 0),
     (1, 1, 1))
    sage: Square.integral_points_count()  # optional - latte_int
    9

.. end of output

Radius and radius square
--------------------------------------------------------------

The radius is the distance from the vertices to the center. All rays and lines
are ignored.

::

    sage: P1.radius()
    sqrt(1/2)
    sage: P1.radius_square()
    1/2

    sage: P2.radius()
    3/8*sqrt(29)
    sage: P2.radius_square()
    261/64

.. end of output

Corresponding linear program
--------------------------------------------------------------

If you would like to use some linear programming on your polyhedron object, use
the :code:`to_linear_program` method to obtain the corresponding linear program object.

::

    sage: P1.to_linear_program()
    Mixed Integer Program  ( maximization, 2 variables, 3 constraints )
    sage: P2.to_linear_program()
    Mixed Integer Program  ( maximization, 2 variables, 4 constraints )
    sage: P3.to_linear_program()
    Mixed Integer Program  ( maximization, 2 variables, 4 constraints )
    sage: P4.to_linear_program()
    Mixed Integer Program  ( maximization, 2 variables, 3 constraints )
    sage: CP.to_linear_program()
    Mixed Integer Program  ( maximization, 5 variables, 30 constraints )

.. end of output

Spaces
--------------------------------------------------------------

There are several spaces related to a polyhedron.

::

    sage: P1.ambient_space()
    Ambient free module of rank 2 over the principal ideal domain Integer Ring
    sage: P1.Hrepresentation_space()
    Ambient free module of rank 3 over the principal ideal domain Integer Ring
    sage: P1.Vrepresentation_space()
    Ambient free module of rank 2 over the principal ideal domain Integer Ring

.. end of output

Notice that the dimension of the :math:`H`-representation space is one more
than the ambient space.

Triangulation
--------------------------------------------------------------

You can triangulate a bounded polyhedron.

::

    sage: T = CP.triangulate()
    sage: for t in T:
    ....:     print t
    (0, 1, 2, 3, 4, 5)
    (0, 1, 2, 3, 5, 6)
    (0, 1, 2, 3, 6, 7)
    (0, 1, 2, 3, 7, 8)
    (0, 1, 3, 4, 5, 6)
    (0, 1, 3, 4, 6, 7)
    (0, 1, 3, 4, 7, 8)
    (0, 1, 4, 5, 6, 7)
    (0, 1, 4, 5, 7, 8)
    (0, 1, 5, 6, 7, 8)
    (1, 2, 3, 4, 5, 6)
    (1, 2, 3, 4, 6, 7)
    (1, 2, 3, 4, 7, 8)
    (1, 2, 4, 5, 6, 7)
    (1, 2, 4, 5, 7, 8)
    (1, 2, 5, 6, 7, 8)
    (2, 3, 4, 5, 6, 7)
    (2, 3, 4, 5, 7, 8)
    (2, 3, 5, 6, 7, 8)
    (3, 4, 5, 6, 7, 8)
    sage: type(T)
    <class 'sage.geometry.triangulation.element.PointConfiguration_with_category.element_class'>

.. end of output

.. note:: 

    If one is interested in studying the triangulations of a polytope, it is
    worth considering the class :ref:`sage.geometry.triangulation.point_configuration`.

Volume
--------------------------------------------------------------

The volume can be computed for full-dimensional bounded polyhedron. Setting
:code:`engine='lrs'` makes it possible to compute volumes of faces without
reducing the dimension of the ambient space.

::

    sage: P4.volume()
    0
    sage: CP.volume()
    1216512
    sage: Square.volume()
    0
    sage: Square.volume(engine='lrs')
    4.0
    sage: Cube.volume()
    8
    sage: Cube.volume(engine='lrs')
    8.0

.. end of output

... combinatorial objects and properties
==============================================================

Face lattice
--------------------------------------------------------------

One of the most important object related to a polyhedron is its *face lattice*
that records faces ordered by inclusion.

::

    sage: S = polytopes.simplex(3)
    sage: FL = S.face_lattice()
    sage: BL = posets.BooleanLattice(4)
    sage: FL.is_isomorphic(BL)
    True

.. end of output

.. note ::

    If one is interested in checking the combinatorial isomorphism of two
    polyhedron objects, one should look at the :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_combinatorially_isomorphic`. 

Facet and Vertex adjacency matrices
--------------------------------------------------------------

In order to know when two facets intersect or two vertices are contained in a
common face, one can looks at adjacency matrices.

::

    sage: Cube.facet_adjacency_matrix()
    [0 1 1 1 0 1]
    [1 0 1 1 1 0]
    [1 1 0 0 1 1]
    [1 1 0 0 1 1]
    [0 1 1 1 0 1]
    [1 0 1 1 1 0]

    sage: Cube.vertex_adjacency_matrix()
    [0 1 1 0 1 0 0 0]
    [1 0 0 1 0 1 0 0]
    [1 0 0 1 0 0 1 0]
    [0 1 1 0 0 0 0 1]
    [1 0 0 0 0 1 1 0]
    [0 1 0 0 1 0 0 1]
    [0 0 1 0 1 0 0 1]
    [0 0 0 1 0 1 1 0]
    sage: Cube.vertex_adjacency_matrix() == Cube.adjacency_matrix()
    True

.. end of output

Graph or 1-skeleton
--------------------------------------------------------------

The graph of a polyhedron consists of its vertices and edges.
For unbounded polyhedron, only the bounded edges are used.
There are two ways to get it.

::

    sage: K4 = graphs.CompleteGraph(4)
    sage: S.graph().is_isomorphic(K4)
    True
    sage: S.vertex_graph().is_isomorphic(K4)
    True

    sage: P1.graph()
    Graph on 2 vertices

.. end of output

Automorphic groups
--------------------------------------------------------------

The first one gives the automorphism group of the vertex graph of the polyhedron.

::

    sage: S.combinatorial_automorphism_group()
    Permutation Group with generators [(3,4), (2,3), (1,2)]

    sage: Octa = polytopes.octahedron()
    sage: Octa.combinatorial_automorphism_group()
    Permutation Group with generators [(3,4), (2,3)(4,5), (1,2)(5,6)]

.. end of output

The second automorphism group is the restricted automorphism group which
contains the affine transformations that preserve the :math:`V`-representation.

::

    sage: P5 = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1], [0, 1]])
    sage: P5.combinatorial_automorphism_group()
    Permutation Group with generators [(2,3)]
    sage: P5.restricted_automorphism_group()
    Permutation Group with generators [()]

.. end of output

Incidence matrix
--------------------------------------------------------------

The entries of the incidence matrix of a polyhedron object are indexed as

 - Rows :math:`\leftrightarrow` Vertices
 - Columns :math:`\leftrightarrow` Facets

There is a 1 when the corresponding vertex belongs to the corresponding facet
and a 0 otherwise.

::

    sage: Cube.incidence_matrix()
    [0 0 0 1 1 1]
    [1 0 0 1 0 1]
    [0 1 0 1 1 0]
    [1 1 0 1 0 0]
    [0 0 1 0 1 1]
    [1 0 1 0 0 1]
    [0 1 1 0 1 0]
    [1 1 1 0 0 0]

.. end of output

Vertex directed graph
--------------------------------------------------------------

Given a linear functional, sometimes also called an *objective function*, one
can give a direction to the edges in the graph of the polyhedron from the
smallest to the biggest value given by the functional (the default setup).

When two vertices have the same value, then two oriented edges are placed
between them. Checkout how :code:`G1` and :code:`G2` look like with the
:code:`plot` method.

::

    sage: G1 = Cube.vertex_digraph(vector([1,1,1]))
    sage: G1.sinks()
    [A vertex at (1, 1, 1)]
    sage: G2 = Cube.vertex_digraph(vector([1,1,0]))
    sage: G2.sinks()
    []
    sage: G2.sources()
    []

.. end of output

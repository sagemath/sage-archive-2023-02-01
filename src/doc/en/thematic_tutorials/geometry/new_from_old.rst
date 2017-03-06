.. -*- coding: utf-8 -*-

.. linkall

.. _new_from_old:

=========================================================
How to obtain ... new polyhedra from old ones
=========================================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>


It is possible to apply various constructions once one has a polyhedron object.
Here is a - not necessarily complete - list of operations.

Minkowski sums
=========================================================

It is possible to do Minkowski sums of polyhedron, using two syntaxes.

::

    sage: P1 = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]])
    sage: P2 = Polyhedron(vertices = [[1/2, 0], [0, 1/2]])
    sage: P1.Minkowski_sum(P2)
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 1 ray

    sage: P1 + P2
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 1 ray

.. end of output

Minkowski differences
=========================================================

After adding, one would like to substract:

::

    sage: Cube = polytopes.cube()
    sage: Square = Polyhedron(vertices = [[1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1]])
    
    sage: Cube.Minkowski_difference(Square)
    A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices
    sage: Square.Minkowski_difference(Cube)
    A 0-dimensional polyhedron in ZZ^3 defined as the convex hull of 1 vertex

    sage: Cube - Square
    A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices
    sage: Square - Cube
    A 0-dimensional polyhedron in ZZ^3 defined as the convex hull of 1 vertex
    
.. end of output

Minkowski decompositions
=========================================================

Given a 2-dimensional bounded polyhedron over :math:`\mathbb{Z}`, you can obtain it Minkowski
decompositions.

::

    sage: Square.Minkowski_decompositions()
    ((A 0-dimensional polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
      A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices),
     (A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
      A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices),
     (A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
      A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices),
     (A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
      A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices),
     (A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
      A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices))
    
.. end of output

Product
=========================================================

It is also possible to multiply polyhedron:

::

    sage: P1.product(P2)
    A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices and 1 ray

    sage: P1 * P2
    A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices and 1 ray

.. end of output

Intersection
=========================================================

Of course, it is possible to intersect two polyhedron objects:

::

    sage: P3 = Polyhedron(vertices = [[3, 0], [4, 1]], rays = [[-1, 1]])
    sage: P1.intersection(P3)
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices

    sage: P1_and_P3 = P1 & P3; P1_and_P3
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices

.. end of output

Union
=========================================================

It is also possible to take the *set theoretic* union of two polyhedron
objects. It does the union of vertices, rays and lines to form the convex hull
of the two objects.

::

    sage: R1 = Polyhedron(rays = [[-1]])
    sage: R2 = Polyhedron(rays = [[1]])
    sage: R1.convex_hull(R2)
    A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex and 1 line

    sage: P1_union_P3 = P1.convex_hull(P3)
    sage: P1_union_P3
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
    and 2 rays
    sage: P1_union_P3.vertices()
    (A vertex at (3, 0), A vertex at (1, 0))
    sage: P1_union_P3.rays()
    (A ray in the direction (-1, 1), A ray in the direction (1, 1))

.. end of output

Making it full-dimensional
=========================================================

Sometimes you deal with a low-dimensional polytope is a high-dimensional space.
There is a function to make an affine transformation and give an affinely
equivalent polyhedron which is full-dimensional.

::

    sage: Hexagon = Polyhedron(vertices = Permutations([1,2,3])); Hexagon
    A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
    sage: Hex_aff = Hexagon.affine_hull(); Hex_aff
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 6 vertices

..

Visualize :code:`Hex_aff` and you will see that it is not a regular hexagon.

Taking a face
=========================================================

It is possible to obtain each face of a polyhedron.

::

    sage: for f in P1.faces(1):
    ....:     print f.ambient_Vrepresentation()
    (A vertex at (0, 1), A ray in the direction (1, 1))
    (A vertex at (0, 1), A vertex at (1, 0))
    (A vertex at (1, 0), A ray in the direction (1, 1))

.. end of output

Faces remember the polyhedron it comes from and can also become a polyhedron
object on its own.

::

    sage: f = P1.faces(1)[0]
    sage: f.polyhedron() is P1
    True

    sage: f.as_polyhedron()
    A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 1 ray

.. end of output

Barycentric subdivision
=========================================================

What is the barycentric subdivision of the simplex?

::

    sage: S = polytopes.simplex(3); S
    A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
    sage: BS = S.barycentric_subdivision(); BS
    A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 14 vertices

.. end of output

Hint: it is the polar dual of a polytope in the library.

Bipyramid
=========================================================

The bipyramid is similar to the suspension in topology. It increases the
dimension of the polytope by 1.

::

    sage: Cube.bipyramid()
    A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 10
    vertices

.. end of output

Dilation
=========================================================

It is possible to dilate a polyhedron by an arbitrary scalar.

::

    sage: D_P1 = P1.dilation(AA(sqrt(2))); D_P1.vertices()
    (A vertex at (0, 1.414213562373095?), A vertex at (1.414213562373095?, 0))

    sage: P4 = Polyhedron(vertices = [[0, 0], [1, 0], [0, 1]])
    sage: 2*P4
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
    sage: P4.dilation(2) == 2*P4
    True

.. end of output

Face-truncation
=========================================================

It is possible to truncate a specific face of a polyhedron. One can also change
the angle of the truncation and how deep the cut is done. 

::

    sage: my_face = P1.faces(0)[0]  # This is a vertex-face
    sage: Trunc1_P1 = P1.face_truncation(P1.faces(0)[0])
    sage: Trunc1_P1.plot()
    Launched png viewer for Graphics object consisting of 6 graphics primitives
    sage: Trunc2_P1 = P1.face_truncation(P1.faces(0)[0],linear_coefficients=(1, 1/2), cut_frac=3/4)
    sage: Trunc2_P1.plot()
    Launched png viewer for Graphics object consisting of 6 graphics primitives

.. end of output

Lattice polytope
=========================================================

This method returns an encompassing lattice polytope.

::

    sage: LP = P2.lattice_polytope(envelope=True)  # envelope=True for rational polytopes
    sage: LP.vertices()
    M(0, 0),
    M(0, 1),
    M(1, 0)
    in 2-d lattice M

.. end of output

Polar
=========================================================

The polar polytope is only defined for compact, or bounded, polyhedron.

::

    sage: P2.polar()
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 1 line
    
    sage: P5 = Polyhedron(vertices = [[1/2, 0, 0], [0, 1/2, 0]],
    ....:                 rays = [[1, 1, 0]],
    ....:                 lines = [[0, 0, 1]])
    sage: P5.polar()
    Traceback (most recent call last):
    ...
    AssertionError: Not a polytope.

.. end of output

Prism
=========================================================

The prism construction is the same as taking the Minkowski sum of the
polyhedron with a segment (a 1-dimensional polytope) in an orthogonal space.

::

    sage: P1.prism()
    A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices and 1 ray

.. end of output

Pyramid
=========================================================

Similar, the pyramid is a join of a vertex with the polyhedron.

::

    sage: (P1_and_P3).pyramid()
    A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 5 vertices

.. end of output

Translation
=========================================================

One can translate a polyhedron by a vector.

::

    sage: (P1_and_P3).vertices()
    (A vertex at (2, 3),
     A vertex at (3, 2),
     A vertex at (2, 1),
     A vertex at (1, 2))
    sage: P1P3_translate = (P1_and_P3).translation([-1, 0])
    sage: P1P3_translate.vertices()
    (A vertex at (0, 2),
     A vertex at (1, 1),
     A vertex at (1, 3),
     A vertex at (2, 2))
    
    sage: P1_and_P3.find_translation(P17_translate)
    (-1, 0)
    sage: P1_and_P3.find_translation(P2)
    Traceback (most recent call last):
    ...
    ValueError: polyhedron is not a translation of self

.. end of output

Truncation
=========================================================

The truncation of a polyhedron is obtained by *chopping* all vertices
simultaneously.

::

    sage: TCube = Cube.truncation()
    sage: TCube2 = polytopes.truncated_cube()
    sage: TCube.is_combinatorially_isomorphic(TCube2)
    True

.. end of output

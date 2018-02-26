.. -*- coding: utf-8 -*-

.. linkall

.. _tips:

=========================
Polyhedra tips and tricks
=========================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>


Operation shortcuts
=================================================

You can obtain different operations using natural symbols:

::

    sage: Cube = polytopes.cube()
    sage: Octahedron = 3/2*Cube.polar()  # Dilation
    sage: Cube + Octahedron   # Minkowski sum
    A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices
    sage: Cube & Octahedron   # Intersection
    A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices
    sage: Cube * Octahedron   # Cartesian product
    A 6-dimensional polyhedron in QQ^6 defined as the convex hull of 48 vertices
    sage: Cube - Polyhedron(vertices=[[-1,0,0],[1,0,0]])  # Minkowski difference
    A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices

.. end of output

Sage input function
==============================================================

If you are working with a polyhedron that was difficult to construct
and you would like to get back the proper Sage input code to reproduce this
object, you can!

::

    sage: Cube = polytopes.cube()
    sage: TCube = Cube.truncation()
    sage: sage_input(TCube)
    Polyhedron(backend='ppl', base_ring=QQ, vertices=[(-QQ(1), -QQ(1), -1/3),
    (-QQ(1), -QQ(1), 1/3), (-QQ(1), -1/3, -QQ(1)), (-QQ(1), -1/3, QQ(1)),
    (-QQ(1), 1/3, -QQ(1)), (-QQ(1), 1/3, QQ(1)), (-QQ(1), QQ(1), -1/3),
    (-QQ(1), QQ(1), 1/3), (-1/3, -QQ(1), -QQ(1)), (-1/3, -QQ(1), QQ(1)), (-1/3,
    QQ(1), -QQ(1)), (-1/3, QQ(1), QQ(1)), (1/3, -QQ(1), -QQ(1)), (1/3, -QQ(1),
    QQ(1)), (1/3, QQ(1), -QQ(1)), (1/3, QQ(1), QQ(1)), (QQ(1), -QQ(1), -1/3),
    (QQ(1), -QQ(1), 1/3), (QQ(1), -1/3, -QQ(1)), (QQ(1), -1/3, QQ(1)), (QQ(1),
    1/3, -QQ(1)), (QQ(1), 1/3, QQ(1)), (QQ(1), QQ(1), -1/3), (QQ(1), QQ(1),
    1/3)])

.. end of output


:code:`repr_pretty_Hrepresentation`
==============================================================

If you would like to visualize the :math:`H`-representation nicely and even get
the latex presentation, there is a method for that!

::

    sage: Nice_repr = TCube.repr_pretty_Hrepresentation(separator='\n')
    sage: print(Nice_repr)
    1 >= x0
    1 >= x1
    3*x1 + 7 >= 3*x0 + 3*x2
    x0 + 1 >= 0
    x1 + 1 >= 0
    3*x0 + 7 >= 3*x1 + 3*x2
    3*x0 + 3*x1 + 7 >= 3*x2
    3*x0 + 3*x2 + 7 >= 3*x1
    3*x0 + 3*x1 + 3*x2 + 7 >= 0
    x2 + 1 >= 0
    1 >= x2
    3*x1 + 3*x2 + 7 >= 3*x0
    3*x2 + 7 >= 3*x0 + 3*x1
    7 >= 3*x0 + 3*x1 + 3*x2

    sage: Latex_repr = LatexExpr(TCube.repr_pretty_Hrepresentation(separator=",\\\\", latex=True))
    sage: view(Latex_repr)  # not tested

.. end of output

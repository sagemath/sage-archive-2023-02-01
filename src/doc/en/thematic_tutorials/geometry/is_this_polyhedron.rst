.. -*- coding: utf-8 -*-

.. linkall

.. _is_this_polyhedron:

==============================================================
Is this polyhedron ...
==============================================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>

Once a polyhedron object is constructed, it is often useful to check if it has certain 
properties.

Here is a list of properties that Sage can check.

.. note::

    The following list may note be complete due to recent additions of features. You can 
    check by making a :code:`tab` completion after typing :code:`is_` to see which methods
    are available.

... combinatorially isomorphic to this one?
==============================================================

Two polyhedra are *combinatorially isomorphic* if their face lattices are isomorphic.


Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_combinatorially_isomorphic`

... compact/a polytope?
==============================================================

A compact polyhedron is also called a *polytope*.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_compact`

... empty?
==============================================================

This function although sounding trivial is very important!

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_empty`

... full-dimensional?
==============================================================

A polyhedron is full dimensional when it does not have any equations in its
:math:`H`-representation.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_full_dimensional`

... a lattice polytope?
==============================================================

This functions checks for compactness and if the base ring is :math:`\mathbb{Z}`.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_lattice_polytope`

... inscribed on a sphere?
==============================================================

This functions checks if the vertices of the polyhedron lie on a sphere.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_inscribed`

... a Minkowski sum of this other one?
==============================================================

This functions checks if the polyhedron can be used to produce another using a
Minkowski sum.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_Minkowski_summand`

... neighborly?
==============================================================

This functions checks if the polyhedron had full skeleton until hafl of the
dimension.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_neighborly`

... reflexive?
==============================================================

This functions checks if the polar of a lattice polytope is also a lattice
polytope.

Reference manual: :meth:`sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.is_reflexive`

... simple? (Aren't they all?)
==============================================================

This functions checks if the degree of all vertices is equal to the dimension
of the polytope.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_simple`

... simplicial?
==============================================================

Simplicial polytopes have only simplices as faces.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_simplicial`

... the simplex?
==============================================================

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_simplex`

... the whole space?
==============================================================

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_universe`

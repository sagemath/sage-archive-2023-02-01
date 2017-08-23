.. -*- coding: utf-8 -*-
.. linkall

.. _polyhedra_quickref:

=====================================
Quick reference for polyhedra in Sage
=====================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>
                  Vincent Delecroix <vincent.delecroix@u-bordeaux.fr>

List of Polyhedron methods
==========================

**H and V-representation**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.base_ring` | ring on which the polyhedron is defined
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.ambient_space` | ambient vector space or free module
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.Hrepresentation_space` | vector space or free module used for the vectors of the H-representation
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.Vrepresentation_space` | vector space or free module used for the vectors of the V-representation
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_Hrepresentation` | number of elements in the H-representation (sum of the number of equations and inequalities)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_Vrepresentation` | number of elements in the V-representation (sum of vertices, rays and lines)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_equations` | number of equations
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_inequalities` | number of inequalities
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_vertices` | number of vertices
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_rays` | number of rays
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_lines` | number of lines
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.n_facets` | number of facets

**Polyhedron boolean properties:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_empty` | test emptyness
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_universe` | test whether a polyhedra is the whole ambient space
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_full_dimensional` | test if the polyhedron has the same dimension has the ambient space
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_combinatorially_isomorphic` | test whether two polyhedra are combinatorially isomorphic
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_compact` | test compactness of a polyhedron
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_lattice_polytope` | test whether a polyhedron is a lattice polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_inscribed` | test whether the polyhedron is inscribed in a sphere
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_Minkowski_summand` | WTH?test if the polyhedron can be used to produce another using a Minkowski sum.
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_neighborly` | test whether the polyhedron has full skeleton until half of the dimension
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.is_reflexive` | test if the polar of a lattice polytope is also a lattice polytope (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_simple` |  checks whether the degree of all vertices is equal to the dimension of the polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_simplex` | test whether a polytope is a simplex
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.is_simplicial` |  check whether a polyhedra has all its faces being simplices

**Enumerative properties**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.ambient_dim` |  the dimension of the ambient vector space
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.dim` |  the dimension of the polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.dimension` |  alias of dim
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.f_vector` |  the `f`-vector (number of faces of each dimension)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.neighborliness` | highest cardinality for which all `k`-subsets of the vertices are faces of the polyhedron

**Transforming polyhedra**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.Minkowski_sum` | Minkowski sum of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.Minkowski_difference` | Minkowski difference of two polyhedra
    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.Minkowski_decompositions` | Minkowski decomposition (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.product` | cartesian product of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.intersection` | intersection of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.convex_hull` | convex hull of the union of two polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.affine_hull` | construct an affinely equivalent full dimensional polyhedra
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.faces` | the list of faces
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.barycentric_subdivision` | the list of faces
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.dilation` |  scalar dilation
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_truncation` | truncate a specific face
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.lattice_polytope` | return an encompassing lattice polytope.
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.polar` | return the polar of a polytope (needs to be compact)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.prism` | the prism of a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient space)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.pyramid` | the pyramid of a polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient space)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.bipyramid` | construct the bipyramid of te polyhedron (increases both the dimension of the polyhedron and the dimension of the ambient space by 1)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.translation` | translation by a given vector
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.truncation` | truncate all vertices

**Combinatorics**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_lattice` | the face lattice
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.combinatorial_automorphism_group` | the automorphism group of the underlying combinatorial polytope
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.graph`, :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_graph` | underlying graph
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_digraph` | digraph (orientation of edges determined by a linear form)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.adjacency_matrix` | adjacency matrix
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix` | incidence matrix
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.facet_adjacency_matrix` | adjacency matrix of the facets
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_adjacency_matrix` | adjacency matrix of the vertices

**Integral points**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ.ehrhart_polynomial` | the Ehrhart polynomial (only for :class:`Polyhedron over ZZ <sage.geometry.polyhedron.base_ZZ.Polyhedron_ZZ>`)
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.integral_points` | list of integral points
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.integral_points_count` | number of integral points

**Other**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.bounded_edges` | generator for bouded edges
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.bounding_box` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.center` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.representative_point` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.interior_contains` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.contains` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.relative_interior_contains` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.gale_transform` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.hyperplane_arrangement` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.radius` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.radius_square` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.to_linear_program` | transform the polyhedra into a Linear Program
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.radius_square` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.triangulate` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.volume` | ?
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.restricted_automorphism_group` | ?
    :meth:`~sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class.lattice_automorphism_group` | Only for :class:`PPL Lattice Polytope <sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class>`


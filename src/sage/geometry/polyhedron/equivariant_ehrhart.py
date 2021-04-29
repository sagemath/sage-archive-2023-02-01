r"""
Script for Equivariant Ehrhart Theory computations.

This script is centered around the computation of the equivariant `H^*` function
for a polytope `P` which is invariant under a linear group action using the
function :func:`Hstar_function`.
The equivariant Ehrhart series is the formal power series
`\sum_{m \geq 0} \chi_{mP}t^m`, where `\chi_{mP}` is
the permutation character of the group action on the lattice points in the
`m`th dilate of `P`. The equivariant Ehrhart series has the rational generating
function
`\sum_{m \geq 0} \chi_{mP}t^m = \frac{H^*(t)} {(1-t)det(Id-\rho(t)}.`
For details, see [Stap2011]_.
There are also functions to test whether the `H^*` series is polynomial,
:func:`is_element_of_base_ring`, and effective, :func:`is_effective`.

One can also use this module to compute fixed subpolytopes of a polytope under
the action of a group with the functions :func:`fixed_subpolytopes` and
:func:`fixed_subpolytope`.`
Finally, the module contains the function :func:`match_perms_to_mats` which
returns a dictionary between elements of the group written as permutations of a
polytope's vertices and written as matrices.

REFERENCES:

    - [Stap2011]_ Alan Stapledon. Equivariant Ehrhart Theory. Advances in Mathematics 226 (2011), no. 4, 3622-3654
    - [ASV2020]_ Federico Ardila, Mariel Supina, and Andr\ ́{e}s R. Vindas-Mel\ ́{e}ndez, The Equivariant Ehrhart Theory of the Permutahedron, 2020.

EXAMPLES:

A first example follows Prop 6.1 of [Stap2011]_  which states that for a simplex,
`H^*` is given by the permutation character of the group on the box points at
each height. The one dimensional simplex [(-1,1),(1,1)] has two box points,
one at (0,0) and the other at (0,1). The group Z/2Z acts linearly on the simplex
by reflection across the y-axis. The box points are fixed under this action,
so we expect an `H^*` polynomial of `\chi_triv + \chi_triv t`::

    sage: simplex = Polyhedron(vertices = [[-1,1],[1,1]], backend = 'normaliz') # optional - pynormaliz
    sage: G = simplex.restricted_automorphism_group(output = 'permutation')     # optional - pynormaliz
    sage: G.order()   # optional - pynormaliz
    2
    sage: Hstar_function(simplex,G) # optional - pynormaliz
    chi_1*t + chi_1
    sage: G.character_table()       # optional - pynormaliz
    [ 1 -1]
    [ 1  1]

A second example is the action of the symmetric group on the 2-dimensional
permutahedron in 3-dimensional space, given by permuting the three basis
vectors. As shown in [ASV2020]_, the corresponding Hstar function is polynomial and
effective, while for all higher dimensional permutahedra under the symmetric
group action, it is not. The following computation agrees with the
`H^*`-polynomial recovered in [ASV2020]_::

    sage: p2 = polytopes.permutahedron(3, backend = 'normaliz')      # optional - pynormaliz
    sage: G = p2.restricted_automorphism_group(output='permutation') # optional - pynormaliz
    sage: H = G.subgroup(gens=[G.gens()[1],G.gens()[2]])             # optional - pynormaliz
    sage: H.order()                                                  # optional - pynormaliz
    6
    sage: [Hstar, Hlin] = [Hstar_function(p2,H), Hstar_function(p2,H, output = 'Hstar_as_lin_comb')]; Hstar # optional - pynormaliz
    chi_0*t^2 + (chi_0 + chi_1 + chi_2)*t + chi_0
    sage: H.character_table()        # optional - pynormaliz
    [ 1  1  1]
    [ 1 -1  1]
    [ 2  0 -1]
    sage: is_element_of_base_ring(Hstar)       # optional - pynormaliz
    True
    sage: is_effective(Hstar,Hlin)   # optional - pynormaliz
    True

Example of a simplex under symmetric group action::

    sage: S = polytopes.simplex(3); S
    A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
    sage: G = S.restricted_automorphism_group(output = 'permutation');
    sage: len(G)
    24
    sage: Hstar_simplex = Hstar_function(S,G); Hstar_simplex # optional - pynormaliz
    chi_4
    sage: is_element_of_base_ring(Hstar_simplex) # optional - pynormaliz
    True

The following is example 7.6 in Stapledon::

    sage: P = Polyhedron(vertices=[[0,0,1],[0,0,-1],[1,0,1],[-1,0,-1],[0,1,1],
    ....: [0,-1,-1],[1,1,1],[-1,-1,-1]],backend='normaliz')           # optional - pynormaliz
    sage: G = P.restricted_automorphism_group(output = 'permutation') # optional - pynormaliz
    sage: H = G.subgroup(gens = [G[6]])                               # optional - pynormaliz
    sage: Hstar = Hstar_function(P,H); Hstar                          # optional - pynormaliz
    (chi_0*t^4 + (3*chi_0 + 3*chi_1)*t^3 + (8*chi_0 + 2*chi_1)*t^2 + (3*chi_0 + 3*chi_1)*t + chi_0)/(t + 1)
    sage: is_element_of_base_ring(Hstar)                  # optional - pynormaliz
    False

AUTHORS:

- Sophia Elia (2021): Initial version
- Jean-Philippe Labbé (2021): Initial version
"""

##############################################################################
#     Copyright (C) 2021 Sophia Elia         <sophiae56 at math.fu-berlin.de>
#                   2021 Jean-Philippe Labbe <labbe at math.fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################










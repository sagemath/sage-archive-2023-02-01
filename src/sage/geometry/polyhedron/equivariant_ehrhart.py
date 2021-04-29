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





def match_perms_to_mats(polytope, conj_class_reps, acting_group=None, additional_elts=None):
    r"""
    An element of ``polytope``'s ``restricted_autormorphism_group`` may
    be represented either as a permutation of the vertices of ``polytope``
    or as a matrix. This function returns a dictionary with permutations
    as keys and matrices as values.

    When ``additional_elts`` is ``None``, the dictionary is returned for the
    generators and conjugacy classes representatives in conj_class_reps of the
    ``restricted_automorphism_group`` or the ``acting_group``.
    When ``additional_elts`` is not ``None``, each element in
    ``additional_elts`` also becomes a key.

    INPUT:

    - ``polytope`` -- polyhedron object. a lattice polytope.

    - ``conj_class_reps`` -- list. A list of representatives of the conjugacy
      classes of the ``restricted_automorphism_group``.

    - ``acting_group`` -- a subgroup of the ``polytope``'s
      ``restricted_automorphism_group``.

    - ``additional_elts`` -- list (default=None). a subset of the
      ``restricted_automorphism_group`` of ``polytope`` expressed as
      permutations.

    OUTPUT:

    A dictionary between elements of ``the restricted_automorphism_group``
    or ``acting_group`` expressed as permutations (keys) and matrices (values).

    EXAMPLES:

    This example shows the dictionary between permutations and matrices
    for the generators of the ``restricted_automorphism_group`` of the
    `\pm 1` 2-dimensional square::

        sage: square = Polyhedron(vertices=[[1,1],[-1,1],[-1,-1],[1,-1]], backend='normaliz') # optional - pynormaliz
        sage: aut_square = square.restricted_automorphism_group(output = 'permutation')       # optional - pynormaliz
        sage: conj_reps = aut_square.conjugacy_classes_representatives()                      # optional - pynormaliz
        sage: gens_dict = match_perms_to_mats(square,conj_reps); gens_dict                    # optional - pynormaliz
        {(): [1 0 0]
         [0 1 0]
         [0 0 1],
         (1,2): [0 1 0]
         [1 0 0]
         [0 0 1],
         (0,1)(2,3): [ 1  0  0]
         [ 0 -1  0]
         [ 0  0  1],
         (0,1,3,2): [ 0  1  0]
         [-1  0  0]
         [ 0  0  1],
         (0,3): [ 0 -1  0]
         [-1  0  0]
         [ 0  0  1],
         (0,3)(1,2): [-1  0  0]
         [ 0 -1  0]
         [ 0  0  1]}
        sage: square.vertices() # optional - pynormaliz
        (A vertex at (-1, -1),
        A vertex at (-1, 1),
        A vertex at (1, -1),
        A vertex at (1, 1))

    This example tests the functionality for extra elements::

        sage: C = polytopes.cross_polytope(2)
        sage: G = C.restricted_automorphism_group(output = 'permutation')
        sage: conj_reps = G.conjugacy_classes_representatives()
        sage: add_elt = [G[6]]; add_elt
        [(0,2,3,1)]
        sage: match_perms_to_mats(C,conj_reps,additional_elts = add_elt)
        {(): [1 0 0]
         [0 1 0]
         [0 0 1],
         (1,2): [ 1  0  0]
         [ 0 -1  0]
         [ 0  0  1],
         (0,1)(2,3): [0 1 0]
         [1 0 0]
         [0 0 1],
         (0,1,3,2): [ 0 -1  0]
         [ 1  0  0]
         [ 0  0  1],
         (0,2,3,1): [ 0  1  0]
         [-1  0  0]
         [ 0  0  1],
         (0,3): [-1  0  0]
         [ 0  1  0]
         [ 0  0  1],
         (0,3)(1,2): [-1  0  0]
         [ 0 -1  0]
         [ 0  0  1]}
    """
    V = [v.homogeneous_vector() for v in polytope.Vrepresentation()]
    Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()
    Vplus = list(matrix(V) * Qplus)
    W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

    G = polytope.restricted_automorphism_group(output='permutation')
    if acting_group is not None:
        G = acting_group

    group_dict = {}

    for perm in G.gens():
        group_dict[perm] = _match_perm(perm, V, Vplus, W)

    for perm in conj_class_reps:
        group_dict[perm] = _match_perm(perm, V, Vplus, W)

    if additional_elts is not None:
        for perm in additional_elts:
            group_dict[perm] = _match_perm(perm, V, Vplus, W)
    return group_dict


def _match_perm(permutation, V, Vplus, W):
    r"""
    Return the matrix representation of a permutation in the
    ``restricted_autormorphism_group`` of ``polytope``.

    INPUT:

    - ``polytope`` -- polyhedron object. A lattice polytope.

    - ``V`` -- list. a list of vectors from the ``match_perms_to_mats`` function.

    - ``Vplus`` -- list. from the ``match_perms_to_mats`` function.

    - ``W`` -- matrix. from the ``match_perms_to_mats`` function.

    OUTPUT:

    A matrix that acts in the same way on the polytope as the ``permutation``.

    EXAMPLES:

    This example shows a reflection across the x-axis::

        sage: cross = polytopes.cross_polytope(2, backend = 'normaliz') # optional - pynormaliz
        sage: cross.vertices() # optional - pynormaliz
        (A vertex at (-1, 0),
         A vertex at (0, -1),
         A vertex at (0, 1),
         A vertex at (1, 0))
        sage: G = cross.restricted_automorphism_group(output = 'permutation') # optional - pynormaliz
        sage: flip = G.gens()[0]; flip   # optional - pynormaliz
        (1,2)
        sage: V = [v.homogeneous_vector() for v in cross.Vrepresentation()]   # optional - pynormaliz
        sage: Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()    # optional - pynormaliz
        sage: Vplus = list(matrix(V) * Qplus)   # optional - pynormaliz
        sage: W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V))) # optional - pynormaliz
        sage: _match_perm(flip, V, Vplus, W)   # optional - pynormaliz
        [ 1  0  0]
        [ 0 -1  0]
        [ 0  0  1]
    """
    A = sum(V[permutation(i)].column() * Vplus[i].row() for i in range(len(V)))
    return A + W




def is_effective(Hstar, Hstar_as_lin_comb):
    r"""
    Check if the `H^*` series is effective.

    The `H^*` series is effective if it is a polynomial and the coefficients
    of each `t^i` are effective representations. The coefficients of each
    irreducible representation must be non-negative integers. The `H^*` series
    must be polynomial in order for it to be effective.

    INPUT:

    - ``Hstar`` -- a rational function in `t` with coefficients in the ring of
                   class functions.

    OUTPUT:

    Boolean. Whether the `H^*` series is effective.

    EXAMPLES:

    The `H^*` series of the two-dimensional permutahedron under the action
    of the symmetric group is effective::

        sage: p2 = polytopes.permutahedron(3, backend = 'normaliz')      # optional - pynormaliz
        sage: G = p2.restricted_automorphism_group(output='permutation') # optional - pynormaliz
        sage: H = G.subgroup(gens=[G.gens()[1],G.gens()[2]])             # optional - pynormaliz
        sage: H.order()                                                  # optional - pynormaliz
        6
        sage: [Hstar, Hlin] = [Hstar_function(p2,H), Hstar_function(p2,H, output = 'Hstar_as_lin_comb')] # optional - pynormaliz
        sage: is_effective(Hstar,Hlin)   # optional - pynormaliz
        True

    The `H^*` series must be a polynomial in order to be effective. If it is
    not polynomial, a value error is returned::

        sage: P = Polyhedron(vertices=[[0,0,1],[0,0,-1],[1,0,1],[-1,0,-1],[0,1,1],
        ....: [0,-1,-1],[1,1,1],[-1,-1,-1]],backend='normaliz')           # optional - pynormaliz
        sage: G = P.restricted_automorphism_group(output = 'permutation') # optional - pynormaliz
        sage: H = G.subgroup(gens = [G[6]])                               # optional - pynormaliz
        sage: Hstar = Hstar_function(P,H); Hstar                          # optional - pynormaliz
        (chi_0*t^4 + (3*chi_0 + 3*chi_1)*t^3 + (8*chi_0 + 2*chi_1)*t^2 + (3*chi_0 + 3*chi_1)*t + chi_0)/(t + 1)
        sage: Hstar_lin = Hstar_function(P,H, output = 'Hstar_as_lin_comb') # optional - pynormaliz
        sage: is_effective(Hstar, Hstar_lin)  # optional - pynormaliz
        Traceback (most recent call last):
        ...
        ValueError: The Hstar vector must be polynomial
    """
    if not is_element_of_base_ring(Hstar):
        raise ValueError("The Hstar vector must be polynomial")
    flag = True
    for irrep in range(len(Hstar_as_lin_comb)):
        coeffs = Hstar_as_lin_comb[irrep].numerator().coefficients()
        for i in coeffs:
            if not i.is_integral() or i < 0:
                flag = False
    return flag

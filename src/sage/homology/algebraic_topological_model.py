# -*- coding: utf-8 -*-
"""
Algebraic topological model for a cell complex

AUTHORS:

- John H. Palmieri (2015-09)
"""

########################################################################
#       Copyright (C) 2015 John H. Palmieri <palmieri@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
########################################################################

# TODO: cythonize this.

from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, zero_matrix
from chain_complex import ChainComplex
from chain_complex_morphism import ChainComplexMorphism
from chain_homotopy import ChainContraction
from sage.rings.rational_field import QQ

def algebraic_topological_model(K, base_ring=None):
    r"""
    Algebraic topological model for cell complex ``K``
    with coefficients in ``base_ring``. ``base_ring`` must be a field.

    INPUT:

    - ``K`` -- either a simplicial complex or a cubical complex.
    - ``base_ring`` -- coefficient ring (optional, default ``QQ``).
      Must be a field.

    OUTPUT: a pair ``(phi, M)`` consisting of

    - the chain contraction ``phi`` associated to `C`, `M`, `\pi: C
      \to M`, and `\iota: M \to C` satisfying `\iota \pi = 1_M` and
      `\pi \iota`.
    - the chain complex `M`

    This construction appears in a paper by Pilarczyk and Réal [PR]_.
    Given a cell complex `K` and a field `F`, there is a chain complex
    `C` associated to `K` with coefficients in `F`. The *algebraic
    topological model* for `K` is a chain complex `M` with trivial
    differential, along with chain maps `\pi: C \to M` and `\iota: M
    \to C` such that

    - `\pi \circ \iota = 1_M`.
    - There is a chain homotopy `\phi` between `1_C` and `\iota \circ \pi`.

    In particular, `\pi` and `\iota` induce isomorphisms on homology,
    and since `M` has trivial differential, it is its own
    homology. Thus `\iota` lifts homology classes to their cycle
    representatives.

    The chain homotopy `\phi` satisfies some additional properties,
    making it a *chain contraction*:

    - `\phi \circ \phi = 0`,
    - `\pi \phi = 0`,
    - `\phi \iota = 0`.

    Given an algebraic topological model for `K`, it is then easy to
    compute cup products and cohomology operations on the cohomology
    of `K`, as described in [G-DR]_ and [PR]_.

    Implementation details: the cell complex `K` must have a
    :meth:`cell_complex.GenericCellComplex.n_cells` method from which
    we can extract a list of cells in each dimension. Combining the
    lists in increasing order of dimension then defines a filtration
    of the complex: a list of cells in which the boundary of each cell
    consists of cells earlier in the list. This is required by
    Pilarczyk and Réal's algorithm.  There must also be a
    :meth:`cell_complex.GenericCellComplex.chain_complex` method, to
    construct the chain complex `C` associated to this chain complex.

    In particular, this should work on simplicial complexes and
    cubical complexes. It doesn't work for delta complexes, though:
    the list of their n-cells has the wrong format.

    Note that from the chain contraction ``\phi``, one can recover the
    chain maps `\pi` and `\iota` via ``phi.pi()`` and
    ``phi.iota()``. Then one can recover `C` and `M` from, for
    example, ``phi.pi()._domain`` and ``phi.pi()._codomain``,
    respectively.

    EXAMPLES::

        sage: from sage.homology.algebraic_topological_model import algebraic_topological_model
        sage: RP2 = simplicial_complexes.RealProjectivePlane()
        sage: phi, M = algebraic_topological_model(RP2, GF(2))
        sage: M.homology()
        {0: Vector space of dimension 1 over Finite Field of size 2,
         1: Vector space of dimension 1 over Finite Field of size 2,
         2: Vector space of dimension 1 over Finite Field of size 2}
        sage: T = cubical_complexes.Torus()
        sage: phi, M = algebraic_topological_model(T, QQ)
        sage: M.homology()
        {0: Vector space of dimension 1 over Rational Field,
         1: Vector space of dimension 2 over Rational Field,
         2: Vector space of dimension 1 over Rational Field}

    If you want to work with cohomology rather than homology, just
    dualize the outputs of this function::

        sage: M.dual().homology()
        {0: Vector space of dimension 1 over Rational Field,
         1: Vector space of dimension 2 over Rational Field,
         2: Vector space of dimension 1 over Rational Field}
        sage: M.dual().degree_of_differential()
        1
        sage: phi.dual()
        Chain homotopy between Chain complex morphism from Chain complex with at most 3 nonzero terms over Rational Field to Chain complex with at most 3 nonzero terms over Rational Field and Chain complex morphism from Chain complex with at most 3 nonzero terms over Rational Field to Chain complex with at most 3 nonzero terms over Rational Field

    In degree 0, the inclusion of the homology `M` into the chain
    complex `C` sends the homology generator to a single vertex::

        sage: K = simplicial_complexes.Simplex(2)
        sage: phi, M = algebraic_topological_model(K, QQ)
        sage: phi.iota().in_degree(0)
        [0]
        [0]
        [1]

    In cohomology, though, one needs the dual of every degree 0 chain
    to detect the degree 0 cohomology generator::

        sage: phi.dual().iota().in_degree(0)
        [1]
        [1]
        [1]

    TESTS::

        sage: T = cubical_complexes.Torus()
        sage: C = T.chain_complex()
        sage: H, M = T.algebraic_topological_model()
        sage: C.differential(1) * H.iota().in_degree(1).column(0) == 0
        True
        sage: C.differential(1) * H.iota().in_degree(1).column(1) == 0
        True
        sage: coC = T.chain_complex(cochain=True)
        sage: coC.differential(1) * H.dual().iota().in_degree(1).column(0) == 0
        True
        sage: coC.differential(1) * H.dual().iota().in_degree(1).column(1) == 0
        True
    """
    if base_ring is None:
        base_ring = QQ
    if not base_ring.is_field():
        raise ValueError('the coefficient ring must be a field.')

    # The following are all dictionaries indexed by dimension.
    # For each n, gens[n] is an ordered list of the n-cells generating the complex M.
    gens = {}
    # For each n, phi_dict[n] is a dictionary of the form {idx:
    # vector}, where idx is the index of an n-cell in the list of
    # n-cells in K, and vector is the image of that n-cell, as an
    # element in the free module of (n+1)-chains for K.
    phi_dict = {}
    # For each n, pi_dict[n] is a dictionary of the same form, except
    # that the target vectors should be elements of the chain complex M.
    pi_dict = {}
    # For each n, iota_dict[n] is a dictionary of the form {cell:
    # vector}, where cell is one of the generators for M and vector is
    # its image in C, as an element in the free module of n-chains.
    iota_dict = {}

    for n in range(K.dimension()+1):
        gens[n] = []
        phi_dict[n] = {}
        pi_dict[n] = {}
        iota_dict[n] = {}

    C = K.chain_complex(base_ring=base_ring)
    # old_cells: cells one dimension lower.
    old_cells = []

    for dim in range(K.dimension()+1):
        n_cells = K.n_cells(dim)
        # No need to set zero values for any of the maps: we will
        # assume any unset values are zero.
        # From the paper: phi_dict[c] = 0.
        diff = C.differential(dim)
        rank = len(n_cells)
        old_rank = len(old_cells)

        for c_idx, c in enumerate(n_cells):
            c_vec = vector(base_ring, rank, {c_idx: 1})

            # c_bar = c - phi(bdry(c))
            c_bar = c_vec
            bdry_c = diff * c_vec
            for (idx, coord) in bdry_c.iteritems():
                try:
                    c_bar -= coord * phi_dict[dim-1][idx]
                except KeyError:
                    pass

            bdry_c_bar = diff * c_bar

            # Evaluate pi(bdry(c_bar)).
            pi_bdry_c_bar = vector(base_ring, old_rank)

            for (idx, coeff) in bdry_c_bar.iteritems():
                try:
                    pi_bdry_c_bar += coeff * pi_dict[dim-1][idx]
                except KeyError:
                    pass

            # One small typo in the published algorithm: it says
            # "if bdry(c_bar) == 0", but should say
            # "if pi(bdry(c_bar)) == 0".
            if not pi_bdry_c_bar:
                # Append c to list of gens.
                gens[dim].append(c)
                # iota(c) = c_bar
                iota_dict[dim][c] = c_bar
                # pi(c) = c
                pi_dict[dim][c_idx] = c_vec
            else:
                # Take any u in gens so that lambda_i = <u, pi(bdry(c_bar))> != 0.
                # u_idx will be the index of the corresponding cell.
                for (u_idx, lambda_i) in pi_bdry_c_bar.iteritems():
                    # Now find the actual cell.
                    u = old_cells[u_idx]
                    if u in gens[dim-1]:
                        break

                # pi(c) = 0: no need to do anything about this.
                for c_j_idx, c_j in enumerate(old_cells):
                    # c_j_vec: representative of c_j in the free
                    # module of chains in dimension dim-1.
                    c_j_vec = vector(base_ring, old_rank, {c_j_idx: 1})
                    # eta_ij = <u, pi(c_j)>.
                    try:
                        eta_ij = pi_dict[dim-1][c_j_idx][u_idx]
                    except (KeyError, IndexError):
                        eta_ij = 0
                    if eta_ij:
                        # Adjust phi(c_j).
                        try:
                            phi_dict[dim-1][c_j_idx] += eta_ij * lambda_i**(-1) * c_bar
                        except KeyError:
                            phi_dict[dim-1][c_j_idx] = eta_ij * lambda_i**(-1) * c_bar
                        # Adjust pi(c_j).
                        try:
                            pi_dict[dim-1][c_j_idx] += -eta_ij * lambda_i**(-1) * pi_bdry_c_bar
                        except KeyError:
                            pi_dict[dim-1][c_j_idx] = -eta_ij * lambda_i**(-1) * pi_bdry_c_bar

                gens[dim-1].remove(u)
                del iota_dict[dim-1][u]
        old_cells = n_cells

    # Now we have constructed the raw data for M, pi, iota, phi, so we
    # have to convert that to data which can be used to construct chain
    # complexes, chain maps, and chain contractions.

    # M_data will contain (trivial) matrices defining the differential
    # on M. Keep track of the sizes using "M_rows" and "M_cols", which are
    # just the ranks of consecutive graded pieces of M.
    M_data = {}
    M_rows = 0
    # pi_data: the matrices defining pi. Similar for iota_data and phi_data.
    pi_data = {}
    iota_data = {}
    phi_data = {}
    for n in range(K.dimension()+1):
        n_cells = K.n_cells(n)
        # Remove zero entries from pi_dict and phi_dict.
        pi_dict[n] = {i: pi_dict[n][i] for i in pi_dict[n] if pi_dict[n][i]}
        phi_dict[n] = {i: phi_dict[n][i] for i in phi_dict[n] if phi_dict[n][i]}
        # Convert gens to data defining the chain complex M with
        # trivial differential.
        M_cols = len(gens[n])
        M_data[n] = zero_matrix(base_ring, M_rows, M_cols)
        M_rows = M_cols
        # Convert the dictionaries for pi, iota, phi to matrices which
        # will define chain maps and chain homotopies.
        pi_cols = []
        phi_cols = []
        for (idx, c) in enumerate(n_cells):
            # First pi:
            if idx in pi_dict[n]:
                column = vector(base_ring, M_rows)
                for (entry, coeff) in pi_dict[n][idx].iteritems():
                    # Translate from cells in n_cells to cells in gens[n].
                    column[gens[n].index(n_cells[entry])] = coeff
            else:
                column = vector(base_ring, M_rows)
            pi_cols.append(column)

            # Now phi:
            try:
                column = phi_dict[n][idx]
            except KeyError:
                column = vector(base_ring, len(K.n_cells(n+1)))
            phi_cols.append(column)
        # Now iota:
        iota_cols = [iota_dict[n][c] for c in gens[n]]

        pi_data[n] = matrix(base_ring, pi_cols).transpose()
        iota_data[n] = matrix(base_ring, len(gens[n]), len(n_cells), iota_cols).transpose()
        phi_data[n] = matrix(base_ring, phi_cols).transpose()

    M = ChainComplex(M_data, base_ring=base_ring, degree=-1)
    pi = ChainComplexMorphism(pi_data, C, M)
    iota = ChainComplexMorphism(iota_data, M, C)
    phi = ChainContraction(phi_data, pi, iota)
    return phi, M

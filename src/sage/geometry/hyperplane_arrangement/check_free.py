r"""
Helper Functions For Freeness Of Hyperplane Arrangements

This contains the algorithms to check for freeness of a hyperplane
arrangement. See
:meth:`sage.geometry.hyperplane_arrangement.HyperplaneArrangementElement.is_free`
for details.

.. NOTE::

    This could be extended to a freeness check for more general modules
    over a polynomial ring.
"""
#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix
import sage.libs.singular.function_factory

def less_generators(X):
    """
    Reduce the generator matrix of the module defined by ``X``.

    This is Algorithm 6.4 in [BC12]_ and relies on the row syzygies of
    the matrix ``X``.

    EXAMPLES::

        sage: from sage.geometry.hyperplane_arrangement.check_free import less_generators
        sage: R.<x,y,z> = QQ[]
        sage: m = matrix([[1, 0, 0], [0, z, -1], [0, 0, 0], [0, y, 1]])
        sage: less_generators(m)
        [ 1  0  0]
        [ 0  z -1]
        [ 0  y  1]
    """
    R = X.base_ring()
    syz = sage.libs.singular.function_factory.ff.syz
    zero = R.zero()
    while True:
        Z = matrix(R, syz(X.transpose()))
        # get_column_independent_unit_positions
        J = range(Z.ncols())
        K = []
        for r in Z.rows():
            for j in J:
                elt = r[j]
                if elt.is_unit():
                    K.append(j)
                    J = [l for l in J if r[l] == zero]
                    break
        if not K: # K is empty
            return X
        Kd = set(range(X.nrows())).difference(K)
        X = X.matrix_from_rows(sorted(Kd))

def construct_free_chain(S, A):
    """
    Construct the free chain for the hyperplanes ``A`` in ``S``.

    ALGORITHM:

    We follow Algorithm 6.5 in [BC12]_.

    INPUT:

    - ``S`` -- set of all hyperplane arrangements
    - ``A`` -- a list of hyperplanes in ``H``

    EXAMPLES::

        sage: from sage.geometry.hyperplane_arrangement.check_free import construct_free_chain
        sage: H.<x,y,z> = HyperplaneArrangements(QQ)
        sage: A = H(z, y+z, x+y+z)
        sage: S = H.ambient_space().symmetric_space()
        sage: construct_free_chain(S, list(A))
        [
        [1 0 0]  [ 1  0  0]  [    0     1     0]
        [0 1 0]  [ 0  z -1]  [y + z     0    -1]
        [0 0 z], [ 0  y  1], [    x     0     1]
        ]
    """
    X = []
    if not A: # Empty arrangement
        return X

    syz = sage.libs.singular.function_factory.ff.syz
    # Compute the morphisms \phi_{H_j} : S^{1xR} \to S / <b_j>
    G = S.gens()
    r = len(G)
    B = [H.to_symmetric_space() for H in A] # As elements in S
    phi = [matrix(S, [[beta.derivative(x)] for x in S.gens()]) for beta in B]
    # Compute the relative row syzygy of phi_1
    mu = phi[0].stack(matrix.diagonal([B[0]]).dense_matrix())
    row_syzygy = matrix(S, syz(mu.transpose())).matrix_from_columns(range(r))
    X.append( less_generators(row_syzygy) )
    if not X[-1].is_square():
        return None
    T = matrix.identity(S, r)
    for i,f in enumerate(phi[1:]):
        T = X[-1] * T
        mu = T * f
        mu = mu.stack(matrix.diagonal([B[i+1]]).dense_matrix())
        row_syzygy = matrix(S, syz(mu.transpose())).matrix_from_columns(range(r))
        X.append( less_generators(row_syzygy) )
        if not X[-1].is_square():
            return None
    return X


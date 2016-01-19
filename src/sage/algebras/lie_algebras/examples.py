"""
Examples of Lie Algebras

There are the following examples of Lie algebras:

- A rather comprehensive family of 3-dimensional Lie
  algebras
- The Lie algebra of affine transformations of the line
- All abelian Lie algebras on free modules
- The Lie algebra of upper triangular matrices
- The Lie algebra of strictly upper triangular matrices

See also
:class:`sage.algebras.lie_algebras.virasoro.LieAlgebraRegularVectorFields`
and
:class:`sage.algebras.lie_algebras.virasoro.VirasoroAlgebra` for
other examples.

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

#from sage.algebras.lie_algebras.classical_lie_algebra import gl, sl, so, sp
from sage.algebras.lie_algebras.virasoro import VirasoroAlgebra # this is used, just not in this file

def three_dimensional(R, a, b, c, d, names=['X', 'Y', 'Z']):
    r"""
    The 3-dimensional Lie algebra over a given commutative ring `R`
    with basis `\{X, Y, Z\}` subject to the relations:

    .. MATH::

        [X, Y] = aZ + dY, \quad [Y, Z] = bX, \quad [Z, X] = cY + dZ

    where `a,b,c,d \in R`.

    This is always a well-defined 3-dimensional Lie algebra, as can
    be easily proven by computation.

    EXAMPLES::

        sage: L = lie_algebras.three_dimensional(QQ, 4, 1, -1, 2)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): 2*Y + 4*Z, ('X', 'Z'): Y - 2*Z, ('Y', 'Z'): X}
        sage: TestSuite(L).run()
        sage: L = lie_algebras.three_dimensional(QQ, 1, 0, 0, 0)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Z}
        sage: L = lie_algebras.three_dimensional(QQ, 0, 0, -1, -1)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): -Y, ('X', 'Z'): Y + Z}
        sage: L = lie_algebras.three_dimensional(QQ, 0, 1, 0, 0)
        sage: L.structure_coefficients()
        Finite family {('Y', 'Z'): X}
        sage: lie_algebras.three_dimensional(QQ, 0, 0, 0, 0)
        Abelian Lie algebra on 3 generators (X, Y, Z) over Rational Field
        sage: Q.<a,b,c,d> = PolynomialRing(QQ)
        sage: L = lie_algebras.three_dimensional(Q, a, b, c, d)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): d*Y + a*Z, ('X', 'Z'): (-c)*Y + (-d)*Z, ('Y', 'Z'): b*X}
        sage: TestSuite(L).run()
    """
    if isinstance(names, str):
        names = names.split(',')
    X = names[0]
    Y = names[1]
    Z = names[2]
    from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
    s_coeff = {(X,Y): {Z:a, Y:d}, (Y,Z): {X:b}, (Z,X): {Y:c, Z:d}}
    return LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))

def cross_product(R, names=['X', 'Y', 'Z']):
    r"""
    The Lie algebra of `\RR^3` defined by the usual cross product
    `\times`.

    EXAMPLES::

        sage: L = lie_algebras.cross_product(QQ)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Z, ('X', 'Z'): -Y, ('Y', 'Z'): X}
        sage: TestSuite(L).run()
    """
    L = three_dimensional(R, 1, 1, 1, 0, names=names)
    L.rename("Lie algebra of RR^3 under cross product over {}".format(R))
    return L

def three_dimensional_by_rank(R, n, a=None, names=['X', 'Y', 'Z']):
    """
    Return a 3-dimensional Lie algebra of rank ``n``, where `0 \leq n \leq 3`.

    Here, the *rank* of a Lie algebra `L` is defined as the dimension
    of its derived subalgebra `[L, L]`. (We are assuming that `R` is
    a field of characteristic `0`; otherwise the Lie algebras
    constructed by this function are still well-defined but no longer
    might have the correct ranks.) This is not to be confused with
    the other standard definition of a rank (namely, as the
    dimension of a Cartan subalgebra, when `L` is semisimple).

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the rank
    - ``a`` -- the deformation parameter (used for `n = 2`); this should
      be a nonzero element of `R` in order for the resulting Lie
      algebra to actually have the right rank(?)
    - ``names`` -- (optional) the generator names

    EXAMPLES::

        sage: lie_algebras.three_dimensional_by_rank(QQ, 0)
        Abelian Lie algebra on 3 generators (X, Y, Z) over Rational Field
        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
        sage: L.structure_coefficients()
        Finite family {('Y', 'Z'): X}
        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 2, 4)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Y, ('X', 'Z'): Y + Z}
        sage: L = lie_algebras.three_dimensional_by_rank(QQ, 2, 0)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Y}
        sage: lie_algebras.three_dimensional_by_rank(QQ, 3)
        sl2 over Rational Field
    """
    if isinstance(names, str):
        names = names.split(',')
    names = tuple(names)

    if n == 0:
        from sage.algebras.lie_algebras.structure_coefficients import AbelianLieAlgebra
        return AbelianLieAlgebra(R, names=names)

    if n == 1:
        L = three_dimensional(R, 0, 1, 0, 0, names=names) # Strictly upper triangular matrices
        L.rename("Lie algebra of 3x3 strictly upper triangular matrices over {}".format(R))
        return L

    if n == 2:
        if a is None:
            raise ValueError("The parameter 'a' must be specified")
        X = names[0]
        Y = names[1]
        Z = names[2]
        from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
        if a == 0:
            s_coeff = {(X,Y): {Y:R.one()}, (X,Z): {Y:R(a)}}
            # Why use R(a) here if R == 0 ? Also this has rank 1.
            L = LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))
            L.rename("Degenerate Lie algebra of dimension 3 and rank 2 over {}".format(R))
        else:
            s_coeff = {(X,Y): {Y:R.one()}, (X,Z): {Y:R.one(), Z:R.one()}}
            # a doesn't appear here :/
            L = LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))
            L.rename("Lie algebra of dimension 3 and rank 2 with parameter {} over {}".format(a, R))
        return L

    if n == 3:
        #return sl(R, 2)
        from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
        E = names[0]
        F = names[1]
        H = names[2]
        s_coeff = { (E,F): {H:R.one()}, (H,E): {E:R(2)}, (H,F): {F:R(-2)} }
        L = LieAlgebraWithStructureCoefficients(R, s_coeff, tuple(names))
        L.rename("sl2 over {}".format(R))
        return L

    raise ValueError("Invalid rank")

# This can probably be replaced (removed) by sl once the classical Lie
#   algebras are implemented
def sl(R, n, representation='bracket'):
    r"""
    Return the Lie algebra `\mathfrak{sl}_n`.

    EXAMPLES::

        sage: sl2 = lie_algebras.sl(QQ, 2); sl2
        sl2 over Rational Field
        sage: E,F,H = sl2.gens()
        sage: E.bracket(F) == H
        True
        sage: H.bracket(E) == 2*E
        True
        sage: H.bracket(F) == -2*F
        True

    TESTS::

        sage: sl2 = lie_algebras.sl(QQ, 2, representation='matrix')
        sage: E,F,H = sl2.gens()
        sage: E.bracket(F) == H
        True
        sage: H.bracket(E) == 2*E
        True
        sage: H.bracket(F) == -2*F
        True
    """
    if n != 2:
        raise NotImplementedError("only n=2 is implemented")

    if representation == 'matrix':
        from sage.matrix.matrix_space import MatrixSpace
        from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative
        MS = MatrixSpace(R, 2)
        E = MS([[0,1],[0,0]])
        F = MS([[0,0],[1,0]])
        H = MS([[1,0],[0,-1]])
        L = LieAlgebraFromAssociative(MS, [E, F, H], ['E', 'F', 'H'])
        L.rename("sl2 as a matrix Lie algebra over {}".format(R))
    elif representation == 'bracket':
        L = three_dimensional_by_rank(R, 3, names=['E', 'F', 'H'])
    else:
        raise ValueError("invalid representation")

    return L

def affine_transformations_line(R, names=['X', 'Y'], representation='bracket'):
    """
    The Lie algebra of affine transformations of the line.

    EXAMPLES::

        sage: L = lie_algebras.affine_transformations_line(QQ)
        sage: L.structure_coefficients()
        Finite family {('X', 'Y'): Y}
        sage: X, Y = L.lie_algebra_generators()
        sage: L[X, Y] == Y
        True
        sage: TestSuite(L).run()
        sage: L = lie_algebras.affine_transformations_line(QQ, representation="matrix")
        sage: X, Y = L.lie_algebra_generators()
        sage: L[X, Y] == Y
        True
        sage: TestSuite(L).run()
    """
    if isinstance(names, str):
        names = names.split(',')
    names = tuple(names)
    if representation == 'matrix':
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(R, 2, sparse=True)
        one = R.one()
        gens = tuple(MS({(0,i):one}) for i in range(2))
        from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative
        return LieAlgebraFromAssociative(MS, gens, names=names)
    X = names[0]
    Y = names[1]
    from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
    s_coeff = {(X,Y): {Y:R.one()}}
    L = LieAlgebraWithStructureCoefficients(R, s_coeff, names=names)
    L.rename("Lie algebra of affine transformations of a line over {}".format(R))
    return L

def abelian(R, names):
    """
    Return the abelian Lie algebra generated by ``names``.

    EXAMPLES::

        sage: lie_algebras.abelian(QQ, 'x, y, z')
        Abelian Lie algebra on 3 generators (x, y, z) over Rational Field
    """
    if isinstance(names, str):
        names = names.split(',')
    from sage.algebras.lie_algebras.structure_coefficients import AbelianLieAlgebra
    return AbelianLieAlgebra(R, tuple(names))

def Heisenberg(R, n, representation="structure"):
    """
    Return the rank ``n`` Heisenberg algebra in the given representation.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the rank (a nonnegative integer or infinity)
    - ``representation`` -- (default: "structure") can be one of the following:

      - ``"structure"`` -- using structure coefficients
      - ``"matrix"`` -- using matrices

    EXAMPLES::

        sage: lie_algebras.Heisenberg(QQ, 3)
        Heisenberg algebra of rank 3 over Rational Field
    """
    from sage.rings.infinity import infinity
    if n == infinity:
        from sage.algebras.lie_algebras.heisenberg import InfiniteHeisenbergAlgebra
        return InfiniteHeisenbergAlgebra(R)
    if representation == "matrix":
        from sage.algebras.lie_algebras.heisenberg import HeisenbergAlgebra_matrix
        return HeisenbergAlgebra_matrix(R, n)
    from sage.algebras.lie_algebras.heisenberg import HeisenbergAlgebra
    return HeisenbergAlgebra(R, n)

def regular_vector_fields(R):
    r"""
    Return the Lie algebra of regular vector fields on `\CC^{\times}`.

    .. SEEALSO::

        :class:`~sage.algebras.lie_algebras.virasoro.LieAlgebraRegularVectorFields`

    EXAMPLES::

        sage: lie_algebras.regular_vector_fields(QQ)
        The Lie algebra of regular vector fields over Rational Field
    """
    from sage.algebras.lie_algebras.virasoro import LieAlgebraRegularVectorFields
    return LieAlgebraRegularVectorFields(R)

def pwitt(R, p):
    r"""
    Return the `p`-Witt Lie algebra over `R`.

    EXAMPLES::

        sage: lie_algebras.pwitt(GF(5), 5)
        The 5-Witt Lie algebra over Finite Field of size 5
    """
    from sage.algebras.lie_algebras.virasoro import WittLieAlgebra_charp
    return WittLieAlgebra_charp(R, p)

def upper_triangular_matrices(R, n):
    r"""
    Return the Lie algebra `\mathfrak{b}_k` of `k \times k` upper
    triangular matrices.

    .. TODO::

        This implementation does not know it is finite-dimensional and
        does not know its basis.

    EXAMPLES::

        sage: L = lie_algebras.upper_triangular_matrices(QQ, 4); L
        Lie algebra of 4-dimensional upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: n0, n1, n2, t0, t1, t2, t3 = L.lie_algebra_generators()
        sage: L[n2, t2] == -n2
        True

    TESTS::

        sage: L = lie_algebras.upper_triangular_matrices(QQ, 1); L
        Lie algebra of 1-dimensional upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: L = lie_algebras.upper_triangular_matrices(QQ, 0); L
        Lie algebra of 0-dimensional upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
    """
    from sage.matrix.matrix_space import MatrixSpace
    from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative
    MS = MatrixSpace(R, n, sparse=True)
    one = R.one()
    names = tuple('n{}'.format(i) for i in range(n-1))
    names += tuple('t{}'.format(i) for i in range(n))
    gens = [MS({(i,i+1):one}) for i in range(n-1)]
    gens += [MS({(i,i):one}) for i in range(n)]
    L = LieAlgebraFromAssociative(MS, gens, names=names)
    L.rename("Lie algebra of {}-dimensional upper triangular matrices over {}".format(n, L.base_ring()))
    return L

def strictly_upper_triangular_matrices(R, n):
    r"""
    Return the Lie algebra `\mathfrak{n}_k` of strictly `k \times k` upper
    triangular matrices.

    .. TODO::

        This implementation does not know it is finite-dimensional and
        does not know its basis.

    EXAMPLES::

        sage: L = lie_algebras.strictly_upper_triangular_matrices(QQ, 4); L
        Lie algebra of 4-dimensional strictly upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: n0, n1, n2 = L.lie_algebra_generators()
        sage: L[n2, n1]
        [ 0  0  0  0]
        [ 0  0  0 -1]
        [ 0  0  0  0]
        [ 0  0  0  0]

    TESTS::

        sage: L = lie_algebras.strictly_upper_triangular_matrices(QQ, 1); L
        Lie algebra of 1-dimensional strictly upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
        sage: L = lie_algebras.strictly_upper_triangular_matrices(QQ, 0); L
        Lie algebra of 0-dimensional strictly upper triangular matrices over Rational Field
        sage: TestSuite(L).run()
    """
    from sage.matrix.matrix_space import MatrixSpace
    from sage.algebras.lie_algebras.lie_algebra import LieAlgebraFromAssociative
    MS = MatrixSpace(R, n, sparse=True)
    one = R.one()
    names = tuple('n{}'.format(i) for i in range(n-1))
    gens = tuple(MS({(i,i+1):one}) for i in range(n-1))
    L = LieAlgebraFromAssociative(MS, gens, names=names)
    L.rename("Lie algebra of {}-dimensional strictly upper triangular matrices over {}".format(n, L.base_ring()))
    return L


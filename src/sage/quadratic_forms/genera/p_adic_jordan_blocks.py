r"""
Algorithms for finite quadratic forms and their relation to genera.

AUTHORS:

- Simon Brandhorst (2017-10-09): initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 Simon Branhdorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import Zp, ZZ, GF
from sage.matrix.constructor import Matrix
from copy import copy
from sage.rings.finite_rings.integer_mod import mod

def p_adic_normal_form(G,p,precision=None,normalize=True):
    r"""
    Return the `p`-adic Jordan form of a symmetric matrix.

    Let `p` be odd and `u` be the smallest non-square modulo `p`.
    The Jordan form is a diagonal matrix with diagonal entries either `p^n` or `u*p^n`.

    If `p=2` is even, then the Jordan form consists of 
    1 x 1 blocks of the form
    `0`, `[2^n]`, `[3*2^n]`, `[5*2^n]`, `[7*2^n]`
    or of 2 x 2 blocks of the form
    `[[2, 1],[1, 2]] * 2^n`
    or  
    `[[0, 1],[1,0]] * 2^n`
    In any case the entries are ordered by their valuation.

    INPUT:

        - ``G`` -- a symmetric n x n matrix in `\QQ`
        - ``p`` -- a prime number -- it is not checked whether ``p`` is prime
        - ``precision`` -- defining precision. If not set, 
          the minimal possible is taken.
        - ``normalize``-- bool (default: ``True``)

    OUTPUT:

        - ``D`` -- the jordan matrix over `\QQ_p`
        - ``U`` -- invertible transformation matrix over `\ZZ_p`, 
          i.e, ``D = U * G * U^T``

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import p_adic_normal_form
        sage: D4 = Matrix(ZZ,4,[2,-1,-1,-1,-1,2,0,0,-1,0,2,0,-1,0,0,2])
        sage: D4
        [ 2 -1 -1 -1]
        [-1  2  0  0]
        [-1  0  2  0]
        [-1  0  0  2]
        sage: D, U = p_adic_normal_form(D4,2)
        sage: D            
        [  2   1   0   0]
        [  1   2   0   0]
        [  0   0 2^2   2]
        [  0   0   2 2^2]
        sage: D == U * D4 * U.T
        True
        sage: A4 = Matrix(ZZ,4,[2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2])
        sage: A4
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
        sage: D, U = p_adic_normal_form(A4,2)
        sage: D        
        [2 1 0 0]
        [1 2 0 0]
        [0 0 0 1]
        [0 0 1 0]

    We can handle degenerate forms as well::

        sage: A4_extended = Matrix(ZZ,5,[2, -1, 0, 0, -1, -1, 2, -1, 0, 0, 0, -1, 2, -1, 0, 0, 0, -1, 2, -1, -1, 0, 0, -1, 2])
        sage: D, U = p_adic_normal_form(A4_extended,5)
        sage: D
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 5 0]
            [0 0 0 0 0]

    Rational matrices as no problem as well::

        sage: A4dual = A4.inverse()
        sage: D, U = p_adic_normal_form(A4dual,5)
        sage: D
        [5^-1    0    0    0]
        [   0    1    0    0]
        [   0    0    1    0]
        [   0    0    0    1]

    TESTS::

        sage: Z = Matrix(ZZ,0,[])
        sage: p_adic_normal_form(Z,3)
        ([], [])
        sage: Z = matrix.zero(10)
        sage: p_adic_normal_form(Z,3)[0] == 0
        True
    """
    # this calculation is local. Thus we can do it over the p-adics
    p = ZZ(p)
    denom = G.denominator()
    G = G * denom
    if precision == None:
        if G.det() != 0:
            # we have to calculate at least mod 8 for 2-adic things to make sense.
            precision = G.det().valuation(p) + 3 
        else:
            # the non-degenerate part
            B = G.row_space().basis_matrix()
            GG = B * G * B.T
            # we have to calculate at least mod 8 for things to make sense.
            precision = GG.det().valuation(p) + 3
    R = Zp(p, prec = precision, type = 'fixed-mod')
    G = G.change_ring(R) # is not changed
    D = copy(G)    # is transformed into jordan form
    n = G.ncols()
    # The trivial case
    if n == 0:
        return G.parent().zero(), G.parent().zero()
    # transformation matrix
    U = Matrix.identity(R,n)

    if(p == 2):
        D, U = _jordan_2_adic(G)
        # we confirm the result
        assert U*G*U.T == Matrix.block_diagonal(_collect_small_blocks(D))
    else:
        D, U = _jordan_odd_adic(G)
        # we confirm the result
        assert U*G*U.T == Matrix.diagonal(D.diagonal())
    assert U.determinant().valuation() == 0     # U is invertible!
    if normalize:
        D, U1 = _normalize(D)
        U = U1 * U
    assert U*G*U.T == Matrix.block_diagonal(_collect_small_blocks(D))
    return D/denom, U

def _jordan_odd_adic(G):
    r"""
    Return the Jordan decomposition over an odd prime.

    INPUT:

        -- a symmetric matrix over ``Zp`` of type ``'fixed-mod'``

    OUTPUT:

        - ``D`` -- a diagonal matrix
        - ``U`` -- a unimodular matrix such that ``D = U*G*U.T``

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _jordan_odd_adic
        sage: R = Zp(3,prec=2,print_mode='series')
        sage: A4 = Matrix(R,4,[2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2])
        sage: A4
        [      2 + O(3^2) 2 + 2*3 + O(3^2)                0                0]
        [2 + 2*3 + O(3^2)       2 + O(3^2) 2 + 2*3 + O(3^2)                0]
        [               0 2 + 2*3 + O(3^2)       2 + O(3^2) 2 + 2*3 + O(3^2)]
        [               0                0 2 + 2*3 + O(3^2)       2 + O(3^2)]
        sage: D, U = _jordan_odd_adic(A4)
        sage: D.change_ring(ZZ)   # just for pretty printing.
        [2 0 0 0]
        [0 2 0 0]
        [0 0 1 0]
        [0 0 0 8]
        sage: D == U*A4*U.T
        True
        sage: U.determinant().valuation() == 0
        True
    """
    R = G.base_ring()
    D = copy(G)
    n = G.ncols()

    # transformation matrix
    U = Matrix.identity(R, n)

    # indices of the diagonal entrys which are already used
    cnt = 0
    while cnt < n:
        pivot = _find_min_p(D, cnt)
        piv1 = pivot[1]
        piv2 = pivot[2]
        minval = pivot[0]
        # the smallest valuation is on the diagonal
        if piv1 == piv2:
            # move pivot to position [cnt,cnt]
            if piv1 != cnt:
                U.swap_rows(cnt, piv1)
                D.swap_rows(cnt, piv1)
                D.swap_columns(cnt, piv1)
            # we are already orthogonal to the part with i < cnt
            # now make the rest orthogonal too
            for i in range(cnt+1,n):
                if D[i, cnt]!= 0:
                    c = D[i, cnt] // D[cnt, cnt]
                    U[i, :] += - c * U[cnt, :]
                    D[i, :] += - c * D[cnt, :]
                    D[:, i] += - c * D[:, cnt]
            cnt = cnt + 1
        else:
            # the smallest valuation is off the diagonal
            row = pivot[1]
            col = pivot[2]
            U[row, :] += U[col, :]
            D[row, :] += D[col, :]
            D[:, row] += D[:, col]
            # the smallest valuation is now on the diagonal
    return D, U

def _jordan_2_adic(G):
    r"""
    Transform a symmetric matrix over the `2`-adic integers into jordan form

    Note that if the precision is too low this method fails.
    The method is only tested for input over `\ZZ_2` of`'type=fixed-mod'`.

    INPUT:

        - ``G`` -- symmetric `n` x `n` matrix in ``Zp``

    OUTPUT:

        - ``D`` -- the jordand matrix
        - ``U`` -- transformation matrix i.e ``D = U * G * U^T``

    The matrix ``D`` is a block diagonal matrix consisting 
    of `1` x `1` and `2` x `2` blocks.
    The 2 x 2 blocks are of the form
    `[2a  b]`
    `[ b 2c] * 2^k`
    with `b` of valuation `0`.

    EXAMPLES:

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _jordan_2_adic
        sage: R = Zp(2,prec=3,print_mode='series')
        sage: A4 = Matrix(R,4,[2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2])
        sage: A4.change_ring(ZZ)   # for pretty printing
        [2 7 0 0]
        [7 2 7 0]
        [0 7 2 7]
        [0 0 7 2]

        sage: D, U = _jordan_2_adic(A4)
        sage: D.change_ring(ZZ)   # for pretty printing.
        [ 2  7  0  0]
        [ 7  2  0  0]
        [ 0  0 12  7]
        [ 0  0  7  2]

        sage: D == U*A4*U.T
        True
        sage: U.determinant().valuation() == 0
        True
    """
    R = G.base_ring()
    D = copy(G)
    n = G.ncols()

    # transformation matrix
    U = Matrix.identity(R, n)

    # indices of the diagonal entrys which are already used
    cnt = 0
    minval = None
    while cnt < n:
            pivot = _find_min_p(D, cnt)
            piv1 = pivot[1]
            piv2 = pivot[2]
            minval_last = minval
            minval = pivot[0]
            # the smallest valuation is on the diagonal
            if piv1 == piv2:
                # move pivot to position [cnt,cnt]
                if piv1 != cnt:
                    U.swap_rows(cnt, piv1)
                    D.swap_rows(cnt, piv1)
                    D.swap_columns(cnt, piv1)
                # we are already orthogonal to the part with i < cnt
                # now make the rest orthogonal too
                for i in range(cnt+1, n):
                    if D[i, cnt] != 0:
                        c = D[i, cnt]//D[cnt, cnt]
                        U[i, :] += -c * U[cnt, :]
                        D[i, :] += -c * D[cnt, :]
                        D[:, i] += -c * D[:, cnt]
                cnt = cnt + 1
            # the smallest valuation is off the diagonal
            else:
                # move this 2 x 2 block to the top left (starting from cnt)
                if piv1 != cnt:
                    U.swap_rows(cnt, piv1)
                    D.swap_rows(cnt, piv1)
                    D.swap_columns(cnt, piv1)
                if piv2 != cnt+1:
                    U.swap_rows(cnt+1, piv2)
                    D.swap_rows(cnt+1, piv2)
                    D.swap_columns(cnt+1, piv2)
                # we split off a 2 x 2 block
                # if it is the last 2 x 2 block, there is nothing to do.
                if cnt != n-2:
                    content = R(2 ** minval)
                    eqn_mat = D[cnt:cnt+2, cnt:cnt+2].list()
                    eqn_mat = Matrix(R, 2, 2, [e // content for e in eqn_mat]) 
                    # calculate the inverse without using division
                    inv = eqn_mat.adjoint() * eqn_mat.det().inverse_of_unit()
                    B = U[cnt:cnt+2, :]
                    C = D[cnt+2:, cnt:cnt+2] * inv
                    for i in range(C.nrows()):
                        for j in range(C.ncols()):
                            C[i, j]=C[i, j] // content 
                    U[cnt+2:, :] -= C * B
                    D[cnt:, cnt:] = U[cnt:, :] * G * U[cnt:, :].transpose()
                # if the previous diagonal entry has the same valuation
                # as the pivot we can split off a 1 x 1 matrix
                #if cnt > 0 and D[cnt-1, cnt-1].valuation() == minval == minval_last:
                #    # a 3 x 3 block of D
                #    # a  0  0
                #    # 0 2b  c
                #    # 0  c 2d
                #    a = D[cnt-1, cnt-1].unit_part()
                #    c = D[cnt, cnt+1].unit_part()
                #    B = Matrix(R, 3, 3, [
                #                    1, 1, 0,
                #                    c, 0,-a,
                #                    0, 0, 1
                #                    ])   # this matrix is invertible
                #    U[cnt-1:cnt+2, :] = B * U[cnt-1:cnt+2, :]
                #    D = U * G * U.transpose()
                #    # at this stage cnt-1 and cnt are orthogonal and both have minimal valuation
                #    cnt -= 3
                cnt += 2
    return D, U

def _get_small_block_indices(G):
    r"""
    Return the indices of the blocks.

    For internal use in :meth:`_collect_small_blocks`

    INPUT:

        - ``G`` -- a block_diagonal matrix consisting of `1 x 1` and `2 x 2` blocks

    OUTPUT:

        - a list of integers

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _get_small_block_indices
        sage: W1 = Matrix([1])
        sage: U = Matrix(ZZ, 2, [0, 1, 1, 0])
        sage: G = Matrix.block_diagonal([W1, U, U, W1, W1, U, W1, U])
        sage: _get_small_block_indices(G)
        [0, 1, 3, 5, 6, 7, 9, 10]
    """
    L = []
    n = G.ncols()
    i = 0
    while i < n-1:
        L.append(i)
        if G[i, i+1]!=0:
            i += 2
        else:
            i += 1
    if i == n-1:
        L.append(i)
    return L[:]

def _collect_small_blocks(G):
    r"""
    Return the blocks as list.

    INPUT:

        - ``G`` -- a block_diagonal matrix consisting of `1` x `1` and `2` x `2` blocks

    OUTPUT:

        - a list of `1` x `1` and `2` x `2` matrices -- the blocks

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _collect_small_blocks
        sage: W1 = Matrix([1])
        sage: V = Matrix(ZZ, 2, [2, 1, 1, 2])
        sage: L = [W1, V, V, W1, W1, V, W1, V]
        sage: G = Matrix.block_diagonal(L)
        sage: L == _collect_small_blocks(G)
        True
    """
    D = copy(G)
    L = _get_small_block_indices(D)[1:]
    D.subdivide(L, L)
    blocks = []
    for i in range(len(L)+1):
        block = copy(D.subdivision(i, i))
        blocks.append(block)
    return blocks

def _get_homogeneous_block_indices(G):
    r"""
    Return the indices of the homogeneous blocks.
    
    EXAMPLES::
    
        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _get_homogeneous_block_indices
        sage: W1 = Matrix(Zp(2),[1])
        sage: V = Matrix(Zp(2), 2, [2, 1, 1, 2])
        sage: G = Matrix.block_diagonal([W1, V, V, 2*W1, 2*W1, 8*V, 8*W1, 16*V])
        sage: _get_homogeneous_block_indices(G)
        [0, 5, 7, 10]
        sage: G = Matrix.block_diagonal([W1, V, V, 2*W1, 2*W1, 8*V, 8*W1, 16*W1])
        sage: _get_homogeneous_block_indices(G)
        [0, 5, 7]
    """
    L = []
    n = G.ncols()
    i = 0
    val = -5
    while i < n-1:
        m = G[i,i].valuation()
        m = min(m,G[i,i+1].valuation())
        if m > val:
            L.append(i)
            val = m
        if G[i, i+1]!=0:
            i += 1
        i += 1
    if i == n-1:
        if m > val:
            L.append(i)
            val = m
    return L
    
def _normalize(G):
    r"""
    Return the partial normal form of ``G``.

    See also :meth:``p_adic_jordan_blocks``.

    INPUT:

        - a symmetric matrix over `Zp` in jordan form --
        the OUTPUT of :meth:``p_adic_normal_form`` or :meth:``_jordan_2_adic``

    OUTPUT:

        - ``(D, U)`` - a pair of matrices such that ``D=U*G*U.T`` is in partial normal form

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import         _normalize
        sage: R = Zp(3, prec = 5, type = 'fixed-mod', print_mode = 'digits')
        sage: G = matrix.diagonal(R,[1,7,3,3*5,3,9,-9,27*13])
        sage: D,U =_normalize(G)
        sage: D
        [...00001        0        0        0        0        0        0        0]
        [       0 ...00001        0        0        0        0        0        0]
        [       0        0 ...00010        0        0        0        0        0]
        [       0        0        0 ...00010        0        0        0        0]
        [       0        0        0        0 ...00020        0        0        0]
        [       0        0        0        0        0 ...00100        0        0]
        [       0        0        0        0        0        0 ...00200        0]
        [       0        0        0        0        0        0        0 ...01000]
    """
    R = G.base_ring()
    D = copy(G)
    p = R.prime()
    n = G.ncols()
    U = copy(G.parent().identity_matrix())
    if p != 2:
        # squareclasses 1, v
        v = _min_nonsquare(p)
        v = R(v)
        non_squares = []
        val = 0
        for i in range(n):
            if D[i, i] != 0:
                if D[i,i].valuation()>val:
                    # a new block starts
                    val = D[i,i].valuation()
                    if len(non_squares)!=0:
                        # move the non-square to
                        # the last entry of the last block
                        j=non_squares.pop()
                        U.swap_rows(j,i-1)
                d = D[i, i].unit_part()
                if d.is_square():
                    D[i, i] = 1
                    U[i, :] *= d.inverse_of_unit().sqrt()
                else:
                    D[i, i] = v
                    U[i, :] *= (v * d.inverse_of_unit()).sqrt()
                    if len(non_squares)!=0:
                        j = non_squares.pop()
                        trafo = _normalize_odd_2x2(D[[i,j],[i,j]])
                        U[[i,j],:] = trafo*U[[i,j],:]
                    else:
                        non_squares.append(i)
        if len(non_squares)!=0:
            j=non_squares.pop()
            U.swap_rows(j,n-1)
    else:
        # find a homogeneous block and treat it
        # squareclasses 1,3,5,7 modulo 8
        for i in range(n):
            d = D[i, i].unit_part()
            if d != 0:
                v = R(mod(d, 8))
                U[i, :] *= (v * d.inverse_of_unit()).sqrt()
        D = U * G * U.T
        for i in range(n-1):
            b = D[i, i+1]
            if b !=0:    # there is a 2 x 2 block here
                block = D[i:i+2, i:i+2]
                trafo = _normalize_2x2(block)
                U[i:i+2, :] = trafo * U[i:i+2, :]
    D = U * G * U.T
    return D, U

def _normalize_2x2(G):
    r"""
    normalize this 2x2 block

    INPUT:

        -  ``G`` - a `2` by `2` matrix over ``Zp`` 
            with ``type = 'fixed-mod'`` of the form
            `[2a  b]`
            `[ b 2c]*2^n`
            with b of valuation 1

    OUTPUT:

        -  a unimodular `2 x 2` matrix ``U`` over ``Zp`` with 
            ``U * G *U.transpose()``
            either 
                `[0 1]`
                `[1 0]`
            or
                `[2 1]`
                `[1 2]`

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import         _normalize_2x2
        sage: R = Zp(2, prec = 10, type = 'fixed-mod', print_mode = 'digits')
        sage: G = Matrix(R,2,[-17*2,3,3,23*2])
        sage: U =_normalize_2x2(G)
        sage: U*G*U.T
        [...0000000010 ...0000000001]
        [...0000000001 ...0000000010]

        sage: G = Matrix(R,2,[-17*4,3,3,23*2])
        sage: U = _normalize_2x2(G)
        sage: U*G*U.T
        [            0 ...0000000001]
        [...0000000001             0]

        sage: G = Matrix(R,2,[-17*2,3,3,23*2])
        sage: U = _normalize_2x2(8*G)
        sage: U*G*U.T
        [...1110000010 ...1010000001]
        [...1010000001 ...1110000010]

    TESTS::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _normalize_2x2
        sage: R = Zp(2, prec = 10, type = 'fixed-mod')
        sage: ref1 = Matrix(R,2,[2,1,1,2])
        sage: ref2 = Matrix(R,2,[0,1,1,0])
        sage: N = _normalize_2x2(G)
        sage: (N*G*N.T == ref1) or (N*G*N.T == ref2)
        True
    """
    from sage.rings.all import PolynomialRing
    from sage.modules.free_module_element import vector
    U = copy(G.parent().identity_matrix())
    R = G.base_ring()
    P = PolynomialRing(R, 'x')
    x = P.gen()

    # The input must be an even block
    odd1 = (G[0, 0].valuation() < G[1, 0].valuation())
    odd2 = (G[1, 1].valuation() < G[1, 0].valuation())
    if  odd1 or odd2:
            raise ValueError("Not a valid 2 x 2 block.")
    scale = 2 ** G[0,1].valuation()
    D = Matrix(R, 2, 2, [d // scale for d in G.list()])
    # now D is of the form
    # [2a b ]
    # [b  2c]
    # where b has valuation 1.
    G = copy(D)

    # Make sure G[1, 1] has valuation 1.
    if D[1, 1].valuation() > D[0, 0].valuation():
        U.swap_columns(0, 1)
        D.swap_columns(0, 1)
        D.swap_rows(0, 1)
    if D[1, 1].valuation() != 1:
        # this works because
        # D[0, 0] has valuation at least 2
        U[1, :] += U[0, :]
        D = U * G * U.transpose()
    assert D[1, 1].valuation() == 1

    if mod(D.det(), 8) == 3:
        #  in this case we can transform D to
        #  2 1
        #  1 2
        # Find a point of norm 2
        # solve: 2 == D[1,1]*x^2 + 2*D[1,0]*x + D[0,0]
        pol = (D[1,1]*x**2 + 2*D[1,0]*x + D[0,0]-2) // 2
        # somehow else pari can get a hickup see `trac`:#24065
        pol = pol // pol.leading_coefficient() 
        sol = pol.roots()[0][0]
        U[0, 1] = sol
        D = U * G * U.transpose()
        # make D[0, 1] = 1
        U[1, :] *= D[1, 0].inverse_of_unit()
        D = U * G * U.transpose()

        # solve: v*D*v == 2 with v = (x,-2*x+1)
        if D[1, 1] != 2:
            v = vector([x, -2*x + 1])
            pol = (v*D*v - 2) // 2
            # somehow else pari can get a hickup `trac`:#24065
            pol = pol // pol.leading_coefficient() 
            sol = pol.roots()[0][0]
            U[1, :] = sol * U[0,:] + (-2*sol + 1)*U[1, :]
            D = U * G * U.transpose()
        # check the result
        assert D == Matrix(G.parent(), 2, 2, [2, 1, 1, 2]), "D1 \n %r" %D
    elif mod(D.det(), 8) == 7:
        # in this case we can transform D to
        #  0 1
        #  1 0
        # Find a point representing 0
        # solve: 0 == D[1,1]*x^2 + 2*D[1,0]*x + D[0,0]
        pol = (D[1,1]*x**2 + 2*D[1,0]*x + D[0,0])//2
        # somehow else pari can get a hickup, see  `trac`:#24065
        pol = pol // pol.leading_coefficient() 
        sol = pol.roots()[0][0]
        U[0,:] += sol*U[1, :]
        D = U * G * U.transpose()
        # make the second basis vector have 0 square as well.
        U[1, :] = U[1, :] - D[1, 1]//(2*D[0, 1])*U[0,:]
        D = U * G * U.transpose()
        # rescale to get D[0,1] = 1
        U[0, :] *= D[1, 0].inverse_of_unit()
        D = U * G * U.transpose()
        # check the result
        assert D == Matrix(G.parent(), 2, 2, [0, 1, 1, 0]), "D2 \n %r" %D
    return U

def _find_min_p(G, cnt):
    r"""
    Find smallest valuation below and right from ``cnt`` prefering the diagonal.

    INPUT:

        - ``G`` -- symmetric `n` x `n` matrix in `\Qp`
        - ``cnt`` -- start search from this index

    OUTPUT:

        - ``min`` -- minimal valuation
        - ``min_i`` -- row of the minimal valuation
        - ``min_j`` -- column of the minimal valuation

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _find_min_p
        sage: G = matrix(Qp(2),3,3,[4,0,1,0,4,2,1,2,1])
        sage: G    
        [2^2 + O(2^22)             0   1 + O(2^20)]
        [            0 2^2 + O(2^22)   2 + O(2^21)]
        [  1 + O(2^20)   2 + O(2^21)   1 + O(2^20)]

        sage: _find_min_p(G,0)
        (0, 2, 2)        

        sage: G = matrix(Qp(3),3,3,[4,0,1,0,4,2,1,2,1])
        sage: G
        [1 + 3 + O(3^20)               0     1 + O(3^20)]
        [              0 1 + 3 + O(3^20)     2 + O(3^20)]
        [    1 + O(3^20)     2 + O(3^20)     1 + O(3^20)]
        sage: _find_min_p(G,0)
        (0, 2, 2)
    """
    n = G.ncols()
    min = G[cnt, cnt].valuation()
    min_i = cnt
    min_j = cnt
    for i in range(cnt, n):
        for j in range(i, n):
            if (G[i, j]!=0) and (G[i, j].valuation() < min):
                min_i = i
                min_j = j
                min = G[i, j].valuation()
        # diagonal has precedence
        if (G[i, i]!=0) and (G[i, i].valuation() <= min):
            min_i = i
            min_j = i
            min = G[i, i].valuation()
    return min, min_i, min_j

def _min_nonsquare(p):
    r"""
    Return the minimal nonsquare modulo the prime ``p``.

    INPUT:

        - ``p`` -- a prime number

    OUTPUT:

        - ``a`` -- the minimal nonsquare mod ``p`` as an element of `\ZZ`.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _min_nonsquare
        sage: _min_nonsquare(2)
        sage: _min_nonsquare(3)
        2
        sage: _min_nonsquare(5)
        2
        sage: _min_nonsquare(7)
        3
    """
    R = GF(p)
    for i in R:
        if not R(i).is_square():
            return i

def _normalize_odd_2x2(G):
    r"""
    normalize this 2x2 block

    INPUT:

        - ``G`` -- a multiple of the 2x2 identity_matrix over the p-adics ``p!=2``

    OUTPUT:

        - A transformation matrix ``U`` such that
          ``U*G*U.T`` is the identiy matrix

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _normalize_odd_2x2
        sage: G = 2*Matrix.identity(Qp(5),2)
        sage: U = _normalize_odd_2x2(G)
        sage: U*G*U.T
        [1 + O(5^20)     O(5^20)]
        [    O(5^20) 1 + O(5^20)]
    """
    assert G[0,0]==G[1,1]
    u = G[0,0]
    y = G.base_ring().zero()
    while not (1/u-y**2).is_square():
        y = y + 1
    x = (1/u-y**2).sqrt()
    U = copy(G.parent().identity_matrix())
    U[0,0] = x
    U[0,1] = y
    U[1,0] = y
    U[1,1] = -x
    return U

def _relations(G,n):
    r"""
    Return relations of 2-adic quadratic forms

    INPUT:

        - ``n`` -- an integer between 1 and 10 the number of the relation
        - ``G`` -- a square matrix - the left side of the relation

    OUTPUT:

        - square matrix ``B`` such that ``B*G*B.T`` is the right side of the relation

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _relations
        sage: R = Zp(2, type = 'fixed-mod')
        sage: U = Matrix(R,2,[0,1,1,0])
        sage: V = Matrix(R,2,[2,1,1,2])
        sage: W1 = Matrix(R,1,[1])
        sage: W3 = Matrix(R,1,[3])
        sage: W5 = Matrix(R,1,[5])
        sage: W7 = Matrix(R,1,[7])
        sage: G = Matrix.block_diagonal([W1,W3])
        sage: b = _relations(G,1)
        sage: b*G*b.T
        [    1 + 2^2 + O(2^20)               O(2^20)]
        [              O(2^20) 1 + 2 + 2^2 + O(2^20)]
        
        sage: G = Matrix.block_diagonal([V,V])
        sage: b = _relations(G,3)
        sage: b*G*b.T    
        [    O(2^20) 1 + O(2^20)     O(2^20)     O(2^20)]
        [1 + O(2^20)     O(2^20)     O(2^20)     O(2^20)]
        [    O(2^20)     O(2^20)     O(2^20) 1 + O(2^20)]
        [    O(2^20)     O(2^20) 1 + O(2^20)     O(2^20)]
    """
    R = G.base_ring()
    if n == 1:
        e1 = G[0,0].unit_part()
        e2 = G[1,1].unit_part()
        B = Matrix(R,2,[1,2,2*e2,-e1])
    if n == 2:
        e1 = G[0,0].unit_part()
        e2 = G[1,1].unit_part()
        e3 = G[2,2].unit_part()
        s1 = e1 + e2 + e3
        s2 = e1*e2 + e1*e3 + e2*e3
        s3 = e1*e2*e3
        B = Matrix(R,3,[1,1,1,e2,-e1,0,e3,0,-e1])
    if n == 3:
        B = Matrix(R,4,[1,1,1,0, 1,1,0,1, 1,0,-1,-1, 0,1,-1,-1])
    if n == 4:
        pass
    if n == 5:
        e1 = G[2,2].unit_part()
        e2 = G[3,3].unit_part()
        B = Matrix(R,4,[ 1,  1,  0, 0,
                            0,  1,  1, 1, 
                           -e1, e1, 1, 0, 
                           -e2, e2, 0, 1])
    if n == 6:
        e1 = G[0,0].unit_part()
        e2 = G[1,1].unit_part()
        B = Matrix(R,2,[1,1,-2*e2,e1])
    if n == 7:
        e = G[0,0]
        B = Matrix(R,3,[ 1, 1, 1,
                        -2, e, 0,
                        -2, 0, e])
    if n == 8:
        e = G[2,2]
        B = Matrix(R,3,[ 2, 2,  1, 
                         e, 0, -1,
                        -2, 0,  e])
    if n == 9:
        B = Matrix(R,3,[ 1, 0,  1,
                        -2, 0, e1,
                         0, 1,  0])
    if n == 10:
        e1 = G[0,0]
        e2 = G[1,1]
        B = Matrix(R,2,[1,1,-4*e2,e1])
    D, U = _normalize(B*G*B.T)
    return U*B

def _partial_normal_form(G):
    r"""
    Return the partial normal form of the homogeneous block ``G``
    
    INPUT:
    
        - ``G`` -- a modular symmetric matrix over the `2`-adic integers
        
    OUTPUT:
    
        - ``B`` -- a transformation matrix such that 
          ``B*G*B.T`` is in homogeneous normal form
    
    EXAMPLES::
        
        sage: from sage.quadratic_forms.genera.p_adic_jordan_blocks import _partial_normal_form
        sage: R = Zp(2,prec=4, type = 'fixed-mod')
        sage: U = Matrix(R,2,[0,1,1,0])
        sage: V = Matrix(R,2,[2,1,1,2])
        sage: W1 = Matrix(R,1,[1])
        sage: W3 = Matrix(R,1,[3])
        sage: W5 = Matrix(R,1,[5])
        sage: W7 = Matrix(R,1,[7])
        sage: G = Matrix.block_diagonal([W1,U,V,W5,V,W3,V,W7])
        sage: B = _partial_normal_form(G)
        sage: (B*G*B.T).change_ring(ZZ)
        [0 1 0 0 0 0 0 0 0 0 0 0]
        [1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 2 1 0 0]
        [0 0 0 0 0 0 0 0 1 2 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 7]
        sage: G = Matrix.block_diagonal([W1,U,V,W1,V,W1,V,W7])
        sage: B = _partial_normal_form(G)
        sage: (B*G*B.T).change_ring(ZZ)
        [0 1 0 0 0 0 0 0 0 0 0 0]
        [1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 3 0]
        [0 0 0 0 0 0 0 0 0 0 0 7]
    """
    D = copy(G)
    R = D.base_ring()
    n = D.ncols()
    B = copy(G.parent().identity_matrix())     # the transformation matrix
    blocks = _get_small_block_indices(D)
    # collect the indices of forms of types U, V and W
    U = []
    W = []
    V = []
    for i in blocks:
        if i+1 in blocks or i==n-1:
            W.append(i)
        else:
            if D[i,i] != 0:
                V += [i,i+1]
            else:
                U += [i,i+1]
        if len(W) == 3:
            # W W W transforms to W U or W V
            T = _relations(D[W,W],2)
            B[W,:] = _relations(D[W,W],2) * B[W,:]
            D = B * G * B.T
            if mod(D[W[1:],W[1:]].det().unit_part(),8) == 3:
                V += W[1:]
            else:
                U += W[1:]
            W = W[:1]
        if len(V) == 4:
            B[V,:] = _relations(D[V,V],3) * B[V,:]
            U += V
            V = []
            D = B * G * B.T
    # put everything into the right order
    UVW = U + V + W
    B = B[UVW,:]
    D = B * G * B.T
    return B
    
def _homogeneous_normal_form(G):
    r"""
    Return the homogeneous normal form of ``G``.
    
    INPUT:
    
        - ``G`` -- a modular symmetric matrix over the `2`-adic integers
    
    OUTPUT:
    
        - ``B`` -- an invertible matrix over the basering of ``G`` 
          such that ``B*G*B.T`` is in homogeneous normal form
    
    EXAMPLES::
    
        - 
    """
    D = copy(G)
    B = _partial_normal_form(G)
    D = B * G * B.T
    return B
from sage.rings.all import Zp
from sage.matrix.constructor import Matrix
from copy import copy
from sage.rings.finite_rings.integer_mod import mod
def diagonalize_p_adic(G,prime,precision=None):
    """
    Diagonalize a symmetric matrix over a p-adic field 
    
    Input:
    
        -- ``G`` - symmetric n x n matrix in `\QQ`
        -- ``p`` - prime number
        -- ``precision`` - defining precision if not set the minimal possible is taken.
    Output:
    
        - ``D`` the diagonalized matrix
        - ``U`` transformation matrix, i.e, D = U * G * U^T
    
    EXAMPLES::
    
        sage: 
    """
    #this calculation is local. Thus we can do it over the p-adics
    p = prime
    denom = G.denominator()
    G = G*denom
    if precision == None:
        precision = G.det().valuation(p) + 10
    R = Zp(p, prec = precision, type = 'fixed-mod', print_mode = 'digits')
    G = G.change_ring(R)
    D = copy(G)
    n = G.ncols()
    
    #transformation matrix
    U = Matrix.identity(R,n)
    
    if(p == 2):
        D, U = _diagonalize_2_adic(G)
        assert U*G*U.T == Matrix.block_diagonal(collect_small_blocks(D))
    else:
        D, U = _diagonalize_odd_adic(G)
        assert D == Matrix.diagonal(D.diagonal())
    assert U.determinant().valuation()==0
    D, U1 = _normalize_blocks(D)
    U = U1 * U
    assert U*G*U.T == Matrix.block_diagonal(collect_small_blocks(D))
    return D, U

def _diagonalize_odd_adic(G):
    """
    """
    R = G.base_ring()
    D = copy(G)
    n = G.ncols()
    
    #transformation matrix
    U = Matrix.identity(R,n)
    
    #indices of the diagonal entrys which are already used
    cnt = 0
    while cnt < n:
        pivot = _find_min_p(D,cnt)
        piv1 = pivot[1]
        piv2 = pivot[2]
        minval = pivot[0]
        assert D==D.T
        #the smallest valuation is on the diagonal
        if piv1 == piv2:
            #move pivot to position [cnt,cnt]
            if piv1 != cnt:
                U.swap_rows(cnt,piv1)
                D.swap_rows(cnt,piv1)
                D.swap_columns(cnt,piv1)
            #we are already orthogonal to the part with i < cnt
            #now make the rest orthogonal too
            for i in range(cnt+1,n):
                if D[i,cnt]!= 0:
                    c = D[i,cnt]//D[cnt,cnt]
                    U[i,:] += -c * U[cnt,:]
                    D[i,:] += -c * D[cnt,:]
                    D[:,i] += -c * D[:,cnt]
            cnt = cnt + 1
        #the smallest valuation is off the diagonal
        else:
            #after this the diagonal entry will have smallest valuation.
            row = pivot[1]
            col = pivot[2]
            U[row,:] += U[col,:]
            D[row,:] += D[col,:]
            D[:,row] += D[:,col]
    return D, U

def _diagonalize_2_adic(G):
    """
    Diagonalize a symmetric matrix over the 2-adics
    
    (This method is called implicitly in diagonalize_p_adic(G,prime) )
    
    Input:

        - ``G`` symmetric nxn matrix in `\Qp`
    
    Output:

        - ``D`` the diagonalized matrix
        - ``U`` transformation matrix i.e D = U * G * U^T
    
    Example:

        sage: G = matrix(QQ,3,3,[1,0,1,0,2,1,1,1,3])

    
    """
    R = G.base_ring()
    D = copy(G)
    n = G.ncols()
    
    #transformation matrix
    U = Matrix.identity(R,n)
    
    #indices of the diagonal entrys which are already used
    cnt = 0
    minval = None
    while cnt < n:
            pivot = _find_min_p(D,cnt)
            piv1 = pivot[1]
            piv2 = pivot[2]
            minval_last = minval
            minval = pivot[0]
            assert D==D.T
            #the smallest valuation is on the diagonal
            if piv1 == piv2:
                #move pivot to position [cnt,cnt]
                if piv1 != cnt:
                    U.swap_rows(cnt,piv1)
                    D.swap_rows(cnt,piv1)
                    D.swap_columns(cnt,piv1)
                #we are already orthogonal to the part with i < cnt
                #now make the rest orthogonal too
                for i in range(cnt+1,n):
                    if D[i,cnt]!= 0:
                        c = D[i,cnt]//D[cnt,cnt]
                        U[i,:] += -c * U[cnt,:]
                        D[i,:] += -c * D[cnt,:]
                        D[:,i] += -c * D[:,cnt]
                cnt = cnt + 1
            #the smallest valuation is off the diagonal
            else:
                #move this 2 x 2 block to the top left (starting from cnt)
                if piv1 != cnt:
                    U.swap_rows(cnt,piv1)
                    D.swap_rows(cnt,piv1)
                    D.swap_columns(cnt,piv1)
                if piv2 != cnt+1:
                    U.swap_rows(cnt+1,piv2)
                    D.swap_rows(cnt+1,piv2)
                    D.swap_columns(cnt+1,piv2)
                #we split off a 2 x 2 block
                if cnt != n-2: 
                #else this is the last 2 x 2 block so there is nothing to do.
                    content = R(2**minval)
                    eqn_mat = D[cnt:cnt+2,cnt:cnt+2].list()
                    eqn_mat = Matrix(R,2,2,[e//content for e in eqn_mat]) #do you know better syntax here?
                    #calculate the inverse without using division
                    inv = eqn_mat.adjoint() * eqn_mat.det().inverse_of_unit()
                    B = U[cnt:cnt+2,:]
                    C = D[cnt+2:,cnt:cnt+2]*inv
                    for i in range(C.nrows()):
                        for j in range(C.ncols()):
                            C[i,j]=C[i,j]//content 
                    U[cnt+2:,:] -= C*B
                    D[cnt:,cnt:] = U[cnt:,:] * G * U[cnt:,:].transpose()
                #if the previous diagonal entry has the same valuation
                #as the pivot we can split off a 1 x 1 matrix
                if cnt > 0 and D[cnt-1,cnt-1].valuation() == minval == minval_last:
                    #a 3 x 3 block
                    # a  0  0
                    # 0 2b  c
                    # 0  c 2d
                    a = D[cnt-1,cnt-1].unit_part()
                    c = D[cnt,cnt+1].unit_part()
                    B = Matrix(R,3,3,[
                                    1,1, 0,
                                    c,0,-a,
                                    0,0, 1
                                    ]) #this matrix is invertible
                    U[cnt-1:cnt+2,:] = B * U[cnt-1:cnt+2,:]
                    D = U * G * U.transpose()
                    #at this stage cnt-1 and cnt are orthogonal and both have minimal valuation
                    cnt -= 3
                cnt += 2
    return D, U

def _get_small_block_indices(G):
    """
    """
    L = []
    n = G.ncols()
    i = 0
    while i < n-1:
        L.append(i)
        if G[i,i+1]!=0:
            i += 2
        else:
            i += 1
    if i == n-1:
        L.append(i)
    return L[1:]

def collect_small_blocks(G):
    """
    """
    L = _get_small_block_indices(G)
    G.subdivide(L,L)
    blocks = []
    for i in range(len(L)+1):
        blocks.append(G.subdivision(i,i))
    return blocks

def get_jordan_blocks():
    """
    """
    pass
    


def _normalize_blocks(G):
    """
    """
    R = G.base_ring()
    D = copy(G)
    p = R.prime()
    n = G.ncols()
    U = copy(G.parent().identity_matrix())
    if p != 2:
        #squareclasses 1,v
        v = _min_nonsquare(p)
        v = R(v)
        for i in range(n):
            d = D[i,i].unit_part()
            if d.is_square():
                D[i,i] = 1
                U[i,:] *= d.inverse_of_unit().sqrt()
            else:
                D[i,i] = v
                U[i,:] *= (v*d.inverse_of_unit()).sqrt()
    else:
        #squareclasses 1,3,5,7 modulo 8
        for i in range(n):
            d = D[i,i].unit_part()
            if d != 0:
                v = R(d.mod(8))
                D[i,i] = v*2**D[i,i].valuation()
                U[i,:] *= (v*d.inverse_of_unit()).sqrt()
        for i in range(n-1):
            b = D[i,i+1]
            if b !=0: #there is a 2 x 2 block here
                block = D[i:i+2,i:i+2]
                trafo = _normalize_2x2(block)
                U[i:i+2,:] = trafo * U[i:i+2,:]
    D = U * G * U.T
    return D, U

def _normalize_2x2(G):
    """
    normalize this 2x2 block
    
    INPUT:
    
        --  ``G`` - a 2 by 2 matrix over Zp with ``type = 'fixed-mod'`` of the form
            [2a  b]
            [ b 2c]
            with b of smallest valuation
    OUTPUT:
    
        -- either 2^k * [0 1]
                        [1 0]
                        
        -- or     2^k * [2 1]
                        [1 2]
    
    EXAMPLES::
        sage: from sage.quadratic_forms.genera.p_adic_jordan_block import _normalize_2x2
        sage: R = Zp(2, prec = 10, type = 'fixed-mod', print_mode = 'digits')
        sage: G = Matrix(R,2,[-17*2,3,3,23*2])
        sage: _normalize_2x2(G)
        sage: G = Matrix(R,2,[-17*4,3,3,23*2])
        sage: _normalize_2x2(G)
        sage: G = Matrix(R,2,[-17*2,3,3,23*2])
        sage: _normalize_2x2(8*G)
    
    TESTS::
    
        sage: from sage.quadratic_forms.genera.p_adic_jordan_block import _normalize_2x2
        sage: R = Zp(2, prec = 10, type = 'fixed-mod', print_mode = 'digits')
        sage: ref1 = Matrix(R,2,[2,1,1,2])
        sage: ref2 = Matrix(R,2,[0,1,1,0])
        sage: N = _normalize_2x2(G)
        sage: (N == ref1) or (N == ref2)
        True
    """
    
    T = copy(G.parent().identity_matrix())
    D = copy(G)
    B = G.base_ring()
    from sage.rings.all import PolynomialRing
    R = PolynomialRing(B,'x')
    x = R.gen()
    
    #The input is an even! block
    assert G[0,0].valuation() > G[1,0].valuation(), G
    assert G[1,1].valuation() > G[1,0].valuation()
    
    #Make sure G[1,1] has the smallest possible valuation.
    if G[0,0].valuation()<G[1,1].valuation():
        T.swap_rows(0,1)
        D.swap_rows(0,1)
        D.swap_columns(0,1)
    if D[1,1].valuation() > D[0,1].valuation()+1:
        T[1,:] += T[0,:]
        D = T*G*T.transpose()
    scale = 2**D[0,1].valuation()

    if mod(D.det().unit_part(),8) == 3:
        #transform it to 2^k *
        # 2 1
        # 1 2
        #Find a point of norm 2
        #solve: 2 == D[1,1]*x^2 + 2*D[1,0]*x + D[0,0]
        a = D[1,1]; b = 2*D[1,0]; c = (D[0,0]-2*scale);
        if a.valuation()>c.valuation():
            t=a
            a=c
            c=t
            T.swap_rows(0,1)
            D.swap_rows(0,1)
            D.swap_columns(0,1)
        pol = (a*x**2 + b*x + c)//(2*scale)
        sol = pol.roots()[0][0]
        T[0,1] = sol
        D = T*G*T.transpose()
        #both should pair with value 1
        T[1,:] *= D[1,0].unit_part().inverse_of_unit()
        D = T*G*T.transpose()
        #solve: v*D*v == 2 with v as below
        from sage.modules.free_module_element import vector
        v = vector([x,-2*x+1])
        pol = (v*D*v-2)//(2*scale)
        sol = pol.roots()[0][0]
        T[1,:] = sol * T[0,:] + (- 2*sol+1)*T[1,:]
        D = T*G*T.transpose()
        assert D == scale*Matrix(G.parent(),2,2,[2,1,1,2]), "D1 \n %r" %D
    elif mod(D.det().unit_part(),8) == 7:
        #transform it to 2^k *
        # 0 1
        # 1 0
        #solve: 0 == D[1,1]*x^2 + 2*D[1,0]*x + D[0,0]
        #Find a point representing 0
        #solve: 0 == D[1,1]*x^2 + 2*D[1,0]*x + D[0,0]
        sol = (D[1,1]*x**2 + 2*D[1,0]*x + D[0,0]).roots()[0][0]
        T[0,:] += sol*T[1,:]
        D = T*G*T.transpose()
        #make the second basis vector have 0 square as well.
        T[1,:] = T[1,:] - D[1,1]//(2*D[0,1])* T[0,:]
        D = T*G*T.transpose()
        #both should pair with value 1
        T[0,:] *= D[1,0].unit_part().inverse_of_unit()
        D = T*G*T.transpose()
        assert D == scale*Matrix(G.parent(),2,2,[0,1,1,0]), "D2 \n %r" %D
    return T
        
def _find_min_p(G,cnt):
    """
    Find the smallest valuation and prefer diagonal entries.
    
    Input:

        -- ``G`` - symmetric n x n matrix in `\Qp`
        -- ``cnt`` - start search from this index
    
    Output:

        - ``min`` minimal valuation
        - ``min_i`` row of the minimal valuation
        - ``min_j`` column of the minimal valuation
    
    EXAMPLES::
    
        sage: from sage.quadratic_forms.genera.p_adic_jordan_block import _find_min_p
        sage: G = matrix(Qp(3),3,3,[4,0,1,0,4,2,1,2,1])
        sage: G
        [1 + 3 + O(3^20)               0     1 + O(3^20)]
        [              0 1 + 3 + O(3^20)     2 + O(3^20)]
        [    1 + O(3^20)     2 + O(3^20)     1 + O(3^20)]
        sage: _find_min_p(G,2,0)
        (0,3,3)
    
    """
    n = G.ncols()
    min = G[cnt,cnt].valuation()
    min_i = cnt
    min_j = cnt
    for i in range(cnt,n):
        for j in range(i,n):
            if (G[i,j]!=0) and (G[i,j].valuation() < min):
                min_i = i
                min_j = j
                min = G[i,j].valuation()
        #diagonal has precedence
        if (G[i,i]!=0) and (G[i,i].valuation() <= min):
            min_i = i
            min_j = i
            min = G[i,i].valuation()
    return min, min_i, min_j

def _min_nonsquare(prime):
    """
    Calculate minimal nonsquare in `\FF_p`
    
    Input:

        - ``p`` prime number for `\Qp`
    
    Output:
    
        - ``a`` minimal nonsquare

    """
    from sage.rings.all import GF
    R = GF(prime)
    for i in R:
        if not R(i).is_square():
            return i

r"""

AUTHORS:

- Simon Brandhorst (2017-09): First created (based on free_module.py)
"""

#*****************************************************************************
#       Copyright (C) 2017 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.modules.fg_pid.fgp_module import FGP_Module_class
from sage.modules.fg_pid.fgp_element import DEBUG, FGP_Element
from sage.modules.free_module import FreeModule
from sage.arith.misc import gcd
from sage.rings.all import ZZ, QQ
from sage.groups.additive_abelian.qmodnz import QmodnZ
from sage.matrix.constructor import matrix


class TorsionQuadraticModuleElement(FGP_Element):
    """
    An element of a torsion quadratic module.
    
    TESTS::

        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: T = TorsionQuadraticModule(ZZ^3,6*ZZ^3)
        sage: loads(dumps(T)) == T
        True
        sage: t = T.0
        sage: loads(dumps(t)) == t
        True    
    """
    def __init__(self, parent, x, check=DEBUG):
        FGP_Element.__init__(self,parent=parent,x=x, check=check)

    def _mul_(self, other):
        r"""
        Compute the inner product of two elements

        OUTPUT:
        
        -- an element of `\Q/m\Z` with `m\Z=b(V,W)` 
        
        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = (1/2)*ZZ^2; W = ZZ^2
            sage: T = TorsionQuadraticModule(V,W)
            sage: x = T.0 
            sage: y = T.0 + T.1
            sage: x
            (1, 0)
            sage: x*y
            1/4
        """
        return self.parent().value_module()(self.lift().inner_product(other.lift()))

    def inner_product(self, other):
        """
        Return the inner product of this element with the other element.
        
        Input:
        
        -- ``other`` an element of the same torsion quadratic module.
        
        Output:
        
        -- an element of `\Q/m\Z` with `m\Z = <V,W>` 

        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = (1/2)*ZZ^2; W = ZZ^2
            sage: T = TorsionQuadraticModule(V,W)
            sage: x = T.0 
            sage: y = T.0 + T.1
            sage: x
            (1, 0)
            sage: x.inner_product(y)
            1/4
        """
        return self * other

    def q(self):
        r"""
        Compute the quadratic_product of self.

        OUTPUT:

        -- an element of `\Q/n\Z`
        
        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: W = FreeQuadraticModule(ZZ,2,2*matrix.identity(2))
            sage: V = (1/2)*W
            sage: T = TorsionQuadraticModule(V,W)
            sage: x = T.0;
            sage: x
            (1, 0)
            sage: x.q()
            1/2
            sage: x.q().parent()
            Q/2Z
            sage: x*x
            1/2
            sage: (x*x).parent()
            Q/Z
        """
        
        return self.parent().value_module_qf()(self.lift().inner_product( self.lift()))

    def quadratic_product(self):
        r"""
        Compute the quadratic_product of self.

        OUTPUT:

        -- an element of `\Q/n\Z` where `n\Z= <V,W> +\Z \{ <w,w> | w \in W \}`

        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: W = FreeQuadraticModule(ZZ,2,2*matrix.identity(2))
            sage: V = (1/2)*W
            sage: T = TorsionQuadraticModule(V,W)
            sage: x = T.0;
            sage: x
            (1, 0)
            sage: x.quadratic_product()
            1/2
            sage: x.quadratic_product().parent()
            Q/2Z
            sage: x*x
            1/2
            sage: (x*x).parent()
            Q/Z
        
        
        """
        return self.q()

class TorsionQuadraticModule(FGP_Module_class):
    r"""
    Let `V` be a symmetric FreeQuadraticModule and `W\subseteq V` a submodule of the same rank as `V`. The quotient `V/W` is a TorsionQuadraticModule. It inherits a bilinear form `b` and a quadratic form `q`.
    ..MATH::
    
    `b: V \times V \rightarrow \Q /m\Z \mbox{ where } m\Z = <V,W> and b(x,y) = <x,y> + m\Z.
    `q: V \rightarrow \Q/n\Z` where `n\Z= <V,W> +\Z \{ <w,w> | w \in W \}`.

    INPUT:

    - ``V`` -- a FreeModule with a symmetric inner_product_matrix
    - ``W`` -- a submodule of V of the same rank as V
    - ``check`` -- bool (default: True)
    - ``modulus`` -- a rational number dividing `m` (default: m). The inner product `b` is defined in `\Q / modulus\Z`
    - ``modulus_qf`` -- a rational number divinding `n`(default: n). The quadratic form `q` is defined in `\Q / modulus_qf\Z`.

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import *
        sage: V = FreeModule(ZZ,3)
        sage: T = TorsionQuadraticModule(V, 5*V)
        sage: T
        Finite quadratic module V/W over Integer Ring with invariants (5, 5, 5).
        Gram matrix of the quadratic form with values in Q/5Z:
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    # The class to be used for creating elements of this
    # module. Should be overridden in derived classes.
    Element = TorsionQuadraticModuleElement
    def __init__(self, V, W, gens=None, check=True, modulus=None,modulus_qf=None):
        r"""
        Constructor for TorsionQuadraticModules
        """
        if check:
            if V.rank() != W.rank():
                raise ValueError("Modules must be of the same rank.")
            if V.base_ring() is not ZZ:
                raise NotImplementedError("Only TorsionQuadraticModules over ZZ are implemented")
            if V.inner_product_matrix() != V.inner_product_matrix().transpose():
                raise ValueError("The cover must have a symmetric inner_product")

            if gens is not None and (V.span(gens) + W != V):
                raise ValueError("Provided gens do not generate quotient")

        FGP_Module_class.__init__(self, V, W, check=check)
        if gens is None:
            self._gens = FGP_Module_class.gens(self)
        else:
            self._gens = [self(v) for v in gens]

        #The inner inner_product of two elements `b(v1+W,v2+W)` is defined `mod <V,W>`
        self._modulus = gcd([x.inner_product(y) for x in V.gens() for y in W.gens()])
        #The quadratic_product of an element `q(v+W)` is defined `\mod 2<V,W> + ZZ\{ <w,w> | w in w\}`
        norm = gcd(self.W().gram_matrix().diagonal())
        self._modulus_qf = gcd(norm, 2 * self._modulus)
        
        if modulus is not None:
            if self._modulus/modulus not in self.base_ring():
                raise ValueError, "The modulus must divide <V,W>."
            self._modulus = modulus
        if modulus_qf is not None:
            if self._modulus_qf/modulus_qf not in self.base_ring():
                raise ValueError, "The modulus must divide <V,W>."
            self._modulus_qf = modulus_qf

    def _repr_(self):
        r"""
        The print representation of this TorsionQuadraticModule.
        """
        return ("Finite quadratic module V/W over %s with invariants %s.\n"%(self.base_ring(),self.invariants()) +
            "Gram matrix of the quadratic form with values in %r:\n%r" %(self.value_module_qf(),self.gram_matrix_quadratic()))
    
    # The method to be used for creating modules - for example submodules.
    # Should be overridden in derived classes.
    def _module_constructor(self, V, W, check=True):
        r"""
        Construct a torsion quadratic module``V/W``.

        INPUT:

        - ``V`` -- an R-module.

        - ``W`` -- an R-submodule of ``V``.

        - ``check`` -- bool (default: True).

        OUTPUT:

        The quotient ``V/W`` as TorsionQuadraticModule.

        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]],ZZ); W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2])
            sage: Q = TorsionQuadraticModule(V,W); Q
            Finite quadratic module V/W over Integer Ring with invariants (4, 12).
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
            sage: Q._module_constructor(V,W)
            Finite quadratic module V/W over Integer Ring with invariants (4, 12).
            Gram matrix of the quadratic form with values in Q/(1/4)Z:
            [0 0]
            [0 0]
        """
        return TorsionQuadraticModule(V, W, check=check)

    def gram_matrix_bilinear(self):
        r"""
        The gram matrix with respect to the generators

        OUTPUT: 
            a rational matrix G with `G_{i,j}` given by the inner product 
            of the `i`-th and `j`-th generator. Its entries are only well defined
            `\mod <V,W>`
        
        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T.gram_matrix_bilinear()
            [1/5   0   0]
            [  0 1/5   0]
            [  0   0 1/5]
        """
        gens = self._gens
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i+1):
                G[i,j] = G[j,i] = (gens[i] * gens[j]).lift()
        return G

    def gram_matrix_quadratic(self):
        r"""
        The gram matrix of the quadratic form with respect to the current generators

        OUTPUT: 
        
            a rational matrix Gq with Gq[i,j]=<gens[i],gens[j]> and G[i,i] = gens[i].q()

        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: D4_gram = Matrix(ZZ,4,4,[2,0,0,-1,0,2,0,-1,0,0,2,-1,-1,-1,-1,2])
            sage: D4 = FreeQuadraticModule(ZZ,4,D4_gram)
            sage: D4dual = D4.span(D4_gram.inverse())
            sage: discrForm = TorsionQuadraticModule(D4dual,D4)
            sage: discrForm.gram_matrix_quadratic()
            [  1 1/2]
            [1/2   1]
            sage: discrForm.gram_matrix_bilinear()
            [  0 1/2]
            [1/2   0]

        """
        gens = self._gens
        n = len(gens)
        Q = self.base_ring().fraction_field()
        G = matrix.zero(Q, n)
        for i in range(n):
            for j in range(i):
                G[i,j] = G[j,i] = (gens[i] * gens[j]).lift()
            G[i,i] = gens[i].q().lift()
        return G
    
    def gens(self):
        """
        Return generators of this module.
        
        There is no assumption on the generators except that they generate the module.
        
        Examples::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ,3)
            sage: T = TorsionQuadraticModule(V, 5*V)
            sage: T.gens()
            ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        """
        return self._gens
    
    def orthogonal_submodule_to(self, S):
        r"""
        Return the submodule orthogonal to ``S``.

        INPUT:

        - ``S`` -- a submodule
        
        EXAMPLES::
        
        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = FreeModule(ZZ,10)
        sage: T = TorsionQuadraticModule(V,3*V)
        sage: S = T.submodule(T.gens()[:5])
        sage: O = T.orthogonal_submodule_to(S)
        sage: O
        Finite quadratic module V/W over Integer Ring with invariants (3, 3, 3, 3, 3).
        Gram matrix of the quadratic form with values in Q/3Z:
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        sage: O.V() + S.V() == T.V()
        True
        """
        if not S.is_submodule(self):
             raise ValueError, "S must be a submodule of this module."


        G = self.V().inner_product_matrix()
        T = self.V().basis_matrix()
        S = S.V().basis_matrix()
        m = self._modulus
        
        Y = (T*G*S.transpose())
        integral = Y.inverse()*T
        orthogonal = m * integral
        orthogonal = self.submodule(orthogonal.rows())
        return orthogonal
    
    def primary_part(self,m):
        """
        Return the m-primary part of this torsion quadraticm module as a submodule.
        
        INPUT:
        
        - ``m`` - an integer
        
        OUTPUT:
        
        - a submodule
        
        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: T = TorsionQuadraticModule((1/6)*ZZ^3,ZZ^3)
            sage: T.primary_part(2)
            Finite quadratic module V/W over Integer Ring with invariants (2, 2, 2).
            Gram matrix of the quadratic form with values in Q/(1/3)Z:
            [1/4   0   0]
            [  0 1/4   0]
            [  0   0 1/4]
    
        TESTS::
            sage: T==T.primary_part(T.annihilator().gen())
            True
        """
        annihilator = self.annihilator().gen()
        a = annihilator.prime_to_m_part(m)
        return self.submodule((a*self.V()).gens())

    def submodule(self, x):
        r"""
        Return the submodule defined by `x` as a TorsionQuadraticModule.

        The modulus of the inner product is inherited from this module.

        INPUT:

        - ``x`` -- list, tuple, or FGP module
        
        EXAMPLES::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T.submodule(T.gens()[:2])
            Finite quadratic module V/W over Integer Ring with invariants (5, 5).
            Gram matrix of the quadratic form with values in Q/Z:
            [1/5   0]
            [  0 1/5]
        """
        T = FGP_Module_class.submodule(self, x)
        T._modulus = self._modulus # is this necessary? Yes else the modulus might increase.
        T._modulus_qf = self._modulus_qf
        return T
    
    def value_module(self):
        """
        Return `\Q/m\Z` with `m= <V,W>`.

        This is where the inner product takes values.
        
        EXAMPLES::
        
            sage: from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
            sage: A2 = Matrix(ZZ,2,2,[2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module V/W over Integer Ring with invariants (2, 6).
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module()
            Q/Z
        """
        return QmodnZ(self._modulus)

    def value_module_qf(self):
        """
        Return `\Q/n\Z` with `n\Z= <V,W> +\Z \{ <w,w> | w \in W \}`.

        This is where the torsion quadratic form takes values.
        
        EXAMPLES::
        
            sage: from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
            sage: A2 = Matrix(ZZ,2,2,[2,-1,-1,2])
            sage: L = IntegralLattice(2*A2)
            sage: D = L.discriminant_group()
            sage: D
            Finite quadratic module V/W over Integer Ring with invariants (2, 6).
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 1/3]
            sage: D.value_module_qf()
            Q/2Z
        """
        return QmodnZ(self._modulus_qf)
    
    def gens_orthogonalize(self):
        """
        Return a free quadratic module such that its generators are mutually orthogonal for `p != 2` or have block size at most `2`.
        """
        gens = []
        for p in self.annihilator().support():
            Qp = self.primary_part(p)
            gens_p = _diagonalize(Qp)
            gens += gens_p
        return gens
    

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
    from sage.rings.all import Zp
    from sage.matrix.constructor import Matrix
    #this calculation is local. Thus we can do it over the p-adics
    p = prime
    denom = G.denominator()
    G = G*denom
    if precision == None:
        precision = G.det().valuation(p) + 10
    R = Zp(p, prec = precision, type = 'fixed-mod', print_mode = 'digits')
    G = G.change_ring(R)
    D = deepcopy(G)
    n = G.ncols()
    
    #transformation matrix
    U = Matrix.identity(R,n)
    
    if(p == 2):
        D, U = _diagonalize_2_adic(G)
        assert U*G*U.T == matrix.block_diagonal(collect_small_blocks(D))
    else:
        D, U = _diagonalize_odd_adic(G)
        assert D == Matrix.diagonal(D.diagonal())
    assert U.determinant().valuation()==0
    D, U1 = _normalize_blocks(D)
    U = U1 * U
    assert U*G*U.T == matrix.block_diagonal(collect_small_blocks(D))
    return D, U

def _diagonalize_odd_adic(G):
    """
    """
    R = G.base_ring()
    D = deepcopy(G)
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
    D = deepcopy(G)
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
    L = _get_small_block_indices(G)
    G.subdivide(L,L)
    blocks = []
    for i in range(len(L)+1):
        blocks.append(G.subdivision(i,i))
    return blocks

def get_jordan_blocks():
    pass
    


def _normalize_blocks(G):
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
    
        -- ``G`` - a 2 by 2 matrix over Qp with (type = 'fixed modulus')
    
    OUTPUT:
    
        -- either 2^k * (0,1,1,0) as a 2 x 2 matrix
        -- or     2^k * (2,1,1,2) as a 2 x 2 matrix
    """
    
    T = copy(G.parent().identity_matrix())
    D = deepcopy(G)
    B = G.base_ring()
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
    
    Example:

        sage: G = matrix(Qp(3),3,3,[4,0,1,0,4,2,1,2,1])
        sage: G
        [4 0 1]
        [0 4 2]
        [1 2 5]
        sage: find_min_p(G,2,0)
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
    R = sage.rings.all.GF(prime)
    for i in R:
        if not R(i).is_square():
            return i

def _something(prime):
    """
    Calculate minimum nonsquare in `\Zp` which satisfies `(1+b^2)` not a square
    
    Input:

        - ``prime`` prime number for `\Zp`
    
    Output:
    
        - ``i`` minimum nonsquare

    """
    R = sage.rings.all.GF(prime)
    for i in R:
        if not (i**2+R(1)).is_square():
            return i

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
    
    def orthogonal_gens(self, normalize = True):
        """
        Return a free quadratic module such that its generators are mutually orthogonal for `p != 2` or have block size at most `2`.
        """
        gens = []
        from sage.quadratic_forms.genera.p_adic_jordan_blocks import jordan_p_adic
        for p in self.annihilator().gen().prime_divisors():
            D_p = self.primary_part(p)
            q_p = D_p.gram_matrix_quadratic()
            q_p = q_p/D_p._modulus
            D,U = jordan_p_adic(q_p,p,normalize=normalize)
            #apply U to the generators
            n = U.ncols()
            U = U.change_ring(ZZ)
            gens_p = []
            for i in range(n):
                g = self.zero()
                for j in range(n):
                    g += self.gens()[j]*U[i,j]
                gens_p.append(g)
            gens += gens_p
            
        return self.submodule_with_gens(gens)
        
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

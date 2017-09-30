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
from sage.rings.all import ZZ,QQ


class TorsionQuadraticModuleElement(FGP_Element):
    def __init__(self, parent, x, check=DEBUG):
        FGP_Element.__init__(self,parent=parent,x=x, check=check)
    
    def inner_product(self,other):
        """
        Compute the inner_product of self and other.
        
        Output:
        
        A rational number - a representative of the equivalence class of <self,other> modulo <V,W> \ZZ.
        """
        b = self.lift().inner_product(other.lift())
        modulus = self.parent()._modulus
        t = b / modulus
        n = t.numerator() % t.denominator()
        x = modulus * n / t.denominator()
        return(x)
    
    def quadratic_product(self):
        """
        Compute the inner_product of self and other.
        
        Output:
        
        An element of the fraction field of R - a representative of the equivalence class of <x,x> modulo 2<V,W> R.
        """
        b = self.lift().inner_product(self.lift())
        modulus = 2 * self.parent()._modulus
        t = b / modulus
        n = t.numerator() % t.denominator()
        x = modulus * n / t.denominator()
        return(x)
    
    
class TorsionQuadraticModule(FGP_Module_class):
    """
    Let $V$ be a symmetric FreeQuadraticModule and $W\subseteq V$ a submodule of the same rank as $V$. The quotient $V/W$ is a TorsionQuadraticModule. It inherits a bilinear form 
    $b: V \times V \rightarrow QQ / m\ZZ$, where $m \ZZ = <V,W>$ and $b(x,y) = <x,y> + m\ZZ$. 
    
    Examples::
    
        sage: from sage.modules.torsion_quadratic_module import *
        sage: V = FreeModule(ZZ,3)
        sage: T = TorsionQuadraticModule(V,5*V)
        sage: T

        Finite quadratic module V/W over Integer Ring with invariants (5, 5, 5).
        Gram matrix:
        [1 0 0]
        [0 1 0]
        [0 0 1]

    """
    
    # The class to be used for creating elements of this
    # module. Should be overridden in derived classes.
    Element = TorsionQuadraticModuleElement
    def __init__(self,V,W,check=True,modulus=None):
        """
        Constructor for TorsionQuadraticModules
        
        Input:
        -- ``V`` a FreeModule with a symmetric inner_product_matrix
        -- ``W`` a submodule of V of the same rank
        -- ``check`` for extra checks on the input
        -- ``modulus`` a rational number. The inner product is defined in QQ / modulus\ZZ.
        
        Output:
        -- a TorsionQuadraticModule
        """
    
        if check:
            if V.rank()!=W.rank():
                raise ValueError, "Modules must be of the same rank."
            if V.base_ring()!=ZZ:
                raise NotImplementedError, "Only TorsionQuadraticModules over ZZ are implemented"
            if V.inner_product_matrix()!=V.inner_product_matrix().transpose():
                raise ValueError, "The cover must have a symmetric inner_product"
        
        FGP_Module_class.__init__(self,V,W,check=check)
        
        self._current_gens = FGP_Module_class.gens(self)
        
        #The inner inner_product of two elements <v1+W,v2+W> is defined mod <V,W> 
        self._modulus = gcd([x.inner_product(y) for x in V.gens() for y in W.gens()])
        if modulus != None:
            if not self._modulus/modulus in self.base_ring():
                raise ValueError, "The modulus must divide <V,W>."
            self._modulus = modulus
            
                
    @property
    def current_gens(self):
        """
        Generators of self specified by the user
        """
        return(self._current_gens)
    
    def _set_gens(self,gens):        
	gens = tuple(self(x) for x in gens)
	if not self.is_submodule(self.submodule(gens)):
            raise ValueError, str(gens) + " does not generate self"
        self._current_gens = gens
    current_gens = current_gens.setter(_set_gens)

    def _repr_(self):
        """
        The printing representation of self.
        """
        s = "Finite quadratic module V/W over %s with invariants %s.\n"%(self.base_ring(),self.invariants()) + \
            "Gram matrix:\n%r" %self.gram_matrix_quadratic() 
        return(s)
    
    def submodule(self,x):
        """
        Return the submodule defined by x as a TorsionQuadraticModule 
        
        The modulus of the inner product is inherited from self.

        INPUT:

        - ``x`` -- list, tuple, or FGP module

        """
        
        #Do all the input checks/conversion in the FGP_Module_class this might be slow
        temp = FGP_Module_class.submodule(self,x)
        return(TorsionQuadraticModule(temp.V(),temp.W(),modulus=self._modulus))

    def gram_matrix_bilinear(self):
        """
        The gram matrix with respect to the current generators
        
        Output:
        -- a rational matrix G with G[i,j]=<gens[i],gens[j]>
        """

        from sage.matrix.special import matrix
        n = len(self.gens())
        G = matrix.zero(self.base_ring().fraction_field(),n)
        for i in range(0,n):
            for j in range(0,n):
                G[i,j] = self.gens()[i].inner_product(self.gens()[j])
        return(G)
        
    def gram_matrix_quadratic(self):
        """
        The gram matrix of the quadratic form with respect to the current generators
        
        Output:
        -- a rational matrix
        """
        
        from sage.matrix.special import matrix
        n = len(self.gens())
        G = matrix.zero(self.base_ring().fraction_field(),n)
        for i in range(0,n):
            for j in range(0,i):
                G[i,j] = self.gens()[i].inner_product(self.gens()[j])
        G = G + G.transpose()
        for i in range(0,n):
            G[i,i] = self.gens()[i].quadratic_product()
        return(G)
        
    def orthogonal_submodule_to(self,S):
        """
        Return the orthogonal subspace of self to S
        
        Input:
        --``S`` - a submodule of self
        
        Output:
        -- A TorsionQuadraticModule
        """
        from sage.matrix.special import matrix
        if not S.is_submodule(self):
             raise ValueError, "S must be a submodule of self."

        n = self.V().rank()
        m = S.V().rank()
        B = self.V().basis_matrix()
        basisS = S.V().basis_matrix()
        s = self._modulus
        
        M = B*self.V().inner_product_matrix()*basisS.transpose()/s
        M,d = M._clear_denom()
        D,U,V = M.smith_form()
        D = D/d
        
        X = matrix.identity(self.base_ring().fraction_field(),n,n)
        for i in range(0,m):
            if D[i,i]!=0:
                X[i,i] = 1 / D[i,i]
        X = X * U * B
        orth = self.V().span(X.rows()).intersection(self.V())
        orth = self.submodule(orth.gens())
        
        return(orth)


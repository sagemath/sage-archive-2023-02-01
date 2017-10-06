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
    def __init__(self, parent, x, check=DEBUG):
        FGP_Element.__init__(self,parent=parent,x=x, check=check)

    def _mul_(self, other):
        r"""
        Compute the inner product of two elements

        INPUT:

        - ``other`` -- 

        OUTPUT:

        a rational number which is a representative of the equivalence
        class of <self,other> modulo <A,B> over `\Z`.
        """
        return self.parent().value_module()(self.lift() * other.lift())

    def inner_product(self, other):
        """
        Return the inner product of this element with the other element.
        
        Input:
        
        -- ``other`` an element of the same torsion quadratic module.
        
        Output:
        
        -- an element of `\Q/n\Z` with `n=` 
        """
        return self * other

    def q(self):
        return self.parent().value_module_qf()(self.lift() * self.lift())

    def quadratic_product(self):
        r"""
        Compute the quadratic_product of self.

        OUTPUT:

        An element of the fraction field of R which is a
        representative of the equivalence class of <x,x> modulo 2<V,W> R.
        """
        return self.q()

class TorsionQuadraticModule(FGP_Module_class):
    r"""
    Let `V` be a symmetric FreeQuadraticModule and `W\subseteq V` a submodule of the same rank as `V`. The quotient `V/W` is a TorsionQuadraticModule. It inherits a bilinear form `b` where
    `b: V \times V \rightarrow \Q /m\Z`$, where m\Z = <V,W> and b(x,y) = <x,y> + m\Z.

    INPUT:

    - ``V`` -- a FreeModule with a symmetric inner_product_matrix
    - ``W`` -- a submodule of V of the same rank as V
    - ``check`` -- bool (default: True)
    - ``modulus`` -- a rational number. The inner product is defined in `\Q / modulus\Z`.

    EXAMPLES::

        sage: from sage.modules.torsion_quadratic_module import *
        sage: V = FreeModule(ZZ,3)
        sage: T = TorsionQuadraticModule(V, 5*V)
        sage: T
        Finite quadratic module V/W over Integer Ring with invariants (5, 5, 5).
        Gram matrix of the quadratic form:
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    # The class to be used for creating elements of this
    # module. Should be overridden in derived classes.
    Element = TorsionQuadraticModuleElement
    def __init__(self, V, W, gens=None, check=True, modulus=None):
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

        #The inner inner_product of two elements <v1+W,v2+W> is defined mod <V,W>
        self._modulus = gcd([x.inner_product(y) for x in V.gens() for y in W.gens()])
        if modulus is not None:
            if self._modulus/modulus not in self.base_ring():
                raise ValueError, "The modulus must divide <V,W>."
            self._modulus = modulus

    def value_module(self):
        """
        This is where the inner product takes values.
        
        Returns `\Q/m\Z` with `m= <V,W>`.
        """
        return QmodnZ(self._modulus)

    def value_module_qf(self):
        """
        This is where the torsion quadratic form takes values.
        
        Returns `\Q/2m\Z` with `m= <V,W>`.
        """
        return QmodnZ(2 * self._modulus)

    def gens(self):
        """
        Return generators of this module.
        
        There is no assumption on the generators except that they generate the module.
        
        Examples::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeModule(ZZ,3)
            sage: T = TorsionQuadraticModule(V, 5*V)
            sage: T.gens()
        """
        return self._gens

    def _repr_(self):
        r"""
        The print representation of this TorsionQuadraticModule.
        """
        return ("Finite quadratic module V/W over %s with invariants %s.\n"%(self.base_ring(),self.invariants()) +
            "Gram matrix of the quadratic form with values in %r:\n%r" %(self.value_module_qf(),self.gram_matrix_quadratic()))

    def _module_constructor(self, V, W, check=True):
        r"""
        Construct a quotient module ``V/W``.

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
            Gram matrix of the quadratic form with values in Q/(1/2)Z:
            [0 0]
            [0 0]
            sage: Q._module_constructor(V,W)
            Finite quadratic module V/W over Integer Ring with invariants (4, 12).
            Gram matrix of the quadratic form with values in Q/(1/2)Z:
            [0 0]
            [0 0]

        """
        return TorsionQuadraticModule(V, W, check=check)

    def submodule(self, x):
        r"""
        Return the submodule defined by `x` as a TorsionQuadraticModule.

        The modulus of the inner product is inherited from this module.

        INPUT:

        - ``x`` -- list, tuple, or FGP module
        
        Examples::
        
            sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
            sage: V = FreeQuadraticModule(ZZ,3,matrix.identity(3)*5)
            sage: T = TorsionQuadraticModule((1/5)*V, V)
            sage: T.submodule(T.gens()[:2])
        """
        T = FGP_Module_class.submodule(self, x)
        T._modulus = self._modulus # is this necessary? Yes the modulus might change.
        return T

    def gram_matrix_bilinear(self):
        r"""
        The gram matrix with respect to the generators

        OUTPUT: a rational matrix G with G[i,j]=<gens[i],gens[j]>
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

        OUTPUT: a rational matrix
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

    def orthogonal_submodule_to(self, S):
        r"""
        Return the submodule orthogonal to ``S``.

        INPUT:

        - ``S`` -- a submodule
        sage: from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        sage: V = FreeModule(ZZ,10)
        sage: T = TorsionQuadraticModule(V,3*V)
        sage: S = T.submodule(T.gens()[:5])
        sage: O = T.orthogonal_submodule_to(S)
        sage: O
        Finite quadratic module V/W over Integer Ring with invariants (3, 3, 3, 3, 3).
        Gram matrix of the quadratic form with values in Q/6Z:
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        sage: O + S
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

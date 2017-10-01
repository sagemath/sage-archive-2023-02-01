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
        return self.parent().value_group()(self.lift() * other.lift())

    def inner_product(self, other):
        return self * other

    def q(self):
        return self.parent().value_group_qform()(self.lift() * self.lift())

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
        Gram matrix:
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

    def value_group(self):
        # Value group of inner product
        return QmodnZ(self._modulus)

    def value_group_qform(self):
        # Value group of quadratic form
        return QmodnZ(2 * self._modulus)

    def gens(self):
        return self._gens

    def _repr_(self):
        r"""
        Print representation.
        """
        return ("Finite quadratic module V/W over %s with invariants %s.\n"%(self.base_ring(),self.invariants()) +
            "Gram matrix:\n%r" %self.gram_matrix_quadratic())

    def _module_constructor(self, V, W, check=True):
        r"""
        Construct a quotient module ``V/W``.

        INPUT:

        - ``V`` -- an R-module.

        - ``W`` -- an R-submodule of ``V``.

        - ``check`` -- bool (default: True).

        OUTPUT:

        The quotient ``V/W``.

        EXAMPLES::

            sage: 
        """
        return TorsionQuadraticModule(V, W, check)

    def submodule(self, x):
        r"""
        Return the submodule defined by `x` as a TorsionQuadraticModule.

        The modulus of the inner product is inherited from this module.

        INPUT:

        - ``x`` -- list, tuple, or FGP module

        """
        T = FGP_Module_class.submodule(self, x)
        T._modulus = self._modulus # is this necessary?
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
                G[i,j] = G[j,i] = gens[i] * gens[j]
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
                G[i,j] = G[j,i] = gens[i] * gens[j]
            G[i,i] = gens[i].q()
        return G

    def orthogonal_submodule_to(self, S):
        r"""
        Return the submodule orthogonal to ``S``.

        INPUT:

        - ``S`` -- a submodule
        """
        if not S.is_submodule(self):
             raise ValueError, "S must be a submodule of this module."

        n = self.V().rank()
        m = S.V().rank()
        B = self.V().basis_matrix()
        basis_S = S.V().basis_matrix()
        s = self._modulus

        M = B*self.V().inner_product_matrix()*basis_S.transpose() / s
        M, d = M._clear_denom()
        D, U, V = M.smith_form()
        D = D / d

        X = matrix.identity(self.base_ring().fraction_field(), n)
        for i in range(m):
            if D[i,i] != 0:
                X[i,i] = 1 / D[i,i]
        X = X * U * B
        orth = self.V().span(X.rows()).intersection(self.V())
        orth = self.submodule(orth.gens())

        return orth


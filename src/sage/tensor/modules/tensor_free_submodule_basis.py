r"""
Free module bases indexed by component indices
"""

#******************************************************************************
#       Copyright (C) 2020 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_basis import Basis_abstract
from sage.tensor.modules.comp import Components

class TensorFreeSubmoduleBasis_comp(Basis_abstract):
    r"""
    Standard basis of a tensor module with prescribed symmetries.

    EXAMPLES::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: T11 = M.tensor_module(1,1)
        sage: e11 = T11.basis('e')
        sage: for a in e11: a.display()
        e_0*e^0
        e_0*e^1
        e_0*e^2
        e_1*e^0
        e_1*e^1
        e_1*e^2
        e_2*e^0
        e_2*e^1
        e_2*e^2

    """

    def __init__(self, tensor_module, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, symbol_dual=None, latex_symbol_dual=None):
        base_module = tensor_module.base_module()
        base_module_basis = base_module.basis(symbol, latex_symbol, indices,
                                              latex_indices, symbol_dual, latex_symbol_dual)
        super().__init__(tensor_module, symbol, latex_symbol, indices, latex_indices)
        self._base_module_basis = base_module_basis
        self._comp = tensor_module._basis_comp()

    def keys(self):
        yield from self._comp.non_redundant_index_generator()

    def __iter__(self):
        r"""
        Generate the basis elements of ``self``.
        """
        for ind in self.keys():
            yield self[ind]

    def __getitem__(self, index):
        r"""
        Return the basis element corresponding to a given index.

        INPUT:

        - ``index`` -- the index of the basis element

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: T11 = M.tensor_module(1,1)
            sage: e11 = T11.basis('e')
            sage: e11[1, 2].display()
            e_1*e^2

            sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
            sage: Sym2M = TensorFreeSubmodule_comp(M, (2, 0), sym=range(2)); Sym2M
            Free module of type-(2,0) tensors
            with Fully symmetric 2-indices components w.r.t. [0, 1, 2]
            on the Rank-3 free module M over the Integer Ring
            sage: eSym2M = Sym2M.basis('e')
            sage: eSym2M[1, 1].display()
            e_1*e_1
            sage: eSym2M[1, 2].display()
            e_1*e_2 + e_2*e_1

        """
        tensor_module = self._fmodule
        base_module = tensor_module.base_module()
        base_module_basis = self._base_module_basis
        element = tensor_module([])
        element.set_comp(base_module_basis)[index] = 1
        return element

# Todo: Make it a Family
#       dual basis
#       add test for dual
# lift/reduce/retract

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
        try:
            # TensorFreeSubmodule_comp
            self._comp = tensor_module._comp
        except AttributeError:
            # TensorFreeModule
            tensor = tensor_module()
            frame = list(base_module.irange())
            self._comp = tensor._new_comp(frame)

    def __iter__(self):
        r"""
        Generate the basis elements of ``self``.
        """
        tensor_module = self._fmodule
        base_module = tensor_module.base_module()
        base_module_basis = self._base_module_basis
        for ind in self._comp.non_redundant_index_generator():
            element = tensor_module.element_class(base_module, tensor_module._tensor_type)
            element.set_comp(base_module_basis)[ind] = 1
            yield element

# Todo: Make it a Family
#       symmetrize/antisymmetrize it
#       dual basis
#       add test for dual
# lift/reduce/retract

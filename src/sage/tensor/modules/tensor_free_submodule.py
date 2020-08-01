r"""
Free submodules of tensor products of free modules
"""

#******************************************************************************
#       Copyright (C) 2020 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from .tensor_free_module import TensorFreeModule
from .finite_rank_free_module import FiniteRankFreeModule

class TensorFreeSubmodule_comp(TensorFreeModule):
    r"""
    Class for free submodules of tensor products of free modules
    that are defined by the symmetries of a
    :class:`~sage.tensor.modules.comp.Components` object.

    EXAMPLES::

        sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: Sym2M = TensorFreeSubmodule_comp(M, (2, 0), sym=range(2)); Sym2M
        Free module of type-(2,0) tensors
        with Fully symmetric 2-indices components w.r.t. [0, 1, 2]
        on the Rank-3 free module M over the Integer Ring
    """
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None,
                 sym=None, antisym=None):
        self._fmodule = fmodule
        self._tensor_type = tuple(tensor_type)
        # Create a tensor only because we need a Components object
        tensor = fmodule.tensor(tensor_type,
                                name=name, latex_name=latex_name,
                                sym=sym, antisym=antisym)
        frame = list(fmodule.irange())
        self._comp = tensor._new_comp(frame)
        rank = len(list(self._comp.non_redundant_index_generator()))
        # Skip TensorFreeModule.__init__
        FiniteRankFreeModule.__init__(self, fmodule._ring, rank, name=name,
                                      latex_name=latex_name,
                                      start_index=fmodule._sindex,
                                    output_formatter=fmodule._output_formatter)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 2, name='M')
            sage: Sym2M = TensorFreeSubmodule_comp(M, (2, 0), sym=range(2)); Sym2M
            Free module of type-(2,0) tensors
            with Fully symmetric 2-indices components w.r.t. [0, 1, 2]
            on the Rank-3 free module M over the Integer Ring

        """
        return "Free module of type-({},{}) tensors with {} on the {}".format(
            self._tensor_type[0], self._tensor_type[1], self._comp, self._fmodule)

    def is_submodule(self, other):
        r"""
        Return ``True`` if ``self`` is a submodule of ``other``.

        """
        raise NotImplementedError

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
from sage.sets.disjoint_set import DisjointSet
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

    Canonical injections from submodules are coercions::

        sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: T60M = M.tensor_module(6, 0)
        sage: Sym0123x45M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3), (4, 5)))
        sage: Sym012x345M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2), (3, 4, 5)))
        sage: Sym012345M  = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3, 4, 5)))
        sage: Sym0123x45M.has_coerce_map_from(Sym012345M)
        True
        sage: T60M.has_coerce_map_from(Sym0123x45M)
        True
        sage: t = e[0]^6
        sage: t.parent()
        FIXME
        sage: Sym012345M(t) is t
        FIXME

    TESTS::

        sage: Sym0123x45M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3), (4, 5)))
        sage: TestSuite(Sym0123x45M).run()

    """
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None,
                 sym=None, antisym=None, *, ambient=None, category=None):
        self._fmodule = fmodule
        self._tensor_type = tuple(tensor_type)
        # Create a tensor only because we need a Components object
        tensor = fmodule.tensor(tensor_type,
                                name=name, latex_name=latex_name,
                                sym=sym, antisym=antisym)
        frame = list(fmodule.irange())
        self._comp = tensor._new_comp(frame)
        rank = len(list(self._comp.non_redundant_index_generator()))
        category = fmodule.category().TensorProducts().FiniteDimensional().Subobjects().or_subcategory(category)
        # Skip TensorFreeModule.__init__
        FiniteRankFreeModule.__init__(self, fmodule._ring, rank, name=name,
                                      latex_name=latex_name,
                                      start_index=fmodule._sindex,
                                      output_formatter=fmodule._output_formatter,
                                      ambient=ambient, category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym2M = TensorFreeSubmodule_comp(M, (2, 0), sym=range(2)); Sym2M
            Free module of type-(2,0) tensors
            with Fully symmetric 2-indices components w.r.t. [0, 1, 2]
            on the Rank-3 free module M over the Integer Ring

        """
        return "Free module of type-({},{}) tensors with {} on the {}".format(
            self._tensor_type[0], self._tensor_type[1], self._comp, self._fmodule)

    def _is_symmetry_coarsening_of(self, coarser_comp, finer_comp):
        self_tensor_type = self.tensor_type()

        def sym_antisym(comp):
            if isinstance(comp, tuple):
                sym, antisym = tuple
                if sym is None:
                    sym = []
                if antisym is None:
                    antisym = []
                return sym, antisym
            # Similar code is in Component.contract, should refactor.
            try:
                return comp._sym, comp._antisym
            except AttributeError:
                return [], []

        def is_coarsening_of(self_sym_list, other_sym_list):
            # Use the union-find data structure
            S = DisjointSet(self_tensor_type[0] + self_tensor_type[1])
            for index_set in self_sym_list:
                i = index_set[0]
                for j in index_set[1:]:
                    S.union(i, j)
            for index_set in other_sym_list:
                i = S.find(index_set[0])
                for j in index_set[1:]:
                    if S.find(j) != i:
                        return False
            return True

        finer_sym, finer_antisym = sym_antisym(finer_comp)
        if not finer_sym and not finer_antisym:
            return True
        coarser_sym, coarser_antisym = sym_antisym(coarser_comp)
        if not is_coarsening_of(coarser_sym, finer_sym):
            return False
        if not is_coarsening_of(coarser_antisym, finer_antisym):
            return False
        return True

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None, sym=None, antisym=None):
        if sym is not None or antisym is not None:
            # Refuse to create a tensor with finer symmetries
            # than those defining the subspace
            if not self._is_symmetry_coarsening_of((sym, antisym), self._comp):
                raise ValueError("cannot create a tensor with symmetries {} as an element of {}".
                                 format((sym, antisym), self))
        try:
            comp_parent = comp.parent()
        except AttributeError:
            comp_parent = None
        if comp_parent == self.ambient_module():
            resu = comp
            # comp is already a tensor.  If its declared symmetries are coarser
            # than the symmetries defining self, we can use it directly.
            if self._is_symmetry_coarsening_of(resu, self._comp):
                return resu
        if sym is None:
            sym = self._comp._sym
        if antisym is None:
            sym = self._comp._antisym
        resu = super()._element_constructor_(comp=comp, basis=basis, name=name,
                                             latex_name=latex_name,
                                             sym=sym, antisym=antisym)
        return resu

    def is_submodule(self, other):
        r"""
        Return ``True`` if ``self`` is a submodule of ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: T60M = M.tensor_module(6, 0)
            sage: Sym0123x45M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3), (4, 5)))
            sage: Sym012x345M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2), (3, 4, 5)))
            sage: Sym012345M  = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3, 4, 5)))
            sage: Sym012345M.is_submodule(Sym012345M)
            True
            sage: Sym012345M.is_submodule(Sym0123x45M)
            True
            sage: Sym0123x45M.is_submodule(Sym012345M)
            False
            sage: Sym012x345M.is_submodule(Sym0123x45M)
            False
            sage: all(S.is_submodule(T60M) for S in (Sym0123x45M, Sym012x345M, Sym012345M))
            True

        """
        if super().is_submodule(other):
            return True
        self_base_module = self.base_module()
        self_tensor_type = self.tensor_type()
        try:
            other_base_module = other.base_module()
            other_tensor_type = other.tensor_type()
        except AttributeError:
            return False
        if self_base_module != other_base_module:
            return False
        if self_tensor_type != other_tensor_type:
            return False

        try:
            other_comp = other._comp
        except AttributeError:
            # other is full tensor module (no symmetry)
            return True
        return self._is_symmetry_coarsening_of(self._comp, other_comp)

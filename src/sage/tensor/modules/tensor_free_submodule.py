r"""
Free submodules of tensor products of free modules
"""

#******************************************************************************
#       Copyright (C) 2020-2022 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

import itertools

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_set import DisjointSet
from sage.typeset.unicode_characters import unicode_otimes

from .comp import CompFullySym, CompFullyAntiSym, CompWithSym
from .tensor_free_module import TensorFreeModule
from .finite_rank_free_module import FiniteRankFreeModule_abstract


class TensorFreeSubmodule_comp(TensorFreeModule):
    r"""
    Class for free submodules of tensor products of free modules
    that are defined by the symmetries of a
    :class:`~sage.tensor.modules.comp.Components` object.

    EXAMPLES::

        sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: T60M = M.tensor_module(6, 0); T60M
        Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring
        sage: T60M._name
        'T^(6, 0)(M)'
        sage: latex(T60M)
        T^{(6, 0)}\left(M\right)
        sage: Sym0123x45M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3), (4, 5))); Sym0123x45M
        Free module of type-(6,0) tensors with 6-indices components w.r.t. (0, 1, 2),
         with symmetry on the index positions (0, 1, 2, 3),
         with symmetry on the index positions (4, 5)
         on the Rank-3 free module M over the Integer Ring
        sage: Sym0123x45M._name
        'Sym^{0,1,2,3}(M)⊗Sym^{4,5}(M)'
        sage: latex(Sym0123x45M)
        \mathrm{Sym}^{\{0,1,2,3\}}\left(M\right) \otimes \mathrm{Sym}^{\{4,5\}}\left(M\right)
        sage: Sym012x345M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2), (3, 4, 5))); Sym012x345M
        Free module of type-(6,0) tensors with 6-indices components w.r.t. (0, 1, 2),
         with symmetry on the index positions (0, 1, 2),
         with symmetry on the index positions (3, 4, 5)
         on the Rank-3 free module M over the Integer Ring
        sage: Sym012x345M._name
        'Sym^{0,1,2}(M)⊗Sym^{3,4,5}(M)'
        sage: latex(Sym012x345M)
        \mathrm{Sym}^{\{0,1,2\}}\left(M\right) \otimes \mathrm{Sym}^{\{3,4,5\}}\left(M\right)
        sage: Sym012345M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3, 4, 5))); Sym012345M
        Free module of type-(6,0) tensors
         with Fully symmetric 6-indices components w.r.t. (0, 1, 2)
         on the Rank-3 free module M over the Integer Ring
        sage: Sym012345M._name
        'Sym^6(M)'
        sage: latex(Sym012345M)
        \mathrm{Sym}^6\left(M\right)

    Canonical injections from submodules are coercions::

        sage: Sym0123x45M.has_coerce_map_from(Sym012345M)
        True
        sage: T60M.has_coerce_map_from(Sym0123x45M)
        True
        sage: t = e[0] * e[0] * e[0] * e[0] * e[0] * e[0]
        sage: t.parent()
        Free module of type-(6,0) tensors on the Rank-3 free module M over the Integer Ring
        sage: Sym012345M(t) is t
        False

    TESTS::

        sage: Sym0123x45M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3), (4, 5)))
        sage: TestSuite(Sym0123x45M).run()
        Traceback (most recent call last):
        ...
        The following tests failed: _test_not_implemented_methods, _test_zero

    """
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None,
                 sym=None, antisym=None, *, category=None, ambient=None):
        self._fmodule = fmodule
        self._tensor_type = tuple(tensor_type)
        if ambient is None:
            ambient = fmodule.tensor_module(*tensor_type)
        self._ambient_module = ambient
        self._sym = sym
        self._antisym = antisym
        basis_comp = self._basis_comp()
        rank = len(list(basis_comp.non_redundant_index_generator()))
        # TODO: Good defaults for name, latex_name for more cases
        if name is None and fmodule._name is not None:
            if isinstance(basis_comp, CompFullySym):
                sym = [tuple(range(tensor_type[0] + tensor_type[1]))]
                antisym = []
            elif isinstance(basis_comp, CompFullyAntiSym):
                sym = []
                antisym = [tuple(range(tensor_type[0] + tensor_type[1]))]
            elif isinstance(basis_comp, CompWithSym):
                sym = basis_comp._sym
                antisym = basis_comp._antisym
            else:
                assert False, "full tensor module"

            def power_name(op, s, latex=False):
                if s[0] < tensor_type[0]:
                    assert all(i < tensor_type[0] for i in s)
                    base = fmodule
                    full = tensor_type[0]
                else:
                    base = fmodule.dual()
                    full = tensor_type[1]
                if len(s) == full:
                    superscript = str(full)
                else:
                    superscript = ','.join(str(i) for i in s)
                    if latex:
                        superscript = r'\{' + superscript + r'\}'
                    else:
                        superscript = '{' + superscript  + '}'
                if latex:
                    if len(superscript) != 1:
                        superscript = '{' + superscript + '}'
                    return r'\mathrm{' + op + '}^' + superscript + \
                        r'\left(' + base._latex_name + r'\right)'
                else:
                    return op + '^' + superscript + '(' + base._name + ')'

            name = unicode_otimes.join(itertools.chain(
                (power_name('Sym', s, latex=False) for s in sym),
                (power_name('ASym', s, latex=False) for s in antisym)))
            latex_name = r' \otimes '.join(itertools.chain(
                (power_name('Sym', s, latex=True) for s in sym),
                (power_name('ASym', s, latex=True) for s in antisym)))

        category = fmodule.category().TensorProducts().FiniteDimensional().Subobjects().or_subcategory(category)
        # Skip TensorFreeModule.__init__
        FiniteRankFreeModule_abstract.__init__(self, fmodule._ring, rank, name=name,
                                               latex_name=latex_name,
                                               category=category, ambient=ambient)

    @cached_method
    def _basis_comp(self):
        frame = tuple(self.base_module().irange())
        # Need to call _element_constructor_ explicitly, or the passed arguments are dropped
        tensor = self.ambient()._element_constructor_(sym=self._sym, antisym=self._antisym)
        return tensor._new_comp(frame)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym2M = TensorFreeSubmodule_comp(M, (2, 0), sym=range(2)); Sym2M
            Free module of type-(2,0) tensors
            with Fully symmetric 2-indices components w.r.t. (0, 1, 2)
            on the Rank-3 free module M over the Integer Ring

        """
        return "Free module of type-({},{}) tensors with {} on the {}".format(
            self._tensor_type[0], self._tensor_type[1], self._basis_comp(), self._fmodule)

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
            if not self._is_symmetry_coarsening_of((sym, antisym), self._basis_comp()):
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
            if self._is_symmetry_coarsening_of(resu, self._basis_comp()):
                return resu
        if sym is None:
            sym = self._basis_comp()._sym
        if antisym is None:
            antisym = self._basis_comp()._antisym
        resu = super()._element_constructor_(comp=comp,
                                             basis=basis, name=name,
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

        other_comp = other._basis_comp()
        return self._is_symmetry_coarsening_of(self._basis_comp(), other_comp)

    @lazy_attribute
    def lift(self):
        r"""
        The lift (embedding) map from ``self`` to the ambient space.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_free_submodule import TensorFreeSubmodule_comp
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: Sym0123x45M = TensorFreeSubmodule_comp(M, (6, 0), sym=((0, 1, 2, 3), (4, 5)))
            sage: Sym0123x45M.lift
            Generic morphism:
              From: Free module of type-(6,0) tensors
                    with 6-indices components w.r.t. (0, 1, 2),
                    with symmetry on the index positions (0, 1, 2, 3),
                    with symmetry on the index positions (4, 5)
                    on the Rank-3 free module M over the Integer Ring
              To:   Free module of type-(6,0) tensors
                    on the Rank-3 free module M over the Integer Ring
         """
        return self.module_morphism(function=lambda x: x, codomain=self.ambient())

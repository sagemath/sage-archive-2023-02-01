r"""
Cell Modules
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.all import ModulesWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.dict_addition import dict_addition, dict_linear_combination
from sage.sets.family import Family

class CellModule(CombinatorialFreeModule):
    r"""
    A cell module.

    Let `R` be a commutative ring. Let `A` be a cellular `R`-algebra
    with cell datum `(\Lambda, i, M, C)`. A *cell module* `W(\lambda)`
    is the `R`-module given by `R\{C_s \mid s \in M(\lambda)\}` with
    an action of `a \in A` given by `a C_s = \sum_{u \in M(\lambda}
    r_a(u, s) C_u`, where `r_a(u, s)` is the same as those
    given by the cellular condition:

    .. MATH::

        a C^{\lambda}_{st} = \sum_{u \ in M(\lambda) r_a(u, s)
        C_{ut}^{\lambda} +
        \sum_{\substack{\mu < \lambda \\ xy \in M(\mu)}} R C^{\mu}_{xy}.

    INPUT:

    - ``A`` -- a cellular algebra
    - ``la`` -- an element of the cellular poset of ``A``
    """
    @staticmethod
    def __classcall_private__(cls, A, la, **kwds):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W1 = S.cell_module([2,1])
            sage: W2 = S.cell_module([2,1], prefix='W')
            sage: W1 is W2
            True
        """
        la = A.cell_poset()(la)
        kwds['prefix'] = kwds.get('prefix', 'W')
        return super(CellModule, cls).__classcall__(cls, A, la, **kwds)

    def __init__(self, A, la, **kwds):
        r"""
        Initialize ``self``.

        TESTS::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: TestSuite(W).run()
        """
        self._algebra = A
        self._la = la
        cat = ModulesWithBasis(A.base_ring()).FiniteDimensional()
        CombinatorialFreeModule.__init__(self, A.base_ring(),
                                         A.cell_module_indices(la),
                                         category=cat, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: S.cell_module([2,1])
            Cell module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field
        """
        return "Cell module indexed by {} of {}".format(self._la, self._algebra)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: latex(W)
            W_{...}\left(...\right)
        """
        from sage.misc.latex import latex
        return "W_{{{}}}\\left({}\\right)".format(latex(self._algebra), latex(self._la))

    def cellular_algebra(self):
        r"""
        Return the cellular algebra of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: W = S.cell_module([2,1])
            sage: W.cellular_algebra() is S.cellular_basis()
            True
            sage: S.has_coerce_map_from(W.cellular_algebra())
            True
        """
        return self._algebra

    @cached_method
    def _action_basis(self, x, s):
        """
        Return the action of the basis element of the cellular algebra
        indexed by ``x`` on the basis element of ``self`` indexed by ``s``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: M = S.cell_module([2,1])
            sage: s = StandardTableaux([2,1])([[1,2],[3]])
            sage: la = Partition([2,1])
            sage: M._action_basis((la, s, s), s)
            W[[1, 2], [3]]
        """
        t = self.basis().keys()[0]
        B = self._algebra.basis()
        elt = B[x] * B[(self._la, s, t)]
        mc = elt.monomial_coefficients(copy=False)
        d = {k[1]: mc[k] for k in mc if k[0] == self._la}
        return self._from_dict(d, remove_zeros=False)

    @cached_method
    def _bilinear_form_on_basis(self, s, t):
        """
        Return the result of the bilinear form on the basis elements
        indexed by ``s`` and ``t``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module([2,1])
            sage: K = M.basis().keys()
            sage: matrix([[M._bilinear_form_on_basis(s, t) for t in K] for s in K])
            [1 0]
            [0 1]
        """
        B = self._algebra.basis()
        elt = B[(self._la, s, s)] * B[(self._la, t, t)]
        return elt[(self._la, s, t)]

    def bilinear_form(self, x, y):
        r"""
        Return the bilinear form on ``x`` and ``y``.

        The cell module `W(\lambda)` has a canonical bilinear form
        `\Phi_{\lambda} : W(\lambda) \times W(\lambda) \to W(\lambda)`
        given by

        .. MATH::

            C_{ss}^{\lambda} C_{tt}^{\lambda} = \Phi_{\lambda}(C_s, C_t)
            C_{st}^{\lambda} +
            \sum_{\substack{\mu < \lambda \\ xy \in M(\mu)}} R C^{\mu}_{xy}.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module([2,1])
            sage: elt = M.an_element(); elt
            2*W[[1, 2], [3]] + 2*W[[1, 3], [2]]
            sage: M.bilinear_form(elt, elt)
            8
        """
        R = self.base_ring()
        return R.sum(self._bilinear_form_on_basis(s, t) * cx * cy
                     for s, cx in x for t, cy in y)

    def bilinear_form_matrix(self, ordering=None):
        """
        Return the matrix corresponding to the bilinear form
        of ``self``.

        INPUT:

        - ``ordering`` -- (optional) an ordering of the indices

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module([2,1])
            sage: M.bilinear_form_matrix()
            [1 0]
            [0 1]
        """
        if ordering is None:
            ordering = sorted(self.basis().keys())
        else:
            sordering = set(ordering)
            if sordering != set(self.basis().keys()) or len(sordering) != len(ordering):
                raise ValueError("not an ordering of the basis indices")
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(self.base_ring(), len(ordering))
        return MS([[self._bilinear_form_on_basis(s, t) for t in ordering]
                   for s in ordering])

    @cached_method
    def nonzero_bilinear_form(self):
        """
        Return ``True`` if the bilinear form of ``self`` is non-zero.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module([2,1])
            sage: M.nonzero_bilinear_form()
            True
        """
        C = list(self.basis().keys())
        # Since the bilinear form is symmetric, it is sufficient
        #   to check on the upper triangular part
        return any(self._bilinear_form_on_basis(s, t)
                   for i,s in enumerate(C) for t in C[i:])

    @cached_method
    def radical_basis(self):
        """
        Return a basis of the radical of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module([2,1])
            sage: M.radical_basis()
            ()
        """
        mat = self.bilinear_form_matrix(self.get_order())
        return tuple([self.from_vector(v) for v in mat.left_kernel().basis()])

    def radical(self):
        """
        Return the radical of ``self``.

        Let `W(\lambda)` denote a cell module. The *radical* of `W(\lambda)`
        is defined as

        .. MATH::

            \operatorname{rad}(\lambda) := \{x \in W(\lambda) \mid
            \Phi_{\lambda}(x, y),

        and note that it is a submodule of `W(\lambda)`.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module([2,1])
            sage: R = M.radical(); R
            Radical of Cell module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field
            sage: R.basis()
            Finite family {}
        """
        from sage.modules.with_basis.subquotient import SubmoduleWithBasis
        radical = self.submodule(self.radical_basis(),
                                category=self.category().Subobjects(),
                                already_echelonized=True)
        radical.rename("Radical of {}".format(self))
        return radical

    def simple_module(self):
        r"""
        Return the corresponding simple module of ``self``.

        Let `W(\lambda)` denote a cell module. The simple module `L(\lambda)`
        is defined as `W(\lambda) / \operatorname{rad}(\lambda)`,
        where `\operatorname{rad}(\lambda)` is the radical of the
        bilinear form `\Phi_{\lambda}`.

        .. SEEALSO::

            :meth:`radical`

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module([2,1])
            sage: L = M.simple_module(); L
            Simple module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field
            sage: L.has_coerce_map_from(M)
            True
        """
        if not self.nonzero_bilinear_form():
            raise ValueError("bilinear form is 0; no corresponding simple module")
        L = self.quotient_module(self.radical(), check=False,
                                 already_echelonized=True)
        L.rename("Simple module indexed by {} of {}".format(self._la, self._algebra))
        # Add the coercion
        self.module_morphism(L.retract, codomain=L).register_as_coercion()
        return L

    class Element(CombinatorialFreeModule.Element):

        def _acted_upon_(self, scalar, self_on_left=False):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: S = SymmetricGroupAlgebra(QQ, 3)
                sage: C = S.cellular_basis()
                sage: M = S.cell_module([2,1])
                sage: elt = M.an_element(); elt
                2*W[[1, 2], [3]] + 2*W[[1, 3], [2]]
                sage: sc = C.an_element()
                sage: sc
                3*C([2, 1], [[1, 3], [2]], [[1, 2], [3]])
                 + 2*C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
                 + 2*C([3], [[1, 2, 3]], [[1, 2, 3]])
                sage: sc * elt
                10*W[[1, 3], [2]]

                sage: 1/2 * elt
                W[[1, 2], [3]] + W[[1, 3], [2]]
            """
            P = self.parent()
            if scalar in P._algebra:
                if self_on_left:
                    raise NotImplementedError
                    #scalar = scalar.cellular_involution()
                mc = self._monomial_coefficients
                scalar_mc = scalar.monomial_coefficients(copy=False)
                D = dict_linear_combination([(P._action_basis(x, k)._monomial_coefficients,
                                              scalar_mc[x] * mc[k])
                                             for k in mc for x in scalar_mc],
                                            factor_on_left=False)

                return P._from_dict(D, remove_zeros=False)

            return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)

        # For backward compatibility
        _lmul_ = _acted_upon_
        _rmul_ = _acted_upon_


"""
Type-(1,1) tensors on free modules

Three derived classes of
:class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor` are devoted
to type-(1,1) tensors:

* :class:`FreeModuleEndomorphismTensor` for endomorphisms viewed as type-(1,1)
  tensors

  * :class:`FreeModuleAutomorphismTensor` for invertible endomorphisms viewed
    as type-(1,1) tensors

    * :class:`FreeModuleIdentityTensor` for the identity map viewed as a 
      type-(1,1) tensor


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

.. TODO::

    Suppress :class:`FreeModuleEndomorphismTensor` ? (since the
    coercion of type-(1,1) tensors to free module endomorphisms is implemented
    now) This would leave only :class:`FreeModuleAutomorphismTensor` and
    :class:`FreeModuleIdentityTensor`.

"""
#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_tensor import FreeModuleTensor

class FreeModuleEndomorphismTensor(FreeModuleTensor):
    r"""
    Endomorphism (considered as a type-`(1,1)` tensor) on a free module.

    INPUT:

    - ``fmodule`` -- free module `M` over a commutative ring `R`
      (must be an instance of :class:`FiniteRankFreeModule`)
    - ``name`` -- (default: ``None``) name given to the endomorphism
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      endomorphism; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Endomorphism tensor on a rank-3 module::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: t = M.endomorphism_tensor('T') ; t
        Endomorphism tensor T on the Rank-3 free module M over the Integer Ring
        sage: t.parent()
        Free module of type-(1,1) tensors on the
         Rank-3 free module M over the Integer Ring
        sage: t.tensor_type()
        (1, 1)
        sage: t.tensor_rank()
        2

    The method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.tensor`
    with the argument ``(1,1)`` can be used as well to create such a tensor::

        sage: t = M.tensor((1,1), name='T') ; t
        Endomorphism tensor T on the Rank-3 free module M over the Integer Ring

    Components of the endomorphism with respect to a given basis::

        sage: e = M.basis('e') ; e
        Basis (e_0,e_1,e_2) on the Rank-3 free module M over the Integer Ring
        sage: t[:] = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: t[:]
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: t.view()
        T = e_0*e^0 + 2 e_0*e^1 + 3 e_0*e^2 + 4 e_1*e^0 + 5 e_1*e^1
            + 6 e_1*e^2 + 7 e_2*e^0 + 8 e_2*e^1 + 9 e_2*e^2

    The matrix of components w.r.t. to a given basis::

        sage: m = matrix(t.components(e)) ; m
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: m.parent()
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

    The endomorphism acting on a module element::

        sage: v = M([1,2,3], basis=e, name='v') ; v
        Element v of the Rank-3 free module M over the Integer Ring
        sage: w = t(v) ; w
        Element T(v) of the Rank-3 free module M over the Integer Ring
        sage: w[:]
        [14, 32, 50]
        sage: for i in M.irange():   # Check:
        ....:     print sum( t[i,j]*v[j] for j in M.irange() ),
        14 32 50

    """
    def __init__(self, fmodule, name=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_tensor_spec import FreeModuleEndomorphismTensor
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: E = FreeModuleEndomorphismTensor(M, name='a')
            sage: E[e,0,1] = -3
            sage: TestSuite(E).run(skip="_test_category") # see below

        In the above test suite, _test_category fails because E is not an
        instance of E.parent().category().element_class.

        """
        FreeModuleTensor.__init__(self, fmodule, (1,1), name=name,
                                  latex_name=latex_name)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: M.endomorphism_tensor()
            Endomorphism tensor on the Rank-3 free module M over the Integer Ring
            sage: M.endomorphism_tensor(name='a')
            Endomorphism tensor a on the Rank-3 free module M over the Integer Ring

        """
        description = "Endomorphism tensor "
        if self._name is not None:
            description += self._name + " "
        description += "on the " + str(self._fmodule)
        return description

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.endomorphism_tensor(name='a')
            sage: a._new_instance()
            Endomorphism tensor on the Rank-3 free module M over the Integer Ring

        """
        return self.__class__(self._fmodule)

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleTensor.__call__` to allow for a single
        argument (module element).

        EXAMPLES:

        Call with a single argument --> return a module element::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.endomorphism_tensor(name='a')
            sage: e = M.basis('e')
            sage: a[0,1], a[1,1], a[2,1] = 2, 4, -5
            sage: v = M([2,1,4], name='v')
            sage: s = a.__call__(v) ; s
            Element a(v) of the Rank-3 free module M over the Integer Ring
            sage: s.view()
            a(v) = 2 e_0 + 4 e_1 - 5 e_2
            sage: s == a(v)
            True
            sage: s == a.contract(v)
            True

        Call with two arguments (:class:`FreeModuleTensor` behaviour)
        --> return a scalar::

            sage: b = M.linear_form(name='b')
            sage: b[:] = 7, 0, 2
            sage: a.__call__(b,v)
            4
            sage: a(b,v) == a.__call__(b,v)
            True
            sage: a(b,v) == s(b)
            True

        """
        from free_module_tensor import FiniteRankFreeModuleElement
        if len(arg) > 1:
            # the endomorphism acting as a type-(1,1) tensor on a pair
            # (linear form, module element), returning a scalar:
            return FreeModuleTensor.__call__(self, *arg)
        # the endomorphism acting as such, on a module element, returning a
        # module element:
        vector = arg[0]
        if not isinstance(vector, FiniteRankFreeModuleElement):
            raise TypeError("the argument must be an element of a free module")
        basis = self.common_basis(vector)
        t = self._components[basis]
        v = vector._components[basis]
        fmodule = self._fmodule
        result = vector._new_instance()
        for i in fmodule.irange():
            res = 0
            for j in fmodule.irange():
                res += t[[i,j]]*v[[j]]
            result.set_comp(basis)[i] = res
        # Name of the output:
        result._name = None
        if self._name is not None and vector._name is not None:
            result._name = self._name + "(" + vector._name + ")"
        # LaTeX symbol for the output:
        result._latex_name = None
        if self._latex_name is not None and vector._latex_name is not None:
            result._latex_name = self._latex_name + r"\left(" + \
                              vector._latex_name + r"\right)"
        return result

#******************************************************************************

class FreeModuleAutomorphismTensor(FreeModuleEndomorphismTensor):
    r"""
    Automorphism (considered as a type-`(1,1)` tensor) on a free module.

    INPUT:

    - ``fmodule`` -- free module `M` over a commutative ring `R`
      (must be an instance of :class:`FiniteRankFreeModule`)
    - ``name`` -- (default: ``None``) name given to the automorphism
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      automorphism; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Automorphism tensor on a rank-2 free module (vector space) on `\QQ`::

        sage: M = FiniteRankFreeModule(QQ, 2, name='M')
        sage: a = M.automorphism_tensor('A') ; a
        Automorphism tensor A on the Rank-2 free module M over the Rational Field

    Automorphisms are tensors of type `(1,1)`::

        sage: a.parent()
        Free module of type-(1,1) tensors on the
         Rank-2 free module M over the Rational Field
        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2

    Setting the components in a basis::

        sage: e = M.basis('e') ; e
        Basis (e_0,e_1) on the Rank-2 free module M over the Rational Field
        sage: a[:] = [[1, 2], [-1, 3]]
        sage: a[:]
        [ 1  2]
        [-1  3]
        sage: a.view(basis=e)
        A = e_0*e^0 + 2 e_0*e^1 - e_1*e^0 + 3 e_1*e^1

    The inverse automorphism is obtained via the method :meth:`inverse`::

        sage: b = a.inverse() ; b
        Automorphism tensor A^(-1) on the Rank-2 free module M over the Rational Field
        sage: b.view(basis=e)
        A^(-1) = 3/5 e_0*e^0 - 2/5 e_0*e^1 + 1/5 e_1*e^0 + 1/5 e_1*e^1
        sage: b[:]
        [ 3/5 -2/5]
        [ 1/5  1/5]
        sage: a[:] * b[:]  # check that b is indeed the inverse of a
        [1 0]
        [0 1]

    """
    def __init__(self, fmodule, name=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_tensor_spec import FreeModuleAutomorphismTensor
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = FreeModuleAutomorphismTensor(M, name='a')
            sage: a[e,:] = [[1,0,1],[0,2,0],[0,0,-3]]
            sage: TestSuite(a).run(skip="_test_category") # see below

        In the above test suite, _test_category fails because a is not an
        instance of a.parent().category().element_class.

        """
        FreeModuleEndomorphismTensor.__init__(self, fmodule, name=name,
                                        latex_name=latex_name)
        self._inverse = None    # inverse automorphism not set yet

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: M.automorphism_tensor()
            Automorphism tensor on the Rank-3 free module M over the Rational Field
            sage: M.automorphism_tensor(name='a')
            Automorphism tensor a on the Rank-3 free module M over the Rational Field

        """
        description = "Automorphism tensor "
        if self._name is not None:
            description += self._name + " "
        description += "on the " + str(self._fmodule)
        return description

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: a = M.automorphism_tensor(name='a')
            sage: e = M.basis('e')
            sage: a[e,:] = [[1,0,-1], [0,3,0], [0,0,2]]
            sage: b = a.inverse()
            sage: a._inverse
            Automorphism tensor a^(-1) on the Rank-3 free module M over the Rational Field
            sage: a._del_derived()
            sage: a._inverse  # has been reset to None

        """
        # First delete the derived quantities pertaining to the mother class:
        FreeModuleEndomorphismTensor._del_derived(self)
        # Then deletes the inverse automorphism:
        self._inverse = None

    def inverse(self):
        r"""
        Return the inverse automorphism.

        OUTPUT:

        - instance of :class:`FreeModuleAutomorphismTensor` representing the
          automorphism that is the inverse of ``self``.

        EXAMPLES:

        Inverse of an automorphism on a rank-3 free module::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: a = M.automorphism_tensor('A')
            sage: e = M.basis('e')
            sage: a[:] = [[1,0,-1], [0,3,0], [0,0,2]]
            sage: b = a.inverse() ; b
            Automorphism tensor A^(-1) on the Rank-3 free module M over the Rational Field
            sage: b[:]
            [  1   0 1/2]
            [  0 1/3   0]
            [  0   0 1/2]

        We may check that ``b`` is the inverse of ``a`` by performing the
        matrix product of the components in the basis ``e``::

            sage: a[:] * b[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]

        Another check is of course::

            sage: b.inverse() == a
            True

        """
        from sage.matrix.constructor import matrix
        from comp import Components
        if self._inverse is None:
            if self._name is None:
                inv_name = None
            else:
                inv_name = self._name  + '^(-1)'
            if self._latex_name is None:
                inv_latex_name = None
            else:
                inv_latex_name = self._latex_name + r'^{-1}'
            fmodule = self._fmodule
            si = fmodule._sindex
            nsi = fmodule._rank + si
            self._inverse = self.__class__(fmodule, inv_name, inv_latex_name)
            for basis in self._components:
                try:
                    mat_self = matrix(
                              [[self.comp(basis)[[i, j]]
                              for j in range(si, nsi)] for i in range(si, nsi)])
                except (KeyError, ValueError):
                    continue
                mat_inv = mat_self.inverse()
                cinv = Components(fmodule._ring, basis, 2, start_index=si,
                                  output_formatter=fmodule._output_formatter)
                for i in range(si, nsi):
                    for j in range(si, nsi):
                        cinv[i, j] = mat_inv[i-si,j-si]
                self._inverse._components[basis] = cinv
        return self._inverse


#******************************************************************************

class FreeModuleIdentityTensor(FreeModuleAutomorphismTensor):
    r"""
    Identity map (considered as a type-(1,1) tensor) on a free module.

    INPUT:

    - ``fmodule`` -- free module `M` over a commutative ring `R`
      (must be an instance of :class:`FiniteRankFreeModule`)
    - ``name`` -- (default: 'Id') name given to the identity tensor.
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the identity
      tensor; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Identity tensor on a rank-3 free module::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: a = M.identity_tensor() ; a
        Identity tensor on the Rank-3 free module M over the Integer Ring

    The LaTeX symbol is set by default to `\mathrm{Id}`, but can be changed::

        sage: latex(a)
        \mathrm{Id}
        sage: a = M.identity_tensor(latex_name=r'\mathrm{1}')
        sage: latex(a)
        \mathrm{1}

    The identity is a tensor of type `(1,1)` on the free module::

        sage: a.parent()
        Free module of type-(1,1) tensors on the
         Rank-3 free module M over the Integer Ring
        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2

    Its components are Kronecker deltas in any basis::

        sage: a[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: a.comp() # components in the module's default basis (e)
        Kronecker delta of size 3x3
        sage: a.view()
        Id = e_0*e^0 + e_1*e^1 + e_2*e^2
        sage: f = M.basis('f')
        sage: a.comp(basis=f)
        Kronecker delta of size 3x3
        sage: a.comp(f)[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]

    The components can be read, but cannot be set::

        sage: a[1,1]
        1
        sage: a[1,1] = 2
        Traceback (most recent call last):
        ...
        TypeError: the components of the identity map cannot be changed

    The identity tensor acting on a module element::

        sage: v = M([2,-3,1], basis=e, name='v')
        sage: v.view()
        v = 2 e_0 - 3 e_1 + e_2
        sage: u = a(v) ; u
        Element v of the Rank-3 free module M over the Integer Ring
        sage: u is v
        True

    The identity tensor acting as a type-`(1,1)` tensor on a pair (linear form,
    module element)::

        sage: w = M.tensor((0,1), name='w') ; w
        Linear form w on the Rank-3 free module M over the Integer Ring
        sage: w[:] = [0, 3, 2]
        sage: s = a(w,v) ; s
        -7
        sage: s == w(v)
        True

    The identity tensor is its own inverse::

        sage: a.inverse() == a
        True
        sage: a.inverse() is a
        True

    """
    def __init__(self, fmodule, name='Id', latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_tensor_spec import FreeModuleIdentityTensor
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: Id = FreeModuleIdentityTensor(M)
            sage: TestSuite(Id).run(skip="_test_category") # see below

        In the above test suite, _test_category fails because Id is not an
        instance of Id.parent().category().element_class.

        """
        if latex_name is None and name == 'Id':
            latex_name = r'\mathrm{Id}'
        FreeModuleAutomorphismTensor.__init__(self, fmodule, name=name,
                                              latex_name=latex_name)
        self._inverse = self    # the identity is its own inverse
        self.comp() # Initializing the components in the module's default basis

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M.identity_tensor()
            Identity tensor on the Rank-3 free module M over the Integer Ring

        """
        description = "Identity tensor "
        if self._name != 'Id':
            description += self._name + " "
        description += "on the " + str(self._fmodule)
        return description

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: id = M.identity_tensor()
            sage: id._del_derived()

        """
        # FreeModuleAutomorphismTensor._del_derived is bypassed:
        FreeModuleEndomorphismTensor._del_derived(self)

    def _new_comp(self, basis):
        r"""
        Create some components in the given basis.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: id = M.identity_tensor()
            sage: id._new_comp(e)
            Kronecker delta of size 3x3
            sage: type(id._new_comp(e))
            <class 'sage.tensor.modules.comp.KroneckerDelta'>

        """
        from comp import KroneckerDelta
        fmodule = self._fmodule  # the base free module
        return KroneckerDelta(fmodule._ring, basis, start_index=fmodule._sindex,
                              output_formatter=fmodule._output_formatter)

    def components(self, basis=None, from_basis=None):
        r"""
        Return the components in a given basis as a Kronecker delta.

        INPUT:

        - ``basis`` -- (default: ``None``) module basis in which the components
          are required; if none is provided, the components are assumed to
          refer to the module's default basis
        - ``from_basis`` -- (default: ``None``) unused (present just for
          ensuring compatibility with ``FreeModuleTensor.comp`` calling list)

        OUTPUT:

        - components in the basis ``basis``, as an instance of the
          class :class:`~sage.tensor.modules.comp.KroneckerDelta`

        EXAMPLES:

        Components of the identity map on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_tensor()
            sage: a.components(basis=e)
            Kronecker delta of size 3x3

        For the module's default basis, the argument ``basis`` can be omitted::

            sage: a.components() is a.components(basis=e)
            True
            sage: a.components()[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]

        A shortcut is ``a.comp()``::

            sage: a.comp() is a.components()
            True

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis not in self._components:
            self._components[basis] = self._new_comp(basis)
        return self._components[basis]

    comp = components

    def set_comp(self, basis=None):
        r"""
        Redefinition of the generic tensor method
        :meth:`FreeModuleTensor.set_comp`: should not be called.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_tensor()
            sage: a.set_comp(e)
            Traceback (most recent call last):
            ...
            TypeError: the components of the identity map cannot be changed

        """
        raise TypeError("the components of the identity map cannot be changed")

    def add_comp(self, basis=None):
        r"""
        Redefinition of the generic tensor method
        :meth:`FreeModuleTensor.add_comp`: should not be called.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_tensor()
            sage: a.add_comp(e)
            Traceback (most recent call last):
            ...
            TypeError: the components of the identity map cannot be changed

        """
        raise TypeError("the components of the identity map cannot be changed")

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleEndomorphismTensor.__call__`.

        EXAMPLES:

        Call with a single argument --> return a module element::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: id = M.identity_tensor()
            sage: v = M([-1,4,3])
            sage: s = id.__call__(v) ; s
            Element of the Rank-3 free module M over the Integer Ring
            sage: s == v
            True
            sage: s == id(v)
            True
            sage: s == id.contract(v)
            True

        Call with two arguments (:class:`FreeModuleTensor` behaviour) -->
        return a scalar::

            sage: b = M.linear_form(name='b')
            sage: b[:] = 7, 0, 2
            sage: id.__call__(b,v)
            -1
            sage: id(b,v) == id.__call__(b,v)
            True
            sage: id(b,v) == b(v)
            True

        """
        from free_module_tensor import FiniteRankFreeModuleElement
        from free_module_alt_form import FreeModuleLinForm
        if len(arg) == 1:
            # the identity map acting as such, on a module element:
            vector = arg[0]
            if not isinstance(vector, FiniteRankFreeModuleElement):
                raise TypeError("the argument must be a module element")
            return vector
            #!# should it be return vector.copy() instead ?
        elif len(arg) == 2:
            # the identity map acting as a type-(1,1) tensor on a pair
            # (1-form, vector), returning a scalar:
            linform = arg[0]
            if not isinstance(linform, FreeModuleLinForm):
                raise TypeError("the first argument must be a linear form")
            vector = arg[1]
            if not isinstance(vector, FiniteRankFreeModuleElement):
                raise TypeError("the second argument must be a module element")
            return linform(vector)
        else:
            raise TypeError("wrong number of arguments")


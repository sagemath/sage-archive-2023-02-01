r"""
Free module bases

The class :class:`FreeModuleBasis` implements bases on a free module `M` of
finite rank over a commutative ring,
while the class :class:`FreeModuleCoBasis` implements the dual bases (i.e. 
bases of the dual module `M^*`). 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version


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

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

class FreeModuleBasis(UniqueRepresentation, SageObject):
    r""" 
    Basis of a free module over a commutative ring `R`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` (must be an instance of 
      :class:`FiniteRankFreeModule`)
    - ``symbol`` -- (string) a letter (of a few letters) to denote a generic 
      element of the basis
    - ``latex_symbol`` -- (string; default: None) symbol to denote a generic 
      element of the basis; if None, the value of ``symbol`` is used. 

    EXAMPLES:
    
    A basis on a rank-3 free module over `\ZZ`::
        
        sage: M0 = FiniteRankFreeModule(ZZ, 3, name='M_0')
        sage: from sage.tensor.modules.free_module_basis import FreeModuleBasis
        sage: e = FreeModuleBasis(M0, 'e') ; e
        basis (e_0,e_1,e_2) on the rank-3 free module M_0 over the Integer Ring

    Instead of importing FreeModuleBasis in the global name space, one can 
    use the module's method :meth:`basis`::
    
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e') ; e
        basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring

    The individual elements constituting the basis are accessed via the 
    square bracket operator::
    
        sage: e[0]
        element e_0 of the rank-3 free module M over the Integer Ring
        sage: e[0] in M
        True

    The LaTeX symbol can be set explicitely, as the second argument of
    :meth:`basis`::
    
        sage: latex(e)
        \left(e_0,e_1,e_2\right)
        sage: eps = M.basis('eps', r'\epsilon') ; eps
        basis (eps_0,eps_1,eps_2) on the rank-3 free module M over the Integer Ring
        sage: latex(eps)
        \left(\epsilon_0,\epsilon_1,\epsilon_2\right)
        
    The individual elements of the basis are labelled according the 
    parameter ``start_index`` provided at the free module construction::
    
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: e = M.basis('e') ; e
        basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring
        sage: e[1]
        element e_1 of the rank-3 free module M over the Integer Ring
        
    """
    def __init__(self, fmodule, symbol, latex_symbol=None):
        r"""
        TESTS::
        
            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: from sage.tensor.modules.free_module_basis import FreeModuleBasis
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: FreeModuleBasis(M, 'e', latex_symbol=r'\epsilon')
            basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
            sage: e = FreeModuleBasis(M, 'e', latex_symbol=r'\epsilon')
            sage: e is FreeModuleBasis(M, 'e', latex_symbol=r'\epsilon') # unique representation
            True
            
        """
        self._fmodule = fmodule
        self._name = "(" + \
          ",".join([symbol + "_" + str(i) for i in fmodule.irange()]) +")"
        if latex_symbol is None:
            latex_symbol = symbol
        self._latex_name = r"\left(" + ",".join([latex_symbol + "_" + str(i) 
                                       for i in fmodule.irange()]) + r"\right)"
        self._symbol = symbol
        self._latex_symbol = latex_symbol
        # The basis is added to the module list of bases
        for other in fmodule._known_bases:
            if symbol == other._symbol:
                raise ValueError("The " + str(other) + " already exist on the " +
                                 str(fmodule))
        fmodule._known_bases.append(self)
        # The individual vectors:
        vl = list()
        for i in fmodule.irange():
            v_name = symbol + "_" + str(i)
            v_symb = latex_symbol + "_" + str(i)
            v = fmodule.element_class(fmodule, name=v_name, latex_name=v_symb)
            for j in fmodule.irange():
                v.set_comp(self)[j] = fmodule._ring.zero()
            v.set_comp(self)[i] = fmodule._ring.one()
            vl.append(v)
        self._vec = tuple(vl)
        # The first defined basis is considered as the default one:
        if fmodule._def_basis is None:
            fmodule._def_basis = self
        # Initialization of the components w.r.t the current basis of the zero 
        # elements of all tensor modules constructed up to now (including the 
        # base module itself, since it is considered as a type-(1,0) tensor 
        # module)        
        for t in fmodule._tensor_modules.itervalues():
            t._zero_element._components[self] = t._zero_element._new_comp(self)
                               # (since new components are initialized to zero)
        # The dual basis:
        self._dual_basis = self._init_dual_basis() 


    ###### Methods to be redefined by derived classes of FreeModuleBasis ######

    def _repr_(self):
        r"""
        String representation of the object.
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._repr_()
            'basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring'
            sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e1 = M1.basis('e')
            sage: e1._repr_()
            'basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring'

        """
        return "basis " + self._name + " on the " + str(self._fmodule)


    def _init_dual_basis(self):
        r""" 
        Construct the basis dual to ``self``.
        
        OUTPUT:
        
        - instance of :class:`FreeModuleCoBasis` representing the dual of
          ``self``
        
        EXAMPLE::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._init_dual_basis()
            dual basis (e^0,e^1,e^2) on the rank-3 free module M over the Integer Ring
                  
        """
        return FreeModuleCoBasis(self, self._symbol, 
                                               latex_symbol=self._latex_symbol)

    def _new_instance(self, symbol, latex_symbol=None):
        r"""
        Construct a new basis on the same module as ``self``. 
        
        INPUT:
        
        - ``symbol`` -- (string) a letter (of a few letters) to denote a
          generic element of the basis
        - ``latex_symbol`` -- (string; default: None) symbol to denote a 
          generic element of the basis; if None, the value of ``symbol`` is 
          used. 

        OUTPUT:
        
        - instance of :class:`FreeModuleBasis`
        
        EXAMPLE::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._new_instance('f')
            basis (f_0,f_1,f_2) on the rank-3 free module M over the Integer Ring
                
        """
        return FreeModuleBasis(self._fmodule, symbol, latex_symbol=latex_symbol)
        
    ###### End of methods to be redefined by derived classes ######


    def dual_basis(self):
        r""" 
        Return the basis dual to ``self``.
        
        OUTPUT:
        
        - instance of :class:`FreeModuleCoBasis` representing the dual of
          ``self``

        EXAMPLES:
        
        Dual basis on a rank-3 free module::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e') ; e
            basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring
            sage: f = e.dual_basis() ; f
            dual basis (e^1,e^2,e^3) on the rank-3 free module M over the Integer Ring
        
        Let us check that the elements of f are tensors of type (0,1) on M::
    
            sage: f[1] in M.tensor_module(0,1)
            True
            sage: f[1]
            linear form e^1 on the rank-3 free module M over the Integer Ring
    
        and that f is indeed the dual of e::
        
            sage: f[1](e[1]), f[1](e[2]), f[1](e[3])
            (1, 0, 0)
            sage: f[2](e[1]), f[2](e[2]), f[2](e[3])
            (0, 1, 0)
            sage: f[3](e[1]), f[3](e[2]), f[3](e[3])
            (0, 0, 1)
        
        """
        return self._dual_basis

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        
        EXAMPLES::

            sage: FiniteRankFreeModule._clear_cache_() # for doctests only
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e._latex_()
            '\\left(e_0,e_1,e_2\\right)'
            sage: latex(e)
            \left(e_0,e_1,e_2\right)
            sage: f = M.basis('eps', latex_symbol=r'\epsilon')
            sage: f._latex_()
            '\\left(\\epsilon_0,\\epsilon_1,\\epsilon_2\\right)'
            sage: latex(f)
            \left(\epsilon_0,\epsilon_1,\epsilon_2\right)
        
        """
        return self._latex_name

    def __hash__(self):
        r"""
        Hash function (since instances of :class:`FreeModuleBasis` are used as
        dictionary keys).
                
        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.__hash__() # random
            140529558663376

        """
        return id(self)

    def __eq__(self, other):
        r"""
        Equality (comparison) operator
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = M.basis('f')
            sage: e.__eq__(f)
            False
            sage: e.__eq__(e)
            True
            sage: e.__eq__(M.basis('e'))
            True
            sage: M2 = FiniteRankFreeModule(ZZ, 2, name='M2')
            sage: e.__eq__(M2.basis('e'))
            False

        """
        return other is self

    def __ne__(self, other):
        r"""
        Non-equality operator.
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = M.basis('f')
            sage: e.__ne__(f)
            True
            sage: e.__ne__(e)
            False

        """
        return not self.__eq__(other)
        
    def __getitem__(self, index):
        r"""
        Returns the basis element corresponding to a given index.
        
        INPUT:
        
        - ``index`` -- the index of the basis element
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.__getitem__(0)
            element e_0 of the rank-3 free module M over the Integer Ring
            sage: e.__getitem__(1)
            element e_1 of the rank-3 free module M over the Integer Ring
            sage: e.__getitem__(2)
            element e_2 of the rank-3 free module M over the Integer Ring
            sage: e[1] is e.__getitem__(1)
            True
            sage: e[1].parent() is M
            True
            sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e1 = M1.basis('e')
            sage: e1.__getitem__(1)
            element e_1 of the rank-3 free module M over the Integer Ring
            sage: e1.__getitem__(2)
            element e_2 of the rank-3 free module M over the Integer Ring
            sage: e1.__getitem__(3)
            element e_3 of the rank-3 free module M over the Integer Ring

        """
        n = self._fmodule._rank
        si = self._fmodule._sindex
        i = index - si
        if i < 0 or i > n-1:
            raise ValueError("Index out of range: " +
                              str(i+si) + " not in [" + str(si) + "," +
                              str(n-1+si) + "]")
        return self._vec[i]

    def __len__(self):
        r"""
        Return the basis length, i.e. the rank of the free module.
        
        NB: the method  __len__()  is required for the basis to act as a 
        "frame" in the class :class:`~sage.tensor.modules.comp.Components`. 
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: e.__len__()
            3
            sage: len(e)
            3

        """
        return self._fmodule._rank

    def new_basis(self, change_of_basis, symbol, latex_symbol=None):
        r"""
        Define a new module basis from the current one. 
        
        The new basis is defined by means of a module automorphism. 
        
        INPUT:
        
        - ``change_of_basis`` -- instance of 
          :class:`~sage.tensor.modules.free_module_tensor_spec.FreeModuleAutomorphism`
          describing the automorphism `P` that relates the current basis 
          `(e_i)` (described by ``self``) to the new basis `(n_i)` according
          to `n_i = P(e_i)`
        - ``symbol`` -- (string) a letter (of a few letters) to denote a 
          generic element of the basis
        - ``latex_symbol`` -- (string; default: None) symbol to denote a 
          generic element of the basis; if None, the value of ``symbol`` is 
          used. 
          
        OUTPUT:
        
        - the new basis `(n_i)`, as an instance of :class:`FreeModuleBasis`
        
        EXAMPLES:
        
        Change of basis on a rank-2 free module::
        
            sage: M = FiniteRankFreeModule(QQ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism()
            sage: a[:] = [[1, 2], [-1, 3]]
            sage: f = e.new_basis(a, 'f') ; f
            basis (f_1,f_2) on the rank-2 free module M over the Rational Field
            sage: f[1].view()
            f_1 = e_1 - e_2
            sage: f[2].view()
            f_2 = 2 e_1 + 3 e_2
            sage: e[1].view(f)
            e_1 = 3/5 f_1 + 1/5 f_2
            sage: e[2].view(f)
            e_2 = -2/5 f_1 + 1/5 f_2

        """
        from free_module_tensor_spec import FreeModuleAutomorphism
        if not isinstance(change_of_basis, FreeModuleAutomorphism):
            raise TypeError("The argument change_of_basis must be some " +
                            "instance of FreeModuleAutomorphism.")
        fmodule = self._fmodule
        # self._new_instance used instead of FreeModuleBasis for a correct
        # construction in case of derived classes:
        the_new_basis = self._new_instance(symbol, latex_symbol=latex_symbol)
        transf = change_of_basis.copy()
        inv_transf = change_of_basis.inverse().copy()
        si = fmodule._sindex
        # Components of the new basis vectors in the old basis: 
        for i in fmodule.irange():
            for j in fmodule.irange():
                the_new_basis._vec[i-si].add_comp(self)[[j]] = \
                                                  transf.comp(self)[[j,i]]
        # Components of the new dual-basis elements in the old dual basis: 
        for i in fmodule.irange():
            for j in fmodule.irange():
                the_new_basis._dual_basis._form[i-si].add_comp(self)[[j]] = \
                                              inv_transf.comp(self)[[i,j]]
        # The components of the transformation and its inverse are the same in 
        # the two bases:
        for i in fmodule.irange():
            for j in fmodule.irange():
                transf.add_comp(the_new_basis)[[i,j]] = transf.comp(self)[[i,j]]
                inv_transf.add_comp(the_new_basis)[[i,j]] = \
                                              inv_transf.comp(self)[[i,j]]
        # Components of the old basis vectors in the new basis: 
        for i in fmodule.irange():
            for j in fmodule.irange():
                self._vec[i-si].add_comp(the_new_basis)[[j]] = \
                                                   inv_transf.comp(self)[[j,i]]
        # Components of the old dual-basis elements in the new cobasis: 
        for i in fmodule.irange():
            for j in fmodule.irange():
                self._dual_basis._form[i-si].add_comp(the_new_basis)[[j]] = \
                                                       transf.comp(self)[[i,j]]
        # The automorphism and its inverse are added to the module's dictionary 
        # of changes of bases:
        fmodule._basis_changes[(self, the_new_basis)] = transf
        fmodule._basis_changes[(the_new_basis, self)] = inv_transf
        #
        return the_new_basis

        
#******************************************************************************

class FreeModuleCoBasis(SageObject):
    r""" 
    Dual basis of a free module over a commutative ring.
    
    INPUT:
    
    - ``basis`` -- basis of a free module `M` of which ``self`` is the dual 
      (must be an instance of :class:`FreeModuleBasis`)
    - ``symbol`` -- a letter (of a few letters) to denote a generic element of
      the cobasis
    - ``latex_symbol`` -- (default: None) symbol to denote a generic element of
      the cobasis; if None, the value of ``symbol`` is used. 

    EXAMPLES:
    
    Dual basis on a rank-3 free module::
    
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
        sage: e = M.basis('e') ; e
        basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring
        sage: from sage.tensor.modules.free_module_basis import FreeModuleCoBasis
        sage: f = FreeModuleCoBasis(e, 'f') ; f
        dual basis (f^1,f^2,f^3) on the rank-3 free module M over the Integer Ring
    
    Let us check that the elements of f are tensors of type (0,1) on M::

        sage: f[1] in M.tensor_module(0,1)
        True
        sage: f[1]
        linear form f^1 on the rank-3 free module M over the Integer Ring

    and that f is indeed the dual of e::
    
        sage: f[1](e[1]), f[1](e[2]), f[1](e[3])
        (1, 0, 0)
        sage: f[2](e[1]), f[2](e[2]), f[2](e[3])
        (0, 1, 0)
        sage: f[3](e[1]), f[3](e[2]), f[3](e[3])
        (0, 0, 1)

    """
    def __init__(self, basis, symbol, latex_symbol=None):
        r"""
        TEST::
        
            sage: from sage.tensor.modules.free_module_basis import FreeModuleCoBasis
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: FreeModuleCoBasis(e, 'f')
            dual basis (f^0,f^1,f^2) on the rank-3 free module M over the Integer Ring
        
        """
        self._basis = basis
        self._fmodule = basis._fmodule
        self._name = "(" + \
          ",".join([symbol + "^" + str(i) for i in self._fmodule.irange()]) +")"
        if latex_symbol is None:
            latex_symbol = symbol
        self._latex_name = r"\left(" + \
          ",".join([latex_symbol + "^" + str(i) 
                    for i in self._fmodule.irange()]) + r"\right)"
        # The individual linear forms:
        vl = list()
        for i in self._fmodule.irange():
            v_name = symbol + "^" + str(i)
            v_symb = latex_symbol + "^" + str(i)
            v = self._fmodule.linear_form(name=v_name, latex_name=v_symb)
            for j in self._fmodule.irange():
                v.set_comp(basis)[j] = 0
            v.set_comp(basis)[i] = 1
            vl.append(v)
        self._form = tuple(vl)
        
    def _repr_(self):
        r"""
        String representation of the object.
        
        EXAMPLE::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = e.dual_basis()
            sage: f._repr_()
            'dual basis (e^0,e^1,e^2) on the rank-3 free module M over the Integer Ring'

        """
        return "dual basis " + self._name + " on the " + str(self._fmodule)

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        
        EXAMPLE::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = e.dual_basis()
            sage: f._latex_()
            '\\left(e^0,e^1,e^2\\right)'
            
        """
        return self._latex_name

    def __getitem__(self, index):
        r"""
        Returns the basis linear form corresponding to a given index.
        
        INPUT:
        
        - ``index`` -- the index of the linear form 
        
        EXAMPLE::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: f = e.dual_basis()
            sage: f.__getitem__(0)
            linear form e^0 on the rank-3 free module M over the Integer Ring
            sage: f.__getitem__(1)
            linear form e^1 on the rank-3 free module M over the Integer Ring
            sage: f.__getitem__(2)
            linear form e^2 on the rank-3 free module M over the Integer Ring
            sage: f[1] is f.__getitem__(1)
            True
            sage: M1 = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: f1 = M1.basis('e').dual_basis()
            sage: f1.__getitem__(1)
            linear form e^1 on the rank-3 free module M over the Integer Ring
            sage: f1.__getitem__(2)
            linear form e^2 on the rank-3 free module M over the Integer Ring
            sage: f1.__getitem__(3)
            linear form e^3 on the rank-3 free module M over the Integer Ring

        """
        n = self._fmodule._rank
        si = self._fmodule._sindex
        i = index - si
        if i < 0 or i > n-1:
            raise ValueError("Index out of range: " +
                              str(i+si) + " not in [" + str(si) + "," +
                              str(n-1+si) + "]")
        return self._form[i]

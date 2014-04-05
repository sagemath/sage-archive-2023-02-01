"""
Specific rank-2 tensors on free modules

Four derived classes of 
:class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor` are devoted 
to rank-2 tensors:

* :class:`FreeModuleEndomorphism` for endomorphisms (type-(1,1) tensors)

  * :class:`FreeModuleAutomorphism` for invertible endomorphisms

    * :class:`FreeModuleIdentityMap` for the identity map on a free module

* :class:`FreeModuleSymBilinForm` for symmetric bilinear forms (symmetric 
  type-(0,2) tensors)
  
Antisymmetric bilinear forms are dealt with by the class
:class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

EXAMPLES:

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

from free_module_tensor import FreeModuleTensor

class FreeModuleEndomorphism(FreeModuleTensor):
    r"""
    Endomorphism (considered as a type-(1,1) tensor) on a free module.

    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the endomorphism
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      endomorphism; if none is provided, the LaTeX symbol is set to ``name``
    
    EXAMPLES:
    
    Endomorphism on a rank-3 module::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: t = M.endomorphism('T') ; t
        endomorphism T on the rank-3 free module M over the Integer Ring

    An endomorphism is type-(1,1) tensor::

        sage: t.parent()
        free module of type-(1,1) tensors on the rank-3 free module M over the Integer Ring
        sage: t.tensor_type
        (1, 1)
        sage: t.tensor_rank
        2

    Consequently, an endomorphism can also be created by the module method 
    :meth:`~sage.tensor.modules.finite_free_module.FiniteFreeModule.tensor`::
    
        sage: t = M.tensor((1,1), name='T') ; t
        endomorphism T on the rank-3 free module M over the Integer Ring
    
    Components of the endomorphism with respect to a given basis::
    
        sage: e = M.basis('e') ; e
        basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
        sage: t[:] = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: t[:]
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: t.view()
        T = e_0*e^0 + 2 e_0*e^1 + 3 e_0*e^2 + 4 e_1*e^0 + 5 e_1*e^1 + 6 e_1*e^2 + 7 e_2*e^0 + 8 e_2*e^1 + 9 e_2*e^2

    The endomorphism acting on a module element::
    
        sage: v = M([1,2,3], basis=e, name='v') ; v
        element v of the rank-3 free module M over the Integer Ring
        sage: w = t(v) ; w
        element T(v) of the rank-3 free module M over the Integer Ring
        sage: w[:]
        [14, 32, 50]
        sage: # check:
        sage: for i in M.irange():
        ....:     print sum( t[i,j]*v[j] for j in M.irange() ),
        ....:     
        14 32 50
        
    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (1,1), name=name, 
                                  latex_name=latex_name)
                                  
    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "endomorphism "
        if self.name is not None:
            description += self.name + " " 
        description += "on the " + str(self.fmodule)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`FreeModuleEndomorphism` instance. 
        
        """
        return FreeModuleEndomorphism(self.fmodule)

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleTensor.__call__` to allow for a single 
        argument (module element). 
        """
        from free_module_tensor import FiniteFreeModuleElement
        if len(arg) > 1:
            # the endomorphism acting as a type (1,1) tensor on a pair 
            # (linear form, module element), returning a scalar:
            return FreeModuleTensor.__call__(self, *arg) 
        # the endomorphism acting as such, on a module element, returning a
        # module element:
        vector = arg[0]
        if not isinstance(vector, FiniteFreeModuleElement):
            raise TypeError("The argument must be an element of a free module.")
        basis = self.common_basis(vector)
        t = self.components[basis]
        v = vector.components[basis]
        fmodule = self.fmodule
        result = vector._new_instance()
        for i in fmodule.irange():
            res = 0
            for j in fmodule.irange():
                res += t[[i,j]]*v[[j]]
            result.set_comp(basis)[i] = res
        # Name of the output:
        result.name = None
        if self.name is not None and vector.name is not None:
            result.name = self.name + "(" + vector.name + ")"
        # LaTeX symbol for the output:
        result.latex_name = None
        if self.latex_name is not None and vector.latex_name is not None:
            result.latex_name = self.latex_name + r"\left(" + \
                              vector.latex_name + r"\right)"
        return result

#******************************************************************************

class FreeModuleAutomorphism(FreeModuleEndomorphism):
    r"""
    Automorphism (considered as a type-(1,1) tensor) on a free module.

    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the automorphism
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      automorphism; if none is provided, the LaTeX symbol is set to ``name``
    
    EXAMPLES:
    
    Automorphism on a rank-2 free module (vector space) on `\QQ`::
    
        sage: M = FiniteFreeModule(QQ, 2, name='M')
        sage: a = M.automorphism('A') ; a
        automorphism A on the rank-2 free module M over the Rational Field

    Automorphisms are tensors of type (1,1)::
    
        sage: a.parent()
        free module of type-(1,1) tensors on the rank-2 free module M over the Rational Field
        sage: a.tensor_type
        (1, 1)
        sage: a.tensor_rank
        2

    Setting the components in a basis::
    
        sage: e = M.basis('e') ; e
        basis (e_0,e_1) on the rank-2 free module M over the Rational Field
        sage: a[:] = [[1, 2], [-1, 3]]
        sage: a[:]
        [ 1  2]
        [-1  3]
        sage: a.view(basis=e)
        A = e_0*e^0 + 2 e_0*e^1 - e_1*e^0 + 3 e_1*e^1

    The inverse automorphism is obtained via the method :meth:`inverse`::
    
        sage: b = a.inverse() ; b
        automorphism A^(-1) on the rank-2 free module M over the Rational Field
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
        FreeModuleEndomorphism.__init__(self, fmodule, name=name, 
                                        latex_name=latex_name)
        self._inverse = None    # inverse automorphism not set yet

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "automorphism "
        if self.name is not None:
            description += self.name + " " 
        description += "on the " + str(self.fmodule)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`FreeModuleAutomorphism` instance. 
        """
        return FreeModuleAutomorphism(self.fmodule)
        
    def _del_derived(self):
        r"""
        Delete the derived quantities

        """
        # First delete the derived quantities pertaining to the mother class:
        FreeModuleEndomorphism._del_derived(self)
        # Then deletes the inverse automorphism:
        self._inverse = None        

    def inverse(self):
        r"""
        Return the inverse automorphism.
        
        OUTPUT:
        
        - instance of :class:`FreeModuleAutomorphism` representing the 
          automorphism that is the inverse of ``self``.
          
        EXAMPLES:
        
        Inverse of an automorphism on a rank-3 free module::
        
            sage: M = FiniteFreeModule(QQ, 3, name='M')
            sage: a = M.automorphism('A')
            sage: e = M.basis('e')
            sage: a[:] = [[1,0,-1], [0,3,0], [0,0,2]]
            sage: b = a.inverse() ; b
            automorphism A^(-1) on the rank-3 free module M over the Rational Field
            sage: b[:]
            [  1   0 1/2]
            [  0 1/3   0]
            [  0   0 1/2]
        
        We may check that b is the inverse of a by performing the matrix 
        product of the components in the basis e::
        
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
            if self.name is None:
                inv_name = None
            else:
                inv_name = self.name  + '^(-1)'
            if self.latex_name is None:
                inv_latex_name = None
            else:
                inv_latex_name = self.latex_name + r'^{-1}'
            fmodule = self.fmodule
            si = fmodule.sindex ; nsi = fmodule._rank + si
            self._inverse = FreeModuleAutomorphism(fmodule, inv_name, 
                                                   inv_latex_name)
            for basis in self.components:
                try:    
                    mat_self = matrix(
                              [[self.comp(basis)[[i, j]]
                              for j in range(si, nsi)] for i in range(si, nsi)])
                except (KeyError, ValueError):
                    continue
                mat_inv = mat_self.inverse()
                cinv = Components(fmodule.ring, basis, 2, start_index=si,
                                  output_formatter=fmodule.output_formatter)
                for i in range(si, nsi):
                    for j in range(si, nsi):   
                        cinv[i, j] = mat_inv[i-si,j-si]
                self._inverse.components[basis] = cinv
        return self._inverse


#******************************************************************************

class FreeModuleIdentityMap(FreeModuleAutomorphism):
    r"""
    Identity map (considered as a type-(1,1) tensor) on a free module.

    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: 'Id') name given to the identity map. 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the identity 
      map; if none is provided, the LaTeX symbol is set to ``name``
    
    EXAMPLES:
    
    Identity map on a rank-3 free module::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: a = M.identity_map() ; a
        identity map on the rank-3 free module M over the Integer Ring

    The LaTeX symbol is set by default to Id, but can be changed::
    
        sage: latex(a)
        \mathrm{Id}
        sage: a = M.identity_map(latex_name=r'\mathrm{1}')
        sage: latex(a)
        \mathrm{1}

    The identity map is a tensor of type (1,1) on the free module::
    
        sage: a.parent()
        free module of type-(1,1) tensors on the rank-3 free module M over the Integer Ring
        sage: a.tensor_type
        (1, 1)
        sage: a.tensor_rank
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
        NotImplementedError: The components of the identity map cannot be changed.

    The identity map acting on a module element::
    
        sage: v = M([2,-3,1], basis=e, name='v')
        sage: v.view()
        v = 2 e_0 - 3 e_1 + e_2
        sage: u = a(v) ; u
        element v of the rank-3 free module M over the Integer Ring
        sage: u is v
        True

    The identity map acting as a type-(1,1) tensor on a pair (linear form, 
    module element)::
    
        sage: w = M.tensor((0,1), name='w') ; w
        linear form w on the rank-3 free module M over the Integer Ring
        sage: w[:] = [0, 3, 2]
        sage: s = a(w,v) ; s
        -7
        sage: s == w(v)
        True
        
    The identity map is its own inverse::
    
        sage: a.inverse() == a
        True
        sage: a.inverse() is a
        True
    
    """
    def __init__(self, fmodule, name='Id', latex_name=None):
        if latex_name is None and name == 'Id':
            latex_name = r'\mathrm{Id}'
        FreeModuleAutomorphism.__init__(self, fmodule, name=name, 
                                        latex_name=latex_name)
        self._inverse = self    # the identity is its own inverse

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "identity map "
        if self.name != 'Id':
            description += self.name + " " 
        description += "on the " + str(self.fmodule)
        return description

    def _del_derived(self):
        r"""
        Delete the derived quantities.
        
        """
        # FreeModuleAutomorphism._del_derived is bypassed:
        FreeModuleEndomorphism._del_derived(self)

    def _new_instance(self):
        r"""
        Create a :class:`FreeModuleIdentityMap` instance. 
        
        """
        return FreeModuleIdentityMap(self.fmodule)
        
    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
                
        """
        from comp import KroneckerDelta
        fmodule = self.fmodule  # the base free module
        return KroneckerDelta(fmodule.ring, basis, start_index=fmodule.sindex,
                              output_formatter=fmodule.output_formatter)

    def comp(self, basis=None):
        r"""
        Return the components in a given basis, as a Kronecker delta.
                
        INPUT:
        
        - ``basis`` -- (default: None) module basis in which the components 
          are required; if none is provided, the components are assumed to 
          refer to the domain's default basis
 
        OUTPUT: 
        
        - components in the basis ``basis``, as an instance of the 
          class :class:`~sage.tensor.modules.comp.KroneckerDelta` 
        
        EXAMPLES:

        Components of the identity map on a rank-3 free module::
    
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_map()
            sage: a.comp(basis=e) 
            Kronecker delta of size 3x3
            
        For the module's default basis, the argument ``basis`` can be omitted::
        
            sage: a.comp() is a.comp(basis=e)
            True
            sage: a.comp()[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]
         
        """
        if basis is None: 
            basis = self.fmodule.def_basis
        if basis not in self.components:
            self.components[basis] = self._new_comp(basis)
        return self.components[basis]

    def set_comp(self, basis=None):
        r"""
        Redefinition of the generic tensor method 
        :meth:`FreeModuleTensor.set_comp`: should not be called.        

        """
        raise NotImplementedError("The components of the identity map " + 
                                  "cannot be changed.")

    def add_comp(self, basis=None):
        r"""
        Redefinition of the generic tensor method 
        :meth:`FreeModuleTensor.add_comp`: should not be called.        

        """
        raise NotImplementedError("The components of the identity map " + 
                                  "cannot be changed.")

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleEndomorphism.__call__`.
        """
        from free_module_tensor import FiniteFreeModuleElement
        from free_module_alt_form import FreeModuleLinForm
        if len(arg) == 1:
            # the identity map acting as such, on a module element:
            vector = arg[0]
            if not isinstance(vector, FiniteFreeModuleElement):
                raise TypeError("The argument must be a module element.")
            return vector
            #!# should it be return vector.copy() instead ? 
        elif len(arg) == 2:
            # the identity map acting as a type (1,1) tensor on a pair 
            # (1-form, vector), returning a scalar:
            linform = arg[0]
            if not isinstance(linform, FreeModuleLinForm):
                raise TypeError("The first argument must be a linear form.")
            vector = arg[1]
            if not isinstance(vector, FiniteFreeModuleElement):
                raise TypeError("The second argument must be a module element.")
            return linform(vector)
        else:
            raise TypeError("Bad number of arguments.")


#******************************************************************************

class FreeModuleSymBilinForm(FreeModuleTensor):
    r"""
    Symmetric bilinear form (considered as a type-(0,2) tensor) on a free 
    module.

    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the symmetric bilinear form
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the symmetric 
      bilinear form; if none is provided, the LaTeX symbol is set to ``name``
    
    EXAMPLES:
    
    Symmetric bilinear form on a rank-3 free module::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: a = M.sym_bilinear_form('A') ; a
        symmetric bilinear form A on the rank-3 free module M over the Integer Ring
        
    A symmetric bilinear form is a type-(0,2) tensor that is symmetric::
    
        sage: a.parent()
        free module of type-(0,2) tensors on the rank-3 free module M over the Integer Ring
        sage: a.tensor_type
        (0, 2)
        sage: a.tensor_rank
        2
        sage: a.symmetries()
        symmetry: (0, 1);  no antisymmetry
    
    Components with respect to a given basis::
    
        sage: e = M.basis('e')
        sage: a[0,0], a[0,1], a[0,2] = 1, 2, 3
        sage: a[1,1], a[1,2] = 4, 5
        sage: a[2,2] = 6
            
    Only independent components have been set; the other ones are deduced by 
    symmetry::
        
        sage: a[1,0], a[2,0], a[2,1]
        (2, 3, 5)
        sage: a[:]
        [1 2 3]
        [2 4 5]
        [3 5 6]
       
    A symmetric bilinear form acts on pairs of module elements::
    
        sage: u = M([2,-1,3]) ; v = M([-2,4,1])
        sage: a(u,v)
        61
        sage: a(v,u) == a(u,v)
        True
    
    The sum of two symmetric bilinear forms is another symmetric bilinear 
    form::

        sage: b = M.sym_bilinear_form('B')
        sage: b[0,0], b[0,1], b[1,2] = -2, 1, -3
        sage: s = a + b ; s
        symmetric bilinear form A+B on the rank-3 free module M over the Integer Ring
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [-2  1  0]  [-1  3  3]
        [2 4 5]  [ 1  0 -3]  [ 3  4  2]
        [3 5 6], [ 0 -3  0], [ 3  2  6]
        )
        
    Adding a symmetric bilinear from with a non-symmetric one results in a 
    generic type-(0,2) tensor::
    
        sage: c = M.tensor((0,2), name='C')
        sage: c[0,1] = 4
        sage: s = a + c ; s
        type-(0,2) tensor A+C on the rank-3 free module M over the Integer Ring
        sage: s.symmetries()
        no symmetry;  no antisymmetry
        sage: s[:]
        [1 6 3]
        [2 4 5]
        [3 5 6]
        


    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (0,2), name=name, 
                                  latex_name=latex_name, sym=(0,1))
                                  
    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "symmetric bilinear form "
        if self.name is not None:
            description += self.name + " " 
        description += "on the " + str(self.fmodule)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`FreeModuleSymBilinForm` instance. 
        
        """
        return FreeModuleSymBilinForm(self.fmodule)

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
                
        """
        from comp import CompFullySym
        fmodule = self.fmodule  # the base free module
        return CompFullySym(fmodule.ring, basis, 2, start_index=fmodule.sindex,
                            output_formatter=fmodule.output_formatter)



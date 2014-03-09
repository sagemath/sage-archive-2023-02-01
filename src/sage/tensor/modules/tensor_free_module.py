r"""
Tensor free modules


The class :class:`TensorFreeModule` implements the tensor products of the type  

.. MATH::

    T^{(k,l)}(M) = \underbrace{M\otimes\cdots\otimes M}_{k\ \; \mbox{times}}
    \otimes \underbrace{M^*\otimes\cdots\otimes M^*}_{l\ \; \mbox{times}}
    
where `M` is a free module of finite rank over a commutative ring `R` and
`M^*=\mathrm{Hom}_R(M,R)` is the dual of `M`. 
`T^{(k,l)}(M)` can be canonically identified with the set of tensors of 
type `(k,l)` acting as multilinear forms on `M`. 
Note that `T^{(1,0)}(M) = M`. 

`T^{(k,l)}(M)` is itself a free module over `R`, of rank `n^{k+l}`, `n`
being the rank of `M`. Accordingly the class :class:`TensorFreeModule` 
inherits from the class :class:`FiniteFreeModule`

Thanks to the canonical isomorphism `M^{**}\simeq M` (which holds because `M` 
is a free module of finite rank), `T^{(k,l)}(M)` can be identified with the 
set of tensors of type `(k,l)` defined as multilinear maps 

.. MATH::

    \underbrace{M^*\times\cdots\times M^*}_{k\ \; \mbox{times}}
    \times \underbrace{M\times\cdots\times M}_{l\ \; \mbox{times}}
    \longrightarrow R
    
Accordingly, :class:`TensorFreeModule` is a Sage *parent* class, whose
*element* class is
:class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`. 

TODO:

- implement more general tensor products, i.e. tensor product of the type 
  `M_1\otimes\cdots\otimes M_n`, where the `M_i`'s are `n` free modules of 
  finite rank over the same ring `R`. 

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

from sage.modules.module import Module
from finite_free_module import FiniteFreeModule
from free_module_tensor import FreeModuleTensor, FiniteFreeModuleElement


# From sage/modules/module.pyx:
# ----------------------------
### The new Module class that should be the base of all Modules
### The derived Module class must implement the element
### constructor:
#
# class MyModule(sage.modules.module.Module):
#     Element = MyElement
#     def _element_constructor_(self, x):
#         return self.element_class(x)
#


class TensorFreeModule(FiniteFreeModule):
    r"""
    Class for the free modules over a commutative ring `R` that are 
    tensor products of a given free module `M` over `R` with itself and its 
    dual `M^*`:
    
    .. MATH::
    
        T^{(k,l)}(M) = \underbrace{M\otimes\cdots\otimes M}_{k\ \; \mbox{times}}
        \otimes \underbrace{M^*\otimes\cdots\otimes M^*}_{l\ \; \mbox{times}}

    As recalled above, `T^{(k,l)}(M)` can be canonically identified with the 
    set of tensors of type `(k,l)` acting as multilinear forms on `M`. 

    INPUT:
    
    - ``fmodule`` -- free module `M` of finite rank (must be an instance of 
      :class:`FiniteFreeModule`)
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank and 
      `l` the covariant rank
    - ``name`` -- (string; default: None) name given to the tensor module
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
      tensor module; if none is provided, it is set to ``name``
    
    EXAMPLES:
    
    Set of tensors of type (1,2) on a free module of rank 3 over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: from sage.tensor.modules.tensor_free_module import TensorFreeModule
        sage: T = TensorFreeModule(M, (1,2)) ; T
        free module of type-(1,2) tensors on the rank-3 free module M over the Integer Ring
    
    Instead of importing TensorFreeModule in the global name space, it is 
    recommended to use the module's method 
    :meth:`~sage.tensor.modules.finite_free_module.FiniteFreeModule.tensor_module`::
    
        sage: T = M.tensor_module(1,2) ; T
        free module of type-(1,2) tensors on the rank-3 free module M over the Integer Ring

    The module `M` itself is considered as the set of tensors of type (1,0)::
    
        sage: M is M.tensor_module(1,0)
        True

    T is a module (actually a free module) over `\ZZ`::
    
        sage: T.category()
        Category of modules over Integer Ring
        sage: T in Modules(ZZ)
        True
        sage: T.rank()
        27
        sage: T.base_ring()
        Integer Ring
        sage: T.base_module()
        rank-3 free module M over the Integer Ring

    T is a *parent* object, whose elements are instances of 
    :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`::

        sage: t = T.an_element() ; t
        type-(1,2) tensor on the rank-3 free module M over the Integer Ring
        sage: from sage.tensor.modules.free_module_tensor import FreeModuleTensor
        sage: isinstance(t, FreeModuleTensor)
        True
        sage: t in T
        True
        sage: T.is_parent_of(t)
        True

    Elements can be constructed by means of the __call__ operator acting 
    on the parent; 0 yields the zero element::
        
        sage: T(0)
        type-(1,2) tensor zero on the rank-3 free module M over the Integer Ring
        sage: T(0) is T.zero()
        True
        
    while non-zero elements are constructed by providing their components in 
    a given basis::
    
        sage: e = M.basis('e') ; e
        basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
        sage: comp = [[[i-j+k for k in range(3)] for j in range(3)] for i in range(3)]
        sage: t = T(comp, basis=e, name='t') ; t
        type-(1,2) tensor t on the rank-3 free module M over the Integer Ring
        sage: t.comp(e)[:]
        [[[0, 1, 2], [-1, 0, 1], [-2, -1, 0]],
         [[1, 2, 3], [0, 1, 2], [-1, 0, 1]],
         [[2, 3, 4], [1, 2, 3], [0, 1, 2]]]
        sage: t.view(e)
        t = e_0*e^0*e^1 + 2 e_0*e^0*e^2 - e_0*e^1*e^0 + e_0*e^1*e^2 - 2 e_0*e^2*e^0 - e_0*e^2*e^1 + e_1*e^0*e^0 + 2 e_1*e^0*e^1 + 3 e_1*e^0*e^2 + e_1*e^1*e^1 + 2 e_1*e^1*e^2 - e_1*e^2*e^0 + e_1*e^2*e^2 + 2 e_2*e^0*e^0 + 3 e_2*e^0*e^1 + 4 e_2*e^0*e^2 + e_2*e^1*e^0 + 2 e_2*e^1*e^1 + 3 e_2*e^1*e^2 + e_2*e^2*e^1 + 2 e_2*e^2*e^2

    An alternative is to construct the tensor from an empty list of components
    and to set the nonzero components afterwards::
    
        sage: t = T([], name='t')
        sage: t.set_comp(e)[0,1,1] = -3
        sage: t.set_comp(e)[2,0,1] = 4
        sage: t.view(e)
        t = -3 e_0*e^1*e^1 + 4 e_2*e^0*e^1
    
    See the documentation of 
    :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor` 
    for the full list of arguments that can be provided to the __call__ 
    operator. For instance, to contruct a tensor symmetric with respect to the 
    last two indices::
    
        sage: t = T([], name='t', sym=(1,2))
        sage: t.set_comp(e)[0,1,1] = -3
        sage: t.set_comp(e)[2,0,1] = 4
        sage: t.view(e)  # notice that t^2_{10} has be set equal to t^2_{01} by symmetry
        t = -3 e_0*e^1*e^1 + 4 e_2*e^0*e^1 + 4 e_2*e^1*e^0
    
    The tensor modules over a given module `M` are unique::
    
        sage: T is M.tensor_module(1,2)
        True
            
    """
    
    Element = FreeModuleTensor
    
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None):
        self.fmodule = fmodule
        self.tensor_type = tuple(tensor_type)
        rank = pow(fmodule._rank, tensor_type[0] + tensor_type[1])
        self._zero_element = 0 # provisory (to avoid infinite recursion in what
                               # follows)
        if tensor_type == (0,1):  # case of the dual
            if name is None and fmodule.name is not None:
                name = fmodule.name + '*'
            if latex_name is None and fmodule.latex_name is not None:
                latex_name = fmodule.latex_name + r'^*'
        FiniteFreeModule.__init__(self, fmodule.ring, rank, name=name, 
                                  latex_name=latex_name, 
                                  start_index=fmodule.sindex,
                                  output_formatter=fmodule.output_formatter)
        # Unique representation:
        if self.tensor_type in self.fmodule._tensor_modules:
            raise TypeError("The module of tensors of type" + 
                            str(self.tensor_type) + 
                            " has already been created.")
        else:
            self.fmodule._tensor_modules[self.tensor_type] = self
        # Zero element 
        self._zero_element = self._element_constructor_(name='zero', 
                                                        latex_name='0')
        def_basis = self.fmodule.def_basis
        if def_basis is not None:
            self._zero_element.components[def_basis] = \
                                        self._zero_element._new_comp(def_basis)
            # (since new components are initialized to zero)
    
    #### Methods required for any Parent 
    def _element_constructor_(self, comp=[], basis=None, name=None, 
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor
        """
        if comp == 0:
            return self._zero_element
        resu = self.element_class(self.fmodule, self.tensor_type, name=name, 
                                  latex_name=latex_name, sym=sym, 
                                  antisym=antisym)
        if comp != []:
            resu.set_comp(basis)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) tensor
        """
        resu = self.element_class(self.fmodule, self.tensor_type)
        if self.fmodule.def_basis is not None:
            sindex = self.fmodule.sindex
            ind = [sindex for i in range(resu.tensor_rank)]
            resu.set_comp()[ind] = self.fmodule.ring.an_element()
        return resu
            
    #### End of methods required for any Parent 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        if self.tensor_type == (0,1):
            return "dual of the " + str(self.fmodule)
        else:
            description = "free module "
            if self.name is not None:
                description += self.name + " "
            description += "of type-(%s,%s)" % \
                           (str(self.tensor_type[0]), str(self.tensor_type[1]))
            description += " tensors on the " + str(self.fmodule)
            return description

    def base_module(self):
        r""" 
        Return the free module on which ``self`` is constructed.
        
        OUTPUT:
        
        - instance of :class:`FiniteFreeModule` representing the free module 
          on which the tensor module is defined. 
        
        EXAMPLE:

        Base module of a type-(1,2) tensor module::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: T = M.tensor_module(1,2)
            sage: T.base_module()
            rank-3 free module M over the Integer Ring

        """
        return self.fmodule


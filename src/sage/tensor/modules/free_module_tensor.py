r"""
Tensors on free modules

The class :class:`FreeModuleTensor` implements tensors over a free module `M`,
i.e. elements of the free module `T^{(k,l)}(M)` of tensors of type `(k,l)`
acting as multilinear forms on `M`. 

A *tensor of type* `(k,l)` is a multilinear map:

.. MATH::

    \underbrace{M^*\times\cdots\times M^*}_{k\ \; \mbox{times}}
    \times \underbrace{M\times\cdots\times M}_{l\ \; \mbox{times}}
    \longrightarrow R
    
where `R` is the commutative ring over which the free module `M` is defined and
`M^*=\mathrm{Hom}_R(M,R)` is the dual of `M`. The integer `k+l` is called the 
*tensor rank*. 

Various derived classes of :class:`FreeModuleTensor` are devoted to specific 
tensors:

* :class:`FiniteFreeModuleElement` for elements of `M`, considered as 
  type-(1,0) tensors thanks to the canonical identification `M^{**}=M`, which
  holds since `M` is a free module of finite rank

* :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm` for 
  fully antisymmetric type-`(0,l)` tensors (alternating forms)

  * :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleLinForm` for 
    type-(0,1) tensors (linear forms)

* :class:`~sage.tensor.modules.free_module_tensor_spec.FreeModuleEndomorphism` 
  for type-(1,1) tensors (endomorphisms)

  * :class:`~sage.tensor.modules.free_module_tensor_spec.FreeModuleAutomorphism` 
    for invertible endomorphisms

    * :class:`~sage.tensor.modules.free_module_tensor_spec.FreeModuleIdentityMap` 
      for the identity map on a free module

* :class:`~sage.tensor.modules.free_module_tensor_spec.FreeModuleSymBilinForm` 
  for symmetric type-(0,2) tensors (symmetric bilinear forms)

:class:`FreeModuleTensor` is a Sage *element* class, the corresponding *parent*
class being :class:`~sage.tensor.modules.tensor_free_module.TensorFreeModule`. 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

EXAMPLES:

    A tensor of type (1,1) on a rank-3 free module over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: t = M.tensor((1,1), name='t') ; t
        endomorphism t on the rank-3 free module M over the Integer Ring
        sage: t.parent()
        free module of type-(1,1) tensors on the rank-3 free module M over the Integer Ring
        sage: t.parent() is M.tensor_module(1,1)
        True
        sage: t in M.tensor_module(1,1)
        True
        
    Setting some component of the tensor in a given basis::
    
        sage: e = M.basis('e') ; e
        basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
        sage: t.set_comp(e)[0,0] = -3  # the component [0,0] w.r.t. basis e is set to -3
        
    The unset components are assumed to be zero::
    
        sage: t.comp(e)[:]  # list of all components w.r.t. basis e
        [-3  0  0]
        [ 0  0  0]
        [ 0  0  0]
        sage: t.view(e)  # expansion of t on the basis e_i*e^j of T^(1,1)(M) 
        t = -3 e_0*e^0

    Since e is M's default basis, shorcuts for the above writings are::
    
        sage: t[0,0] = -3
        sage: t[:]
        [-3  0  0]
        [ 0  0  0]
        [ 0  0  0]
    
    Tensor components can be modified (reset) at any time::
    
        sage: t[0,0] = 0
        sage: t[:]
        [0 0 0]
        [0 0 0]
        [0 0 0]

    Checking that t is zero::
    
        sage: t.is_zero()
        True
        sage: t == 0
        True
        sage: t == M.tensor_module(1,1).zero()  # the zero element of the module of all type-(1,1) tensors on M
        True

    The components are managed by the class :class:`~sage.tensor.modules.comp.Components`::
    
        sage: type(t.comp(e))
        <class 'sage.tensor.modules.comp.Components'>

    Only non-zero components are actually stored, in the dictionary :attr:`_comp`
    of class :class:`~sage.tensor.modules.comp.Components`, whose keys are the indices::
    
        sage: t.comp(e)._comp
        {}
        sage: t.set_comp(e)[0,0] = -3 ; t.set_comp(e)[1,2] = 2
        sage: t.comp(e)._comp  # random output order (dictionary)
        {(0, 0): -3, (1, 2): 2}
        sage: t.view(e)
        t = -3 e_0*e^0 + 2 e_1*e^2

    Further tests of the comparison operator::
    
        sage: t.is_zero()
        False
        sage: t == 0
        False
        sage: t == M.tensor_module(1,1).zero()
        False
        sage: t1 = t.copy()
        sage: t1 == t
        True
        sage: t1[2,0] = 4
        sage: t1 == t
        False

    As a multilinear map `M^*\times M \rightarrow \ZZ`, the type-(1,1) tensor t 
    acts on pairs formed by a linear form and a module element::
    
        sage: a = M.linear_form(name='a') ; a[:] = (2, 1, -3) ; a
        linear form a on the rank-3 free module M over the Integer Ring
        sage: b = M([1,-6,2], name='b') ; b
        element b of the rank-3 free module M over the Integer Ring
        sage: t(a,b)
        -2

    
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

from sage.rings.integer import Integer
from sage.structure.element import ModuleElement  
#!# or from sage.structure.element import Element
# to avoid arithmetics defined in ModuleElement ??

from comp import Components, CompWithSym, CompFullySym, CompFullyAntiSym

class FreeModuleTensor(ModuleElement):
    r"""
    Tensor over a free module of finite rank over a commutative ring.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` (must be an 
      instance of :class:`FiniteFreeModule`)
    - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
      the covariant rank
    - ``name`` -- (default: None) name given to the tensor
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor; 
      if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: None) a symmetry or a list of symmetries among the 
      tensor arguments: each symmetry is described by a tuple containing 
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:

      * sym=(0,1) for a symmetry between the 1st and 2nd arguments 
      * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
         arguments and a symmetry between the 2nd, 4th and 5th arguments.

    - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
      among the arguments, with the same convention as for ``sym``. 

    EXAMPLES:

    A tensor of type (1,1) on a rank-3 free module over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: t = M.tensor((1,1), name='t') ; t
        endomorphism t on the rank-3 free module M over the Integer Ring

    Tensors are *Element* objects whose parents are tensor free modules::
    
        sage: t.parent()
        free module of type-(1,1) tensors on the rank-3 free module M over the Integer Ring
        sage: t.parent() is M.tensor_module(1,1)
        True
        
    """
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None,
                 sym=None, antisym=None):
        ModuleElement.__init__(self, fmodule.tensor_module(*tensor_type))
        self.fmodule = fmodule
        self.tensor_type = tuple(tensor_type)
        self.tensor_rank = self.tensor_type[0] + self.tensor_type[1]
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.components = {}    # components on various bases (not set yet)
        # Treatment of symmetry declarations:
        self.sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):  
                # a single symmetry is provided as a tuple -> 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                 " not in [0," + str(self.tensor_rank-1) + "]")
                    self.sym.append(tuple(isym))       
        self.antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):  
                # a single antisymmetry is provided as a tuple -> 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                " not in [0," + str(self.tensor_rank-1) + "]")
                    self.antisym.append(tuple(isym))
        # Final consistency check:
        index_list = []
        for isym in self.sym:
            index_list += isym
        for isym in self.antisym:
            index_list += isym
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("Incompatible lists of symmetries: the same " + 
                             "position appears more than once.")
        # Initialization of derived quantities:
        FreeModuleTensor._init_derived(self) 

    ####### Required methods for ModuleElement (beside arithmetic) #######
    
    def __nonzero__(self):
        r"""
        Return True if ``self`` is nonzero and False otherwise. 
        
        This method is called by self.is_zero(). 
        """
        basis = self.pick_a_basis()
        return not self.components[basis].is_zero()
        
    ####### End of required methods for ModuleElement (beside arithmetic) #######
    
    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "type-(%s,%s) tensor" % \
                           (str(self.tensor_type[0]), str(self.tensor_type[1]))
        if self.name is not None:
            description += " " + self.name
        description += " on the " + str(self.fmodule)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self.latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        pass # no derived quantities

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        pass # no derived quantities

    def view(self, basis=None, format_spec=None):
        r"""
        Displays the tensor in terms of its expansion onto a given basis.
        
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        INPUT:
                
        - ``basis`` -- (default: None) basis of the free module with respect to 
          which the tensor is expanded; if none is provided, the module's 
          default basis is assumed
        - ``format_spec`` -- (default: None) format specification passed to 
          ``self.fmodule.output_formatter`` to format the output.

        EXAMPLES:
        
        Display of a module element (type-(1,0) tensor)::
        
            sage: M = FiniteFreeModule(QQ, 2, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: v = M([1/3,-2], name='v')
            sage: v.view()
            v = 1/3 e_1 - 2 e_2
            sage: latex(v.view())  # display in the notebook
            v = \frac{1}{3} e_1 -2 e_2

        Display of a linear form (type-(0,1) tensor)::
        
            sage: de = e.dual_basis()
            sage: w = - 3/4 * de[1] + de[2] ; w
            linear form on the rank-2 free module M over the Rational Field
            sage: w.set_name('w', latex_name='\omega')
            sage: w.view()
            w = -3/4 e^1 + e^2
            sage: latex(w.view())  # display in the notebook
            \omega = -\frac{3}{4} e^1 +e^2

        Display of a type-(1,1) tensor::
        
            sage: t = v*w ; t  # the type-(1,1) is formed as the tensor product of v by w
            endomorphism v*w on the rank-2 free module M over the Rational Field
            sage: t.view()
            v*w = -1/4 e_1*e^1 + 1/3 e_1*e^2 + 3/2 e_2*e^1 - 2 e_2*e^2
            sage: latex(t.view())  # display in the notebook
            v\otimes \omega = -\frac{1}{4} e_1\otimes e^1 + \frac{1}{3} e_1\otimes e^2 + \frac{3}{2} e_2\otimes e^1 -2 e_2\otimes e^2

        Display in a basis which is not the default one::
        
            sage: a = M.automorphism()
            sage: a[:] = [[1,2],[3,4]]
            sage: f = e.new_basis(a, 'f')
            sage: v.view(f) # the components w.r.t basis f are first computed via the change-of-basis formula defined by a
            v = -8/3 f_1 + 3/2 f_2
            sage: w.view(f)
            w = 9/4 f^1 + 5/2 f^2
            sage: t.view(f)
            v*w = -6 f_1*f^1 - 20/3 f_1*f^2 + 27/8 f_2*f^1 + 15/4 f_2*f^2

        The output format can be set via the argument ``output_formatter`` 
        passed at the module construction::

            sage: N = FiniteFreeModule(QQ, 2, name='N', start_index=1, output_formatter=Rational.numerical_approx)
            sage: e = N.basis('e')
            sage: v = N([1/3,-2], name='v')
            sage: v.view()  # default format (53 bits of precision)
            v = 0.333333333333333 e_1 - 2.00000000000000 e_2
            sage: latex(v.view()) 
            v = 0.333333333333333 e_1 -2.00000000000000 e_2
            
        The output format is then controled by the argument ``format_spec`` of
        the method :meth:`view`::
        
            sage: v.view(format_spec=10)  # 10 bits of precision
            v = 0.33 e_1 - 2.0 e_2
                    
        """
        from sage.misc.latex import latex
        from format_utilities import is_atomic, FormattedExpansion
        if basis is None:
            basis = self.fmodule.def_basis
        cobasis = basis.dual_basis()
        comp = self.comp(basis)
        terms_txt = []
        terms_latex = []
        n_con = self.tensor_type[0]
        for ind in comp.index_generator():
            ind_arg = ind + (format_spec,)
            coef = comp[ind_arg]
            if coef != 0:
                bases_txt = []
                bases_latex = []
                for k in range(n_con):
                    bases_txt.append(basis[ind[k]].name)
                    bases_latex.append(latex(basis[ind[k]]))
                for k in range(n_con, self.tensor_rank):
                    bases_txt.append(cobasis[ind[k]].name)
                    bases_latex.append(latex(cobasis[ind[k]]))
                basis_term_txt = "*".join(bases_txt)    
                basis_term_latex = r"\otimes ".join(bases_latex)    
                if coef == 1:
                    terms_txt.append(basis_term_txt)
                    terms_latex.append(basis_term_latex)
                elif coef == -1:
                    terms_txt.append("-" + basis_term_txt)
                    terms_latex.append("-" + basis_term_latex)
                else:
                    coef_txt = repr(coef)
                    coef_latex = latex(coef)
                    if is_atomic(coef_txt):
                        terms_txt.append(coef_txt + " " + basis_term_txt)
                    else:
                        terms_txt.append("(" + coef_txt + ") " + 
                                         basis_term_txt)
                    if is_atomic(coef_latex):
                        terms_latex.append(coef_latex + basis_term_latex)
                    else:
                        terms_latex.append(r"\left(" + coef_latex + r"\right)" + 
                                           basis_term_latex)
        if terms_txt == []:
            expansion_txt = "0"
        else:
            expansion_txt = terms_txt[0]
            for term in terms_txt[1:]:
                if term[0] == "-":
                    expansion_txt += " - " + term[1:]
                else:
                    expansion_txt += " + " + term
        if terms_latex == []:
            expansion_latex = "0"
        else:
            expansion_latex = terms_latex[0]
            for term in terms_latex[1:]:
                if term[0] == "-":
                    expansion_latex += term
                else:
                    expansion_latex += "+" + term
        result = FormattedExpansion(self)            
        if self.name is None:
            result.txt = expansion_txt
        else:
            result.txt = self.name + " = " + expansion_txt
        if self.latex_name is None:
            result.latex = expansion_latex
        else:
            result.latex = latex(self) + " = " + expansion_latex
        return result
    
    def symmetries(self):
        r"""
        Print the list of symmetries and antisymmetries.
        
        EXAMPLES:
        
        Various symmetries / antisymmetries for a rank-4 tensor::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: t = M.tensor((4,0), name='T') # no symmetry declared
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: t = M.tensor((4,0), name='T', sym=(0,1))
            sage: t.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: t = M.tensor((4,0), name='T', sym=[(0,1), (2,3)])
            sage: t.symmetries()
            symmetries: [(0, 1), (2, 3)];  no antisymmetry
            sage: t = M.tensor((4,0), name='T', sym=(0,1), antisym=(2,3))
            sage: t.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
            
        """
        if len(self.sym) == 0:
            s = "no symmetry; "
        elif len(self.sym) == 1:
            s = "symmetry: " + str(self.sym[0]) + "; "
        else:
            s = "symmetries: " + str(self.sym) + "; " 
        if len(self.antisym) == 0:
            a = "no antisymmetry"
        elif len(self.antisym) == 1:
            a = "antisymmetry: " + str(self.antisym[0])
        else:
            a = "antisymmetries: " + str(self.antisym)   
        print s, a
         
    def set_name(self, name, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of the tensor.

        INPUT:
        
        - ``name`` -- name given to the tensor
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor; 
          if none is provided, the LaTeX symbol is set to ``name``

        """
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
       
    def _new_instance(self):
        r"""
        Create a :class:`FreeModuleTensor` instance of the same tensor type and 
        with the same symmetries.

        This method must be redefined by derived classes of 
        :class:`FreeModuleTensor`.
        
        """
        return FreeModuleTensor(self.fmodule, self.tensor_type, sym=self.sym, 
                                antisym=self.antisym)

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
        
        This method, to be called by :meth:`comp`, must be redefined by derived 
        classes to adapt the output to the relevant subclass of 
        :class:`~sage.tensor.modules.comp.Components`.
        
        OUTPUT:
        
        - an instance of :class:`~sage.tensor.modules.comp.Components` (or of one of its subclass)
        
        """
        fmodule = self.fmodule  # the base free module
        if self.sym == [] and self.antisym == []:
            return Components(fmodule.ring, basis, self.tensor_rank,
                              start_index=fmodule.sindex,
                              output_formatter=fmodule.output_formatter)
        for isym in self.sym:
            if len(isym) == self.tensor_rank:
                return CompFullySym(fmodule.ring, basis, self.tensor_rank,
                                    start_index=fmodule.sindex,
                                    output_formatter=fmodule.output_formatter)
        for isym in self.antisym:
            if len(isym) == self.tensor_rank:
                return CompFullyAntiSym(fmodule.ring, basis, self.tensor_rank, 
                                        start_index=fmodule.sindex,
                                     output_formatter=fmodule.output_formatter)
        return CompWithSym(fmodule.ring, basis, self.tensor_rank, 
                           start_index=fmodule.sindex, 
                           output_formatter=fmodule.output_formatter,
                           sym=self.sym, antisym=self.antisym)        

    def comp(self, basis=None, from_basis=None):
        r"""
        Return the components in a given basis.
        
        If the components are not known already, they are computed by the tensor
        change-of-basis formula from components in another basis. 
        
        INPUT:
        
        - ``basis`` -- (default: None) basis in which the components are 
          required; if none is provided, the components are assumed to refer to
          the module's default basis
        - ``from_basis`` -- (default: None) basis from which the
          required components are computed, via the tensor change-of-basis 
          formula, if they are not known already in the basis ``basis``; 
          if none, a basis is picked in ``self.components``.
 
        OUTPUT: 
        
        - components in the basis ``basis``, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components` 
        
        EXAMPLES:
        
        Components of a tensor of type-(1,1)::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t[1,2] = -3 ; t[3,3] = 2
            sage: t.comp()
            2-indices components w.r.t. basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring
            sage: t.comp() is t.comp(e)  # since e is M's default basis
            True
            sage: t.comp()[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]

        A direct access to the components w.r.t. the module's default basis is
        provided by the square brackets applied to the tensor itself::

            sage: t[1,2] is t.comp(e)[1,2]
            True
            sage: t[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]

        Components computed via a change-of-basis formula::
        
            sage: a = M.automorphism()
            sage: a[:] = [[0,0,1], [1,0,0], [0,-1,0]]
            sage: f = e.new_basis(a, 'f')
            sage: t.comp(f)
            2-indices components w.r.t. basis (f_1,f_2,f_3) on the rank-3 free module M over the Integer Ring
            sage: t.comp(f)[:]
            [ 0  0  0]
            [ 0  2  0]
            [-3  0  0]
            
        """
        fmodule = self.fmodule
        if basis is None: 
            basis = fmodule.def_basis
        if basis not in self.components:
            # The components must be computed from 
            # those in the basis from_basis
            if from_basis is None: 
                for known_basis in self.components:
                    if (known_basis, basis) in self.fmodule.basis_changes \
                        and (basis, known_basis) in self.fmodule.basis_changes:
                        from_basis = known_basis
                        break
                if from_basis is None:
                    raise ValueError("No basis could be found for computing " + 
                                     "the components in the " + str(basis))
            elif from_basis not in self.components:
                raise ValueError("The tensor components are not known in the " +
                                 "basis "+ str(from_basis))
            (n_con, n_cov) = self.tensor_type
            if n_cov > 0:
                if (from_basis, basis) not in fmodule.basis_changes:
                    raise ValueError("The change-of-basis matrix from the " + 
                                     str(from_basis) + " to the " + str(basis) 
                                     + " has not been set.")
                pp = \
                  fmodule.basis_changes[(from_basis, basis)].comp(from_basis)
                # pp not used if n_cov = 0 (pure contravariant tensor)
            if n_con > 0:
                if (basis, from_basis) not in fmodule.basis_changes:
                    raise ValueError("The change-of-basis matrix from the " + 
                                     str(basis) + " to the " + str(from_basis) +
                                     " has not been set.")
                ppinv = \
                  fmodule.basis_changes[(basis, from_basis)].comp(from_basis)
                # ppinv not used if n_con = 0 (pure covariant tensor)
            old_comp = self.components[from_basis]
            new_comp = self._new_comp(basis)
            rank = self.tensor_rank
            # loop on the new components:
            for ind_new in new_comp.non_redundant_index_generator(): 
                # Summation on the old components multiplied by the proper 
                # change-of-basis matrix elements (tensor formula): 
                res = 0 
                for ind_old in old_comp.index_generator(): 
                    t = old_comp[[ind_old]]
                    for i in range(n_con): # loop on contravariant indices
                        t *= ppinv[[ind_new[i], ind_old[i]]]
                    for i in range(n_con,rank):  # loop on covariant indices
                        t *= pp[[ind_old[i], ind_new[i]]]
                    res += t
                new_comp[ind_new] = res
            self.components[basis] = new_comp
            # end of case where the computation was necessary
        return self.components[basis]

    def set_comp(self, basis=None):
        r"""
        Return the components in a given basis for assignment.
        
        The components with respect to other bases are deleted, in order to 
        avoid any inconsistency. To keep them, use the method :meth:`add_comp` 
        instead.
        
        INPUT:
        
        - ``basis`` -- (default: None) basis in which the components are
          defined; if none is provided, the components are assumed to refer to 
          the module's default basis.
         
        OUTPUT: 
        
        - components in the given basis, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components`; if such components did not exist
          previously, they are created.  
        
        EXAMPLES:
        
        Setting components of a type-(1,1) tensor::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t.set_comp()[0,1] = -3
            sage: t.view()
            t = -3 e_0*e^1
            sage: t.set_comp()[1,2] = 2
            sage: t.view()
            t = -3 e_0*e^1 + 2 e_1*e^2
            sage: t.set_comp(e)
            2-indices components w.r.t. basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
            
        Setting components in a new basis::
        
            sage: f =  M.basis('f')
            sage: t.set_comp(f)[0,1] = 4
            sage: t.components.keys() # the components w.r.t. basis e have been deleted
            [basis (f_0,f_1,f_2) on the rank-3 free module M over the Integer Ring]
            sage: t.view(f)
            t = 4 f_0*f^1
            
        The components w.r.t. basis e can be deduced from those w.r.t. basis f,
        once a relation between the two bases has been set::
        
            sage: a = M.automorphism()
            sage: a[:] = [[0,0,1], [1,0,0], [0,-1,0]]
            sage: M.basis_changes[(e,f)] = a
            sage: M.basis_changes[(f,e)] = a.inverse()
            sage: t.view(e)
            t = -4 e_1*e^2
            sage: t.components.keys()  # random output (dictionary keys)
            [basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring,
             basis (f_0,f_1,f_2) on the rank-3 free module M over the Integer Ring]

        """
        if self is self.parent()._zero_element: #!# this is maybe not very efficient
            raise ValueError("The zero element cannot be changed.")
        if basis is None: 
            basis = self.fmodule.def_basis
        if basis not in self.components:
            if basis not in self.fmodule.known_bases:
                raise ValueError("The " + str(basis) + " has not been " +
                                 "defined on the " + str(self.fmodule))
            self.components[basis] = self._new_comp(basis)
        self._del_derived() # deletes the derived quantities
        self.del_other_comp(basis)
        return self.components[basis]

    def add_comp(self, basis=None):
        r"""
        Return the components in a given basis for assignment, keeping the
        components in other bases. 
        
        To delete the components in other bases, use the method 
        :meth:`set_comp` instead. 
        
        INPUT:
        
        - ``basis`` -- (default: None) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis.
          
        .. WARNING::
        
            If the tensor has already components in other bases, it 
            is the user's responsability to make sure that the components
            to be added are consistent with them. 
         
        OUTPUT: 
        
        - components in the given basis, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components`; if such components did not exist
          previously, they are created.  
        
        EXAMPLES:
        
        Setting components of a type-(1,1) tensor::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t.add_comp()[0,1] = -3
            sage: t.view()
            t = -3 e_0*e^1
            sage: t.add_comp()[1,2] = 2
            sage: t.view()
            t = -3 e_0*e^1 + 2 e_1*e^2
            sage: t.add_comp(e)
            2-indices components w.r.t. basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
            
        Adding components in a new basis::
        
            sage: f =  M.basis('f')
            sage: t.add_comp(f)[0,1] = 4
            
        The components w.r.t. basis e have been kept::
        
            sage: t.components.keys() # # random output (dictionary keys) 
            [basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring,
             basis (f_0,f_1,f_2) on the rank-3 free module M over the Integer Ring]
            sage: t.view(f)
            t = 4 f_0*f^1
            sage: t.view(e)
            t = -3 e_0*e^1 + 2 e_1*e^2

        """
        if basis is None: basis = self.fmodule.def_basis
        if basis not in self.components:
            if basis not in self.fmodule.known_bases:
                raise ValueError("The " + str(basis) + " has not been " +
                                 "defined on the " + str(self.fmodule))
            self.components[basis] = self._new_comp(basis)
        self._del_derived() # deletes the derived quantities
        return self.components[basis]


    def del_other_comp(self, basis=None):
        r"""
        Delete all the components but those corresponding to ``basis``.
        
        INPUT:
        
        - ``basis`` -- (default: None) basis in which the components are
          kept; if none the module's default basis is assumed

        EXAMPLE:
        
        Deleting components of a module element::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: u = M([2,1,-5])
            sage: f = M.basis('f')
            sage: u.add_comp(f)[:] = [0,4,2]
            sage: u.components.keys() # random output (dictionary keys)
            [basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring,
             basis (f_1,f_2,f_3) on the rank-3 free module M over the Integer Ring]
            sage: u.del_other_comp(f)
            sage: u.components.keys()
            [basis (f_1,f_2,f_3) on the rank-3 free module M over the Integer Ring]

        Let us restore the components w.r.t. e and delete those w.r.t. f::
        
            sage: u.add_comp(e)[:] = [2,1,-5]
            sage: u.components.keys()  # random output (dictionary keys)
            [basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring,
             basis (f_1,f_2,f_3) on the rank-3 free module M over the Integer Ring]
            sage: u.del_other_comp()  # default argument: basis = e
            sage: u.components.keys()
            [basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring]

        """
        if basis is None: basis = self.fmodule.def_basis
        if basis not in self.components:
            raise ValueError("The components w.r.t. the " + 
                             str(basis) + " have not been defined.")
        to_be_deleted = []
        for other_basis in self.components:
            if other_basis != basis:
                to_be_deleted.append(other_basis)
        for other_basis in to_be_deleted:
            del self.components[other_basis]

    def __getitem__(self, indices):
        r"""
        Return a component w.r.t. the free module's default basis.

        INPUT:
        
        - ``indices`` -- list of indices defining the component
    
        """
        return self.comp()[indices]
        
    def __setitem__(self, indices, value):
        r"""
        Set a component w.r.t. the free module's default basis.

        INPUT:
        
        - ``indices`` -- list of indices defining the component
    
        """        
        self.set_comp()[indices] = value


    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        The name and the derived quantities are not copied. 
        
        EXAMPLES:
        
        Copy of a tensor of type (1,1)::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: t = M.tensor((1,1), name='t')
            sage: t[1,2] = -3 ; t[3,3] = 2
            sage: t1 = t.copy()
            sage: t1[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]
            sage: t1 == t
            True

        If the original tensor is modified, the copy is not::
        
            sage: t[2,2] = 4
            sage: t1[:]
            [ 0 -3  0]
            [ 0  0  0]
            [ 0  0  2]
            sage: t1 == t
            False

        """
        resu = self._new_instance()
        for basis, comp in self.components.items():
             resu.components[basis] = comp.copy()
        return resu

    def common_basis(self, other):
        r"""
        Find a common basis for the components of ``self`` and ``other``. 
        
        In case of multiple common bases, the free module's default basis is 
        privileged. 
        If the current components of ``self`` and ``other`` are all relative to
        different bases, a common basis is searched by performing a component
        transformation, via the transformations listed in 
        ``self.fmodule.basis_changes``, still privileging transformations to 
        the free module's default basis.
        
        INPUT:
        
        - ``other`` -- a tensor (instance of :class:`FreeModuleTensor`)
        
        OUPUT:
        
        - instance of 
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing the common basis; if no common basis is found, None is 
          returned. 
          
        EXAMPLES:
        
        Common basis for the components of two module elements::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: u = M([2,1,-5])
            sage: f = M.basis('f')
            sage: M.basis_changes.clear() # to ensure that bases e and f are unrelated at this stage
            sage: v = M([0,4,2], basis=f)
            sage: u.common_basis(v) 
            
        The above result is None since u and v have been defined on different
        bases and no connection between these bases have been set::
        
            sage: u.components.keys()
            [basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring]
            sage: v.components.keys()
            [basis (f_1,f_2,f_3) on the rank-3 free module M over the Integer Ring]

        Linking bases e and f changes the result::
        
            sage: a = M.automorphism()
            sage: a[:] = [[0,0,1], [1,0,0], [0,-1,0]]
            sage: M.basis_changes[(e,f)] = a
            sage: M.basis_changes[(f,e)] = a.inverse()
            sage: u.common_basis(v)
            basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring

        Indeed, v is now known in basis e::
        
            sage: v.components.keys() # random output (dictionary keys)
            [basis (f_1,f_2,f_3) on the rank-3 free module M over the Integer Ring,
             basis (e_1,e_2,e_3) on the rank-3 free module M over the Integer Ring]

        """
        # Compatibility checks:
        if not isinstance(other, FreeModuleTensor):
            raise TypeError("The argument must be a tensor on a free module.")
        fmodule = self.fmodule
        if other.fmodule != fmodule:
            raise TypeError("The two tensors are not defined on the same " +
                            "free module.")
        def_basis = fmodule.def_basis
        #
        # 1/ Search for a common basis among the existing components, i.e. 
        #    without performing any component transformation. 
        #    -------------------------------------------------------------
        if def_basis in self.components and def_basis in other.components:
            return def_basis # the module's default basis is privileged
        for basis1 in self.components:
            if basis1 in other.components:
                return basis1
        # 2/ Search for a common basis via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least 
        # one component transformation to get a common basis
        if def_basis in self.components:
            for obasis in other.components:
                if (obasis, def_basis) in fmodule.basis_changes:
                    other.comp(def_basis, from_basis=obasis)
                    return def_basis
        if def_basis in other.components:
            for sbasis in self.components:
                if (sbasis, def_basis) in fmodule.basis_changes:
                    self.comp(def_basis, from_basis=sbasis)
                    return def_basis
        # If this point is reached, then def_basis cannot be a common basis
        # via a single component transformation
        for sbasis in self.components:
            for obasis in other.components:
                if (obasis, sbasis) in fmodule.basis_changes:
                    other.comp(sbasis, from_basis=obasis)
                    return sbasis
                if (sbasis, obasis) in fmodule.basis_changes:
                    self.comp(obasis, from_basis=sbasis)
                    return obasis
        #
        # 3/ Search for a common basis via two component transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at two
        # component transformation to get a common basis
        for sbasis in self.components:
            for obasis in other.components:
                if (sbasis, def_basis) in fmodule.basis_changes and \
                   (obasis, def_basis) in fmodule.basis_changes:
                    self.comp(def_basis, from_basis=sbasis)
                    other.comp(def_basis, from_basis=obasis)
                    return def_basis
                for basis in fmodule.known_bases:
                    if (sbasis, basis) in fmodule.basis_changes and \
                       (obasis, basis) in fmodule.basis_changes:
                        self.comp(basis, from_basis=sbasis)
                        other.comp(basis, from_basis=obasis)
                        return basis
        #
        # If this point is reached, no common basis could be found, even at 
        # the price of component transformations:
        return None
    
    def pick_a_basis(self):
        r"""
        Return a basis in which the tensor components are defined. 
        
        The free module's default basis is privileged. 

        OUTPUT:
        
        - instance of 
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis` 
          representing the basis 

        """
        if self.fmodule.def_basis in self.components:
            return self.fmodule.def_basis  # the default basis is privileged
        else:
            # a basis is picked arbitrarily:
            return self.components.items()[0][0]  

    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a tensor or 0
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other`` and False otherwise
        
        """
        if self.tensor_rank == 0:
            raise NotImplementedError("Scalar comparison not implemented.")
        if isinstance(other, (int, Integer)): # other should be 0
            if other == 0:
                return self.is_zero()
            else:
                return False
        elif not isinstance(other, FreeModuleTensor):
            return False
        else: # other is another tensor
            if other.fmodule != self.fmodule:
                return False
            if other.tensor_type != self.tensor_type:
                return False
            basis = self.common_basis(other)
            if basis is None:
                raise ValueError("No common basis for the comparison.")
            return bool(self.components[basis] == other.components[basis])

    def __ne__(self, other):
        r"""
        Inequality operator. 
        
        INPUT:
        
        - ``other`` -- a tensor or 0
        
        OUTPUT:
        
        - True if ``self`` is different from ``other`` and False otherwise
        
        """
        return not self.__eq__(other)

    def __pos__(self):
        r"""
        Unary plus operator. 
        
        OUTPUT:
        
        - an exact copy of ``self``
    
        """
        result = self._new_instance()
        for basis in self.components:
            result.components[basis] = + self.components[basis]
        if self.name is not None:
            result.name = '+' + self.name 
        if self.latex_name is not None:
            result.latex_name = '+' + self.latex_name
        return result

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - the tensor `-T`, where `T` is ``self``
    
        """
        result = self._new_instance()
        for basis in self.components:
            result.components[basis] = - self.components[basis]
        if self.name is not None:
            result.name = '-' + self.name 
        if self.latex_name is not None:
            result.latex_name = '-' + self.latex_name
        return result

    ######### ModuleElement arithmetic operators ########
    
    def _add_(self, other):
        r"""
        Tensor addition. 
        
        INPUT:
        
        - ``other`` -- a tensor, of the same type as ``self``
        
        OUPUT:
        
        - the tensor resulting from the addition of ``self`` and ``other``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, FreeModuleTensor):
            raise TypeError("For the addition, other must be a tensor.")
        if other.tensor_type != self.tensor_type:
            raise TypeError("The two tensors are not of the same type.")
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("No common basis for the addition.")
        comp_result = self.components[basis] + other.components[basis]
        result = self.fmodule.tensor_from_comp(self.tensor_type, comp_result)
        if self.name is not None and other.name is not None:
            result.name = self.name + '+' + other.name
        if self.latex_name is not None and other.latex_name is not None:
            result.latex_name = self.latex_name + '+' + other.latex_name
        return result

    def _sub_(self, other):
        r"""
        Tensor subtraction. 
        
        INPUT:
        
        - ``other`` -- a tensor, of the same type as ``self``
        
        OUPUT:
        
        - the tensor resulting from the subtraction of ``other`` from ``self``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, FreeModuleTensor):
            raise TypeError("For the subtraction, other must be a tensor.")
        if other.tensor_type != self.tensor_type:
            raise TypeError("The two tensors are not of the same type.")
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("No common basis for the subtraction.")
        comp_result = self.components[basis] - other.components[basis]
        result = self.fmodule.tensor_from_comp(self.tensor_type, comp_result)
        if self.name is not None and other.name is not None:
            result.name = self.name + '-' + other.name
        if self.latex_name is not None and other.latex_name is not None:
            result.latex_name = self.latex_name + '-' + other.latex_name
        return result

    def _rmul_(self, other):
        r"""
        Multiplication on the left by ``other``. 
        
        """
        if isinstance(other, FreeModuleTensor):
            raise NotImplementedError("Left tensor product not implemented.")
        # Left multiplication by a scalar: 
        result = self._new_instance()
        for basis in self.components:
            result.components[basis] = other * self.components[basis]
        return result

    ######### End of ModuleElement arithmetic operators ########
    
    def __radd__(self, other):
        r"""
        Addition on the left with ``other``. 
        
        """
        return self.__add__(other)

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``. 
        
        """
        return (-self).__add__(other)

    def __mul__(self, other):
        r"""
        Tensor product. 
        """
        from format_utilities import format_mul_txt, format_mul_latex
        if isinstance(other, FreeModuleTensor):
            basis = self.common_basis(other)
            if basis is None:
                raise ValueError("No common basis for the tensor product.")
            comp_prov = self.components[basis] * other.components[basis]
            # Reordering of the contravariant and covariant indices:
            k1, l1 = self.tensor_type
            k2, l2 = other.tensor_type
            if l1 != 0:
                comp_result = comp_prov.swap_adjacent_indices(k1, 
                                                              self.tensor_rank, 
                                                              self.tensor_rank+k2)
            else:
                comp_result = comp_prov  # no reordering is necessary
            result = self.fmodule.tensor_from_comp((k1+k2, l1+l2), comp_result)
            result.name = format_mul_txt(self.name, '*', other.name)
            result.latex_name = format_mul_latex(self.latex_name, r'\otimes ', 
                                                 other.latex_name)
            return result
        else:
            # multiplication by a scalar: 
            result = self._new_instance()
            for basis in self.components:
                result.components[basis] = other * self.components[basis]
            return result


    def __div__(self, other):
        r"""
        Division (by a scalar)
        """
        result = self._new_instance()
        for basis in self.components:
            result.components[basis] = self.components[basis] / other
        return result
        

    def __call__(self, *args):
        r"""
        The tensor acting on linear forms and module elements as a multilinear 
        map.
        
        INPUT:
        
        - ``*args`` -- list of k linear forms and l module elements, ``self`` 
          being a tensor of type (k,l). 
          
        """
        from free_module_alt_form import FreeModuleLinForm
        # Consistency checks:
        p = len(args)
        if p != self.tensor_rank:
            raise TypeError(str(self.tensor_rank) + 
                            " arguments must be provided.")
        for i in range(self.tensor_type[0]):
            if not isinstance(args[i], FreeModuleLinForm):
                raise TypeError("The argument no. " + str(i+1) + 
                                " must be a linear form.")
        for i in range(self.tensor_type[0],p):
            if not isinstance(args[i], FiniteFreeModuleElement):
                raise TypeError("The argument no. " + str(i+1) + 
                                " must be a module element.")
        fmodule = self.fmodule
        # Search for a common basis
        basis = None
        # First try with the module's default basis 
        def_basis = fmodule.def_basis
        if def_basis in self.components:
            basis = def_basis
            for arg in args:
                if def_basis not in arg.components:
                    basis = None
                    break
        if basis is None:
            # Search for another basis:
            for bas in self.components:
                basis = bas
                for arg in args:
                    if bas not in arg.components:
                        basis = None
                        break
                if basis is not None: # common basis found ! 
                    break
        if basis is None:
            raise ValueError("No common basis for the components.")
        t = self.components[basis]
        v = [args[i].components[basis] for i in range(p)]
        
        res = 0
        for ind in t.index_generator():
            prod = t[[ind]]
            for i in range(p):
                prod *= v[i][[ind[i]]]
            res += prod
        # Name of the output:
        if hasattr(res, 'name'): 
            res_name = None
            if self.name is not None:
                res_name = self.name + "("
                for i in range(p-1):
                    if args[i].name is not None:
                        res_name += args[i].name + ","
                    else:
                        res_name = None
                        break
                if res_name is not None:
                    if args[p-1].name is not None:
                        res_name += args[p-1].name + ")"
                    else:
                        res_name = None
            res.name = res_name       
        # LaTeX symbol of the output:
        if hasattr(res, 'latex_name'): 
            res_latex = None
            if self.latex_name is not None:
                res_latex = self.latex_name + r"\left("
                for i in range(p-1):
                    if args[i].latex_name is not None:
                        res_latex += args[i].latex_name + ","
                    else:
                        res_latex = None
                        break
                if res_latex is not None:
                    if args[p-1].latex_name is not None:
                        res_latex += args[p-1].latex_name + r"\right)"
                    else:
                        res_latex = None
            res.latex_name = res_latex
        return res

    def self_contract(self, pos1, pos2):
        r""" 
        Contraction on two slots of the tensor. 
        
        INPUT:
            
        - ``pos1`` -- position of the first index for the contraction, with the
          convention ``pos1=0`` for the first slot
        - ``pos2`` -- position of the second index for the contraction, with 
          the same convention as for ``pos1``. 
          
        OUTPUT:
        
        - tensor resulting from the (pos1, pos2) contraction
       
        EXAMPLES:
        
        Contraction on the two slots of a type-(1,1) tensor::

            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e') ; e
            basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
            sage: a = M.tensor((1,1), name='a') ; a
            endomorphism a on the rank-3 free module M over the Integer Ring
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: a.self_contract(0,1)  # contraction of slot 0 with slot 1
            15
            sage: a.self_contract(1,0)  # the order of the slots does not matter
            15

        The contraction on two slots having the same tensor type cannot occur::
        
            sage: b =  M.tensor((2,0), name='b') ; b
            type-(2,0) tensor b on the rank-3 free module M over the Integer Ring
            sage: b[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b.self_contract(0,1)
            Traceback (most recent call last):
            ...
            IndexError: Contraction on two contravariant indices is not allowed.

        The contraction either preserves or destroys the symmetries::
        
            sage: b = M.alternating_form(2, 'b') ; b
            alternating form b of degree 2 on the rank-3 free module M over the Integer Ring
            sage: b[0,1], b[0,2], b[1,2] = 3, 2, 1
            sage: t = a*b ; t
            type-(1,3) tensor a*b on the rank-3 free module M over the Integer Ring
            sage: # by construction, t is a tensor field antisymmetric w.r.t. its last two slots:
            sage: t.symmetries()
            no symmetry;  antisymmetry: (2, 3)
            sage: s = t.self_contract(0,1) ; s   # contraction on the first two slots
            alternating form of degree 2 on the rank-3 free module M over the Integer Ring
            sage: s.symmetries()    # the antisymmetry is preserved
            no symmetry;  antisymmetry: (0, 1)
            sage: s[:]
            [  0  45  30]
            [-45   0  15]
            [-30 -15   0]
            sage: s == 15*b  # check
            True
            sage: s = t.self_contract(0,2) ; s   # contraction on the first and third slots
            type-(0,2) tensor on the rank-3 free module M over the Integer Ring
            sage: s.symmetries()  # the antisymmetry has been destroyed by the above contraction:
            no symmetry;  no antisymmetry
            sage: s[:]  # indeed:
            [-26  -4   6]
            [-31  -2   9]
            [-36   0  12]
            sage: s[:] == matrix( [[sum(t[k,i,k,j] for k in M.irange()) for j in M.irange()] for i in M.irange()] )  # check
            True
            
        """
        # The indices at pos1 and pos2 must be of different types: 
        k_con = self.tensor_type[0]
        l_cov = self.tensor_type[1]
        if pos1 < k_con and pos2 < k_con:
            raise IndexError("Contraction on two contravariant indices is " +
                             "not allowed.")
        if pos1 >= k_con and pos2 >= k_con:
            raise IndexError("Contraction on two covariant indices is " +
                             "not allowed.")
        # Frame selection for the computation: 
        if self.fmodule.def_basis in self.components:
            basis = self.fmodule.def_basis
        else: # a basis is picked arbitrarily:
            basis = self.pick_a_basis()     
        resu_comp = self.components[basis].self_contract(pos1, pos2)
        if self.tensor_rank == 2:  # result is a scalar
            return resu_comp
        else:
            return self.fmodule.tensor_from_comp((k_con-1, l_cov-1), resu_comp)

    def contract(self, *args):
        r""" 
        Contraction with another tensor. 
        
        INPUT:
            
        - ``pos1`` -- position of the first index (in ``self``) for the 
          contraction; if not given, the last index position is assumed
        - ``other`` -- the tensor to contract with
        - ``pos2`` -- position of the second index (in ``other``) for the 
          contraction; if not given, the first index position is assumed
          
        OUTPUT:
        
        - tensor resulting from the (pos1, pos2) contraction of ``self`` with 
          ``other``
       
        EXAMPLES:
        
        Contraction of a tensor of type (0,1) with a tensor of type (1,0)::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.linear_form()  # tensor of type (0,1) is a linear form
            sage: a[:] = [-3,2,1]
            sage: b = M([2,5,-2])  # tensor of type (1,0) is a module element
            sage: s = a.contract(b) ; s
            2
            sage: s in M.base_ring()
            True
            sage: s == a[0]*b[0] + a[1]*b[1] + a[2]*b[2]  # check of the computation
            True

        The positions of the contraction indices can be set explicitely::
        
            sage: s == a.contract(0, b, 0)
            True
            sage: s == a.contract(0, b)
            True
            sage: s == a.contract(b, 0)
            True
           
        In the present case, performing the contraction is identical to 
        applying the linear form to the module element::
        
            sage: a.contract(b) == a(b)
            True

        or to applying the module element, considered as a tensor of type (1,0),
        to the linear form::
        
            sage: a.contract(b) == b(a)
            True

        We have also::

            sage: a.contract(b) == b.contract(a)
            True

        Contraction of a tensor of type (1,1) with a tensor of type (1,0)::
        
            sage: a = M.endomorphism()  # tensor of type (1,1)
            sage: a[:] = [[-1,2,3],[4,-5,6],[7,8,9]]
            sage: s = a.contract(b) ; s
            element of the rank-3 free module M over the Integer Ring
            sage: s.view()
            2 e_0 - 29 e_1 + 36 e_2
    
        Since the index positions have not been specified, the contraction
        takes place on the last position of a (i.e. no. 1) and the first
        position of b (i.e. no. 0)::
        
            sage: a.contract(b) == a.contract(1, b, 0)
            True
            sage: a.contract(b) == b.contract(0, a, 1)
            True
            sage: a.contract(b) == b.contract(a, 1) 
            True
            
        Contraction is possible only between a contravariant index and a 
        covariant one::
        
            sage: a.contract(0, b)
            Traceback (most recent call last):
            ...
            TypeError: Contraction not possible: the two index positions are both contravariant.

        In the present case, performing the contraction is identical to 
        applying the endomorphism to the module element::
        
            sage: a.contract(b) == a(b)
            True

        Contraction of a tensor of type (2,1) with a tensor of type (0,2)::
        
            sage: a = a*b ; a
            type-(2,1) tensor on the rank-3 free module M over the Integer Ring
            sage: b = M.tensor((0,2))
            sage: b[:] = [[-2,3,1], [0,-2,3], [4,-7,6]]
            sage: s = a.contract(1, b, 1) ; s
            type-(1,2) tensor on the rank-3 free module M over the Integer Ring
            sage: s[:]
            [[[-9, 16, 39], [18, -32, -78], [27, -48, -117]],
             [[36, -64, -156], [-45, 80, 195], [54, -96, -234]],
             [[63, -112, -273], [72, -128, -312], [81, -144, -351]]]
            sage: # check of the computation:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == a[i,0,j]*b[k,0]+a[i,1,j]*b[k,1]+a[i,2,j]*b[k,2],
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True

        The two tensors do not have to be defined on the same basis for the 
        contraction to take place, reflecting the fact that the contraction is 
        basis-independent::
        
            sage: A = M.automorphism()
            sage: A[:] =  [[0,0,1], [1,0,0], [0,-1,0]]
            sage: f = e.new_basis(A, 'f')
            sage: b.comp(f)[:]  # forces the computation of b's components w.r.t. basis f
            [-2 -3  0]
            [ 7  6 -4]
            [ 3 -1 -2]
            sage: b.del_other_comp(f)  # deletes components w.r.t. basis e
            sage: b.components.keys()  # indeed:
            [basis (f_0,f_1,f_2) on the rank-3 free module M over the Integer Ring]
            sage: a.components.keys()  # while a is known only in basis e:
            [basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring]
            sage: s1 = a.contract(1, b, 1) ; s1  # yet the computation is possible
            type-(1,2) tensor on the rank-3 free module M over the Integer Ring
            sage: s1 == s  # ... and yields the same result as previously:
            True

        """
        nargs = len(args)
        if nargs == 1:
            pos1 = self.tensor_rank - 1
            other = args[0]
            pos2 = 0
        elif nargs == 2:
            if isinstance(args[0], FreeModuleTensor):
                pos1 = self.tensor_rank - 1
                other = args[0]
                pos2 = args[1]
            else:
                pos1 = args[0]
                other = args[1]
                pos2 = 0
        elif nargs == 3:
            pos1 = args[0]
            other = args[1]
            pos2 = args[2]
        else:
            raise TypeError("Wrong number of arguments in contract(): " + 
                str(nargs) + 
                " arguments provided, while between 1 and 3 are expected.")
        if not isinstance(other, FreeModuleTensor):
            raise TypeError("For the contraction, other must be a tensor " + 
                            "field.")
        k1, l1 = self.tensor_type
        k2, l2 = other.tensor_type
        if pos1 < k1 and pos2 < k2:
            raise TypeError("Contraction not possible: the two index " + 
                            "positions are both contravariant.")
        if pos1 >= k1 and pos2 >= k2:
            raise TypeError("Contraction not possible: the two index " + 
                            "positions are both covavariant.")
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("No common basis for the contraction.")
        cmp_res = self.components[basis].contract(pos1, 
                                            other.components[basis], pos2)
        # reordering of the indices to have all contravariant indices first:
        if k2 > 1:
            if pos1 < k1:
                cmp_res = cmp_res.swap_adjacent_indices(k1-1, k1+l1-1, k1+l1+k2-1)
            else:
                cmp_res = cmp_res.swap_adjacent_indices(k1, k1+l1-1, k1+l1+k2-2)
        type_res = (k1+k2-1, l1+l2-1)
        if type_res == (0, 0):
            return cmp_res  # scalar case
        else:
            return self.fmodule.tensor_from_comp(type_res, cmp_res)


    def symmetrize(self, pos=None, basis=None):
        r"""
        Symmetrization over some arguments.
        
        INPUT:
        
        - ``pos`` -- (default: None) list of argument positions involved in the 
          symmetrization (with the convention position=0 for the first 
          argument); if none, the symmetrization is performed over all the 
          arguments
        - ``basis`` -- (default: None) module basis with respect to which the 
          component computation is to be performed; if none, the module's 
          default basis is used if the tensor field has already components
          in it; otherwise another basis w.r.t. which the tensor has 
          components will be picked
                  
        OUTPUT:
        
        - the symmetrized tensor (instance of :class:`FreeModuleTensor`)
          
        EXAMPLES:
        
        Symmetrization of a tensor of type (2,0)::
        
            sage: M = FiniteFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((2,0))
            sage: t[:] = [[2,1,-3],[0,-4,5],[-1,4,2]]
            sage: s = t.symmetrize() ; s
            type-(2,0) tensor on the rank-3 free module M over the Rational Field
            sage: t[:], s[:]
            (
            [ 2  1 -3]  [  2 1/2  -2]
            [ 0 -4  5]  [1/2  -4 9/2]
            [-1  4  2], [ -2 9/2   2]
            )
            sage: s.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         print s[i,j] == 1/2*(t[i,j]+t[j,i]),
            ....:         
            True True True True True True True True True
    
        Symmetrization of a tensor of type (0,3) on the first two arguments::
        
            sage: t = M.tensor((0,3))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]], [[10,-11,12], [13,14,-15], [16,17,18]], [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.symmetrize((0,1)) ; s  # (0,1) = the first two arguments
            type-(0,3) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: s[:]
            [[[1, 2, 3], [3, -3, 9], [13, -6, -15]],
             [[3, -3, 9], [13, 14, -15], [-3, 20, 21]],
             [[13, -6, -15], [-3, 20, 21], [25, 26, -27]]]
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/2*(t[i,j,k]+t[j,i,k]), 
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.symmetrize((0,1)) == s  # another test
            True

        Symmetrization of a tensor of type (0,3) on the first and last arguments::

            sage: s = t.symmetrize((0,2)) ; s  # (0,2) = first and last arguments
            type-(0,3) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            symmetry: (0, 2);  no antisymmetry
            sage: s[:]
            [[[1, 6, 11], [-4, 9, -8], [7, 12, 8]],
             [[6, -11, -4], [9, 14, 4], [12, 17, 22]],
             [[11, -4, -21], [-8, 4, 24], [8, 22, -27]]]
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/2*(t[i,j,k]+t[k,j,i]),
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.symmetrize((0,2)) == s  # another test
            True

        Symmetrization of a tensor of type (0,3) on the last two arguments::
        
            sage: s = t.symmetrize((1,2)) ; s  # (1,2) = the last two arguments
            type-(0,3) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            symmetry: (1, 2);  no antisymmetry
            sage: s[:]
            [[[1, -1, 5], [-1, 5, 7], [5, 7, -9]],
             [[10, 1, 14], [1, 14, 1], [14, 1, 18]],
             [[19, -21, 2], [-21, 23, 25], [2, 25, -27]]]
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/2*(t[i,j,k]+t[i,k,j]),  
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.symmetrize((1,2)) == s  # another test
            True
    
        Full symmetrization of a tensor of type (0,3)::
        
            sage: s = t.symmetrize() ; s
            type-(0,3) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            symmetry: (0, 1, 2);  no antisymmetry
            sage: s[:]
            [[[1, 8/3, 29/3], [8/3, 7/3, 0], [29/3, 0, -5/3]],
             [[8/3, 7/3, 0], [7/3, 14, 25/3], [0, 25/3, 68/3]],
             [[29/3, 0, -5/3], [0, 25/3, 68/3], [-5/3, 68/3, -27]]]
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/6*(t[i,j,k]+t[i,k,j]+t[j,k,i]+t[j,i,k]+t[k,i,j]+t[k,j,i]),
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.symmetrize() == s  # another test
            True

        Symmetrization can be performed only on arguments on the same type::
        
            sage: t = M.tensor((1,2))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]], [[10,-11,12], [13,14,-15], [16,17,18]], [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.symmetrize((0,1)) 
            Traceback (most recent call last):
            ...
            TypeError: 0 is a contravariant position, while 1 is a covariant position; 
            symmetrization is meaningfull only on tensor arguments of the same type.
            sage: s = t.symmetrize((1,2)) # OK: both 1 and 2 are covariant positions

        The order of positions does not matter::
        
            sage: t.symmetrize((2,1)) == t.symmetrize((1,2))
            True
            
        """
        if pos is None:
            pos = range(self.tensor_rank)
        # check whether the symmetrization is possible:
        pos_cov = self.tensor_type[0]   # first covariant position 
        pos0 = pos[0]
        if pos0 < pos_cov:  # pos0 is a contravariant position
            for k in range(1,len(pos)):
                if pos[k] >= pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a contravariant position, while " + 
                        str(pos[k]) + " is a covariant position; \n"
                        "symmetrization is meaningfull only on tensor " + 
                        "arguments of the same type.")
        else:  # pos0 is a covariant position
            for k in range(1,len(pos)):
                if pos[k] < pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a covariant position, while " + \
                        str(pos[k]) + " is a contravariant position; \n"
                        "symmetrization is meaningfull only on tensor " + 
                        "arguments of the same type.")                
        if basis is None:
            basis = self.pick_a_basis()
        res_comp = self.components[basis].symmetrize(pos)
        return self.fmodule.tensor_from_comp(self.tensor_type, res_comp)

        
    def antisymmetrize(self, pos=None, basis=None):
        r"""
        Antisymmetrization over some arguments.
        
        INPUT:
        
        - ``pos`` -- (default: None) list of argument positions involved in the 
          antisymmetrization (with the convention position=0 for the first 
          argument); if none, the antisymmetrization is performed over all the 
          arguments
        - ``basis`` -- (default: None) module basis with respect to which the 
          component computation is to be performed; if none, the module's 
          default basis is used if the tensor field has already components
          in it; otherwise another basis w.r.t. which the tensor has 
          components will be picked
                  
        OUTPUT:
        
        - the antisymmetrized tensor (instance of :class:`FreeModuleTensor`)
          
        EXAMPLES:
                
        Antisymmetrization of a tensor of type (2,0)::
        
            sage: M = FiniteFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: t = M.tensor((2,0))
            sage: t[:] = [[1,-2,3], [4,5,6], [7,8,-9]]
            sage: s = t.antisymmetrize() ; s
            type-(2,0) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 1)
            sage: t[:], s[:]
            (
            [ 1 -2  3]  [ 0 -3 -2]
            [ 4  5  6]  [ 3  0 -1]
            [ 7  8 -9], [ 2  1  0]
            )
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         print s[i,j] == 1/2*(t[i,j]-t[j,i]),
            ....:         
            True True True True True True True True True
            sage: s.antisymmetrize() == s  # another test
            True
            sage: t.antisymmetrize() == t.antisymmetrize((0,1))
            True

        Antisymmetrization of a tensor of type (0,3) over the first two 
        arguments::

            sage: t = M.tensor((0,3))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]], [[10,-11,12], [13,14,-15], [16,17,18]], [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.antisymmetrize((0,1)) ; s  # (0,1) = the first two arguments
            type-(0,3) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 1)
            sage: s[:]
            [[[0, 0, 0], [-7, 8, -3], [-6, 14, 6]],
             [[7, -8, 3], [0, 0, 0], [19, -3, -3]],
             [[6, -14, -6], [-19, 3, 3], [0, 0, 0]]]
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/2*(t[i,j,k]-t[j,i,k]),
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.antisymmetrize((0,1)) == s  # another test
            True
            sage: s.symmetrize((0,1)) == 0  # of course
            True
        
        Antisymmetrization of a tensor of type (0,3) over the first and last 
        arguments::

            sage: s = t.antisymmetrize((0,2)) ; s  # (0,2) = first and last arguments
            type-(0,3) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 2)
            sage: s[:]
            [[[0, -4, -8], [0, -4, 14], [0, -4, -17]],
             [[4, 0, 16], [4, 0, -19], [4, 0, -4]],
             [[8, -16, 0], [-14, 19, 0], [17, 4, 0]]]
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/2*(t[i,j,k]-t[k,j,i]),
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.antisymmetrize((0,2)) == s  # another test
            True
            sage: s.symmetrize((0,2)) == 0  # of course
            True
            sage: s.symmetrize((0,1)) == 0  # no reason for this to hold
            False
        
        Antisymmetrization of a tensor of type (0,3) over the last two 
        arguments::

            sage: s = t.antisymmetrize((1,2)) ; s  # (1,2) = the last two arguments
            type-(0,3) tensor on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (1, 2)
            sage: s[:]
            [[[0, 3, -2], [-3, 0, -1], [2, 1, 0]],
             [[0, -12, -2], [12, 0, -16], [2, 16, 0]],
             [[0, 1, -23], [-1, 0, -1], [23, 1, 0]]]
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/2*(t[i,j,k]-t[i,k,j]),
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.antisymmetrize((1,2)) == s  # another test
            True
            sage: s.symmetrize((1,2)) == 0  # of course
            True

        Full antisymmetrization of a tensor of type (0,3)::
        
            sage: s = t.antisymmetrize() ; s
            alternating form of degree 3 on the rank-3 free module M over the Rational Field
            sage: s.symmetries()
            no symmetry;  antisymmetry: (0, 1, 2)
            sage: s[:]
            [[[0, 0, 0], [0, 0, 2/3], [0, -2/3, 0]],
             [[0, 0, -2/3], [0, 0, 0], [2/3, 0, 0]],
             [[0, 2/3, 0], [-2/3, 0, 0], [0, 0, 0]]]
            sage: # Check:
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         for k in range(3):
            ....:             print s[i,j,k] == 1/6*(t[i,j,k]-t[i,k,j]+t[j,k,i]-t[j,i,k]+t[k,i,j]-t[k,j,i]),
            ....:             
            True True True True True True True True True True True True True True True True True True True True True True True True True True True
            sage: s.antisymmetrize() == s  # another test
            True
            sage: s.symmetrize((0,1)) == 0  # of course
            True
            sage: s.symmetrize((0,2)) == 0  # of course
            True
            sage: s.symmetrize((1,2)) == 0  # of course
            True
            sage: t.antisymmetrize() == t.antisymmetrize((0,1,2))
            True

        Antisymmetrization can be performed only on arguments on the same type::
        
            sage: t = M.tensor((1,2))
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]], [[10,-11,12], [13,14,-15], [16,17,18]], [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: s = t.antisymmetrize((0,1)) 
            Traceback (most recent call last):
            ...
            TypeError: 0 is a contravariant position, while 1 is a covariant position; 
            antisymmetrization is meaningfull only on tensor arguments of the same type.
            sage: s = t.antisymmetrize((1,2)) # OK: both 1 and 2 are covariant positions

        The order of positions does not matter::
        
            sage: t.antisymmetrize((2,1)) == t.antisymmetrize((1,2))
            True
            
        """
        if pos is None:
            pos = range(self.tensor_rank)
        # check whether the antisymmetrization is possible:
        pos_cov = self.tensor_type[0]   # first covariant position 
        pos0 = pos[0]
        if pos0 < pos_cov:  # pos0 is a contravariant position
            for k in range(1,len(pos)):
                if pos[k] >= pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a contravariant position, while " + 
                        str(pos[k]) + " is a covariant position; \n"
                        "antisymmetrization is meaningfull only on tensor " + 
                        "arguments of the same type.")
        else:  # pos0 is a covariant position
            for k in range(1,len(pos)):
                if pos[k] < pos_cov:
                    raise TypeError(
                        str(pos[0]) + " is a covariant position, while " + \
                        str(pos[k]) + " is a contravariant position; \n"
                        "antisymmetrization is meaningfull only on tensor " + 
                        "arguments of the same type.")                
        if basis is None:
            basis = self.pick_a_basis()
        res_comp = self.components[basis].antisymmetrize(pos)
        return self.fmodule.tensor_from_comp(self.tensor_type, res_comp)
        


#******************************************************************************

# From sage/modules/module.pyx:
#-----------------------------
### The Element should also implement _rmul_ (or _lmul_)
#
# class MyElement(sage.structure.element.ModuleElement):
#     def _rmul_(self, c):
#         ...


class FiniteFreeModuleElement(FreeModuleTensor):
    r"""
    Element of a free module of finite rank over a commutative ring.
    
    The class :class:`FiniteFreeModuleElement` inherits from 
    :class:`FreeModuleTensor` because the elements of a free module `M` of 
    finite rank over a commutative ring `R` are identified with tensors of 
    type (1,0) on `M` via the canonical map
    
    .. MATH::

        \begin{array}{lllllll}
        \Phi: & M & \longrightarrow & M^{**}  &  &   & \\
              & v & \longmapsto & \bar v : & M^* & \longrightarrow & R \\
              &   &             &          & a   & \longmapsto & a(v) 
        \end{array}
    
    Note that for free modules of finite rank, this map is actually an 
    isomorphism, enabling the canonical identification: `M^{**}= M`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` (must be an 
      instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the element
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the element; 
      if none is provided, the LaTeX symbol is set to ``name``
    
    EXAMPLES:
    
    Let us consider a rank-3 free module over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e') ; e
        basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
        
    There are four ways to construct an element of the free module M: the first 
    one (recommended) is via the operator __call__ acting on the free module::
    
        sage: v = M([2,0,-1], basis=e, name='v') ; v
        element v of the rank-3 free module M over the Integer Ring
        sage: v.view()  # expansion on the default basis (e)
        v = 2 e_0 - e_2
        sage: v.parent() is M
        True

    The second way is by a direct call to the class constructor::
    
        sage: from sage.tensor.modules.free_module_tensor import FiniteFreeModuleElement
        sage: v2 = FiniteFreeModuleElement(M, name='v')
        sage: v2[0], v2[2] = 2, -1 # setting the nonzero components in the default basis (e)
        sage: v2
        element v of the rank-3 free module M over the Integer Ring
        sage: v2.view()
        v = 2 e_0 - e_2
        sage: v2 == v
        True

    The third way is to construct a tensor of type (1,0) on `M` (cf. the
    canonical identification `M^{**}=M` recalled above)::
    
        sage: v3 = M.tensor((1,0), name='v')
        sage: v3[0], v3[2] = 2, -1 ; v3
        element v of the rank-3 free module M over the Integer Ring
        sage: v3.view()
        v = 2 e_0 - e_2
        sage: v3 == v
        True
        
    Finally, the fourth way is via some linear combination of the basis 
    elements::
    
        sage: v4 = 2*e[0] - e[2]
        sage: v4.set_name('v') ; v4 # in this case, the name has to be set separately
        element v of the rank-3 free module M over the Integer Ring
        sage: v4.view()
        v = 2 e_0 - e_2
        sage: v4 == v
        True
    
    The canonical identification `M^{**}=M` is implemented by letting the
    module elements act on linear forms, providing the same result as the
    reverse operation (cf. the map `\Phi` defined above)::
    
        sage: a = M.linear_form(name='a')
        sage: a[:] = (2, 1, -3) ; a
        linear form a on the rank-3 free module M over the Integer Ring
        sage: v(a)
        7
        sage: a(v)
        7
        sage: a(v) == v(a)
        True
    
    ARITHMETIC EXAMPLES:
    
    Addition::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e') ; e
        basis (e_0,e_1,e_2) on the rank-3 free module M over the Integer Ring
        sage: a = M([0,1,3], name='a') ; a
        element a of the rank-3 free module M over the Integer Ring
        sage: a.view()
        a = e_1 + 3 e_2
        sage: b = M([2,-2,1], name='b') ; b
        element b of the rank-3 free module M over the Integer Ring
        sage: b.view()
        b = 2 e_0 - 2 e_1 + e_2
        sage: s = a + b ; s
        element a+b of the rank-3 free module M over the Integer Ring
        sage: s.view()
        a+b = 2 e_0 - e_1 + 4 e_2
        sage: # Test of the addition:
        sage: for i in M.irange(): print s[i] == a[i] + b[i],
        True True True

    Subtraction::
    
        sage: s = a - b ; s
        element a-b of the rank-3 free module M over the Integer Ring
        sage: s.view()
        a-b = -2 e_0 + 3 e_1 + 2 e_2
        sage: # Test of the substraction:
        sage: for i in M.irange(): print s[i] == a[i] - b[i],
        True True True
        
    Multiplication by a scalar::
    
        sage: s = 2*a ; s
        element of the rank-3 free module M over the Integer Ring
        sage: s.view()
        2 e_1 + 6 e_2
        sage: a.view()
        a = e_1 + 3 e_2

    Tensor product::
    
        sage: s = a*b ; s
        type-(2,0) tensor a*b on the rank-3 free module M over the Integer Ring
        sage: s.symmetries()
        no symmetry;  no antisymmetry
        sage: s[:]
        [ 0  0  0]
        [ 2 -2  1]
        [ 6 -6  3]
        sage: s = a*s ; s
        type-(3,0) tensor a*a*b on the rank-3 free module M over the Integer Ring
        sage: s[:]
        [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
         [[0, 0, 0], [2, -2, 1], [6, -6, 3]],
         [[0, 0, 0], [6, -6, 3], [18, -18, 9]]]

    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (1,0), name=name, 
                                  latex_name=latex_name)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "element "
        if self.name is not None:
            description += self.name + " " 
        description += "of the " + str(self.fmodule)
        return description

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
              
        This method, which is already implemented in 
        :meth:`FreeModuleTensor._new_comp`, is redefined here for efficiency
        """
        fmodule = self.fmodule  # the base free module
        return Components(fmodule.ring, basis, 1, start_index=fmodule.sindex,
                          output_formatter=fmodule.output_formatter)


    def _new_instance(self):
        r"""
        Create a :class:`FiniteFreeModuleElement` instance.
        
        """
        return FiniteFreeModuleElement(self.fmodule)

        

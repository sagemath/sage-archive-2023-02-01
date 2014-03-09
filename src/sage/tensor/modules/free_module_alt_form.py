r"""
Alternating forms on free modules

The class :class:`FreeModuleAltForm` implement alternating forms on a free 
module of finite rank over a commutative ring. 

It is a subclass of 
:class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`, alternating 
forms being a special type of tensors. 

A subclass of :class:`FreeModuleAltForm` is :class:`FreeModuleLinForm` for 
alternating forms of degree 1, i.e. linear forms. 

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

from free_module_tensor import FreeModuleTensor, FiniteFreeModuleElement
from comp import Components, CompFullyAntiSym

class FreeModuleAltForm(FreeModuleTensor):
    r"""
    Alternating form over a free module `M`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` of finite rank over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``degree`` -- the degree of the alternating form (i.e. its tensor rank)
    - ``name`` -- (default: None) name given to the alternating form
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the alternating 
      form; if none is provided, the LaTeX symbol is set to ``name``

    """
    def __init__(self, fmodule, degree, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (0,degree), name=name, 
                                  latex_name=latex_name, antisym=range(degree))
        FreeModuleAltForm._init_derived(self) # initialization of derived quantities

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "alternating form "
        if self.name is not None:
            description += self.name + " " 
        description += "of degree " + str(self.tensor_rank) + " on the " + \
                       str(self.fmodule)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        FreeModuleTensor._init_derived(self)  

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        FreeModuleTensor._del_derived(self)  

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 

        This method, which is already implemented in 
        :meth:`FreeModuleTensor._new_comp`, is redefined here for efficiency
        """
        fmodule = self.fmodule  # the base free module
        if self.tensor_rank == 1: 
            return Components(fmodule.ring, basis, 1,
                              start_index=fmodule.sindex,
                              output_formatter=fmodule.output_formatter)
        else:
            return CompFullyAntiSym(fmodule.ring, basis, self.tensor_rank, 
                                    start_index=fmodule.sindex,
                                    output_formatter=fmodule.output_formatter)

    def degree(self):
        r"""
        Return the degree of the alternating form. 
        """
        return self.tensor_rank


    def view(self, basis=None, format_spec=None):
        r"""
        Displays the alternating form in terms of its expansion onto a given 
        cobasis.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        INPUT:        
        
        - ``basis`` -- (default: None) basis of the free module with respect to 
          which the alternating form is expanded; if none is provided, the 
          module's default basis is assumed
        - ``format_spec`` -- (default: None) format specification passed to 
          ``self.fmodule.output_formatter`` to format the output.
          
        EXAMPLES:
        
        Display of an alternating form of degree 1 (linear form) on a rank-3 
        free module::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.linear_form('a', latex_name=r'\alpha')
            sage: a[:] = [1,-3,4]
            sage: a.view()
            a = e^0 - 3 e^1 + 4 e^2
            sage: latex(a.view())  # display in the notebook
            \alpha = e^0 -3 e^1 + 4 e^2

        Display of an alternating form of degree 2 on a rank-3 free module::
        
            sage: b = M.alternating_form(2, 'b', latex_name=r'\beta')
            sage: b[0,1], b[0,2], b[1,2] = 3, 2, -1
            sage: b.view()
            b = 3 e^0/\e^1 + 2 e^0/\e^2 - e^1/\e^2
            sage: latex(b.view())  # display in the notebook
            \beta = 3 e^0\wedge e^1 + 2 e^0\wedge e^2 -e^1\wedge e^2

        Display of an alternating form of degree 3 on a rank-3 free module::
        
            sage: c = M.alternating_form(3, 'c')
            sage: c[0,1,2] = 4
            sage: c.view()
            c = 4 e^0/\e^1/\e^2
            sage: latex(c.view())
            c = 4 e^0\wedge e^1\wedge e^2

        Display of a vanishing alternating form::
        
            sage: c[0,1,2] = 0  # the only independent component set to zero
            sage: c.is_zero()
            True
            sage: c.view()
            c = 0
            sage: latex(c.view())
            c = 0
            sage: c[0,1,2] = 4  # value restored for what follows
            
        Display in a basis which is not the default one::
        
            sage: aut = M.automorphism()
            sage: aut[:] = [[0,1,0], [0,0,-1], [1,0,0]]
            sage: f = e.new_basis(aut, 'f')
            sage: a.view(f)
            a = 4 f^0 + f^1 + 3 f^2
            sage: b.view(f)
            b = -2 f^0/\f^1 - f^0/\f^2 - 3 f^1/\f^2
            sage: c.view(f)
            c = -4 f^0/\f^1/\f^2

        The output format can be set via the argument ``output_formatter`` 
        passed at the module construction::

            sage: N = FiniteFreeModule(QQ, 3, name='N', start_index=1, output_formatter=Rational.numerical_approx)
            sage: e = N.basis('e')
            sage: b = N.alternating_form(2, 'b')
            sage: b[1,2], b[1,3], b[2,3] = 1/3, 5/2, 4
            sage: b.view()  # default format (53 bits of precision)
            b = 0.333333333333333 e^1/\e^2 + 2.50000000000000 e^1/\e^3 + 4.00000000000000 e^2/\e^3
            
        The output format is then controled by the argument ``format_spec`` of
        the method :meth:`view`::
        
            sage: b.view(format_spec=10)  # 10 bits of precision
            b = 0.33 e^1/\e^2 + 2.5 e^1/\e^3 + 4.0 e^2/\e^3
    
        """
        from sage.misc.latex import latex
        from format_utilities import is_atomic, FormattedExpansion
        if basis is None:
            basis = self.fmodule.def_basis
        cobasis = basis.dual_basis()
        comp = self.comp(basis)
        terms_txt = []
        terms_latex = []
        for ind in comp.non_redundant_index_generator():
            ind_arg = ind + (format_spec,)
            coef = comp[ind_arg]
            if coef != 0:
                bases_txt = []
                bases_latex = []
                for k in range(self.tensor_rank):
                    bases_txt.append(cobasis[ind[k]].name)
                    bases_latex.append(latex(cobasis[ind[k]]))
                basis_term_txt = "/\\".join(bases_txt)    
                basis_term_latex = r"\wedge ".join(bases_latex)    
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


    def wedge(self, other):
        r"""
        Exterior product with another alternating form. 
        
        INPUT:
        
        - ``other``: another alternating form
        
        OUTPUT:
        
        - instance of :class:`FreeModuleAltForm` representing the exterior 
          product self/\\other. 
        
        EXAMPLES:
        
        Exterior product of two linear forms::
        
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.linear_form('A')
            sage: a[:] = [1,-3,4]
            sage: b = M.linear_form('B')
            sage: b[:] = [2,-1,2]
            sage: c = a.wedge(b) ; c
            alternating form A/\B of degree 2 on the rank-3 free module M over the Integer Ring
            sage: c.view()
            A/\B = 5 e^0/\e^1 - 6 e^0/\e^2 - 2 e^1/\e^2
            sage: latex(c)
            A\wedge B
            sage: latex(c.view())
            A\wedge B = 5 e^0\wedge e^1 -6 e^0\wedge e^2 -2 e^1\wedge e^2

        Test of the computation::
        
            sage: a.wedge(b) == a*b - b*a
            True

        Exterior product of a linear form and an alternating form of degree 2::
        
            sage: d = M.linear_form('D')
            sage: d[:] = [-1,2,4]
            sage: s = d.wedge(c) ; s
            alternating form D/\A/\B of degree 3 on the rank-3 free module M over the Integer Ring
            sage: s.view()
            D/\A/\B = 34 e^0/\e^1/\e^2

        Test of the computation::
  
            sage: s[0,1,2] == d[0]*c[1,2] + d[1]*c[2,0] + d[2]*c[0,1]
            True
            
        Let us check that the exterior product is associative::
        
            sage: d.wedge(a.wedge(b)) == (d.wedge(a)).wedge(b)
            True
            
        and that it is graded anticommutative::
            
            sage: a.wedge(b) == - b.wedge(a)  
            True
            sage: d.wedge(c) == c.wedge(d)
            True
        
        """
        from format_utilities import is_atomic
        if not isinstance(other, FreeModuleAltForm):
            raise TypeError("The second argument for the exterior product " + 
                            "must be an alternating form.")
        if other.tensor_rank == 0:
            return other*self
        if self.tensor_rank == 0:
            return self*other
        fmodule = self.fmodule
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("No common basis for the exterior product.")
        rank_r = self.tensor_rank + other.tensor_rank
        cmp_s = self.components[basis]
        cmp_o = other.components[basis]
        cmp_r = CompFullyAntiSym(fmodule.ring, basis, rank_r, 
                                 start_index=fmodule.sindex,
                                 output_formatter=fmodule.output_formatter)
        for ind_s, val_s in cmp_s._comp.items():
            for ind_o, val_o in cmp_o._comp.items():
                ind_r = ind_s + ind_o
                if len(ind_r) == len(set(ind_r)): # all indices are different
                    cmp_r[ind_r] += val_s * val_o
        result = FreeModuleAltForm(fmodule, rank_r)
        result.components[basis] = cmp_r
        if self.name is not None and other.name is not None:
            sname = self.name
            oname = other.name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            result.name = sname + '/\\' + oname
        if self.latex_name is not None and other.latex_name is not None:
            slname = self.latex_name
            olname = other.latex_name
            if not is_atomic(slname):
                slname = '(' + slname + ')'
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            result.latex_name = slname + r'\wedge ' + olname
        return result


#******************************************************************************

class FreeModuleLinForm(FreeModuleAltForm):
    r"""
    Linear form on a free module `M` over a commutative ring `R`.

    A *linear form* is a map `M\rightarrow R` that is linear.

    INPUT:
    
    - ``fmodule`` -- free module `M` of finite rank over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the linear form
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the linear 
      form; if none is provided, the LaTeX symbol is set to ``name``
      
    EXAMPLES:
    
    Linear form on a rank-3 free module::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: a = M.linear_form('A') ; a
        linear form A on the rank-3 free module M over the Integer Ring
        sage: a[:] = [2,-1,3]  # components w.r.t. the module's default basis (e)

    The members of a dual basis are linear forms::
    
        sage: e.dual_basis()[0]
        linear form e^0 on the rank-3 free module M over the Integer Ring
        sage: e.dual_basis()[1]
        linear form e^1 on the rank-3 free module M over the Integer Ring
        sage: e.dual_basis()[2]
        linear form e^2 on the rank-3 free module M over the Integer Ring

    Any linear form is expanded onto them::
    
        sage: a.view(basis=e)  # e being the default basis, it is equivalent to write a.view()
        A = 2 e^0 - e^1 + 3 e^2
        
    A linear form maps module elements to ring elements::

        sage: v = M([1,1,1])
        sage: a(v)
        4
        sage: a(v) in M.base_ring()
        True

    Test of linearity::
        
        sage: u = M([-5,-2,7])
        sage: a(3*u - 4*v) == 3*a(u) - 4*a(v)
        True

    A linear form is an element of the dual module::
        
        sage: a.parent()
        dual of the rank-3 free module M over the Integer Ring

    As such, it is a tensor of type (0,1)::

        sage: a.tensor_type
        (0, 1)

    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleAltForm.__init__(self, fmodule, 1, name=name, 
                                   latex_name=latex_name)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "linear form "
        if self.name is not None:
            description += self.name + " " 
        description += "on the " + str(self.fmodule)
        return description

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
              
        This method, which is already implemented in 
        :meth:`FreeModuleAltForm._new_comp`, is redefined here for efficiency
        """
        fmodule = self.fmodule  # the base free module
        return Components(fmodule.ring, basis, 1, start_index=fmodule.sindex,
                          output_formatter=fmodule.output_formatter)

    def __call__(self, vector):
        r"""
        The linear form acting on an element of the module.
        
        INPUT:
        
        - ``vector`` -- an element of the module (instance of 
          :class:`FiniteFreeModuleElement`)
        
        OUTPUT:
        
        - ring element `\langle \omega, v \rangle`
          
        """
        if not isinstance(vector, FiniteFreeModuleElement):
            raise TypeError("The argument must be a free module element.")
        basis = self.common_basis(vector)
        if basis is None:
            raise ValueError("No common basis for the components.")
        omega = self.components[basis]
        vv = vector.components[basis]
        resu = 0
        for i in self.fmodule.irange():
            resu += omega[[i]]*vv[[i]]
        # Name and LaTeX symbol of the output:
        if hasattr(resu, 'name'): 
            if self.name is not None and vector.name is not None:
                resu.name = self.name + "(" + vector.name + ")"
        if hasattr(resu, 'latex_name'): 
            if self.latex_name is not None and vector.latex_name is not None:
                resu.latex_name = self.latex_name + r"\left(" + \
                                  vector.latex_name + r"\right)"
        return resu

















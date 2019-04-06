r"""
Mixed Differential Forms

Let `M` and `N` be differentiable manifolds and `\Phi : M \longrightarrow N` a
differentiable map. A *mixed differential form along* `\Phi` is an element of
the graded algebra represented by
:class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`.
Its homogeneous components consist of differential forms along `\Phi`. Mixed
forms are useful to represent characteristic classes and perform computations
of such.

AUTHORS:

- Michael Jung (2019) : initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.element import AlgebraElement
from sage.rings.integer import Integer

class MixedForm(AlgebraElement):
    r"""
    An instance of this class is a mixed form along some differentiable map
    `\varphi: M \to N` between two differentiable manifolds `M` and `N`. More
    precisely, a mixed form `a` along `\varphi: M \to N` can be considered as a
    differentiable map
    
    .. MATH::

        a: M \longrightarrow \bigoplus^n_{k=0} T^{(0,k)}N,
    
    where `T^{(0,k)}` denotes the tensor bundle of type `(0,k)`, `\bigoplus`
    the Whitney sum and `n` the dimension of `N`, such that
    
    .. MATH:

        \forall x\in M, \quad a(x) \in \bigoplus^n_{k=0} \Lambda^k\left( T_{\varphi(x)}^* N \right),
    
    where `\Lambda^k(T^*_{\varphi(x)} N)` is the `k`-th exterior power of the
    dual of the tangent space `T_{\varphi(x)} N`.

    The standard case of a mixed form *on* `M` corresponds to `M=N` with
    `\varphi = \mathrm{Id}_M`.

    INPUT:

    - ``parent`` -- graded algebra of mixed forms represented by
      :class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`
      where the mixed form ``self`` shall belong to
    - ``comp```-- (default: ``None``) homogeneous components of the mixed form;
      if none is provided, the components are set to innocent unnamed
      differential forms
    - ``name`` -- (default: ``None``) name given to the mixed form
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      mixed form; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Initialize a mixed form on a 2-dimensional parallelizable differentiable
    manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: F = M.mixed_form(name='F'); F
        Mixed differential form F on the 2-dimensional differentiable manifold M
        sage: F.parent()
        Graded algebra Omega^*(M) of mixed differential forms on the
         2-dimensional differentiable manifold M

    To define the homogenous components of a mixed form, it is convenient to
    define some differential forms first::

        sage: f = M.scalar_field(x, name='f'); f
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: omega = M.diff_form(1, name='omega', latex_name=r'\omega'); omega
        1-form omega on the 2-dimensional differentiable manifold M
        sage: omega[X.frame(),0] = y*x; omega.disp()
        omega = x*y dx
        sage: eta = M.diff_form(2, name='eta', latex_name=r'\eta'); eta
        2-form eta on the 2-dimensional differentiable manifold M
        sage: eta[X.frame(),0,1] = y^2*x; eta.disp()
        eta = x*y^2 dx/\dy

    The components of the mixed form ``F`` can be manipulated very easily::

        sage: F[:] = [f, omega, eta]; F.disp()
        F = f + omega + eta # display names
        sage: F.disp(X.frame()) # display in coordinates
        F = [x] + [x*y dx] + [x*y^2 dx/\dy]
        sage: F[0]
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: F[0] is f
        True
        sage: F[1]
        1-form omega on the 2-dimensional differentiable manifold M
        sage: F[1] is omega
        True
        sage: F[2]
        2-form eta on the 2-dimensional differentiable manifold M
        sage: F[2] is eta
        True

    Moreover, it is possible to compute the exterior derivative of a
    mixed form::

        sage: dF = F.exterior_derivative(); dF.disp()
        dF = zero + df + domega
        sage: dF.disp(X.frame())
        dF = [0] + [dx] + [-x dx/\dy]

    Mixed forms can also be added and multiplied (via wedge product)::

        sage: F+dF
        Mixed differential form F+dF on the 2-dimensional differentiable
         manifold M
        sage: (F+dF).disp(X.frame())
        F+dF = [x] + [(x*y + 1) dx] + [(x*y^2 - x) dx/\dy]
        sage: (F*dF).disp(X.frame())
        F/\dF = [0] + [x dx] + [-x^2 dx/\dy]
        sage: (3*F).disp(X.frame())
        [3*x] + [3*x*y dx] + [3*x*y^2 dx/\dy]

    And thanks to the coercion system, these computations can be done with
    simple differential forms as well::

        sage: (omega*F).disp(X.frame())
        omega/\F = [0] + [x^2*y dx] + [0]
        sage: (omega+F).disp(X.frame())
        omega+F = [x] + [2*x*y dx] + [x*y^2 dx/\dy]

    Initialize a mixed form on a 2-dimensional non-parallelizable differentiable
    manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
        ....:                   intersection_name='W', restrictions1= x>0,
        ....:                   restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: e_xy = c_xy.frame(); e_uv = c_uv.frame() # define frames
        sage: omega = M.diff_form(1, name='omega', latex_name=r'\omega')
        sage: omega[c_xy.frame(),0] = y*x; omega.disp(e_xy)
        omega = x*y dx
        sage: eta = M.diff_form(2, name='eta', latex_name=r'\eta')
        sage: eta[c_uv.frame(),0,1] = u*v^2; eta.disp(e_uv)
        eta = u*v^2 du/\dv
        sage: F = M.mixed_form(name='F'); F
        Mixed differential form F on the 2-dimensional differentiable manifold M
        sage: F[:] = [x, omega, eta]
        sage: F.add_comp_by_continuation(e_uv, V.intersection(U), c_uv)
        sage: F.disp(e_uv)
        F = [1/2*u + 1/2*v] + [(1/8*u^2 - 1/8*v^2) du + (1/8*u^2 - 1/8*v^2) dv]
         + [u*v^2 du/\dv]
        sage: F.add_comp_by_continuation(e_xy, V.intersection(U), c_xy)
        sage: F.disp(e_xy)
        F = [x] + [x*y dx] + [(-2*x^3 + 2*x^2*y + 2*x*y^2 - 2*y^3) dx/\dy]


    """
    def __init__(self, parent, comp=None, name=None, latex_name=None):
        r"""
        Construct a mixed form.

        TESTS:

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
        ....:                   intersection_name='W', restrictions1= x>0,
        ....:                   restrictions2= u+v>0)
        sage: inv = transf.inverse()
        sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
        sage: omega = M.diff_form(1, name='omega')
        sage: omega[c_xy.frame(),0] = y*x
        sage: eta = M.diff_form(2, name='eta')
        sage: eta[c_uv.frame(),0,1] = u*v^2
        sage: A = M.mixed_form_algebra()
        sage: f = A([x, omega, eta], name='F')
        sage: f.add_comp_by_continuation(e_uv, V.intersection(U), c_uv)
        sage: TestSuite(f).run(skip='_test_pickling')
        """
        if parent is None:
            raise ValueError("a parent must be provided")

        # Add this element (instance) to the set of parent algebra:
        AlgebraElement.__init__(self, parent)

        vmodule = parent._vmodule

        # Define attributes:
        self._vmodule = vmodule
        self._dest_map = vmodule._dest_map
        self._domain = vmodule._domain
        self._ambient_domain = vmodule._ambient_domain
        self._max_deg = vmodule._ambient_domain.dim()

        # Set components:
        if comp is None:
            self._comp = [self._domain.diff_form(j)
                          for j in range(0, self._max_deg + 1)]
        else:
            self._comp = comp

        # Set names:
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

    def _repr_(self):
        r"""
        """
        description = "Mixed differential form "
        if self._name is not None:
            description += self._name + " "
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} with values on the {} ".format(
                self._domain, self._ambient_domain)
            if self._dest_map._name is None:
                dm_name = "unnamed map"
            else:
                dm_name = self._dest_map._name
            description += "via the map " + dm_name
        return description

    def _latex_(self):
        if self._name is None:
            return r'\mbox{' + repr(self) + r'}'
        else:
            return self._latex_name

    #####
    # display
    #
    # Return mixed differential form with all its components nicely.
    #####

    def display(self, basis=None, chart=None, from_chart=None):

        from sage.misc.latex import latex, LatexExpr
        from sage.tensor.modules.format_utilities import FormattedExpansion

        #####
        # display._display_form
        #
        # Return coord expression of a single form.
        #####

        def _display_form(self, rst, basis, chart):
            from sage.tensor.modules.format_utilities import is_atomic

            cobasis = basis.dual_basis()
            comp = rst.comp(basis)
            terms_txt = []
            terms_latex = []
            for ind in comp.non_redundant_index_generator():
                ind_arg = ind + (chart,)
                coef = comp[ind_arg]
                # Check whether the coefficient is zero, preferably via
                # the fast method is_trivial_zero():
                if hasattr(coef, 'is_trivial_zero'):
                    zero_coef = coef.is_trivial_zero()
                else:
                    zero_coef = coef == 0
                if not zero_coef:
                    bases_txt = []
                    bases_latex = []
                    for k in range(rst._tensor_rank):
                        bases_txt.append(cobasis[ind[k]]._name)
                        bases_latex.append(latex(cobasis[ind[k]]))
                    basis_term_txt = "/\\".join(bases_txt)
                    basis_term_latex = r"\wedge ".join(bases_latex)
                    coef_txt = repr(coef)
                    if coef_txt == "1":
                        terms_txt.append(basis_term_txt)
                        terms_latex.append(basis_term_latex)
                    elif coef_txt == "-1":
                        terms_txt.append("-" + basis_term_txt)
                        terms_latex.append("-" + basis_term_latex)
                    else:
                        coef_latex = latex(coef)
                        if is_atomic(coef_txt):
                            terms_txt.append(coef_txt + " " + basis_term_txt)
                        else:
                            terms_txt.append(
                                "(" + coef_txt + ") " + basis_term_txt)
                        if is_atomic(coef_latex):
                            terms_latex.append(coef_latex + basis_term_latex)
                        else:
                            terms_latex.append(r"\left(" + coef_latex +
                                               r"\right)" + basis_term_latex)
            if not terms_txt:
                resu_txt = "0"
            else:
                resu_txt = terms_txt[0]
                for term in terms_txt[1:]:
                    if term[0] == "-":
                        resu_txt += " - " + term[1:]
                    else:
                        resu_txt += " + " + term
            if not terms_latex:
                resu_latex = r"0"
            else:
                resu_latex = terms_latex[0]
                for term in terms_latex[1:]:
                    if term[0] == "-":
                        resu_latex += term
                    else:
                        resu_latex += "+" + term
            return FormattedExpansion(resu_txt, resu_latex)

        plain_txt = ""
        latex_txt = r""

        if self._name is not None:
            plain_txt += self._name + " = "
        if self._latex_name is not None:
            latex_txt += latex(self) + LatexExpr(r"=")

        if basis is None:
            if self._comp[0]._name is None:
                plain_txt += "(unnamed scalar field) "
            else:
                plain_txt += self._comp[0]._name
            if self._comp[0]._latex_name is None:
                latex_txt += LatexExpr(r"\mbox{(unnamed scalar field)}")
            else:
                latex_txt += latex(self._comp[0])

            for j in range(1, self._max_deg + 1):
                if self._comp[j]._name is None:
                    plain_txt += " + (unnamed " + str(j) + "-form)"
                else:
                    plain_txt += " + " + self._comp[j]._name
                if self._comp[j]._latex_name is None:
                    latex_txt += LatexExpr(
                        r"+\mbox{(unnamed " + str(j) + r"-form)}")
                else:
                    latex_txt += LatexExpr(r"+") + latex(self._comp[j])
            return FormattedExpansion(plain_txt, latex_txt)

        if chart is None:
            chart = basis._chart

        plain_txt += "[" + repr(self._comp[0].expr(chart, from_chart)) + "]"
        latex_txt += LatexExpr(r"\left[") + latex(
            self._comp[0].expr(chart, from_chart)) \
                     + LatexExpr(r"\right]_0")

        for j in range(1, self._max_deg + 1):
            rst = self._comp[j].restrict(basis._domain,
                                         dest_map=basis._dest_map)
            format_txt = _display_form(self, rst, basis, chart)
            plain_txt += " + [" + repr(format_txt) + "]"
            latex_txt += LatexExpr(r"+ \left[") + latex(format_txt) \
                         + LatexExpr(r"\right]_{}".format(j))
        return FormattedExpansion(plain_txt, latex_txt)

    disp = display

    def _new_instance(self):
        r"""
        Create an instance of the same class, on the same domain.

        TESTS::



        """
        return type(self)(self.parent())

    def set_name(self, name=None, latex_name=None):
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name

    #####
    # __bool__
    #
    # Return whether ``self`` is nonzero
    #####

    def __bool__(self):
        return any(bool(form) for form in self._comp)

    __nonzero__ = __bool__  # For Python2 compatibility

    #####
    # _richcmp_
    #
    # Return comparison
    #####

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import richcmp

        # Compare all elements separately:
        for j in range(0, self._max_deg + 1):
            if not richcmp(self._comp[j], other._comp[j], op):
                return False
        return True

    #####
    # _add_
    #
    # Ring addition
    #####

    def _add_(self, other):
        resu_comp = [self._comp[j] + other._comp[j]
                     for j in range(0, self._max_deg + 1)]
        resu = self.__class__(self.parent(), comp=resu_comp)
        # Compose name:
        if self._name is not None and other._name is not None:
            resu._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            resu._latex_name = self._latex_name + '+' + other._latex_name

        return resu

    #####
    # _sub_
    #
    # Ring substraction
    #####

    def _sub_(self, other):
        resu_comp = [self._comp[j] - other._comp[j]
                     for j in range(0, self._max_deg + 1)]
        resu = self.__class__(self.parent(), comp=resu_comp)
        # Compose name:
        from sage.tensor.modules.format_utilities import is_atomic

        if self._name is not None and other._name is not None:
            sname = self._name
            oname = other._name
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            resu._name = sname + '-' + oname
        if self._latex_name is not None and other._latex_name is not None:
            slname = self._latex_name
            olname = other._latex_name
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            resu._latex_name = slname + '-' + olname

        return resu

        #####

    # _mul_
    #
    # Ring multiplication
    #####

    def _mul_(self, other):
        return self.wedge(other)

    #####
    # wedge
    #
    # Return the wedge product of two mixed forms.
    #####

    def wedge(self, other):
        resu_comp = [None] * (self._max_deg + 1)
        resu_comp[0] = self._comp[0] * other._comp[0]
        for j in range(1, self._max_deg + 1):
            resu_comp[j] = self._comp[j] * other._comp[0]
            resu_comp[j] = resu_comp[j] + self._comp[0] * other._comp[j]
            for k in range(1, j):
                resu_comp[j] = resu_comp[j] + \
                               self._comp[k].wedge(other._comp[j - k])
        resu = self.__class__(self.parent(), comp=resu_comp)
        # Compose name:
        from sage.tensor.modules.format_utilities import format_mul_txt, \
            format_mul_latex

        resu._name = format_mul_txt(self._name, '/\\', other._name)
        resu._latex_name = format_mul_latex(self._latex_name, r'\wedge ',
                                            other._latex_name)
        return resu

    #####
    # exterior_derivative
    #
    # Return the exterior derivative of ``self``
    #####

    @cached_method
    def exterior_derivative(self):
        resu_comp = list()
        resu_comp.append(self._domain.scalar_field_algebra().zero())
        resu_comp.append(self._comp[0].differential())
        resu_comp.extend([self._comp[j].exterior_derivative()
                          for j in range(1, self._max_deg)])
        resu = self.__class__(self.parent(), comp=resu_comp)
        # Compose name:
        from sage.tensor.modules.format_utilities import (format_unop_txt,
                                                          format_unop_latex)
        resu._name = format_unop_txt('d', self._name)
        resu._latex_name = format_unop_latex(r'\mathrm{d}', self._latex_name)

        return resu

    #####
    # copy
    #
    # Return an exact copy of ``self``.
    #####

    def copy(self):
        return self.__class__(self.parent(), comp=self._comp, name=self._name,
                              latex_name=self._latex_name)

    #####
    # __setitem__
    #
    # Set the differential forms via components.
    #####

    def __setitem__(self, index, values):
        if isinstance(index, (int, Integer)):
            start = index
            stop = index + 1
            step = 1
        elif isinstance(index, slice):
            start, stop, step = index.indices(len(self._comp))
        else:
            raise TypeError("index must be int, Integer or slice")
        if isinstance(values, list):
            form_list = values
        else:
            form_list = [values]

        if len(form_list) != len(range(start, stop, step)):
            raise IndexError("either input or index out of range")

        for deg, j in zip(range(start, stop, step), range(0, len(form_list))):
            if form_list[j] in self._domain.diff_form_module(deg,
                                                             self._dest_map):
                self._comp[deg] = \
                    self._domain.diff_form_module(deg,
                                            self._dest_map).coerce(form_list[j])
            else:
                # A very datailed error message:
                error_msg = "input must be a " \
                            "differential form of degree {}".format(deg)
                if self._dest_map is self._domain.identity_map():
                    error_msg += " on the {}".format(self._domain)
                else:
                    error_msg += " along {} mapped " \
                                 "into {}".format(self._domain,
                                                   self._ambient_domain)
                    if self._dest_map._name is None:
                        dm_name = "unnamed map"
                    else:
                        dm_name = self._dest_map._name
                    error_msg += " via " + dm_name
                raise TypeError(error_msg)

    #####
    # __setitem__
    #
    # Get the differential form of degree ``deg``.
    #####

    def __getitem__(self, deg):
        return self._comp[deg]

    #####
    # set_restriction
    #
    # Define a restriction of the components of ``self`` to some subdomain.
    #####

    def set_restriction(self, rst):
        if not isinstance(rst, MixedForm):
            raise TypeError("the argument must be a mixed form")
        for j in range(0, self._max_deg + 1):
            self._comp[j].restrict(rst._comp[j])

    #####
    # restrict
    #
    # Return the restriction of ``self`` to some subdomain.
    #####

    def restrict(self, subdomain, dest_map=None):
        resu_comp = list()
        resu_comp.append(self._comp[0].restrict(subdomain))
        resu_comp.extend([self._comp[j].restrict(subdomain, dest_map)
                          for j in range(1, self._max_deg + 1)])
        return self.__class__(subdomain.mixed_form_algebra(dest_map=dest_map),
                              resu_comp, self._name, self._latex_name)

    #####
    # add_comp_by_continuation
    #
    # Set components with respect to a vector frame by continuation of the
    # coordinate expression of the components in a subframe.
    #####

    def add_comp_by_continuation(self, frame, subdomain, chart=None):
        if chart is None:
            chart = frame._chart
        self._comp[0].add_expr_by_continuation(chart, subdomain)
        for j in range(1, self._max_deg + 1):
            self._comp[j].add_comp_by_continuation(frame, subdomain, chart)

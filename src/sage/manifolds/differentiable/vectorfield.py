r"""
Vector Fields

Given two differentiable manifolds `U` and `M` over the same topological field
`K` and a differentiable map

.. MATH::

    \Phi:\ U \longrightarrow  M,

we define a *vector field along* `U` *with values on* `M` to be a
differentiable map

.. MATH::

    v:\ U  \longrightarrow TM

(`TM` being the tangent bundle of `M`) such that

.. MATH::

    \forall p \in U,\ v(p) \in T_{\Phi(p)}M.

The standard case of vector fields *on* a differentiable manifold corresponds
to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
being an immersion and `\Phi` being a curve in `M` (`U` is then an open
interval of `\RR`).

Vector fields are implemented via two classes: :class:`VectorFieldParal` and
:class:`VectorField`, depending respectively whether the manifold `M`
is parallelizable or not, i.e. whether the bundle `TM` is trivial or not.


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [KN1963]_
- [Lee2013]_
- [ONe1983]_
- [BG1988]_

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_tensor import FiniteRankFreeModuleElement
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal

class VectorField(TensorField):
    r"""
    Vector field along a differentiable manifold.

    An instance of this class is a vector field along a differentiable
    manifold `U` with values on a differentiable manifold `M`, via a
    differentiable map `U \rightarrow M`. More precisely, given a
    differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *vector field along* `U` *with values on* `M` is a differentiable map

    .. MATH::

        v:\ U  \longrightarrow TM

    (`TM` being the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M.

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are
    `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then an
    open interval of `\RR`).

    .. NOTE::

        If `M` is parallelizable, then
        :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal`
        *must* be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` with values on `M\supset\Phi(U)`
    - ``name`` -- (default: ``None``) name given to the vector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the vector
      field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a non-parallelizable 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_tu.<t,u> = V.chart()
        sage: transf = c_xy.transition_map(c_tu, (x+y, x-y), intersection_name='W',
        ....:                              restrictions1= x>0, restrictions2= t+u>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_tu.frame()
        sage: c_tuW = c_tu.restrict(W) ; eVW = c_tuW.frame()
        sage: v = M.vector_field('v') ; v
        Vector field v on the 2-dimensional differentiable manifold M
        sage: v.parent()
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M

    The vector field is first defined on the domain `U` by means of its
    components with respect to the frame ``eU``::

        sage: v[eU,:] = [-y, 1+x]

    The components with respect to the frame ``eV`` are then deduced
    by continuation of the components with respect to the frame ``eVW``
    on the domain `W = U \cap V`, expressed in terms on the coordinates
    covering `V`::

        sage: v[eV,0] = v[eVW,0,c_tuW].expr()
        sage: v[eV,1] = v[eVW,1,c_tuW].expr()

    At this stage, the vector field is fully defined on the whole manifold::

        sage: v.display(eU)
        v = -y d/dx + (x + 1) d/dy
        sage: v.display(eV)
        v = (u + 1) d/dt + (-t - 1) d/du

    The vector field acting on scalar fields::

        sage: f = M.scalar_field({c_xy: (x+y)^2, c_tu: t^2}, name='f')
        sage: s = v(f) ; s
        Scalar field v(f) on the 2-dimensional differentiable manifold M
        sage: s.display()
        v(f): M --> R
        on U: (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
        on V: (t, u) |--> 2*t*u + 2*t

    Some checks::

        sage: v(f) == f.differential()(v)
        True
        sage: v(f) == f.lie_der(v)
        True

    The result is defined on the intersection of the vector field's
    domain and the scalar field's one::

        sage: s = v(f.restrict(U)) ; s
        Scalar field v(f) on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s == v(f).restrict(U)
        True
        sage: s = v(f.restrict(W)) ; s
        Scalar field v(f) on the Open subset W of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): W --> R
           (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
           (t, u) |--> 2*t*u + 2*t
        sage: s = v.restrict(U)(f) ; s
        Scalar field v(f) on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): U --> R
           (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
        on W: (t, u) |--> 2*t*u + 2*t
        sage: s = v.restrict(U)(f.restrict(V)) ; s
        Scalar field v(f) on the Open subset W of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): W --> R
           (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
           (t, u) |--> 2*t*u + 2*t

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a vector field with values on a non-parallelizable manifold.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``VectorField``, to fit with the category framework::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: XM = M.vector_field_module()
            sage: a = XM.element_class(XM, name='a'); a
            Vector field a on the 2-dimensional differentiable manifold M
            sage: a[c_xy.frame(),:] = [x, y]
            sage: a[c_uv.frame(),:] = [-u, -v]
            sage: TestSuite(a).run(skip='_test_pickling')

        Construction with ``DifferentiableManifold.vector_field``::

            sage: a1 = M.vector_field(name='a'); a1
            Vector field a on the 2-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True

        .. TODO::

            Fix ``_test_pickling`` (in the superclass :class:`TensorField`).

        """
        TensorField.__init__(self, vector_field_module, (1,0), name=name,
                             latex_name=latex_name)
        # Initialization of derived quantities:
        TensorField._init_derived(self)
        # Initialization of list of quantities depending on self:
        self._init_dependencies()

    def _repr_(self) :
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._repr_()
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: repr(v)  # indirect doctest
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: v  # indirect doctest
            Vector field v on the 2-dimensional differentiable manifold M

        """
        description = "Vector field "
        if self._name is not None:
            description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same module.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: u = v._new_instance(); u
            Vector field on the 2-dimensional differentiable manifold M
            sage: u.parent() is v.parent()
            True

        """
        return type(self)(self._vmodule)

    def _init_dependencies(self):
        r"""
        Initialize list of quantities that depend on ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._init_dependencies()

        """
        self._lie_der_along_self = {}

    def _del_dependencies(self):
        r"""
        Clear list of quantities that depend on ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._del_dependencies()

        """
        if self._lie_der_along_self != {}:
            for idtens, tens in self._lie_der_along_self.items():
                del tens._lie_derivatives[id(self)]
            self._lie_der_along_self.clear()

    def __call__(self, scalar):
        r"""
        Action on a scalar field (or on a 1-form).

        INPUT:

        - ``scalar`` -- scalar field `f`

        OUTPUT:

        - scalar field representing the derivative of `f` along the vector
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: a = M.vector_field(name='a')
            sage: a[c_xy.frame(),:] = [x, y]
            sage: a[c_uv.frame(),:] = [-u, -v]
            sage: f = M.scalar_field({c_xy: atan(x^2+y^2), c_uv: pi/2-atan(u^2+v^2)},
            ....:                    name='f')
            sage: s = a.__call__(f); s
            Scalar field a(f) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(f): M --> R
            on U: (x, y) |--> 2*(x^2 + y^2)/(x^4 + 2*x^2*y^2 + y^4 + 1)
            on V: (u, v) |--> 2*(u^2 + v^2)/(u^4 + 2*u^2*v^2 + v^4 + 1)
            sage: s == f.differential()(a)
            True

        """
        if scalar._tensor_type == (0,1):
            # This is actually the action of the vector field on a 1-form,
            # as a tensor field of type (1,0):
            return scalar(self)
        if scalar._tensor_type != (0,0):
            raise TypeError("the argument must be a scalar field")
        #!# Could it be simply
        # return scalar.differential()(self)
        # ?
        dom_resu = self._domain.intersection(scalar._domain)
        self_r = self.restrict(dom_resu)
        scalar_r = scalar.restrict(dom_resu)
        if scalar_r._is_zero:
            return dom_resu._zero_scalar_field
        if isinstance(self_r, VectorFieldParal):
            return self_r(scalar_r)
        # Creation of the result:
        if self._name is not None and scalar._name is not None:
            resu_name = "{}({})".format(self._name, scalar._name)
        else:
            resu_name = None
        if self._latex_name is not None and scalar._latex_name is not None:
            resu_latex = r"{}\left({}\right)".format(self._latex_name , scalar._latex_name)
        else:
            resu_latex = None
        resu = dom_resu.scalar_field(name=resu_name, latex_name=resu_latex)
        for dom, rst in self_r._restrictions.items():
            resu_rst = rst(scalar_r.restrict(dom))
            for chart, funct in resu_rst._express.items():
                resu._express[chart] = funct
        return resu

#******************************************************************************

class VectorFieldParal(FiniteRankFreeModuleElement, TensorFieldParal, VectorField):
    r"""
    Vector field along a differentiable manifold, with values on a
    parallelizable manifold.

    An instance of this class is a vector field along a differentiable
    manifold `U` with values on a parallelizable manifold `M`, via a
    differentiable map `\Phi: U \rightarrow M`. More precisely, given
    a differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *vector field along* `U` *with values on* `M` is a differentiable map

    .. MATH::

        v:\ U  \longrightarrow TM

    (`TM` being the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M.

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases
    are `\Phi` being an immersion and `\Phi` being a curve in `M` (`U`
    is then an open interval of `\RR`).

    .. NOTE::

        If `M` is not parallelizable, then
        :class:`~sage.manifolds.differentiable.vectorfield.VectorField`
        *must* be used instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` with values on `M\supset\Phi(U)`
    - ``name`` -- (default: ``None``) name given to the vector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the vector
      field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a parallelizable 3-dimensional manifold::

        sage: M = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart()
        sage: v = M.vector_field('V') ; v
        Vector field V on the 3-dimensional differentiable manifold M
        sage: latex(v)
        V

    Vector fields are considered as elements of a module over the ring
    (algebra) of scalar fields on `M`::

        sage: v.parent()
        Free module X(M) of vector fields on the 3-dimensional differentiable
         manifold M
        sage: v.parent().base_ring()
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M
        sage: v.parent() is M.vector_field_module()
        True

    A vector field is a tensor field of rank 1 and of type `(1,0)`::

        sage: v.tensor_rank()
        1
        sage: v.tensor_type()
        (1, 0)

    Components of a vector field with respect to a given frame::

        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v[0], v[1], v[2] = (1, 4, 9)  # components on M's default frame (e)
        sage: v.comp()
        1-index components w.r.t. Vector frame (M, (e_0,e_1,e_2))

    The totality of the components are accessed via the operator ``[:]``::

        sage: v[:] = (1, 4, 9)  # equivalent to v[0], v[1], v[2] = (1, 4, 9)
        sage: v[:]
        [1, 4, 9]

    The components are also read on the expansion on the frame ``e``,
    as provided by the method
    :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.display`::

        sage: v.display()  # displays the expansion in the default frame
        V = e_0 + 4 e_1 + 9 e_2

    A subset of the components can be accessed by using slice notation::

        sage: v[1:] = (-2, -3)
        sage: v[:]
        [1, -2, -3]
        sage: v[:2]
        [1, -2]

    Components in another frame::

        sage: f = M.vector_frame('f')
        sage: for i in range(3):
        ....:     v.set_comp(f)[i] = (i+1)**3
        ....:
        sage: v.comp(f)[2]
        27
        sage: v[f, 2]  # equivalent to above
        27
        sage: v.display(f)
        V = f_0 + 8 f_1 + 27 f_2

    The range of the indices depends on the convention set for the manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v = M.vector_field('V')
        sage: (v[1], v[2], v[3]) = (1, 4, 9)
        sage: v[0]
        Traceback (most recent call last):
        ...
        IndexError: index out of range: 0 not in [1, 3]

    A vector field acts on scalar fields (derivation along the vector field)::

        sage: M = Manifold(2, 'M')
        sage: c_cart.<x,y> = M.chart()
        sage: f = M.scalar_field(x*y^2, name='f')
        sage: v = M.vector_field('v')
        sage: v[:] = (-y, x)
        sage: v.display()
        v = -y d/dx + x d/dy
        sage: v(f)
        Scalar field v(f) on the 2-dimensional differentiable manifold M
        sage: v(f).expr()
        2*x^2*y - y^3
        sage: latex(v(f))
        v\left(f\right)

    Example of a vector field associated with a non-trivial map `\Phi`;
    a vector field along a curve in `M`::

        sage: R = Manifold(1, 'R')
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: Phi = R.diff_map(M, [cos(t), sin(t)], name='Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold R
         to the 2-dimensional differentiable manifold M
        sage: Phi.display()
        Phi: R --> M
           t |--> (x, y) = (cos(t), sin(t))
        sage: w = R.vector_field('w', dest_map=Phi) ; w
        Vector field w along the 1-dimensional differentiable manifold R with
         values on the 2-dimensional differentiable manifold M
        sage: w.parent()
        Free module X(R,Phi) of vector fields along the 1-dimensional
         differentiable manifold R mapped into the 2-dimensional differentiable
         manifold M
        sage: w[:] = (-sin(t), cos(t))
        sage: w.display()
        w = -sin(t) d/dx + cos(t) d/dy

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a vector field with values on a parallelizable manifold.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``VectorFieldParal``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()  # the parent
            sage: v = XM.element_class(XM, name='v'); v
            Vector field v on the 2-dimensional differentiable manifold M
            sage: v[:] = (-y, x)
            sage: v.display()
            v = -y d/dx + x d/dy
            sage: TestSuite(v).run()

        Construction via ``DifferentiableManifold.vector_field``::

            sage: u = M.vector_field(name='u'); u
            Vector field u on the 2-dimensional differentiable manifold M
            sage: type(u) == type(v)
            True
            sage: u.parent() is v.parent()
            True
            sage: u[:] = (1+x, 1-y)
            sage: TestSuite(u).run()

        """
        FiniteRankFreeModuleElement.__init__(self, vector_field_module,
                                             name=name, latex_name=latex_name)
        # TensorFieldParal attributes:
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # VectorField attributes:
        self._vmodule = vector_field_module
        # Initialization of derived quantities:
        TensorFieldParal._init_derived(self)
        VectorField._init_derived(self)
        # Initialization of list of quantities depending on self:
        self._init_dependencies()

    def _repr_(self) :
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: v._repr_()
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: repr(v)  # indirect doctest
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: v  # indirect doctest
            Vector field v on the 2-dimensional differentiable manifold M

        """
        return VectorField._repr_(self)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same module.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: u = v._new_instance(); u
            Vector field on the 2-dimensional differentiable manifold M
            sage: u.parent() is v.parent()
            True

        """
        return type(self)(self._fmodule)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- (default: ``True``) determines whether
          the restrictions of ``self`` to subdomains are deleted

        TEST::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: v._del_derived()

        """
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)
        VectorField._del_derived(self)
        self._del_dependencies()

    def __call__(self, scalar):
        r"""
        Action on a scalar field.

        INPUT:

        - ``scalar`` -- scalar field `f`

        OUTPUT:

        - scalar field representing the derivative of `f` along the vector
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`

        EXAMPLES:

        Action of a vector field on a scalar field on a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_cart.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: v = M.vector_field()
            sage: v[:] = (-y, x)
            sage: v(f)
            Scalar field on the 2-dimensional differentiable manifold M
            sage: v(f).display()
            M --> R
            (x, y) |--> 2*x^2*y - y^3

        """
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        if scalar._tensor_type == (0,1):
            # This is actually the action of the vector field on a 1-form,
            # as a tensor field of type (1,0):
            return scalar(self)
        if scalar._tensor_type != (0,0):
            raise TypeError("the argument must be a scalar field")
        #!# Could it be simply
        # return scalar.differential()(self)
        # ?
        dom_resu = self._domain.intersection(scalar._domain)
        self_r = self.restrict(dom_resu)
        scalar_r = scalar.restrict(dom_resu)
        if scalar_r._is_zero:
            return dom_resu._zero_scalar_field
        # Creation of the result:
        if self._name is not None and scalar._name is not None:
            resu_name = self._name + "(" + scalar._name + ")"
        else:
            resu_name = None
        if self._latex_name is not None and scalar._latex_name is not None:
            resu_latex = (self._latex_name + r"\left(" + scalar._latex_name
                          + r"\right)")
        else:
            resu_latex = None
        resu = dom_resu.scalar_field(name=resu_name, latex_name=resu_latex)
        # Search for common charts for the computation:
        common_charts = set()
        for chart in scalar_r._express:
            try:
                self_r.comp(chart._frame)
                common_charts.add(chart)
            except ValueError:
                pass
        for frame in self_r._components:
            if isinstance(frame, CoordFrame):
                chart = frame._chart
                try:
                    scalar_r.coord_function(chart)
                    common_charts.add(chart)
                except ValueError:
                    pass
        if not common_charts:
            raise ValueError("no common chart found")
        # The computation:
        manif = scalar._manifold
        for chart in common_charts:
            v = self_r.comp(chart._frame)
            f = scalar_r.coord_function(chart)
            res = 0
            for i in manif.irange():
                res += v[i, chart] * f.diff(i)
            resu._express[chart] = res
        return resu


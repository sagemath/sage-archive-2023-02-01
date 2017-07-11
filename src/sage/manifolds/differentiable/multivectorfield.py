r"""
Multivector Fields

Let `U` and `M` be two differentiable manifolds.
Given a positive integer `p` and a differentiable map `\Phi: U \rightarrow M`,
a *multivector field of degree* `p`, or `p`-*vector field*,
*along* `U` *with values on* `M` is a field along `U` of alternating
contravariant tensors of rank `p` in the tangent spaces to `M`.
The standard case of a multivector field *on* a differentiable manifold
corresponds to `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are
`\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then an open
interval of `\RR`).

Two classes implement multivector fields, depending whether the manifold
`M` is parallelizable:

* :class:`MultivectorFieldParal` when `M` is parallelizable
* :class:`MultivectorField` when `M` is not assumed parallelizable.

AUTHORS:

- Eric Gourgoulhon (2017): initial version

REFERENCES:

- C.-M. Marle (1997)

"""

#******************************************************************************
#       Copyright (C) 2017 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.tensor.modules.alternating_contr_tensor import AlternatingContrTensor
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal

class MultivectorField(TensorField):
    r"""
    Multivector field with values on a generic (i.e. a priori not
    parallelizable) differentiable manifold.

    Given a differentiable manifold `U`, a differentiable map
    `\Phi: U \rightarrow M` to a differentiable manifold `M` and a positive
    integer `p`, a *multivector field of degree* `p` (or `p`-*vector field*)
    *along* `U` *with values on* `M\supset\Phi(U)` is a differentiable map

    .. MATH::

        a:\ U  \longrightarrow T^{(p,0)}M

    (`T^{(p,0)}M` being the tensor bundle of type `(p,0)` over `M`) such that

    .. MATH::

        \forall x \in U,\quad a(x) \in \Lambda^p(T_{\Phi(x)} M) ,

    where `T_{\Phi(x)} M` is the vector space tangent to `M` at `\Phi(x)` and
    `\Lambda^p` stands for the exterior power of degree `p` (cf.
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule`).
    In other words, `a(x)` is an alternating contravariant tensor of degree `p`
    of the tangent vector space `T_{\Phi(x)} M`.

    The standard case of a multivector field *on* a manifold `M` corresponds to
    `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval of
    `\RR`).

    .. NOTE::

        If `M` is parallelizable, the class :class:`MultivectorFieldParal`
        must be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector fields
      along `U` with values on `M` via the map `\Phi`
    - ``degree`` -- the degree of the multivector field (i.e. its tensor rank)
    - ``name`` -- (default: ``None``) name given to the multivector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      multivector field; if none is provided, the LaTeX symbol is set to
      ``name``

    EXAMPLES:

    Multivector field of degree 2 on a non-parallelizable 2-dimensional
    manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
        ....:                         intersection_name='W',
        ....:                         restrictions1= x>0, restrictions2= u+v>0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_uv.frame()
        sage: a = M.multivector_field(2, name='a') ; a
        2-vector field a on the 2-dimensional differentiable manifold M
        sage: a.parent()
        Module A^2(M) of 2-vector fields on the 2-dimensional differentiable
         manifold M
        sage: a.degree()
        2

    Setting the components of ``a``::

        sage: a[eU,0,1] = x*y^2 + 2*x
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: a.display(eU)
        a = (x*y^2 + 2*x) d/dx/\d/dy
        sage: a.display(eV)
        a = (-1/4*u^3 + 1/4*u*v^2 - 1/4*v^3 + 1/4*(u^2 - 8)*v - 2*u) d/du/\d/dv

    The exterior product of two vector fields is a 2-vector field::

        sage: a = M.vector_field(name='a')
        sage: a[eU,:] = [-y, x]
        sage: a.add_comp_by_continuation(eV, W, c_uv)
        sage: b = M.vector_field(name='b')
        sage: b[eU,:] = [1+x*y, x^2]
        sage: b.add_comp_by_continuation(eV, W, c_uv)
        sage: s = a.wedge(b) ; s
        2-vector field a/\b on the 2-dimensional differentiable manifold M
        sage: s.display(eU)
        a/\b = (-2*x^2*y - x) d/dx/\d/dy
        sage: s.display(eV)
        a/\b = (1/2*u^3 - 1/2*u*v^2 - 1/2*v^3 + 1/2*(u^2 + 2)*v + u) d/du/\d/dv

    Multiplying a 2-vector field by a scalar field results in another
    2-vector field::

        sage: f = M.scalar_field({c_xy: (x+y)^2, c_uv: u^2}, name='f')
        sage: s = f*s ; s
        2-vector field on the 2-dimensional differentiable manifold M
        sage: s.display(eU)
        (-2*x^2*y^3 - x^3 - (4*x^3 + x)*y^2 - 2*(x^4 + x^2)*y) d/dx/\d/dy
        sage: s.display(eV)
        (1/2*u^5 - 1/2*u^3*v^2 - 1/2*u^2*v^3 + u^3 + 1/2*(u^4 + 2*u^2)*v)
          d/du/\d/dv

    """
    def __init__(self, vector_field_module, degree, name=None, latex_name=None):
        r"""
        Construct a multivector field.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``MultivectorField`, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: A = M.multivector_module(2)
            sage: XM = M.vector_field_module()
            sage: a = A.element_class(XM, 2, name='a'); a
            2-vector field a on the 2-dimensional differentiable manifold M
            sage: a[e_xy,0,1] = x+y
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: TestSuite(a).run(skip='_test_pickling')

        Construction with ``DifferentiableManifold.multivector_field``::

            sage: a1 = M.multivector_field(2, name='a'); a1
            2-vector field a on the 2-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True

        .. TODO::

            Fix ``_test_pickling`` (in the superclass :class:`TensorField`).

        """
        TensorField.__init__(self, vector_field_module, (degree, 0), name=name,
                             latex_name=latex_name, antisym=range(degree),
                             parent=vector_field_module.exterior_power(degree))
        self._init_derived() # initialization of derived quantities

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: a = M.multivector_field(2, name='a')
            sage: a._repr_()
            '2-vector field a on the 3-dimensional differentiable manifold M'
            sage: repr(a)  # indirect doctest
            '2-vector field a on the 3-dimensional differentiable manifold M'
            sage: a  # indirect doctest
            2-vector field a on the 3-dimensional differentiable manifold M
            sage: b = M.multivector_field(2)
            sage: b._repr_()
            '2-vector field on the 3-dimensional differentiable manifold M'

        """
        description = "{}-vector field ".format(self._tensor_rank)
        if self._name is not None:
            description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class, of the same degree and on the
        same domain.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: a = M.multivector_field(2, name='a')
            sage: a1 = a._new_instance(); a1
            2-vector field on the 3-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True

        """
        return type(self)(self._vmodule, self._tensor_rank)

    def wedge(self, other):
        r"""
        Exterior product with another multivector field.

        INPUT:

        - ``other`` -- another multivector field (on the same manifold)

        OUTPUT:

        - instance of :class:`MultivectorField` representing the exterior
          product ``self/\other``

        EXAMPLES:

        Exterior product of two vector fields on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1) # the sphere S^2
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() # stereographic coord. North
            sage: c_uv.<u,v> = V.chart() # stereographic coord. South
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: a = M.vector_field(name='a')
            sage: a[e_xy,:] = y, x
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: b = M.vector_field(name='b')
            sage: b[e_xy,:] = x^2 + y^2, y
            sage: b.add_comp_by_continuation(e_uv, W, c_uv)
            sage: c = a.wedge(b); c
            2-vector field a/\b on the 2-dimensional differentiable
             manifold S^2
            sage: c.display(e_xy)
            a/\b = (-x^3 - (x - 1)*y^2) d/dx/\d/dy
            sage: c.display(e_uv)
            a/\b = (-v^2 + u) d/du/\d/dv

        """
        from sage.tensor.modules.format_utilities import is_atomic
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior " +
                                 "product")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior " +
                                 "product")
        dom_resu = self._domain.intersection(other._domain)
        ambient_dom_resu = self._ambient_domain.intersection(other._ambient_domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the AlternatingContrTensor version:
            return AlternatingContrTensor.wedge(self_r, other_r)
        # otherwise, the result is created here:
        if self._name is not None and other._name is not None:
            sname = self._name
            oname = other._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            resu_name = sname + '/\\' + oname
        if self._latex_name is not None and other._latex_name is not None:
            slname = self._latex_name
            olname = other._latex_name
            if not is_atomic(slname):
                slname = '(' + slname + ')'
            if not is_atomic(olname):
                olname = '(' + olname + ')'
            resu_latex_name = slname + r'\wedge ' + olname
        dest_map = self._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu,
                                          subcodomain=ambient_dom_resu)
        vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)
        resu_degree = self._tensor_rank + other._tensor_rank
        resu = vmodule.alternating_contravariant_tensor(resu_degree,
                                    name=resu_name, latex_name=resu_latex_name)
        for dom in self_r._restrictions:
            if dom in other_r._restrictions:
                resu._restrictions[dom] = self_r._restrictions[dom].wedge(
                                          other_r._restrictions[dom])
        return resu

    def interior_product(self, form):
        r"""
        Interior product with a differential form.

        If ``self`` is a multivector field `A` of degree `p` and `B` is a
        differential form of degree `q\geq p` on the same manifold as `A`, the
        interior product of `A` by `B` is the differential form `\iota_A B` of
        degree `q-p` defined by

        .. MATH::

            (\iota_A B)_{i_1\ldots i_{q-p}} = A^{k_1\ldots k_p}
                            B_{k_1\ldots k_p i_1\ldots i_{q-p}}

        .. NOTE::

            ``A.interior_product(B)`` yields the same result as
            ``A.contract(0,..., p-1, B, 0,..., p-1)`` (cf.
            :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.contract`),
            but ``interior_product`` is more efficient, the alternating
            character of `A` being not used to reduce the computation in
            :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.contract`

        INPUT:

        - ``form`` -- differential form `B` (instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`);
          the degree of `B` must be at least equal to the degree of ``self``

        OUTPUT:

        - scalar field (case `p=q`) or
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
          (case `p<q`) representing the interior product `\iota_A B`, where `A`
          is ``self``

        .. SEEALSO::

            :meth:`~sage.manifolds.differentiable.diff_form.DiffForm.interior_product`
            for the interior product of an differential form by a multivector
            field

        EXAMPLES:

        Interior product of a vector field (`p=1`) by a 2-form (`q=2`) on the
        2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1) # the sphere S^2
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() # stereographic coord. North
            sage: c_uv.<u,v> = V.chart() # stereographic coord. South
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: a = M.vector_field(name='a')
            sage: a[e_xy,:] = y, x
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: b = M.diff_form(2, name='b')
            sage: b[e_xy,1,2] = 4/(x^2+y^2+1)^2   # the standard area 2-form
            sage: b.add_comp_by_continuation(e_uv, W, c_uv)
            sage: b.display(e_xy)
            b = 4/(x^2 + y^2 + 1)^2 dx/\dy
            sage: b.display(e_uv)
            b = -4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) du/\dv
            sage: s = a.interior_product(b); s
            1-form i_a b on the 2-dimensional differentiable manifold S^2
            sage: s.display(e_xy)
            i_a b = -4*x/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) dx
             + 4*y/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1) dy
            sage: s.display(e_uv)
            i_a b = 4*(u^3 - 3*u*v^2)/(u^6 + v^6 + (3*u^2 + 2)*v^4 + 2*u^4
             + (3*u^4 + 4*u^2 + 1)*v^2 + u^2) du
             + 4*(3*u^2*v - v^3)/(u^6 + v^6 + (3*u^2 + 2)*v^4 + 2*u^4
             + (3*u^4 + 4*u^2 + 1)*v^2 + u^2) dv
            sage: s == a.contract(b)
            True

        Example with `p=2` and `q=2`::

            sage: a = M.multivector_field(2, name='a')
            sage: a[e_xy,1,2] = x*y
            sage: a.add_comp_by_continuation(e_uv, W, c_uv)
            sage: a.display(e_xy)
            a = x*y d/dx/\d/dy
            sage: a.display(e_uv)
            a = -u*v d/du/\d/dv
            sage: s = a.interior_product(b); s
            Scalar field i_a b on the 2-dimensional differentiable manifold S^2
            sage: s.display()
            i_a b: S^2 --> R
            on U: (x, y) |--> 8*x*y/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1)
            on V: (u, v) |--> 8*u*v/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1)

        Some checks::

            sage: s == a.contract(0, 1, b, 0, 1)
            True
            sage: s.restrict(U) == 2 * a[[e_xy,1,2]] * b[[e_xy,1,2]]
            True
            sage: s.restrict(V) == 2 * a[[e_uv,1,2]] * b[[e_uv,1,2]]
            True

        """
        from sage.tensor.modules.format_utilities import is_atomic
        if self._domain.is_subset(form._domain):
            if not self._ambient_domain.is_subset(form._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        elif form._domain.is_subset(self._domain):
            if not form._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        dom_resu = self._domain.intersection(form._domain)
        ambient_dom_resu = self._ambient_domain.intersection(form._ambient_domain)
        self_r = self.restrict(dom_resu)
        form_r = form.restrict(dom_resu)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the AlternatingContrTensor version:
            return AlternatingContrTensor.interior_product(self_r, form_r)
        # Otherwise, the result is created here:
        # Name of the result
        resu_name = None
        if self._name is not None and form._name is not None:
            sname = self._name
            oname = form._name
            if not is_atomic(sname):
                sname = '(' + sname + ')'
            if not is_atomic(oname):
                oname = '(' + oname + ')'
            resu_name = 'i_' + sname + ' ' + oname
        resu_latex_name = None
        if self._latex_name is not None and form._latex_name is not None:
            slname = self._latex_name
            olname = form._latex_name
            if not is_atomic(olname):
                olname = r'\left(' + olname + r'\right)'
            resu_latex_name = r'\iota_{' + slname + '} ' + olname
        # Domain and computation of the result
        dest_map = self._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu,
                                          subcodomain=ambient_dom_resu)
        vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)
        resu_degree = form._tensor_rank - self._tensor_rank
        resu = vmodule.alternating_form(resu_degree,
                                    name=resu_name, latex_name=resu_latex_name)
        for dom in self_r._restrictions:
            if dom in form_r._restrictions:
                resu._restrictions[dom] = \
                    self_r._restrictions[dom].interior_product(
                                                     form_r._restrictions[dom])
        if resu_degree == 0:
            if not resu._express: # only the restrictions to subdomains have
                                  # been initialized
                for chart in dom_resu.top_charts():
                    resu._express[chart] = \
                            resu.restrict(chart.domain()).coord_function(chart)
        return resu

    def degree(self):
        r"""
        Return the degree of ``self``.

        OUTPUT:

        - integer `p` such that ``self`` is a `p`-vector field

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: a = M.multivector_field(2); a
            2-vector field on the 3-dimensional differentiable manifold M
            sage: a.degree()
            2
            sage: b = M.vector_field(); b
            Vector field on the 3-dimensional differentiable manifold M
            sage: b.degree()
            1

        """
        return self._tensor_rank


#******************************************************************************

class MultivectorFieldParal(AlternatingContrTensor, TensorFieldParal):
    r"""
    Multivector field with values on a parallelizable manifold.

    Given a differentiable manifold `U`,  a differentiable map
    `\Phi: U \rightarrow M` to a parallelizable manifold `M` and a positive
    integer `p`, a *multivector field of degree* `p` (or `p`-*vector field*)
    *along* `U` *with values on* `M\supset\Phi(U)` is a differentiable map

    .. MATH::

        a:\ U  \longrightarrow T^{(p,0)}M

    (`T^{(p,0)}M` being the tensor bundle of type `(p,0)` over `M`) such that

    .. MATH::

        \forall x \in U,\quad a(x) \in \Lambda^p(T_{\Phi(x)} M) ,

    where `T_{\Phi(x)} M` is the vector space tangent to `M` at `\Phi(x)` and
    `\Lambda^p` stands for the exterior power of degree `p` (cf.
    :class:`~sage.tensor.modules.ext_pow_free_module.ExtPowerFreeModule`).
    In other words, `a(x)` is an alternating contravariant tensor of degree `p`
    of the tangent vector space `T_{\Phi(x)} M`.

    The standard case of a multivector field *on* a manifold `M` corresponds to
    `U = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval of
    `\RR`).

    .. NOTE::

        If `M` is not parallelizable, the class :class:`MultivectorField` must
        be used instead.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` with values on `M` via the map `\Phi`
    - ``degree`` -- the degree of the multivector field (i.e. its tensor rank)
    - ``name`` -- (default: ``None``) name given to the multivector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      multivector field; if none is provided, the LaTeX symbol is set to
      ``name``

    EXAMPLES:

    A 2-vector field on a 4-dimensional manifold::

        sage: M = Manifold(4, 'M')
        sage: c_txyz.<t,x,y,z> = M.chart()
        sage: a = M.multivector_field(2, 'a') ; a
        2-vector field a on the 4-dimensional differentiable manifold M
        sage: a.parent()
        Free module A^2(M) of 2-vector fields on the 4-dimensional
         differentiable manifold M

    A multivector field is a tensor field of purely contravariant type::

        sage: a.tensor_type()
        (2, 0)

    It is antisymmetric, its components being
    :class:`~sage.tensor.modules.comp.CompFullyAntiSym`::

        sage: a.symmetries()
        no symmetry;  antisymmetry: (0, 1)
        sage: a[0,1] = 2
        sage: a[1,0]
        -2
        sage: a.comp()
        Fully antisymmetric 2-indices components w.r.t. Coordinate frame
         (M, (d/dt,d/dx,d/dy,d/dz))
        sage: type(a.comp())
        <class 'sage.tensor.modules.comp.CompFullyAntiSym'>

    Setting a component with repeated indices to a non-zero value results in
    an error::

        sage: a[1,1] = 3
        Traceback (most recent call last):
        ...
        ValueError: by antisymmetry, the component cannot have a nonzero value
         for the indices (1, 1)
        sage: a[1,1] = 0  # OK, albeit useless
        sage: a[1,2] = 3  # OK

    The expansion of a multivector field with respect to a given frame is
    displayed via the method
    :meth:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor.display`::

        sage: a.display() # expansion w.r.t. the default frame
        a = 2 d/dt/\d/dx + 3 d/dx/\d/dy
        sage: latex(a.display()) # output for the notebook
        a = 2 \frac{\partial}{\partial t }\wedge \frac{\partial}{\partial x }
         + 3 \frac{\partial}{\partial x }\wedge \frac{\partial}{\partial y }

    Multivector fields can be added or subtracted::

        sage: b = M.multivector_field(2)
        sage: b[0,1], b[0,2], b[0,3] = (1,2,3)
        sage: s = a + b ; s
        2-vector field on the 4-dimensional differentiable manifold M
        sage: a[:], b[:], s[:]
        (
        [ 0  2  0  0]  [ 0  1  2  3]  [ 0  3  2  3]
        [-2  0  3  0]  [-1  0  0  0]  [-3  0  3  0]
        [ 0 -3  0  0]  [-2  0  0  0]  [-2 -3  0  0]
        [ 0  0  0  0], [-3  0  0  0], [-3  0  0  0]
        )
        sage: s = a - b ; s
        2-vector field on the 4-dimensional differentiable manifold M
        sage: s[:]
        [ 0  1 -2 -3]
        [-1  0  3  0]
        [ 2 -3  0  0]
        [ 3  0  0  0]

    An example of 3-vector field in `\RR^3` with Cartesian coordinates::

        sage: M = Manifold(3, 'R3', '\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart()
        sage: a = M.multivector_field(3, name='a')
        sage: a[1,2,3] = x^2+y^2+z^2  # the only independent component
        sage: a[:] # all the components are set from the previous line:
        [[[0, 0, 0], [0, 0, x^2 + y^2 + z^2], [0, -x^2 - y^2 - z^2, 0]],
         [[0, 0, -x^2 - y^2 - z^2], [0, 0, 0], [x^2 + y^2 + z^2, 0, 0]],
         [[0, x^2 + y^2 + z^2, 0], [-x^2 - y^2 - z^2, 0, 0], [0, 0, 0]]]
        sage: a.display()
        a = (x^2 + y^2 + z^2) d/dx/\d/dy/\d/dz

    Spherical components from the tensorial change-of-frame formula::

        sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: spher_to_cart = c_spher.transition_map(c_cart,
        ....:                [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
        sage: cart_to_spher = spher_to_cart.set_inverse(sqrt(x^2+y^2+z^2),
        ....:                              atan2(sqrt(x^2+y^2),z), atan2(y, x))
        sage: a.comp(c_spher.frame()) # computation of the components in the spherical frame
        Fully antisymmetric 3-indices components w.r.t. Coordinate frame
         (R3, (d/dr,d/dth,d/dph))
        sage: a.comp(c_spher.frame())[1,2,3, c_spher]
        1/sin(th)
        sage: a.display(c_spher.frame())
        a = sqrt(x^2 + y^2 + z^2)/sqrt(x^2 + y^2) d/dr/\d/dth/\d/dph
        sage: a.display(c_spher.frame(), c_spher)
        a = 1/sin(th) d/dr/\d/dth/\d/dph

    The exterior product of two multivector fields is performed via the method
    :meth:`~sage.tensor.modules.alternating_contr_tensor.AlternatingContrTensor.wedge`::

        sage: a = M.vector_field(name='A')
        sage: a[:] = (x*y*z, -z*x, y*z)
        sage: b = M.vector_field(name='B')
        sage: b[:] = (y, z+y^2, x^2-z^2)
        sage: ab = a.wedge(b) ; ab
        2-vector field A/\B on the 3-dimensional differentiable manifold R3
        sage: ab[:]
        [                             0      x*y*z^2 + (x*y^3 + x*y)*z     -x*y*z^3 + (x^3*y - y^2)*z]
        [    -x*y*z^2 - (x*y^3 + x*y)*z                              0  x*z^3 - y*z^2 - (x^3 + y^3)*z]
        [     x*y*z^3 - (x^3*y - y^2)*z -x*z^3 + y*z^2 + (x^3 + y^3)*z                              0]
        sage: ab.display()
        A/\B = (x*y*z^2 + (x*y^3 + x*y)*z) d/dx/\d/dy + (-x*y*z^3
         + (x^3*y - y^2)*z) d/dx/\d/dz + (x*z^3 - y*z^2
         - (x^3 + y^3)*z) d/dy/\d/dz

    Let us check the formula relating the exterior product to the tensor
    product for vector fields::

        sage: a.wedge(b) == a*b - b*a
        True

    The tensor product of a vector field and a 2-vector field is not a 3-vector
    field but a tensor field of type `(3,0)` with less symmetries::

        sage: c = a*ab ; c
        Tensor field A*(A/\B) of type (3,0) on the 3-dimensional differentiable
         manifold R3
        sage: c.symmetries()  # the antisymmetry is only w.r.t. the last 2 arguments:
        no symmetry;  antisymmetry: (1, 2)
        sage: d = ab*a ; d
        Tensor field (A/\B)*A of type (3,0) on the 3-dimensional differentiable
         manifold R3
        sage: d.symmetries()  # the antisymmetry is only w.r.t. the first 2 arguments:
        no symmetry;  antisymmetry: (0, 1)

    The Lie derivative of a 2-vector field is a 2-vector field::

        sage: v = M.vector_field(name='v')
        sage: v[:] = (y*z, -x*z, x*y)
        sage: ab.lie_der(v)
        2-vector field on the 3-dimensional differentiable manifold R3

    """
    def __init__(self, vector_field_module, degree, name=None,
                 latex_name=None):
        r"""
        Construct a multivector field.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``MultivectorFieldParal``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: A = M.multivector_module(2)
            sage: XM = M.vector_field_module()
            sage: a = A.element_class(XM, 2, name='a'); a
            2-vector field a on the 2-dimensional differentiable manifold M
            sage: a[0,1] = x*y
            sage: TestSuite(a).run()

        Construction via ``DifferentiableManifold.multivector_field``::

            sage: a1 = M.multivector_field(2, name='a'); a1
            2-vector field a on the 2-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True

        """
        AlternatingContrTensor.__init__(self, vector_field_module, degree,
                                        name=name, latex_name=latex_name)
        # TensorFieldParal attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.multivector_field(2, name='a')
            sage: a._repr_()
            '2-vector field a on the 3-dimensional differentiable manifold M'
            sage: repr(a)  # indirect doctest
            '2-vector field a on the 3-dimensional differentiable manifold M'
            sage: a  # indirect doctest
            2-vector field a on the 3-dimensional differentiable manifold M
            sage: b = M.multivector_field(2)
            sage: b._repr_()
            '2-vector field on the 3-dimensional differentiable manifold M'

        """
        description = "{}-vector field ".format(self._tensor_rank)
        if self._name is not None:
            description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class, of the same degree and on the
        same domain.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.multivector_field(2, name='a')
            sage: a1 = a._new_instance(); a1
            2-vector field on the 3-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True
            sage: a1.parent() is a.parent()
            True

        """
        return type(self)(self._fmodule, self._tensor_rank)

    # This method is needed to redirect to the correct class (TensorFieldParal)
    def _init_derived(self):
        r"""
        Initialize the derived quantities of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.multivector_field(2, name='a')
            sage: a._init_derived()

        """
        TensorFieldParal._init_derived(self)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- (default: ``True``) determines whether the
          restrictions of ``self`` to subdomains are deleted

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()  # makes M parallelizable
            sage: a = M.multivector_field(2, name='a')
            sage: a._del_derived()

        """
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)

    def __call__(self, *args):
        r"""
        Redefinition of
        :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.__call__`
        to allow for domain treatment.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.multivector_field(2, name='a')
            sage: a[0,1] = x*y
            sage: a.display()
            a = x*y d/dx/\d/dy
            sage: b = M.one_form(name='b')
            sage: b[:] = [1+x, 2-y]
            sage: c = M.one_form(name='c')
            sage: c[:] = [-y, x]
            sage: s = a.__call__(b,c); s
            Scalar field a(b,c) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(b,c): M --> R
               (x, y) |--> -x*y^3 + 2*x*y^2 + (x^3 + x^2)*y
            sage: s == a[[0,1]]*(b[[0]]*c[[1]] - b[[1]]*c[[0]])
            True
            sage: s == a(b,c)  # indirect doctest
            True

        """
        return TensorFieldParal.__call__(self, *args)

    def wedge(self, other):
        r"""
        Exterior product of ``self`` with another multivector field.

        INPUT:

        - ``other`` -- another multivector field

        OUTPUT:

        - instance of :class:`MultivectorFieldParal` representing the exterior
          product ``self/\other``

        EXAMPLES:

        Exterior product of a vector field and a 2-vector field on a
        3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: a = M.vector_field(name='a')
            sage: a[:] = [2, 1+x, y*z]
            sage: b = M.multivector_field(2, name='b')
            sage: b[1,2], b[1,3], b[2,3] = y^2, z+x, z^2
            sage: a.display()
            a = 2 d/dx + (x + 1) d/dy + y*z d/dz
            sage: b.display()
            b = y^2 d/dx/\d/dy + (x + z) d/dx/\d/dz + z^2 d/dy/\d/dz
            sage: s = a.wedge(b); s
            3-vector field a/\b on the 3-dimensional differentiable manifold M
            sage: s.display()
            a/\b = (-x^2 + (y^3 - x - 1)*z + 2*z^2 - x) d/dx/\d/dy/\d/dz

        Check::

            sage: s[1,2,3] == a[1]*b[2,3] + a[2]*b[3,1] + a[3]*b[1,2]
            True

        """
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior " +
                                 "product")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for exterior " +
                                 "product")
        dom_resu = self._domain.intersection(other._domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        return AlternatingContrTensor.wedge(self_r, other_r)

    def interior_product(self, form):
        r"""
        Interior product with a differential form.

        If ``self`` is a multivector field `A` of degree `p` and `B` is a
        differential form of degree `q\geq p` on the same manifold as `A`, the
        interior product of `A` by `B` is the differential form `\iota_A B` of
        degree `q-p` defined by

        .. MATH::

            (\iota_A B)_{i_1\ldots i_{q-p}} = A^{k_1\ldots k_p}
                            B_{k_1\ldots k_p i_1\ldots i_{q-p}}

        .. NOTE::

            ``A.interior_product(B)`` yields the same result as
            ``A.contract(0,..., p-1, B, 0,..., p-1)`` (cf.
            :meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.contract`),
            but ``interior_product`` is more efficient, the alternating
            character of `A` being not used to reduce the computation in
            :meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.contract`

        INPUT:

        - ``form`` -- differential form `B` (instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`);
          the degree of `B` must be at least equal to the degree of ``self``

        OUTPUT:

        - scalar field (case `p=q`) or
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`
          (case `p<q`) representing the interior product `\iota_A B`, where `A`
          is ``self``

        .. SEEALSO::

            :meth:`~sage.manifolds.differentiable.diff_form.DiffFormParal.interior_product`
            for the interior product of an differential form by a multivector
            field

        EXAMPLES:

        Interior product with `p=1` and `q=1` on 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: a = M.vector_field(name='a')
            sage: a[:] = [x, 1+t^2, x*z, y-3]
            sage: b = M.one_form(name='b')
            sage: b[:] = [-z^2, 2, x*y, x-y]
            sage: s = a.interior_product(b); s
            Scalar field i_a b on the 4-dimensional differentiable manifold M
            sage: s.display()
            i_a b: M --> R
               (t, x, y, z) |--> x^2*y*z - x*z^2 + 2*t^2 + (x + 3)*y - y^2 - 3*x + 2

        In this case, we have `\iota_a b = a^i b_i = a(b) = b(a)`::

            sage: all([s == a.contract(b), s == a(b), s == b(a)])
            True

        Case `p=1` and `q=3`::

            sage: c = M.diff_form(3, name='c')
            sage: c[0,1,2], c[0,1,3] = x*y - z, -3*t
            sage: c[0,2,3], c[1,2,3] = t^2+x, y^2
            sage: s = a.interior_product(c); s
            2-form i_a c on the 4-dimensional differentiable manifold M
            sage: s.display()
            i_a c = (x^2*y*z - x*z^2 - 3*t*y + 9*t) dt/\dx + (-3*t^2
             - (t^2*x - t^2)*y + (t^2 + 1)*z - 3*x) dt/\dy + (3*t^3
             - (t^2*x + x^2)*z + 3*t) dt/\dz + (x^2*y + y^3 - 3*y^2
             - x*z) dx/\dy + (-x*y^2*z - 3*t*x) dx/\dz + (t^2*x
             + (t^2 + 1)*y^2 + x^2) dy/\dz
            sage: s == a.contract(c)
            True

        Case `p=2` and `q=3`::

            sage: d = M.multivector_field(2, name='d')
            sage: d[0,1], d[0,2], d[0,3] = t-x, 2*z, y-1
            sage: d[1,2], d[1,3], d[2,3] = z^2, y*t, 4
            sage: s = d.interior_product(c); s
            1-form i_d c on the 4-dimensional differentiable manifold M
            sage: s.display()
            i_d c = (2*x*y*z^2 - 6*t^2*y - 2*z^3 + 8*t^2 + 8*x) dt
             + (-4*x*y*z + 6*t*y + 8*y^2 + 4*z^2 - 6*t) dx
             + (-2*t*y^3 + 2*t^2 - 2*(t^2 - (t - 1)*x + x^2)*y - 2*(t - x)*z
             + 2*x) dy + (2*y^2*z^2 - 6*t^2 + 6*t*x + 4*(t^2 + x)*z) dz
            sage: s == d.contract(0, 1, c, 0, 1)
            True

        """
        if self._domain.is_subset(form._domain):
            if not self._ambient_domain.is_subset(form._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        elif form._domain.is_subset(self._domain):
            if not form._ambient_domain.is_subset(self._ambient_domain):
                raise ValueError("incompatible ambient domains for interior " +
                                 "product")
        dom_resu = self._domain.intersection(form._domain)
        self_r = self.restrict(dom_resu)
        form_r = form.restrict(dom_resu)
        return AlternatingContrTensor.interior_product(self_r, form_r)

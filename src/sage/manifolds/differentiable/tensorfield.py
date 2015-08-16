r"""
Tensor fields

The class :class:`TensorField` implements tensor fields on differentiable
manifolds.
The derived class :class:`TensorFieldParal` is devoted to tensor fields with
values on parallelizable open subsets.

Various derived classes of :class:`TensorField` are devoted to specific tensor
fields:

* :class:`~sage.manifolds.differentiable.vectorfield.VectorField` for vector
  fields (rank-1 contravariant tensor fields)

* :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
  for fields of tangent-space automorphisms

* :class:`~sage.manifolds.differentiable.diff_form.DiffForm` for differential
  forms (fully antisymmetric covariant tensor fields)


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

REFERENCES:

- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
  Interscience Publishers (New York) (1963)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)
- B O'Neill : *Semi-Riemannian Geometry*, Academic Press (San Diego) (1983)

EXAMPLES:

A tensor field of type (1,1) on a 2-dimensional differentiable manifold::

    sage: M = DiffManifold(2, 'M', start_index=1)
    sage: c_xy.<x,y> = M.chart()
    sage: t = M.tensor_field(1, 1, 'T') ; t
    Tensor field T of type (1,1) on the 2-dimensional differentiable manifold M
    sage: t.tensor_type()
    (1, 1)
    sage: t.tensor_rank()
    2

Components w.r.t. the manifold's default frame are created by providing the
relevant indices inside square brackets::

    sage: t[1,1] = x^2

Unset components are initialized to zero::

    sage: t[:]  # list of components w.r.t. the manifold's default vector frame
    [x^2   0]
    [  0   0]

The full set of components w.r.t. a given vector frame is returned by the
method
:meth:`~sage.manifolds.differentiable.tensorfield.TensorFieldParal.comp`; it is
an instance of the class :class:`~sage.tensor.modules.comp.Components`::

    sage: t.comp(c_xy.frame())
    2-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
    sage: type(t.comp(c_xy.frame()))
    <class 'sage.tensor.modules.comp.Components'>

If no vector frame is mentionned in the argument of
:meth:`~sage.manifolds.differentiable.tensorfield.TensorFieldParal.comp`, it is
assumed to be the manifold's default frame::

    sage: M.default_frame()
    Coordinate frame (M, (d/dx,d/dy))
    sage: t.comp() is t.comp(c_xy.frame())
    True

Individual components w.r.t. the manifold's default frame are accessed by
listing their indices inside double square brackets; they are scalar
fields on the manifold, and therefore instances of the class
:class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`::

    sage: t[[1,1]]
    Scalar field on the 2-dimensional differentiable manifold M
    sage: t[[1,1]].display()
    M --> R
    (x, y) |--> x^2
    sage: t[[1,2]]
    Scalar field zero on the 2-dimensional differentiable manifold M
    sage: t[[1,2]].display()
    zero: M --> R
       (x, y) |--> 0

A direct access to the coordinate expression of some component is obtained
via the single square brackets::

    sage: t[1,1]
    x^2
    sage: t[1,1] is t[[1,1]].coord_function() # the coordinate function
    True
    sage: t[1,1] is t[[1,1]].coord_function(c_xy)
    True
    sage: t[1,1].expr() is t[[1,1]].expr() # the symbolic expression
    True

In other words, the single square brackets return an instance of
:class:`~sage.manifolds.coord_func.CoordFunction` that is the
coordinate function representing the component in some chart (by default,
the manifold's default chart)::

    sage: type(t[1,1])    # single bracket --> coordinate function
    <class 'sage.manifolds.coord_func_symb.CoordFunctionSymb'>
    sage: type(t[[1,1]])  # double bracket --> scalar field
    <class 'sage.manifolds.differentiable.scalarfield.DiffScalarFieldAlgebra_with_category.element_class'>

Expressions in a chart different from the manifold's default one are
obtained by specifying the chart as the last argument inside the
single square brackets::

    sage: c_uv.<u,v> = M.chart()
    sage: xy_to_uv = c_xy.transition_map(c_uv, [x+y, x-y])
    sage: uv_to_xy = xy_to_uv.inverse()
    sage: t[1,1, c_uv]
    1/4*u^2 + 1/2*u*v + 1/4*v^2

Note that ``t[1,1, c_uv]`` is the component of the tensor t w.r.t. to
the coordinate frame associated to the chart (x,y) expressed in terms of
the coordinates (u,v). Indeed, ``t[1,1, c_uv]`` is a shortcut for
``t.comp(c_xy.frame())[[1,1]].coord_function(c_uv)``::

    sage: t[1,1, c_uv] is t.comp(c_xy.frame())[[1,1]].coord_function(c_uv)
    True

Similarly, ``t[1,1]`` is a shortcut for
``t.comp(c_xy.frame())[[1,1]].coord_function(c_xy)``::

    sage: t[1,1] is t.comp(c_xy.frame())[[1,1]].coord_function(c_xy)
    True
    sage: t[1,1] is t.comp()[[1,1]].coord_function()  # since c_xy.frame() and c_xy are the manifold's default values
    True

Internally, the components are stored as a dictionary (attribute
:attr:`_comp` of the class
:class:`~sage.tensor.modules.comp.Components`) whose
keys are the indices. Only the non-zero components and non-redundant
components (in case of symmetries) are stored::

    sage: t.comp()._comp
    {(1, 1): Scalar field on the 2-dimensional differentiable manifold M}

All the components can be set at once via [:]::

    sage: t[:] = [[1, -x], [x*y, 2]]
    sage: t[:]
    [  1  -x]
    [x*y   2]

To set the components in a vector frame different from the manifold's
default one, the method
:meth:`~sage.manifolds.differentiable.tensorfield.TensorFieldParal.set_comp`
can be employed::

    sage: e = M.vector_frame('e')
    sage: t.set_comp(e)[1,1] = x+y
    sage: t.set_comp(e)[2,1], t.set_comp(e)[2,2] = y, -3*x

but, as a shortcut, one may simply specify the frame as the first argument
of the square brackets::

    sage: t[e,1,1] = x+y
    sage: t[e,2,1], t[e,2,2] = y, -3*x
    sage: t.comp(e)
    2-indices components w.r.t. Vector frame (M, (e_1,e_2))
    sage: t.comp(e)[:]
    [x + y     0]
    [    y  -3*x]
    sage: t[e,:]  # a shortcut of the above
    [x + y     0]
    [    y  -3*x]

All the components in some frame can be set at once, via the operator
[:]::

    sage: t[e,:] = [[x+y, 0], [y, -3*x]]
    sage: t[e,:]  # same as above:
    [x + y     0]
    [    y  -3*x]

To avoid any insconstency between the various components, the method
:meth:`~sage.manifolds.differentiable.tensorfield.TensorFieldParal.set_comp`
clears the components in other frames.
To keep the other components, one must use the method
:meth:`~sage.manifolds.differentiable.tensorfield.TensorFieldParal.add_comp`::

    sage: t = M.tensor_field(1, 1, 'T')  # Let us restart
    sage: t[:] = [[1, -x], [x*y, 2]]  # by first setting the components in the frame c_xy.frame()
    sage: # We now set the components in the frame e with add_comp:
    sage: t.add_comp(e)[:] = [[x+y, 0], [y, -3*x]]

The expansion of the tensor field in a given frame is obtained via the
method
:meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.display` (the
symbol * stands for tensor product)::

    sage: t.display()  # expansion in the manifold's default frame
    T = d/dx*dx - x d/dx*dy + x*y d/dy*dx + 2 d/dy*dy
    sage: t.display(e)
    T = (x + y) e_1*e^1 + y e_2*e^1 - 3*x e_2*e^2

By definition, a tensor field acts as a multilinear map on 1-forms and vector
fields; in the present case, T being of type (1,1), it acts on pairs
(1-form, vector field)::

    sage: a = M.one_form('a')
    sage: a[:] = (1, x)
    sage: v = M.vector_field('V')
    sage: v[:] = (y, 2)
    sage: t(a,v)
    Scalar field T(a,V) on the 2-dimensional differentiable manifold M
    sage: t(a,v).display()
    T(a,V): M --> R
       (x, y) |--> x^2*y^2 + 2*x + y
       (u, v) |--> 1/16*u^4 - 1/8*u^2*v^2 + 1/16*v^4 + 3/2*u + 1/2*v
    sage: latex(t(a,v))
    T\left(a,V\right)

Check by means of the component expression of t(a,v)::

    sage: t(a,v).expr() - t[1,1]*a[1]*v[1] - t[1,2]*a[1]*v[2] \
    ....: - t[2,1]*a[2]*v[1] - t[2,2]*a[2]*v[2]
    0

A scalar field (rank-0 tensor field)::

    sage: f = M.scalar_field(x*y + 2, name='f') ; f
    Scalar field f on the 2-dimensional differentiable manifold M
    sage: f.tensor_type()
    (0, 0)

A scalar field acts on points on the manifold::

    sage: p = M.point((1,2))
    sage: f(p)
    4

See :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField` for
more details on scalar fields.

A vector field (rank-1 contravariant tensor field)::

    sage: v = M.vector_field('v') ; v
    Vector field v on the 2-dimensional differentiable manifold M
    sage: v.tensor_type()
    (1, 0)
    sage: v[1], v[2] = -x, y
    sage: v.display()
    v = -x d/dx + y d/dy

A field of symmetric bilinear forms::

    sage: q = M.sym_bilin_form_field('Q') ; q
    Field of symmetric bilinear forms Q on the 2-dimensional differentiable
     manifold M
    sage: q.tensor_type()
    (0, 2)

The components of a symmetric bilinear form are dealt by the subclass
:class:`~sage.tensor.modules.comp.CompFullySym` of the class
:class:`~sage.tensor.modules.comp.Components`, which takes into
account the symmetry between the two indices::

    sage: q[1,1], q[1,2], q[2,2] = (0, -x, y) # no need to set the component (2,1)
    sage: type(q.comp())
    <class 'sage.tensor.modules.comp.CompFullySym'>
    sage: q[:] # note that the component (2,1) is equal to the component (1,2)
    [ 0 -x]
    [-x  y]
    sage: q.display()
    Q = -x dx*dy - x dy*dx + y dy*dy

Internally (dictionary :attr:`_comp` of the class
:class:`~sage.tensor.modules.comp.Components`), only
the non-zero and non-redundant components are stored::

    sage: q.comp()._comp  # random (dictionary output)
    {(1, 2): Scalar field on the 2-dimensional differentiable manifold M,
     (2, 2): Scalar field on the 2-dimensional differentiable manifold M}
    sage: q.comp()._comp[(1,2)].expr()
    -x
    sage: q.comp()._comp[(2,2)].expr()
    y

More generally, tensor symmetries or antisymmetries can be specified via
the keywords ``sym`` and ``antisym``. For instance a rank-4 covariant
tensor symmetric with respect to its first two arguments (no. 0 and no. 1) and
antisymmetric with respect to its last two ones (no. 2 and no. 3) is declared
as follows::

    sage: t = M.tensor_field(0, 4, 'T', sym=(0,1), antisym=(2,3))
    sage: t[1,2,1,2] = 3
    sage: t[2,1,1,2] # check of the symmetry with respect to the first 2 indices
    3
    sage: t[1,2,2,1] # check of the antisymmetry with respect to the last 2 indices
    -3

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
from sage.tensor.modules.free_module_tensor import FreeModuleTensor
from sage.tensor.modules.tensor_with_indices import TensorWithIndices

class TensorField(ModuleElement):
    r"""
    Tensor field along an open set of a differentiable manifold.

    An instance of this class is a tensor field along an open subset `U`
    of some differentiable manifold `S` with values in an open subset `V`
    of a differentiable manifold `M`, via a differentiable map
    `\Phi: U \rightarrow V`. The standard case of a tensor field *on* a
    differentiable manifold corresponds to `S=M`, `U=V` and
    `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `V` (`U` is then an open interval
    of `\RR`).

    If `V` is parallelizable, the class :class:`TensorFieldParal` should be
    used instead.

    A tensor field of type `(k,l)` is a field `t` on `U`, such that at each
    point `p\in U`, `t(p)` is a multilinear map

    .. MATH::

        t(p):\ \underbrace{T_q^*M\times\cdots\times T_q^*M}_{k\ \; \mbox{times}}
        \times \underbrace{T_q M\times\cdots\times T_q M}_{l\ \; \mbox{times}}
        \longrightarrow \RR

    where `T_q M` stands for the tangent space to the manifold `M` at the point
    `q=\Phi(p)` and `T_q^* M` for its dual vector space. The integer `k+l`
    is called the tensor rank.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the map `\Phi: U \rightarrow V` (cf.
      :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`)
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank
      and `l` the covariant rank
    - ``name`` -- (default: ``None``) name given to the tensor field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the tensor
      field; if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the tensor arguments: each symmetry is described by a tuple containing
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:

      * sym=(0,1) for a symmetry between the 1st and 2nd arguments
      * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
        arguments and a symmetry between the 2nd, 4th and 5th arguments.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the arguments, with the same convention as for ``sym``.
    - ``parent`` -- (default: ``None``) some specific parent (e.g. exterior
      power for differential forms); if ``None``,
      ``vector_field_module.tensor_module(k,l)`` is used

    EXAMPLES:

    Tensor field of type (0,2) on the sphere `S^2`::

        sage: M = DiffManifold(2, 'S^2') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: W = U.intersection(V)
        sage: t = M.tensor_field(0,2, name='t') ; t
        Tensor field t of type (0,2) on the 2-dimensional differentiable
         manifold S^2
        sage: t.parent()
        Module T^(0,2)(S^2) of type-(0,2) tensors fields on the 2-dimensional
         differentiable manifold S^2
        sage: t.parent().category()
        Category of modules over Algebra of differentiable scalar fields on the
         2-dimensional differentiable manifold S^2

    The parent of `t` is not a free module, for the sphere `S^2` is not
    parallelizable::

        sage: isinstance(t.parent(), FiniteRankFreeModule)
        False

    To fully define `t`, we have to specify its components in some vector
    frames defined on subsets of `S^2`; let us start by the open subset `U`::

        sage: eU = c_xy.frame()
        sage: t = M.tensor_field(0,2, name='t')
        sage: t[eU,:] = [[1,0], [-2,3]]
        sage: t.display(eU)
        t = dx*dx - 2 dy*dx + 3 dy*dy

    To set the components of `t` on `V` consistently, we copy the expressions
    of the components in the common subset `W`::

        sage: eV = c_uv.frame()
        sage: eVW = eV.restrict(W)
        sage: c_uvW = c_uv.restrict(W)
        sage: t[eV,0,0] = t[eVW,0,0,c_uvW].expr()
        sage: t[eV,0,1] = t[eVW,0,1,c_uvW].expr()
        sage: t[eV,1,0] = t[eVW,1,0,c_uvW].expr()
        sage: t[eV,1,1] = t[eVW,1,1,c_uvW].expr()

    Actually, the above operation can by performed in a single line by means
    of the method
    :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.add_comp_by_continuation`::

        sage: t.add_comp_by_continuation(eV, W, chart=c_uv)

    At this stage, `t` is fully defined, having components in frames eU and eV
    and the union of the domains of eU and eV being whole manifold::

        sage: t.display(eV)
        t = (u^4 - 4*u^3*v + 10*u^2*v^2 + 4*u*v^3 + v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du*du
         - 4*(u^3*v + 2*u^2*v^2 - u*v^3)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) du*dv
         + 2*(u^4 - 2*u^3*v - 2*u^2*v^2 + 2*u*v^3 + v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv*du
         + (3*u^4 + 4*u^3*v - 2*u^2*v^2 - 4*u*v^3 + 3*v^4)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) dv*dv

    Let us consider two vector fields, `a` and `b`, on `S^2`::

        sage: a = M.vector_field(name='a')
        sage: a[eU,:] = [1,x]
        sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: a.display(eV)
        a = -(u^4 - v^4 + 2*u^2*v)/(u^2 + v^2) d/du
         - (2*u^3*v + 2*u*v^3 - u^3 + u*v^2)/(u^2 + v^2) d/dv
        sage: b = M.vector_field(name='b')
        sage: b[eU,:] = [y,-1]
        sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
        sage: b.display(eV)
        b = ((2*u + 1)*v^3 + (2*u^3 - u^2)*v)/(u^2 + v^2) d/du
         - (u^4 - v^4 + 2*u*v^2)/(u^2 + v^2) d/dv

    As a tensor field of type (0,2), `t` acts on the pair `(a,b)`, resulting in
    a scalar field::

        sage: f = t(a,b) ; f
        Scalar field t(a,b) on the 2-dimensional differentiable manifold S^2
        sage: f.display()
        t(a,b): S^2 --> R
        on U: (x, y) |--> -(2*x - 1)*y - 3*x
        on V: (u, v) |--> -(3*u^3 + 3*u*v^2 - v^3 - (u^2 - 2*u)*v)/(u^4 + 2*u^2*v^2 + v^4)

    The vectors can be defined only on subsets of `S^2`, the domain of the
    result is then the common subset::

        sage: s = t(a.restrict(U), b) ; s
        Scalar field t(a,b) on the Open subset U of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()
        t(a,b): U --> R
           (x, y) |--> -(2*x - 1)*y - 3*x
        on W: (u, v) |--> -(3*u^3 + 3*u*v^2 - v^3 - (u^2 - 2*u)*v)/(u^4 + 2*u^2*v^2 + v^4)
        sage: s = t(a.restrict(U), b.restrict(W)) ; s
        Scalar field t(a,b) on the Open subset W of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()
        t(a,b): W --> R
           (x, y) |--> -(2*x - 1)*y - 3*x
           (u, v) |--> -(3*u^3 + 3*u*v^2 - v^3 - (u^2 - 2*u)*v)/(u^4 + 2*u^2*v^2 + v^4)

    The tensor itself can be defined only on some open subset of `S^2`,
    yielding a result whose domain is this subset::

        sage: s = t.restrict(V)(a,b) ; s
        Scalar field t(a,b) on the Open subset V of the 2-dimensional
         differentiable manifold S^2
        sage: s.display()
        t(a,b): V --> R
           (u, v) |--> -(3*u^3 + 3*u*v^2 - v^3 - (u^2 - 2*u)*v)/(u^4 + 2*u^2*v^2 + v^4)
        on W: (x, y) |--> -(2*x - 1)*y - 3*x

    Tests regarding the multiplication by a scalar field::

        sage: t.parent().base_ring() is f.parent()
        True
        sage: s = f*t ; s
        Tensor field of type (0,2) on the 2-dimensional differentiable
         manifold S^2
        sage: s[[0,0]] == f*t[[0,0]]
        True
        sage: s.restrict(U) == f.restrict(U)*t.restrict(U)
        True
        sage: s = f*t.restrict(U) ; s
        Tensor field of type (0,2) on the Open subset U of the 2-dimensional
         differentiable manifold S^2
        sage: s.restrict(U) == f.restrict(U)*t.restrict(U)
        True

    """
    def __init__(self, vector_field_module, tensor_type, name=None,
                 latex_name=None, sym=None, antisym=None, parent=None):
        r"""
        Construct a tensor field.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``TensorField``, to fit with the category framework::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:          intersection_name='W', restrictions1= x>0,
            ....:          restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame() ; e_uv = c_uv.frame()
            sage: XM = M.vector_field_module()
            sage: T02 = M.tensor_field_module((0,2))
            sage: t = T02.element_class(XM, (0,2), name='t'); t
            Tensor field t of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: t[e_xy,:] = [[1+x^2, x*y], [0, 1+y^2]]
            sage: t.add_comp_by_continuation(e_uv, W, c_uv)
            sage: t.display(e_xy)
            t = (x^2 + 1) dx*dx + x*y dx*dy + (y^2 + 1) dy*dy
            sage: t.display(e_uv)
            t = (3/16*u^2 + 1/16*v^2 + 1/2) du*du
             + (-1/16*u^2 + 1/4*u*v + 1/16*v^2) du*dv
             + (1/16*u^2 + 1/4*u*v - 1/16*v^2) dv*du
             + (1/16*u^2 + 3/16*v^2 + 1/2) dv*dv
            sage: TestSuite(t).run(skip='_test_pickling')

        Construction with ``DiffManifold.tensor_field``::

            sage: t1 = M.tensor_field(0, 2, name='t'); t1
            Tensor field t of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: type(t1) == type(t)
            True

        .. TODO::

            fix _test_pickling (in the superclass TensorField)

        """
        if parent is None:
            parent = vector_field_module.tensor_module(*tensor_type)
        ModuleElement.__init__(self, parent)
        self._vmodule = vector_field_module
        self._tensor_type = tuple(tensor_type)
        self._tensor_rank = self._tensor_type[0] + self._tensor_type[1]
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        self._restrictions = {} # dict. of restrictions of self on subdomains
                                # of self._domain, with the subdomains as keys
        # Treatment of symmetry declarations:
        self._sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):
                # a single symmetry is provided as a tuple -> 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self._tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                 " not in [0," + str(self._tensor_rank-1) + "]")
                    self._sym.append(tuple(isym))
        self._antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):
                # a single antisymmetry is provided as a tuple -> 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self._tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                " not in [0," + str(self._tensor_rank-1) + "]")
                    self._antisym.append(tuple(isym))
        # Final consistency check:
        index_list = []
        for isym in self._sym:
            index_list += isym
        for isym in self._antisym:
            index_list += isym
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("Incompatible lists of symmetries: the same " +
                             "position appears more than once.")
        # Initialization of derived quantities:
        self._init_derived()

    ####### Required methods for ModuleElement (beside arithmetic) #######

    def __nonzero__(self):
        r"""
        Return True if ``self`` is nonzero and False otherwise.

        This method is called by self.is_zero().

        EXAMPLE:

        Tensor field defined by parts on a 2-dimensional manifold::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1, 2, name='t')
            sage: tu = U.tensor_field(1, 2, name='t')
            sage: tv = V.tensor_field(1, 2, name='t')
            sage: tu[0,0,0] = 0
            sage: tv[0,0,0] = 0
            sage: t.set_restriction(tv)
            sage: t.set_restriction(tu)
            sage: t.__nonzero__()
            False
            sage: t.is_zero()  # indirect doctest
            True
            sage: tv[0,0,0] = 1
            sage: t.set_restriction(tv)
            sage: t.__nonzero__()
            True
            sage: t.is_zero()  # indirect doctest
            False

        """
        resu = False
        for rst in self._restrictions.itervalues():
            resu = resu or rst.__nonzero__()
        return resu

    ####### End of required methods for ModuleElement (beside arithmetic) #######

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = DiffManifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._repr_()
            'Tensor field t of type (1,3) on the 2-dimensional differentiable manifold M'
            sage: repr(t)  # indirect doctest
            'Tensor field t of type (1,3) on the 2-dimensional differentiable manifold M'
            sage: t  # indirect doctest
            Tensor field t of type (1,3) on the 2-dimensional differentiable
             manifold M

        """
        # Special cases
        if self._tensor_type == (0,2) and self._sym == [(0,1)]:
            description = "Field of symmetric bilinear forms "
            if self._name is not None:
                description += self._name + " "
        else:
        # Generic case
            description = "Tensor field "
            if self._name is not None:
                description += self._name + " "
            description += "of type ({},{}) ".format(
                                    self._tensor_type[0], self._tensor_type[1])
        return self._final_repr(description)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = DiffManifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._latex_()
            't'
            sage: t = M.tensor_field(1, 3, name='t', latex_name=r'\tau')
            sage: t._latex_()
            '\\tau'
            sage: latex(t)  # indirect doctest
            \tau

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def set_name(self, name=None, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of the tensor field.

        INPUT:

        - ``name`` -- (string; default: ``None``) name given to the tensor
          field
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the tensor field; if ``None`` while ``name`` is provided, the LaTeX
          symbol is set to ``name``.

        EXAMPLES::

            sage: M = DiffManifold(2, 'M')
            sage: t = M.tensor_field(1, 3); t
            Tensor field of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: t.set_name(name='t')
            sage: t
            Tensor field t of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: latex(t)
            t
            sage: t.set_name(latex_name=r'\tau')
            sage: latex(t)
            \tau
            sage: t.set_name(name='a')
            sage: t
            Tensor field a of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: latex(t)
            a

        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name
        for rst in self._restrictions.itervalues():
            rst.set_name(name=name, latex_name=latex_name)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same
        vector field module, with the same tensor type and same symmetries

        TESTS::

            sage: M = DiffManifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t1 = t._new_instance(); t1
            Tensor field of type (1,3) on the 2-dimensional differentiable
             manifold M
            sage: type(t1) == type(t)
            True
            sage: t1.parent() is t.parent()
            True

        """
        return self.__class__(self._vmodule, self._tensor_type, sym=self._sym,
                              antisym=self._antisym)

    def _final_repr(self, description):
        r"""
        Part of string representation common to all derived classes of
        :class:`TensorField`.

        TEST::

            sage: M = DiffManifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._final_repr('Tensor field t ')
            'Tensor field t on the 2-dimensional differentiable manifold M'

        """
        if self._domain == self._ambient_domain:
            description += "on the {}".format(self._domain)
        else:
            description += "along the {}".format(self._domain) + \
                           " with values on the {}".format(self._ambient_domain)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TEST::

            sage: M = DiffManifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._init_derived()

        """
        self._lie_derivatives = {} # dict. of Lie derivatives of self
                                   # (keys: id(vector))

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TEST::

            sage: M = DiffManifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._del_derived()

        """
        # First deletes any reference to self in the vectors' dictionaries:
        for vid, val in self._lie_derivatives.iteritems():
            del val[0]._lie_der_along_self[id(self)]
        # Then clears the dictionary of Lie derivatives
        self._lie_derivatives.clear()

    #### Simple accessors ####

    def domain(self):
        r"""
        Return the manifold on which the tensor field is defined.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.differentiable.manifold.DiffManifold`

        EXAMPLES::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: t = M.tensor_field(1,2)
            sage: t.domain()
            2-dimensional differentiable manifold M
            sage: U = M.open_subset('U', coord_def={c_xy: x<0})
            sage: h = t.restrict(U)
            sage: h.domain()
            Open subset U of the 2-dimensional differentiable manifold M

        """
        return self._domain


    def base_module(self):
        r"""
        Return the vector field module on which the current tensor field acts
        as a tensor.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`

        EXAMPLES:

        The module of vector fields on the 2-sphere as a "base module"::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'S^2') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: t = M.tensor_field(0,2)
            sage: t.base_module()
            Module X(S^2) of vector fields on the 2-dimensional differentiable
             manifold S^2
            sage: t.base_module() is M.vector_field_module()
            True
            sage: XM = M.vector_field_module()
            sage: XM.an_element().base_module() is XM
            True

        """
        return self._vmodule


    def tensor_type(self):
        r"""
        Return the tensor type of the object.

        OUTPUT:

        - pair (k,l), where k is the contravariant rank and l is the covariant
          rank

        EXAMPLES::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'S^2') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: t = M.tensor_field(1,2)
            sage: t.tensor_type()
            (1, 2)
            sage: v = M.vector_field()
            sage: v.tensor_type()
            (1, 0)

        """
        return self._tensor_type

    def tensor_rank(self):
        r"""
        Return the tensor rank of the object.

        OUTPUT:

        - integer k+l, where k is the contravariant rank and l is the covariant
          rank

        EXAMPLES::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'S^2') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: t = M.tensor_field(1,2)
            sage: t.tensor_rank()
            3
            sage: v = M.vector_field()
            sage: v.tensor_rank()
            1

        """
        return self._tensor_rank

    def symmetries(self):
        r"""
        Print the list of symmetries and antisymmetries.

        EXAMPLES::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'S^2') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: t = M.tensor_field(1,2)
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: t = M.tensor_field(1,2, sym=(1,2))
            sage: t.symmetries()
            symmetry: (1, 2);  no antisymmetry
            sage: t = M.tensor_field(2,2, sym=(0,1), antisym=(2,3))
            sage: t.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
            sage: t = M.tensor_field(2,2, antisym=[(0,1),(2,3)])
            sage: t.symmetries()
            no symmetry;  antisymmetries: [(0, 1), (2, 3)]

        """
        if len(self._sym) == 0:
            s = "no symmetry; "
        elif len(self._sym) == 1:
            s = "symmetry: " + str(self._sym[0]) + "; "
        else:
            s = "symmetries: " + str(self._sym) + "; "
        if len(self._antisym) == 0:
            a = "no antisymmetry"
        elif len(self._antisym) == 1:
            a = "antisymmetry: " + str(self._antisym[0])
        else:
            a = "antisymmetries: " + str(self._antisym)
        print s, a

    #### End of simple accessors #####


    def set_restriction(self, rst):
        r"""
        Define a restriction of the tensor field to some subdomain.

        INPUT:

        - ``rst`` -- tensor field of the same type and symmetries as ``self``,
          defined on a subdomain of ``self._domain`` (must be an instance of
          :class:`TensorField`)

        EXAMPLE::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: t = M.tensor_field(1, 2, name='t')
            sage: s = U.tensor_field(1, 2)
            sage: s[0,0,1] = x+y
            sage: t.set_restriction(s)
            sage: t.display(c_xy.frame())
            t = (x + y) d/dx*dx*dy
            sage: t.restrict(U) == s
            True

        """
        if not isinstance(rst, TensorField):
            raise TypeError("The argument must be a tensor field.")
        if not rst._domain.is_subset(self._domain):
            raise ValueError("The domain of the declared restriction is not " +
                             "a subset of the field's domain.")
        if not rst._ambient_domain.is_subset(self._ambient_domain):
            raise ValueError("The ambient domain of the declared " +
                             "restriction is not a subset of the " +
                             "field's ambient domain.")
        if rst._tensor_type != self._tensor_type:
            raise ValueError("The declared restriction has not the same " +
                             "tensor type as the current tensor field.")
        if rst._tensor_type != self._tensor_type:
            raise ValueError("The declared restriction has not the same " +
                             "tensor type as the current tensor field.")
        if rst._sym != self._sym:
            raise ValueError("The declared restriction has not the same " +
                             "symmetries as the current tensor field.")
        if rst._antisym != self._antisym:
            raise ValueError("The declared restriction has not the same " +
                             "antisymmetries as the current tensor field.")
        self._restrictions[rst._domain] = rst.copy()
        self._restrictions[rst._domain].set_name(name=self._name,
                                                 latex_name=self._latex_name)

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of the tensor field to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of ``self._domain`` (must be an
          instance of :class:`~sage.manifolds.differentiable.manifold.DiffManifold`)
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `V` is a subdomain of
          ``self._codomain``
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          If ``None``, the restriction of ``self._vmodule._dest_map`` to `U` is
          used.

        OUTPUT:

        - instance of :class:`TensorField` representing the restriction.

        EXAMPLES:

        Restrictions of a vector field on the 2-sphere::

            sage: M = DiffManifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of the North pole
            sage: stereoN.<x,y> = U.chart()  # stereographic coordinates from the North pole
            sage: eN = stereoN.frame() # the associated vector frame
            sage: V =  M.open_subset('V') # the complement of the South pole
            sage: stereoS.<u,v> = V.chart()  # stereographic coordinates from the South pole
            sage: eS = stereoS.frame() # the associated vector frame
            sage: transf = stereoN.transition_map(stereoS, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:               intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:               restrictions2= u^2+v^2!=0)
            sage: inv = transf.inverse() # transformation from stereoS to stereoN
            sage: W = U.intersection(V) # the complement of the North and South poles
            sage: stereoN_W = W.atlas()[0]  # restriction of stereographic coord. from North pole to W
            sage: stereoS_W = W.atlas()[1]  # restriction of stereographic coord. from South pole to W
            sage: eN_W = stereoN_W.frame() ; eS_W = stereoS_W.frame()
            sage: v = M.vector_field('v')
            sage: v.set_comp(eN)[1] = 1  # given the default settings, this can be abriged to v[1] = 1
            sage: v.display()
            v = d/dx
            sage: vU = v.restrict(U) ; vU
            Vector field v on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: vU.display()
            v = d/dx
            sage: vU == eN[1]
            True
            sage: vW = v.restrict(W) ; vW
            Vector field v on the Open subset W of the 2-dimensional
             differentiable manifold S^2
            sage: vW.display()
            v = d/dx
            sage: vW.display(eS_W, stereoS_W)
            v = (-u^2 + v^2) d/du - 2*u*v d/dv
            sage: vW == eN_W[1]
            True

        At this stage, defining the restriction of v to the open subset V fully
        specifies v::

            sage: v.restrict(V)[1] = vW[eS_W, 1, stereoS_W].expr()  # note that eS is the default frame on V
            sage: v.restrict(V)[2] = vW[eS_W, 2, stereoS_W].expr()
            sage: v.display(eS, stereoS)
            v = (-u^2 + v^2) d/du - 2*u*v d/dv
            sage: v.restrict(U).display()
            v = d/dx
            sage: v.restrict(V).display()
            v = (-u^2 + v^2) d/du - 2*u*v d/dv

        The restriction of the vector field to its own domain is of course
        itself::

            sage: v.restrict(M) is v
            True
            sage: vU.restrict(U) is vU
            True

        """
        if subdomain == self._domain and \
                    (dest_map is None or dest_map == self._vmodule._dest_map) :
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("The provided domain is not a subset of " +
                                 "the field's domain.")
            if dest_map is None:
                dest_map = self._vmodule._dest_map.restrict(subdomain)
            elif not dest_map._codomain.is_subset(self._ambient_domain):
                raise ValueError("Argument dest_map not compatible with " +
                                 "self._ambient_domain")
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.iteritems():
                if subdomain.is_subset(dom):
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    break
            # If this fails, the restriction is created from scratch:
            else:
                smodule = subdomain.vector_field_module(dest_map=dest_map)
                self._restrictions[subdomain] = smodule.tensor(self._tensor_type,
                                                    name=self._name,
                                                    latex_name=self._latex_name,
                                                    sym=self._sym,
                                                    antisym=self._antisym,
                                                    specific_type=self.__class__)
        return self._restrictions[subdomain]

    def set_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame
        for assignment.

        The components with respect to other frames having the same domain
        as the provided vector frame are deleted, in order to avoid any
        inconsistency. To keep them, use the method :meth:`add_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are defined; if none is provided, the components are assumed to refer
          to the tensor field domain's default frame.

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        EXAMPLE::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 2, name='t')
            sage: t.set_comp(e_uv)
            3-indices components w.r.t. Coordinate frame (V, (d/du,d/dv))
            sage: t.set_comp(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) d/dv*du*dv

        Setting the components in a new frame (``e``)::

            sage: e = V.vector_frame('e')
            sage: t.set_comp(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t.set_comp(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0*e^1*e^1

        Since the frames ``e`` and ``e_uv`` are defined on the same domain, the
        components w.r.t. ``e_uv`` have been erased::

            sage: t.display(c_uv.frame())
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (V, (d/du,d/dv))

        """
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.set_comp(basis)

    def add_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame
        for assignment.

        The components with respect to other frames having the same domain
        as the provided vector frame are kept. To delete them them, use the
        method :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are defined; if none is provided, the components are assumed to refer
          to the tensor field domain's default frame.

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.


        EXAMPLE::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 2, name='t')
            sage: t.add_comp(e_uv)
            3-indices components w.r.t. Coordinate frame (V, (d/du,d/dv))
            sage: t.add_comp(e_uv)[1,0,1] = u+v
            sage: t.display(e_uv)
            t = (u + v) d/dv*du*dv

        Setting the components in a new frame::

            sage: e = V.vector_frame('e')
            sage: t.add_comp(e)
            3-indices components w.r.t. Vector frame (V, (e_0,e_1))
            sage: t.add_comp(e)[0,1,1] = u*v
            sage: t.display(e)
            t = u*v e_0*e^1*e^1

        The components w.r.t. ``e_uv`` are kept::

            sage: t.display(e_uv)
            t = (u + v) d/dv*du*dv

        """
        if basis is None:
            basis = self._domain._def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.add_comp(basis)

    def add_comp_by_continuation(self, frame, subdomain, chart=None):
        r"""
        Set components w.r.t to a vector frame by continuation of the
        coordinate expression of the components in a subframe.

        The continuation is performed by demanding that the components have
        the same coordinate expression as those on the restriction of the frame
        to a given subdomain.

        INPUT:

        - ``frame`` -- vector frame `e` in which the components are to be set
        - ``subdomain`` -- open subset of `e`'s domain in which the
          components are known or can be evaluated from other components
        - ``chart`` -- (default: ``None``) coordinate chart on `e`'s domain in
          which the extension of the expression of the components is to be
          performed; if ``None``, the default's chart of `e`'s domain is
          assumed

        EXAMPLES:

        Components of a vector field on the sphere `S^2`::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'S^2', start_index=1)
            sage: # The two open subsets covered by stereographic coordinates (North and South):
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coordinates
            sage: transf = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), intersection_name='W', restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.vector_field('a')
            sage: a[eU,:] = [x, 2+y]

        At this stage, the vector field has been defined only on the open
        subset U (through its components in the frame eU)::

            sage: a.display(eU)
            a = x d/dx + (y + 2) d/dy

        The components with respect to the restriction of eV to the common
        subdomain W, in terms of the (u,v) coordinates, are obtained by a
        change-of-frame formula on W::

            sage: a.display(eV.restrict(W), c_uv.restrict(W))
            a = (-4*u*v - u) d/du + (2*u^2 - 2*v^2 - v) d/dv

        The continuation consists in extending the definition of the vector
        field to the whole open subset V by demanding that the components in
        the frame eV have the same coordinate expression as the above one::

            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)

        We have then::

            sage: a.display(eV)
            a = (-4*u*v - u) d/du + (2*u^2 - 2*v^2 - v) d/dv

        and `a` is defined on the entire manifold `S^2`.

        """
        dom = frame._domain
        if not dom.is_subset(self._domain):
            raise ValueError("The vector frame is not defined on a subset" +
                             " of the tensor field domain.")
        if chart is None:
            chart = dom._def_chart
        sframe = frame.restrict(subdomain)
        schart = chart.restrict(subdomain)
        scomp = self.comp(sframe)
        resu = self.add_comp(frame) # _del_derived is performed here
        for ind in resu.non_redundant_index_generator():
            resu[[ind]] = dom.scalar_field({chart: scomp[[ind]].expr(schart)})

    def comp(self, basis=None, from_basis=None):
        r"""
        Return the components in a given vector frame.

        If the components are not known already, they are computed by the
        tensor change-of-basis formula from components in another vector frame.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are required; if none is provided, the components are assumed to
          refer to the tensor field domain's default frame
        - ``from_basis`` -- (default: ``None``) vector frame from which the
          required components are computed, via the tensor change-of-basis
          formula, if they are not known already in the basis ``basis``

        OUTPUT:

        - components in the vector frame ``basis``, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES:

        Components of a type-(1,1) tensor field defined on two open subsets::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: e = U.default_frame() ; e
            Coordinate frame (U, (d/dx,d/dy))
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: f = V.default_frame() ; f
            Coordinate frame (V, (d/du,d/dv))
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e,0,0] = - x + y^3
            sage: t[e,0,1] = 2+x
            sage: t[f,1,1] = - u*v
            sage: t.comp(e)
            2-indices components w.r.t. Coordinate frame (U, (d/dx,d/dy))
            sage: t.comp(e)[:]
            [y^3 - x   x + 2]
            [      0       0]
            sage: t.comp(f)
            2-indices components w.r.t. Coordinate frame (V, (d/du,d/dv))
            sage: t.comp(f)[:]
            [   0    0]
            [   0 -u*v]

        Since e is M's default frame, the argument e can be omitted::

            sage: e is M.default_frame()
            True
            sage: t.comp() is t.comp(e)
            True

        Example of computation of the components via a change of frame::

            sage: a = V.automorphism_field()
            sage: a[:] = [[1+v, -u^2], [0, 1-u]]
            sage: h = f.new_frame(a, 'h')
            sage: t.comp(h)
            2-indices components w.r.t. Vector frame (V, (h_0,h_1))
            sage: t.comp(h)[:]
            [             0 -u^3*v/(v + 1)]
            [             0           -u*v]

        """
        if basis is None:
            basis = self._domain._def_frame
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.comp(basis=basis, from_basis=from_basis)

    def display(self, basis=None, chart=None):
        r"""
        Display the tensor field in terms of its expansion w.r.t to a given
        vector frame.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame with respect to
          which the tensor is expanded; if none is provided, the default frame
          of the domain of definition of the tensor field is assumed.
        - ``chart`` -- (default: ``None``) chart with respect to which the
          components of the tensor field in the selected frame are expressed;
          if none is provided, the default chart of the vector frame domain
          is assumed.

        EXAMPLES:

        Display of a type-(1,1) tensor field defined on two open subsets::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: e = U.default_frame() ; e
            Coordinate frame (U, (d/dx,d/dy))
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: f = V.default_frame() ; f
            Coordinate frame (V, (d/du,d/dv))
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e,0,0] = - x + y^3
            sage: t[e,0,1] = 2+x
            sage: t[f,1,1] = - u*v
            sage: t.display(e)
            t = (y^3 - x) d/dx*dx + (x + 2) d/dx*dy
            sage: t.display(f)
            t = -u*v d/dv*dv

        Since e is M's default frame, the argument e can be omitted::

            sage: e is M.default_frame()
            True
            sage: t.display()
            t = (y^3 - x) d/dx*dx + (x + 2) d/dx*dy

        Similarly, since f is V's default frame, the argument f can be omitted
        when considering the restriction of t to V::

            sage: t.restrict(V).display()
            t = -u*v d/dv*dv

        Display w.r.t a frame in which t has not been initialized (automatic
        use of a change-of-frame formula)::

            sage: a = V.automorphism_field()
            sage: a[:] = [[1+v, -u^2], [0, 1-u]]
            sage: h = f.new_frame(a, 'h')
            sage: t.display(h)
            t = -u^3*v/(v + 1) h_0*h^1 - u*v h_1*h^1

        A shortcut of ``display()`` is ``disp()``::

            sage: t.disp(h)
            t = -u^3*v/(v + 1) h_0*h^1 - u*v h_1*h^1

        """
        if basis is None:
            if self._vmodule._dest_map.is_identity():
                basis = self._domain._def_frame
            else:
                for rst in self._restrictions.values():
                    try:
                        return rst.display()
                    except ValueError:
                        pass
            if basis is None:  # should be "is still None" ;-)
                raise ValueError("a frame must be provided for the display")
        rst = self.restrict(basis._domain, dest_map=basis._dest_map)
        return rst.display(basis, chart)

    disp = display

    def display_comp(self, frame=None, chart=None, coordinate_labels=True,
                     only_nonzero=True, only_nonredundant=False):
        r"""
        Display the tensor components w.r.t. a given frame, one per
        line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to which
          the tensor field components are defined; if ``None``, then

          - if ``chart`` is not ``None``, the coordinate frame associated to
            ``chart`` is used
          - otherwise, the default basis of the vector field module on which
            the tensor field is defined is used

        - ``chart`` -- (default: ``None``) chart specifying the coordinate
          expression of the components; if ``None``, the default chart of the
          tensor field domain is used
        - ``coordinate_labels`` -- (default: ``True``) boolean; if ``True``,
          coordinate symbols are used by default (instead of integers) as
          index labels whenever ``frame`` is a coordinate frame
        - ``only_nonzero`` -- (default: ``True``) boolean; if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- (default: ``False``) boolean; if ``True``,
          only nonredundant components are displayed in case of symmetries

        EXAMPLES:

        Display of the components of a type-(1,1) tensor field defined on two
        open subsets::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x, y> = U.chart()
            sage: e = U.default_frame()
            sage: V = M.open_subset('V')
            sage: c_uv.<u, v> = V.chart()
            sage: f = V.default_frame()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e,0,0] = - x + y^3
            sage: t[e,0,1] = 2+x
            sage: t[f,1,1] = - u*v
            sage: t.display_comp(e)
            t^x_x = y^3 - x
            t^x_y = x + 2
            sage: t.display_comp(f)
            t^v_v = -u*v

        Components in a chart frame::

            sage: t.display_comp(chart=c_xy)
            t^x_x = y^3 - x
            t^x_y = x + 2
            sage: t.display_comp(chart=c_uv)
            t^v_v = -u*v

        See documentation of :meth:`TensorFieldParal.display_comp` for more
        options.

        """
        if frame is None:
            if chart is not None:
                frame = chart.frame()
            else:
                if self._vmodule._dest_map.is_identity():
                    frame = self._domain.default_frame()
                else:
                    for rst in self._restrictions.values():
                        try:
                            return rst.display_comp(chart=chart,
                                       coordinate_labels=coordinate_labels,
                                       only_nonzero=only_nonzero,
                                       only_nonredundant=only_nonredundant)
                        except ValueError:
                            pass
                if frame is None:  # should be "is still None" ;-)
                    raise ValueError("a frame must be provided for the display")
        rst = self.restrict(frame.domain(), dest_map=frame._dest_map)
        return rst.display_comp(frame=frame, chart=chart,
                                coordinate_labels=coordinate_labels,
                                only_nonzero=only_nonzero,
                                only_nonredundant=only_nonredundant)


    def __getitem__(self, args):
        r"""
        Return a component w.r.t. some frame.

        NB: if ``args`` is a string, this method acts as a shortcut for
        tensor contractions and symmetrizations, the string containing
        abstract indices.

        INPUT:

        - ``args`` -- list of indices defining the component; if [:] is
          provided, all the components are returned. The frame can be passed
          as the first item of ``args``; if not, the default frame of the
          tensor field's domain is assumed.

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_xy = c_xy.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy, :] = [[x+y, -2], [3*y^2, x*y]]
            sage: t.__getitem__((1,0))
            3*y^2
            sage: t.__getitem__((1,1))
            x*y
            sage: t.__getitem__((e_xy,1,0))
            3*y^2
            sage: t.__getitem__(slice(None))
            [x + y    -2]
            [3*y^2   x*y]
            sage: t.__getitem__((e_xy,slice(None)))
            [x + y    -2]
            [3*y^2   x*y]
            sage: t.__getitem__('^a_a')  # trace
            Scalar field on the 2-dimensional differentiable manifold M
            sage: t.__getitem__('^a_a').display()
            M --> R
            on U: (x, y) |--> (x + 1)*y + x

        """
        if isinstance(args, str): # tensor with specified indices
            return TensorWithIndices(self, args).update()
        if isinstance(args, list):  # case of [[...]] syntax
            if not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._domain._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        return self.comp(frame)[args]

    def __setitem__(self, args, value):
        r"""
        Sets a component w.r.t to some vector frame.

        INPUT:

       - ``args`` -- list of indices; if [:] is provided, all the components
          are set. The frame can be passed as the first item of ``args``; if
          not, the default frame of the tensor field's domain is assumed.
        - ``value`` -- the value to be set or a list of values if ``args``
          == ``[:]``

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: e_xy = c_xy.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t.__setitem__((e_xy, 0, 1), x+y^2)
            sage: t.display(e_xy)
            t = (y^2 + x) d/dx*dy
            sage: t.__setitem__((0, 1), x+y^2)  # same as above since e_xy is the default frame on M
            sage: t.display()
            t = (y^2 + x) d/dx*dy
            sage: t.__setitem__(slice(None), [[x+y, -2], [3*y^2, x*y]])
            sage: t.display()
            t = (x + y) d/dx*dx - 2 d/dx*dy + 3*y^2 d/dy*dx + x*y d/dy*dy

        """
        if isinstance(args, list):  # case of [[...]] syntax
            if not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._domain._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._domain._def_frame
        self.set_comp(frame)[args] = value


    def copy(self):
        r"""
        Return an exact copy of ``self``.

        The name and the derived quantities are not copied.

        EXAMPLE:

        Copy of a type-(1,1) tensor field on the 2-sphere::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = t.copy(); s
            Tensor field of type (1,1) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s == t
            True

        If the original tensor field is modified, the copy is not:

            sage: t[e_xy,0,0] = -1
            sage: t.display(e_xy)
            t = -d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s.display(e_xy)
            (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s == t
            False

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.iteritems():
            resu._restrictions[dom] = rst.copy()
        return resu

    def _common_subdomains(self, other):
        r"""
        Return the list of subdomains of self._domain on which both ``self``
        and ``other`` have known restrictions.

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: t._common_subdomains(t)  # random
            [Open subset V of the 2-dimensional differentiable manifold M,
             Open subset U of the 2-dimensional differentiable manifold M,
             Open subset W of the 2-dimensional differentiable manifold M]
            sage: a = M.tensor_field(1, 1, name='a')
            sage: t._common_subdomains(a)
            []
            sage: a[e_xy, 0, 1] = 0
            sage: t._common_subdomains(a)
            [Open subset U of the 2-dimensional differentiable manifold M]
            sage: a[e_uv, 0, 0] = 0
            sage: t._common_subdomains(a)  # random
            [Open subset V of the 2-dimensional differentiable manifold M,
             Open subset U of the 2-dimensional differentiable manifold M]

        """
        resu = []
        for dom in self._restrictions:
            if dom in other._restrictions:
                resu.append(dom)
        return resu

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a tensor field or 0

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: t.__eq__(t)
            True
            sage: t.__eq__(t.copy())
            True
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a.set_restriction(t.restrict(U))
            sage: a.set_restriction(t.restrict(V))
            sage: t.__eq__(a)
            True
            sage: a[e_xy, 0, 0] = -1
            sage: t.__eq__(a)
            False
            sage: t.parent().zero().__eq__(0)
            True

        """
        if isinstance(other, (int, Integer)): # other should be 0
            if other == 0:
                return self.is_zero()
            else:
                return False
        elif not isinstance(other, TensorField):
            return False
        else: # other is another tensor field
            if other._vmodule != self._vmodule:
                return False
            if other._tensor_type != self._tensor_type:
                return False
            resu = True
            for dom, rst in self._restrictions.iteritems():
                if dom in other._restrictions:
                    resu = resu and bool(rst == other._restrictions[dom])
            return resu

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- a tensor field or 0

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        TEST::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.__ne__(t)
            False
            sage: t.__ne__(t.copy())
            False

        """
        return not self.__eq__(other)

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = t.__pos__(); s
            Tensor field +t of type (1,1) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            +t = (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.iteritems():
            resu._restrictions[dom] = + rst
        if self._name is not None:
            resu._name = '+' + self._name
        if self._latex_name is not None:
            resu._latex_name = '+' + self._latex_name
        return resu

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the tensor field `-T`, where `T` is ``self``

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy, :] = [[x, -x], [y, -y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: t.display(e_xy)
            t = x d/dx*dx - x d/dx*dy + y d/dy*dx - y d/dy*dy
            sage: t.display(e_uv)
            t = (u^3 - 2*u^2*v - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*du
             + (u^3 + 2*u^2*v - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*dv
             + (u^2*v - 2*u*v^2 - v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*du
             + (u^2*v + 2*u*v^2 - v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*dv
            sage: s = t.__neg__(); s
            Tensor field -t of type (1,1) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            -t = -x d/dx*dx + x d/dx*dy - y d/dy*dx + y d/dy*dy
            sage: s.display(e_uv)
            -t = -(u^3 - 2*u^2*v - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*du
             - (u^3 + 2*u^2*v - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*dv
             - (u^2*v - 2*u*v^2 - v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*du
             - (u^2*v + 2*u*v^2 - v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*dv
            sage: s == -t  # indirect doctest
            True

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.iteritems():
            resu._restrictions[dom] = - rst
        if self._name is not None:
            resu._name = '-' + self._name
        if self._latex_name is not None:
            resu._latex_name = '-' + self._latex_name
        return resu

    ######### ModuleElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Tensor field addition.

        INPUT:

        - ``other`` -- a tensor field, in the same tensor module as ``self``

        OUPUT:

        - the tensor field resulting from the addition of ``self`` and ``other``

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.tensor_field(1, 1, name='b')
            sage: b[e_xy,:] = [[2, y], [x, -x]]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a._add_(b); s
            Tensor field a+b of type (1,1) on the 2-dimensional differentiable manifold M
            sage: a.display(e_xy)
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: b.display(e_xy)
            b = 2 d/dx*dx + y d/dx*dy + x d/dy*dx - x d/dy*dy
            sage: s.display(e_xy)
            a+b = (x + 2) d/dx*dx + (y + 1) d/dx*dy + (x + y) d/dy*dx - x d/dy*dy
            sage: s == a + b  # indirect doctest
            True
            sage: z = a.parent().zero(); z
            Tensor field zero of type (1,1) on the 2-dimensional differentiable manifold M
            sage: a._add_(z) == a
            True
            sage: z._add_(a) == a
            True

        """
        if other == 0:
            return +self
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] + other._restrictions[dom]
        some_rst = resu_rst.values()[0]
        resu_sym = some_rst._sym
        resu_antisym = some_rst._antisym
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_sym,
                                    antisym=resu_antisym)
        resu._restrictions = resu_rst
        if self._name is not None and other._name is not None:
            resu._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            resu._latex_name = self._latex_name + '+' + other._latex_name
        return resu

    def _sub_(self, other):
        r"""
        Tensor field subtraction.

        INPUT:

        - ``other`` -- a tensor field, in the same tensor module as ``self``

        OUPUT:

        - the tensor field resulting from the subtraction of ``other`` from
          ``self``

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.tensor_field(1, 1, name='b')
            sage: b[e_xy,:] = [[2, y], [x, -x]]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a._sub_(b); s
            Tensor field a-b of type (1,1) on the 2-dimensional differentiable manifold M
            sage: a.display(e_xy)
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: b.display(e_xy)
            b = 2 d/dx*dx + y d/dx*dy + x d/dy*dx - x d/dy*dy
            sage: s.display(e_xy)
            a-b = (x - 2) d/dx*dx + (-y + 1) d/dx*dy + (-x + y) d/dy*dx + x d/dy*dy
            sage: s == a - b
            True
            sage: z = a.parent().zero()
            sage: a._sub_(z) == a
            True
            sage: z._sub_(a) == -a
            True

        """
        if other == 0:
            return +self
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] - other._restrictions[dom]
        some_rst = resu_rst.values()[0]
        resu_sym = some_rst._sym
        resu_antisym = some_rst._antisym
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_sym,
                                   antisym=resu_antisym)
        resu._restrictions = resu_rst
        if self._name is not None and other._name is not None:
            resu._name = self._name + '-' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            resu._latex_name = self._latex_name + '-' + other._latex_name
        return resu

    def _rmul_(self, scalar):
        r"""
        Reflected multiplication operator: performs ``scalar * self``

        This is actually the multiplication by an element of the ring over
        which the tensor field module is constructed.

        INPUT:

        - ``scalar`` -- a scalar field, in the scalar field algebra over which
          the module containing ``self`` is defined

        OUPUT:

        - the tensor field ``scalar * self``

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2)}, name='f')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: f.display()
            f: M --> R
            on U: (x, y) |--> 1/(x^2 + y^2 + 1)
            on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)
            sage: s = a._rmul_(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable manifold M
            sage: a.display(e_xy)
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: s.display(e_xy)
            x/(x^2 + y^2 + 1) d/dx*dx + 1/(x^2 + y^2 + 1) d/dx*dy
             + y/(x^2 + y^2 + 1) d/dy*dx
            sage: a.display(e_uv)
            a = (2*u^3*v - 2*u*v^3 + u^3 - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*du
             - (u^4 - 2*u^2*v^2 + v^4 - 2*u^2*v)/(u^4 + 2*u^2*v^2 + v^4) d/du*dv
             + (4*u^2*v^2 + u^2*v - v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*du
             - 2*(u^3*v - u*v^3 - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/dv*dv
            sage: s.display(e_uv)
            (2*u^3*v - 2*u*v^3 + u^3 - u*v^2)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/du*du
             - (u^4 - 2*u^2*v^2 + v^4 - 2*u^2*v)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/du*dv
             + (4*u^2*v^2 + u^2*v - v^3)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/dv*du
             - 2*(u^3*v - u*v^3 - u*v^2)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/dv*dv
            sage: s == f*a  # indirect doctest
            True
            sage: z = a.parent().zero(); z
            Tensor field zero of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a._rmul_(M.zero_scalar_field()) == z
            True
            sage: z._rmul_(f) == z
            True

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.iteritems():
            resu._restrictions[dom] = scalar.restrict(dom) * rst
        return resu

    ######### End of ModuleElement arithmetic operators ########

    def __mul__(self, other):
        r"""
        Tensor product (or multiplication of the right by a scalar).

        INPUT:

        - ``other`` -- a tensor field, on the same manifold as ``self`` (or an
          object that can be coerced to a scalar field on the same manifold
          as ``self``)

        OUPUT:

        - the tensor field resulting from the tensor product of ``self`` with
          ``other`` (or from the product ``other * self`` if ``other`` is a
          scalar)

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)

        Tensor product with another tensor field::

            sage: b = M.vector_field(name='b')
            sage: b[e_xy, :] = [x, y]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a.__mul__(b); s
            Tensor field of type (2,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            a*b = x^2 d/dx*d/dx*dx + x d/dx*d/dx*dy + x*y d/dx*d/dy*dx
             + y d/dx*d/dy*dy + x*y d/dy*d/dx*dx + y^2 d/dy*d/dy*dx
            sage: s.display(e_uv)
            a*b = -(2*u^4*v - 2*u^2*v^3 + u^4 - u^2*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*d/du*du
             + (u^5 - 2*u^3*v^2 + u*v^4 - 2*u^3*v)/(u^4 + 2*u^2*v^2 + v^4) d/du*d/du*dv
             - (2*u^3*v^2 - 2*u*v^4 + u^3*v - u*v^3)/(u^4 + 2*u^2*v^2 + v^4) d/du*d/dv*du
             + (u^4*v - 2*u^2*v^3 + v^5 - 2*u^2*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*d/dv*dv
             - (4*u^3*v^2 + u^3*v - u*v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*d/du*du
             + 2*(u^4*v - u^2*v^3 - u^2*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/dv*d/du*dv
             - (4*u^2*v^3 + u^2*v^2 - v^4)/(u^4 + 2*u^2*v^2 + v^4) d/dv*d/dv*du
             + 2*(u^3*v^2 - u*v^4 - u*v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*d/dv*dv

        Multiplication on the right by a scalar field::

            sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2)}, name='f')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: f.display()
            f: M --> R
            on U: (x, y) |--> 1/(x^2 + y^2 + 1)
            on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)
            sage: s = a.__mul__(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            x/(x^2 + y^2 + 1) d/dx*dx + 1/(x^2 + y^2 + 1) d/dx*dy
             + y/(x^2 + y^2 + 1) d/dy*dx
            sage: s.display(e_uv)
            (2*u^3*v - 2*u*v^3 + u^3 - u*v^2)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/du*du
             - (u^4 - 2*u^2*v^2 + v^4 - 2*u^2*v)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/du*dv
             + (4*u^2*v^2 + u^2*v - v^3)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/dv*du
             - 2*(u^3*v - u*v^3 - u*v^2)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2) d/dv*dv
            sage: s == f*a
            True

        Multiplication on the right by a number::

            sage: s = a.__mul__(2); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: s.display(e_xy)
            2*x d/dx*dx + 2 d/dx*dy + 2*y d/dy*dx
            sage: s.display(e_uv)
            2*(2*u^3*v - 2*u*v^3 + u^3 - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*du
             - 2*(u^4 - 2*u^2*v^2 + v^4 - 2*u^2*v)/(u^4 + 2*u^2*v^2 + v^4) d/du*dv
             + 2*(4*u^2*v^2 + u^2*v - v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*du
             - 4*(u^3*v - u*v^3 - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/dv*dv
            sage: s.restrict(U) == 2*a.restrict(U)
            True
            sage: s.restrict(V) == 2*a.restrict(V)
            True
            sage: s == 2*a
            True

        """
        if not isinstance(other, TensorField):
            # Multiplication by a scalar field or a number
            return other * self
        # Tensor product:
        dom_resu = self._domain.intersection(other._domain)
        ambient_dom_resu = self._ambient_domain.intersection(
                                                         other._ambient_domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the FreeModuleTensor version:
            return FreeModuleTensor.__mul__(self_r, other_r)
        dest_map = self._vmodule._dest_map
        dest_map_resu = dest_map.restrict(dom_resu,
                                          subcodomain=ambient_dom_resu)
        vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)
        com_dom = []
        for dom in self_r._restrictions:
            if dom in other_r._restrictions:
                com_dom.append(dom)
        resu_rst = []
        for dom in com_dom:
            self_rr = self_r._restrictions[dom]
            other_rr = other_r._restrictions[dom]
            resu_rst.append(self_rr * other_rr)
        k1, l1 = self._tensor_type
        k2, l2 = other._tensor_type
        resu = vmodule.tensor((k1+k2, l1+l2), sym=resu_rst[0]._sym,
                              antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def __div__(self, scalar):
        r"""
        Division by a scalar field.

        INPUT:

        - ``scalar`` -- a scalar field, in the scalar field algebra over which
          the module containing ``self`` is defined

        OUPUT:

        - the tensor field ``scalar * self``

        TESTS::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: f = M.scalar_field({c_xy: (x^2 + y^2 + 2)/(x^2 + y^2 + 1)}, name='f')
            sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
            sage: f.display()
            f: M --> R
            on U: (x, y) |--> (x^2 + y^2 + 2)/(x^2 + y^2 + 1)
            on V: (u, v) |--> (2*u^2 + 2*v^2 + 1)/(u^2 + v^2 + 1)
            sage: s = a.__div__(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            (x^3 + x*y^2 + x)/(x^2 + y^2 + 2) d/dx*dx
             + (x^2 + y^2 + 1)/(x^2 + y^2 + 2) d/dx*dy
             + (y^3 + (x^2 + 1)*y)/(x^2 + y^2 + 2) d/dy*dx
            sage: f*s == a
            True

        Division by a number::

            sage: s = a.__div__(2); s
            Tensor field of type (1,1) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            1/2*x d/dx*dx + 1/2 d/dx*dy + 1/2*y d/dy*dx
            sage: s.display(e_uv)
            1/2*(2*u^3*v - 2*u*v^3 + u^3 - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du*du - 1/2*(u^4 - 2*u^2*v^2 + v^4 - 2*u^2*v)/(u^4 + 2*u^2*v^2 + v^4) d/du*dv + 1/2*(4*u^2*v^2 + u^2*v - v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv*du - (u^3*v - u*v^3 - u*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/dv*dv
            sage: s == a/2
            True
            sage: 2*s == a
            True

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.iteritems():
            resu._restrictions[dom] = rst / scalar
        return resu

    def __call__(self, *args):
        r"""
        The tensor field acting on 1-forms and vector fields as a multilinear
        map.

        In the particular case of tensor field of type (1,1), the action can
        be on a single vector field, the tensor field being identified to a
        field of tangent-space endomorphisms. The output is then a vector
        field.

        INPUT:

        - ``*args`` -- list of k 1-forms and l vector fields, ``self``
          being a tensor of type (k,l).

        OUTPUT:

        - either the scalar field resulting from the action of ``self`` on
          the 1-forms and vector fields passed as arguments or the vector
          field resulting from the action of ``self`` as a field of
          tangent-space endomorphisms (case of a type-(1,1) tensor field)

        TESTS:

        Action of a tensor field of type (1,1) on the 2-sphere::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e_xy,:] = [[x, 1], [y, 0]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: w = M.vector_field(name='w')
            sage: w[e_xy,:] = [y^2, x^2]
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: a = M.one_form(name='a')
            sage: a[e_xy,:] = [-1+y, x*y]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)

        The tensor field acting on a pair (1-form, vector field)::

            sage: s = t.__call__(a,w); s
            Scalar field t(a,w) on the 2-dimensional differentiable manifold M
            sage: s.display()
            t(a,w): M --> R
            on U: (x, y) |--> x*y^4 + x*y^3 + x^2*y - x*y^2 - x^2
            on V: (u, v) |--> -(u^8 - u^6*v + (u^2 + u)*v^6 - (u^2 + u)*v^5 + (3*u^4 + 2*u^3 - u)*v^4 - (2*u^4 + u^3)*v^3 + (3*u^6 + u^5)*v^2)/(u^10 + 5*u^8*v^2 + 10*u^6*v^4 + 10*u^4*v^6 + 5*u^2*v^8 + v^10)
            sage: s.restrict(U) == t.restrict(U)(a.restrict(U), w.restrict(U))
            True
            sage: s.restrict(V) == t.restrict(V)(a.restrict(V), w.restrict(V))
            True

        The tensor field acting on vector field, as a field of tangent-space
        endomorphisms::

            sage: s = t.__call__(w); s
            Vector field t(w) on the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            t(w) = (x*y^2 + x^2) d/dx + y^3 d/dy
            sage: s.display(e_uv)
            t(w) = -(u^4 - (u^2 - u)*v^2)/(u^4 + 2*u^2*v^2 + v^4) d/du
             - (2*u^3*v + v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv
            sage: s.restrict(U) == t.restrict(U)(w.restrict(U))
            True
            sage: s.restrict(V) == t.restrict(V)(w.restrict(V))
            True

        """
        p = len(args)
        if p == 1 and self._tensor_type == (1,1):
            # type-(1,1) tensor acting as a a field of tangent-space
            # endomorphisms:
            vector = args[0]
            if vector._tensor_type != (1,0):
                raise TypeError("The argument must be a vector field.")
            dom_resu = self._domain.intersection(vector._domain)
            if dom_resu.is_manifestly_parallelizable():
                # call of the TensorFieldParal version:
                return self.restrict(dom_resu)(vector.restrict(dom_resu))
            if self._name is not None and vector._name is not None:
                name_resu = self._name + "(" + vector._name + ")"
            else:
                name_resu = None
            if self._latex_name is not None and vector._latex_name is not None:
                latex_name_resu = self._latex_name + r"\left(" + \
                                  vector._latex_name + r"\right)"
            else:
                latex_name_resu = None
            dest_map = vector._vmodule._dest_map
            dest_map_resu = dest_map.restrict(dom_resu)
            resu = dom_resu.vector_field(name=name_resu,
                                         latex_name=latex_name_resu,
                                         dest_map=dest_map_resu)
            for dom in self._common_subdomains(vector):
                if dom.is_subset(dom_resu):
                    resu._restrictions[dom] = \
                        self._restrictions[dom](vector._restrictions[dom])
            return resu
        # Generic case
        if p != self._tensor_rank:
            raise TypeError(str(self._tensor_rank) +
                            " arguments must be provided.")
        # Domain of the result
        dom_resu = self._domain
        ambient_dom = self._ambient_domain
        for arg in args:
            dom_resu = dom_resu.intersection(arg._domain)
            ambient_dom = ambient_dom.intersection(arg._ambient_domain)
        self_r = self.restrict(dom_resu)
        args_r = [args[i].restrict(dom_resu) for i in range(p)]
        if ambient_dom.is_manifestly_parallelizable():
            # TensorFieldParal version
            return self_r(*args_r)
        else:
            resu = dom_resu.scalar_field()
            com_dom = []
            for dom in self_r._restrictions:
                for arg in args_r:
                    if dom not in arg._restrictions:
                        break
                else:
                    com_dom.append(dom)
            for dom in com_dom:
                self_rr = self_r._restrictions[dom]
                args_rr = [args_r[i]._restrictions[dom] for i in range(p)]
                resu_rr = self_rr(*args_rr)
                if resu_rr == 0:
                    for chart in resu_rr._domain._atlas:
                        resu._express[chart] = chart._zero_function
                else:
                    for chart, expr in resu_rr._express.iteritems():
                        resu._express[chart] = expr
            if resu == 0:
                return dom_resu._zero_scalar_field
            # Name of the output:
            res_name = None
            if self._name is not None:
                res_name = self._name + "("
                for i in range(p-1):
                    if args[i]._name is not None:
                        res_name += args[i]._name + ","
                    else:
                        res_name = None
                        break
                if res_name is not None:
                    if args[p-1]._name is not None:
                        res_name += args[p-1]._name + ")"
                    else:
                        res_name = None
            resu._name = res_name
            # LaTeX symbol of the output:
            res_latex = None
            if self._latex_name is not None:
                res_latex = self._latex_name + r"\left("
                for i in range(p-1):
                    if args[i]._latex_name is not None:
                        res_latex += args[i]._latex_name + ","
                    else:
                        res_latex = None
                        break
                if res_latex is not None:
                    if args[p-1]._latex_name is not None:
                        res_latex += args[p-1]._latex_name + r"\right)"
                    else:
                        res_latex = None
            resu._latex_name = res_latex
            return resu

    def trace(self, pos1=0, pos2=1):
        r"""
        Trace (contraction) on two slots of the tensor field.

        INPUT:

        - ``pos1`` -- (default: 0) position of the first index for the
          contraction, with the convention ``pos1=0`` for the first slot
        - ``pos2`` -- (default: 1) position of the second index for the
          contraction, with the same convention as for ``pos1``. The variance
          type of ``pos2`` must be opposite to that of ``pos1``

        OUTPUT:

        - tensor field resulting from the (pos1, pos2) contraction

        EXAMPLES:

        Trace of a type-(1,1) tensor field on a 2-dimensional
        non-parallelizable manifold::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(1,1, name='a')
            sage: a[eU,:] = [[1,x], [2,y]]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = a.trace() ; s
            Scalar field on the 2-dimensional differentiable manifold M
            sage: s.display()
            M --> R
            on U: (x, y) |--> y + 1
            on V: (u, v) |--> 1/2*u - 1/2*v + 1
            sage: s == a.trace(0,1) # explicit mention of the positions
            True

        Instead of the explicit call to the method :meth:`trace`, one
        may use the index notation with Einstein convention (summation over
        repeated indices); it suffices to pass the indices as a string inside
        square brackets::

            sage: a['^i_i']
            Scalar field on the 2-dimensional differentiable manifold M
            sage: a['^i_i'] == s
            True

        Any letter can be used to denote the repeated index::

            sage: a['^b_b'] == s
            True

        Trace of a type-(1,2) tensor field::

            sage: b = M.tensor_field(1,2, name='b') ; b
            Tensor field b of type (1,2) on the 2-dimensional differentiable
             manifold M
            sage: b[eU,:] = [[[1,x], [2,y]], [[0,y-3], [x*y,y^2]]]
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = b.trace(0,1) ; s # contraction on first and second slots
            1-form on the 2-dimensional differentiable manifold M
            sage: s.display(eU)
            (x*y + 1) dx + (y^2 + x) dy
            sage: s.display(eV)
            (1/4*u^2 - 1/4*(u - 1)*v + 1/4*u + 1/2) du
             + (1/4*(u - 1)*v - 1/4*v^2 - 1/4*u + 1/2) dv

        Use of the index notation::

            sage: b['^k_ki']
            1-form on the 2-dimensional differentiable manifold M
            sage: b['^k_ki'] == s
            True

        Indices not involved in the contraction may be replaced by dots::

            sage: b['^k_k.'] == s
            True

        The symbol '^' may be omitted::

            sage: b['k_k.'] == s
            True

        LaTeX notations are allowed::

            sage: b['^{k}_{ki}'] == s
            True

        Contraction on first and third slots::

            sage: s = b.trace(0,2) ; s
            1-form on the 2-dimensional differentiable manifold M
            sage: s.display(eU)
            (y - 2) dx + (y^2 + 2) dy
            sage: s.display(eV)
            (1/8*u^2 - 1/4*(u + 1)*v + 1/8*v^2 + 1/4*u) du
             + (-1/8*u^2 + 1/4*(u - 1)*v - 1/8*v^2 + 1/4*u - 2) dv

        Use of index notation::

            sage: b['^k_.k'] == s
            True

        """
        # The indices at pos1 and pos2 must be of different types:
        k_con = self._tensor_type[0]
        l_cov = self._tensor_type[1]
        if pos1 < k_con and pos2 < k_con:
            raise IndexError("Contraction on two contravariant indices is " +
                             "not allowed.")
        if pos1 >= k_con and pos2 >= k_con:
            raise IndexError("Contraction on two covariant indices is " +
                             "not allowed.")
        resu_rst = []
        for rst in self._restrictions.itervalues():
            resu_rst.append(rst.trace(pos1, pos2))
        if (k_con, l_cov) == (1,1):
            # scalar field result
            resu = self._domain.scalar_field()
            all_zero = True
            for rst in resu_rst:
                if rst == 0:
                    for chart in rst._domain._atlas:
                        resu._express[chart] = 0
                else:
                    all_zero = False
                    for chart, funct in rst._express.iteritems():
                        resu._express[chart] = funct
            if all_zero:
                resu = self._domain._zero_scalar_field
        else:
            # tensor field result
            resu = self._vmodule.tensor((k_con-1, l_cov-1),
                            sym=resu_rst[0]._sym, antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def contract(self, *args):
        r"""
        Contraction with another tensor field, on one or more indices.

        INPUT:

        - ``pos1`` -- positions of the indices in ``self`` involved in the
          contraction; ``pos1`` must be a sequence of integers, with 0 standing
          for the first index position, 1 for the second one, etc. If ``pos1``
          is not provided, a single contraction on the last index position of
          ``self`` is assumed
        - ``other`` -- the tensor field to contract with
        - ``pos2`` -- positions of the indices in ``other`` involved in the
          contraction, with the same conventions as for ``pos1``. If ``pos2``
          is not provided, a single contraction on the first index position of
          ``other`` is assumed

        OUTPUT:

        - tensor field resulting from the contraction at the positions
          ``pos1`` and ``pos2`` of ``self`` with ``other``

        EXAMPLES:

        Contractions of a type-(1,1) tensor field with a type-(2,0) one on
        a 2-dimensional non-parallelizable manifold::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(1,1, name='a')
            sage: a[eU,:] = [[1,x], [0,2]]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: b = M.tensor_field(2,0, name='b')
            sage: b[eU,:] = [[y,-1], [x+y,2]]
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = a.contract(b) ; s   # contraction on last index of a and first one of b
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M

        Check 1: components w.r.t. the manifold's default frame (eU)::

            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         print bool(s[i,j] == sum(a[i,k]*b[k,j] for k in M.irange())),
            ....:
            True True True True

        Check 2: components w.r.t. frame eV::

            sage: for i in M.irange():
            ....:     for j in M.irange():
            ....:         print bool(s[eV,i,j] == sum(a[eV,i,k]*b[eV,k,j] for k in M.irange())),
            ....:
            True True True True

        Instead of the explicit call to the method :meth:`contract`, one
        may use the index notation with Einstein convention (summation over
        repeated indices); it suffices to pass the indices as a string inside
        square brackets::

            sage: a['^i_k']*b['^kj'] == s
            True

        Indices not involved in the contraction may be replaced by dots::

            sage: a['^._k']*b['^k.'] == s
            True

        LaTeX notation may be used::

            sage: a['^{i}_{k}']*b['^{kj}'] == s
            True

        Contraction on the last index of a and last index of b::

            sage: s = a.contract(b, 1) ; s
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: a['^i_k']*b['^jk'] == s
            True

        Contraction on the first index of b and the last index of a::

            sage: s = b.contract(0,a,1) ; s
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: b['^ki']*a['^j_k'] == s
            True

        The domain of the result is the intersection of the domains of the two
        tensor fields::

            sage: aU = a.restrict(U) ; bV = b.restrict(V)
            sage: s = aU.contract(b) ; s
            Tensor field of type (2,0) on the Open subset U of the
             2-dimensional differentiable manifold M
            sage: s = a.contract(bV) ; s
            Tensor field of type (2,0) on the Open subset V of the
             2-dimensional differentiable manifold M
            sage: s = aU.contract(bV) ; s
            Tensor field of type (2,0) on the Open subset W of the
             2-dimensional differentiable manifold M
            sage: s0 = a.contract(b)
            sage: s == s0.restrict(W)
            True

        The contraction can be performed on more than one index: c being a
        type-(2,2) tensor, contracting the indices in positions 2 and 3 of c
        with respectively those in positions 0 and 1 of b is::

            sage: c = a*a ; c
            Tensor field of type (2,2) on the 2-dimensional differentiable
             manifold M
            sage: s = c.contract(2,3, b, 0,1) ; s
            Tensor field of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: s == c['^.._kl']*b['^kl']  # the same double contraction in index notation
            True

        The symmetries are either conserved or destroyed by the contraction::

            sage: c = c.symmetrize(0,1).antisymmetrize(2,3)
            sage: c.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
            sage: s = b.contract(0, c, 2) ; s
            Tensor field of type (3,1) on the 2-dimensional differentiable
             manifold M
            sage: s.symmetries()
            symmetry: (1, 2);  no antisymmetry

        Case of a scalar field result::

            sage: a = M.one_form('a')
            sage: a[eU,:] = [y, 1+x]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: b = M.vector_field('b')
            sage: b[eU,:] = [x, y^2]
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a.display(eU)
            a = y dx + (x + 1) dy
            sage: b.display(eU)
            b = x d/dx + y^2 d/dy
            sage: s = a.contract(b) ; s
            Scalar field on the 2-dimensional differentiable manifold M
            sage: s.display()
            M --> R
            on U: (x, y) |--> (x + 1)*y^2 + x*y
            on V: (u, v) |--> 1/8*u^3 - 1/8*u*v^2 + 1/8*v^3 + 1/2*u^2 - 1/8*(u^2 + 4*u)*v
            sage: s == a['_i']*b['^i'] # use of index notation
            True
            sage: s == b.contract(a)
            True

        Case of a vanishing scalar field result::

            sage: b = M.vector_field('b')
            sage: b[eU,:] = [1+x, -y]
            sage: b.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: s = a.contract(b) ; s
            Scalar field zero on the 2-dimensional differentiable manifold M
            sage: s.display()
            zero: M --> R
            on U: (x, y) |--> 0
            on V: (u, v) |--> 0

        """
        nargs = len(args)
        for i, arg in enumerate(args):
            if isinstance(arg, TensorField):
                other = arg
                it = i
                break
        else:
            raise TypeError("A tensor field must be provided in the " +
                            "argument list.")
        if it == 0:
            pos1 = (self._tensor_rank - 1,)
        else:
            pos1 = args[:it]
        if it == nargs-1:
            pos2 = (0,)
        else:
            pos2 = args[it+1:]
        ncontr = len(pos1) # number of contractions
        if len(pos2) != ncontr:
            raise TypeError("Different number of indices for the contraction.")
        if self._domain.is_subset(other._domain):
            if not self._ambient_domain.is_subset(other._ambient_domain):
                raise TypeError("Incompatible ambient domains for contraction.")
        elif other._domain.is_subset(self._domain):
            if not other._ambient_domain.is_subset(self._ambient_domain):
                raise TypeError("Incompatible ambient domains for contraction.")
        dom_resu = self._domain.intersection(other._domain)
        ambient_dom_resu = self._ambient_domain.intersection(
                                                         other._ambient_domain)
        self_r = self.restrict(dom_resu)
        other_r = other.restrict(dom_resu)
        k1, l1 = self._tensor_type
        k2, l2 = other._tensor_type
        tensor_type_resu = (k1+k2-ncontr, l1+l2-ncontr)
        if ambient_dom_resu.is_manifestly_parallelizable():
            # call of the FreeModuleTensor version:
            args = pos1 + (other_r,) + pos2
            return FreeModuleTensor.contract(self_r, *args)
        com_dom = []
        for dom in self_r._restrictions:
            if dom in other_r._restrictions:
                com_dom.append(dom)
        resu_rst = []
        for dom in com_dom:
            self_rr = self_r._restrictions[dom]
            other_rr = other_r._restrictions[dom]
            args = pos1 + (other_rr,) + pos2
            resu_rst.append(self_rr.contract(*args))
        if tensor_type_resu == (0,0):
            # scalar field result
            resu = dom_resu.scalar_field()
            all_zero = True
            for rst in resu_rst:
                if rst == 0:
                    for chart in rst._domain._atlas:
                        resu._express[chart] = 0
                else:
                    all_zero = False
                    for chart, funct in rst._express.iteritems():
                        resu._express[chart] = funct
            if all_zero:
                resu = dom_resu._zero_scalar_field
        else:
            # tensor field result
            dest_map = self._vmodule._dest_map
            dest_map_resu = dest_map.restrict(dom_resu,
                                              subcodomain=ambient_dom_resu)
            vmodule = dom_resu.vector_field_module(dest_map=dest_map_resu)

            resu = vmodule.tensor(tensor_type_resu, sym=resu_rst[0]._sym,
                                  antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def symmetrize(self, *pos):
        r"""
        Symmetrization over some arguments.

        INPUT:

        - ``pos`` -- (default: ``None``) list of argument positions involved in
          the symmetrization (with the convention position=0 for the first
          argument); if none, the symmetrization is performed over all the
          arguments

        OUTPUT:

        - the symmetrized tensor field (instance of :class:`TensorField`)

        EXAMPLES:

        Symmetrization of a type-(0,2) tensor field on a 2-dimensional
        non-parallelizable manifold::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(0,2, name='a')
            sage: a[eU,:] = [[1,x], [2,y]]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a[eV,:]
            [ 1/4*u + 3/4 -1/4*u + 3/4]
            [ 1/4*v - 1/4 -1/4*v - 1/4]
            sage: s = a.symmetrize() ; s
            Field of symmetric bilinear forms on the 2-dimensional
             differentiable manifold M
            sage: s[eU,:]
            [        1 1/2*x + 1]
            [1/2*x + 1         y]
            sage: s[eV,:]
            [         1/4*u + 3/4 -1/8*u + 1/8*v + 1/4]
            [-1/8*u + 1/8*v + 1/4         -1/4*v - 1/4]
            sage: s == a.symmetrize(0,1)  # explicit positions
            True

        See
        :meth:`sage.tensor.modules.free_module_tensor.FreeModuleTensor.symmetrize`
        for more details and examples.

        """
        resu_rst = []
        for rst in self._restrictions.itervalues():
            resu_rst.append(rst.symmetrize(*pos))
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_rst[0]._sym,
                                    antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def antisymmetrize(self, *pos):
        r"""
        Antisymmetrization over some arguments.

        INPUT:

        - ``pos`` -- (default: ``None``) list of argument positions involved in
          the antisymmetrization (with the convention position=0 for the first
          argument); if none, the antisymmetrization is performed over all the
          arguments

        OUTPUT:

        - the antisymmetrized tensor field (instance of :class:`TensorField`)

        EXAMPLES:

        Antisymmetrization of a type-(0,2) tensor field on a 2-dimensional
        non-parallelizable manifold::

            sage: M = DiffManifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y), intersection_name='W',
            ....:                              restrictions1= x>0, restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(0,2, name='a')
            sage: a[eU,:] = [[1,x], [2,y]]
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a[eV,:]
            [ 1/4*u + 3/4 -1/4*u + 3/4]
            [ 1/4*v - 1/4 -1/4*v - 1/4]
            sage: s = a.antisymmetrize() ; s
            2-form on the 2-dimensional differentiable manifold M
            sage: s[eU,:]
            [         0  1/2*x - 1]
            [-1/2*x + 1          0]
            sage: s[eV,:]
            [                   0 -1/8*u - 1/8*v + 1/2]
            [ 1/8*u + 1/8*v - 1/2                    0]
            sage: s == a.antisymmetrize(0,1)  # explicit positions
            True
            sage: s == a.antisymmetrize(1,0)  # the order of positions does not matter
            True

        See
        :meth:`sage.tensor.modules.free_module_tensor.FreeModuleTensor.antisymmetrize`
        for more details and examples.

        """
        resu_rst = []
        for rst in self._restrictions.itervalues():
            resu_rst.append(rst.antisymmetrize(*pos))
        resu = self._vmodule.tensor(self._tensor_type, sym=resu_rst[0]._sym,
                                    antisym=resu_rst[0]._antisym)
        for rst in resu_rst:
            resu._restrictions[rst._domain] = rst
        return resu

    def lie_der(self, vector):
        r"""
        Compute the Lie derivative with respect to a vector field.

        The Lie derivative is stored in the dictionary
        :attr:`_lie_derivatives`, so that there is no need to
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile.

        INPUT:

        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken

        OUTPUT:

        - the tensor field that is the Lie derivative of ``self`` with respect
          to ``vector``

        EXAMPLES:

        Lie derivative of a type-(1,1) tensor field along a vector field on
        the 2-sphere::

            sage: M = DiffManifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e_xy,:] = [[x, 1], [y, 0]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: w = M.vector_field(name='w')
            sage: w[e_xy,:] = [y^2, x^2]
            sage: w.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: lt = t.lie_der(w); lt
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: lt.display(e_xy)
            (-y^2 + 2*x) d/dx*dx + 2*x*y d/dx*dy - x^2 d/dy*dx
             + (2*y^2 - 2*x) d/dy*dy
            sage: lt.display(e_uv)
            (2*u^7 + (2*u - 1)*v^6 - 2*u^5*v + 2*u^3*v^3 - 2*(5*u^3 - 3*u^2)*v^4 - (10*u^5 - 3*u^4)*v^2)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) d/du*du
             - 2*(2*u^4*v^2 + u^3*v^3 + 2*(2*u^2 - u)*v^5 - (4*u^6 - u^5)*v)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) d/du*dv
             + (8*u^6*v + u^6 - 2*u^4*v^2 + 2*u^3*v^3 + u^2*v^4 - 2*(4*u^2 - 3*u)*v^5)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) d/dv*du
             - 2*(u^7 + (u - 1)*v^6 - u^5*v + u^3*v^3 - (5*u^3 - 2*u^2)*v^4 - (5*u^5 - u^4)*v^2)/(u^8 + 4*u^6*v^2 + 6*u^4*v^4 + 4*u^2*v^6 + v^8) d/dv*dv

        The result is cached::

            sage: t.lie_der(w) is lt
            True

        Lie derivative of a vector field::

            sage: a = M.vector_field(name='a')
            sage: a[e_xy,:] = [1-x, x-y]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: a.lie_der(w)
            Vector field on the 2-dimensional differentiable manifold M
            sage: a.lie_der(w).display(e_xy)
            (-2*x*y + y^2) d/dx + (x^2 + y^2 - 2*x) d/dy
            sage: a.lie_der(w).display(e_uv)
            (4*u^4*v - u^2*v^2 + 4*(u^2 - u)*v^3 + v^4)/(u^4 + 2*u^2*v^2 + v^4) d/du - (2*u^5 - (2*u - 1)*v^4 - u^4 - 4*u^2*v^2 + 2*u*v^3)/(u^4 + 2*u^2*v^2 + v^4) d/dv

        The Lie derivative is antisymmetric:

            sage: a.lie_der(w) == - w.lie_der(a)
            True

        and it coincides with the commutator of the two vector fields::

            sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2), c_uv: (u^2+v^2)/(1+u^2+v^2)})
            sage: a.lie_der(w)(f) == w(a(f)) - a(w(f))  # long time
            True

        """
        if vector._tensor_type != (1,0):
            raise TypeError("The argument must be a vector field.")
        if id(vector) not in self._lie_derivatives:
            # the computation must be performed:
            resu_rst = []
            for dom, rst in self._restrictions.iteritems():
                resu_rst.append(rst.lie_der(vector.restrict(dom)))
            resu = self._vmodule.tensor(self._tensor_type,
                                        sym=resu_rst[0]._sym,
                                        antisym=resu_rst[0]._antisym)
            for rst in resu_rst:
                resu._restrictions[rst._domain] = rst
            self._lie_derivatives[id(vector)] = (vector, resu)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]


#******************************************************************************

class TensorFieldParal(FreeModuleTensor, TensorField):
    r"""
    Tensor field along an open set of a differentiable manifold,
    with values on a parallelizable open set of a differentiable manifold.

    An instance of this class is a tensor field along an open subset `U`
    of some differentiable manifold `S` with values in a parallelizable open
    subset `V` of a differentiable manifold `M`, via a differentiable map
    `\Phi: U \rightarrow V`.
    The standard case of a tensor field *on* a differentiable manifold
    corresponds to `S=M`, `U=V` and `\Phi = \mathrm{Id}_U`. Other common cases
    are `\Phi` being an immersion and `\Phi` being a curve in `V` (`U` is then
    an open interval of `\RR`).

    A tensor field of type `(k,l)` is a field `t` on `U`, such that at each
    point `p\in U`, `t(p)` is a multilinear map

    .. MATH::

        t(p):\ \underbrace{T_q^*M\times\cdots\times T_q^*M}_{k\ \; \mbox{times}}
        \times \underbrace{T_q M\times\cdots\times T_q M}_{l\ \; \mbox{times}}
        \longrightarrow \RR

    where `T_q M` stands for the tangent space to the manifold `M` at the point
    `q=\Phi(p)` and `T_q^* M` for its dual vector space. The integer `k+l`
    is called the tensor rank.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldFreeModule`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the map `\Phi: U \rightarrow V` (cf.
      :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`)
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank
      and `l` the covariant rank
    - ``name`` -- (default: ``None``) name given to the tensor field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the tensor
      field; if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
      the tensor arguments: each symmetry is described by a tuple containing
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:

      * sym=(0,1) for a symmetry between the 1st and 2nd arguments
      * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
        arguments and a symmetry between the 2nd, 4th and 5th arguments.

    - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
      among the arguments, with the same convention as for ``sym``.


    EXAMPLES:

    A tensor field of type (2,0) on a 3-dimensional manifold::

        sage: M = DiffManifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart()
        sage: t = M.tensor_field(2, 0, 'T') ; t
        Tensor field T of type (2,0) on the 3-dimensional differentiable
         manifold M

    Tensor fields are considered as elements of a module over the ring
    `C^k(M)` of scalar fields on `M`::

        sage: t.parent()
        Free module T^(2,0)(M) of type-(2,0) tensors fields on the
         3-dimensional differentiable manifold M
        sage: t.parent().base_ring()
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M

    The components with respect to the manifold's default frame are set or read
    by means of square brackets::

        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: for i in M.irange():
        ...       for j in M.irange():
        ...           t[i,j] = (i+1)**(j+1)
        ...
        sage: [[ t[i,j] for j in M.irange()] for i in M.irange()]
        [[1, 1, 1], [2, 4, 8], [3, 9, 27]]

    A shortcut for the above is using [:]::

        sage: t[:]
        [ 1  1  1]
        [ 2  4  8]
        [ 3  9 27]

    The components with respect to another frame are set via the method
    :meth:`set_comp` and read via the method :meth:`comp`; both return an
    instance of :class:`~sage.tensor.modules.comp.Components`::

        sage: f = M.vector_frame('f')  # a new frame defined on M, in addition to e
        sage: t.set_comp(f)[0,0] = -3
        sage: t.comp(f)
        2-indices components w.r.t. Vector frame (M, (f_0,f_1,f_2))
        sage: t.comp(f)[0,0]
        -3
        sage: t.comp(f)[:]  # the full list of components
        [-3  0  0]
        [ 0  0  0]
        [ 0  0  0]

    To avoid any insconstency between the various components, the method
    :meth:`set_comp` deletes the components in other frames.
    Accordingly, the components in the frame e have been deleted::

        sage: t._components
        {Vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. Vector
         frame (M, (f_0,f_1,f_2))}

    To keep the other components, one must use the method :meth:`add_comp`::

        sage: t = M.tensor_field(2, 0, 'T')  # Let us restart
        sage: t[0,0] = 2                     # sets the components in the frame e
        sage: # We now set the components in the frame f with add_comp:
        sage: t.add_comp(f)[0,0] = -3
        sage: # The components w.r.t. frame e have been kept:
        sage: t._components  # random (dictionary output)
        {Vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. Vector frame (M, (e_0,e_1,e_2)),
         Vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. Vector frame (M, (f_0,f_1,f_2))}

    The basic properties of a tensor field are::

        sage: t.domain()
        3-dimensional differentiable manifold M
        sage: t.tensor_type()
        (2, 0)

    Internally, the components w.r.t. various vector frames are stored in the
    dictionary :attr:`_components`::

        sage: t._components  # random (dictionary output)
        {Vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. Vector frame (M, (e_0,e_1,e_2)),
         Vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. Vector frame (M, (f_0,f_1,f_2))}

    Symmetries and antisymmetries are declared via the keywords ``sym`` and
    ``antisym``. For instance, a rank-6 covariant tensor that is symmetric with
    respect to its 1st and 3rd arguments and antisymmetric with respect to the
    2nd, 5th and 6th arguments is set up as follows::

        sage: a = M.tensor_field(0, 6, 'T', sym=(0,2), antisym=(1,4,5))
        sage: a[0,0,1,0,1,2] = 3
        sage: a[1,0,0,0,1,2] # check of the symmetry
        3
        sage: a[0,1,1,0,0,2], a[0,1,1,0,2,0] # check of the antisymmetry
        (-3, 3)

    Multiple symmetries or antisymmetries are allowed; they must then be
    declared as a list. For instance, a rank-4 covariant tensor that is
    antisymmetric with respect to its 1st and 2nd arguments and with respect to
    its 3rd and 4th argument must be declared as::

        sage: r = M.tensor_field(0, 4, 'T', antisym=[(0,1), (2,3)])
        sage: r[0,1,2,0] = 3
        sage: r[1,0,2,0] # first antisymmetry
        -3
        sage: r[0,1,0,2] # second antisymmetry
        -3
        sage: r[1,0,0,2] # both antisymmetries acting
        3

    Tensor fields of the same type can be added and subtracted::

        sage: a = M.tensor_field(2, 0)
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: b = M.tensor_field(2, 0)
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + 2*b ; s
        Tensor field of type (2,0) on the 3-dimensional differentiable
         manifold M
        sage: a[:], (2*b)[:], s[:]
        (
        [1 2 3]  [ 8  0 14]  [ 9  2 17]
        [0 0 0]  [ 0 10  0]  [ 0 10  0]
        [0 0 0], [ 0  0 12], [ 0  0 12]
        )
        sage: s = a - b ; s
        Tensor field of type (2,0) on the 3-dimensional differentiable
         manifold M
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [-3  2 -4]
        [0 0 0]  [0 5 0]  [ 0 -5  0]
        [0 0 0], [0 0 6], [ 0  0 -6]
        )

    Symmetries are preserved by the addition whenever it is possible::

        sage: a = M.tensor_field(2, 0, sym=(0,1))
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [0 0 6], [ 3  0  6]
        )
        sage: a.symmetries()
        symmetry: (0, 1);  no antisymmetry
        sage: b.symmetries()
        no symmetry;  no antisymmetry
        sage: s.symmetries()
        no symmetry;  no antisymmetry
        sage: # let us now make b symmetric:
        sage: b = M.tensor_field(2, 0, sym=(0,1))
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [7 0 6], [10  0  6]
        )
        sage: s.symmetries()  # s is symmetric because both a and b are
        symmetry: (0, 1);  no antisymmetry

    The tensor product is taken with the operator \*::

        sage: c = a*b ; c
        Tensor field of type (4,0) on the 3-dimensional differentiable
         manifold M
        sage: c.symmetries()  # since a and b are both symmetric, a*b has two symmetries:
        symmetries: [(0, 1), (2, 3)];  no antisymmetry

    The tensor product of two fully contravariant tensors is not symmetric in
    general::

        sage: a*b == b*a
        False

    The tensor product of a fully contravariant tensor by a fully covariant one
    is symmetric::

        sage: d = M.diff_form(2)  # a fully covariant tensor field
        sage: d[0,1], d[0,2], d[1,2] = (3, 2, 1)
        sage: s = a*d ; s
        Tensor field of type (2,2) on the 3-dimensional differentiable
         manifold M
        sage: s.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: s1 = d*a ; s1
        Tensor field of type (2,2) on the 3-dimensional differentiable
         manifold M
        sage: s1.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: d*a == a*d
        True

    Example of tensor field associated with a non-trivial differentiable
    map `\Phi`: tensor field along a curve in `M`::

        sage: R = DiffManifold(1, 'R')  # R as a 1-dimensional manifold
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: Phi = R.diff_map(M, [cos(t), sin(t), t], name='Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold R
         to the 3-dimensional differentiable manifold M
        sage: h = R.tensor_field(2, 0, name='h', dest_map=Phi) ; h
        Tensor field h of type (2,0) along the 1-dimensional differentiable
         manifold R with values on the 3-dimensional differentiable manifold M
        sage: h.parent()
        Free module T^(2,0)(R,Phi) of type-(2,0) tensors fields along the
         1-dimensional differentiable manifold R mapped into the 3-dimensional
         differentiable manifold M
        sage: h[0,0], h[0,1], h[2,0] = 1+t, t^2, sin(t)
        sage: h.display()
        h = (t + 1) d/dx*d/dx + t^2 d/dx*d/dy + sin(t) d/dz*d/dx

    """
    def __init__(self, vector_field_module, tensor_type, name=None,
                 latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor field.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``TensorFieldParal``, to fit with the category framework::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: XM = M.vector_field_module()
            sage: T02 = M.tensor_field_module((0,2))
            sage: t = T02.element_class(XM, (0,2), name='t'); t
            Tensor field t of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: t[:] = [[1+x^2, x*y], [0, 1+y^2]]
            sage: t.display()
            t = (x^2 + 1) dx*dx + x*y dx*dy + (y^2 + 1) dy*dy
            sage: t.parent()
            Free module T^(0,2)(M) of type-(0,2) tensors fields on the
             2-dimensional differentiable manifold M
            sage: TestSuite(t).run()

        """
        FreeModuleTensor.__init__(self, vector_field_module, tensor_type,
                                  name=name, latex_name=latex_name,
                                  sym=sym, antisym=antisym)
        # TensorField attributes:
        self._vmodule = vector_field_module
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: the TensorField attribute self._restrictions is considered as a
        #     derived quantity in the present case (the primary attribute
        #     being self._components, which is initialized by
        #     FreeModuleTensor.__init__ ); accordingly self._restrictions is
        #     initialized by _init_derived() and cleared by _del_derived().
        # Initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._repr_()
            'Tensor field t of type (1,1) on the 2-dimensional differentiable manifold M'
            sage: repr(t)  # indirect doctest
            'Tensor field t of type (1,1) on the 2-dimensional differentiable manifold M'
            sage: t  # indirect doctest
            Tensor field t of type (1,1) on the 2-dimensional differentiable
             manifold M

        """
        return TensorField._repr_(self)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same
        vector field module, with the same tensor type and same symmetries.

        TESTS::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._new_instance()
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: type(t._new_instance()) is type(t)
            True

        """
        return self.__class__(self._fmodule, self._tensor_type, sym=self._sym,
                                antisym=self._antisym)

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TEST::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._init_derived()

        """
        FreeModuleTensor._init_derived(self)
        TensorField._init_derived(self)
        self._restrictions = {} # dict. of restrictions of self on subdomains
                                # of self._domain, with the subdomains as keys

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities

        INPUT:

        - ``del_restrictions`` -- (default: ``True``) determines whether the
          restrictions of ``self`` to subdomains are deleted.

        TEST::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._del_derived()

        """
        FreeModuleTensor._del_derived(self)
        TensorField._del_derived(self)
        if del_restrictions:
            self._restrictions.clear()

    def set_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame
        for assignment.

        The components with respect to other frames on the same domain are
        deleted, in order to avoid any inconsistency. To keep them, use the
        method :meth:`add_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are defined; if none is provided, the components are assumed to refer
          to the tensor field domain's default frame.

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        EXAMPLES::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e_xy = X.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t.set_comp(e_xy)
            2-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.set_comp(e_xy)[1,0] = 2
            sage: t.display(e_xy)
            t = 2 d/dy*dx

        Setting components in a new frame (``e``)::

            sage: e = M.vector_frame('e')
            sage: t.set_comp(e)
            2-indices components w.r.t. Vector frame (M, (e_0,e_1))
            sage: t.set_comp(e)[0,1] = x
            sage: t.display(e)
            t = x e_0*e^1

        The components w.r.t. the frame ``e_xy`` have be erased::

            sage: t.display(e_xy)
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (M, (d/dx,d/dy))

        Setting components in a frame defined on a subdomain deletes previously
        defined components as well::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: f = U.vector_frame('f')
            sage: t.set_comp(f)
            2-indices components w.r.t. Vector frame (U, (f_0,f_1))
            sage: t.set_comp(f)[0,1] = 1+y
            sage: t.display(f)
            t = (y + 1) f_0*f^1
            sage: t.display(e)
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Vector frame (M, (e_0,e_1))

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis._domain == self._domain:
            # Setting components on the tensor field domain:
            return FreeModuleTensor.set_comp(self, basis=basis)
        else:
            # Setting components on a subdomain:
            #
            # Creating or saving the restriction to the subdomain:
            rst = self.restrict(basis._domain, dest_map=basis._dest_map)
            # Deleting all the components on self._domain and the derived
            # quantities:
            self._components.clear()
            self._del_derived()
            # Restoring the restriction to the subdomain (which has been
            # deleted by _del_derived):
            self._restrictions[basis._domain] = rst
            # The set_comp operation is performed on the subdomain:
            return rst.set_comp(basis=basis)

    def add_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame
        for assignment.

        The components with respect to other frames on the same domain are
        kept. To delete them them, use the method :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are defined; if none is provided, the components are assumed to refer
          to the tensor field domain's default frame.

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created.

        EXAMPLES::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e_xy = X.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t.add_comp(e_xy)
            2-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.add_comp(e_xy)[1,0] = 2
            sage: t.display(e_xy)
            t = 2 d/dy*dx

        Adding components w.r.t. a new frame (``e``)::

            sage: e = M.vector_frame('e')
            sage: t.add_comp(e)
            2-indices components w.r.t. Vector frame (M, (e_0,e_1))
            sage: t.add_comp(e)[0,1] = x
            sage: t.display(e)
            t = x e_0*e^1

        The components w.r.t. frame ``e_xy`` are kept::

            sage: t.display(e_xy)
            t = 2 d/dy*dx

        Adding components in a frame defined on a subdomain::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: f = U.vector_frame('f')
            sage: t.add_comp(f)
            2-indices components w.r.t. Vector frame (U, (f_0,f_1))
            sage: t.add_comp(f)[0,1] = 1+y
            sage: t.display(f)
            t = (y + 1) f_0*f^1

        The components previously defined are kept::

            sage: t.display(e_xy)
            t = 2 d/dy*dx
            sage: t.display(e)
            t = x e_0*e^1

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis._domain == self._domain:
            # Adding components on the tensor field domain:
            return FreeModuleTensor.add_comp(self, basis=basis)
        else:
            # Adding components on a subdomain:
            #
            # Creating or saving the restriction to the subdomain:
            rst = self.restrict(basis._domain, dest_map=basis._dest_map)
            # Deleting the derived quantities except for the restrictions to
            # subdomains:
            self._del_derived(del_restrictions=False)
            # The add_comp operation is performed on the subdomain:
            return rst.add_comp(basis=basis)

    def comp(self, basis=None, from_basis=None):
        r"""
        Return the components in a given vector frame.

        If the components are not known already, they are computed by the
        tensor change-of-basis formula from components in another vector frame.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the components
          are required; if none is provided, the components are assumed to
          refer to the tensor field domain's default frame
        - ``from_basis`` -- (default: ``None``) vector frame from which the
          required components are computed, via the tensor change-of-basis
          formula, if they are not known already in the basis ``basis``

        OUTPUT:

        - components in the vector frame ``basis``, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES::

            sage: M = DiffManifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(1,2, name='t')
            sage: t[1,2,1] = x*y
            sage: t.comp(X.frame())
            3-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.comp()  # the default frame is X.frame()
            3-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.comp()[:]
            [[[0, 0], [x*y, 0]], [[0, 0], [0, 0]]]
            sage: e = M.vector_frame('e')
            sage: t[e, 2,1,1] = x-3
            sage: t.comp(e)
            3-indices components w.r.t. Vector frame (M, (e_1,e_2))
            sage: t.comp(e)[:]
            [[[0, 0], [0, 0]], [[x - 3, 0], [0, 0]]]

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis._domain == self._domain:
            # components on the tensor field domain:
            return FreeModuleTensor.comp(self, basis=basis,
                                         from_basis=from_basis)
        else:
            # components on a subdomain:
            rst = self.restrict(basis._domain, dest_map=basis._dest_map)
            return rst.comp(basis=basis, from_basis=from_basis)


    def common_coord_frame(self, other):
        r"""
        Find a common coordinate frame for the components of ``self`` and
        ``other``.

        In case of multiple common bases, the domain's default coordinate
        basis is privileged.
        If the current components of ``self`` and ``other`` are all relative to
        different frames, a common frame is searched by performing a component
        transformation, via the transformations listed in
        ``self._domain._frame_changes``, still privileging transformations to
        the domain's default frame.

        INPUT:

        - ``other`` -- a tensor field (instance of :class:`TensorFieldParal`)

        OUPUT:

        - common coordinate frame; if no common basis is found, None is
          returned.

        EXAMPLES::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.tensor_field(1,2, name='a')
            sage: a[0,1,0] = 2
            sage: b = M.vector_field(name='b')
            sage: b[:] = [-y, x]
            sage: a.common_coord_frame(b)
            Coordinate frame (M, (d/dx,d/dy))

        Vector field defined on a new chart::

            sage: Y.<u,v> = M.chart()
            sage: c = M.vector_field(name='c')
            sage: c[Y.frame(), :, Y] = (1+u, u*v)
            sage: c.display(Y.frame(), Y)
            c = (u + 1) d/du + u*v d/dv

        There is no common ccordinate frame::

            sage: a.common_coord_frame(c)

        Connecting the two coordinate charts enables to find a common frame::

            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: Y_to_X = X_to_Y.inverse()
            sage: a.common_coord_frame(c)
            Coordinate frame (M, (d/dx,d/dy))

        Indeed, the components of ``c`` w.r.t. the frame ``(M, (d/dx,d/dy))``
        have been computed via the change-of-coordinate formulas::

            sage: c.display(a.common_coord_frame(c))
            c = (1/2*x^2 - 1/2*y^2 + 1/2*x + 1/2*y + 1/2) d/dx
             + (-1/2*x^2 + 1/2*y^2 + 1/2*x + 1/2*y + 1/2) d/dy

        """
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        # Compatibility checks:
        if not isinstance(other, TensorFieldParal):
            raise TypeError("The argument must be of type TensorFieldParal.")
        dom = self._domain
        def_frame = dom._def_frame
        #
        # 1/ Search for a common frame among the existing components, i.e.
        #    without performing any component transformation.
        #    -------------------------------------------------------------
        # 1a/ Direct search
        if def_frame in self._components and \
           def_frame in other._components and \
           isinstance(dom._def_frame, CoordFrame):
            return def_frame # the domain's default frame is privileged
        for frame1 in self._components:
            if frame1 in other._components and \
               isinstance(frame1, CoordFrame):
                return frame1
        # 1b/ Search involving subframes
        dom2 = other._domain
        for frame1 in self._components:
            if not isinstance(frame1, CoordFrame):
                continue
            for frame2 in other._components:
                if not isinstance(frame2, CoordFrame):
                    continue
                if frame2 in frame1._subframes:
                    self.comp(frame2)
                    return frame2
                if frame1 in frame2._subframes:
                    other.comp(frame1)
                    return frame1
        #
        # 2/ Search for a common frame via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least
        # one component transformation to get a common frame
        if isinstance(dom._def_frame, CoordFrame):
            if def_frame in self._components:
                for oframe in other._components:
                    if (oframe, def_frame) in dom._frame_changes:
                        other.comp(def_frame, from_basis=oframe)
                        return def_frame
            if def_frame in other._components:
                for sframe in self._components:
                    if (sframe, def_frame) in dom._frame_changes:
                        self.comp(def_frame, from_basis=sframe)
                        return def_frame
        # If this point is reached, then def_frame cannot be a common frame
        # via a single component transformation
        for sframe in self._components:
            if not isinstance(sframe, CoordFrame):
                continue
            for oframe in other._components:
                if not isinstance(oframe, CoordFrame):
                    continue
                if (oframe, sframe) in dom._frame_changes:
                    other.comp(sframe, from_basis=oframe)
                    return sframe
                if (sframe, oframe) in dom._frame_changes:
                    self.comp(oframe, from_basis=sframe)
                    return oframe
        #
        # 3/ Search for a common frame via two component transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least
        # two component transformations to get a common frame
        for sframe in self._components:
            for oframe in other._components:
                if (sframe, def_frame) in dom._frame_changes and \
                   (oframe, def_frame) in dom._frame_changes and \
                   isinstance(def_frame, CoordFrame):
                    self.comp(def_frame, from_basis=sframe)
                    other.comp(def_frame, from_basis=oframe)
                    return def_frame
                for frame in dom._frames:
                    if (sframe, frame) in dom._frame_changes and \
                       (oframe, frame) in dom._frame_changes and \
                       isinstance(frame, CoordFrame):
                        self.comp(frame, from_basis=sframe)
                        other.comp(frame, from_basis=oframe)
                        return frame
        #
        # If this point is reached, no common frame could be found, even at
        # the price of component transformations:
        return None


    def lie_der(self, vector):
        r"""
        Compute the Lie derivative with respect to a vector field.

        The Lie derivative is stored in the dictionary
        :attr:`_lie_derivatives`, so that there is no need to
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile.

        INPUT:

        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken

        OUTPUT:

        - the tensor field that is the Lie derivative of ``self`` with respect
          to ``vector``

        EXAMPLES:

        Lie derivative of a vector::

            sage: M = DiffManifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: v = M.vector_field('v')
            sage: v[:] = (-y, x)
            sage: w = M.vector_field()
            sage: w[:] = (2*x+y, x*y)
            sage: w.lie_der(v)
            Vector field on the 2-dimensional differentiable manifold M
            sage: w.lie_der(v).display()
            ((x - 2)*y + x) d/dx + (x^2 - y^2 - 2*x - y) d/dy

        The Lie derivative is antisymmetric::

            sage: w.lie_der(v) == -v.lie_der(w)
            True

        For vectors, it coincides with the commutator::

            sage: f = M.scalar_field(x^3 + x*y^2)
            sage: w.lie_der(v)(f).display()
            M --> R
            (x, y) |--> -(x + 2)*y^3 + 3*x^3 - x*y^2 + 5*(x^3 - 2*x^2)*y
            sage: w.lie_der(v)(f) == v(w(f)) - w(v(f))  # rhs = commutator [v,w] acting on f
            True

        Lie derivative of a 1-form::

            sage: om = M.one_form()
            sage: om[:] = (y^2*sin(x), x^3*cos(y))
            sage: om.lie_der(v)
            1-form on the 2-dimensional differentiable manifold M
            sage: om.lie_der(v).display()
            (-y^3*cos(x) + x^3*cos(y) + 2*x*y*sin(x)) dx
             + (-x^4*sin(y) - 3*x^2*y*cos(y) - y^2*sin(x)) dy

        Check of Cartan identity::

            sage: om.lie_der(v) == v.contract(0,om.exterior_der(),0) + (om(v)).exterior_der()
            True

        """
        if vector._tensor_type != (1,0):
            raise TypeError("The argument must be a vector field.")
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            #
            # 1/ Search for a common coordinate frame:
            coord_frame = self.common_coord_frame(vector)
            if coord_frame is None:
                raise TypeError("No common coordinate frame found.")
            chart = coord_frame._chart
            #
            # 2/ Component computation:
            tc = self._components[coord_frame]
            vc = vector._components[coord_frame]
            # the result has the same tensor type and same symmetries as self:
            resc = self._new_comp(coord_frame)
            n_con = self._tensor_type[0]
            vf_module = vector._fmodule
            for ind in resc.non_redundant_index_generator():
                rsum = 0
                for i in vf_module.irange():
                    rsum += vc[[i]].coord_function(chart) * \
                           tc[[ind]].coord_function(chart).diff(i)
                # loop on contravariant indices:
                for k in range(n_con):
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i
                        rsum -= tc[[indk]].coord_function(chart) * \
                                vc[[ind[k]]].coord_function(chart).diff(i)
                # loop on covariant indices:
                for k in range(n_con, self._tensor_rank):
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i
                        rsum += tc[[indk]].coord_function(chart) * \
                                vc[[i]].coord_function(chart).diff(ind[k])
                resc[[ind]] = rsum.scalar_field()
            #
            # 3/ Final result (the tensor)
            res = vf_module.tensor_from_comp(self._tensor_type, resc)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]

    def restrict(self, subdomain, dest_map=None):
        r"""
        Return the restriction of ``self`` to some subset of its domain.

        If such restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `U` of ``self._domain`` (must be an
          instance of :class:`~sage.manifolds.differentiable.manifold.DiffManifold`)
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `V` is a subset of
          ``self._codomain``
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          If None, the restriction of ``self._vmodule._dest_map`` to `U` is
          used.

        OUTPUT:

        - instance of :class:`TensorFieldParal` representing the restriction.

        EXAMPLES:

        Restriction of a vector field defined on `\RR^2` to a disk::

            sage: M = DiffManifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: v = M.vector_field('v')
            sage: v[:] = [x+y, -1+x^2]
            sage: D = M.open_subset('D') # the unit open disc
            sage: c_cart_D = c_cart.restrict(D, x^2+y^2<1)
            sage: v_D = v.restrict(D) ; v_D
            Vector field v on the Open subset D of the 2-dimensional
             differentiable manifold R^2
            sage: v_D.display()
            v = (x + y) d/dx + (x^2 - 1) d/dy

        The symbolic expressions of the components w.r.t. Cartesian coordinates
        are equal::

            sage: bool( v_D[1].expr() == v[1].expr() )
            True

        but neither the chart functions representing the components (they are
        defined on different charts)::

            sage: v_D[1] == v[1]
            False

        nor the scalar fields representing the components (they are
        defined on different open subsets)::

            sage: v_D[[1]] == v[[1]]
            False

        The restriction of the vector field to its own domain is of course
        itself::

            sage: v.restrict(M) is v
            True

        """
        if subdomain == self._domain and \
                    (dest_map is None or dest_map == self._vmodule._dest_map) :
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("The provided domain is not a subset of " +
                                 "the field's domain.")
            if dest_map is None:
                dest_map = self._fmodule._dest_map.restrict(subdomain)
            elif not dest_map._codomain.is_subset(self._ambient_domain):
                raise ValueError("Argument dest_map not compatible with " +
                                 "self._ambient_domain")
            #!# First one tries to derive the restriction from a tighter domain:
            #for dom, rst in self._restrictions.iteritems():
            #    if subdomain.is_subset(dom):
            #        self._restrictions[subdomain] = rst.restrict(subdomain)
            #        break
            # If this fails, the restriction is created from scratch:
            #else:
            smodule = subdomain.vector_field_module(dest_map=dest_map)
            resu = smodule.tensor(self._tensor_type, name=self._name,
                                  latex_name=self._latex_name, sym=self._sym,
                                  antisym=self._antisym,
                                  specific_type=self.__class__)
            for frame in self._components:
                for sframe in subdomain._covering_frames:
                    if sframe in frame._subframes:
                        comp_store = self._components[frame]._comp
                        scomp = resu._new_comp(sframe)
                        scomp_store = scomp._comp
                        # the components of the restriction are evaluated
                        # index by index:
                        for ind, value in comp_store.iteritems():
                            scomp_store[ind] = value.restrict(subdomain)
                        resu._components[sframe] = scomp
            self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def __call__(self, *args):
        r"""
        The tensor field acting on 1-forms and vector fields as a multilinear
        map.

        In the particular case of tensor field of type (1,1), the action can
        be on a single vector field, the tensor field being identified to a
        field of tangent-space endomorphisms. The output is then a vector
        field.

        INPUT:

        - ``*args`` -- list of k 1-forms and l vector fields, ``self``
          being a tensor of type (k,l).

        OUTPUT:

        - either the scalar field resulting from the action of ``self`` on
          the 1-forms and vector fields passed as arguments or the vector
          field resulting from the action of ``self`` as a field of
          tangent-space endomorphisms (case of a type-(1,1) tensor field)

        TESTS:

        Action of a tensor field of type (1,1)::

            sage: DiffManifold._clear_cache_() # for doctests only
            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[:] = [[1+x, 2], [y, -x^2]]
            sage: v = M.vector_field(name='v')
            sage: v[:] = [-y, x]
            sage: a = M.one_form(name='a')
            sage: a[:] = [3, 1-y]
            sage: s = t.__call__(a,v); s
            Scalar field t(a,v) on the 2-dimensional differentiable manifold M
            sage: s.display()
            t(a,v): M --> R
               (x, y) |--> -x^3 + y^3 + (x^3 - 3*x - 3)*y - y^2 + 6*x
            sage: s.coord_function() == sum(sum(t[i,j]*a[i]*v[j] for j in [0..1])
            ....:                           for i in [0..1])
            True
            sage: s == t(a,v)  # indirect doctest
            True

        The tensor field acting on vector field, as a field of tangent-space
        endomorphisms::

            sage: s = t.__call__(v); s
            Vector field t(v) on the 2-dimensional differentiable manifold M
            sage: s.display()
            t(v) = (-(x + 1)*y + 2*x) d/dx + (-x^3 - y^2) d/dy
            sage: s[0] == t[0,0]*v[0] + t[0,1]*v[1]
            True
            sage: s[1] == t[1,0]*v[0] + t[1,1]*v[1]
            True
            sage: s == t(v)  # indirect doctest
            True

        """
        from sage.categories.homset import End
        p = len(args)
        if p == 1 and self._tensor_type == (1,1):
            # type-(1,1) tensor acting as an endomorphism:
            vector = args[0]
            if vector._tensor_type != (1,0):
                raise TypeError("The argument must be a vector field.")
            dom = self._domain.intersection(vector._domain)
            sd = self.restrict(dom)
            vd = vector.restrict(dom)
            endom = End(vd.parent())(sd)
            return endom(vd)
        # Generic case
        if p != self._tensor_rank:
            raise TypeError(str(self._tensor_rank) +
                            " arguments must be provided.")
        # Domain of the result
        dom_resu = self._domain
        for arg in args:
            dom_resu = dom_resu.intersection(arg._domain)
        # Restriction to the result domain
        self_r = self.restrict(dom_resu)
        args_r = [args[i].restrict(dom_resu) for i in range(p)]
        # Call of the FreeModuleTensor version
        return FreeModuleTensor.__call__(self_r, *args_r)

    def contract(self, *args):
        r"""
        Contraction with another tensor field, on one or more indices.

        INPUT:

        - ``pos1`` -- positions of the indices in ``self`` involved in the
          contraction; ``pos1`` must be a sequence of integers, with 0 standing
          for the first index position, 1 for the second one, etc. If ``pos1``
          is not provided, a single contraction on the last index position of
          ``self`` is assumed
        - ``other`` -- the tensor field to contract with
        - ``pos2`` -- positions of the indices in ``other`` involved in the
          contraction, with the same conventions as for ``pos1``. If ``pos2``
          is not provided, a single contraction on the first index position of
          ``other`` is assumed

        OUTPUT:

        - tensor field resulting from the contraction at the positions
          ``pos1`` and ``pos2`` of ``self`` with ``other``

        EXAMPLES:

        Contraction of a tensor field of type (2,0) with a tensor field of
        type (1,1)::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.tensor_field(2,0, name='a')
            sage: a[:] = [[1+x, 2], [y, -x^2]]
            sage: b = M.tensor_field(1,1, name='b')
            sage: b[:] = [[-y, 1], [x, x+y]]
            sage: s = a.contract(0, b, 1); s
            Tensor field of type (2,0) on the 2-dimensional differentiable manifold M
            sage: s.display()
            -x*y d/dx*d/dx + (x^2 + x*y + y^2 + x) d/dx*d/dy
             + (-x^2 - 2*y) d/dy*d/dx + (-x^3 - x^2*y + 2*x) d/dy*d/dy

        Check::

            sage: all([s[ind] == sum(a[k, ind[0]]*b[ind[1], k] for k in [0..1])
            ....:      for ind in M.index_generator(2)])
            True

        The same contraction with repeated index notation::

            sage: s == a['^ki']*b['^j_k']
            True

        Contraction on the second index of ``a``::

            sage: s = a.contract(1, b, 1); s
            Tensor field of type (2,0) on the 2-dimensional differentiable manifold M
            sage: s.display()
            (-(x + 1)*y + 2) d/dx*d/dx + (x^2 + 3*x + 2*y) d/dx*d/dy
             + (-x^2 - y^2) d/dy*d/dx + (-x^3 - (x^2 - x)*y) d/dy*d/dy

        Check::

            sage: all([s[ind] == sum(a[ind[0], k]*b[ind[1], k] for k in [0..1])
            ....:      for ind in M.index_generator(2)])
            True

        The same contraction with repeated index notation::

            sage: s == a['^ik']*b['^j_k']
            True

        See
        :meth:`sage.manifolds.differentiable.tensorfield.TensorField.contract`
        for more examples.

        """
        # This is to ensure the call to the TensorField version instead of
        # the FreeModuleTensor one
        return TensorField.contract(self, *args)

    def __mul__(self, other):
        r"""
        Tensor product (or multiplication of the right by a scalar).

        INPUT:

        - ``other`` -- a tensor field, on the same manifold as ``self`` (or an
          object that can be coerced to a scalar field on the same manifold
          as ``self``)

        OUPUT:

        - the tensor field resulting from the tensor product of ``self`` with
          ``other`` (or from the product ``other * self`` if ``other`` is a
          scalar)

        TESTS::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.tensor_field(0,2, name='a')
            sage: a[:] = [[1+x, 2], [y, -x^2]]

        Tensor product with another tensor field::

            sage: v = M.vector_field(name='v')
            sage: v[:] = [-y, x]
            sage: s = a.__mul__(v); s
            Tensor field a*v of type (1,2) on the 2-dimensional differentiable
             manifold M
            sage: s.display()
            a*v = -(x + 1)*y d/dx*dx*dx - 2*y d/dx*dx*dy - y^2 d/dx*dy*dx
             + x^2*y d/dx*dy*dy + (x^2 + x) d/dy*dx*dx + 2*x d/dy*dx*dy
             + x*y d/dy*dy*dx - x^3 d/dy*dy*dy
            sage: all([s[ind] == v[ind[0]] * a[ind[1],ind[2]]
            ....:      for ind in M.index_generator(3)])
            True

        Multiplication on the right by a scalar field::

            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: s = a.__mul__(f); s
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: s.display()
            (x^2 + (x + 1)*y + x) dx*dx + (2*x + 2*y) dx*dy + (x*y + y^2) dy*dx
             + (-x^3 - x^2*y) dy*dy
            sage: s == f*a
            True

        """
        # This is to ensure the call to the TensorField version instead of
        # the FreeModuleTensor one
        return TensorField.__mul__(self, other)


    def display_comp(self, frame=None, chart=None, coordinate_labels=True,
                     only_nonzero=True, only_nonredundant=False):
        r"""
        Display the tensor components w.r.t. a given frame, one per
        line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to which
          the tensor field components are defined; if ``None``, then

          - if ``chart`` is not ``None``, the coordinate frame associated to
            ``chart`` is used
          - otherwise, the default basis of the vector field module on which
            the tensor field is defined is used

        - ``chart`` -- (default: ``None``) chart specifying the coordinate
          expression of the components; if ``None``, the default chart of the
          tensor field domain is used
        - ``coordinate_labels`` -- (default: ``True``) boolean; if ``True``,
          coordinate symbols are used by default (instead of integers) as
          index labels whenever ``frame`` is a coordinate frame
        - ``only_nonzero`` -- (default: ``True``) boolean; if ``True``, only
          nonzero components are displayed
        - ``only_nonredundant`` -- (default: ``False``) boolean; if ``True``,
          only nonredundant components are displayed in case of symmetries

        EXAMPLES:

        Display of the components of a type-(2,1) tensor field on a
        2-dimensional manifold::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(2, 1, name='t', sym=(0,1))
            sage: t[0,0,0], t[0,1,0], t[1,1,1] = x+y, x*y, -3
            sage: t.display_comp()
            t^xx_x = x + y
            t^xy_x = x*y
            t^yx_x = x*y
            t^yy_y = -3

        By default, only the non-vanishing components are displayed; to see
        all the components, the argument ``only_nonzero`` must be set to
        ``False``::

            sage: t.display_comp(only_nonzero=False)
            t^xx_x = x + y
            t^xx_y = 0
            t^xy_x = x*y
            t^xy_y = 0
            t^yx_x = x*y
            t^yx_y = 0
            t^yy_x = 0
            t^yy_y = -3

        ``t`` being symmetric w.r.t. to its first two indices, one may ask to
        skip the components that can be deduced by symmetry::

            sage: t.display_comp(only_nonredundant=True)
            t^xx_x = x + y
            t^xy_x = x*y
            t^yy_y = -3

        Instead of coordinate labels, one may ask for integers::

            sage: t.display_comp(coordinate_labels=False)
            t^00_0 = x + y
            t^01_0 = x*y
            t^10_0 = x*y
            t^11_1 = -3

        Display in a frame different from the default one (note that since
        ``f`` is not a coordinate frame, integer are used to label the
        indices)::

            sage: a = M.automorphism_field()
            sage: a[:] = [[1+y^2, 0], [0, 2+x^2]]
            sage: f = X.frame().new_frame(a, 'f')
            sage: t.display_comp(frame=f)
            t^00_0 = (x + y)/(y^2 + 1)
            t^01_0 = x*y/(x^2 + 2)
            t^10_0 = x*y/(x^2 + 2)
            t^11_1 = -3/(x^2 + 2)

        Display with respect to a chart different from the default one::

            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: Y_to_X = X_to_Y.inverse()
            sage: t.display_comp(chart=Y)
            t^uu_u = 1/4*u^2 - 1/4*v^2 + 1/2*u - 3/2
            t^uu_v = 1/4*u^2 - 1/4*v^2 + 1/2*u + 3/2
            t^uv_u = 1/2*u + 3/2
            t^uv_v = 1/2*u - 3/2
            t^vu_u = 1/2*u + 3/2
            t^vu_v = 1/2*u - 3/2
            t^vv_u = -1/4*u^2 + 1/4*v^2 + 1/2*u - 3/2
            t^vv_v = -1/4*u^2 + 1/4*v^2 + 1/2*u + 3/2

        Note that the frame defining the components is the coordinate frame
        associated with chart ``Y``, i.e. we have::

            sage: str(t.display_comp(chart=Y)) == str(t.display_comp(frame=Y.frame(), chart=Y))
            True

        Display of the components with respect to a specific frame, expressed
        in terms of a specific chart::

            sage: t.display_comp(frame=f, chart=Y)
            t^00_0 = 4*u/(u^2 - 2*u*v + v^2 + 4)
            t^01_0 = (u^2 - v^2)/(u^2 + 2*u*v + v^2 + 8)
            t^10_0 = (u^2 - v^2)/(u^2 + 2*u*v + v^2 + 8)
            t^11_1 = -12/(u^2 + 2*u*v + v^2 + 8)

        """
        from sage.misc.latex import latex
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        if frame is None:
            if chart is not None:
                frame = chart.frame()
            else:
                frame = self._fmodule.default_basis()
        if chart is None:
            chart = self._domain.default_chart()
        index_labels = None
        index_latex_labels = None
        if isinstance(frame, CoordFrame) and coordinate_labels:
            ch = frame.chart()
            index_labels = map(str, ch[:])
            index_latex_labels = map(latex, ch[:])
        return FreeModuleTensor.display_comp(self, basis=frame,
                                  format_spec=chart, index_labels=index_labels,
                                  index_latex_labels=index_latex_labels,
                                  only_nonzero=only_nonzero,
                                  only_nonredundant=only_nonredundant)

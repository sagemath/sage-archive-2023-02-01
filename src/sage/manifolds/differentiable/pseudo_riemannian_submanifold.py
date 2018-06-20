r"""
Pseudo-Riemannian submanifolds

An *embedded (resp. immersed) submanifold of a pseudo-Riemannian manifold*
`(M,g)` is an embedded (resp. immersed) submanifold `N` of `M` as a
differentiable manifold such that pull back of the metric tensor `g` via the
embedding (resp. immersion) endows `N` with the structure of a
pseudo-Riemannian manifold.

An actual limitation of this implementation is that a foliation is required to
perform nearly all the calculations (except the induced metric). This is because
the normal vector is easily computed with a foliation, but otherwise requires
some operations which are not yet implemented in Sage (contraction over
different domains).

To correctly compute the normal vector, the submanifold must be declared either
Riemannian or Lorentzian.


The following example explains how to compute the various quantities associated
with the hyperbolic slicing of the 3-dimensional Minkowski space (4D is
similar but 10 time slower to run).


The manifolds must first be declared::

    sage: M = Manifold(3, 'M', structure="Lorentzian")
    sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")

Because the slice considered are spacelike hypersurface, they are Riemannian,
despite being embedded in a Lorentzian manifold.

Let's continue with chart declaration and various free variables::

    sage: E.<w,x,y> = M.chart()
    sage: C.<rh,th> = N.chart(\
    ....:       r'rh:(0,+oo):\rho th:(0,pi):\theta')
    sage: b = var('b',domain='real')
    sage: assume(b>0)
    sage: t = var('t',domain='real')

Here b is the hyperbola semi major axis, and t is the parameter of the
foliation.

One must then define the embedding, as well as the inverse embedding and the
inverse concerning the foliation parameter::

    sage: phi = N.diff_map(M,{(C,E): [b*cosh(rh)+t,
    ....:                             b*sinh(rh)*cos(th),
    ....:                             b*sinh(rh)*sin(th)]})
    sage: phi_inv = M.diff_map(N,{(E,C):[log(sqrt(x**2+y**2+b**2)/b+
    ....:                                sqrt((x**2+y**2+b**2)/b^2-1)),
    ....:                      atan2(y,x)]})
    sage: phi_inv_t = M.scalar_field({E:w-sqrt(x**2+y**2+b**2)})

One can check that the inverse is correct with::

    sage: (phi*phi_inv).display()
    M --> M
       (w, x, y) |--> ((b^2 + x^2 + y^2 + sqrt(b^2 + x^2 + y^2)*(t + sqrt(x^2 +
     y^2)) + sqrt(x^2 + y^2)*t)/(sqrt(b^2 + x^2 + y^2) + sqrt(x^2 + y^2)), x, y)

The first parameter cannot be evaluated yet, because the inverse for t is not
taken into account. To prove that it is correct, one can temporarily inject it
in the result::

    sage: assume(w-t>0)
    sage: (phi*phi_inv).expr()[0].subs({b**2:(w-t)**2-x**2-y**2})\
    ....:           .simplify().expand().simplify_full()
    w
    sage: forget(w-t>0)

The immersion can then be declared::

    sage: N.set_immersion(phi, inverse=phi_inv, var=t,
    ....:                 t_inverse = {t: phi_inv_t})

This line doesn't do any calculation yet. It just check the coherence of the
arguments, but not the inverse, the user is trusted on this point. The user can
also declare that his immersion is in fact an embedding::

    sage: N.declare_embedding()

Don't forget to specify the metric of the Minkowski space::

    sage: g = M.metric('g')
    sage: g[0,0], g[1,1], g[2,2] = -1, 1, 1

With this the declaration of our manifolds is finished, and calculation can be
performed.

The first step is always to find a chart adapted to the foliation. This is done
by the method "adapted_chart"::

    sage: T = N.adapted_chart(); T
    [Chart (M, (rh_M, th_M, t_M))]

T now contains the new charts defined on M. By default, the name of a variable
will be the name of the variable in the submanifold chart indexed by the name
of the manifold.

One can check that the coordinates
changes have been defined correctly::

    sage: len(M.coord_changes())
    2

Let's now compute some quantities:

The induced metric (or first fundamental form)::

    sage: N.induced_metric()[:] # long time
    [           b^2              0]
    [             0 b^2*sinh(rh)^2]

The normal vector::

    sage: N.normal().display() # long time
    n = sqrt(b^2 + x^2 + y^2)/b d/dw + x/b d/dx + y/b d/dy

Check that the hypersurface is indeed spacelike::

    sage: N.ambient_metric()(N.normal(), N.normal()).display() # long time
    g(n,n): M --> R
       (w, x, y) |--> -1
       (rh_M, th_M, t_M) |--> -1

Lapse function::

    sage: N.lapse().display() # long time
    N: M --> R
       (w, x, y) |--> sqrt(b^2 + x^2 + y^2)/b
       (rh_M, th_M, t_M) |--> cosh(rh_M)

Shift vector::

    sage: N.shift().display() # long time
    beta = -(x^2 + y^2)/b^2 d/dw - sqrt(b^2 + x^2 + y^2)*x/b^2 d/dx
     - sqrt(b^2 + x^2 + y^2)*y/b^2 d/dy

The extrinsic curvature (second fundamental form) as a tensor of the ambient
manifold::

    sage: N.ambient_extrinsic_curvature()[:] # long time
    [                                 -(x^2 + y^2)/b^3 (b^2*x + x^3 + x*y^2)/(sqrt(b^2 + x^2 + y^2)*b^3) (y^3 + (b^2 + x^2)*y)/(sqrt(b^2 + x^2 + y^2)*b^3)]
    [                      sqrt(b^2 + x^2 + y^2)*x/b^3                                  -(b^2 + x^2)/b^3                                          -x*y/b^3]
    [                      sqrt(b^2 + x^2 + y^2)*y/b^3                                          -x*y/b^3                                  -(b^2 + y^2)/b^3]

The extrinsic curvature (second fundamental form) as a tensor of the
submanifold (can be quite long because of all the simplifications, about 50
seconds)::

    sage: N.extrinsic_curvature()[:] # long time
    [           -b             0]
    [            0 -b*sinh(rh)^2]



AUTHORS:

- Florentin Jaffredo (2018): initial version

REFERENCES:

- \B. O'Neill : *Semi-Riemannian Geometry* [ONe1983]_
- \J. M. Lee : *Riemannian Manifolds* [Lee1997]_

"""

# *****************************************************************************
#  Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.differentiable.pseudo_riemannian import \
    PseudoRiemannianManifold
from sage.manifolds.differentiable.differentiable_submanifold import \
    DifferentiableSubmanifold
from sage.rings.infinity import infinity
from sage.matrix.constructor import matrix
from sage.functions.other import factorial
from sage.symbolic.ring import SR
from Queue import Queue
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer


class PseudoRiemannianSubmanifold(PseudoRiemannianManifold,
                                  DifferentiableSubmanifold):
    r"""
    Pseudo-Riemannian submanifold.

    An *embedded (resp. immersed) submanifold of a pseudo-Riemannian manifold*
    `(M,g)` is an embedded (resp. immersed) submanifold `N` of `M` as a
    differentiable manifold such that pull back of the metric tensor `g` via
    the embedding (resp. immersion) endows `N` with the structure of a
    pseudo-Riemannian manifold.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``field`` -- field `K` on which the manifold is
      defined; allowed values are

        - ``'real'`` or an object of type ``RealField`` (e.g., ``RR``) for
           a manifold over `\RR`
        - ``'complex'`` or an object of type ``ComplexField`` (e.g., ``CC``)
           for a manifold over `\CC`
        - an object in the category of topological fields (see
          :class:`~sage.categories.fields.Fields` and
          :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
          for other types of manifolds

    - ``structure`` -- manifold structure (see
      :class:`~sage.manifolds.structure.TopologicalStructure` or
      :class:`~sage.manifolds.structure.RealTopologicalStructure`)
    - ``ambient`` -- (default: ``None``) manifold of destination
      of the immersion. If ``None``, set to ``self``
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      topological manifold; the created object is then an open subset of
      ``base_manifold``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none are provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g., coordinates
      in a chart
      - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``Manifolds(field)`` is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    Let `N` be a 2-dimensional submanifold of a 3-dimensional manifold `M`::

        sage: M = Manifold(3, 'M', structure ="pseudo-Riemannian")
        sage: N = Manifold(2, 'N', ambient=M, structure="pseudo-Riemannian")
        sage: N
        2-dimensional pseudo-Riemannian submanifold N embedded in 3-dimensional
         differentiable manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()

    Let us define a 1-dimension foliation indexed by `t`. The inverse map is
    needed in order to compute the adapted chart in the ambient manifold::

        sage: t = var('t')
        sage: phi = N.diff_map(M, {(CN,CM):[u, v, t+u**2+v**2]}); phi
        Differentiable map from the 2-dimensional pseudo-Riemannian submanifold
         N embedded in 3-dimensional differentiable manifold M to the
         3-dimensional Riemannian manifold M
        sage: phi_inv = M.diff_map(N,{(CM, CN): [x,y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x**2-y**2})

    `\phi` can then be declared as an embedding `N\to M`::

        sage: N.set_embedding(phi, inverse=phi_inv, var=t,
        ....:                 t_inverse={t: phi_inv_t})

    The foliation can also be used to find new charts on the ambient manifold
    that are adapted to the foliation, ie in which the expression of the
    immersion is trivial. At the same time, the appropriate coordinate changes
    are computed::

        sage: N.adapted_chart()
        [Chart (M, (u_M, v_M, t_M))]
        sage: len(M.coord_changes())
        2

    .. SEEALSO::

        :mod:`~sage.manifolds.manifold` and
        :mod:`~sage.manifolds.differentiable.differentiable_submanifold`

   """
    def __init__(self, n, name, ambient=None, metric_name='g', signature=None,
                 base_manifold=None, diff_degree=infinity, latex_name=None,
                 metric_latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a pseudo-Riemannian submanifold.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="pseudo-Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="pseudo-Riemannian")
            sage: N
            2-dimensional pseudo-Riemannian submanifold N embedded in
             3-dimensional differentiable manifold M

        """

        PseudoRiemannianManifold.__init__(self, n, name=name,
                                          metric_name=metric_name,
                                          signature=signature,
                                          base_manifold=base_manifold,
                                          diff_degree=diff_degree,
                                          latex_name=latex_name,
                                          metric_latex_name=metric_latex_name,
                                          start_index=start_index,
                                          category=category)
        DifferentiableSubmanifold.__init__(self, n, name, self._field,
                                           self._structure, ambient=ambient,
                                           base_manifold=base_manifold,
                                           latex_name=latex_name,
                                           start_index=start_index,
                                           category=category)

        self._difft = None
        self._gradt = None
        self._normal = None
        self._lapse = None
        self._shift = None
        self._first_fundamental_form = None
        self._ambient_first_fundamental_form = None
        self._second_fundamental_form = None
        self._ambient_second_fundamental_form = None
        self._ambient_metric = None
        self._projector = None
        self._gauss_curvature = None
        self._principal_directions = {}
        self._principal_curvatures = {}
        self._mean_curvature = None
        self._shape_operator = None
        self._sgn = 1 if ambient._structure.name == "Riemannian" else -1

    def _repr_(self):
        r"""
        Return a string representation of the submanifold.

        If no ambient manifold is specified, the submanifold is considered
        as a manifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure="pseudo-Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="pseudo-Riemannian")
            sage: N._repr_()
            '2-dimensional pseudo-Riemannian submanifold N embedded in
             3-dimensional differentiable manifold M'

        """
        if self._ambient is None:
            return super(PseudoRiemannianManifold, self).__repr__()
        return "{}-dimensional pseudo-Riemannian submanifold {} embedded " \
               "in {}-dimensional differentiable " \
               "manifold {}".format(self._dim, self._name, self._ambient._dim,
                                    self._ambient._name)

    @cached_method
    def ambient_metric(self):
        r"""
        Return the metric of the ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the metric of the ambient manifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.ambient_metric()[:]) # long time
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if not self._embedded or not isinstance(self._ambient,
                                                PseudoRiemannianManifold):
            raise ValueError("Submanifold must be "
                             "embedded in a pseudo-Riemannian manifold")
        self._ambient_metric = self._ambient.metric()
        self._ambient_metric.set_name("g", r"g")
        return self._ambient_metric

    @cached_method
    def first_fundamental_form(self):
        r"""
        Return the first fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.induced_metric`
        is an alias.

        OUTPUT:

        - (0,2) tensor field on the submanifold equal to the induced metric

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.first_fundamental_form()[:]) # long time
            [          r^2             0]
            [            0 r^2*sin(th)^2]

        """
        self._first_fundamental_form = self.metric()
        self._first_fundamental_form\
            .set(self._immersion.pullback(self.ambient_metric()))
        self._first_fundamental_form.set_name("gamma", r"\gamma")
        return self._first_fundamental_form

    @cached_method
    def induced_metric(self):
        r"""
        Return the first fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        Alias of
        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.first_fundamental_form`.

        OUTPUT:

        - (0,2) tensor field on the submanifold equal to the induced metric

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.first_fundamental_form()[:]) # long time
            [          r^2             0]
            [            0 r^2*sin(th)^2]

        """
        return self.first_fundamental_form()

    @cached_method
    def difft(self):
        r"""
        Return the differential of the first scalar field defining the
        submanifold

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - 1-form field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.difft().display()) # long time
            dr = x/sqrt(x^2 + y^2 + z^2) dx + y/sqrt(x^2 + y^2 + z^2) dy +
             z/sqrt(x^2 + y^2 + z^2) dz

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to "
                             "perform this calculation")
        self._difft = self._t_inverse[self._var[0]].differential()
        self._difft.set_name("d" + self._var[0]._repr_(),
                             r"\mathrm{d}" + self._var[0]._latex_())
        return self._difft

    @cached_method
    def gradt(self):
        r"""
        Return the gradient of the first scalar field defining the
        submanifold

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - vector field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.gradt().display()) # long time
            grad_r = x/sqrt(x^2 + y^2 + z^2) d/dx + y/sqrt(x^2 + y^2 + z^2) d/dy
             + z/sqrt(x^2 + y^2 + z^2) d/dz

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform "
                             "this calculation")
        self._gradt = self.ambient_metric().inverse()\
            .contract(self.difft())
        self._gradt.set_name("grad_" + self._var[0]._repr_(),
                             r"\nabla " + self._var[0]._latex_())
        return self._gradt

    @cached_method
    def normal(self):
        r"""
        Return a normal unit vector to the submanifold.

        If a foliation is defined, it is used to compute the gradient and then
        the normal vector. If not, the normal vector is computed using the
        following formula:

        .. MATH::

            n = \vec{*}(\mathrm{d}x_0\wedge\mathrm{d}x_1\wedge...
            \wedge\mathrm{d}x_{n-1})

        where the star is the hodge dual operator and de wedge the product on
        the exterior algebra.

        This formula does not always define a proper vector field when multiple
        charts overlap, because of the arbitrariness of the direction of the
        normal vector. To avoid this problem, this function considers the graph
        defined by the atlas of the submanifold and the changes of coordinates,
        and only calculate the normal vector once by connected component. The
        expression is then propagate by restriction, continuation, or change of
        coordinates using a breadth-first exploration of the graph.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - vector field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidian space foliated on the radius::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.normal().display()) # long time
            n = x/sqrt(x^2 + y^2 + z^2) d/dx + y/sqrt(x^2 + y^2 + z^2) d/dy
             + z/sqrt(x^2 + y^2 + z^2) d/dz

        Or in spherical coordinates::

            sage: print(N.normal().display(T[0].frame(),T[0])) # long time
            n = d/dr_M

        The same sphere of constant radius, ie not foliated in stereographic
        coordinates::

            sage: M = Manifold(3, 'M', structure="Riemannian", start_index=1)
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U, V)
            sage: stereoN.<x,y> = U.chart()
            sage: stereoS.<xp,yp> = V.chart("xp:x' yp:y'")
            sage: stereoN_to_S = stereoN.transition_map(stereoS,
            ....:                                 (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                 intersection_name='W',
            ....:                                 restrictions1= x^2+y^2!=0,
            ....:                                 restrictions2= xp^2+yp^2!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: W = U.intersection(V)
            sage: stereoN_W = stereoN.restrict(W)
            sage: stereoS_W = stereoS.restrict(W)
            sage: A = W.open_subset('A', coord_def={stereoN_W: (y!=0, x<0),
            ....:                                  stereoS_W: (yp!=0, xp<0)})
            sage: spher.<the,phi> = A.chart(r'the:(0,pi):\theta phi:(0,2*pi):\phi')
            sage: stereoN_A = stereoN_W.restrict(A)
            sage: spher_to_stereoN = spher.transition_map(stereoN_A,
            ....:                              (sin(the)*cos(phi)/(1-cos(the)),
            ....:                               sin(the)*sin(phi)/(1-cos(the))))
            sage: spher_to_stereoN.set_inverse(2*atan(1/sqrt(x^2+y^2)),
            ....:                                    atan2(-y,-x)+pi)
            sage: stereoN_to_S_A = stereoN_to_S.restrict(A)
            sage: spher_to_stereoS = stereoN_to_S_A * spher_to_stereoN
            sage: stereoS_to_N_A = stereoN_to_S.inverse().restrict(A)
            sage: stereoS_to_spher = spher_to_stereoN.inverse() * stereoS_to_N_A
            sage: E.<X,Y,Z> = M.chart()
            sage: r=1
            sage: phi = N.diff_map(M, {(stereoN, E): [2*x/(1+x^2+y^2),
            ....:                                    2*y/(1+x^2+y^2),
            ....:                                    (x^2+y^2-1)/(1+x^2+y^2)],
            ....:                   (stereoS, E): [2*xp/(1+xp^2+yp^2),
            ....:                                  2*yp/(1+xp^2+yp^2),
            ....:                                 (1-xp^2-yp^2)/(1+xp^2+yp^2)]},
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: N.set_embedding(phi)
            sage: g = M.metric('g')
            sage: g[3,3],g[1,1],g[2,2]=1,1,1
            sage: N.ambient_metric()[:] # long time
            [1 0 0]
            [0 1 0]
            [0 0 1]

        The normal vector is computed the same way, but now returns a
        tensorfield along N::

            sage: n = N.normal() # long time
            sage: print(n) # long time
            Vector field n along the 2-dimensional pseudo-Riemannian submanifold
             N embedded in 3-dimensional differentiable manifold M with values
             on the 3-dimensional Riemannian manifold M

        Let's check that the choice of orientation is coherent on the two top
         frames::

            sage: n.restrict(V).display(format_spec = spher) # long time
            n = -cos(phi)*sin(the) d/dX - sin(phi)*sin(the) d/dY - cos(the) d/dZ
            sage: n.restrict(U).display(format_spec = spher) # long time
            n = -cos(phi)*sin(the) d/dX - sin(phi)*sin(the) d/dY - cos(the) d/dZ

        """
        if self._dim_foliation != 0:    # case foliation
            self._normal = self._sgn*self.lapse()*self.gradt()
            self._normal.set_name("n", r"n")
            return self._normal

        # case no foliation:
        max_frame = self._ambient.default_frame().along(self._immersion)
        self._normal = self.vector_field("n", r"n", self._immersion)

        # an auxiliary functions:
        def calc_normal(chart):
            """
            Calculate the normal vector field according to the formula in the
            documentation in a given chart.
            """
            eps = self.ambient_metric().volume_form(self._dim).along(
                self._immersion).restrict(chart.domain())
            args = list(range(self._dim)) + [eps] + list(range(self._dim))
            r = self.irange()
            n_form = self._immersion.restrict(chart.domain()).pushforward(
                chart.frame()[next(r)]).down(
                self.ambient_metric().along(self._immersion).restrict(
                    chart.domain()))
            for i in r:
                n_form = n_form.wedge(
                    self._immersion.restrict(chart.domain()).pushforward(
                        chart.frame()[i]).down(
                        self.ambient_metric().along(
                            self._immersion).restrict(
                            chart.domain())))
            n_comp = (n_form.contract(*args) / factorial(
                self._dim)).contract(
                self.ambient_metric().inverse().along(self._immersion))
            n_comp = n_comp / n_comp.norm(
                self.ambient_metric().along(self._immersion))

            norm_rst = self._normal.restrict(chart.domain())
            norm_rst.add_comp(max_frame.restrict(chart.domain()))[:] = n_comp[:]
            self._normal.add_comp_by_continuation(max_frame, chart.domain(),
                                                  chart)

        # start breadth-first graph exploration
        marked = set()
        f = Queue()

        for v in self.top_charts():
            if v not in marked:
                f.put(v)
                calc_normal(v)  # initial calculus
                marked.add(v)
                while not f.empty():
                    v = f.get()
                    # for each neighbor:
                    for vp in self.atlas():
                        # case restriction
                        if vp in v._subcharts and vp not in marked:
                            f.put(vp)
                            self._normal.restrict(vp.domain())
                            marked.add(vp)

                        # case continuation
                        if vp in v._supercharts and vp not in marked:
                            f.put(vp)
                            self._normal.add_comp_by_continuation(
                                max_frame.restrict(vp.domain()), v.domain(), vp)
                            marked.add(vp)

                        # case coordinates change
                        if (v, vp) in self.coord_changes() and vp not in marked:
                            f.put(vp)
                            self._normal.comp(max_frame, vp)
                            marked.add(vp)

        # Going up from each top_chart to the full manifold :
        for v in self.top_charts():
            self._normal.add_expr_from_subdomain(max_frame, v.domain())

        return self._normal

    @cached_method
    def ambient_first_fundamental_form(self):
        r"""
        Return the first fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.ambient_induced_metric`
        is an alias.

        OUTPUT:

        - (0,2) tensor field on the ambient manifold describing the induced
          metric before projection on the submanifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.ambient_first_fundamental_form().\
            ....:   display(T[0].frame(),T[0])) # long time
            gamma = r_M^2 dth_M*dth_M + r_M^2*sin(th_M)^2 dph_M*dph_M

        """
        g = self.ambient_metric()

        if self._dim_foliation == 0:  # case no foliation
            g = g.along(self._immersion)

        self._ambient_first_fundamental_form = g - self._sgn * g.contract(
            self.normal()) * g.contract(self.normal())

        self._ambient_first_fundamental_form.set_name("gamma", r"\gamma")
        return self._ambient_first_fundamental_form

    @cached_method
    def ambient_induced_metric(self):
        r"""
        Return the first fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        Alias of
        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.ambient_first_fundamental_form`.

        OUTPUT:

        - (0,2) tensor field on the ambient manifold describing the induced
          metric before projection on the submanifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.ambient_first_fundamental_form().\
            ....:   display(T[0].frame(),T[0])) # long time
            gamma = r_M^2 dth_M*dth_M + r_M^2*sin(th_M)^2 dph_M*dph_M

        """
        return self.ambient_first_fundamental_form()

    @cached_method
    def lapse(self):
        r"""
        Return the lapse function of the foliation

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - scalar field on the ambient manifold equal to the lapse function

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.lapse().display()) # long time
            N: M --> R
               (x, y, z) |--> 1
               (th_M, ph_M, r_M) |--> 1

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        self._lapse = 1 / (self._sgn * self.ambient_metric()(
            self.gradt(), self.gradt())).sqrt()
        self._lapse.set_name("N", r"N")
        return self._lapse

    @cached_method
    def shift(self):
        r"""
        Return the shift function of the foliation

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - shift vector field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.shift().display()) # long time
            beta = 0

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        self._shift = self._adapted_charts[0].frame()[self._dim]\
            - self.lapse() * self.normal()
        self._shift.set_name("beta", r"\beta")
        return self._shift

    @cached_method
    def ambient_second_fundamental_form(self):
        r"""
        Return the second fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.ambient_second_fundamental_form`
        is an alias.

        OUTPUT:

        - (0,2) tensor field on the ambient manifold equal to the extrinsic
          curvature before projection on the submanifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.ambient_second_fundamental_form()\
            ....:       .display(T[0].frame(),T[0])) #long time
            K = -r_M dth_M*dth_M - r_M*sin(th_M)^2 dph_M*dph_M

        """
        if self._dim_foliation == 0:
            self._ambient_second_fundamental_form = \
                           self.tensor_field(0, 2, sym=[(0, 1)], antisym=[],
                                             dest_map=self._immersion)
            k = self.second_fundamental_form()
            g = self.ambient_metric().along(self._immersion)
            max_frame = self._ambient.default_frame().along(self._immersion)
            for chart in self.top_charts():
                pf = [self._immersion.restrict(chart.domain()).pushforward(
                    chart.frame()[i]) for i in self.irange()]
                for i in range(self._dim):
                    pf[i] = pf[i]/g(pf[i], pf[i])
                gam_rst = sum(
                    g.restrict(chart.domain()).contract(pf[i]) *
                    g.restrict(chart.domain()).contract(pf[j]) *
                    self.scalar_field({chart: k.comp(chart.frame())[:][i, j]})
                    for i in range(self._dim) for j in range(self._dim))
                gam_rst._sym = [(0, 1)]
                self._ambient_second_fundamental_form.set_restriction(gam_rst)

            charts = iter(self.top_charts())
            self._ambient_second_fundamental_form.add_comp_by_continuation(
                max_frame, next(charts).domain())
            for chart in charts:
                self._ambient_second_fundamental_form.add_expr_from_subdomain(
                    max_frame, chart.domain())
        else:
            nab = self.ambient_metric().connection('nabla', r'\nabla')
            self._ambient_second_fundamental_form = \
                -self.ambient_metric().contract(nab(self.normal())) \
                - nab(self.normal()).contract(self.normal())\
                .contract(self.ambient_metric())\
                * self.normal().contract(self.ambient_metric())
        self._ambient_second_fundamental_form.set_name("K", r"K")
        return self._ambient_second_fundamental_form

    @cached_method
    def ambient_extrinsic_curvature(self):
        r"""
        Return the second fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        Alias of
        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.ambient_second_fundamental_form`.

        OUTPUT:

        - (0,2) tensor field on the ambient manifold equal to the extrinsic
          curvature before projection on the submanifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.ambient_second_fundamental_form()\
            ....:       .display(T[0].frame(),T[0])) #long time
            K = -r_M dth_M*dth_M - r_M*sin(th_M)^2 dph_M*dph_M

        """
        return self.ambient_second_fundamental_form()

    @cached_method
    def second_fundamental_form(self):
        r"""
        Return the second fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.extrinsic_curvature`.
        is an alias.

        OUTPUT:

        - (0,2) tensor field on the submanifold equal to the extrinsic curvature

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.second_fundamental_form().display()) # long time
            K = -r dth*dth - r*sin(th)^2 dph*dph

        """
        resu = self.vector_field_module() \
            .tensor((0, 2), name='K', latex_name='K', sym=[(0, 1)], antisym=[])
        if self._dim_foliation != 0:
            if self._adapted_charts is None:
                self.adapted_chart()
            inverse_subs = {v: k for k, v in self._subs[0].items()}
            self.ambient_extrinsic_curvature()
            r = list(self._ambient.irange())
            for i in self.irange():
                for j in self.irange():
                    resu[i, j] = self.ambient_extrinsic_curvature()[
                        self._adapted_charts[0].frame(), [r[i], r[j]]].expr(
                        self._adapted_charts[0]).subs(inverse_subs)
        else:
            nab = self.ambient_metric().connection('nabla', r'\nabla')
            n = self.normal()

            for chart in self.atlas():
                gamma_n = matrix(SR, self._dim+1, self._dim+1)
                for i in range(self._dim+1):
                    for j in range(self._dim+1):
                        gamma_n[i, j] = sum(
                            nab[self._ambient.frames()[0], :][i][j][k].expr() *
                            n.restrict(chart.domain()).comp(
                                n.restrict(chart.domain())._fmodule.bases()[0])
                            [:][k].expr() for k in
                            range(self._dim + 1))
                dXdu = self._immersion.differential_functions(chart)
                dNdu = matrix(SR, self._dim+1, self._dim)
                for i in range(self._dim+1):
                    for j in range(self._dim):
                        dNdu[i, j] = n.restrict(chart.domain()).comp(
                            n.restrict(chart.domain())._fmodule.bases()[0])[:,
                            chart][i].diff(chart[:][j]).expr()
                g = self.ambient_metric().along(
                    self._immersion.restrict(chart.domain())).restrict(
                    chart.domain())[:, chart]
                K = dXdu.transpose()*g*(dNdu+gamma_n*dXdu)
                resu[chart.frame(), :] = K

        self._second_fundamental_form = resu
        return self._second_fundamental_form

    @cached_method
    def extrinsic_curvature(self):
        r"""
        Return the second fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        Alias of
        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.second_fundamental_form`.

        OUTPUT:

        - (0,2) tensor field on the submanifold equal to the extrinsic curvature

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.second_fundamental_form().display()) # long time
            K = -r dth*dth - r*sin(th)^2 dph*dph

        """
        return self.second_fundamental_form()

    @cached_method
    def projector(self):
        r"""
        Return the projector on the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the (1,1) tensor field on the ambient manifold corresponding to the
          projector along the normal of the submanifold

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1

        Print the projector::

            sage: print(N.projector()) # long time
            Tensor field gamma of type (1,1) on the 3-dimensional Riemannian
             manifold M

        Check that the projector applied to the normal vector is zero::

            sage: N.projector().contract(N.normal()).display() # long time
            0

        """

        g = self.ambient_metric().inverse()
        if self._dim_foliation == 0:
            g = g.along(self._immersion)

        self._projector = self.ambient_first_fundamental_form().contract(0, g)
        self._projector.set_name("gamma", r"\vec{\gamma}")
        return self._projector

    def project(self, tensor):
        r"""
        Return the projection of a tensor on the submanifold.

        INPUT:

        - ``tensor`` -- any tensor field to be projected on the manifold. If no
          foliation is provided, must be a tensorfield along the submanifold.

        OUTPUT:

        - tensor field on the ambient manifold, projection of the input tensor
          along the normal of the submanifold.

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi,inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1

        Let's check that the first fundamental form is the projection of the
        metric on the submanifold::

            sage: N.ambient_first_fundamental_form()==N.project(M.metric()) # long time
            True

        The result is not cached :

            sage: N.ambient_first_fundamental_form() is N.project(M.metric()) # long time
            False

        """
        resu = tensor.copy()
        resu.set_name(tensor._name + "_" + self._name,
                      r"{" + tensor._latex_() + r"}_{" + self._latex_() + r"}")
        for i in range(tensor.tensor_type()[0]):
            resu = self.projector().contract(1, resu, i)
        for i in range(tensor.tensor_type()[1]):
            resu = self.projector().contract(0, resu, i)
        return resu

    def mixed_projection(self, tensor, indices=0):
        r"""
        Return de n+1 decomposition of a tensor on the submanifold and the
        normal vector.

        The n+1 decomposition of a tensor of rank `k` can be obtained by
        contracting each indices either with the normal vector or the projection
        operator of the submanifold (see
        :meth:`sage.manifold.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.projector`).

        INPUT:

        - ``tensor`` -- any tensor field, eventually along the submanifold if
          no foliation is provided.
        - ``indices`` -- (default: ``0``) list of integers containing the
          indices on which the projection is made on the normal vector.
          By default, all projection are made on the submanifold. If
          an integer `n` is provided, the `n` first contractions are made with
          the normal vector, all the other with the projection operator.

        OUTPUT:

        - tensorfield of rank `k`-``len(indices)``.

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi,inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1

        Let's check that the first fundamental form is the projection of the
        metric on the submanifold::

            sage: N.ambient_first_fundamental_form()==N.mixed_projection(M.metric()) # long time
            True

        The other non redundant projection are::

            sage: N.mixed_projection(M.metric(), [0]) # long time
            1-form on the 3-dimensional Riemannian manifold M

        and::
            sage: N.mixed_projection(M.metric(), [0,1]) # long time
            Scalar field on the 3-dimensional Riemannian manifold M

        Which is constant and equal to 1::

            sage: N.mixed_projection(M.metric(), [0,1]).display() # long time
            M --> R
            (x, y, z) |--> 1
            (th_M, ph_M, r_M) |--> 1

        """
        if isinstance(indices, (Integer, int)):
            indices = range(indices)

        if len(indices)>tensor.tensor_rank():
            raise ValueError("Too much contractions")

        g = self.ambient_metric()
        if self._dim_foliation == 0:
            g = g.along(self._immersion)

        multiprojector = 1
        k = tensor.tensor_rank()    # order of the tensor
        kp = 2*k-len(indices)       # order of the multiprojector
        for i in range(tensor.tensor_type()[1]):
            if i in indices:
                multiprojector = multiprojector * self.normal()
            else:
                multiprojector = multiprojector * self.projector()
        for i in range(tensor.tensor_type()[0]):
            if i in indices:
                multiprojector = multiprojector * self.normal().contract(g)
            else:
                multiprojector = multiprojector * self.projector()
        args = range(kp - tensor.tensor_type()[0], kp) + range(
            tensor.tensor_type()[1]) + [tensor] + range(k)
        return multiprojector.contract(*args)

    @cached_method
    def gauss_curvature(self):
        r"""
        Return the gauss curvature of the submanifold.

        The Gauss curvature is the product or the principal curvatures, or
        equivalently the determinant of the projection operator.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - scalar field on the submanifold equal to the Gauss curvature

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: N.gauss_curvature().display() # long time
            N --> R
            (th, ph) |--> r^(-2)

        """
        a = self.shape_operator()
        self._gauss_curvature = self.scalar_field(
            {chart: a[chart.frame(), :, chart].determinant()
             for chart in self.top_charts()})
        return self._gauss_curvature

    @cached_method
    def principal_directions(self, chart):
        r"""
        Return the principal directions of the submanifold.

        The principal directions are the eigenvectors of the projection
        operator. The result is formatted as a list of couples
        (eigenvector, eigenvalue).

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``chart`` --  chart in which the principal directions are to be
          computed

        OUTPUT:

        - List of couples  (vectorfields,scalarfield) representing the
          principal directions and the principal curvature associated

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.principal_directions(C)[0][0].display()) # long time
            e_0 = d/dth

        """
        a = self.shape_operator()
        pr_d = matrix(
            [[a[chart.frame(), :, chart][i, j].expr() for i in self.irange()]
             for j in self.irange()]).eigenvectors_right()
        res = []
        v = self.vector_field()
        counter = self.irange()
        for eigen_space in pr_d:
            for eigen_vector in eigen_space[1]:
                v[chart.frame(), :] = eigen_vector
                res.append((v.copy(), eigen_space[0]))
                res[-1][0].set_name("e_{}".format(next(counter)))
        self._principal_directions[chart] = res
        return res

    @cached_method
    def principal_curvatures(self, chart):
        r"""
        Return the principal curvatures of the submanifold.

        The principal curvatures are the eigenvalues of the projection operator.
        The resulting scalarfield are named "k_i" with i index of the
        submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``chart`` --  chart in which the principal curvatures are to be
          computed

        OUTPUT:

        - list of scalar field on the submanifold equal to the principal
          curvatures

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.principal_curvatures(C)[0].display()) # long time
            k_0: N --> R
               (th, ph) |--> -1/r
        """
        a = self.shape_operator()
        res = matrix(
            [[a[chart.frame(), :, chart][i, j].expr() for i in self.irange()]
             for j in self.irange()]).eigenvalues()
        counter = self.irange()
        for i in range(self._dim):
            res[i] = self.scalar_field({chart: res[i]},
                                       name="k_{}".format(next(counter)))
        self._principal_curvatures[chart] = res
        return res

    @cached_method
    def mean_curvature(self):
        r"""
        Return the mean curvature of the submanifold.

        The mean curvature is the arithmetic mean of the principal curvatures,
        or equivalently the trace of the projection operator.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - scalar field on the submanifold equal to the mean curvature

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.mean_curvature().display()) # long time
            N --> R
            (th, ph) |--> -1/r


        """
        self._shape_operator = self.scalar_field({chart: self._sgn * sum(
            self.principal_curvatures(chart)).expr(chart) / self._dim
                                                  for chart in
                                                  self.top_charts()})
        return self._shape_operator

    @cached_method
    def shape_operator(self):
        r"""
        Return the shape operator of the submanifold.

        The shape operator is equal to the second fundamental form with one of
        the indices upped.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - (1,1) tensor field on the submanifold equal to the shape operator

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: print(N.shape_operator()[:]) # long time
            [-1/r    0]
            [   0 -1/r]

        """
        self._shape_operator = self.second_fundamental_form().contract(
            self.induced_metric().inverse())
        return self._shape_operator

    def clear_cache(self):
        """
        Reset all the cached functions and the derived quantities.

        Use this function if you modified the immersion (or embedding) of the
        submanifold. Note that when calling calculus function after clearing,
        new pythons objects will be created.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E): [r*sin(th)*cos(ph),
            ....:                             r*sin(th)*sin(ph),
            ....:                             r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_embedding(phi, inverse = phi_inv, var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[2,2] = 1, 1, 1
            sage: n = N.normal()
            sage: n is N.normal()
            True
            sage: N.clear_cache()
            sage: n is N.normal()
            False

        """
        self.ambient_metric.clear_cache()
        self.first_fundamental_form.clear_cache()
        self.induced_metric.clear_cache()
        self.difft.clear_cache()
        self.gradt.clear_cache()
        self.normal.clear_cache()
        self.ambient_first_fundamental_form.clear_cache()
        self.ambient_induced_metric.clear_cache()
        self.lapse.clear_cache()
        self.shift.clear_cache()
        self.ambient_second_fundamental_form.clear_cache()
        self.ambient_extrinsic_curvature.clear_cache()
        self.second_fundamental_form.clear_cache()
        self.extrinsic_curvature.clear_cache()
        self.projector.clear_cache()
        self.gauss_curvature.clear_cache()
        self.principal_directions.clear_cache()
        self.principal_curvatures.clear_cache()
        self.shape_operator.clear_cache()
        self._difft = None
        self._gradt = None
        self._normal = None
        self._lapse = None
        self._shift = None
        self._first_fundamental_form = None
        self._ambient_first_fundamental_form = None
        self._second_fundamental_form = None
        self._ambient_second_fundamental_form = None
        self._ambient_metric = None
        self._projector = None
        self._gauss_curvature = None
        self._principal_directions = {}
        self._principal_curvatures = {}
        self._mean_curvature = None
        self._shape_operator = None


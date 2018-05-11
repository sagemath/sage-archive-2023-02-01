r"""
Pseudo-Riemannian submanifold of a differentiable manifold

A pseudo-Riemannian submanifold of a differentiable manifold is a differentiable
submanifold which is also pseudo-Riemannian.

An actual limitation of this implementation is that a foliation is required to
perform nearly all the calculations (except the induced metric). This is because
the normal vector is easily computed with a foliation, but otherwise requires
some operations which are not yet implemented in Sage (contraction over
different domains).

To correctly compute the normal vector, the submanifold must be declared either
Riemannian or Lorentzian.


The following example explains how to compute the various quantities associated
with the hyperbolic slicing of the Minkowski space.

The manifolds must first be declared::

    sage: M = Manifold(4,'M',structure = "Lorentzian")
    sage: N = Manifold(3,'N',ambient = M,structure = "Riemannian")

Because the slice considered are spacelike hypersurface, they are Riemannian,
despite being embedded in a Lorentzian manifold.

Let's continue with chart declaration and various free variables::

    sage: E.<w,x,y,z> = M.chart()
    sage: C.<rh,th,ph> = N.chart(r'rh:(0,+oo):\rho th:(0,pi):\theta ph:(-pi,pi):\phi')
    sage: b = var('b',domain = 'real')
    sage: assume(b>0)
    sage: t = var('t',domain = 'real')

Here b is the hyperbola semi major axis, and t is the parameter of the
foliation.

One must then define the embedding, as well as the inverse embedding and the
inverse concerning the foliation parameter::

    sage: phi = N.diff_map(M,{(C,E):[b*cosh(rh)+t,
    ....:                            b*sinh(rh)*sin(th)*cos(ph),
    ....:                            b*sinh(rh)*sin(th)*sin(ph),
    ....:                            b*sinh(rh)*cos(th)]})
    sage: phi_inv = M.diff_map(N,{(E,C):[log(sqrt(x**2+y**2+z**2+b**2)\
    ....:                      /b+sqrt((x**2+y**2+z**2+b**2)/b^2-1)),
    ....:                      arccos(z/(b*sqrt((x**2+y**2+z**2+b**2)/b^2-1))),
    ....:                      atan2(y,x)]})
    sage: phi_inv_t = M.scalar_field({E:w-sqrt(x**2+y**2+z**2+b**2)})

One can check that the inverse is correct with::

    sage: (phi*phi_inv).display()
    M --> M
    (w, x, y, z) |--> ((b^2 + x^2 + y^2 + z^2 + sqrt(b^2 + x^2 + y^2 + z^2)*(t + sqrt(x^2 + y^2 + z^2)) + sqrt(x^2 + y^2 + z^2)*t)/(sqrt(b^2 + x^2 + y^2 + z^2) + sqrt(x^2 + y^2 + z^2)), x, y, z)

The first parameter cannot be evaluated yet, because the inverse for t is not
taken into account. To prove that it is correct, one can temporarily inject it
in the result::

    sage: assume(w-t>0)
    sage: (phi*phi_inv).expr()[0].subs({b**2:(w-t)**2-x**2-y**2-z**2})\
    ....:           .simplify().expand().simplify_full()
    w
    sage: forget(w-t>0)

The immersion can then be declared::

    sage: N.set_immersion(phi,phi_inverse = phi_inv,var = t,t_inverse = {t:phi_inv_t})

This line doesn't do any calculation yet. It just check the coherence of the
arguments, but not the inverse, the user is trusted on this point. The user can
also declare that his immersion is in fact an embedding::

    sage: N.declare_embedding()

Don't forget to specify the metric of the Minkowski space::

    sage: g = M.metric('g')
    sage: g[0,0], g[1,1], g[2,2], g[3,3]=-1, 1, 1, 1

With this the declaration of our manifolds is finished, and calculation can be
performed.

The first step is always to find a chart adapted to the foliation. This is done
by the method "adapted_chart"::

    sage: T = N.adapted_chart(); T
    [Chart (M, (rh_M, th_M, ph_M, t_M))]

T now contains the new charts defined on M. By default, the name of a variable
will be the name of the variable in the submanifold chart indexed by the name
of the manifold.

One can check that the coordinates
changes have been defined correctly::

    sage: len(M.coord_changes())
    2

Let's now compute some quantities:

The induced metric (or first fundamental form)::

    sage: N.induced_metric()[:]
    [                     b^2                        0                        0]
    [                       0           b^2*sinh(rh)^2                        0]
    [                       0                        0 b^2*sin(th)^2*sinh(rh)^2]

The normal vector::

    sage: N.normal().display()
    sqrt(b^2 + x^2 + y^2 + z^2)/b d/dw + x/b d/dx + y/b d/dy + z/b d/dz

Check that the hypersurface is indeed spacelike::

    sage: N.ambient_metric()(N.normal(),N.normal()).display()
    M --> R
    (w, x, y, z) |--> -1
    (rh_M, th_M, ph_M, t_M) |--> -1

Lapse function::

    sage: N.lapse().display()
    N: M --> R
       (w, x, y, z) |--> sqrt(b^2 + x^2 + y^2 + z^2)/b
       (rh_M, th_M, ph_M, t_M) |--> cosh(rh_M)

Shift vector::

    sage: N.shift().display()
    beta = -(x^2 + y^2 + z^2)/b^2 d/dw - sqrt(b^2 + x^2 + y^2 + z^2)*x/b^2 d/dx
     - sqrt(b^2 + x^2 + y^2 + z^2)*y/b^2 d/dy - sqrt(b^2 + x^2 + y^2 + z^2)*z/b^2 d/dz

The extrinsic curvature (second fundamental form) as a tensor of the ambient
manifold::

    sage: N.ambient_extrinsic_curvature()[:]
    [           -(x^2 + y^2 + z^2)/b^3 sqrt(b^2 + x^2 + y^2 + z^2)*x/b^3 sqrt(b^2 + x^2 + y^2 + z^2)*y/b^3 sqrt(b^2 + x^2 + y^2 + z^2)*z/b^3]
    [sqrt(b^2 + x^2 + y^2 + z^2)*x/b^3                  -(b^2 + x^2)/b^3                          -x*y/b^3                          -x*z/b^3]
    [sqrt(b^2 + x^2 + y^2 + z^2)*y/b^3                          -x*y/b^3                  -(b^2 + y^2)/b^3                          -y*z/b^3]
    [sqrt(b^2 + x^2 + y^2 + z^2)*z/b^3                          -x*z/b^3                          -y*z/b^3                  -(b^2 + z^2)/b^3]

The extrinsic curvature (second fundamental form) as a tensor of the
submanifold (can be quite long because of all the simplifications, about 50
seconds)::

    sage: N.extrinsic_curvature()[:] # long time
    [                     -b                       0                       0]
    [                      0           -b*sinh(rh)^2                       0]
    [                      0                       0 -b*sin(th)^2*sinh(rh)^2]



AUTHORS:

- Florentin Jaffredo

"""

# *****************************************************************************
#   Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.differentiable.pseudo_riemannian import \
    PseudoRiemannianManifold
from sage.manifolds.submanifold.differentiable_submanifold import \
    DifferentiableSubmanifold
from sage.rings.infinity import infinity
from sage.matrix.constructor import matrix
from sage.functions.other import factorial
from sage.symbolic.ring import SR
from Queue import Queue

class PseudoRiemannianSubmanifold(PseudoRiemannianManifold,
                                  DifferentiableSubmanifold):
    r"""
    Pseudo-Riemannian submanifold of a differentiable manifold

    A pseudo-Riemannian submanifold of a differentiable manifold is a
    differentiable submanifold which is also pseudo-Riemannian.

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

    Let N be a 2-dimensional submanifold of M, 3-dimensional manifold::

        sage: M = Manifold(3, 'M', structure ="pseudo-Riemannian")
        sage: N = Manifold(2, 'N', ambient = M, structure ="pseudo-Riemannian")
        sage: N
        2-dimensional pseudo-Riemannian submanifold N embedded in 3-dimensional
         differentiable manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()

    Let's define a 1-dimension foliation indexed by t. The inverse map is needed
    in order to compute the adapted chart in the ambient manifold::

        sage: t = var('t')
        sage: phi = N.diff_map(M, {(CN,CM):[u, v, t+u**2+v**2]}); phi
        Differentiable map from the 2-dimensional pseudo-Riemannian submanifold
         N embedded in 3-dimensional differentiable manifold M to the
         3-dimensional Riemannian manifold M
        sage: phi_inv = M.diff_map(N,{(CM, CN): [x,y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x**2-y**2})

    \phi can then be declared as an embedding from N to M::

        sage: N.set_immersion(phi, phi_inverse = phi_inv, var = t,\
        ....:                 t_inverse = {t: phi_inv_t})
        sage: N.declare_embedding()

    The foliation can also be used to find new charts on the ambient manifold
    that are adapted to the foliation, ie in which the expression of the
    immersion is trivial. At the same time coordinates changes or computed::

        sage: N.adapted_chart()
        [Chart (M, (u_M, v_M, t_M))]
        sage: len(M._coord_changes)
        2

    .. SEEALSO::

        :mod:`sage.manifolds.manifold`
        :mod:`sage.manifolds.submanifold.differentiable_submanifold`
   """
    def __init__(self, n, name, ambient=None, metric_name='g', signature=None,
                 base_manifold=None, diff_degree=infinity, latex_name=None,
                 metric_latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct an immersion of a given differentiable manifold.

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
        self._principal_directions = None
        self._principal_curvatures = None
        self._shape_operator = None
        self._sgn = 1 if ambient._structure.name == "Riemannian" else -1

    def _repr_(self):
        r"""
        Return a string representation of the submanifold. If no ambient
        manifold is specified, the submanifold is considered as a manifold

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

    def ambient_metric(self, recache=False):
        r"""
        Return the metric of the ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - the metric of the ambient manifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.ambient_metric()[:])
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if not self._embedded or not isinstance(self._ambient,
                                                PseudoRiemannianManifold):
            raise ValueError("Submanifold must be "
                             "embedded in a pseudo-Riemannian manifold")
        if self._ambient_metric is not None and not recache:
            return self._ambient_metric
        self._ambient_metric = self._ambient.metric()
        self._ambient_metric.set_name("g", r"g")
        return self._ambient_metric

    def first_fundamental_form(self, recache=False):
        r"""
        Return the first fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - (0,2) tensor field on the submanifold equal to the induced metric

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.first_fundamental_form()[:])
            [          r^2             0]
            [            0 r^2*sin(th)^2]

        """
        if self._first_fundamental_form is not None and not recache:
            return self._first_fundamental_form
        self._first_fundamental_form = self.metric()
        self._first_fundamental_form\
            .set(self._immersion.pullback(self.ambient_metric(recache)))
        self._first_fundamental_form.set_name("gamma", r"\gamma")
        return self._first_fundamental_form

    induced_metric = first_fundamental_form

    def difft(self, recache=False):
        r"""
        Return the differential of the first scalar field defining the
        submanifold

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - 1-form field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.difft().display())
            dr = x/sqrt(x^2 + y^2 + z^2) dx + y/sqrt(x^2 + y^2 + z^2) dy +
             z/sqrt(x^2 + y^2 + z^2) dz

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to "
                             "perform this calculation")
        if self._difft is not None and not recache:
            return self._difft
        self._difft = self._t_inverse[self._var[0]].differential()
        self._difft.set_name("d" + self._var[0]._repr_(),
                             r"\mathrm{d}" + self._var[0]._latex_())
        return self._difft

    def gradt(self, recache=False):
        r"""
        Return the gradient of the first scalar field defining the
        submanifold

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - vector field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.gradt().display())
            grad_r = x/sqrt(x^2 + y^2 + z^2) d/dx + y/sqrt(x^2 + y^2 + z^2) d/dy
             + z/sqrt(x^2 + y^2 + z^2) d/dz

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform "
                             "this calculation")
        if self._gradt is not None and not recache:
            return self._gradt
        self._gradt = self.ambient_metric(recache).inverse()\
            .contract(self.difft(recache))
        self._gradt.set_name("grad_" + self._var[0]._repr_(),
                             r"\nabla " + self._var[0]._latex_())
        return self._gradt

    def normal(self, recache=False):
        r"""
        Return a normal unit vector to the submanifold.

        If a foliation is defined, it is used to compute the gradient and then
        the normal vector. If not, the normal vector is computed using the
        following formula:

        .. MATH::

            n = \overrightarrow{*}(\mathrm{d}x_0\wedge\mathrm{d}x_1\wedge...
            \wedge\mathrm{d}x_{n-1})

        where the star is the hodge dual operator and de wedge the product on
        the exterior algebra.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - vector field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.normal().display())
            x/sqrt(x^2 + y^2 + z^2) d/dx + y/sqrt(x^2 + y^2 + z^2) d/dy
             + z/sqrt(x^2 + y^2 + z^2) d/dz

        Or in spherical coordinates::

            sage: print(N.normal().display(T[0].frame(),T[0]))
            d/dr_M


        """
        if self._normal is not None and not recache:
            return self._normal
        if self._dim_foliation != 0:    # case foliation
            self._normal = self._sgn*self.lapse(recache)*self.gradt(recache)
            self._normal.set_name("n", r"n")
            return self._normal
        else:                           # case no foliation
            self._normal = self._ambient.vector_field().along(self._immersion)
            self._normal.set_name("n", r"n")

            for chart in self.atlas():
                eps = self.ambient_metric(recache).volume_form(self._dim).along(
                    self._immersion).restrict(chart.domain())
                args = list(range(self._dim)) + [eps] + list(range(self._dim))
                r = self.irange()
                n_form = self._immersion.restrict(chart.domain()).pushforward(
                    chart.frame()[r.next()]).down(
                    self.ambient_metric().along(self._immersion).restrict(
                        chart.domain()))
                for i in r:
                    n_form = n_form.wedge(
                        self._immersion.restrict(chart.domain()).pushforward(
                            chart.frame()[i]).down(
                            self.ambient_metric().along(self._immersion).restrict(
                                chart.domain())))
                n_comp = (n_form.contract(*args) / factorial(
                        self._dim)).contract(
                        self.ambient_metric().inverse().along(self._immersion))
                n_comp = n_comp / n_comp.norm(
                    self.ambient_metric().along(self._immersion))

                frame = next(iter(n_comp._components))
                self._normal.add_comp(frame)[:] = n_comp[frame,:]

            return self._normal

    def normal_by_breadth_first_search(self,recache=False):
        if self._normal is not None and not recache:
            return self._normal
        if self._dim_foliation != 0:
            return self.normal(recache)
        else:

            def calc_normal(chart):
                eps = self.ambient_metric(recache).volume_form(self._dim).along(
                    self._immersion).restrict(chart.domain())
                args = list(range(self._dim)) + [eps] + list(range(self._dim))
                r = self.irange()
                n_form = self._immersion.restrict(chart.domain()).pushforward(
                    chart.frame()[r.next()]).down(
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

                frame = next(iter(n_comp._components))
                self._normal.add_comp(frame)[:] = n_comp[frame, :]

            def calc_by_restrict(chart1, chart2):
                frame = next(iter(self._normal.restrict(chart1.domain())._components))
                self._normal.add_comp(frame.restrict(chart2.domain()))[:] = self._normal[frame, :]

            self._normal = self._ambient.vector_field().along(self._immersion)
            self._normal.set_name("n", r"n")

            G = self.atlas()
            marked = set()
            f = Queue()

            for v in G:
                if v not in marked:
                    f.put(v)
                    calc_normal(v) # first calculus
                    marked.add(v)
                    while not f.empty():
                        v = f.get()
                        # for each neighbor:
                        print("first chart : ",v)
                        for vp in G:
                            if vp in v._subcharts and vp not in marked:
                                print("    second chart : "+ vp._repr_()+" (restriction)")
                                f.put(vp)
                                calc_by_restrict(v,vp)
                                marked.add(vp)
                            if vp in v._supercharts and vp not in marked:
                                print("    second chart : "+ vp._repr_()+ " (continuation)")
                                f.put(vp)
                                self._normal.add_comp_by_continuation(vp.frame(),v.domain())
                                marked.add(vp)
                            if (v,vp) in self.coord_changes() and vp not in marked:
                                print("    second chart : "+ vp._repr_()+ " (coord change)")
                                f.put(vp)
                                calc_normal(vp)
                                if self._normal.restrict(v.domain())[vp,:] == -self._normal.restrict(vp.domain)[vp,:]:
                                    frame = next(iter(self._normal.restrict(vp.domain)._components))
                                    self._normal.restrict(v.domain())._component[frame] = -self._normal._component[frame]
                                marked.add(vp)
            return self._normal

    def ambient_first_fundamental_form(self, recache=False):
        r"""
        Return the first fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - (0,2) tensor field on the ambient manifold describing the induced
          metric before projection on the submanifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.ambient_first_fundamental_form().\
            ....:   display(T[0].frame(),T[0]))
            gamma = r_M^2 dth_M*dth_M + r_M^2*sin(th_M)^2 dph_M*dph_M

        """
        if self._ambient_first_fundamental_form is not None and not recache:
            return self._ambient_first_fundamental_form
        self.ambient_metric(recache)
        self.normal(recache)
        self._ambient_first_fundamental_form = \
            self.ambient_metric() - self._sgn\
            * self.ambient_metric().contract(self.normal())\
            * self.ambient_metric().contract(self.normal())
        self._ambient_first_fundamental_form.set_name("gamma", r"\gamma")
        return self._ambient_first_fundamental_form

    ambient_induced_metric = ambient_first_fundamental_form

    def lapse(self, recache=False):
        r"""
        Return the lapse function of the foliation

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - scalar field on the ambient manifold equal to the lapse function

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.lapse().display())
            N: M --> R
               (x, y, z) |--> 1
               (th_M, ph_M, r_M) |--> 1

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        if self._lapse is not None and not recache:
            return self._lapse
        self._lapse = 1 / (self._sgn * self.ambient_metric(recache)(
            self.gradt(recache), self.gradt(recache))).sqrt()
        self._lapse.set_name("N", r"N")
        return self._lapse

    def shift(self, recache=False):
        r"""
        Return the shift function of the foliation

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - shift vector field on the ambient manifold.

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.shift().display())
            beta = 0

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        if self._shift is not None and not recache:
            return self._shift
        self._shift = self._adapted_charts[0].frame()[self._dim]\
            - self.lapse(recache) * self.normal(recache)
        self._shift.set_name("beta", r"\beta")
        return self._shift

    def ambient_second_fundamental_form(self, recache=False):
        r"""
        Return the second fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - (0,2) tensor field on the ambient manifold equal to the extrinsic
          curvature before projection on the submanifold

        EXAMPLES:

        A sphere embedded in euclidian space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.ambient_second_fundamental_form()\
            ....:       .display(T[0].frame(),T[0])) #long time
            K = -r_M dth_M*dth_M - r_M*sin(th_M)^2 dph_M*dph_M

        """
        if self._ambient_second_fundamental_form is not None and not recache:
            return self._ambient_second_fundamental_form
        nab = self.ambient_metric().connection('nabla', r'\nabla')
        self.normal(recache)
        if self._dim_foliation == 0:
            # g = self.ambient_metric(recache).along(self._immersion)
            # self._ambient_second_fundamental_form = -g.contract(nab(self.normal())) - nab(self.normal()).contract(self.normal()).contract(g) * self.normal().contract(g)
            pass
        else:
            self._ambient_second_fundamental_form = \
                -self.ambient_metric().contract(nab(self.normal())) \
                - nab(self.normal()).contract(self.normal())\
                .contract(self.ambient_metric())\
                * self.normal().contract(self.ambient_metric())
        self._ambient_second_fundamental_form.set_name("K", r"K")
        return self._ambient_second_fundamental_form

    ambient_extrinsic_curvature = ambient_second_fundamental_form

    def second_fundamental_form(self, recache=False):
        r"""
        Return the second fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - (0,2) tensor field on the submanifold equal to the extrinsic curvature

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.second_fundamental_form().display()) # long time
            K = -r dth*dth - r*sin(th)^2 dph*dph

        """
        if self._second_fundamental_form is not None and not recache:
            return self._second_fundamental_form
        resu = self._immersion._domain.vector_field_module() \
            .tensor((0, 2), name='K', latex_name='K', sym=[(0, 1)], antisym=[])
        if self._dim_foliation != 0:
            inverse_subs = {v: k for k, v in self._subs[0].items()}
            self.ambient_extrinsic_curvature(recache)
            r = list(self._ambient.irange())
            for i in self.irange():
                for j in self.irange():
                    resu[i, j] = self.ambient_extrinsic_curvature()[
                        self._adapted_charts[0].frame(), [r[i], r[j]]].expr(
                        self._adapted_charts[0]).subs(inverse_subs)
        else:
            nab = self.ambient_metric().connection('nabla', r'\nabla')
            n = self.normal(recache)

            gamma_n = matrix(self._dim+1,self._dim+1)
            for i in range(self._dim+1):
                for j in range(self._dim+1):
                    gamma_n[i, j] = sum(
                        nab[self._ambient.frames()[0], :][i][j][k].expr() *
                        n.comp(n._fmodule.bases()[0])[:][k].expr() for k in
                        range(self._dim + 1))
            dXdu = self._immersion.differential(
                self(self.atlas()[0][:])).matrix()
            dNdu = matrix(SR,self._dim+1,self._dim)
            for i in range(self._dim+1):
                for j in range(self._dim):
                    dNdu[i,j] = n.comp(n._fmodule.bases()[0])[:][i].diff(self.atlas()[0][:][j]).expr()
            g = self.ambient_metric().along(self._immersion)[:]
            K = -dXdu.transpose()*g*(dNdu+gamma_n*dXdu)
            resu[self.atlas()[0].frame(),:] = K
        self._second_fundamental_form = resu
        return self._second_fundamental_form

    extrinsic_curvature = second_fundamental_form

    def projector(self, recache=False):
        r"""
        Return the projector on the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - the (1,1) tensor field on the ambient manifold corresponding to the
          projector along the normal of the submanifold

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1

        Print the projector::

            sage: print(N.projector()) # long time
            Tensor field gamma of type (1,1) on the 3-dimensional Riemannian
             manifold M

        Check that the projector applied to the normal vector is zero::

            sage: N.projector().contract(N.normal()).display()
            0

        """
        if self._projector is not None and not recache:
            return self._projector
        self._projector = self.ambient_first_fundamental_form(recache).up(
            self.ambient_metric(recache), pos=0)
        self._projector.set_name("gamma", r"\overrightarrow{\gamma}")
        return self._projector

    def project(self, tensor):
        r"""
        Return the projection of a tensor on the submanifold.

        INPUT:

        - ``tensor`` -- any tensor field to be projected on the manifold

        OUTPUT:

        - tensor field on the ambient manifold, projection of the input tensor
          along the normal of the submanifold.

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1

        Let's check that the first fundamental form is the projection of the
        metric on the submanifold::

            sage: N.ambient_first_fundamental_form()==N.project(M.metric())
            True

        The result is not cached :

            sage: N.ambient_first_fundamental_form() is N.project(M.metric())
            False

        """
        resu = tensor.copy()
        resu.set_name(tensor._name + "_" + self._name,
                      r"{" + tensor._latex_() + r"}_" + self._latex_())
        for i in range(tensor.tensor_type()[0]):
            resu = self.projector().contract(1, resu, i)
        for i in range(tensor.tensor_type()[1]):
            resu = self.projector().contract(0, resu, i)
        return resu

    def gauss_curvature(self, recache=False):
        r"""
        Return the gauss curvature of the submanifold.

        The Gauss curvature is the product or the principal curvatures, or
        equivalently the determinant of the projection operator.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - scalar field on the submanifold equal to the Gauss curvature

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.gauss_curvature().display())
            N --> R
            (th, ph) |--> r^(-2)

        """
        if self._gauss_curvature is not None and not recache:
            return self._gauss_curvature
        a = self.shape_operator(recache)
        e = matrix([[a[i, j].expr() for i in self.irange()] for j in
                    self.irange()]).determinant()
        self._gauss_curvature = self.scalar_field({self.default_chart(): e})
        return self._gauss_curvature

    def principal_directions(self, recache=False):
        r"""
        Return the principal directions of the submanifold.

        The principal directions are the eigenvectors of the projection
        operator. The result is formatted as a list of couples
        (eigenvector, eigenvalue).

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - List of couples  (vectorfields,scalarfield) representing the
          principal directions and the principal curvature associated

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.principal_directions()[0][0].display())
            e_0 = d/dth

        """
        if self._principal_directions is not None and not recache:
            return self._principal_directions
        a = self.shape_operator(recache)
        pr_d = matrix(
            [[a[i, j].expr() for i in self.irange()] for j in
             self.irange()]).eigenvectors_right()
        self._principal_directions = []
        v = self.vector_field()
        counter = self.irange()
        for eigen_space in pr_d:
            for eigen_vector in eigen_space[1]:
                v[self.default_frame(), :] = eigen_vector
                self._principal_directions.append((v.copy(), eigen_space[0]))
                self._principal_directions[-1][0].set_name(
                    "e_%i" % counter.next())
        return self._principal_directions

    def principal_curvatures(self, recache=False):
        r"""
        Return the principal curvatures of the submanifold.

        The principal curvatures are the eigenvalues of the projection operator.
        The resulting scalarfield are named "k_i" with i index of the
        submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - list of scalar field on the submanifold equal to the principal
          curvatures

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.principal_curvatures()[0].display())
            k_0: N --> R
               (th, ph) |--> -1/r
        """
        if self._principal_curvatures is not None and not recache:
            return self._principal_curvatures
        a = self.shape_operator(recache)
        self._principal_curvatures = matrix(
            [[a[i, j].expr() for i in self.irange()] for j in
             self.irange()]).eigenvalues()
        counter = self.irange()
        for i in range(self._dim):
            self._principal_curvatures[i] = self.scalar_field(
                {self.default_chart(): self._principal_curvatures[i]},
                name="k_%i" % counter.next())
        return self._principal_curvatures

    def mean_curvature(self, recache=False):
        r"""
        Return the mean curvature of the submanifold.

        The mean curvature is the arithmetic mean of the principal curvatures,
        or equivalently the trace of the projection operator.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - scalar field on the submanifold equal to the mean curvature

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.mean_curvature().display())
            N --> R
            (th, ph) |--> -1/r


        """
        return self._sgn*sum(self.principal_curvatures(recache)) / self._dim

    def shape_operator(self, recache=False):
        r"""
        Return the shape opeator of the submanifold.

        The shape operator is equal to the second fundamental form with one of
        the indices upped.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``recache`` -- (default: ``False``) if True, the cached value will be
          ignored and all the functions this function depends on will be
          reevaluated (potentially long). Use only after a modification of the
          submanifold

        OUTPUT:

        - (1,1) tensor field on the submanifold equal to the shape operator

        EXAMPLES:

        A sphere embedded in euclidan space::

            sage: M = Manifold(3,'M',structure = "Riemannian")
            sage: N = Manifold(2,'N',ambient = M,structure = "Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r')
            sage: assume(r>0)
            sage: E.<x,y,z> = M.chart()
            sage: phi = N.diff_map(M,{(C,E):[r*sin(th)*cos(ph),
            ....:                            r*sin(th)*sin(ph),r*cos(th)]})
            sage: phi_inv = M.diff_map(N,{(E,C):[arccos(z/r),atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E:sqrt(x**2+y**2+z**2)})
            sage: N.set_immersion(phi,phi_inverse = phi_inv,var = r,
            ....:                 t_inverse = {r:phi_inv_r})
            sage: N.declare_embedding()
            sage: T = N.adapted_chart()
            sage: g = M.metric('g')
            sage: g[0,0],g[1,1],g[2,2]=1,1,1
            sage: print(N.shape_operator()[:])
            [-1/r    0]
            [   0 -1/r]

        """
        if self._shape_operator is not None and not recache:
            return self._shape_operator
        self._shape_operator = self.second_fundamental_form(recache).contract(
            self.induced_metric(recache).inverse())
        return self._shape_operator

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
    M --> R
    (w, x, y, z) |--> sqrt(b^2 + x^2 + y^2 + z^2)/b
    (rh_M, th_M, ph_M, t_M) |--> cosh(rh_M)

Shift vector::

    sage: N.shift().display()
    -(x^2 + y^2 + z^2)/b^2 d/dw - sqrt(b^2 + x^2 + y^2 + z^2)*x/b^2 d/dx -
     sqrt(b^2 + x^2 + y^2 + z^2)*y/b^2 d/dy - sqrt(b^2 + x^2 + y^2 + z^2)*z/b^2 d/dz

The extrinsic curvature (second fundamental form) as a tensor of the ambient
manifold:

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

    def ambient_metric(self):
        r"""
        Return the metric of the ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
        if self._ambient_metric is not None:
            return self._ambient_metric
        self._ambient_metric = self._ambient.metric()
        self._ambient_metric.set_name("g",r"g")
        return self._ambient_metric

    def first_fundamental_form(self):
        r"""
        Return the first fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
        if self._first_fundamental_form is not None:
            return self._first_fundamental_form
        self._first_fundamental_form = self.metric()
        self._first_fundamental_form\
            .set(self._immersion.pullback(self.ambient_metric()))
        self._first_fundamental_form.set_name("gamma",r"\gamma")
        return self._first_fundamental_form

    induced_metric = first_fundamental_form

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
            x/sqrt(x^2 + y^2 + z^2) dx + y/sqrt(x^2 + y^2 + z^2) dy + z/sqrt(x^2 + y^2 + z^2) dz

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to "
                             "perform this calculation")
        if self._difft is not None:
            return self._difft
        self._difft = self._t_inverse[self._var[0]].differential()
        self._difft.set_name("d" + self._var[0]._repr_(),
                             r"\mathrm{d}" + self._var[0]._latex_())
        return self._difft

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
            x/sqrt(x^2 + y^2 + z^2) d/dx + y/sqrt(x^2 + y^2 + z^2) d/dy + z/sqrt(x^2 + y^2 + z^2) d/dz

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform "
                             "this calculation")
        if self._gradt is not None:
            return self._gradt
        self._gradt = self.ambient_metric().inverse().contract(self.difft())
        self._gradt.set_name("grad_" + self._var[0]._repr_(),
                             r"\nabla " + self._var[0]._latex_())
        return self._gradt

    def normal(self):
        r"""
        Return a normal unit vector to the submanifold.

        For now, it can only do by computing gradt() first, i.e. a foliation
        is needed to perform this computation, although it is possible to
        proceed without, using the formula:

        .. MATH::

            n = \overrightarrow{*}(\mathrm{d}x_0\wedge\mathrm{d}x_1\wedge...\wedge\mathrm{d}x_{n-1})

        where the star is the hodge dual operator and de wedge the product on
        the exterior algebra.

        This is not currently possible to implement in sagemanifold because the
        tensors are defined on different domains, despite having the same
        codomain, which should make the contraction possible


        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
            x/sqrt(x^2 + y^2 + z^2) d/dx + y/sqrt(x^2 + y^2 + z^2) d/dy + z/sqrt(x^2 + y^2 + z^2) d/dz

        Or in spherical coordinates::

            sage: print(N.normal().display(T[0].frame(),T[0]))
            d/dr_M

        .. TODO:

            Implement normal() in a way that doesn't need a foliation.

        """
        if self._normal is not None:
            return self._normal
        self._normal = self._sgn*self.lapse()*self.gradt()
        self._gradt.set_name("n",r"n")
        return self._normal

    def ambient_first_fundamental_form(self):
        r"""
        Return the first fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
            r_M^2 dth_M*dth_M + r_M^2*sin(th_M)^2 dph_M*dph_M + 2 dr_M*dr_M

        """
        if self._ambient_first_fundamental_form is not None:
            return self._ambient_first_fundamental_form
        self._ambient_first_fundamental_form = \
            self.ambient_metric() \
            + self.ambient_metric().contract(self.normal())\
            * self.ambient_metric().contract(self.normal())
        self._ambient_first_fundamental_form.set_name("gamma",r"\gamma")
        return self._ambient_first_fundamental_form

    ambient_induced_metric = ambient_first_fundamental_form

    def lapse(self):
        r"""
        Return the lapse function of the foliation

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
            M --> R
            (x, y, z) |--> 1
            (th_M, ph_M, r_M) |--> 1

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        if self._lapse is not None:
            return self._lapse
        self._lapse = 1/(self._sgn * self.ambient_metric()(self.gradt(),
                                                           self.gradt())).sqrt()
        self._lapse.set_name("N", r"N")
        return self._lapse

    def shift(self):
        r"""
        Return the shift function of the foliation

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
            0

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        if self._shift is not None:
            return self._shift
        self._shift = self._adapted_charts[0].frame()[self._dim]\
            - self.lapse() * self.normal()
        self._shift.set_name("beta",r"\beta")
        return self._shift

    def ambient_second_fundamental_form(self):
        r"""
        Return the second fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
            -r_M dth_M*dth_M - r_M*sin(th_M)^2 dph_M*dph_M

        """
        if self._ambient_second_fundamental_form is not None:
            return self._ambient_second_fundamental_form
        nab = self.ambient_metric().connection('nabla', r'\nabla')
        self._ambient_second_fundamental_form = \
            -self.ambient_metric().contract(nab(self.normal())) \
            - nab(self.normal()).contract(self.normal())\
            .contract(self.ambient_metric())\
            * self.normal().contract(self.ambient_metric())
        self._ambient_second_fundamental_form.set_name("K",r"K")
        return self._ambient_second_fundamental_form

    ambient_extrinsic_curvature = ambient_second_fundamental_form

    def second_fundamental_form(self):
        r"""
        Return the second fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

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
            -r dth*dth - r*sin(th)^2 dph*dph

        """
        if self._second_fundamental_form is not None:
            return self._second_fundamental_form
        # self._second_fundamental_form = \
        #    self._immersion.pullback(self.ambient_second_fundamental_form())

        inverse_subs = {v: k for k, v in self._subs[0].items()}
        resu = self._immersion._domain.vector_field_module().tensor((0, 2),
                                                                    name='K',
                                                                    latex_name='K',
                                                                    sym=[
                                                                        (0, 1)],
                                                                    antisym=[])
        for i in self.irange():
            for j in self.irange():
                resu[i, j] = self.ambient_extrinsic_curvature()[
                    self._adapted_charts[0].frame(), [i, j]].expr(
                    self._adapted_charts[0]).subs(inverse_subs)

        self._second_fundamental_form = resu
        return self._second_fundamental_form

    extrinsic_curvature = second_fundamental_form

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
            Tensor field gamma of type (1,1) on the 3-dimensional Riemannian manifold M

        Check that the projector applied to the normal vector is zero::

            sage: N.projector().contract(N.normal()).display()
            0

        """
        if self._projector is not None:
            return self._projector
        self._projector = self._ambient_first_fundamental_form.up(
            self.ambient_metric(), pos=0)
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

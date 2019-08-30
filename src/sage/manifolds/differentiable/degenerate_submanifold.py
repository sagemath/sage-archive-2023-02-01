r"""
Degenerate submanifolds

An *embedded (resp. immersed) degenerate submanifold of a proper 
pseudo-Riemannian manifold* `(M,g)` is an embedded (resp. immersed) 
submanifold `H` of `M` as a differentiable manifold such that pull 
back of the metric tensor `g` via the embedding (resp. immersion) 
endows `H` with the structure of a degenerate manifold.

Degenerate submanifolds are study in many fields of mathematics and physics,
for instance in Differential Geometry (especially in geometry of
lightlike submanifold) and in General Relativity. In geometry of lightlike 
submanifolds, according to the dimension `r` of the radical distribution
(see below for definition of radical distribution), degenerate submanifolds 
have been classify into 4 subgroups: `r-`lightlike submanifolds, Coisotropic 
submanifolds, Isotropic submanifolds and Totally lightlike submanifolds.
(See the book of Krishan L. Duggal and Aurel Bejancu in *REFERENCES*.) 

In the present module, you can definie any of the 4 types but most of the
methods are implemented only for degenerate hypersurfaces who belong to
`r-`lightlike submanifolds. However, their might be generalized to
`1-`lightlike submanifolds. In the litterature there is a new approach
(the rigging technique) for studying `1-`lightlike submanifolds but
here we we the method of Krishan L. Duggal and Aurel Bejancu base on
the screen distribution.

Let `H` be a lightlike hypersurface of a pseudo-Riemannian manifold
`(M,g)`. Then the normal bundle `T^\perp H` intersect the tangent 
bundle `TH`. The radical distribution is defined as
'Rad(TH)=TH\cap T^\perp H'. In case of hypersurfaces, and more
generally `1-`lightlike submanifolds, this is a rank 1 distribution. 
A screen distribution `S(TH)` is a complementary of `Rad(TH)` in `TH`.

Giving a screen distribution `S(TH)` and a null vector field `\xi` 
locally defined and spanning `Rad(TH)`, there exists a unique 
transversal null vector field 'N' locally defined and such that 
`g(N,\xi)=1`. From a transverse vector 'v', the null rigging 'N'
is giving by the formula

.. MATH::

    N = \frac{1}{g(\xi, v)}\left(v-\frac{g(v,v)}{2g(xi, v)}\xi\right)

Tensors on the ambient manifold 'M' are projected on 'H' along 'N'
to obtain induced objects. For instance, induced connection is the
linear connexion defined on H through the Levi-Civitta connection of
'g' along `N`.

To work on a degenerate submanifold, after defining `H` as an instance
of :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`,
with the keyword *structure='degenerate_metric'*, you have to set a transvervector
`v` and screen distribution together with the radical distribution.

An example of degenerate submanifold from General Relativity is the
horizon of the Shawrzschild black hole. Allow us to recall that 
Shawrzschild black hole is the first non-trivial solution of Einstein's
equations. It describes the metric inside a star of radius `R = 2m`,
being `m` the inertial mass of the star. It can be seen as an open
ball in a Lorentzian manifold structure on `\RR^4`::

    sage: M = Manifold(4, 'M', structure="Lorentzian")
    sage: X_M.<t, r, th, ph> = \
    ....: M.chart(r"t r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi") 
    sage: var('m'); assume(m>0)
    m
    sage: g = M.metric()
    sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = \
    ....: -1+2*m/r, 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2

Let us define the horizon as a degenerate hypersurface::

    sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
    sage: H
    degenerate hypersurface H embedded in 4-dimensional differentiable 
    manifold M
    sage: X_H.<ht,hth,hph> = \
    ....: H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta hph:(0,2*pi):\phi")
    sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, \
    ....: name='Phi', latex_name=r'\Phi')
    sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, \
    ....: name='Phi_inv', latex_name=r'\Phi^{-1}')
    sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()

The induced metric is the pullback of `g` on `H`::

    sage: h = H.induced_metric()
    sage: h
    degenerate metric gamma on the degenerate hypersurface H embedded 
    in 4-dimensional differentiable manifold M 
    sage: h.disp()
    gamma = 4*m^2 dhth*dhth + 4*m^2*sin(hth)^2 dhph*dhph

Let us now set a transversal direction and a screen ::

    sage: xi = M.vector_field(); v = M.vector_field(); \
    ....: e1 = M.vector_field(); e2 = M.vector_field(); 
    sage: xi[0] =-1; v[0] = r; v[1] = -r; e1[2] = 1; e2[3] = 1
    sage: H.set_transverse(rigging=v)
    sage: S = H.screen('S', [e1, e2], (xi)); S
    screen distribution S along the degenerate hypersurface H embedded 
    in 4-dimensional differentiable manifold M mapped into the 
    4-dimensional Lorentzian manifold M

In the field of geometry of degenerate submanifolds, it is common to
work in a field `(v_0, v_1, v_2, v_3)` where `span(v_0, v_1)=S(TM)`,
`v_2=\xi` is a null vector field spanning the radical distribution and
`v_3=N` is the corresponding null transverse vector field::

    sage: T = H.adapted_frame(); T
    Vector frame (H, (v_0,v_1,v_2,v_3)) with values on the 4-dimensional 
    Lorentzian manifold M
    sage: T[3].disp()
    N = d/dt - d/dr

The Weingarten map::

    sage: W = H.weingarten_map(); W
    Tensor field nabla_g(xi)|X(H) of type (1,1) along the degenerate 
    hypersurface H embedded in 4-dimensional differentiable manifold M 
    with values on the 4-dimensional Lorentzian manifold M
    sage: W.display()
    nabla_g(xi)|X(H) = -1/4/m d/dt*dt - 1/4/m d/dt*dr
    sage: W.display(T)
    nabla_g(xi)|X(H) = -1/4/m xi*xi^b
    sage: W(xi.along(Phi))
    Vector field along the degenerate hypersurface H embedded in 
    4-dimensional differentiable manifold M with values on the 
    4-dimensional Lorentzian manifold M

The shape operator::

    sage: SO = H.shape_operator(); SO
    Tensor field A^* of type (1,1) along the degenerate hypersurface 
    H embedded in 4-dimensional differentiable manifold M with values 
    on the 4-dimensional Lorentzian manifold M
    sage: SO.disp()
    A^* = 0

The rotation one form::

    sage: tau = H.rotation_one_form(); tau
    Tensor field tau of type (0,1) along the degenerate hypersurface 
    H embedded in 4-dimensional differentiable manifold M with values on 
    the 4-dimensional Lorentzian manifold M
    sage: tau.disp()
    tau = -1/4/m dt - 1/4/m dr
    sage: tau.display(T)
    tau = 1/4/m xi^b

Total umbilicity::

    sage: H.is_maximal()
    True
    sage: H.is_totally_umbilical()
    True
    sage: H.is_totally_geodesic()
    True

REFERENCES:

- [DB1996]_
- [DS2010]_
- [FNO2019]_

"""

from sage.manifolds.differentiable.pseudo_riemannian import \
    PseudoRiemannianManifold
from sage.manifolds.differentiable.degenerate import DegenerateManifold
from sage.manifolds.differentiable.differentiable_submanifold import \
    DifferentiableSubmanifold
from sage.rings.infinity import infinity
from sage.matrix.constructor import matrix
from sage.functions.other import factorial
from sage.symbolic.ring import SR
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from queue import Queue
from sage.symbolic.expression import Expression


class DegenerateSubmanifold(DegenerateManifold, DifferentiableSubmanifold):
    r"""
    Degenerate submanifolds

    An *embedded (resp. immersed) degenerate submanifold of a proper 
    pseudo-Riemannian manifold* `(M,g)` is an embedded (resp. immersed) 
    submanifold `H` of `M` as a differentiable manifold such that pull 
    back of the metric tensor `g` via the embedding (resp. immersion) 
    endows `H` with the structure of a degenerate manifold.
    
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

    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      tuple: `S = (n_+, n_-, n_0)`, where `n_+` (resp. `n_-`, resp. `n_0`) is the
      number of positive terms (resp. negative terms, resp. zero tems) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to `(ndim-1, 0, 1)`, being `ndim` the manifold's dimension
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

    A `2-`dimensional degenerate submanifold of a Lorentzian manifold::

        sage: M = Manifold(4, 'M', structure="Lorentzian")
        sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric') 
        sage: S
        2-dimensional degenerate submanifold S embedded in 4-dimensional 
        differentiable manifold M

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

        A `2-`dimensional degenerate submanifold of a Lorentzian manifold::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric') 
            sage: S
            2-dimensional degenerate submanifold S embedded in 4-dimensional 
            differentiable manifold M

        """

        DegenerateManifold.__init__(self, n, name=name,
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
        self._normal = None
        self._first_fundamental_form = None
        self._induced_metric = None
        self._second_fundamental_form = {}
        self._ambient_metric = None
        self._gauss_curvature = {}
        self._principal_directions = {}
        self._principal_curvatures = {}
        self._mean_curvature = {}
        self._screens = {}
        self._default_screen = None
        self._adapted_frame = {}
        self._shape_operator = {}
        self._rotation_one_form = {}
        signature = self._ambient.metric().signature()
        ndim = self._ambient._dim
        try:
           if signature[0]==ndim or signature[1]==ndim:
            raise ValueError("ambient must be a proper pseudo-Riemannian"+
                              " or a degenerate manifold")
        except TypeError:
          if signature==ndim or signature==-ndim:
            raise ValueError("ambient must be a proper pseudo-Riemannian"+
                              " or a degenerate manifold")
        self._transverse = {}
        
    def _repr_(self):
        r"""
        Return a string representation of the submanifold.

        If no ambient manifold is specified, the submanifold is considered
        as a manifold.

        TESTS::

        A `2-`dimensional degenerate submanifold of a Lorentzian manifold::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric') 
            sage: S.__repr__()
            '2-dimensional degenerate submanifold S embedded in 4-dimensional 
            differentiable manifold M'

        """
        if self._ambient is None:
            return super(DegenerateManifold, self).__repr__()
        if self._ambient._dim-self._dim==1:
          return "degenerate hypersurface {} embedded " \
               "in {}-dimensional differentiable " \
               "manifold {}".format(self._name, self._ambient._dim,
                                    self._ambient._name)  
        return "{}-dimensional degenerate submanifold {} embedded " \
               "in {}-dimensional differentiable " \
               "manifold {}".format(self._dim, self._name, self._ambient._dim,
                                    self._ambient._name)
  
    def ambient_metric(self):
        r"""
        Return the metric of the ambient manifold. The submanifold has to be
        embedded

        OUTPUT:

        - the metric of the ambient manifold

        EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane 
        passing through the origin::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), 0, u, v]}, 
            ....:           name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:               latex_name=r'\Phi^{-1}'); 
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: S.ambient_metric()
            Lorentzian metric g on the 4-dimensional Lorentzian manifold M

        """
        if self._ambient_metric is None:
            if not self._embedded or not isinstance(self._ambient,
                            (PseudoRiemannianManifold, DegenerateManifold)):
                raise ValueError("degenerate submanifold must be embedded in a "
                                 "pseudo-Riemannian or degenerate manifold")
            self._ambient_metric = self._ambient.metric()
        return self._ambient_metric
    
    def default_screen(self):
        r"""
        Return the default screen distribution

        OUTPUT:

        - an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:                   name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                           latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y 
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi)
            sage: S.default_screen()
            screen distribution Sc along the 2-dimensional degenerate submanifold S 
            embedded in 4-dimensional differentiable manifold M mapped into the 
            4-dimensional Lorentzian manifold M
            sage: S.default_screen()[0].disp()
            u d/dt + sqrt(u^2 + v^2) d/dx

        """
        return self._default_screen

    def set_transverse(self, rigging=None, normal=None):
        r"""
        For setting a transversal disttribution of the degenerate submanifold.
        according to the type of the submanifold amoung the 4 possible types,
        one must enter a list of normal transversal vector fields and/or a
        list of transversal and not normal vector fields spanning a transverse
        distribution. 

        INPUT:

        - ``rigging`` -- list or tuple (default: ``None``); list of vector fields
          of the ambient manifold or chart function; of the ambient manifold in 
          the latter case, the corresponding gradient vector field with respect to 
          the ambient metric is calculated; the vectors must be linearly independent, 
          transversal to the submanifold but not normal
        - ``normal`` -- list or tuple (default: ``None``); list of vector fields
          of the ambient manifold or chart function; of the ambient manifold in 
          the latter case, the corresponding gradient vector field with respect to 
          the ambient metric is calculated; the vectors must be linearly independent, 
          transversal and normal to the submanifold

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:                     name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                             latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: S._transverse
            {'normal': (Vector field on the 4-dimensional Lorentzian manifold M,),
            'rigging': (Vector field on the 4-dimensional Lorentzian manifold M,)}

        """
        if isinstance(rigging, (list, tuple)):
            rigging = [elt for elt in rigging]
        else:
            if rigging is not None:
                rigging = [rigging]
        if isinstance(normal, (list, tuple)):
            normal = [elt for elt in normal]
        else:
            if normal is not None:
                normal = [normal]
        frame = self.default_frame()
        im = self.immersion()
        g = self.ambient_metric().along(im)
        nor = []
        rig = []
        l1, l2 = 0, 0
        if normal is not None:
            for u in normal:
                if isinstance(u, Expression):
                    u = self._ambient.scalar_field(u).gradient()
                for v in frame:
                    v = im.pushforward(v)
                    if not g(u.along(im), v).is_zero():
                        raise ValueError("{} is not tangent to {}".format(u.disp(), self._name))
                nor.append(u)
                l1 += 1
        if rigging is not None:
            for u in rigging:
                if isinstance(u, Expression):
                    u = self._ambient.scalar_field(u).gradient()
                rigg = False
                for v in frame:
                    v = im.pushforward(v)
                    if not g(u.along(im), v).is_zero():
                        rigg = True
                if not rigg:
                    raise ValueError("{} is normal to {}".format(u._name, self._name))
                rig.append(u)
                l2 += 1
        if l1+l2!=self._codim:
            raise ValueError("lenght of the transverse must be {}".format(self._codim))
        self._transverse['normal'] = tuple(nor)
        self._transverse['rigging'] = tuple(rig)

    def screen(self, name, screen, rad, latex_name=None):
        r"""
        For setting a screen distribution and vector fields of the radical distribution
        that will be used for computations

        INPUT:

        - ``name`` -- string (default: ``None``); name given to the screen
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to denote
          the screen; if ``None``, the LaTeX symbol is set to ``name``
        - ``screen`` -- list or tuple  of vector fields
          of the ambient manifold or chart function; of the ambient manifold in 
          the latter case, the corresponding gradient vector field with respect to 
          the ambient metric is calculated; the vectors must be linearly independent, 
          tangent to the submanifold but not normal
        - ``rad`` -- -- list or tuple  of vector fields
          of the ambient manifold or chart function; of the ambient manifold in 
          the latter case, the corresponding gradient vector field with respect to 
          the ambient metric is calculated; the vectors must be linearly independent, 
          tangent and normal to the submanifold

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi); Sc
            screen distribution Sc along the 2-dimensional degenerate submanifold S 
            embedded in 4-dimensional differentiable manifold M mapped into the 
            4-dimensional Lorentzian manifold M

        """
        if isinstance(screen, (list, tuple)):
            screen = [elt for elt in screen]
        else:
            screen = [screen]
        if isinstance(rad, (list, tuple)):
            rad = [elt for elt in rad]
        else:
            rad = [rad]
        if name in self._screens:
            if list(screen)==self._screens[name]._screen and list(rad)==self._screens[name]._rad:
                return self._screens[name]
            else:
                raise ValueError("a different screen distribution with the same name \
                  had already been set")
        if len(screen)+len(rad)!=self._dim:
            raise ValueError("total lenght screen+rad must be {}".format(self._dim))
        frame = self.default_frame()
        im = self.immersion()
        g = self.ambient_metric().along(im)
        for (i, u) in enumerate(screen):
            if isinstance(u, Expression):
                u = self._ambient.scalar_field(u).gradient()
            screen[i] = u
            sc = False
            for v in frame:
                v = im.pushforward(v)
                if not g(u.along(im), v).is_zero():
                    sc = True
            if not sc:
                raise ValueError("{} cannot belong to a screen distribution".format(u))
        for (i, u) in enumerate(rad):
            if isinstance(u, Expression):
                u = self._ambient.scalar_field(u).gradient()
            rad[i] = u
            normal = False
            for v in frame:
                v = im.pushforward(v)
                if not g(u.along(im), v).is_zero():
                    raise ValueError("{} is not normal to {}".format(u, self._name))
        self._screens[name] = Screen(self, name, tuple(screen), tuple(rad), latex_name=latex_name)
        self._default_screen = self._screens[name]  
        return self._screens[name]

    def induced_metric(self):
        r"""
        Return the pullback of the ambient metric.

        OUTPUT:

        - induced mettric, as an instance of
          :class:`~sage.manifolds.differentiable.metric.DegenerateMetric`

        EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane 
        passing through the origin::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:             name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:           latex_name=r'\Phi^{-1}'); 
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: h = S.induced_metric(); h
            degenerate metric gamma on the 2-dimensional degenerate 
            submanifold S embedded in 4-dimensional differentiable manifold M
            sage: h[:]
            [ v^2/(u^2 + v^2) -u*v/(u^2 + v^2)]
            [-u*v/(u^2 + v^2)  u^2/(u^2 + v^2)] 

        """
        if self._induced_metric is None or self._induced_metric._components=={}:
            self._induced_metric = self.metric()
            self._induced_metric.set(
                               self.immersion().pullback(self.ambient_metric()))
            self._induced_metric.set_name("gamma", r"\gamma")
        return self._induced_metric

    def first_fundamental_form(self):
        r"""
        Return the restriction of the ambient metric on vector field
        along the submanifold and tangent to it. It is difference from 
        induced metric who gives the pullback of the ambient metric
        on the submanifold. 

        OUTPUT:

        - the first fundamental form, as an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.TangentTensor`

        EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane 
        passing through the origin::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:         name='Phi', latex_name=r'\Phi'); 
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:           latex_name=r'\Phi^{-1}'); 
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: h = S.first_fundamental_form()

        """
        if self._first_fundamental_form is None:
            g = self.ambient_metric()
            h = g.copy()
            h.set_name(g._name+"|"+self._name, g._latex_name+r"|_"+self._latex_name)
            h = TangentTensor(h, self.immersion())
            self._first_fundamental_form = h
        return self._first_fundamental_form
        
    def _ambient_decomposition(self, screen=None):
        r"""
        Return a list ``[screen, rad, normal, rig]`` where `screen`
        is a list a vector fields on the ambient manifoldspanning 
        the giving screen distribution, `rad` a list of vector fields 
        spanning radical distribution, `normal` list of normal transversal 
        vector fields, and `rig` list of rigging vector fields 
        corresponding to the giving screen. 

        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a list

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:             name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                     latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y 
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi)
            sage: S._ambient_decomposition()
            [[Vector field on the 4-dimensional Lorentzian manifold M],
             [Vector field on the 4-dimensional Lorentzian manifold M],
             (Vector field on the 4-dimensional Lorentzian manifold M,),
             (Vector field on the 4-dimensional Lorentzian manifold M,)]

        """
        try:
            normal = self._transverse['normal']
        except KeyError:
            raise ValueError("Set first transverse by using set_transverse")
        if screen is None:
            screen = self.default_screen()
        if isinstance(screen, Screen):
            rad = screen._rad
            screen = screen._screen
            rig = self._transverse['rigging']
        else:
            raise ValueError("set first a screen")
        if self._codim==1:
            xi = rad[0]
            v = rig[0]
            g = self.ambient_metric()
            N = (1/g(xi, v))*(v-(g(v,v)/(2*g(xi, v)))*xi)
            N.set_name(name='N')
            rig = [N]
        return [screen, rad, normal, rig]
    
    def _adapted_frame_(self, screen=None):
        r"""
        Return a frame `(e_1,\ldots,e_p, xi_1,\ldots, xi_r, v_1,\ldots, v_q, N_1, \ldots, N_n)`
        of the ambient manifold, being `e_i` vector fields
        spanning the giving screen distribution, `xi_i` vector fields spanning
        radical distribution, `v_i` normal transversal vector fields, `N_i`
        rigging vector fields corresponding to the giving screen. 

        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a frame on the ambient manifold

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:             name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                     latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y 
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi)
            sage: T = S._adapted_frame_(); T
            Vector frame (M, (v_0,v_1,v_2,v_3))

        """

        if screen is None:
            screen = self.default_screen()
        if screen is None:
            raise ValueError("set first a screen distribution")
        if screen._name in self._adapted_frame:
            return self._adapted_frame[screen._name]
        decomposition = self._ambient_decomposition(screen)
        sc, rad = decomposition[0], decomposition[1]
        normal, rigging = decomposition[2], decomposition[3]
        A = self._ambient.automorphism_field()
        i = self._ambient._sindex
        for u in sc:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        for u in rad:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        for u in normal:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        for u in rigging:
            for j in self._ambient.irange():
                A[j,i] = u[j]
            i += 1
        f = self._ambient.default_frame()
        GLHPhi = f.along(self.immersion())[0].parent().general_linear_group()
        e = f.new_frame(A, 'v')
        self.set_change_of_frame(f.along(self.immersion()), e.along(
                  self.immersion()), GLHPhi(A.along(self.immersion())))
        b = e.dual_basis()
        if self._codim==1:
            e[self._dim-self._sindex].set_name('N')
            e[self._dim-self._sindex-1].set_name('xi', latex_name=r'\xi')
            b[self._dim-self._sindex].set_name('N^b', latex_name=r'N^\flat')
            b[self._dim-self._sindex-1].set_name('xi^b', latex_name=r'\xi^\flat')
        self._adapted_frame[screen._name] = e
        return e

    def adapted_frame(self, screen=None):
        r"""
        Return a frame `(e_1,\ldots,e_p, xi_1,\ldots, xi_r, v_1,\ldots, v_q, N_1, \ldots, N_n)`
        of the ambient manifold along the submanifold, being `e_i` vector fields
        spanning the giving screen distribution, `xi_i` vector fields spanning
        radical distribution, `v_i` normal transversal vector fields, `N_i`
        rigging vector fields corresponding to the giving screen. 

        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a frame on the ambient manifold along the submanifold

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:             name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                     latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y 
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi)
            sage: T = S.adapted_frame(); T
            Vector frame (S, (v_0,v_1,v_2,v_3)) with values on the 4-dimensional 
            Lorentzian manifold M
            sage: T[0].disp()
            v_0 = u d/dt + sqrt(u^2 + v^2) d/dx
            sage: T[0].disp(T)
            v_0 = v_0
            sage: T[3].disp()
            v_3 = -d/dt

        """

        f = self._adapted_frame_(screen).along(self.immersion())
        b = f.dual_basis()
        if self._codim==1:
            f[self._dim-self._sindex].set_name('N')
            f[self._dim-self._sindex-1].set_name('xi', latex_name=r'\xi')
            b[self._dim-self._sindex].set_name('N^b', latex_name=r'N^\flat')
            b[self._dim-self._sindex-1].set_name('xi^b', latex_name=r'\xi^\flat')
        return self._adapted_frame_(screen).along(self.immersion())

    def second_fundamental_form(self, screen=None):
        r"""

        This method is implemented only for null hypersurfaces. The method
        returns a tensor `B` of type `(o,2)` instance of
        :class:`~sage.manifolds.differentiable.degenerate_submanifold.TangentTensor` 
        such that for two vector fields `U, V` on the ambient manifold along
        the null hypersurface, one has:

        .. MATH::

        \nabla_UV=D(U,V)+B(U,V)N

        being `\nabla` the ambient connection, `D` the induced connection
        and `N` the choosen rigging.
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - an instance of
        :class:`TangentTensor`

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....:                     ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \ 
            ....:                     hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:                     name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....:                           1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi = M.vector_field(); v = M.vector_field(); 
            sage: e1 = M.vector_field(); e2 = M.vector_field()
            sage: xi[0] =-1; v[0] = r; v[1] = -r; e1[2] = 1; e2[3] = 1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi))
            sage: B = H.second_fundamental_form(); B
            Field of symmetric bilinear forms B along the degenerate hypersurface 
            H embedded in 4-dimensional differentiable manifold M with values on 
            the 4-dimensional Lorentzian manifold M
            sage: B.disp()
            B = 0

        """
        if self._ambient._dim-self._dim != 1:
            raise ValueError("'second_fundamental_form' is defined"+
                                      " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name not in self._second_fundamental_form:
            resu = self._ambient.vector_field_module() \
                .tensor((0, 2), name='B', latex_name='B', sym=[(0, 1)], antisym=[]) \
                .along(self.immersion())
            f = self.adapted_frame(screen)
            A = self.screen_projection(self._weingarten_map())
            resu[f, :] = A[f, :]
            im = self.immersion()
            self._second_fundamental_form[screen._name] = TangentTensor(resu, im)
        return self._second_fundamental_form[screen._name]

    extrinsic_curvature = second_fundamental_form

    def projection(self, tensor, screen=None):
        r"""

        For a giving tensor `T` of type `(r, 1)` on the ambient manifold, this 
        method returns the tensor `T'` of type `(r,1)` such that for `r`
        vector fields `v_1,\ldots,v_r`, `T'(v_1,\ldots,v_r)` is the projection
        of  `T(v_1,\ldots,v_r)` on ``self`` along the bundle spans by the
        transversal vector fields provided by the 
        :method:`~sage.manifolds.differentiable.degenerate_submanifold.set_transverse`

        INPUT:

        - ``tensor`` -- a tensor of type `(r,1)` on the ambient manifold;

        OUTPUT:

        - a tensor of type `(r,1)` on the ambient manifold along ``self``;

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y 
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi)
            sage: U1 = S.projection(U); U1
            Vector field along the 2-dimensional degenerate submanifold S embedded in 
            4-dimensional differentiable manifold M with values on the 4-dimensional 
            Lorentzian manifold M
            sage: U1.disp()
            u d/dt + sqrt(u^2 + v^2) d/dx

        """
        if tensor.tensor_type()[0]!=1:
            raise NotImplementedError("``projection`` is implemented only for "+
                                      "tensors with 1 as contravariant order")
        #TensorField.display(T, basis=self._ambient.default_frame())
        frame = self.adapted_frame(screen)
        T = tensor.copy()
        try:
            T = T.along(self.immersion())
        except ValueError:
            pass
        T.disp(frame)       
        for i in self._ambient.index_generator(T.tensor_rank()):
            if i[0] in range(self._dim-self._sindex,self._ambient._dim-self._sindex):
                T[frame, i] = 0
        if tensor._latex_name is None:
            T.set_name(tensor._name)
        else:
            T.set_name(tensor._name, latex_name=tensor._latex_name)
        return T

    def screen_projection(self, tensor, screen=None):
        r"""

        For a giving tensor `T` of type `(r, 1)` on the ambient manifold, this 
        method returns the tensor `T'` of type `(r,1)` such that for `r`
        vector fields `v_1,\ldots,v_r`, `T'(v_1,\ldots,v_r)` is the projection
        of  `T(v_1,\ldots,v_r)` on the bundle spans by``screen`` along the 
        bundle spans by the transversal plus the radical vector fields provided.

        INPUT:

        - ``tensor`` -- a tensor of type `(r,1)` on the ambient manifold;

        OUTPUT:

        - a tensor of type `(r,1)` on the ambient manifold;

        EXAMPLES:

        The lightcone of the 4-dimensional Minkowski space `\RR^4_1`::
            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:     name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:     latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y 
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi)
            sage: U1 = S.screen_projection(U); U1
            Vector field along the 2-dimensional degenerate submanifold S embedded 
            in 4-dimensional differentiable manifold M with values on the 4-dimensional 
            Lorentzian manifold M
            sage: U1.disp()
            u d/dt + sqrt(u^2 + v^2) d/dx

        """
        if tensor.tensor_type()[0]!=1:
            raise NotImplementedError("``projection`` is implemented only for "+
                                      "tensors with 1 as contravariant order")
        frame = self.adapted_frame(screen)
        T = tensor.copy()
        try:
            T = T.along(self.immersion())
        except ValueError:
            pass
        T.disp(frame)
        for i in self._ambient.index_generator(T.tensor_rank()):
            if i[0] in range(self._dim-self._sindex-1,self._ambient._dim-self._sindex):
                T[frame, i] = 0
        if tensor._latex_name is None:
            T.set_name(tensor._name)
        else:
            T.set_name("P"+tensor._name, latex_name=r'P'+tensor._latex_name)
        return T
     
    def gauss_curvature(self, screen=None):
        r"""
        Gauss curvature is the product of all  eigenfunctions of the shape operator.
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a scalar function on ``self``

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....:       1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi = M.vector_field(); v = M.vector_field(); e1 = M.vector_field() 
            sage: e2 = M.vector_field(); xi[0]=-1; v[0]=r; v[1]=-r; e1[2]=1; e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi)); T = H.adapted_frame()
            sage: K = H.gauss_curvature(); K
            Scalar field on the degenerate hypersurface H embedded in 4-dimensional
            differentiable manifold M
            sage: K.disp(T)
            H --> R

        """
        if self._ambient._dim-self._dim != 1:
            raise ValueError("'gauss_curvature' is defined"+
                                      " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name not in self._gauss_curvature:
            f = self.adapted_frame()
            A = self.shape_operator()
            self._gauss_curvature[screen._name] = self.scalar_field(
              {chart: A[f,:,chart].determinant()
                for chart in self.top_charts()})
        return self._gauss_curvature[screen._name]

    def principal_directions(self, screen=None):
        r"""

        Principal directions are eigenvectors of the shape operator. This
        method is implemented only for hypersurfaces.
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a dictionnary

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....:       1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi = M.vector_field(); v = M.vector_field(); e1 = M.vector_field() 
            sage: e2 = M.vector_field(); xi[0]=-1; v[0]=r; v[1]=-r; e1[2]=1; e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi)); T = H.adapted_frame()
            sage: V = H.principal_directions()
            sage: V[2][0].disp(T)
            e_2 = xi

        """
        if self._codim != 1:
            raise ValueError("'principal directions' is defined"+
                                      " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._principal_directions:
            return self._principal_directions[screen._name]
        a = self.shape_operator(screen)
        frame = self.adapted_frame(screen)
        pr_d = matrix(
            [[a[frame, :][i, j].expr() for i in self.irange()]
             for j in self.irange()]).eigenvectors_right()
        res = []
        counter = self.irange()
        for eigen_space in pr_d:
            for eigen_vector in eigen_space[1]:
                v = self._ambient.vector_field(name="e_{}".format(next(counter))).along(self.immersion())
                v[frame, :] = [elt for elt in eigen_vector]+ [0]
                res.append((TangentTensor(v, self.immersion()), eigen_space[0]))
                #res[-1][0].set_name("e_{}".format(next(counter)))
        self._principal_directions[screen._name] = res
        return res

    def principal_curvatures(self, screen=None):
        r"""

        Principal curvatures are eigenfunctions of the shape operator. This
        method is implemented only for hypersurfaces.
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a dictionnary

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field(); 
            sage: e2=M.vector_field();xi[0]=-1;v[0] =r; v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi)); 
            sage: rho = H.principal_curvatures()
            sage: rho[2].disp()
            k_2: H --> R
                (ht, hth, hph) |--> 0

        """
        if self._codim != 1:
            raise ValueError("'principal_curvatures' is defined"+
                                      " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._principal_curvatures:
            return self._principal_curvatures[screen._name]
        a = self.shape_operator(screen)
        frame = self.adapted_frame()
        res = matrix([[a[frame, :][i, j].expr() for i in self.irange()] 
                      for j in self.irange()]).eigenvalues()
        counter = self.irange()
        for i in range(self._dim):
            res[i] = self.scalar_field(res[i], name="k_{}".format(next(counter)))
        self._principal_curvatures[screen._name] = res
        return res

    def mean_curvature(self, screen=None):
        r"""

        Mean curvature is the sum of principal curvatures. This
        method is implemented only for hypersurfaces.
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a dictionnary

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1; v[0] =r; v[1]=-r; e1[2]=1; e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi)); 
            sage: m = H.mean_curvature(); m
            Scalar field on the degenerate hypersurface H embedded in 4-dimensional 
            differentiable manifold M
            sage: m.disp()
            H --> R
              (ht, hth, hph) |--> 0

        """
        if self._codim != 1:
            raise ValueError("'mean_curvature' is defined"+
                    " only for hypersurfaces.")
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._mean_curvature:
            return self._mean_curvature[screen._name]
        self._mean_curvature[screen._name] = self.scalar_field({chart: sum(
            self.principal_curvatures(screen)).expr(chart) / self._dim
                                                  for chart in
                                                  self.top_charts()})
        return self._mean_curvature[screen._name]

    def is_tangent(self, v):
        r"""
        This method determine whether a vector field on the ambient
        manifold along ``self`` is tangent to ``self`` or not
        
        INPUT:

        - ``v`` -- field on the ambient manifold along ``self``

        OUTPUT:

        - True if ``v`` is everywhere tangent to ``self`` or False if not

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1; v[0] =r; v[1]=-r; e1[2]=1; e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi)); 
            sage: H.is_tangent(xi.along(Phi))
            True
            sage: H.is_tangent(v.along(Phi))
            False

        """
        g = self.ambient_metric()
        im = self.immersion()
        decomposition = self._ambient_decomposition()
        rad, normal = decomposition[1], decomposition[2]
        for u in rad:
            if not g.along(im)(u.along(im),v).is_zero():
                return False
        for u in normal:
            if not g.along(im)(u.along(im),v).is_zero():
                return False
        return True 

    def _weingarten_map(self, screen=None):
        r"""

        This method is implemented only for hypersurfaces. 
        *Weigarten map* is the `1-`form `W` defined for a vector field
        `U` tangent to ``self`` by 

        .. MATH::

        W(U)=\nabla_U\xi

        being `\nabla` the Levi-Civita connection of the ambient manifold
        and `\xi` the choosen vector field spanning the radical distribution.
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - a dictionnary

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r; v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi)); 
            sage: H._weingarten_map()
            Tensor field of type (1,1) on the 4-dimensional Lorentzian manifold M

        """
        if self._ambient._dim-self._dim != 1:
            raise NotImplementedError("'shape_operator' is implemented"+
                                      " only for hypersurfaces.")
        rad = self._ambient_decomposition(screen)[1]
        nab = self.ambient_metric().connection()
        xi = rad[0]
        return nab(xi)

    def weingarten_map(self, screen=None):
        r"""

        This method is implemented only for hypersurfaces. 
        *Weigarten map* is the `1-`form `W` defined for a vector field
        `U` tangent to ``self`` by 

        .. MATH::

        W(U)=\nabla_U\xi

        being `\nabla` the Levi-Civita connection of the ambient manifold
        and `\xi` the choosen vector field spanning the radical distribution.
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - tensor of type `(1,1)` instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.TangentTensor`

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, \
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1; v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi));  T = H.adapted_frame()
            sage: W = H.weingarten_map(); W
            Tensor field nabla_g(xi)|X(H) of type (1,1) along the degenerate 
            hypersurface H embedded in 4-dimensional differentiable manifold 
            M with values on the 4-dimensional Lorentzian manifold M
            sage: W.display()
            nabla_g(xi)|X(H) = -1/4/m d/dt*dt - 1/4/m d/dt*dr
            sage: W.display(T)
            nabla_g(xi)|X(H) = -1/4/m xi*xi^b
            sage: W(xi.along(Phi)).disp(T)
            -1/4/m xi

        """

        im = self.immersion()
        T = self.projection(self._weingarten_map(screen=screen))
        T.set_name("nabla_g(xi)|X("+self._name+")", 
          latex_name=r'\nabla_g(\xi)|_{\mathfrak{X}('+self._latex_name+r')}')
        return TangentTensor(T, im)

    def shape_operator(self, screen=None):
        r"""

        This method is implemented only for hypersurfaces. 
        *shape operator* is the projection of the Weingarten map
        on the screen distribution along the radical distribution. 
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - tensor of type `(1,1)` instance of
          :class:`~sage.manifolds.differentiable.degenerate_submanifold.TangentTensor`

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, \
            ....: 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi)); 
            sage: SO = H.shape_operator(); SO
            Tensor field A^* of type (1,1) along the degenerate hypersurface 
            H embedded in 4-dimensional differentiable manifold M with values 
            on the 4-dimensional Lorentzian manifold M
            sage: SO.disp()
            A^* = 0

        """
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._shape_operator:
            return self._shape_operator[screen._name]
        im = self.immersion()
        T = -self.screen_projection(self. \
                  _weingarten_map(screen=screen), screen=screen)
        T.set_name("A^*", latex_name =  r'A^\star')
        A = TangentTensor(T, im)
        self._shape_operator[screen._name] = A
        return A

    def rotation_one_form(self, screen=None):
        r"""

        This method is implemented only for hypersurfaces. 
        *rotation one form* is the `1-`form `\tau` such that for any
        vector field `U` tangent to ``self``, the Weingarten map `W`
        and the shape operator `A^*` are related by

        .. MATH::

        W(U) = -A^*(U)-\tau(U)\xi

        being `\xi` the choosen vector field spanning the radical distribution.  
        
        INPUT:

        - ``screen`` -- (default: ``None``); an instance of
          :class:`Screen`
          if ``None`` default screen is used.
        OUTPUT:

        - tensor of type `(0,1)` instance of
          :class:`TangentTensor`

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi));  T = H.adapted_frame()
            sage: tau = H.rotation_one_form(); tau
            Tensor field tau of type (0,1) along the degenerate hypersurface H 
            embedded in 4-dimensional differentiable manifold M with values on 
            the 4-dimensional Lorentzian manifold M
            sage: tau.disp()
            tau = -1/4/m dt - 1/4/m dr
            sage: tau.display(T)
            tau = 1/4/m xi^b

        """
        if screen is None:
            screen = self.default_screen()
        if screen._name in self._rotation_one_form:
            return self._rotation_one_form[screen._name]
        im = self.immersion()
        tau = self._ambient.diff_form(1, name='tau', latex_name=r'\tau').along(self.immersion())
        frame = self.adapted_frame(screen)
        W = -self._weingarten_map(screen).along(self.immersion())
        i = self._dim-self._sindex-1
        tau[frame, i] = W[frame, i, i]
        tau.set_name("tau", latex_name = r'\tau')
        tau = TangentTensor(tau, im)
        self._rotation_one_form[screen._name] = tau
        return tau

    def is_totally_geodesic(self):
        r"""

        This method is implemented only for hypersurfaces. 
        An hypersurface is said to be *totally geodesic* when all
        principal curvatures are `0`.

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi));
            sage: H.is_totally_geodesic()
            True

        """
        cur = self.principal_curvatures()
        for elt in cur:
            if not elt.is_zero():
                return False
        return True

    def is_totally_umbilical(self):
        r"""

        This method is implemented only for hypersurfaces. 
        An hypersurface is said to be *totally umbilical* when apart from
        eigenfunctions associated to vector of the radical distribution,
        all the others are equal.

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi));
            sage: H.is_totally_umbilical()
            True

        """
        screen = self.default_screen()
        B = self.second_fundamental_form(screen)
        g = self.first_fundamental_form()
        #elt = screen[0]
        #rho = B(elt, elt)/g(elt, elt)
        for elt in screen:
            try:
                elt = elt.along(self.immersion())
            except ValueError:
                pass
            if not g(elt, elt)==0:
                rho = B(elt, elt)/g(elt, elt)
                break
        for elt in screen:
            try:
                elt = elt.along(self.immersion())
            except ValueError:
                pass
            if not B(elt, elt)==rho*g(elt, elt):
                return False
        return True

    def is_maximal(self):
        r"""

        This method is implemented only for hypersurfaces. 
        An hypersurface is said to be *totally umbilical* when the
        sum of all eigenfunctions is `0`.

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = M.chart(r"t r:(0,oo) th:(0,pi):\theta \
            ....: ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta \
            ....: hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, 
            ....:       name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, 
            ....:       name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = -1+2*m/r, 2*m/r, \
            ....: 1+2*m/r, r^2, r^2*sin(th)^2
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1
            sage: H.set_transverse(rigging=[v])
            sage: S = H.screen('S', [e1, e2], (xi));
            sage: H.is_maximal()
            True

        """
        m = self.mean_curvature()
        if m.is_zero():
            return True
        return False
#**************************************************************************************

from sage.manifolds.differentiable.vectorfield_module import VectorFieldModule

class Screen(VectorFieldModule):
    r"""
    Let `H` be a lightlike submanifold embedded in a pseudo-Riemannian 
    manifold `(M,g)` with `Phi` the embedding map. A screen distribution
    is a complementary `S` of the radical distribution `Rad(TM)=TH\cap 
    T^\perp H` in `TH`. One then has 

    .. MATH::

        TH=S\oplus_{orth}Rad(TH)

    INPUT:

    - ``submanifold`` -- a lightlike submanifold instance of 
      :class:`DegenerateSubmanifold`
    - ``name`` -- name given to the screen distribution
    - ``screen`` -- vector fields of the ambient manifold which
      span the screen distribution
    - ``rad`` -- vector fields of the ambient manifold which
      span the radical distribution
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the 
      screen distribution; if ``None``, it is formed from ``name``

    OUTPUT:

    - an instance of
      :class:`sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`

    EXAMPLES:

    The horizon of the Shawrzschild black hole::

        sage: M = Manifold(4, 'M', structure="Lorentzian")
        sage: X_M.<t, r, th, ph> = \
        ....: M.chart(r"t r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi") 
        sage: var('m'); assume(m>0)
        m
        sage: g = M.metric()
        sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = \
        ....: -1+2*m/r, 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2
        sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
        sage: X_H.<ht,hth,hph> = \
        ....: H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta hph:(0,2*pi):\phi")
        sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, \
        ....: name='Phi', latex_name=r'\Phi')
        sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, \
        ....: name='Phi_inv', latex_name=r'\Phi^{-1}')
        sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
        sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
        sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1

    A screen distribution for the Shawrzschild black hole horizon::

        sage: H.set_transverse(rigging=v)
        sage: S = H.screen('S', [e1, e2], (xi)); S
        screen distribution S along the degenerate hypersurface H embedded 
        in 4-dimensional differentiable manifold M mapped into the 
        4-dimensional Lorentzian manifold M

    The corresponding normal tangent null vector field and null
    transversal vector field::

        sage: xi = S.normal_tangent_vector(); xi.disp()
        xi = -d/dt
        sage: N = S.rigging(); N.disp()
        N = d/dt - d/dr

    Those vector fields are normalized by `g(xi,N)=1`::

        sage: g.along(Phi)(xi, N).disp()
        g(xi,N): H --> R
        (ht, hth, hph) |--> 1

    Section of the lightcone of the Minkowski space with a hyperplane 
    passing through the origin::

        sage: M = Manifold(4, 'M', structure="Lorentzian")
        sage: X.<t,x,y,z> = M.chart()
        sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
        sage: X_S.<u,v> = S.chart()
        sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
        ....:               name='Phi', latex_name=r'\Phi')
        sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
        ....:                       latex_name=r'\Phi^{-1}')
        sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
        sage: g = M.metric()
        sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            
    A screen distribution for this section of the lightcone::

        sage: V = M.vector_field(); V[3] = 1
        sage: S.set_transverse(rigging=t, normal=V)
        sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
        sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
        sage: Sc = S.screen('Sc', U, xi); Sc
        screen distribution Sc along the 2-dimensional degenerate submanifold S 
        embedded in 4-dimensional differentiable manifold M mapped into the 
        4-dimensional Lorentzian manifold M

    Different lists of vector fields spanning respectively the radical
    distribution and the complementary of the normal distribution in
    the transverse bundle::

        sage:: Rad = Sc.normal_tangent_vector(); Rad[0].disp()
        sqrt(u^2 + v^2) d/dt + u d/dx + v d/dy
        sage: rig = Sc.rigging(); rig[0].disp()
        -d/dt        

        """
    
    def __init__(self, submanifold, name, screen, rad, latex_name=None):
        r"""

        TESTS::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi); Sc
            screen distribution Sc along the 2-dimensional degenerate submanifold S 
            embedded in 4-dimensional differentiable manifold M mapped into the 
            4-dimensional Lorentzian manifold M

        """
        if not isinstance(submanifold, DegenerateSubmanifold):
            raise TypeError("the first argument must be a null submanifold")
        VectorFieldModule.__init__(self, submanifold, dest_map=submanifold.immersion())
        try:
            self._rigging = submanifold._transverse['rigging']
        except KeyError:
            raise ValueError("set first transverse by using `set_transverse` method")
        self._screen = list(screen)
        self._rad = list(rad)
        self._name = name
        if not latex_name:
            latex_name = name
        self._latex_name = latex_name

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi); Sc._repr_()
            'screen distribution Sc along the 2-dimensional degenerate submanifold S 
            embedded in 4-dimensional differentiable manifold M mapped into the 
            4-dimensional Lorentzian manifold M'

        """
        description = "screen distribution "+self._name
        if self._dest_map is self._domain.identity_map():
            description += " on the {}".format(self._domain)
        else:
            description += (" along the {}".format(self._domain)
                            + " mapped into the {}".format(self._ambient_domain))
        return description

    def __getitem__(self, i):
        r"""
        Access vector fields spanning the screen distribution.

        INPUT:

        - ``i`` -- index of the coordinate; if the slice ``[:]``, then all
            the coordinates are returned

        OUTPUT:

        - a vector field on the ambient manifold that belong to 
        the screen distribution

        TESTS::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi)
            sage: Sc.__getitem__(0)
            Vector field along the 2-dimensional degenerate submanifold S 
            embedded in 4-dimensional differentiable manifold M with 
            values on the 4-dimensional Lorentzian manifold M

        """
        sc = [elt.along(self._domain.immersion()) for elt in self._screen]
        return sc[i-self._domain._sindex]
  
    def normal_tangent_vector(self):
        r"""

        return either a list ``Rad`` of vector fields spanning the radical 
        distribution or (in case of hypersurface) a normal tangent null 
        vector field spanning the radical distribution.

        OUTPUT:

        - either a list made by vector fields or a vector field in 
          case of hypersurface

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = \
            ....: M.chart(r"t r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = \
            ....: -1+2*m/r, 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = \
            ....: H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, \
            ....: name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, \
            ....: name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1

        A screen distribution for the Shawrzschild black hole horizon::

            sage: H.set_transverse(rigging=v)
            sage: S = H.screen('S', [e1, e2], (xi));

        The corresponding normal tangent null vector field::

            sage: Rad = S.normal_tangent_vector(); Rad.disp()
            xi = -d/dt

        """
        rad = [elt.along(self._domain.immersion()) for elt in self._rad]
        if self._domain._codim==1:
            xi = rad[0]
            xi.set_name(name='xi', latex_name=r'\xi')
            return xi
        else:
            return rad

    def rigging(self):
        r"""

        Return either a list ``Rad`` of vector fields spanning the 
        complementary of the normal distribution `T^\perp H` in the 
        transverse bundle or (When `H` is a null hypersurface) the
        null transversal vector field defined in `\cite{DB1996]`.

        OUTPUT:

        - either a list made by vector fields or a vector field in 
          case of hypersurface

        EXAMPLES:

        The horizon of the Shawrzschild black hole::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X_M.<t, r, th, ph> = \
            ....: M.chart(r"t r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi") 
            sage: var('m'); assume(m>0)
            m
            sage: g = M.metric()
            sage: g[0,0], g[0,1], g[1,1], g[2,2], g[3,3] = \
            ....: -1+2*m/r, 2*m/r, 1+2*m/r, r^2, r^2*sin(th)^2
            sage: H = Manifold(3, 'H', ambient=M, structure='degenerate_metric')
            sage: X_H.<ht,hth,hph> = \
            ....: H.chart(r"ht:(-oo,oo):t hth:(0,pi):\theta hph:(0,2*pi):\phi")
            sage: Phi = H.diff_map(M, {(X_H, X_M): [ht, 2*m,hth, hph]}, \
            ....: name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(H, {(X_M, X_H): [t,th, ph]}, \
            ....: name='Phi_inv', latex_name=r'\Phi^{-1}')
            sage: H.set_immersion(Phi, inverse=Phi_inv); H.declare_embedding()
            sage: xi=M.vector_field(); v=M.vector_field(); e1=M.vector_field()
            sage: e2=M.vector_field();xi[0]=-1;v[0]=r;v[1]=-r;e1[2]=1;e2[3]=1

        A screen distribution for the Shawrzschild black hole horizon::

            sage: H.set_transverse(rigging=v)
            sage: S = H.screen('S', [e1, e2], (xi));

        The corresponding null transversal vector field::

            sage: N = S.rigging(); N.disp()
            N = d/dt - d/dr

        """
        im = self._domain.immersion()
        rig = [elt.along(im) for elt in self._domain._transverse['rigging']]
        if self._domain._codim!=1:
            return rig
        xi = self.normal_tangent_vector()
        v = rig[0]
        g = self._domain.ambient_metric().along(im)
        N = (1/g(xi, v))*(v-(g(v,v)/(2*g(xi, v)))*xi)
        N.set_name(name='N')
        return N

    
#*******************************************************************************************

from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.tensor.modules.tensor_with_indices import TensorWithIndices

class TangentTensor(TensorFieldParal):
    r"""
    Let `S` be a lightlike submanifold embedded in a pseudo-Riemannian 
    manifold `(M,g)` with `Phi` the embedding map. Let `T1` be a tensor 
    on `M` along `S` or not. `TangentTensor(T1,Phi)` returns restriction
    `T2` of `T1` along `S` that in addition can be applied only on vector
    fields tangent to `S`, when `T1` has a covariant part.

    INPUT:

    - ``tensor`` -- a tensor field on the ambient manifold

    OUTPUT:

    - a tensor field on the ambient manifold along the submanifold
        
    EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane 
        passing through the origin::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi); 
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: V1 = M.vector_field().along(Phi); V1[0] = 1; V1[1]=1
            sage: T1(V1).disp()
            d/dt
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi); V2 = S.projection(V1)
            sage: T2(V2).disp()
            u/sqrt(u^2 + v^2) d/dt

        Of course `T1` and `T2` give the same output on vector fields tangent to S::

            sage: T1(xi.along(Phi)).disp()
            sqrt(u^2 + v^2) d/dt
            sage: T2(xi.along(Phi)).disp()
            sqrt(u^2 + v^2) d/dt



        """
    
    def __init__(self, tensor, embedding):
        r"""

        TESTS::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi); 
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: V1 = M.vector_field().along(Phi); V1[0] = 1; V1[1]=1
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi); V2 = S.projection(V1)
            sage: T1(xi.along(Phi)).disp()
            sqrt(u^2 + v^2) d/dt
            sage: T2(xi.along(Phi)).disp()
            sqrt(u^2 + v^2) d/dt
        
        """
        if not isinstance(tensor, TensorField):
            raise TypeError("the second argument must be a tensor field")
        self._tensor = tensor
        self._tensor.set_name(tensor._name, latex_name=tensor._latex_name)
        self._embedding = embedding
        try:
            tensor = tensor.along(embedding)
        except ValueError:
            pass
        if isinstance(tensor, TensorFieldParal):
            TensorFieldParal.__init__(self, tensor._vmodule, tensor._tensor_type, name=tensor._name,
                   latex_name=tensor._latex_name, sym=tensor._sym, antisym=tensor._antisym)
        else:
            TensorField.__init__(self, tensor._vmodule, tensor._tensor_type, name=tensor._name,
                   latex_name=tensor._latex_name, sym=tensor._sym, antisym=tensor._antisym)
        f = tensor._domain._ambient.default_frame().along(embedding)
        self[f, :] = tensor[f, :]
        
    def __call__(self, *args):
        r"""

        TESTS::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi); 
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: V1 = M.vector_field().along(Phi); V1[0] = 1; V1[1]=1
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi); V2 = S.projection(V1)
            sage: T1(xi.along(Phi)).disp()
            sqrt(u^2 + v^2) d/dt
            sage: T2(xi.along(Phi)).disp()
            sqrt(u^2 + v^2) d/dt

        """
        for vector in args:
            try:
                vector = vector.along(self._embedding)
            except ValueError:
                pass
            if not self._domain.is_tangent(vector):
                raise ValueError("The provided vector field is not "+
                        "tangent to {}".format(self._domain._name))
        try:
            return TensorField.__call__(self._tensor.along(self._embedding), *args)
        except ValueError:    
            return TensorField.__call__(self._tensor, *args)

    def extension(self):
        r"""

        Return initial tensor

        EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane 
        passing through the origin::

            sage: M = Manifold(4, 'M', structure="Lorentzian")
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]}, 
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv', 
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = sqrt(x^2+y^2+z^2); U[0] = x
            sage: Sc = S.screen('Sc', U, xi); 
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi); T3 = T2.extension() 
            sage: T3 is T2
            False
            sage: T3 is T1
            True
        
        """
        return self._tensor
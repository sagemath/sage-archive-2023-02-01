r"""
Local Sections

The class :class:`~sage.manifolds.section.Section` implements local sections on
vector bundles. The derived class :class:`~sage.manifolds.section.TrivialSection`
is devoted to sections on trivial parts of a vector bundle.

AUTHORS:

- Michael Jung (2019): initial version

EXAMPLES:

TODO

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.structure.element import ModuleElement
from sage.tensor.modules.free_module_element import FiniteRankFreeModuleElement
from sage.tensor.modules.tensor_with_indices import TensorWithIndices
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

class Section(ModuleElement):
    r"""
    Section in a vector bundle.

    An instance of this class is a local section in a vector bundle `E \to M` of
    class `C^k`. More precisely, a *section* on a subset `U \in M` is a map of
    class `C^k`

    .. MATH::

        s: U \longrightarrow E

    such that

    .. MATH::

        \forall p \in U,\ s(p) \in E_p

    where `E_p` denotes the vector bundle fiber of `E` over the point `p`.

    If `E|_U` is trivial, the class
    :class:`~sage.manifolds.section.TrivialSection` should be used instead.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.section_module.SectionModule`.

    EXAMPLES:

        TODO

    """
    def __init__(self, section_module, name=None, latex_name=None):
        r"""
        Construct a local section.

        TESTS::

            sage: M = Manifold(1, 'S^1', structure='top', start_index=1)
            sage: U = M.open_subset('U')  # the complement of one point
            sage: c_t.<t> =  U.chart('t:(0,2*pi)') # the standard angle coordinate
            sage: V = M.open_subset('V') # the complement of the point t=pi
            sage: M.declare_union(U,V)   # S^1 is the union of U and V
            sage: c_u.<u> = V.chart('u:(0,2*pi)') # the angle t-pi
            sage: t_to_u = c_t.transition_map(c_u, (t-pi,), intersection_name='W',
            ....:                     restrictions1 = t!=pi, restrictions2 = u!=pi)
            sage: u_to_t = t_to_u.inverse()
            sage: W = U.intersection(V)
            sage: E = M.vector_bundle(1, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: transf = phi_U.transition_map(phi_V, [[-1]])
            sage: C0 = E.section_module()
            sage: s = E.section([t], frame=phi_U.frame(), name='s')
            sage: s.add_comp_by_continuation(phi_V.frame(), W)
            sage: s in C0
            True
            sage: TestSuite(s).run()

        """
        ModuleElement.__init__(self, section_module)
        self._smodule = section_module
        self._domain = section_module.domain()
        self._base_space = section_module.base_space()
        self._vbundle = section_module.vector_bundle()
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._init_derived()

    ####### Required methods for ModuleElement (beside arithmetic) #######

    def __bool__(self):
        r"""
        Return ``True`` if ``self`` is nonzero and ``False`` otherwise.

        This method is called by :meth:`is_zero`.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U,V)
            sage: XU = X.restrict(U); XV = X.restrict(V)
            sage: XUV = X.restrict(U.intersection(V))
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U)
            sage: f = E.local_frame('f', domain=V)
            sage: s = E.section(name='s')
            sage: su = E.section(name='s', domain=U)
            sage: sv = E.section(name='s', domain=V)
            sage: su[0] = 0
            sage: sv[0] = 0
            sage: s.set_restriction(su); s.set_restriction(sv)
            sage: bool(s)
            False
            sage: s.is_zero()  # indirect doctest
            True
            sage: sv[0] = 1
            sage: s.set_restriction(sv)
            sage: bool(s)
            True
            sage: s.is_zero()  # indirect doctest
            False

        """
        return any(bool(rst) for rst in self._restrictions.values())

    __nonzero__ = __bool__  # For Python2 compatibility

    ##### End of required methods for ModuleElement (beside arithmetic) #####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: s = E.section(name='s')
            sage: s._repr_()
            'Section s on the 3-dimensional differentiable manifold M with
             values in the real vector bundle E of rank 2'
            sage: repr(s) # indirect doctest
            'Section s on the 3-dimensional differentiable manifold M with
             values in the real vector bundle E of rank 2'
            sage: s # indirect doctest
            Section s on the 3-dimensional differentiable manifold M with values
             in the real vector bundle E of rank 2

        """
        desc = "Section "
        if self._name is not None:
            desc += self._name + " "
        desc += "on the {} ".format(self._domain)
        desc += "with values in the " + self._vbundle.base_field_type() + " "
        desc += "vector bundle " + self._vbundle._name + " "
        desc += "of rank {}".format(self._vbundle.rank())
        return desc

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: omega = E.section(name='\omega')
            sage: omega._latex_()
            '\\omega'
            sage: latex(omega) # indirect doctest
            \omega

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: t = M.tensor_field(1, 3, name='t')
            sage: t._init_derived()

        """
        self._restrictions = {} # dict. of restrictions of self on subdomains
                                # of self._domain, with the subdomains as keys
        self._extensions_graph = {self._domain: self}
        self._restrictions_graph = {self._domain: self}

    def _del_derived(self, del_restrictions=False):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: s = E.section(name='s')
            sage: s._del_derived()

        """
        if del_restrictions:
            self._restrictions.clear()
            self._extensions_graph = {self._domain: self}
            self._restrictions_graph = {self._domain: self}

    def set_name(self, name=None, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``name`` -- string (default: ``None``); name given to the section
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to denote
          the section; if ``None`` while ``name`` is provided, the LaTeX
          symbol is set to ``name``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: s = E.section(); s
            Section on the 3-dimensional topological manifold M with values in
             the real vector bundle E of rank 2
            sage: s.set_name(name='s')
            sage: s
            Section s on the 3-dimensional topological manifold M with values in
             the real vector bundle E of rank 2
            sage: latex(s)
            s
            sage: s.set_name(latex_name=r'\sigma')
            sage: latex(s)
            \sigma
            sage: s.set_name(name='a')
            sage: s
            Section a on the 3-dimensional topological manifold M with values in
             the real vector bundle E of rank 2
            sage: latex(s)
            a

        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name
        for rst in self._restrictions.values():
            rst.set_name(name=name, latex_name=latex_name)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same
        section module.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: s = E.section(name='s', domain=U); s
            Section s on the Open subset U of the 3-dimensional differentiable
             manifold M with values in the real vector bundle E of rank 2
            sage: s1 = s._new_instance(); s1
            Section on the Open subset U of the 3-dimensional differentiable
             manifold M with values in the real vector bundle E of rank 2
            sage: type(s1) == type(s)
            True
            sage: s1.parent() is s.parent()
            True

        """
        return type(self)(self._smodule)

    def _init_components(self, *comp, **kwargs):
        r"""
        Initialize the section components in some given local frames.

        INPUT:

        - ``comp`` -- either the components of the section with respect
          to the local frame specified by the argument ``frame`` or a
          dictionary of components, the keys of which are local frames or
          pairs ``(f,c)`` where ``f`` is a local frame and ``c`` a chart
        - ``frame`` -- (default: ``None``; unused if ``comp`` is a dictionary)
          local frame in which the components are given; if ``None``, the
          default local frame on the domain of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is a dictionary)
          coordinate chart in which the components are expressed; if ``None``,
          the default chart on some subdomain of ``frame`` is assumed

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: s = E.section(name='s')
            sage: s._init_components([1+x, y^2])
            sage: s.display()
            s = (x + 1) e_0 + y^2 e_1
            sage: f = E.local_frame('f')
            sage: Y.<u,v,t> = M.chart()
            sage: s._init_components([2*v, -u], frame=f)
            sage: s.display(f, Y)
            s = 2*v f_0 - u f_1
            sage: s._init_components({e: [1-y, x]})
            sage: s.display()
            s = (-y + 1) e_0 + x e_1
            sage: s._init_components({(f, Y): [t, v^3]})
            sage: s.display(f, Y)
            s = t f_0 + v^3 f_1

        """
        comp0 = comp[0]
        if isinstance(comp0, dict):
            for frame, components in comp0.items():
                chart = None
                if isinstance(frame, tuple):
                    # frame is actually a pair (frame, chart):
                    frame, chart = frame
                self.add_comp(frame)[:, chart] = components
        else:
            if hasattr(comp0, '__getitem__'):
                # comp0 is a list/vector of components
                # otherwise comp is the tuple of components in a specific frame
                comp = comp0
            frame = kwargs.get('frame')
            chart = kwargs.get('chart')
            self.add_comp(frame)[:, chart] = comp

    def domain(self):
        r"""
        Return the manifold on which ``self`` is defined.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.manifold.TopologicalManifold`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0_U = E.section_module(domain=U, force_free=True)
            sage: z = C0_U.zero()
            sage: z.domain()
            Open subset U of the 3-dimensional topological manifold M

        """
        return self._domain

    def base_module(self):
        r"""
        Return the section module on which ``self`` acts as a section.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.section_module.SectionModule`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: s = E.section(domain=U)
            sage: s.base_module()
            Module C^0(U;E) of sections on the Open subset U of the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2

        """
        return self._smodule

    def set_restriction(self, rst):
        r"""
        Define a restriction of ``self`` to some subdomain.

        INPUT:

        - ``rst`` -- :class:`Section` defined on a subdomain of the domain of
          ``self``

        EXAMPLES::

            sage: S2 = Manifold(2, 'S^2', structure='top')
            sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
            sage: S2.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                 restrictions2= u^2+v^2!=0)
            sage: W = U.intersection(V)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: E = S2.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: s = E.section(name='s')
            sage: su = E.section(domain=U, name='s')
            sage: su[:] = x+y, x
            sage: s.set_restriction(su)
            sage: s.display(phi_U.frame())
            s = (x + y) (phi_U^*e_1) + x (phi_U^*e_2)
            sage: s.restrict(U) == su
            True

        """
        self._restrictions[rst._domain] = rst.copy()
        self._restrictions[rst._domain].set_name(name=self._name,
                                                 latex_name=self._latex_name)

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` --
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
          open subset `U` of the section domain `S`

        OUTPUT:

        - :class:`Section` representing the restriction

        EXAMPLES:

        Restrictions of a vector field on the 2-sphere::

            sage: M = Manifold(2, 'S^2', start_index=1)
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
            sage: v = M.vector_field({eN: [1, 0]}, name='v')
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

        At this stage, defining the restriction of ``v`` to the open
        subset ``V`` fully specifies ``v``::

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
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subset of " +
                                 "the field's domain")

            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and subdomain in rst._restrictions:
                    res = rst._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and dom is not self._domain:
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    self._restrictions_graph[subdomain] = rst.restrict(subdomain)
                    return self._restrictions[subdomain]

            # Secondly one tries to get the restriction from one previously
            # defined on a larger domain:
            for dom, ext in self._extensions_graph.items():
                if subdomain in ext._restrictions:
                    res = ext._restrictions_graph[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

            # If this fails, the restriction is created from scratch:
            smodule = self._vbundle.section_module(domain=subdomain)
            res = smodule.element_class(smodule, name=self._name,
                                        latex_name=self._latex_name)
            res._extensions_graph.update(self._extensions_graph)

            for dom, ext in self._extensions_graph.items():
                ext._restrictions[subdomain] = res
                ext._restrictions_graph[subdomain] = res
            for dom, rst in self._restrictions.items():
                if dom.is_subset(subdomain):
                    if rst is not res:
                        res._restrictions.update(rst._restrictions)
                    res._restrictions_graph.update(rst._restrictions_graph)
                    rst._extensions_graph.update(res._extensions_graph)
            self._restrictions[subdomain] = res
            self._restrictions_graph[subdomain] = res
            res._extensions_graph.update(self._extensions_graph)

        return self._restrictions[subdomain]

    def set_comp(self, basis=None):
        r"""
        Return the components of ``self`` in a given vector frame
        for assignment.

        The components with respect to other frames having the same domain
        as the provided vector frame are deleted, in order to avoid any
        inconsistency. To keep them, use the method :meth:`add_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if none is provided, the components are
          assumed to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as a
          :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        EXAMPLES::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
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
            basis = self._smodule._def_frame
            if basis is None: # should be "is still None" ;-)
                raise ValueError("a frame must be provided for the display")
        rst = self.restrict(basis._domain)
        return rst.set_comp(basis)

    def add_comp(self, basis=None):
        r"""
        Return the components of ``self`` in a given vector frame
        for assignment.

        The components with respect to other frames having the same domain
        as the provided vector frame are kept. To delete them, use the
        method :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if ``None``, the components are assumed
          to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as a
          :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        EXAMPLES::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
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

        The components with respect to ``e_uv`` are kept::

            sage: t.display(e_uv)
            t = (u + v) d/dv*du*dv

        """
        if basis is None:
            basis = self._smodule._def_frame
            if basis is None: # should be "is still None" ;-)
                raise ValueError("a frame must be provided for the display")
        rst = self.restrict(basis._domain)
        return rst.add_comp(basis)

    def add_comp_by_continuation(self, frame, subdomain, chart=None):
        r"""
        Set components with respect to a vector frame by continuation of the
        coordinate expression of the components in a subframe.

        The continuation is performed by demanding that the components have
        the same coordinate expression as those on the restriction of the
        frame to a given subdomain.

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

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: # The two open subsets covered by stereographic coordinates (North and South):
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coordinates
            sage: transf = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:             intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:             restrictions2= u^2+v^2!=0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V) # The complement of the two poles
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.vector_field({eU: [x, 2+y]}, name='a')

        At this stage, the vector field has been defined only on the open
        subset ``U`` (through its components in the frame ``eU``)::

            sage: a.display(eU)
            a = x d/dx + (y + 2) d/dy

        The components with respect to the restriction of ``eV`` to the common
        subdomain ``W``, in terms of the ``(u,v)`` coordinates, are obtained
        by a change-of-frame formula on ``W``::

            sage: a.display(eV.restrict(W), c_uv.restrict(W))
            a = (-4*u*v - u) d/du + (2*u^2 - 2*v^2 - v) d/dv

        The continuation consists in extending the definition of the vector
        field to the whole open subset ``V`` by demanding that the components
        in the frame eV have the same coordinate expression as the above one::

            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)

        We have then::

            sage: a.display(eV)
            a = (-4*u*v - u) d/du + (2*u^2 - 2*v^2 - v) d/dv

        and `a` is defined on the entire manifold `S^2`.

        """
        dom = frame._domain
        if not dom.is_subset(self._domain):
            raise ValueError("the vector frame is not defined on a subset " +
                             "of the tensor field domain")
        if chart is None:
            chart = dom._def_chart
        sframe = frame.restrict(subdomain)
        schart = chart.restrict(subdomain)
        scomp = self.comp(sframe)
        resu = self.add_comp(frame) # _del_derived is performed here
        for ind in resu.non_redundant_index_generator():
            resu[[ind]] = dom.scalar_field({chart: scomp[[ind]].expr(schart)})

    def add_expr_from_subdomain(self, frame, subdomain):
        r"""
        Add an expression to an existing component from a subdomain.

        INPUT:

        - ``frame`` -- vector frame `e` in which the components are to be set
        - ``subdomain`` -- open subset of `e`'s domain in which the
          components have additional expressions.

        EXAMPLES:

        We are going to consider a vector field in `\RR^3` along the 2-sphere::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: S = Manifold(2, 'S', structure="Riemannian")
            sage: E.<X,Y,Z> = M.chart()

        Let us define ``S`` in terms of stereographic charts::

            sage: U = S.open_subset('U')
            sage: V = S.open_subset('V')
            sage: S.declare_union(U,V)
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

        The embedding of `S^2` in `\RR^3`::

            sage: phi = S.diff_map(M, {(stereoN, E): [2*x/(1+x^2+y^2),
            ....:                                     2*y/(1+x^2+y^2),
            ....:                                     (x^2+y^2-1)/(1+x^2+y^2)],
            ....:                        (stereoS, E): [2*xp/(1+xp^2+yp^2),
            ....:                                       2*yp/(1+xp^2+yp^2),
            ....:                               (1-xp^2-yp^2)/(1+xp^2+yp^2)]},
            ....:                   name='Phi', latex_name=r'\Phi')

        To define a vector field ``v`` along ``S`` taking its values in ``M``,
        we first set the components on ``U``::

            sage: v = M.vector_field(name='v').along(phi)
            sage: vU = v.restrict(U)
            sage: vU[:] = [x,y,x**2+y**2]

        But because ``M`` is parallelizable, these components can be extended
        to ``S`` itself::

            sage: v.add_comp_by_continuation(E.frame().along(phi), U)

        One can see that ``v`` is not yet fully defined: the components
        (scalar fields) do not have values on the whole manifold::

            sage: sorted(v._components.values())[0]._comp[(0,)].display()
            S --> R
            on U: (x, y) |--> x

        To fix that, we first extend the components from ``W`` to ``V`` using
        :meth:`add_comp_by_continuation`::

            sage: v.add_comp_by_continuation(E.frame().along(phi).restrict(V),
            ....:                            W, stereoS)

        Then, the expression on the subdomain ``V`` is added to the
        already known components on ``S`` by::

            sage: v.add_expr_from_subdomain(E.frame().along(phi), V)

        The definition of ``v`` is now complete::

            sage: sorted(v._components.values())[0]._comp[(2,)].display()
            S --> R
            on U: (x, y) |--> x^2 + y^2
            on V: (xp, yp) |--> 1/(xp^2 + yp^2)

        """
        dom = frame._domain
        if not dom.is_subset(self._domain):
            raise ValueError("the vector frame is not defined on a subset " +
                             "of the tensor field domain")
        if frame not in self.restrict(frame.domain())._components:
            raise ValueError("the tensor doesn't have an expression in "
                             "the frame"+frame._repr_())
        comp = self.comp(frame)
        scomp = self.restrict(subdomain).comp(frame.restrict(subdomain))
        for ind in comp.non_redundant_index_generator():
            comp[[ind]]._express.update(scomp[[ind]]._express)

        rst = self._restrictions.copy()
        self._restrictions = rst

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

        - components in the vector frame ``basis``, as a
          :class:`~sage.tensor.modules.comp.Components`

        EXAMPLES:

        Components of a type-`(1,1)` tensor field defined on two
        open subsets::

            sage: M = Manifold(2, 'M')
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

        Since ``e`` is ``M``'s default frame, the argument ``e`` can
        be omitted::

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
            basis = self._smodule._def_frame
            if basis is None: # should be "is still None" ;-)
                raise ValueError("a frame must be provided for the display")

        rst = self.restrict(basis._domain)
        return rst.comp(basis=basis, from_basis=from_basis)

    def display(self, frame=None, chart=None):
        r"""
        Display the tensor field in terms of its expansion with respect
        to a given vector frame.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to
          which the tensor is expanded; if ``frame`` is ``None`` and ``chart``
          is not ``None``, the coordinate frame associated with ``chart`` is
          assumed; if both ``frame`` and ``chart`` are ``None``, the default
          frame of the domain of definition of the tensor field is assumed
        - ``chart`` -- (default: ``None``) chart with respect to which the
          components of the tensor field in the selected frame are expressed;
          if ``None``, the default chart of the vector frame domain is assumed

        EXAMPLES:

        Display of a type-`(1,1)` tensor field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t[e_xy,:] = [[x, 1], [y, 0]]
            sage: t.add_comp_by_continuation(e_uv, W, c_uv)
            sage: t.display(e_xy)
            t = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: t.display(e_uv)
            t = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv
              + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv

        Since ``e_xy`` is ``M``'s default frame, the argument ``e_xy`` can
        be omitted::

            sage: e_xy is M.default_frame()
            True
            sage: t.display()
            t = x d/dx*dx + d/dx*dy + y d/dy*dx

        Similarly, since ``e_uv`` is ``V``'s default frame, the argument ``e_uv``
        can be omitted when considering the restriction of ``t`` to ``V``::

            sage: t.restrict(V).display()
            t = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv
              + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv

        If the coordinate expression of the components are to be displayed in
        a chart distinct from the default one on the considered domain, then
        the chart has to be passed as the second argument of ``display``.
        For instance, on `W = U \cap V`, two charts are available:
        ``c_xy.restrict(W)`` (the default one) and ``c_uv.restrict(W)``.
        Accordingly, one can have two views of the expansion of ``t`` in the
        *same* vector frame ``e_uv.restrict(W)``::

            sage: t.display(e_uv.restrict(W))  # W's default chart assumed
            t = (1/2*x + 1/2*y + 1/2) d/du*du + (1/2*x + 1/2*y - 1/2) d/du*dv
              + (1/2*x - 1/2*y + 1/2) d/dv*du + (1/2*x - 1/2*y - 1/2) d/dv*dv
            sage: t.display(e_uv.restrict(W), c_uv.restrict(W))
            t = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv
              + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv

        As a shortcut, one can pass just a chart to ``display``. It is then
        understood that the expansion is to be performed with respect to the
        coordinate frame associated with this chart. Therefore the above
        command can be abridged to::

            sage: t.display(c_uv.restrict(W))
            t = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv
              + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv

        and one has::

            sage: t.display(c_xy)
            t = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: t.display(c_uv)
            t = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv
              + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv
            sage: t.display(c_xy.restrict(W))
            t = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: t.restrict(W).display(c_uv.restrict(W))
            t = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv
              + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv

        One can ask for the display with respect to a frame in which ``t`` has
        not been initialized yet (this will automatically trigger the use of
        the change-of-frame formula for tensors)::

            sage: a = V.automorphism_field()
            sage: a[:] = [[1+v, -u^2], [0, 1-u]]
            sage: f = e_uv.new_frame(a, 'f')
            sage: [f[i].display() for i in M.irange()]
            [f_0 = (v + 1) d/du, f_1 = -u^2 d/du + (-u + 1) d/dv]
            sage: t.display(f)
            t = -1/2*(u^2*v + 1)/(u - 1) f_0*f^0
              - 1/2*(2*u^3 - 5*u^2 - (u^4 + u^3 - u^2)*v + 3*u - 1)/((u - 1)*v + u - 1) f_0*f^1
              - 1/2*(v^2 + 2*v + 1)/(u - 1) f_1*f^0
              + 1/2*(u^2 + (u^2 + u - 1)*v - u + 1)/(u - 1) f_1*f^1

        A shortcut of ``display()`` is ``disp()``::

            sage: t.disp(e_uv)
            t = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv
              + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv

        """
        if frame is None:
            frame = self._smodule._def_frame
            if frame is None:  # should be "is still None" ;-)
                raise ValueError("a frame must be provided for the display")
        else:
            try:
                frame0 = frame.frame()
                # if this succeeds, frame is actually not a local frame, but
                # a trivialization
                if chart is None:
                    chart = frame._domain._def_chart
                frame = frame0
            except AttributeError:
                # case of a genuine local frame
                pass
        rst = self.restrict(frame._domain)
        return rst.display(frame, chart)

    disp = display

    def display_comp(self, frame=None, chart=None, only_nonzero=True):
        r"""
        Display the tensor components with respect to a given frame,
        one per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to which
          the tensor field components are defined; if ``None``, then

          * if ``chart`` is not ``None``, the coordinate frame associated to
            ``chart`` is used
          * otherwise, the default basis of the vector field module on which
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

        Display of the components of a type-`(1,1)` tensor field defined
        on two open subsets::

            sage: M = Manifold(2, 'M')
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

        See documentation of
        :meth:`sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.display_comp`
        for more options.

        """
        if frame is None:
            frame = self._smodule._def_frame
            if frame is None:  # should be "is still None" ;-)
                raise ValueError("a frame must be provided for the display")
        rst = self.restrict(frame.domain())
        return rst.display_comp(frame=frame, chart=chart,
                                only_nonzero=only_nonzero)

    def at(self, point):
        r"""
        Value of ``self`` at a point of its domain.

        If the current tensor field is

        .. MATH::

            t:\ U  \longrightarrow T^{(k,l)} M

        associated with the differentiable map

        .. MATH::

            \Phi:\ U \longrightarrow M,

        where `U` and `M` are two manifolds (possibly `U = M` and
        `\Phi = \mathrm{Id}_M`), then for any point `p \in U`, `t(p)`
        is a tensor on the tangent space to `M` at the point `\Phi(p)`.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` in the domain of the tensor field `U`

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          representing the tensor `t(p)` on the tangent vector space
          `T_{\Phi(p)} M`

        EXAMPLES:

        Tensor on a tangent space of a non-parallelizable 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: W = U.intersection(V)
            sage: eU = c_xy.frame() ; eV = c_uv.frame()
            sage: a = M.tensor_field(1, 1, {eU: [[1+y,x], [0,x+y]]}, name='a')
            sage: a.add_comp_by_continuation(eV, W, chart=c_uv)
            sage: a.display(eU)
            a = (y + 1) d/dx*dx + x d/dx*dy + (x + y) d/dy*dy
            sage: a.display(eV)
            a = (u + 1/2) d/du*du + (-1/2*u - 1/2*v + 1/2) d/du*dv
             + 1/2 d/dv*du + (1/2*u - 1/2*v + 1/2) d/dv*dv
            sage: p = M.point((2,3), chart=c_xy, name='p')
            sage: ap = a.at(p) ; ap
            Type-(1,1) tensor a on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: ap.parent()
            Free module of type-(1,1) tensors on the Tangent space at Point p
             on the 2-dimensional differentiable manifold M
            sage: ap.display(eU.at(p))
            a = 4 d/dx*dx + 2 d/dx*dy + 5 d/dy*dy
            sage: ap.display(eV.at(p))
            a = 11/2 d/du*du - 3/2 d/du*dv + 1/2 d/dv*du + 7/2 d/dv*dv
            sage: p.coord(c_uv) # to check the above expression
            (5, -1)

        """
        if point not in self._domain:
            raise ValueError("the {} is not a point in the ".format(point) +
                             "domain of {}".format(self))
        for dom, rst in self._restrictions.items():
            if point in dom:
                return rst.at(point)

    def __getitem__(self, args):
        r"""
        Return a component with respect to some frame.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are returned

        The frame can be passed as the first item of ``args``. If not, the
        default frame of the tensor field's domain is assumed. If ``args``
        is a string, this method acts as a shortcut for tensor contractions
        and symmetrizations, the string containing abstract indices.

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
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
                frame = self._smodule._def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._smodule._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._smodule._def_frame
        return self.comp(frame)[args]

    def __setitem__(self, args, value):
        r"""
        Sets a component with respect to some vector frame.

        INPUT:

        - ``args`` -- list of indices; if ``[:]`` is provided, all the
          components are set; the frame can be passed as the first item
          of ``args``; if not, the default frame of the tensor field's
          domain is assumed
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
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
                frame = self._smodule._def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self._smodule._def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self._smodule._def_frame
        self.set_comp(frame)[args] = value

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        .. NOTE::

            The name and the derived quantities are not copied.

        EXAMPLES:

        Copy of a type-`(1,1)` tensor field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = t.copy(); s
            Tensor field of type (1,1) on
             the 2-dimensional differentiable manifold M
            sage: s.display(e_xy)
            (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s == t
            True

        If the original tensor field is modified, the copy is not::

            sage: t[e_xy,0,0] = -1
            sage: t.display(e_xy)
            t = -d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s.display(e_xy)
            (x + y) d/dx*dx + 2 d/dy*dx + (-y + 1) d/dy*dy
            sage: s == t
            False

        """
        resu = self._new_instance()
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = rst.copy()
        return resu

    def _common_subdomains(self, other):
        r"""
        Return the list of subdomains of ``self._domain`` on which
        both ``self`` and ``other`` have known restrictions.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 0]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: sorted(t._common_subdomains(t), key=str)
            [Open subset U of the 2-dimensional differentiable manifold M,
             Open subset V of the 2-dimensional differentiable manifold M,
             Open subset W of the 2-dimensional differentiable manifold M]
            sage: a = M.tensor_field(1, 1, name='a')
            sage: t._common_subdomains(a)
            []
            sage: a[e_xy, 0, 1] = 0
            sage: t._common_subdomains(a)
            [Open subset U of the 2-dimensional differentiable manifold M]
            sage: a[e_uv, 0, 0] = 0
            sage: sorted(t._common_subdomains(a), key=str)
            [Open subset U of the 2-dimensional differentiable manifold M,
             Open subset V of the 2-dimensional differentiable manifold M]

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

        - ``other`` -- a section or 0

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        TESTS::

            TODO

        """
        if other is self:
            return True
        if other in ZZ: # to compare with 0
            if other == 0:
                return self.is_zero()
            return False
        elif not isinstance(other, Section):
            return False
        else: # other is another tensor field
            if other._smodule != self._smodule:
                return False
            # Non-trivial open covers of the domain:
            open_covers = self._domain.open_covers()[1:]  # the open cover 0
                                                          # is trivial
            for oc in open_covers:
                resu = True
                for dom in oc:
                    try:
                        resu = resu and \
                                bool(self.restrict(dom) == other.restrict(dom))
                    except ValueError:
                        break
                else:
                    # If this point is reached, no exception has occured; hence
                    # the result is valid and can be returned:
                    return resu
            # If this point is reached, the comparison has not been possible
            # on any open cover; we then compare the restrictions to
            # subdomains:
            if not self._restrictions:
                return False  # self is not initialized
            if len(self._restrictions) != len(other._restrictions):
                return False  # the restrictions are not on the same subdomains
            resu = True
            for dom, rst in self._restrictions.items():
                if dom in other._restrictions:
                    resu = resu and bool(rst == other._restrictions[dom])
                else:
                    return False  # the restrictions are not on the same
                                  # subdomains
            return resu

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- a tensor field or 0

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[e_xy,:] = [[x+y, 0], [2, 1-y]]
            sage: t.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: t != t
            False
            sage: t != t.copy()
            False
            sage: t != 0
            True

        """
        return not (self == other)

    ######### ModuleElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Section addition.

        INPUT:

        - ``other`` -- a section, in the same section module as ``self``

        OUTPUT:

        - the section resulting from the addition of ``self`` and ``other``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.tensor_field(1, 1, name='b')
            sage: b[e_xy,:] = [[2, y], [x, -x]]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a._add_(b); s
            Tensor field a+b of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a.display(e_xy)
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: b.display(e_xy)
            b = 2 d/dx*dx + y d/dx*dy + x d/dy*dx - x d/dy*dy
            sage: s.display(e_xy)
            a+b = (x + 2) d/dx*dx + (y + 1) d/dx*dy + (x + y) d/dy*dx - x d/dy*dy
            sage: s == a + b  # indirect doctest
            True
            sage: z = a.parent().zero(); z
            Tensor field zero of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a._add_(z) == a
            True
            sage: z._add_(a) == a
            True

        """
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] + other._restrictions[dom]
        resu = self._vbundle.section(domain=self._domain)
        resu._restrictions = resu_rst
        if self._name is not None and other._name is not None:
            resu._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            resu._latex_name = self._latex_name + '+' + other._latex_name
        return resu

    def _sub_(self, other):
        r"""
        Section subtraction.

        INPUT:

        - ``other`` -- a section in the same section module as ``self``

        OUTPUT:

        - the section resulting from the subtraction of ``other`` from ``self``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: e_xy = c_xy.frame(); e_uv = c_uv.frame()
            sage: a = M.tensor_field(1, 1, name='a')
            sage: a[e_xy,:] = [[x, 1], [y, 0]]
            sage: a.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: b = M.tensor_field(1, 1, name='b')
            sage: b[e_xy,:] = [[2, y], [x, -x]]
            sage: b.add_comp_by_continuation(e_uv, U.intersection(V), c_uv)
            sage: s = a._sub_(b); s
            Tensor field a-b of type (1,1) on the 2-dimensional differentiable
             manifold M
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
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self._restrictions[dom] - other._restrictions[dom]
        resu = self._vbundle.section(domain=self._domain)
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

        - ``scalar`` -- scalar field in the scalar field algebra over which
          the module containing ``self`` is defined

        OUTPUT:

        - the tensor field ``scalar * self``

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:                    intersection_name='W', restrictions1= x>0,
            ....:                    restrictions2= u+v>0)
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
            on V: (u, v) |--> 2/(u^2 + v^2 + 2)
            sage: s = a._rmul_(f); s
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: a.display(e_xy)
            a = x d/dx*dx + d/dx*dy + y d/dy*dx
            sage: s.display(e_xy)
            x/(x^2 + y^2 + 1) d/dx*dx + 1/(x^2 + y^2 + 1) d/dx*dy + y/(x^2 + y^2 + 1) d/dy*dx
            sage: a.display(e_uv)
            a = (1/2*u + 1/2) d/du*du + (1/2*u - 1/2) d/du*dv + (1/2*v + 1/2) d/dv*du + (1/2*v - 1/2) d/dv*dv
            sage: s.display(e_uv)
            (u + 1)/(u^2 + v^2 + 2) d/du*du + (u - 1)/(u^2 + v^2 + 2) d/du*dv + (v + 1)/(u^2 + v^2 + 2) d/dv*du + (v - 1)/(u^2 + v^2 + 2) d/dv*dv
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
        for dom, rst in self._restrictions.items():
            resu._restrictions[dom] = scalar.restrict(dom) * rst
        return resu

    ######### End of ModuleElement arithmetic operators ########

#******************************************************************************

class TrivialSection(FiniteRankFreeModuleElement, Section):
    r"""
    Trivial sections
    """
    def __init__(self, section_module, name=None, latex_name=None):
        r"""

        """
        FiniteRankFreeModuleElement.__init__(self, section_module,
                                             name=name, latex_name=latex_name)
        self._domain = section_module.domain()
        self._vbundle = section_module.vector_bundle()
        self._base_space = section_module.base_space()
        # Initialization of derived quantities:
        self._init_derived()

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._init_derived()

        """
        FiniteRankFreeModuleElement._init_derived(self)
        Section._init_derived(self)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities.

        INPUT:

        - ``del_restrictions`` -- (default: ``True``) determines whether the
          restrictions of ``self`` to subdomains are deleted

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._del_derived()

        """
        FiniteRankFreeModuleElement._del_derived(self)
        Section._del_derived(self, del_restrictions=del_restrictions)

    def _repr_(self) :
        r"""
        String representation of ``self``.

        TESTS::

            TODO

        """
        return Section._repr_(self)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same
        vector field module, with the same tensor type and same symmetries.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: t = M.tensor_field(1,1, name='t')
            sage: t._new_instance()
            Tensor field of type (1,1) on the 2-dimensional differentiable
             manifold M
            sage: type(t._new_instance()) is type(t)
            True

        """
        return type(self)(self._fmodule)

    def set_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame
        for assignment.

        The components with respect to other frames on the same domain are
        deleted, in order to avoid any inconsistency. To keep them, use the
        method :meth:`add_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if none is provided, the components are
          assumed to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        EXAMPLES::

            sage: M = Manifold(2, 'M')
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

        The components with respect to the frame ``e_xy`` have be erased::

            sage: t.display(e_xy)
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components
             in the Coordinate frame (M, (d/dx,d/dy))

        Setting components in a frame defined on a subdomain deletes
        previously defined components as well::

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
            # Setting components on the section domain:
            return FiniteRankFreeModuleElement.set_comp(self, basis=basis)

        # Setting components on a subdomain:
        #
        # Creating or saving the restriction to the subdomain:
        rst = self.restrict(basis._domain)
        # Deleting all the components on self._domain and the derived
        # quantities:
        self._components.clear()
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
        kept. To delete them, use the method :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) vector frame in which the
          components are defined; if none is provided, the components are
          assumed to refer to the tensor field domain's default frame

        OUTPUT:

        - components in the given frame, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such
          components did not exist previously, they are created

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e_xy = X.frame()
            sage: t = M.tensor_field(1,1, name='t')
            sage: t.add_comp(e_xy)
            2-indices components w.r.t. Coordinate frame (M, (d/dx,d/dy))
            sage: t.add_comp(e_xy)[1,0] = 2
            sage: t.display(e_xy)
            t = 2 d/dy*dx

        Adding components with respect to a new frame (``e``)::

            sage: e = M.vector_frame('e')
            sage: t.add_comp(e)
            2-indices components w.r.t. Vector frame (M, (e_0,e_1))
            sage: t.add_comp(e)[0,1] = x
            sage: t.display(e)
            t = x e_0*e^1

        The components with respect to the frame ``e_xy`` are kept::

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
            # We perform a backup of the restrictions, since
            # they are deleted by FreeModuleTensor.add_comp (which
            # invokes del_derived()), and restore them afterwards
            restrictions_save = self._restrictions.copy()
            comp = FiniteRankFreeModuleElement.add_comp(self, basis=basis)
            self._restrictions = restrictions_save
            return comp

        # Adding components on a subdomain:
        #
        # Creating or saving the restriction to the subdomain:
        rst = self.restrict(basis._domain)
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

            sage: M = Manifold(2, 'M', start_index=1)
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
            # components on the local section domain:
            return FiniteRankFreeModuleElement.comp(self, basis=basis,
                                         from_basis=from_basis)

        # components on a subdomain:
        rst = self.restrict(basis._domain)
        return rst.comp(basis=basis, from_basis=from_basis)

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some subdomain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` --
          :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
          open subset `U` of the tensor field domain `S`
        - ``dest_map`` --
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          (default: ``None``); destination map
          `\Psi:\ U \rightarrow V`, where `V` is an open subset
          of the manifold `M` where the tensor field takes it values;
          if ``None``, the restriction of `\Phi` to `U` is used, `\Phi`
          being the differentiable map `S \rightarrow M` associated
          with the tensor field

        OUTPUT:

        - instance of :class:`TensorFieldParal` representing the restriction

        EXAMPLES:

        Restriction of a vector field defined on `\RR^2` to a disk::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: v = M.vector_field(x+y, -1+x^2, name='v')
            sage: D = M.open_subset('D') # the unit open disc
            sage: c_cart_D = c_cart.restrict(D, x^2+y^2<1)
            sage: v_D = v.restrict(D) ; v_D
            Vector field v on the Open subset D of the 2-dimensional
             differentiable manifold R^2
            sage: v_D.display()
            v = (x + y) d/dx + (x^2 - 1) d/dy

        The symbolic expressions of the components with respect to
        Cartesian coordinates are equal::

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

        The restriction of the vector field to its own domain is of
        course itself::

            sage: v.restrict(M) is v
            True

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subset of " +
                                 "the field's domain")
            # First one tries to derive the restriction from a tighter domain:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and subdomain in rst._restrictions:
                    res = rst._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and dom is not self._domain:
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    self._restrictions_graph[subdomain] = rst.restrict(subdomain)
                    return self._restrictions[subdomain]

            # Secondly one tries to get the restriction from one previously
            # defined on a larger domain:
            for dom, ext in self._extensions_graph.items():
                if subdomain in ext._restrictions_graph:
                    res = ext._restrictions_graph[subdomain]
                    self._restrictions[subdomain] = res
                    self._restrictions_graph[subdomain] = res
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions_graph[subdomain] = res
                    return self._restrictions[subdomain]

            # If this fails, the restriction is created from scratch:
            smodule = self._vbundle.section_module(domain=subdomain)
            res = smodule.element_class(smodule, name=self._name,
                                        latex_name=self._latex_name)

            for frame in self._components:
                for sframe in self._vbundle._frames:
                    if (sframe.domain() is subdomain and
                            sframe in frame._subframes):
                        comp_store = self._components[frame]._comp
                        scomp = res._new_comp(sframe)
                        scomp_store = scomp._comp
                        # the components of the restriction are evaluated
                        # index by index:
                        for ind, value in comp_store.items():
                            scomp_store[ind] = value.restrict(subdomain)
                        res._components[sframe] = scomp

            res._extensions_graph.update(self._extensions_graph)
            for dom, ext in self._extensions_graph.items():
                ext._restrictions[subdomain] = res
                ext._restrictions_graph[subdomain] = res

            for dom, rst in self._restrictions.items():
                if dom.is_subset(subdomain):
                    if rst is not res:
                        res._restrictions.update(rst._restrictions)
                    res._restrictions_graph.update(rst._restrictions_graph)
                    rst._extensions_graph.update(res._extensions_graph)

            self._restrictions[subdomain] = res
            self._restrictions_graph[subdomain] = res

        return self._restrictions[subdomain]

    def display_comp(self, frame=None, chart=None, only_nonzero=True):
        r"""
        Display the tensor components with respect to a given frame,
        one per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) vector frame with respect to which
          the tensor field components are defined; if ``None``, then

          * if ``chart`` is not ``None``, the coordinate frame associated to
            ``chart`` is used
          * otherwise, the default basis of the vector field module on which
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

        Display of the components of a type-`(2,1)` tensor field on a
        2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(2, 1, name='t', sym=(0,1))
            sage: t[0,0,0], t[0,1,0], t[1,1,1] = x+y, x*y, -3
            sage: t.display_comp()
            t^xx_x = x + y
            t^xy_x = x*y
            t^yx_x = x*y
            t^yy_y = -3

        By default, only the non-vanishing components are displayed;
        to see all the components, the argument ``only_nonzero`` must
        be set to ``False``::

            sage: t.display_comp(only_nonzero=False)
            t^xx_x = x + y
            t^xx_y = 0
            t^xy_x = x*y
            t^xy_y = 0
            t^yx_x = x*y
            t^yx_y = 0
            t^yy_x = 0
            t^yy_y = -3

        ``t`` being symmetric with respect to its first two indices, one
        may ask to skip the components that can be deduced by symmetry::

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

        Display in a frame different from the default one (note that
        since ``f`` is not a coordinate frame, integer are used to
        label the indices)::

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
                frame = self._fmodule.default_basis()
        if chart is None:
            chart = self._domain.default_chart()
        return FiniteRankFreeModuleElement.display_comp(self, basis=frame,
                                  format_spec=chart,
                                  only_nonzero=only_nonzero)

    def at(self, point):
        r"""
        Value of ``self`` at a point of its domain.

        If the current tensor field is

        .. MATH::

            t:\ U  \longrightarrow T^{(k,l)} M

        associated with the differentiable map

        .. MATH::

            \Phi:\ U \longrightarrow M,

        where `U` and `M` are two manifolds (possibly `U = M` and
        `\Phi = \mathrm{Id}_M`), then for any point `p\in U`, `t(p)` is
        a tensor on the tangent space to `M` at the point `\Phi(p)`.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`
          point `p` in the domain of the tensor field `U`

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor`
          representing the tensor `t(p)` on the tangent vector space
          `T_{\Phi(p)} M`

        EXAMPLES:

        Vector in a tangent space of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: p = M.point((-2,3), name='p')
            sage: v = M.vector_field(y, x^2, name='v')
            sage: v.display()
            v = y d/dx + x^2 d/dy
            sage: vp = v.at(p) ; vp
            Tangent vector v at Point p on the 2-dimensional differentiable
             manifold M
            sage: vp.parent()
            Tangent space at Point p on the 2-dimensional differentiable
             manifold M
            sage: vp.display()
            v = 3 d/dx + 4 d/dy

        A 1-form gives birth to a linear form in the tangent space::

            sage: w = M.one_form(-x, 1+y, name='w')
            sage: w.display()
            w = -x dx + (y + 1) dy
            sage: wp = w.at(p) ; wp
            Linear form w on the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: wp.parent()
            Dual of the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: wp.display()
            w = 2 dx + 4 dy

        A tensor field of type `(1,1)` yields a tensor of type `(1,1)`
        in the tangent space::

            sage: t = M.tensor_field(1, 1, name='t')
            sage: t[0,0], t[0,1], t[1,1] = 1+x, x*y, 1-y
            sage: t.display()
            t = (x + 1) d/dx*dx + x*y d/dx*dy + (-y + 1) d/dy*dy
            sage: tp = t.at(p) ; tp
            Type-(1,1) tensor t on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: tp.parent()
            Free module of type-(1,1) tensors on the Tangent space at Point p
             on the 2-dimensional differentiable manifold M
            sage: tp.display()
            t = -d/dx*dx - 6 d/dx*dy - 2 d/dy*dy

        A 2-form yields an alternating form of degree 2 in the tangent space::

            sage: a = M.diff_form(2, name='a')
            sage: a[0,1] = x*y
            sage: a.display()
            a = x*y dx/\dy
            sage: ap = a.at(p) ; ap
            Alternating form a of degree 2 on the Tangent space at Point p on
             the 2-dimensional differentiable manifold M
            sage: ap.parent()
            2nd exterior power of the dual of the Tangent space at Point p on
             the 2-dimensional differentiable manifold M
            sage: ap.display()
            a = -6 dx/\dy

        Example with a non trivial map `\Phi`::

            sage: U = Manifold(1, 'U')  # (0,2*pi) as a 1-dimensional manifold
            sage: T.<t> = U.chart(r't:(0,2*pi)')  # canonical chart on U
            sage: Phi = U.diff_map(M, [cos(t), sin(t)], name='Phi',
            ....:                  latex_name=r'\Phi')
            sage: v = U.vector_field(1+t, t^2, name='v', dest_map=Phi) ; v
            Vector field v along the 1-dimensional differentiable manifold U
             with values on the 2-dimensional differentiable manifold M
            sage: v.display()
            v = (t + 1) d/dx + t^2 d/dy
            sage: p = U((pi/6,))
            sage: vp = v.at(p) ; vp
            Tangent vector v at Point on the 2-dimensional differentiable
             manifold M
            sage: vp.parent() is M.tangent_space(Phi(p))
            True
            sage: vp.display()
            v = (1/6*pi + 1) d/dx + 1/36*pi^2 d/dy

        """
        if point not in self._domain:
            raise ValueError("the {} is not in the domain of ".format(point) +
                             "the {}".format(self))
        vbf = self._vbundle.fiber(point)
        resu = vbf.tensor((1,0), name=self._name,
                         latex_name=self._latex_name)
        for frame, comp in self._components.items():
            comp_resu = resu.add_comp(frame.at(point))
            for ind, val in comp._comp.items():
                comp_resu._comp[ind] = val(point)
        return resu
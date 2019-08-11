r"""
Local Frames

The class :class:`LocalFrame` implements local frames on topological vector
bundles (see :class:`sage.manifolds.vector_bundle.TopologicalVectorBundle`).

A *local frame* on a topological vector bundle `E \to M` is a continuous local
section `(e_1,\dots,e_n):U \to E defined on some subset `U` of the base space
`M`, such that `e(p)` is a basis of the fiber `E_p` for any `p \in U`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Travis Scrimshaw (2016): review tweaks
- Eric Gourgoulhon (2018): some refactoring and more functionalities in the
  choice of symbols for vector frame elements (:trac:`24792`)
- Michael Jung (2019): Generalization to vector bundles

EXAMPLES:

    TODO

"""

#******************************************************************************
#       Copyright (C) 2013-2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_basis import (FreeModuleBasis,
                                                   FreeModuleCoBasis)
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule

class LocalCoFrame(FreeModuleCoBasis):
    r"""
    Local coframe on a vector bundle.

    A *local coframe* on a topological vector bundle `E \to M` is a continuous
    local section `e^*: U \to E` on some subset `U` of the base space `M`, such
    that `e^*(p)` is a basis of the fiber `E^*_p` of the dual bundle for any
    `p \in U`.

    INPUT:

    - ``frame`` -- the local frame dual to the coframe
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the 1-forms constituting the coframe, or a tuple of strings,
      representing the individual symbols of the 1-forms
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the 1-forms constituting the
      coframe, or a tuple of strings, representing the individual LaTeX symbols
      of the 1-forms; if ``None``, ``symbol`` is used in place of
      ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the 1-forms
      of the coframe; if ``None``, the indices will be generated as integers
      within the range declared on the coframe's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the 1-forms of the coframe; if
      ``None``, ``indices`` is used instead

    EXAMPLES:

        TODO

    """
    def __init__(self, frame, symbol, latex_symbol=None, indices=None,
                 latex_indices=None):
        r"""
        Construct a local coframe, dual to a given local frame.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: from sage.manifolds.local_frame import LocalCoFrame
            sage: f = LocalCoFrame(e, 'f'); f
            Local coframe (E|_M, (f^0,f^1))
            sage: TestSuite(f).run()

        """
        self._domain = frame.domain()
        self._base_space = frame.base_space()
        self._vbundle = frame.vector_bundle()
        FreeModuleCoBasis.__init__(self, frame, symbol,
                                   latex_symbol=latex_symbol, indices=indices,
                                   latex_indices=latex_indices)
        # The coframe is added to the vector bundle's set of coframes
        self._vbundle._coframes.append(self)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: f = e.coframe()
            sage: f._repr_()
            'Local coframe (E|_M, (e^0,e^1))'
            sage: repr(f)  # indirect doctest
            'Local coframe (E|_M, (e^0,e^1))'
            sage: f  # indirect doctest
            Local coframe (E|_M, (e^0,e^1))

        """
        desc = "Local coframe " + self._name
        return desc

    def at(self, point):
        r"""
        Return the value of ``self`` at a given point on the base space, this
        value being a basis of the dual of the vector bundle at the point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` in the domain `U` of the coframe (denoted `f` hereafter)

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_basis.FreeModuleCoBasis`
          representing the basis `f(p)` of the vector space `E^*_p`,
          dual to the vector bundle fiber `E_p`

        EXAMPLES:

            TODO

        """
        return self._basis.at(point).dual_basis()

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='up',
                 include_domain=True):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the 1-forms constituting the coframe, or a list/tuple of
          strings, representing the individual symbols of the 1-forms
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the 1-forms constituting
          the coframe, or a list/tuple of strings, representing the individual
          LaTeX symbols of the 1-forms; if ``None``, ``symbol`` is used in
          place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the 1-forms of the coframe; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the 1-forms;
          if ``None``, ``indices`` is used instead
        - ``index_position`` -- (default: ``'up'``) determines the position
          of the indices labelling the 1-forms of the coframe; can be
          either ``'down'`` or ``'up'``
        - ``include_domain`` -- (default: ``True``) boolean determining whether
          the name of the domain is included in the beginning of the coframe
          name

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e').coframe(); e
            Local coframe (E|_M, (e^0,e^1))
            sage: e.set_name('f'); e
            Local coframe (E|_M, (f^0,f^1))
            sage: e.set_name('e', latex_symbol=r'\epsilon')
            sage: latex(e)
            \left(E|_{M}, \left(\epsilon^{0},\epsilon^{1}\right)\right)
            sage: e.set_name('e', include_domain=False); e
            Local coframe (e^0,e^1)
            sage: e.set_name(['a', 'b'], latex_symbol=[r'\alpha', r'\beta']); e
            Local coframe (E|_M, (a,b))
            sage: latex(e)
            \left(E|_{M}, \left(\alpha,\beta\right)\right)
            sage: e.set_name('e', indices=['x','y'],
            ....:            latex_indices=[r'\xi', r'\zeta']); e
            Local coframe (E|_M, (e^x,e^y))
            sage: latex(e)
            \left(E|_{M}, \left(e^{\xi},e^{\zeta}\right)\right)

        """
        super(LocalCoFrame, self).set_name(symbol, latex_symbol=latex_symbol,
                                      indices=indices,
                                      latex_indices=latex_indices,
                                      index_position=index_position)
        if include_domain:
            # Redefinition of the name and the LaTeX name to include the domain
            self._name = "({}|_{}, {})".format(self._vbundle._name,
                                                self._domain._name, self._name)
            self._latex_name = r"\left({}|_{{{}}}, {}\right)".format(
                                                    self._vbundle._latex_name,
                                                    self._domain._latex_name,
                                                    self._latex_name)

#******************************************************************************

class LocalFrame(FreeModuleBasis):
    r"""
    Local frame on a vector bundle.

    A *local frame* on a topological vector bundle `E \to M` is a continuous
    local section `(e_1,\dots,e_n):U \to E defined on some subset `U` of the
    base space `M`, such that `e(p)` is a basis of the fiber `E_p` for any
    `p \in U`.

    For each instanciation of a local frame, a local coframe is automatically
    created, as an instance of the class :class:`LocalCoFrame`. It is returned
    by the method :meth:`coframe`.

    INPUT:

    - ``section_module`` -- free module of local sections over `U` in the given
      vector bundle `E \to M`
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the local sections constituting the local frame, or a tuple
      of strings, representing the individual symbols of the local sections
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the local sections constituting
      the local frame, or a tuple of strings, representing the individual
      LaTeX symbols of the local sections; if ``None``, ``symbol`` is used in
      place of ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the local
      sections of the frame; if ``None``, the indices will be generated as
      integers within the range declared on the local frame's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the  local sections; if
      ``None``, ``indices`` is used instead
    - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
      dual coframe; if ``None``, ``symbol`` must be a string and is used
      for the common base of the symbols of the elements of the dual coframe
    - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
      but for the dual coframe

    EXAMPLES:

    Defining a local frame on a 3-dimensional vector bundle over a 3-dimensional
    manifold::

        sage: M = Manifold(3, 'M', start_index=1, structure='top')
        sage: E = M.vector_bundle(3, 'E')
        sage: e = E.local_frame('e') ; e
        Local frame (E|_M, (e_1,e_2,e_3))
        sage: latex(e)
        \left(E|_{M}, \left(e_{1},e_{2},e_{3}\right)\right)

    TODO

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = LocalCoFrame

    @staticmethod
    def __classcall_private__(cls, section_module, symbol,
                              latex_symbol=None, indices=None,
                              latex_indices=None, symbol_dual=None,
                              latex_symbol_dual=None):
        """
        Transform input lists into tuples for the unique representation of
        LocalFrame.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True)
            sage: from sage.manifolds.local_frame import LocalFrame
            sage: e = LocalFrame(C0, ['a', 'b'], symbol_dual=['A', 'B']); e
            Local frame (E|_M, (a,b))
            sage: e.dual_basis()
            Local coframe (E|_M, (A,B))
            sage: e is LocalFrame(C0, ('a', 'b'), symbol_dual=('A', 'B'))
            True

        """
        if isinstance(symbol, list):
            symbol = tuple(symbol)
        if isinstance(latex_symbol, list):
            latex_symbol = tuple(latex_symbol)
        if isinstance(indices, list):
            indices = tuple(indices)
        if isinstance(latex_indices, list):
            latex_indices = tuple(latex_indices)
        if isinstance(symbol_dual, list):
            symbol_dual = tuple(symbol_dual)
        if isinstance(latex_symbol_dual, list):
            latex_symbol_dual = tuple(latex_symbol_dual)
        return super(LocalFrame, cls).__classcall__(cls, section_module,
                                        symbol, latex_symbol=latex_symbol,
                                        indices=indices,
                                        latex_indices=latex_indices,
                                        symbol_dual=symbol_dual,
                                        latex_symbol_dual=latex_symbol_dual)

    def __init__(self, section_module, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, symbol_dual=None, latex_symbol_dual=None):
        r"""
        Construct a local frame on a vector bundle.

        TESTS:

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True)
            sage: from sage.manifolds.local_frame import LocalFrame
            sage: e = LocalFrame(C0, 'e', latex_symbol=r'\epsilon'); e
            Local frame (E|_M, (e_0,e_1))
            sage: TestSuite(e).run()

        """
        from sage.tensor.modules.finite_rank_free_module import \
            FiniteRankFreeModule
        if not isinstance(section_module, FiniteRankFreeModule):
            raise ValueError("'section_module' must be a free module")
        self._domain = section_module.domain()
        self._base_space = section_module.base_space()
        self._vbundle = section_module.vector_bundle()
        FreeModuleBasis.__init__(self, section_module,
                                 symbol, latex_symbol=latex_symbol,
                                 indices=indices, latex_indices=latex_indices,
                                 symbol_dual=symbol_dual,
                                 latex_symbol_dual=latex_symbol_dual)
        self._coframe = self.dual_basis()  # Shortcut for self._dual_basis
        ###
        # Add this frame to the list of frames of the overlying vector bundle:
        self._vbundle._add_local_frame(self)
        ###
        # Frame restrictions:
        self._subframes = set([self]) # Set of frames which are just a
                        # restriction of self
        self._superframes = set([self]) # Set of frames for which self is a
                        # restriction of
        self._restrictions = {} # Key: subdomain of self._domain; value:
                        # restriction of self on this subdomain

    ###### Methods that must be redefined by derived classes of ######
    ###### FreeModuleBasis                                      ######

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.local_frame('e')
            'Local frame (E|_M, (e_0,e_1))'
            sage: repr(e)  # indirect doctest
            'Local frame (E|_M, , (e_0,e_1))'
            sage: e  # indirect doctest
            Local frame (E|_M, , (e_0,e_1))

        """
        desc = "Local frame " + self._name
        return desc

    def _new_instance(self, symbol, latex_symbol=None, indices=None,
                      latex_indices=None, symbol_dual=None,
                      latex_symbol_dual=None):
        r"""
        Construct a new local frame on the same section module as ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the sections constituting the local frame, or a
          tuple of strings, representing the individual symbols of the local
          sections
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the local sections
          constituting the local frame, or a tuple of strings, representing
          the individual LaTeX symbols of the local sections; if ``None``,
          ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the local sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on the local frame's
          domain
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the local sections;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe

        OUTPUT:

        - instance of :class:`LocalFrame`

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: e._new_instance('f')
            Local frame (E|_M, (f_0,f_1))

        """
        return LocalFrame(self._fmodule, symbol, latex_symbol=latex_symbol,
                           indices=indices, latex_indices=latex_indices,
                           symbol_dual=symbol_dual,
                           latex_symbol_dual=latex_symbol_dual)

    ###### End of methods to be redefined by derived classes ######

    def domain(self):
        r"""
        Return the domain on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U); e
            Local frame (E|_U, (e_0,e_1))
            sage: e.domain()
            Open subset U of the 3-dimensional topological manifold M

        """
        return self._domain

    def base_space(self):
        r"""
        Return the base space on which the overlying vector bundle is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U)
            sage: e.base_space()
            3-dimensional topological manifold M

        """
        return self._base_space

    def vector_bundle(self):
        r"""
        Return the vector bundle on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U)
            sage: e.vector_bundle()
            Topological real vector bundle E -> M of rank 2 over the base space
            3-dimensional topological manifold M
            sage: e.vector_bundle() is E
            True

        """
        return self._vbundle

    def coframe(self):
        r"""
        Return the coframe of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: e.coframe()
            Local coframe (E|_M, (e^0,e^1))

        """
        return self._coframe

    def new_frame(self, change_of_frame, symbol, latex_symbol=None,
                  indices=None, latex_indices=None, symbol_dual=None,
                  latex_symbol_dual=None):
        r"""
        Define a new local frame from ``self``.

        The new local frame is defined from vector bundle automorphisms; its
        module is the same as that of the current frame.

        INPUT:

        - ``change_of_frame`` --
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`;
          vector bundle automorphisms `P` that relates
          the current frame `(e_i)` to the new frame `(f_i)` according
          to `f_i = P(e_i)`
        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the sections constituting the local frame, or a
          list/tuple of strings, representing the individual symbols of the
          sections
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the sections
          constituting the local frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the sections;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the sections;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe

        OUTPUT:

        - the new frame `(f_i)`, as an instance of :class:`LocalFrame`

        EXAMPLES:

            TODO

        """
        the_new_frame = self.new_basis(change_of_frame, symbol,
                                       latex_symbol=latex_symbol,
                                       indices=indices,
                                       latex_indices=latex_indices,
                                       symbol_dual=symbol_dual,
                                       latex_symbol_dual=latex_symbol_dual)
        self._vbundle._frame_changes[(self, the_new_frame)] = \
                            self._fmodule._basis_changes[(self, the_new_frame)]
        self._vbundle._frame_changes[(the_new_frame, self)] = \
                            self._fmodule._basis_changes[(the_new_frame, self)]
        return the_new_frame

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some open subset of its domain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `V` of the current frame domain `U`

        OUTPUT:

        - the restriction of the current frame to `V` as a :class:`LocalFrame`

        EXAMPLES:

            TODO

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subdomain of " +
                                 "the current frame's domain")
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and subdomain in rst._restrictions:
                    res = rst._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    res._superframes.update(self._superframes)
                    for sframe2 in self._superframes:
                        sframe2._subframes.add(res)
                    return self._restrictions[subdomain]
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom):
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    return self._restrictions[subdomain]
            # Secondly one tries to get the restriction from one previously
            # defined on a larger domain:
            for sframe in self._superframes:
                if subdomain in sframe._restrictions:
                    res = sframe._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    res._superframes.update(self._superframes)
                    for sframe2 in self._superframes:
                        sframe2._subframes.add(res)
                    return self._restrictions[subdomain]
            # If this point is reached, the restriction has to be created
            # from scratch:
            resmodule = self._vbundle.section_module(subdomain, force_free=True)
            res = LocalFrame(resmodule,
                             self._symbol, latex_symbol=self._latex_symbol,
                             indices=self._indices,
                             latex_indices=self._latex_indices,
                             symbol_dual=self._symbol_dual,
                             latex_symbol_dual=self._latex_symbol_dual)
            #new_vectors = list()
            #for i in self._fmodule.irange():
            #    vrest = self[i].restrict(subdomain)
            #    for j in self._fmodule.irange():
            #        vrest.add_comp(res)[j] = 0
            #    vrest.add_comp(res)[i] = 1
            #    new_vectors.append(vrest)
            #res._vec = tuple(new_vectors)
            # Update of superframes and subframes:
            res._superframes.update(self._superframes)
            for sframe in self._superframes:
                sframe._subframes.add(res)
                sframe._restrictions[subdomain] = res # includes sframe = self
            for dom, rst in self._restrictions.items():
                if dom.is_subset(subdomain):
                    res._restrictions.update(rst._restrictions)
                    res._subframes.update(rst._subframes)
                    rst._superframes.update(res._superframes)
        return self._restrictions[subdomain]

    def at(self, point):
        r"""
        Return the value of ``self`` at a given point, this value being
        a basis of the vector bundle fiber at the point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`; point
          `p` in the domain `U` of the local frame (denoted `e` hereafter)

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing the basis `e(p)` of the vector bundle fiber
          `E_p`

        EXAMPLES:

            TODO

        """
        # Determination of the vector bundle fiber:
        if point not in self._domain:
            raise ValueError("the {} is not a point in the ".format(point) +
                             "domain of {}".format(self))
        vbf = self._vbundle.fiber(point)
        # If the basis has already been constructed, it is simply returned:
        vbf_frame_bases = vbf._frame_bases
        if self in vbf_frame_bases:
            return vbf_frame_bases[self]
        for frame in vbf_frame_bases:
            if self in frame._subframes or self in frame._superframes:
                return vbf_frame_bases[frame]
        # If this point is reached, the basis has to be constructed from
        # scratch.
        # The names of the basis vectors set to those of the frame vector
        # fields:
        basis = vbf.basis(self._symbol, latex_symbol=self._latex_symbol,
                         indices=self._indices,
                         latex_indices=self._latex_indices,
                         symbol_dual=self._symbol_dual,
                         latex_symbol_dual=self._latex_symbol_dual)
        vbf_frame_bases[self] = basis
        # Update of the change of bases in the fiber:
        for frame_pair, automorph in self._vbundle._frame_changes.items():
            frame1 = frame_pair[0]; frame2 = frame_pair[1]
            if frame1 is self:
                fr2 = None
                for frame in vbf_frame_bases:
                    if frame2 in frame._subframes:
                        fr2 = frame
                        break
                if fr2 is not None:
                    basis1 = basis
                    basis2 = vbf_frame_bases[fr2]
                    auto = vbf.automorphism()
                    for frame, comp in automorph._components.items():
                        bas = None
                        if frame is frame1:
                            bas = basis1
                        if frame is frame2:
                            bas = basis2
                        if bas is not None:
                            cauto = auto.add_comp(bas)
                            for ind, val in comp._comp.items():
                                cauto._comp[ind] = val(point)
                    vbf._basis_changes[(basis1, basis2)] = auto
            if frame2 is self:
                fr1 = None
                for frame in vbf_frame_bases:
                    if frame1 in frame._subframes:
                        fr1 = frame
                        break
                if fr1 is not None:
                    basis1 = vbf_frame_bases[fr1]
                    basis2 = basis
                    auto = vbf.automorphism()
                    for frame, comp in automorph._components.items():
                        bas = None
                        if frame is frame1:
                            bas = basis1
                        if frame is frame2:
                            bas = basis2
                        if bas is not None:
                            cauto = auto.add_comp(bas)
                            for ind, val in comp._comp.items():
                                cauto._comp[ind] = val(point)
                    vbf._basis_changes[(basis1, basis2)] = auto
        return basis

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='down',
                 include_domain=True):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the local sections constituting the local frame, or a
          list/tuple of strings, representing the individual symbols of the
          local sections
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the local sections
          constituting the local frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the local sections;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the local sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the local sections;
          if ``None``, ``indices`` is used instead
        - ``index_position`` -- (default: ``'down'``) determines the position
          of the indices labelling the local sections of the frame; can be
          either ``'down'`` or ``'up'``
        - ``include_domain`` -- (default: ``True``) boolean determining whether
          the name of the domain is included in the beginning of the vector
          frame name

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: e.set_name('f'); e
            Local frame (E|_M, (f_0,f_1))
            sage: e.set_name('e', include_domain=False); e
            Local frame (e_0,e_1)
            sage: e.set_name(['a', 'b']); e
            Local frame (E|_M, (a,b))
            sage: e.set_name('e', indices=['x', 'y']); e
            Local frame (E|_M, (e_x,e_y))
            sage: e.set_name('e', latex_symbol=r'\epsilon')
            sage: latex(e)
            \left(E|_{M}, \left(\epsilon_{0},\epsilon_{1}\right)\right)
            sage: e.set_name('e', latex_symbol=[r'\alpha', r'\beta'])
            sage: latex(e)
            \left(E|_{M}, \left(\alpha,\beta\right)\right)
            sage: e.set_name('e', latex_symbol='E',
            ....:            latex_indices=[r'\alpha', r'\beta'])
            sage: latex(e)
            \left(E|_{M}, \left(E_{\alpha},E_{\beta}\right)\right)

        """
        super(LocalFrame, self).set_name(symbol, latex_symbol=latex_symbol,
                                          indices=indices,
                                          latex_indices=latex_indices,
                                          index_position=index_position)
        if include_domain:
            # Redefinition of the name and the LaTeX name to include the domain
            self._name = "({}|_{}, {})".format(self._vbundle._name,
                                                self._domain._name, self._name)
            self._latex_name = r"\left({}|_{{{}}}, {}\right)".format(
                                                    self._vbundle._latex_name,
                                                    self._domain._latex_name,
                                                    self._latex_name)

#******************************************************************************

class TrivializationCoFrame(LocalCoFrame):
    r"""
    Trivialization coframe on a topological vector bundle.

    A *trivialization coframe* is the coframe of the trivialization frame
    induced by a trivialization (see: :class:`~sage.manifolds.local_frame.TrivializationFrame`).

    More precisely, a *trivialization frame* on a topological vector bundle
    `E \to M` of rank `n` over the topological field `K` and over a topological
    manifold `M` is a local coframe induced by a local trivialization
    `\varphi:E|_U \to U \times K^n` of the domain `U \in M`. Namely, the local
    dual sections

    .. MATH::

        \varphi^*e^i := \varphi(\;\cdot\;, e^i)

    on `U` induce a local frame `(\varphi^*e^1, \dots, \varphi^*e^n)`, where
    `(e^1, \dots, e^n)` is the dual of the standard basis of `K^n`.

    INPUT:

    - ``triv_frame`` -- trivialization frame dual to the trivialization coframe
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the dual sections constituting the coframe, or a tuple of
      strings, representing the individual symbols of the dual sections
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the dual sections constituting
      the coframe, or a tuple of strings, representing the individual LaTeX
      symbols of the dual sections; if ``None``, ``symbol`` is used in place of
      ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the dual
      sections of the coframe; if ``None``, the indices will be generated as
      integers within the range declared on the local frame's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the dual sections of the coframe; if
      ``None``, ``indices`` is used instead

    EXAMPLES:

    TODO

    """
    def __init__(self, triv_frame, symbol, latex_symbol=None,
                 indices=None, latex_indices=None):
        r"""
        Construct a local coframe from a local trivialization.

        TESTS::

            TODO

        """
        if not isinstance(triv_frame, TrivializationFrame):
            raise TypeError("the first argument must be a local trivialization "
                            "frame")
        LocalCoFrame.__init__(self, triv_frame, symbol,
                              latex_symbol=latex_symbol, indices=indices,
                              latex_indices=latex_indices)
        self._trivialization = triv_frame._trivialization

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: e = phi.frame().coframe()
            sage: e._repr_()
            'Trivialization coframe (E|_M, (phi^*e^1,phi^*e^2))'
            sage: repr(e)  # indirect doctest
            'Trivialization coframe (E|_M, (phi^*e^1,phi^*e^2))'
            sage: e  # indirect doctest
            Trivialization coframe (E|_M, (phi^*e^1,phi^*e^2))

        """
        return "Trivialization coframe " + self._name

#******************************************************************************

class TrivializationFrame(LocalFrame):
    r"""
    Trivialization frame on a topological vector bundle.

    A *trivialization frame* on a topological vector bundle `E \to M` of rank
    `n` over the topological field `K` and over a topological manifold `M` is a
    local frame induced by a local trivialization `\varphi:E|_U \to U \times K^n`
    of the domain `U \in M`. More precisely, the local sections

    .. MATH::

        \varphi^*e_i := \varphi(\;\cdot\;, e_i)

    on `U` induce a local frame `(\varphi^*e_1, \dots, \varphi^*e_n)`, where
    `(e_1, \dots, e_n)` is the standard basis of `K^n`.

    INPUT:

    - ``trivialization`` -- the trivialization defined on the vector bundle

    EXAMPLES::

        sage: M = Manifold(3, 'M')
        sage: U = M.open_subset('U')
        sage: E = M.vector_bundle(2, 'E')
        sage: phi_U = E.trivialization('phi_U', domain=U)
        sage: phi_U.frame()
        Trivialization frame (E|_U, (phi_U^*e_1,phi_U^*e_2))
        sage: latex(phi_U.frame())
        \left(E|_{U}, \left(\left(phi_U^* e_{ 1 }\right),\left(phi_U^* e_{ 2 }
         \right)\right)\right)

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = TrivializationCoFrame

    def __init__(self, trivialization):
        r"""
        Construct a trivialization frame.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: e = phi.frame()
            sage: TestSuite(e).run()

        """
        from sage.misc.latex import latex
        from .trivialization import Trivialization
        if not isinstance(trivialization, Trivialization):
            raise TypeError("the first argument must be a trivialization")
        triv = trivialization
        domain = triv.domain()
        self._trivialization = triv
        ###
        # Define trivialization names
        vbundle = triv.vector_bundle()
        rank = vbundle.rank()
        if triv._name is None:
            symbol = tuple("e_" + str(i) +
                "({},{})".format(triv.domain()._name, vbundle._name)
                for i in range(1, rank + 1))
            symbol_dual = tuple("e^" + str(i) +
                "({},{})".format(triv.domain()._name, vbundle._name)
                for i in range(1, rank + 1))
        else:
            symbol = tuple(triv._name + "^*" + "e_" + str(i)
                for i in range(1, rank + 1))
            symbol_dual = tuple(triv._name + "^*" + "e^" + str(i)
                for i in range(1, rank + 1))
        if triv._latex_name is None:
            latex_symbol = tuple(r'e_{' + latex(i) + r'}(' +
                                 triv.domain()._latex_name + r',' +
                                 vbundle._latex_name + r')'
                                 for i in range(1, rank + 1))
            latex_symbol_dual = tuple(r'e^{' + latex(i) + r'}(' +
                                 triv.domain()._latex_name + r',' +
                                 vbundle._latex_name + r')'
                                 for i in range(1, rank + 1))
        else:
            latex_symbol = tuple(r'\left(' + triv._latex_name + r'^* e_{' +
                                 latex(i) + r'}\right)'
                                 for i in range(1, rank + 1))
            latex_symbol_dual = tuple(r'\left(' + triv._latex_name + r'^* e^{' +
                                 latex(i) + r'}\right)'
                                 for i in range(1, rank + 1))
        LocalFrame.__init__(self,
                        vbundle.section_module(domain=domain, force_free=True),
                        symbol=symbol, latex_symbol=latex_symbol,
                        symbol_dual=symbol_dual,
                        latex_symbol_dual=latex_symbol_dual)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: e = phi.frame()
            sage: e._repr_()
            'Trivialization frame (E|_M, (phi^*e_1,phi^*e_2))'
            sage: repr(e)  # indirect doctest
            'Trivialization frame (E|_M, (phi^*e_1,phi^*e_2))'
            sage: e  # indirect doctest
            Trivialization frame (E|_M, (phi^*e_1,phi^*e_2))

        """
        return "Trivialization frame " + self._name

    def trivialization(self):
        r"""
        Return the underlying trivialization of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: e = phi_U.frame()
            sage: e.trivialization()
            Trivialization (phi_U:E|_M -> M)
            sage: e.trivialization() is phi_U
            True

        """
        return self._trivialization
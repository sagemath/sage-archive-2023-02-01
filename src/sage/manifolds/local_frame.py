r"""
Local Frames

AUTHORS:

- Michael Jung (2019): initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
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

    INPUT:

    - ``frame`` -- the vector frame dual to the coframe

    EXAMPLES:



    """
    def __init__(self, frame, symbol, latex_symbol=None, indices=None,
                 latex_indices=None):
        r"""
        Construct a local coframe, dual to a given local frame.

        TESTS::



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

        """
        desc = "Local coframe " + self._name + " "
        desc += "on {}".format(self._basis.module())
        return desc

    def at(self, point):
        r"""

        """
        return self._basis.at(point).dual_basis()

    def domain(self):
        r"""

        """
        return self._domain

    def base_space(self):
        r"""

        """
        return self._base_space

    def vector_bundle(self):
        r"""

        """
        return self._vbundle

#******************************************************************************

class LocalFrame(FreeModuleBasis):
    r"""

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = LocalCoFrame

    @staticmethod
    def __classcall_private__(cls, section_module, symbol,
                              latex_symbol=None, from_frame=None, indices=None,
                              latex_indices=None, symbol_dual=None,
                              latex_symbol_dual=None):
        """
        Transform input lists into tuples for the unique representation of
        LocalFrame.

        TESTS::



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
                                        from_frame=from_frame, indices=indices,
                                        latex_indices=latex_indices,
                                        symbol_dual=symbol_dual,
                                        latex_symbol_dual=latex_symbol_dual)

    def __init__(self, section_module, symbol, latex_symbol=None,
                 from_frame=None, indices=None, latex_indices=None,
                 symbol_dual=None, latex_symbol_dual=None):
        r"""

        """
        from sage.tensor.modules.finite_rank_free_module import \
            FiniteRankFreeModule
        if not isinstance(section_module, FiniteRankFreeModule):
            raise ValueError("'section_module' must be a free module")
        self._domain = section_module.domain()
        self._base_space = section_module.base_space()
        self._vbundle = section_module.vector_bundle()
        self._from_frame = from_frame
        if from_frame is not None:
            if not from_frame.domain().is_subset(self._domain):
                raise ValueError("the domain of the frame 'from_frame' is " +
                                 "not included in the given domain")
        if symbol is None:
            if from_frame is None:
                raise TypeError("some frame symbol must be provided")
            symbol = from_frame._symbol
            latex_symbol = from_frame._latex_symbol
            indices = from_frame._indices
            latex_indices = from_frame._latex_indices
            symbol_dual = from_frame._symbol_dual
            latex_symbol_dual = from_frame._latex_symbol_dual
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

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: e._repr_()
            'Vector frame (M, (e_0,e_1))'
            sage: repr(e)  # indirect doctest
            'Vector frame (M, (e_0,e_1))'
            sage: e  # indirect doctest
            Vector frame (M, (e_0,e_1))

        Test with a nontrivial destination map::

            sage: N = Manifold(3, 'N', start_index=1)
            sage: phi = M.diff_map(N)
            sage: h = M.vector_frame('h', dest_map=phi)
            sage: h._repr_()
            'Vector frame (M, (h_1,h_2,h_3)) with values on the 3-dimensional differentiable manifold N'

        """
        desc = "Local frame " + self._name + " on {}".format(self._vbundle)
        return desc

    def _new_instance(self, symbol, latex_symbol=None, indices=None,
                      latex_indices=None, symbol_dual=None,
                      latex_symbol_dual=None):
        pass

    ###### End of methods to be redefined by derived classes ######

    def coframe(self):
        r"""

        """
        return self._coframe

    def new_frame(self, change_of_frame, symbol, latex_symbol=None,
                  indices=None, latex_indices=None, symbol_dual=None,
                  latex_symbol_dual=None):
        r"""

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
            res._from_frame = self._from_frame
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

    def domain(self):
        r"""

        """
        return self._domain

    def base_space(self):
        r"""

        """
        return self._base_space

    def vector_bundle(self):
        r"""

        """
        return self._vbundle


#******************************************************************************

class TrivializationCoFrame(LocalCoFrame):
    r"""

    """
    def __init__(self, trivialization_frame, symbol, latex_symbol=None,
                 indices=None, latex_indices=None):
        r"""
        Construct a local coframe from a local trivialization.

        TESTS::


        """
        if not isinstance(trivialization_frame, TrivializationFrame):
            raise TypeError("the first argument must be a local trivialization "
                            "frame")
        LocalCoFrame.__init__(self, trivialization_frame, symbol,
                              latex_symbol=latex_symbol, indices=indices,
                              latex_indices=latex_indices)
        self._trivialization = trivialization_frame._trivialization

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::



        """
        return "Trivialization coframe " + self._name

#******************************************************************************

class TrivializationFrame(LocalFrame):
    r"""

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = TrivializationCoFrame

    def __init__(self, trivialization):
        r"""

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
        along = triv._name
        if along is None:
            along = triv.domain()._name
        along_latex = triv._latex_name
        if along_latex is None:
            along_latex = triv.domain()._latex_name
        symbol = tuple("e_" + str(i) + "({},{})".format(along, vbundle._name)
                       for i in range(1, rank + 1))
        latex_symbol = tuple(r'e_{' + latex(i) + r'}(' + along_latex + r','
                             + vbundle._latex_name + r')'
                             for i in range(1, rank + 1))
        symbol_dual = tuple("e^" + str(i) + "({},{})".format(along, vbundle._name)
                       for i in range(1, rank + 1))
        latex_symbol_dual = tuple(r'e^{' + latex(i) + r'}(' + along_latex + r','
                                  + vbundle._latex_name + r')'
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



        """
        return "Trivialization frame " + self._name

    def trivialization(self):
        r"""
        Return the underlying trivialization of ``self``.

        OUTPUT:

        -

        EXAMPLES::



        """
        return self._trivialization
r"""
Vector Bundle Fibers

The class :class:`VectorBundleFiber` implements fibers over a vector bundle.

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

from sage.symbolic.ring import SR
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.manifolds.vector_bundle_fiber_element import VectorBundleFiberElement

class VectorBundleFiber(FiniteRankFreeModule):
    r"""
    Fiber of a given vector bundle at a given point.

    Let `\pi: E \to M` be a topological vector bundle of rank `n` over the field
    `K` (see :class:`~sage.manifolds.vector_bundle.TopologicalVectorBundle`) and
    `p \in M`. The fiber `E_p` at `p` is defined via `E_p := \pi^{-1}(p)` and
    takes the structure of an `n`-dimensional vector space over the field `K`.

    INPUT:

    - ``vector_bundle`` -- :class:`~sage.manifolds.vector_bundle.VectorBundle`;
      vector bundle `E` on which the fiber is defined
    - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
      point `p` at which the fiber is defined

    EXAMPLES::

        sage: M = Manifold(4, 'M', structure='top')
        sage: X.<x,y,z,t> = M.chart()
        sage: p = M((0,0,0,0), name='p')
        sage: E = M.vector_bundle(2, 'E')
        sage: E.fiber(p)
        Fiber of E at Point p on the 4-dimensional topological manifold M

    """
    Element = VectorBundleFiberElement

    def __init__(self, vector_bundle, point):
        r"""

        """
        if point._manifold is not vector_bundle._base_space:
            raise ValueError("Point must be an element "
                             "of {}".format(vector_bundle._manifold))
        name = "{}_{}".format(vector_bundle._name, point._name)
        latex_name = r'{}_{{{}}}'.format(vector_bundle._latex_name,
                                         point._latex_name)
        self._rank = vector_bundle._rank
        self._vbundle = vector_bundle
        self._point = point
        self._base_space = point._manifold
        FiniteRankFreeModule.__init__(self, SR, self._rank, name=name,
                                      latex_name=latex_name,
                                      start_index=self._base_space._sindex)
        ###
        # Construct basis
        self._frame_bases = {} # dictionary of bases of the vector bundle fiber
                        # derived from local frames around the point
                        # (keys: local frames)
        for frame in vector_bundle._frames:
            # the frame is used to construct a basis of the vector bundle fiber
            # only if it is a frame for the given point:
            if point in frame.domain():
                coframe = frame.coframe()
                basis = self.basis(frame._symbol,
                                   latex_symbol=frame._latex_symbol,
                                   indices=frame._indices,
                                   latex_indices=frame._latex_indices,
                                   symbol_dual=coframe._symbol,
                                   latex_symbol_dual=coframe._latex_symbol)
                self._frame_bases[frame] = basis
        def_frame = vector_bundle._def_frame
        if def_frame in self._frame_bases:
            self._def_basis = self._frame_bases[def_frame]
        # The basis induced by the default frame of the manifold subset
        # in which the point has been created is declared the default
        # basis of self:
        def_frame = point.parent()._def_frame
        if def_frame in self._frame_bases:
            self._def_basis = self._frame_bases[def_frame]
        # Initialization of the changes of bases from the existing changes of
        # frames around the point:
        for frame_pair, automorph in self._vbundle._frame_changes.items():
            if point in automorph.domain():
                frame1, frame2 = frame_pair[0], frame_pair[1]
                fr1, fr2 = None, None
                for frame in self._frame_bases:
                    if frame1 in frame._subframes:
                        fr1 = frame
                        break
                for frame in self._frame_bases:
                    if frame2 in frame._subframes:
                        fr2 = frame
                        break
                if fr1 is not None and fr2 is not None:
                    basis1 = self._frame_bases[fr1]
                    basis2 = self._frame_bases[fr2]
                    auto = self.automorphism()
                    for frame, comp in automorph._components.items():
                        try:
                            basis = None
                            if frame is frame1:
                                basis = basis1
                            if frame is frame2:
                                basis = basis2
                            if basis is not None:
                                cauto = auto.add_comp(basis)
                                for ind, val in comp._comp.items():
                                    cauto._comp[ind] = val(point)
                        except ValueError:
                            pass
                    self._basis_changes[(basis1, basis2)] = auto

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,0,0), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.fiber(p)._repr_()
            'Fiber of E at Point p on the 3-dimensional topological manifold M'

        """
        return "Fiber of {} at {}".format(self._vbundle._name,
                                          self._point)

    def dimension(self):
        r"""
        Return the vector space dimension of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,0,0), name='p')
            sage: E = M.vector_bundle(2, 'E')
            sage: Ep = E.fiber(p)
            sage: Ep.dim()
            2

        """
        # The dimension is the rank of self as a free module:
        return self._rank

    dim = dimension

    def _an_element_(self):
        r"""
        Construct some (unnamed) vector in ``self``.

        EXAMPLES::



        """
        resu = self.element_class(self)
        if self._def_basis is not None:
            resu.set_comp()[:] = range(1, self._rank + 1)
        return resu

    def base_point(self):
        r"""
        Return the manifold point over which ``self`` is defined.

        EXAMPLES::



        """
        return self._point

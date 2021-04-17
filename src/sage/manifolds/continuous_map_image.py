r"""
Images of Continuous Maps as Manifold Subsets

:class:`ImageManifoldSubset` implements the image of a continuous map `\Phi`
from a manifold `M` to some manifold `N` as a subset `\Phi(M)` of `N`.

"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.manifolds.subset import ManifoldSubset

class ImageManifoldSubset(ManifoldSubset):
    r"""
    Subset of a topological manifold that is the image of a continuous map.

    INPUT:

    - ``map`` -- continuous map `\Phi`
    - ``inverse`` -- (default: ``None``) continuous map from
      ``map.codomain()`` to ``map.domain()``, which once restricted to the image
      of `\Phi` is the inverse of `\Phi` onto its image if the latter
      exists (NB: no check of this is performed)

    """

    def __init__(self, map, inverse=None, name=None, latex_name=None):
        r"""
        Construct a manifold subset that is the image of a continuous map.

        TESTS::

        """
        self._map = map
        self._inverse = inverse
        base_manifold = map.codomain()
        map_name = map._name or 'f'
        map_latex_name = map._latex_name or map_name
        if latex_name is None:
            if name is None:
                latex_name = map_latex_name + r'(' + map.domain()._latex_name + ')'
            else:
                latex_name = name
        if name is None:
            name = map_name + '_' + map.domain()._name
        ManifoldSubset.__init__(self, base_manifold, name, latex_name=latex_name)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

        """
        return "Image of the {}".format(self._map)

    def _an_element_(self):
        r"""
        Construct some point in the subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure="topological")
            sage: N = Manifold(1, 'N', ambient=M, structure="topological")
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart()
            sage: CN.add_restrictions([u > -1, u < 1])
            sage: Phi = N.continuous_map(M, {(CN,CM): [u, 1 + u^2]}, name='Phi')
            sage: Phi_N = Phi.image()
            sage: p = Phi_N.an_element(); p   # indirect doctest
            Point on the 2-dimensional topological manifold M
            sage: p.coordinates()
            (0, 1)
        """
        return self._map(self._map.domain().an_element())

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure="topological")
            sage: N = Manifold(1, 'N', ambient=M, structure="topological")
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart()
            sage: CN.add_restrictions([u > -1, u < 1])
            sage: Phi = N.continuous_map(M, {(CN,CM): [u, 1 + u^2]}, name='Phi')
            sage: Phi_inv = M.continuous_map(N, {(CM, CN): [x]}, name='Phi_inv')
            sage: Phi_N = Phi.image(inverse=Phi_inv)
            sage: M((0, 0)) in Phi_N
            False
            sage: M((0, 1)) in Phi_N
            True

        """
        if super().__contains__(point):
            return True
        if point not in self._map.codomain():
            return False
        if self._inverse is not None:
            preimage = self._inverse(point)
            if preimage not in self._map.domain():
                False
            return self._map(preimage) == point
        raise NotImplementedError

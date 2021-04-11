r"""
Images of Continuous Maps as Manifold Subsets
"""

from sage.manifolds.subset import ManifoldSubset

class ImageManifoldSubset(ManifoldSubset):

    def __init__(self, map, name=None, latex_name=None):

        self._map = map
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

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        TESTS::
        """
        if super().__contains__(point):
            return True
        if point not in self._map.codomain():
            return False
        raise NotImplementedError

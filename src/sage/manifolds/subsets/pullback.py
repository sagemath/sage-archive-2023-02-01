r"""
Manifold Subsets Defined as Pullbacks of Subsets under Continuous Maps
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

class ManifoldSubsetPullback(ManifoldSubset):

    """
    Manifold subset defined as a pullback of a subset under a continuous map.

    EXAMPLES::

        sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
        sage: M = Manifold(2, 'R^2', structure='topological')
        sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
        sage: r_squared = M.scalar_field(x^2+y^2)
        sage: r_squared.set_immutable()
        sage: I = RealSet((1, 4)); I
        (1, 4)
        sage: O = ManifoldSubsetPullback(r_squared, None, I); O
        Subset f_inv_(1, 4) of the 2-dimensional topological manifold R^2
        sage: M.point((1, 0)) in O
        False
        sage: M.point((1, 1)) in O
        True

    Using the embedding map of a submanifold::

        sage: M = Manifold(3, 'M', structure="topological")
        sage: N = Manifold(2, 'N', ambient=M, structure="topological")
        sage: N
        2-dimensional topological submanifold N immersed in the 3-dimensional topological manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()
        sage: t = var('t')
        sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
        sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
        sage: N.set_immersion(phi, inverse=phi_inv, var=t,
        ....:                 t_inverse={t: phi_inv_t})
        sage: N.declare_embedding()

        sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
        sage: S = M.open_subset('S', coord_def={CM: z<1})
        sage: D = ManifoldSubsetPullback(phi, codomain_subset=S)
        sage: N.point((2,0)) in D   ### BUG - the foliation parameters are in the way!
        True

    """

    def __init__(self, map, inverse=None, codomain_subset=None, name=None, latex_name=None):

        self._map = map
        self._inverse = inverse
        if codomain_subset is None:
            codomain_subset = map.codomain()
        self._codomain_subset = codomain_subset
        base_manifold = map.domain()
        map_name = map._name or 'f'
        map_latex_name = map._latex_name or map_name
        try:
            codomain_subset_latex_name = codomain_subset._latex_name
            codomain_subset_name = codomain_subset._name
        except AttributeError:
            from sage.misc.latex import latex
            codomain_subset_latex_name = str(latex(codomain_subset))
            codomain_subset_name = repr(codomain_subset)
        if latex_name is None:
            if name is None:
                latex_name = map_latex_name + r'^{-1}(' + codomain_subset_latex_name + ')'
            else:
                latex_name = name
        if name is None:
            name = map_name + '_inv_' + codomain_subset_name
        ManifoldSubset.__init__(self, base_manifold, name, latex_name=latex_name)

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        """
        if super().__contains__(point):
            return True
        return self._map(point) in self._codomain_subset

    def is_open(self):
        return self._codomain_subset.is_open()

    def is_closed(self):
        return self._codomain_subset.is_closed()

    def closure(self, name=None, latex_name=None):
        """
        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: r_squared = M.scalar_field(x^2+y^2)
            sage: r_squared.set_immutable()
            sage: I = RealSet((1, 2)); I
            (1, 2)
            sage: O = ManifoldSubsetPullback(r_squared, None, I); O
            Subset f_inv_(1, 2) of the 2-dimensional topological manifold R^2
            sage: latex(O)
            f^{-1}((1, 2))
            sage: O.closure()
            Subset f_inv_[1, 2] of the 2-dimensional topological manifold R^2
        """
        if self.is_closed():
            return self
        try:
            codomain_subset_closure = self._codomain_subset.closure()
        except AttributeError:
            return super().closure()
        closure = ManifoldSubsetPullback(self._map, self._inverse,
                                         codomain_subset_closure,
                                         name=name, latex_name=latex_name)
        closure.declare_superset(self)
        return closure

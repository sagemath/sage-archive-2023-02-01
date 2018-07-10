# *****************************************************************************
#  Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.manifold import Manifold

class ManifoldGenerator():
    @staticmethod
    def minkowski(names, sgn = 1):
        M = Manifold(4, 'M', structure = 'Lorentzian')
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        g = M.metric('g')
        g[0,0] = -sgn
        g[1,1], g[2,2], g[3,3] = sgn, sgn, sgn
        return M

    @staticmethod
    def sphere(dim=2, radius = 1, names=('th', 'phi'), embedded=False):
        M = Manifold(dim, 'S', structure='Riemannian')
        raise NotImplementedError("compliqué pour l'instant")

    @staticmethod
    def kerr(a, m, names, coordinates = "standard"):
        raise NotImplementedError()

    @staticmethod
    def euclidian(names):
        n = len(names)
        M = Manifold(n, 'M', structure="Riemannian")
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        g = M.metric('g')
        for i in range(n):
            g[i,i] = 1
        return M
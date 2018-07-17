r"""
Manifold Generator

the class :class:`ManifoldGenerator` implements some shortcuts to rapidly
create various simple manifolds.

AUTHORS:

- Florentin Jaffredo (2018) : initial version

"""


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
    def minkowski(names=None, sgn=1):
        """
        Generate a Minkowski space of dimension 4.

        By default the signature is set to `(- + + +)`, but can be changed to
        `(+ - - -)` by setting the optionnal argument to -1. The shortcut
        operator ``.<,>`` can be used to specify the coordinates.

        INPUT:

        - ``names`` -- (default: ``None``) name of the coordinates,
          automatically set by the shortcut operator.
        - ``sng`` -- (default: ``1``) set to `-1` to switch sign convention.

        OUTPUT:

        - Riemannian manifold of dimension 4 with Minkowskian metric.

        EXAMPLES::

            sage: M.<t, x, y, z> = ManifoldGenerator.minkowski()

        """
        M = Manifold(4, 'M', structure = 'Lorentzian')
        if names is None:
            names = ("t", "x", "y", "z")
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        g = M.metric('g')
        g[0,0] = -sgn
        g[1,1], g[2,2], g[3,3] = sgn, sgn, sgn
        return M

    @staticmethod
    def sphere(dim=2, radius=1, names=None):
        """
        Generate a sphere embedded in Euclidean space.

        The shortcut operator ``.<,>`` can be used to specify the coordinates.

        INPUT:

        - ``dim`` -- (default: ``2``) dimension of the sphere
        - ``radius`` -- (default: ``1``) radius of the sphere
        - ``names`` -- (default: ``None``) name of the coordinates,
          automatically set by the shortcut operator.

        OUTPUT:

        - Riemannian manifold.

        EXAMPLES::

            sage: S.<th, ph> = ManifoldGenerator.sphere()
            sage: S.plot({},srange(-pi,pi*21/20,pi/20),
            ....:   srange(-pi/2,pi/2*21/20,pi/20), viewer = 'threejs')
            Graphics3d Object

        """
        from functools import reduce
        from sage.functions.trig import cos, sin
        import operator

        def prod(iterable):
            return reduce(operator.mul, iterable, 1)

        if names is not None:
            dim = len(names)

        xnames = tuple(["x_{}".format(i) for i in range(dim+1)])
        E = ManifoldGenerator.euclidean(xnames)
        M = Manifold(dim, 'S', ambient=E, structure='Riemannian')
        if names is None:
            names = tuple(
                ["phi_{}:(0,pi)".format(i) for i in range(dim-1)] +
                ["phi_{}:(-pi,pi)".format(dim-1)])
        else:
            names = tuple([names[i]+":(0,pi)"for i in range(dim - 1)]+
                          [names[dim-1]+":(-pi,pi)"])
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        phi = M._first_ngens(dim)[:]
        coordfunc = [radius * cos(phi[i]) * prod(sin(phi[j]) for j in range(i))
                     for i in range(dim)] + [
                        radius * prod(sin(phi[j]) for j in range(dim))]
        imm = M.diff_map(E, coordfunc)
        M.set_embedding(imm)
        M.induced_metric()
        return M

    @staticmethod
    def kerr(m=1, a=0, names=None, coordinates = "BL"):
        """
        Generate a Kerr spacetime.

        A Kerr spacetime is a 4 dimensional manifold describing a rotating black
        hole. Two coordinates system are implemented : Boyer-Lindquist and ADM.
        The first only has a single non-diagonal term in he metric, making its
        Ricci tensor faster to compute, but is divergent on the event horizon.
        The second is only divergent at the center.

        The shortcut operator ``.<,>`` can be used to specify the coordinates.

        INPUT:

        - ``m`` -- (default: ``1``) Mass of the black hole in natural units
          (`c=1`, `G=1`)
        - ``a`` -- (default: ``0``) angular momentum in natural units. If set to
          0, the resulting spacetime corresponds to a Schwarzschild black hole.
        - ``names`` -- (default: ``None``) name of the coordinates,
          automatically set by the shortcut operator.
        - ``coordinates`` -- (default: ``"BL"``) either ``"BL"`` or ``"ADM"``.

        OUTPUT:

        - Lorentzian manifold.

        EXAMPLES::

            sage: K = ManifoldGenerator.kerr(2,1)

        """
        from sage.functions.all import sqrt, cos, sin
        M = Manifold(4, 'M', structure="Lorentzian")
        if coordinates == "ADM":
            if names is None:
                names = r't:(-oo,+oo) r:(0,+oo) th:(0,pi):\theta ph:(-pi,pi):\phi'
            C = M.chart(names)
            M._first_ngens = C._first_ngens
            g = M.metric('g')
            t, r, th, ph = C[:]
            rho = sqrt(r**2+a**2*cos(th)**2)
            g[0, 0], g[1, 1], g[2, 2], g[3, 3] = -(1-2*m*r/rho**2), 1+2*m*r/rho**2,\
                    rho**2, (r**2+a**2+2*a**2*m*r*sin(th)**2/rho**2)*sin(th)**2
            g[0, 1] = 2*m*r/rho**2
            g[0, 3] = -2*a*m*r/rho**2*sin(th)**2
            g[1, 3] = -a*sin(th)**2*(1+2*m*r/rho**2)
            return M
        if coordinates == "BL":
            if names is None:
                names = r't:(-oo,+oo) r:(0,+oo) th:(0,pi):\theta ph:(-pi,pi):\phi'
            C = M.chart(names)
            M._first_ngens = C._first_ngens
            g = M.metric('g')
            t, r, th, ph = C[:]
            rho = sqrt(r**2+a**2*cos(th)**2)
            g[0, 0], g[1, 1], g[2, 2], g[3, 3] = -(1-2*m*r/rho**2), \
                rho**2/(r**2-2*m*r+a**2), rho**2, \
                (r**2+a**2+2*m*r/rho**2*sin(th)**2)*sin(th)**2
            g[0, 3] = 2*m*r*a*sin(th)**2/rho**2
            return M
        raise NotImplementedError("Coordinates system not implemented, see help"
                                  " for details")

    @staticmethod
    def euclidean(names):
        """
        Generate a Euclidean space.

        The shortcut operator ``.<,>`` can be used to specify the coordinates.

        The function
        :func:`~sage.manifolds.differentiable.euclidean.EuclideanSpace` should
        be preferred in most cases.

        INPUT:

        - ``names`` -- (default: ``None``) name of the coordinates,
          automatically set by the shortcut operator.

        OUTPUT:

        - Riemannian manifold.

        EXAMPLES:

            sage: E.<x, y, z> = ManifoldGenerator.euclidean()

        """
        n = len(names)
        M = Manifold(n, 'M', structure="Riemannian")
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        g = M.metric('g')
        for i in range(n):
            g[i,i] = 1
        return M

    @staticmethod
    def torus(R=2, r=1, names=None):
        """
        Generate a 2-dimensional torus embedded in Euclidean space.

        The shortcut operator ``.<,>`` can be used to specify the coordinates.

        INPUT:

        - ``R`` -- (default: ``2``) Distance form the center to the center of
          the tube
        - ``r`` -- (default: ``1``) radius of the tube
        - ``names`` -- (default: ``None``) name of the coordinates,
          automatically set by the shortcut operator.

        OUTPUT:

        - Riemannian manifold.

        EXAMPLES::

            sage: T.<theta, phi> = ManifoldGenerator.torus(3,1)
            sage: T.plot({},srange(-pi,pi*21/20,pi/20),
            ....:   srange(-pi,pi*21/20,pi/20), viewer = 'threejs')
            Graphics3d Object

        """
        from sage.functions.all import cos, sin
        xnames = tuple(["x_{}".format(i) for i in range(3)])
        E = ManifoldGenerator.euclidean(xnames)
        M = Manifold(2, 'M', ambient=E, structure="Riemannian")
        if names is None:
            names = ("th", "ph")
        names = tuple([names[i] + ":(-pi,pi)" for i in range(2)])
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        th, ph = C[:]
        coordfunc = [(R+r*cos(th))*cos(ph), (R+r*cos(th))*sin(ph), r*sin(th)]
        imm = M.diff_map(E, coordfunc)
        M.set_embedding(imm)
        M.induced_metric()
        return M


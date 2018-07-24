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

def Minkowski(names=None, sgn=1):
    """
    Generate a Minkowski space of dimension 4.

    By default the signature is set to `(- + + +)`, but can be changed to
    `(+ - - -)` by setting the optionnal argument to -1, a full signature
    can also be provided. The shortcut operator ``.<,>`` can be used to
    specify the coordinates.

    INPUT:

    - ``names`` -- (default: ``None``) name of the coordinates,
      automatically set by the shortcut operator.
    - ``sgn`` -- (default: ``1``) set to `-1` to switch sign convention. A
      length 4 tuple is also accepted.

    OUTPUT:

    - pseudo-Riemannian manifold of dimension 4 with Minkowskian metric.

    EXAMPLES::

        sage: M.<t, x, y, z> = manifolds.Minkowski(sgn=(-1, 1, 1, 1))
        sage: M.metric()[:]
        [-1  0  0  0]
        [ 0  1  0  0]
        [ 0  0  1  0]
        [ 0  0  0  1]

    """
    from sage.manifolds.manifold import Manifold
    M = Manifold(4, 'M', structure = 'pseudo-Riemannian')
    if names is None:
        names = ("t", "x", "y", "z")
    C = M.chart(names=names)
    M._first_ngens = C._first_ngens
    g = M.metric('g')

    if hasattr(sgn, '__len__'):
        if len(sgn)!=4:
            raise ValueError("'sgn' must be a tuple of length 4")
        for i in range(4):
            g[i, i] = sgn[i]
    else:
        g[0,0] = -sgn
        g[1,1], g[2,2], g[3,3] = sgn, sgn, sgn
    return M

def Sphere(dim=2, radius=1, names=None, stereo2d=False, stereo_lim=None):
    """
    Generate a sphere embedded in Euclidean space.

    The shortcut operator ``.<,>`` can be used to specify the coordinates.

    INPUT:

    - ``dim`` -- (default: ``2``) dimension of the sphere
    - ``radius`` -- (default: ``1``) radius of the sphere
    - ``names`` -- (default: ``None``) name of the coordinates,
      automatically set by the shortcut operator.
    - ``stereo2d`` -- (default: ``False``) if True, defines only the
      stereographic charts, only implemented in 2d.
    - ``stereo_lim`` -- (default: ``None``) Parameter used to restrict the span
      of the stereographic charts, so that they don't cover the whole sphere.


    OUTPUT:

    - Riemannian manifold.

    EXAMPLES::

        sage: S.<th, ph> = manifolds.Sphere()
        sage: S.plot({},srange(-pi,pi*21/20,pi/20),
        ....:   srange(-pi/2,pi/2*21/20,pi/20), viewer = 'threejs').show()

    """
    from sage.manifolds.manifold import Manifold
    from functools import reduce
    from sage.functions.trig import cos, sin, atan, atan2
    from sage.functions.other import sqrt
    from sage.symbolic.constants import pi
    import operator

    def prod(iterable):
        return reduce(operator.mul, iterable, 1)

    if stereo2d:
        dim = 2
        xnames = ("X","Y","Z")
        E = Euclidean(names=xnames)
        S2 = Manifold(dim, 'S', ambient=E, structure='Riemannian')
        U = S2.open_subset('U')
        V = S2.open_subset('V')
        stereoN = U.chart(names=("x", "y"))
        x, y = stereoN[:]
        stereoS = V.chart(names=("xp", "yp"))
        xp, yp = stereoS[:]
        if stereo_lim is not None:
            stereoN.add_restrictions(x**2+y**2<stereo_lim**2)
            stereoS.add_restrictions(xp**2+yp**2<stereo_lim**2)

        stereoN_to_S = stereoN.transition_map(stereoS,
          (x / (x**2 + y**2), y / (x**2 + y**2)), intersection_name='W',
          restrictions1=x**2 + y**2 != 0, restrictions2=xp**2+yp**2!=0)
        stereoS_to_N = stereoS.transition_map(stereoN,
          (xp / (xp**2 + yp**2), yp / (xp**2 + yp**2)), intersection_name='W',
          restrictions1=x**2 + y**2 != 0, restrictions2=xp**2+yp**2!=0)
        stereoN_to_S.set_inverse(xp / (xp**2 + yp**2), yp / (xp**2 + yp**2))
        stereoS_to_N.set_inverse(x / (x**2 + y**2), y / (x**2 + y**2))
        W = U.intersection(V)
        stereoN_W = stereoN.restrict(W)
        stereoS_W = stereoS.restrict(W)
        A = W.open_subset('A', coord_def={stereoN_W: (y != 0, x < 0),
                                          stereoS_W: (yp != 0, xp < 0)})
        stereoN_A = stereoN_W.restrict(A)
        spher = A.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        th, ph = spher[:]
        spher_to_stereoN = spher.transition_map(stereoN_A, (sin(th)*cos(ph)
                /(1-cos(th)), sin(th)*sin(ph)/(1-cos(th))))
        spher_to_stereoN.set_inverse(2*atan(1/sqrt(x**2+y**2)), atan2(-y, -x)+pi)
        stereoN_to_S_A = stereoN_to_S.restrict(A)
        spher_to_stereoS = stereoN_to_S_A * spher_to_stereoN
        stereoS_to_N_A = stereoN_to_S.inverse().restrict(A)
        stereoS_to_spher = spher_to_stereoN.inverse() * stereoS_to_N_A

        coordfunc = [cos(th)*cos(ph), cos(th)*sin(ph), cos(ph)]
        imm = S2.diff_map(E, coordfunc)
        S2.set_embedding(imm)
        S2.induced_metric()

        return S2


    if names is not None:
        dim = len(names)

    xnames = tuple(["x_{}".format(i) for i in range(dim+1)])
    E = Euclidean(names=xnames)
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
    coordfunc = [radius * prod(sin(phi[j]) for j in range(dim))] +\
                [radius * cos(phi[i]) * prod(sin(phi[j]) for j in range(i))
                 for i in range(dim)]
    imm = M.diff_map(E, coordfunc)
    M.set_embedding(imm)
    M.induced_metric()
    return M

def Kerr(m=1, a=0, names=None, coordinates="BL"):
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

        sage: m, a = var('m, a')
        sage: K = manifolds.Kerr(m, a)

    """
    from sage.manifolds.manifold import Manifold
    from sage.functions.other import sqrt
    from sage.functions.trig import cos, sin
    M = Manifold(4, 'M', structure="Lorentzian")
    if coordinates == "ADM":
        if names is None:
            names = (r't:(-oo,+oo)', r'r:(0,+oo)', r'th:(0,pi):\theta', r'ph:(-pi,pi):\phi')
        C = M.chart(names=names)
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
            names = (r't:(-oo,+oo)', r'r:(0,+oo)', r'th:(0,pi):\theta', r'ph:(-pi,pi):\phi')
        C = M.chart(names=names)
        M._first_ngens = C._first_ngens
        g = M.metric('g')
        t, r, th, ph = C[:]
        rho = sqrt(r**2+a**2*cos(th)**2)
        g[0, 0], g[1, 1], g[2, 2], g[3, 3] = -(1-2*m*r/rho**2), \
            rho**2/(r**2-2*m*r+a**2), rho**2, \
            (r**2+a**2+2*m*r*a**2/rho**2*sin(th)**2)*sin(th)**2
        g[0, 3] = 2*m*r*a*sin(th)**2/rho**2
        return M
    raise NotImplementedError("Coordinates system not implemented, see help"
                              " for details")

def Euclidean(*args, **kwargs):
    """
    Generate a Euclidean space.

    Alias of :func:`~sage.manifolds.differentiable.euclidean.EuclideanSpace`.

    OUTPUT:

    - Riemannian manifold.

    EXAMPLES::

        sage: E.<x, y, z> = manifolds.Euclidean()

    """
    from sage.manifolds.differentiable.euclidean import EuclideanSpace
    return EuclideanSpace(*args, **kwargs)

def Torus(R=2, r=1, names=None):
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

        sage: T.<theta, phi> = manifolds.Torus(3,1)
        sage: T.plot({},srange(-pi,pi*21/20,pi/20),
        ....:   srange(-pi,pi*21/20,pi/20), viewer = 'threejs').show()

    """
    from sage.manifolds.manifold import Manifold
    from sage.functions.all import cos, sin
    xnames = tuple(["x_{}".format(i) for i in range(3)])
    E = Euclidean(names=xnames)
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


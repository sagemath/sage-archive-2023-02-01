r"""

"""

from sage.manifolds.differentiable.pseudo_riemannian_submanifold import \
                                                     PseudoRiemannianSubmanifold
from sage.categories.metric_spaces import MetricSpaces
from sage.categories.manifolds import Manifolds
from sage.rings.real_mpfr import RR
from sage.manifolds.differentiable.euclidean import EuclideanSpace

class StandardSphere(PseudoRiemannianSubmanifold):
    r"""

    """

    def __init__(self, n, radius=1, name=None, latex_name=None,
                 coordinates='spherical', category=None,
                 init_coord_methods=None, unique_tag=None):
        r"""

        """
        if radius <= 0:
            raise ValueError('radius must be greater than zero')
        E = EuclideanSpace(n+1)
        if name is None:
            name = 'S^{}'.format(n)
            if radius > 1:
                name += '({})'.format(radius)
            if latex_name is None:
                latex_name = r'\mathbb{S}^{' + str(n) + r'}'
                if radius > 1:
                    latex_name += '({})'.format(radius)
        if category is None:
            category = Manifolds(RR).Smooth() & MetricSpaces().Complete()
        PseudoRiemannianSubmanifold.__init__(self, n, name, ambient=E,
                                             signature=n, latex_name=latex_name,
                                             start_index=1, category=category)
        self._radius = radius
        self._coordinates = {}  # established coordinates
        self._init_coordinates = {'spherical': self._init_spherical,
                                  'stereographic': self._init_stereographic}
                                 # predefined coordinates
        if init_coord_methods:
            self._init_coordinates.update(init_coord_methods)
        # up here, the actual initialization begins:
        self._init_chart_domains()
        self._init_embedding()
        self._init_coordinates[coordinates]()

    def _init_embedding(self):
        r"""

        """
        name = 'iota'
        latex_name = r'\iota'
        iota = self.diff_map(self._ambient, name=name, latex_name=latex_name)
        self.set_embedding(iota)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage:

        """
        s = "Standard {}-sphere {} of radius {}".format(self._dim, self._name,
                                                        self._radius)
        return s

    def _init_chart_domains(self):
        r"""

        """
        # without south pole:
        name = self._name + '-{SP}'
        latex_name = self._latex_name + r'-\{\mathrm{SP}\}'
        self._stereoN_dom = self.open_subset(name, latex_name)
        # without north pole:
        name = self._name + '-{NP}'
        latex_name = self._latex_name + r'-\{\mathrm{NP}\}'
        self._stereoS_dom = self.open_subset(name, latex_name)
        # intersection:
        int = self._stereoN_dom.intersection(self._stereoS_dom)
        int._name = self._name + '-{NP,SP}'
        int._latex_name = self._latex_name + r'-\{\mathrm{NP}, \mathrm{SP}\}'
        # without half circle:
        self._spher_dom = int.open_subset('A')
        # declare union:
        self.declare_union(self._stereoN_dom, self._stereoS_dom)

    def _init_spherical(self):
        r"""
        Construct the chart of spherical coordinates.

        """
        # speed-up via simplification method
        self.set_simplify_function(lambda expr: expr.simplify_trig())

        # get domain
        A = self._spher_dom

        # initialize coordinates
        n = self._dim
        names = tuple(["phi_1:(-pi,pi):periodic"] +
                      ["phi_{}:(0,pi)".format(i) for i in range(2,n+1)])
        spher = A.chart(names=names)
        coord = spher[:]

        # manage embedding
        from sage.misc.misc_c import prod
        from sage.functions.trig import cos, sin

        R = self._radius

        coordfunc = [R*prod(sin(coord[i]) for i in range(n))]
        for k in range(n):
            c = R*cos(coord[k])*prod(sin(coord[i]) for i in range(k+1,n))
            coordfunc.append(c)
        cart = self._ambient.cartesian_coordinates()
        self._immersion.add_expr(spher, cart, coordfunc)

        # finish process
        self._coordinates['spherical'] = spher

        # adapt other coordinates
        if 'stereographic' in self._coordinates:
            self._transition_spher_stereo()

    def coordinates(self, type):
        r"""

        """
        if type not in self._coordinates:
            self._init_coordinates[type]()
        return self._coordinates[type]

    def _init_stereographic(self):
        r"""
        Construct the chart of stereographic coordinates.

        """
        # speed-up via simplification method
        self.set_simplify_function(lambda expr: expr.simplify_rational())

        # TODO: More speed-up?

        # get domains
        U = self._stereoN_dom
        V = self._stereoS_dom

        # initialize coordinates
        symbols_N = ''
        symbols_S = ''
        for i in self.irange():
            symbols_N += "y{}".format(i) + r":y_{" + str(i) + r"} "
            symbols_S += "yp{}".format(i) + r":y'_{" + str(i) + r"} "
        symbols_N = symbols_N[:-1]
        symbols_S = symbols_S[:-1]
        stereoN = U.chart(coordinates=symbols_N)
        stereoS = V.chart(coordinates=symbols_S)
        coordN = stereoN[:]
        coordS = stereoS[:]

        # predefine variables
        r2_N = sum(y ** 2 for y in coordN)
        r2_S = sum(yp ** 2 for yp in coordS)
        R = self._radius
        R2 = R**2

        # define transition map
        coordN_to_S = tuple(R*y/r2_N for y in coordN)
        coordS_to_N = tuple(R*yp/r2_S for yp in coordS)
        stereoN_to_S = stereoN.transition_map(stereoS, coordN_to_S,
                                              restrictions1=r2_N != 0,
                                              restrictions2=r2_S != 0)
        stereoN_to_S.set_inverse(*coordS_to_N, check=False)

        # manage embedding
        coordfuncN = [2*y*R2 / (R2+r2_N) for y in coordN]
        coordfuncN += [(R*r2_N-R*R2)/(R2+r2_N)]
        coordfuncS = [2*yp*R2 / (R2+r2_S) for yp in coordS]
        coordfuncS += [(R*R2-R*r2_S)/(R2+r2_S)]
        cart = self._ambient.cartesian_coordinates()
        self._immersion.add_expr(stereoN, cart, coordfuncN)
        self._immersion.add_expr(stereoS, cart, coordfuncS)
        self.clear_cache()  # clear cache of extrinsic information

        # define orientation
        eN = stereoN.frame()
        eS = stereoS.frame()  # oriented w.r.t. embedding
        frame_comp = list(eN[:])
        frame_comp[0] = -frame_comp[0]  # reverse orientation
        f = U.vector_frame('f', frame_comp)
        self.set_orientation([eS, f])

        # finish process
        self._coordinates['stereographic'] = [stereoN, stereoS]

        # adapt other coordinates
        if 'spherical' in self._coordinates:
            self._transition_spher_stereo()

    def _transition_spher_stereo(self):
        r"""

        """
        self.set_simplify_function(lambda expr: expr.simplify())

        # configure preexisting charts
        W = self._stereoN_dom.intersection(self._stereoS_dom)
        A = self._spher_dom
        stereoN, stereoS = self._coordinates['stereographic'][:]
        coordN = stereoN[:]
        coordS = stereoS[:]
        rstN = (coordN[0] != 0,)
        rstS = (coordS[0] != 0,)
        if self._dim > 1:
            rstN += (coordN[1] > 0,)
            rstS += (coordS[1] > 0,)
        stereoN_A = stereoN.restrict(A, rstN)
        stereoS_A = stereoS.restrict(A, rstS)
        self._coord_changes[(stereoN.restrict(W),
                             stereoS.restrict(W))].restrict(A)
        self._coord_changes[(stereoS.restrict(W),
                             stereoN.restrict(W))].restrict(A)
        spher = self._coordinates['spherical']

        R = self._radius
        n = self._dim

        # transition: spher to stereoN
        imm = self.embedding()
        cart = self._ambient.cartesian_coordinates()
        x = imm.expr(spher, cart)
        coordfunc = [(R*x[i])/(R-x[-1]) for i in range(n)]
        spher_to_stereoN = spher.transition_map(stereoN_A, coordfunc)

        # transition: stereoN to spher
        from sage.functions.trig import acos, atan2
        from sage.functions.special import sqrt
        from sage.symbolic.constants import pi
        x = imm.expr(stereoN, cart)
        coordfunc = [atan2(x[0],x[1])]
        for k in range(2, n+1):
            c = acos(x[k]/sqrt(sum(x[i]**2 for i in range(k+1))))
            coordfunc.append(c)
        spher_to_stereoN.set_inverse(*coordfunc, check=False)

        # transition spher <-> stereoS
        stereoN_to_S_A = self.coord_change(stereoN_A, stereoS_A)
        stereoN_to_S_A * spher_to_stereoN  # generates spher_to_stereoS
        stereoS_to_N_A = self.coord_change(stereoS_A, stereoN_A)
        spher_to_stereoN.inverse() * stereoS_to_N_A  # generates stereoS_to_spher

    def dist(self, p, q):
        r"""
        Return the geodesic distance between point ``p`` and ``q`` on ``self``.

        """
        from sage.functions.trig import acos
        # get Euclidean points:
        x = self._immersion(p)
        y = self._immersion(q)
        cart = self._ambient.cartesian_coordinates()
        x_coord = x.coord(chart=cart)
        y_coord = y.coord(chart=cart)
        n = len(x_coord)
        r = self._radius
        inv_angle = sum(x_coord[i]*y_coord[i] for i in range(n)) / r**2
        return r * acos(inv_angle)

    def radius(self):
        r"""
        Return the radius of ``self``.

        """
        return self._radius

    def simplicial_complex(self):
        r"""

        """
        from sage.homology.examples import Sphere
        return Sphere(self._dim)

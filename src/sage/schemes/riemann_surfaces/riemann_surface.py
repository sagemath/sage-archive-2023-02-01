r"""
Computation of Riemann matrices and endomorphism rings of algebraic Riemann surfaces.

This module provides a class, RiemannSurface, to model the Riemann surface
determined by a plane algebraic curve over a subfield of the complex numbers.

A homology basis is derived from the edges of a voronoi cell decomposition based on
the branch locus. The pull-back of these edges to the Riemann surface provides
a graph on it that contains a homology basis.

The class provides methods for computing the Riemann period matrix of the
surface numerically, using a certified homotopy continuation method due to
[Stefan Kranich, An epsilon-delta bound for plane algebraic curves and its use
for certified homotopy continuation of systems of plane algebraic curves,
arXiv:1505.03432, 2016]

The class also provides facilities for computing the endomorphism ring of the
period lattice numerically, by determining integer (near) solutions to the relevant
approximate linear equations.

AUTHORS:

- Alexandre Zotine, Nils Bruin (2017-06-10): initial version

EXAMPLES:

We compute the Riemann matrix of a genus 3 curve::

    sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
    sage: R.<x,y>=QQ[]
    sage: f=x^4-x^3*y+2*x^3+2*x^2*y+2*x^2-2*x*y^2+4*x*y-y^3+3*y^2+2*y+1
    sage: S=RiemannSurface(f,prec=100)
    sage: M=S.riemann_matrix()

We test the usual properties, i.e., that the period matrix is symmetric and that
the imaginary part is positive definite::

    sage: all(abs(a) < 1e-20 for a in (M-M.T).list())
    True
    sage: iM=Matrix(RDF,3,3,[a.imag_part() for a in M.list()])
    sage: iM.is_positive_definite()
    True

We compute the endomorphism ring and check it has ZZ-rank 6::

    sage: A=S.endomorphism_basis(80,8)
    sage: len(A) == 6
    True

In fact it is an order in a number field. We check this.

    sage: T.<t>=QQ[]
    sage: K.<a>=NumberField(t^6 - t^5 + 2*t^4 + 8*t^3 - t^2 - 5*t + 7)
    sage: all(len(a.minpoly().roots(K)) == a.minpoly().degree() for a in A)
    True

"""

#*****************************************************************************
#       Copyright (C) 2017 Alexandre Zotine, Nils Bruin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from scipy.spatial import Voronoi, voronoi_plot_2d
from sage.misc.cachefunc import cached_method
from sage.rings.complex_field import ComplexField
from sage.rings.real_mpfr import RealField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.arith.srange import srange
from sage.ext.fast_callable import fast_callable
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.interfaces.singular import singular
from sage.graphs.graph import Graph
from sage.matrix.constructor import Matrix
from sage.modules.free_module import VectorSpace
from sage.numerical.gauss_legendre import integrate_vector
from sage.misc.misc_c import prod

def voronoi_ghost(cpoints):
    """
    Converts a set of complex points to a list of real tuples `(x,y)`,
    as well as 'ghost' points to surround the complex points.

    In particular, given some list of complex points to surround by Voronoi
    cells, creates a list of additional points so that the voronoi cells
    surrounding the points are bounded in the plane.

    INPUT:

    - ``cpoints`` -- a list of complex numbers.

    OUTPUT:

    Returns a list of real tuples `(x,y)` consisting of the original points and
    a set of points which surround them.

    EXAMPLES:

    Create a set of complex points, then runs the function::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import voronoi_ghost
        sage: L = [1 + 1*I, 1 - 1*I, -1 + 1*I, -1 - 1*I]
        sage: voronoi_ghost(L)
        [(1, 1), (1, -1), (-1, 1), (-1, -1), (-2, -2), (-2, 2), (2, -2), (2, 2), (4, 0), (-4, 0), (0, 4), (0, -4)]

    NOTES:

    - Mainly for internal use.

    - Box formation is somewhat naive. The box is not necessarily
      centered at the average of the points, but is just constructed so
      that it contains all of the points.

    TODO:

    It may be useful to add an option to compute a circle instead of
    a box.
    """
    points = [(z.real_part(), z.imag_part()) for z in cpoints]
    # If there's only a single point, return unit box.
    if len(points) == 1:
        return [(points[0][0] + 0.5, points[0][1] + 0.5), (points[0][0] + 0.5, points[0][1] - 0.5),
                (points[0][0] - 0.5, points[0][1] + 0.5), (points[0][0] - 0.5, points[0][1] - 0.5)]
    # Initialize maxes and mins. these form the initial box [x1,x2] x [y1, y2].
    maxx = 2*max(x for x,y in points)
    minx = 2*min(x for x,y in points)
    maxy = 2*max(y for x,y in points)
    miny = 2*min(y for x,y in points)
    # Because of the structure of the initial box, these conditionals are
    # added to ensure that none of the original points lie on the boundary
    # of the box.
    if abs(maxx) < 1e-8:
        maxx = -minx
    if abs(minx) < 1e-8:
        minx = -maxx
    if abs(maxy) < 1e-8:
        maxy = -miny
    if abs(miny) < 1e-8:
        miny = -maxy
    # These conditions are to ensure that the box is not trivial-
    # i.e. is not of the form [x1,x1] x [y1, y2], etc.
    if (abs(minx - maxx) < 10**(-8)):
        minx = points[0][0] - (abs(maxy) + abs(miny))/2
        maxx = points[0][0] + (abs(maxy) + abs(miny))/2
    if (abs(miny - maxy) < 10**(-8)):
        miny = points[0][1] - (abs(maxx) + abs(minx))/2
        maxy = points[0][1] + (abs(maxx) + abs(minx))/2
    return points + [(x,y) for x in [minx,maxx] for y in [miny,maxy]] \
                  + [(2*maxx,0),(2*minx,0),(0,2*maxy),(0,2*miny)]

def bisect(L,t):
    """
    Find position in a sorted list using bisection.

    Given a list `L = [(t_0,...),(t_1,...),...(t_n,...)]` with
    increasing t_i, find the index i such that `t_i <= t < t_{i+1}` using bisection.
    The rest of the tuple is available for whatever use required.

    INPUT:

    - ``L`` -- A list of tuples such that the first term of each tuple is a real
      number between 0 and 1. These real numbers must be increasing.

    - ``t`` -- A real number between `t_0` and `t_n`.

    OUTPUT:

    An integer i, giving the position in L where t would be in

    EXAMPLES:

    Form a list of the desired form, and pick a real number between 0 and 1::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import bisect
        sage: L = [(0.0, 'a'), (0.3, 'b'), (0.7, 'c'), (0.8, 'd'), (0.9, 'e'), (1.0, 'f')]
        sage: t = 0.5
        sage: bisect(L,t)
        1

    Another example which demonstrates that if t is equal to one of the t_i, it
    returns that index::

        sage: L = [(0.0, 'a'), (0.1, 'b'), (0.45, 'c'), (0.5, 'd'), (0.65, 'e'), (1.0, 'f')]
        sage: t = 0.5
        sage: bisect(L,t)
        3

    """
    # Defining starting indices for the loop.
    min = 0
    max = len(L) - 1
    # If the input t is not between 0 and 1, raise an error.
    if t < L[min][0] or t > L[max][0]:
        raise ValueError("Value for t out of range")
    # Main loop.
    while (min < max-1):
        # Bisect.
        mid = (max+min)//2
        # If it's equal, return the index we bisected to.
        if t == L[mid][0]:
            return mid
        # If it's smaller, then we're on the left side.
        elif t < L[mid][0]:
            max = mid
        # Otherwise we're on the right side.
        else:
            min = mid
    # Once the loop terminates, we return what the indices converged to.
    return min

class ConvergenceError(ValueError):
    r"""
    Error object suitable for raising and catching when Newton iteration fails.

    EXAMPLES::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import ConvergenceError
        sage: raise ConvergenceError("test")
        Traceback (most recent call last):
        ...
        ConvergenceError: test
        sage: isinstance(ConvergenceError(),ValueError)
        True
    """
    pass

class RiemannSurface(object):
    """
    Constructs a Riemann Surface. This is specified by the zeroes of a
    bivariate polynomial with rational coefficients `f(z,w) = 0`.

    INPUT:

    - ``f`` -- a bivariate polynomial with rational coefficients.
      The surface is interpreted as the covering space of the
      coordinate plane in the first variable.

    - ``prec`` -- the desired precision of computations on the surface
      in bits. Default 53.

    - ``certification`` -- a boolean value indicating whether homotopy
      continuation is certified or not. Certified results have a higher
      chance of being accurate at the cost of computation time.
      Default True.

    EXAMPLES:

    Constructing a Riemann surface.::

        sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
        sage: R.<z,w> = QQ[]
        sage: f = w^2 - z^3 + 1
        sage: RiemannSurface(f)
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 53 bits of precision

    Constructing another Riemann surface with 100 bits of precision.::

        sage: S = RiemannSurface(f, prec=100); S
        Riemann surface defined by polynomial f = -z^3 + w^2 + 1 = 0, with 100 bits of precision
        sage: S.riemann_matrix() #abs tol 0.00000001
        [0.500000000000000000000000... + 0.866025403784438646763723...*I]

    TODO: have an example where self.genus < 0 to raise error

    """
    def __init__(self, f, prec=53, certification=True):
        # Initializations.
        self._prec = prec
        self._certification = certification
        self._R = f.parent()
        if len(self._R.gens()) != 2:
            raise ValueError('Only bivariate polynomials supported.')
        z, w = self._R.gen(0), self._R.gen(1)
        self._CC = ComplexField(self._prec)
        self._RR = RealField(self._prec)
        self._CCz = PolynomialRing(self._CC, [self._R.gen(0)])
        self._CCw = PolynomialRing(self._CC, [self._R.gen(1)])
        self.f = f
        # Using singular to compute the genus quickly. Later, use rank/2 of
        # gram matrix?
        self.genus = self._R.ideal(self.f).genus()
        # TODO: can we get rid of the comments below?
        # I'm not sure if this check will work, but could be useful?
        # I think singular's genus computation returns a negative genus
        # if it's reducible. Either way, if the genus is negative we should
        # terminate.
        if self.genus < 0:
            raise AttributeError('Error: polynomial is geometricially reducible')
        self.degree = self.f.degree(w)
        self._dfdw = self.f.derivative(w)
        self._dfdz = self.f.derivative(z)
        self._discriminant = self.f.resultant(self._dfdw,w)
        # Coefficients of the polynomial for use in homotopy continuation.
        self._a0 = self._CCz(self.f.coefficient({w:self.degree})(self._CCz.gen(),0))
        self._a0roots = self._a0.roots(multiplicities=False)
        self._aks = [self._CCz(self.f.coefficient({w:self.degree - k - 1})
                               (self._CCz.gen(),0)) for k in range(self.degree)]
        # Compute the branch locus. Takes the square-free part of the discriminant
        # because of numerical issues.
        self.branch_locus = []
        for x in self._discriminant.factor():
            self.branch_locus += self._CCz(x[0](self._CCz.gen(),0)).roots(multiplicities=False)
        # Voronoi diagram and the important points associated with it.
        self.voronoi_diagram = Voronoi(voronoi_ghost(self.branch_locus))
        self._vertices = [self._CC(x0,y0) for x0,y0 in self.voronoi_diagram.vertices]
        self._wvalues = [self.w_values(z0) for z0 in self._vertices]
        self._Sn = SymmetricGroup(srange(self.degree))
        self._L = dict()
        self._PM = None
        self._RM = None
        self._differentials = None
        self._fastcall_f = fast_callable(f,domain=self._CC)
        self._fastcall_dfdw = fast_callable(self._dfdw,domain=self._CC)
        self._fastcall_dfdz = fast_callable(self._dfdz,domain=self._CC)

    def __repr__(self):
        """
        Return a string representation of the Riemann surface class.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: RiemannSurface(f)
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision

        """
        s = 'Riemann surface defined by polynomial f = %s = 0, with %s bits of precision'%(self.f, self._prec)
        return s

    def w_values(self, z0):
        """
        Returns the points lying on the surface above `z0`.

        INPUT:

        - ``z0`` -- (complex) a point in the complex z-plane.

        OUTPUT:

        A set of complex numbers corresponding to solutions of `f(z0,w) = 0`.

        EXAMPLES:

        Setting up a Riemann Surface.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Find the w-values above the origin, i.e. the solutions of `w^2 + 1 = 0`.::

            sage: S.w_values(0)
            [-1.00000000000000*I, 1.00000000000000*I]

        """
        return self.f(z0,self._CCw.gen(0)).roots(multiplicities=False)

    @cached_method
    def downstairs_edges(self):
        """
        Procedure to compute the edgeset of the Voronoi diagram. For use in
        forming a graph.

        OUTPUT:

        A list of integer tuples corresponding to edges between vertices
        in the Voronoi diagram.

        EXAMPLES:

        Form a Riemann surface, one with a particularly simple branch locus.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^3 - z^2
            sage: S = RiemannSurface(f)

        Compute the edges::

            sage: S.downstairs_edges()
            [(0, 1), (0, 3), (1, 2), (2, 6), (3, 4), (3, 6), (4, 5), (5, 6)]

        This now gives an edgeset which one could use to form a graph.

        .. NOTE::

            The numbering of the vertices is given by the Voronoi package.
        """
        # Because of how we constructed the Voronoi diagram, the first n points
        # correspond to the branch locus points.
        # The regions of these points are all of the edges which don't go off
        # to infinity, which are exactly the ones we want.
        n = len(self.branch_locus)
        desired_edges = [self.voronoi_diagram.regions[self.voronoi_diagram.point_region[i]] for i in range(n)]
        # First I construct the edges as a set- because the regions will overlap
        # and I don't want to have two of the same edge.
        edges1 = set()
        for c in desired_edges:
            for j in range(len(c)-1):
                edges1.add(frozenset((c[j],c[j+1])))
            edges1.add(frozenset((c[0],c[-1])))
        # Then I just make it into a list and sort it.
        # The sorting is important - it will make computing the monodromy group
        # MUCH easier.
        # We orient all the edges so that we go from lower to higher
        # numbered vertex for the continuation.

        #TODO: Why not make a list in the first place with the edges?
        edges = []
        for e in edges1:
            i0, i1 = e
            if i0 > i1:
                i1,i0 = i0,i1
            edges += [(i0,i1)]
        edges.sort()
        return edges

    def _compute_delta(self, z1, epsilon):
        """
        Computes a delta for homotopy continuation when moving along a path.

        INPUT:

        - ``z1`` -- a complex number in the z-plane.

        - ``epsilon`` -- a real number, which is the minimum distance between
          the w-values above z1.

        OUTPUT:

        A real number, which is a step size for moving along a path.

        EXAMPLES:

        Form a Riemann Surface.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)

        Pick a point which lies on the voronoi diagram, and compute an
        appropriate epsilon.::

            sage: z1 = S._vertices[0]
            sage: currw = S.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: S._compute_delta(z1, epsilon)
            0.40790892870...

        If the Riemann surface doesn't have certified homotopy continuation,
        then the delta will just be the minimum distance away from a branch
        point.::

            sage: T = RiemannSurface(f, certification=0)
            sage: z1 = T._vertices[0]
            sage: currw = T.w_values(z1)
            sage: n = len(currw)
            sage: epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            sage: T._compute_delta(z1, epsilon)
            0.79056941504...

        """
        # If we don't want certified computations, then doing this is
        # unneccessary.
        if (self._certification == False):
            # Instead, we just compute the minimum distance between branch
            # points and the point in question.
            return min([abs(z1 - self.branch_locus[i]) for i in range(len(self.branch_locus))])/2
        else:
            wvalues = self.w_values(z1)
            # For computation of rho. Need the branch locus + roots of a0.
            badpoints = self.branch_locus + self._a0roots
            # Compute rho.
            rho = min(abs(z1 - z) for z in badpoints)/2
            # Compute Y.
            Y = max(abs(self._fastcall_dfdz(z1,wi)/self._fastcall_dfdw(z1,wi)) for wi in wvalues)
            # Compute M.
            upperbounds = [sum(ak[k]*(abs(z1) + rho)**k for k in range(ak.degree())) for ak in self._aks]
            upperbounds.reverse()
            # If a0 is a constant polynomial, it is obviously bounded below.
            if self._a0roots == []:
                lowerbound = self._CC(self._a0)/2
            else:
                lowerbound = self._a0[self._a0.degree()]*prod(abs((zk - z1) - rho) for zk in self._a0roots)/2
            M = 2*max(abs((upperbounds[k]/lowerbound))**(1/(k+1)) for k in range(self.degree-1))
            return rho*( ((rho*Y - epsilon)**2 + 4*epsilon*M).sqrt() - (rho*Y + epsilon))/(2*M - 2*rho*Y)

    def homotopy_continuation(self, edge):
        """
        Performs homotopy continuation along an edge of the Voronoi diagram
        using Newton iteration.

        INPUT:

        - ``edge`` -- A tuple of integers indicating an edge of the Voronoi
          diagram.

        OUTPUT:

        Returns a list of complex numbers corresponding to the points which are
        reached when traversing along the direction of the edge. The ordering of
        these points indicates how they have been permuted due to the weaving of
        the curve.

        EXAMPLES:

        Here we consider two examples::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)

        First we traverse the first edge on the Voronoi diagram, which is
        `(0,1)`.::

            sage: edge1 = (0,1)
            sage: homotopyvalues = S.homotopy_continuation(edge1); homotopyvalues #abs tol 0.00000001
            [-5.83317387179... + 1.6687560771...*I, -0.08373298599... - 0.0320669890287...*I, 5.91690685778... - 1.63668908809...*I]

        Compare this with the w-values at the index 1 vertex::

            sage: originalwvalues = S.w_values(S._vertices[1])
            sage: abs(originalwvalues[0] - homotopyvalues[0]) #abs tol 0.00000001
            0.000000000000...
            sage: abs(originalwvalues[1] - homotopyvalues[1]) #abs tol 0.00000001
            0.000000000000...
            sage: abs(originalwvalues[2] - homotopyvalues[2]) #abs tol 0.00000001
            0.000000000000...

        This indicates that traversing the path of the edge along the surface
        does not permute the points. On the other hand, consider the following::

            sage: edge2 = (9,11)
            sage: homotopyvalues = S.homotopy_continuation(edge2); homotopyvalues #abs tol 0.00000001
            [-0.47094999230... - 0.8602729692...*I, -0.4725503696... + 0.59255942108...*I, 0.94350036195... + 0.267713548165...*I]

        In this situation, the first and second roots are interchanged, which
        indicates that the layers of the surface are interchanging when
        traversing along the direction of the edge `(9,11)`.::

            sage: originalwvalues = S.w_values(S._vertices[11])
            sage: abs(originalwvalues[0] - homotopyvalues[1]) #abs tol 0.00000001
            0.00000000000...
            sage: abs(originalwvalues[1] - homotopyvalues[0]) #abs tol 0.00000001
            0.00000000000...

        """
        i0, i1 = edge
        ZERO=self._RR.zero()
        ONE=self._RR.one()
        datastorage = []
        z_start = self._CC(self._vertices[i0])
        z_end = self._CC(self._vertices[i1])
        def path(t):
            return z_start*(1-t) + z_end*t
        # Primary procedure.
        T = ZERO
        currw = self.w_values(path(T))
        n = len(currw)
        epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
        datastorage += [(T,currw,epsilon)]
        while T < ONE:
            epsilon = min([abs(currw[i] - currw[n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            delta = self._compute_delta(path(T), epsilon)/abs(z_end - z_start)
            # Move along the path by delta.
            T += delta
            # If T exceeds 1, just set it to 1 and compute.
            if T > ONE:
                delta -= (T-ONE)
                T = ONE
            while True:
                try:
                    neww = self._determine_new_w(path(T),currw,epsilon)
                except ConvergenceError:
                    delta /= 2
                    T -= delta
                else:
                    break
            currw = neww
            datastorage += [(T,currw,epsilon)]
        self._L[edge] = datastorage
        return currw

    def _determine_new_w(self,z0,oldw,epsilon):
        """
        A procedure to Newton iterate a list of w-values simultaneously.
        Used primarily for moving along the surface for integration or
        homotopy continuation.

        INPUT:

        - ``z0`` -- A complex number.

        - ``oldw`` -- A list of w-values which are presumed to be guesses of
          the w-values above ``z0``.

        - ``epsilon`` -- ('The minimum distance between the points of oldw')/3

        OUTPUT:

        A list of points the same length as oldw corresponding to the
        new newton iterated points.

        However, if the newton iteration exceeded the alloted attempts, or
        exited the epsilon ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the roots are.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)
            sage: neww = S._determine_new_w(z0,oldw,epsilon); neww #abs tol 0.00000001
            [-5.9288247295... - 2.53001238598...*I, 5.9288247295... + 2.53001238598...*I]

        Which should be exactly the same as the w-values we started with.::

            sage: abs(neww[0] - oldw[0]) #abs tol 0.00000001
            0.000000000000...
            sage: abs(neww[1] - oldw[1]) #abs tol 0.00000001
            0.000000000000...

        Here is an example where we exit the epsilon bound. This approach
        is based on the homotopy continuation procedure which traverses along
        a path and attempts newton iteration.::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) - T._vertices[15]*(0.1)
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])
            sage: T._determine_new_w(z0,oldw,epsilon)
            Traceback (most recent call last):
            ...
            ConvergenceError: Newton iteration escaped neighbourhood

        """
        # Tools of newton iteration.
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        neww = []
        # Iterate over all roots.
        for i in range(len(oldw)):
            delta = F(z0,oldw[i])/dF(z0,oldw[i])
            Ndelta = delta.norm()
            wi = oldw[i]-delta
            for j in range(100):
                # If we exceed the epsilon bound from homotopy continuation,
                # terminate.
                if abs(wi - oldw[i]) >= epsilon:
                    raise ConvergenceError("Newton iteration escaped neighbourhood")
                new_delta = F(z0,wi)/dF(z0,wi)
                Nnew_delta = new_delta.norm()
                # If we found the root exactly, or our deltas stop getting
                # smaller, terminate and change flag.
                if (new_delta == 0) or (Nnew_delta>=Ndelta):
                    neww.append(wi)
                    break
                delta=new_delta
                Ndelta=Nnew_delta
                wi-=delta
            # If we run 100 iterations without a result, terminate.
            else:
                raise ConvergenceError("Newton interation fails to converge after %s iterations"%(j,))
        return neww

    def _newton_iteration(self,z0,oldw,epsilon):
        """
        A non-vectorized newton iteration procedure used for integration.

        INPUT:

        - ``z0`` -- A complex number.

        - ``oldw`` -- A w-value which is presumed to be a guess of one of
          the w-values above ``z0``.

        - ``epsilon`` -- the minimum distance between the w-values divided by 3.

        OUTPUT:

        A complex number, which should be a w-value above ``z0``.

        However, if the newton iteration exceeded the alloted attempts, or
        exited the epsilon ball, raises a convergence error.

        EXAMPLES:

        First, a trivial example where we guess exactly what the root is.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: z0 = S._vertices[0]
            sage: epsilon = 0.1
            sage: oldw = S.w_values(z0)[0]
            sage: neww = S._newton_iteration(z0,oldw,epsilon); neww #abs tol 0.00000001
            -5.9288247295... - 2.53001238598...*I

        Which should be exactly the same as the w-value we started with.::

            sage: oldw - neww #abs tol 0.00000001
            0.000000000000000

        Here is an example where we exit the epsilon bound. This approach
        is based on the homotopy continuation procedure which traverses along
        a path and attempts newton iteration.::

            sage: g = z^3*w + w^3 + z
            sage: T = RiemannSurface(g)
            sage: z0 = T._vertices[2]*(0.9) - T._vertices[15]*(0.1)
            sage: epsilon = 0.5
            sage: oldw = T.w_values(T._vertices[2])[0]
            sage: T._newton_iteration(z0, oldw, epsilon)
            Traceback (most recent call last):
            ...
            ConvergenceError: Newton iteration escaped neighbourhood

        """
        F = self._fastcall_f
        dF = self._fastcall_dfdw
        delta = F(z0,oldw)/dF(z0,oldw)
        Ndelta = delta.norm()
        neww = oldw-delta
        eps_squared = epsilon**2
        for j in range(100):
            if (neww-oldw).norm() > eps_squared:
                raise ConvergenceError("Newton iteration escaped neighbourhood")
            new_delta = F(z0,neww)/dF(z0,neww)
            Nnew_delta = new_delta.norm()
            if (new_delta == 0) or (Nnew_delta>=Ndelta):
                return neww
            delta = new_delta
            Ndelta = Nnew_delta
            neww-=delta
        raise ConvergenceError("Newton iteration fails to converge")

    @cached_method
    def upstairs_edges(self):
        """
        Computes the edgeset of the lift of the downstairs graph onto the
        Riemann surface.

        OUTPUT:

        An edgeset between vertices (i, j), where i corresponds to the i-th
        point in the Voronoi diagram vertices, and j is the j-th w-value
        associated with that point.

        EXAMPLES:

        Computing the upstairs edges.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 + z^3 - z^2
            sage: S = RiemannSurface(f)
            sage: edgeset = S.upstairs_edges(); edgeset
            [[(0, 0), (1, 1)], [(0, 1), (1, 0)], [(0, 0), (3, 0)], [(0, 1), (3, 1)], [(1, 0), (2, 1)], [(1, 1), (2, 0)], [(2, 0), (6, 0)], [(2, 1), (6, 1)], [(3, 0), (4, 0)], [(3, 1), (4, 1)], [(3, 0), (6, 0)], [(3, 1), (6, 1)], [(4, 0), (5, 1)], [(4, 1), (5, 0)], [(5, 0), (6, 0)], [(5, 1), (6, 1)]]

        """
        edgeset = []
        n = len(self._wvalues[0])
        # Lifts each edge individually.
        for e in self.downstairs_edges():
            i0, i1 = e
            # Epsilon for checking w-value later.
            epsilon = min([abs(self._wvalues[i1][i] - self._wvalues[i1][n-j-1]) for i in range(n) for j in range(n-i-1)])/3
            # Homotopy continuation along e.
            homotopycont = self.homotopy_continuation(e)
            for i in range(len(homotopycont)):
                # Checks over the w-values of the next point to check which it is.
                for j in range(len(self._wvalues[i1])):
                    if abs(homotopycont[i] - self._wvalues[i1][j]) < epsilon:
                        # Once it finds the appropriate w-value, adds the edge.
                        edgeset = edgeset + [[(i0,i),(i1,j)]]
                        continue
        return edgeset

    def _edge_permutation(self, edge):
        """
        Computes the permutation of the w-values above a point in the z-plane
        when moving along an edge in the Voronoi diagram.

        INPUT:

        - ``edge`` -- an edge on the Voronoi diagram. Raises an error
          if input isn't on the Voronoi diagram.

        OUTPUT:

        A permutation corresponding to how the roots interchange when
        moving along the edge.

        EXAMPLES:

        Compute the edge permutation of (5,16) on the Voronoi diagram.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: S._edge_permutation((5,16))
            (1,2)

        This indicates that while traversing along the direction of `(5,16)`,
        the 2nd and 3rd layers of the Riemann surface are interchanging.
        """
        if edge in self.downstairs_edges():
            # For every edge downstairs, I get self.degree = d many edges
            # upstairs. This means I count in multiples of d.
            # Find the position downstairs.
            n = self.downstairs_edges().index(edge)
            tempPerm = []
            for i in range(self.degree):
                # Count up to the n*d-th edge, and then take the resulting d
                # edges. Because of how I computed these edges to begin with,
                # the first vertex in the edge is always ordered 1..d. This
                # means that the second vertex in the edge completely encodes
                # the permutation.
                tempPerm += [self.upstairs_edges()[(self.degree)*n+i][1][1]]
            return self._Sn(tempPerm)
        else:
            raise ValueError('edge not in Voronoi diagram.')

    @cached_method
    def edge_permutations(self):
        """
        Computes all of the edge permutations of the Voronoi diagram and
        packages them together.

        OUTPUT:

        Returns a list indexed by the edges of the Voronoi diagram
        of the permutation associated with moving along the corresponding edge.

        EXAMPLES:

        Compute the set of edge permutations.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: 2 + 2
            4
            sage: f = w^2 + z^3 - z^2
            sage: S = RiemannSurface(f)
            sage: S.edge_permutations()
            [[(0, 1), (0,1)], [(0, 3), ()], [(1, 2), (0,1)], [(2, 6), ()], [(3, 4), ()], [(3, 6), ()], [(4, 5), (0,1)], [(5, 6), ()]]

        We see that if we start at the branch labeled 0 above the vertex 0 and continue along
        the edge (0, 1), we arrive at the branch labeled 1 above vertex 1.
        On the other hand, if we do the same thing along the edge (0,3), we end up at the branch labeled 0
        above 3.
        """
        datastorage = []
        for e in self.downstairs_edges():
            tempPerm = self._edge_permutation(e)
            datastorage += [[e, tempPerm]]
        return datastorage

    @cached_method
    def monodromy_group(self):
        """
        Computes the monodromy group of a Riemannian surface.

        OUTPUT:

        A list of the local monodromy associated with each branch point,
        encoded as a permutation. The path is positively oriented, but we presently
        do not control the base point, so the permutation is only well-defined up to
        conjugation. One can read off the ramification type from the length of the
        disjoint cycles.

        EXAMPLES:

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z, w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: S.monodromy_group()
            [(0,2,1), (), (1,2), (0,2,1), (0,1,2), (0,2,1), (0,2), (1,2)]
            sage: zip(S.branch_locus, S.monodromy_group()) #abs tol 0.0000001
            [(0.000000000000000, (0,2,1)),
             (-1.31362670141929, ()),
             (-0.819032851784253 - 1.02703471138023*I, (1,2)),
             (-0.819032851784253 + 1.02703471138023*I, (0,2,1)),
             (0.292309440469772 - 1.28069133740100*I, (0,1,2)),
             (0.292309440469772 + 1.28069133740100*I, (0,2,1)),
             (1.18353676202412 - 0.569961265016465*I, (0,2)),
             (1.18353676202412 + 0.569961265016465*I, (1,2))]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: S.monodromy_group()
            [(0,1), (), (), ()]

        ..NOTE::

            The i-th generator corresponds to the i-th branch point.
            Later, we can implement a call which asks what
            the monodromy around a point is and return the corresponding
            generator.
        """
        # Computing the set of permutations associated with each edge.
        edgesd = self.downstairs_edges()
        edgesu = self.upstairs_edges()
        PermSet = []
        for e in self.edge_permutations():
            PermSet += [e[1]]
        n = len(self.branch_locus)
        cycbasisd = [self.voronoi_diagram.regions[self.voronoi_diagram.point_region[i]] for i in range(n)]
        monodromygens = []
        perm = self._Sn(())
        for center_point, c in zip(self.branch_locus,cycbasisd):
            for i in range(len(c)-1):
                # My permutations are 'directed'- they depend on the direction
                # which you're moving. Since my PermSet was constructed so that
                # the permutations only correspond to 'lower going to higher',
                # I add this conditional checking which case I'm in.
                # If I have an edge like (5,3), it takes the inverse instead.
                if (c[i] < c[i+1]):
                    perm = perm*PermSet[c[i]]
                else:
                    perm = perm*(PermSet[c[i]])**(-1)
            if (c[-1] < c[0]):
                perm = perm*PermSet[c[-1]]
            else:
                perm = perm*(PermSet[c[-1]])**(-1)
            if Matrix([list(self._vertices[c[0]]-center_point), list(self._vertices[c[1]]-center_point)]).determinant() < 0:
                perm = perm**(-1)
            monodromygens += [perm]
        return monodromygens

    @cached_method
    def homology_basis(self):
        """
        Computes the homology basis of the Riemannian surface.

        OUTPUT:

        A list L, with entries being lists of paths `[P_1 ... P_n]`.

        Each path `P_i` is of the form `(k, [p_1 ... p_m, p_1])`, where
        k is the number of times to traverse the path (if negative, to
        traverse it backwards), and the `p_i` are vertices of the
        upstairs graph.

        EXAMPLES:

        Computing a homology basis.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: g = w^2 - z^4 + 1
            sage: S = RiemannSurface(g)
            sage: S.homology_basis()
            [[(1,
               [(8, 1),
                (1, 0),
                (0, 0),
                (4, 0),
                (12, 0),
                (5, 0),
                (2, 1),
                (3, 1),
                (9, 1),
                (12, 1),
                (8, 1)])],
             [(1,
               [(5, 1),
                (6, 1),
                (7, 1),
                (4, 0),
                (12, 0),
                (5, 0),
                (2, 1),
                (3, 1),
                (9, 1),
                (12, 1),
                (5, 1)])]]

        Hence there are two paths that form the homology basis for the surface.
        """
        # Fringe case.
        if self.genus == 0:
            return []
        edgesu = self.upstairs_edges()
        cycles = Graph(edgesu).cycle_basis()
        # Computing the Gram matrix.
        cn = len(cycles)
        # Forming a list of lists of zeroes. Later this will be converted into a
        # matrix.
        intersectionprod = [[0 for c in cycles] for c in cycles]
        # This loop will start at the entry (0,1), and proceed along the row up
        # til (0,cn-1).
        # Then it will go to entry (1,2), and proceed along the row, etc.
        for i in range(cn-1):
            for j in range(cn-1-i):
                # Initializing the intersection product value.
                intsum = 0
                # Intersection of the edges
                intsec = set(cycles[i]).intersection(set(cycles[j+i+1]))
                for v in intsec:
                    # Get indices of the vertex in the cycles.
                    i0 = cycles[i].index(v)
                    i1 = cycles[j+i+1].index(v)
                    # Get the complex value of the vertex v.
                    vd = self._vertices[cycles[i][i0][0]]
                    # Get the two adjacent vertices to v in each cycle. There is
                    #some modular arithmetic here for if v is the last vertex in
                    # the cycle, so that the next adjacent vertex is the first
                    # element in the list. Also, I just skip straight to the
                    # downstairs vertices since that's what we care about.
                    vds = [cycles[i][i0-1][0],
                           cycles[i][(i0+1)%len(cycles[i])][0],
                           cycles[j+i+1][i1-1][0],
                           cycles[j+i+1][(i1+1)%len(cycles[j+i+1])][0]]
                    # The angles of the vectors.
                    cds = [(self._vertices[vds[0]] - vd).argument(),
                           (self._vertices[vds[1]] - vd).argument(),
                           (self._vertices[vds[2]] - vd).argument(),
                           (self._vertices[vds[3]] - vd).argument()]
                    # Above, a_in = cds[0], a_out = cds[1]
                    #        b_in = cds[2], b_out = cds[3].
                    # Lots of cases.
                    if cds[0] < cds[1]:
                        # 1
                        if (cds[0] < cds[2] < cds[1]):
                            intsum += 1
                        # 2
                        if ((cds[2] < cds[0]) or (cds[1] < cds[2])):
                            intsum -= 1
                        # 3
                        if (cds[0] < cds[3] < cds[1]):
                            intsum -= 1
                        # 4
                        if ((cds[3] < cds[0]) or (cds[1] < cds[3])):
                            intsum += 1
                    else:
                        assert(cds[1] < cds[0])
                        # 5
                        if (cds[1] < cds[2] < cds[0]):
                            intsum -= 1
                        # 6
                        if ((cds[2] < cds[1]) or (cds[0] < cds[2])):
                            intsum += 1
                        # 7
                        if (cds[1] < cds[3] < cds[0]):
                            intsum += 1
                        # 8
                        if ((cds[3] < cds[1]) or (cds[0] < cds[3])):
                            intsum -= 1
                assert (intsum%2) == 0
                intsum = intsum//2
                intersectionprod[i][j+i+1] = intsum
                # Skew Symmetry
                intersectionprod[j+i+1][i] = -intsum
        Gmatrix = Matrix(intersectionprod)
        P = Gmatrix.symplectic_form()[1]
        # Define the cycle sets.
        acycles = [[] for i in range(self.genus)]
        bcycles = [[] for i in range(self.genus)]
        # There are g a and b cycles.
        for i in range(self.genus):
            # Range over the size of the Gram matrix.
            for j in range(cn):
                # Forms the acycles and bcycles. If the entry in the
                # transformation matrix is non-zero, it adds the coefficient at
                # that entry, and the corresponding cycle. (also, forms it
                # into a loop)
                if P[i][j] != 0:
                    acycles[i] += [(P[i][j],[x for x in cycles[j]]+[cycles[j][0]])]
                if P[self.genus + i][j] != 0:
                    bcycles[i] += [(P[self.genus + i][j],[x for x in cycles[j]]+[cycles[j][0]])]
        return acycles + bcycles

    def make_zw_interpolator(self,upstairs_edge):
        """
        Given an upstairs edge for which continuation data has been stored,
        return a function that computes `z(t),w(t)` , where t in `[0,1]` is a
        parametrization of the edge.

        Also returns the Delta_z of the path as a second return argument.

        INPUT:

        - ``upstairs_edge`` -- A pair of integer tuples indicating an edge
          on the upstairs graph of the surface. Of the form [(a,b), (c,d)]

        OUTPUT:

        A tuple (g, d), where g is the function that computes the interpolation
        along the edge and d is the difference of the z-values of the end and start point.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f)
            sage: _ = S.homology_basis()
            sage: g,d = S.make_zw_interpolator([(0,0),(1,0)]);
            sage: all(f(*g(i*0.1)).abs() < 1e-13for i in range(10))
            True
            sage: abs((g(1)[0]-g(0)[0]) - d) < 1e-13
            True

        """
        eindex = tuple(u[0] for u in upstairs_edge)
        i0, i1 = eindex
        z_start = self._vertices[i0]
        z_end = self._vertices[i1]
        currL = self._L[eindex]
        windex = upstairs_edge[0][1]
        def w_interpolate(t):
            if t < 0 or t > 1:
                raise ValueError("t outside path range")
            if t == 0:
                return z_start,currL[0][1][windex]
            elif t == 1:
                return z_end,currL[-1][1][windex]
            while True:
                i = bisect(currL,t)
                t1, w1 ,epsilon = currL[i]
                w1 = w1[windex]
                t2, w2, _ = currL[i+1]
                w2 = w2[windex]
                z0 = (1-t)*z_start+t*z_end
                w0 = self._CC(((t2-t)*w1+(t-t1)*w2)/(t2-t1))
                try:
                    desired_result = self._newton_iteration(z0,w0,epsilon)
                except ConvergenceError:
                    pass
                else:
                    return z0,desired_result
                #If we did not succeed, we insert a new point in our interpolation list
                tnew=t
                while True:
                    tnew = (t1 + tnew)/2
                    znew = (1-tnew)*self._vertices[i0]+tnew*self._vertices[i1]
                    try:
                        neww1 = self._determine_new_w(znew,currL[i][1],epsilon)
                    except ConvergenceError:
                        pass
                    else:
                        #When *no* ConvergenceError is raised, we have succeeded and we can exit
                        break
                #once the loop has succeeded we insert our new value
                t1 = tnew
                self._L[eindex].insert(i+1,(t1,neww1,epsilon))
        return w_interpolate,(z_end-z_start)

    def simple_vector_line_integral(self,upstairs_edge,differentials):
        """
        Perfoms vectorized integration along a straight path.

        Suppose the Riemann surface is given by the polynomial equation
        f(z,w)=0. We consider differentials of the form g(z,w)/(df/dw)*dz,
        and the differentials are given by specifying the polynomial g.

        INPUT:

        - ``upstairs_edge`` -- a pair of integer tuples corresponding to an
          edge of the upstairs graph.

        - ``differentials`` -- a list of polynomials.

        OUTPUT:

        A complex number from evaluating the line integral.

        EXAMPLES:

        Setting up a Riemann surface.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = w^2 - z^4 + 1
            sage: S = RiemannSurface(f); S
            Riemann surface defined by polynomial f = -z^4 + w^2 + 1 = 0, with 53 bits of precision

        Since we make use of data from homotopy continuation, we need to compute
        the necessary data.::

            sage: M = S.riemann_matrix()
            sage: differentials = S.cohomology_basis()
            sage: S.simple_vector_line_integral([(0,0),(1,0)], differentials) #abs tol 0.00000001
            (3.25260651745651e-19 + 0.0777066972231166*I)

        ..NOTE:

            This integration makes use of data which was previously computed
            during homotopy continuation.
        """
        w_of_t,Delta_z = self.make_zw_interpolator(upstairs_edge)
        V = VectorSpace(self._CC,self.genus)
        def integrand(t):
            zt,wt = w_of_t(t)
            dfdwt = self._fastcall_dfdw(zt,wt)
            return V([omega(zt,wt)/dfdwt for omega in differentials])

        I=integrate_vector(integrand,self._prec)*Delta_z
        return I

    def cohomology_basis(self,option=1):
        """
        Computes the cohomology basis of the surface.

        INPUT:

        - ``option`` -- Presently, this routine uses Singular's adjointIdeal
            and passes the ``option`` parameter on. Legal values are 1, 2, 3 ,4,
            where 1 is the default. See the Singular documentation for the meaning.
            The backend for this function may change, and support for this parameter may
            disappear.

        OUTPUT:

        Returns a list of polynomials g representing the holomorphic differentials g/(df/dw)*dz,
        where f(z,w)=0 is the equation specifying the Riemann surface.

        EXAMPLES:

        Computing the cohomology basis.::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f)
            sage: S.cohomology_basis()
            [1, w, z]
        """
        if self.genus == 0:
            self._differentials = [[],option]
            return self._differentials[0]
        if (self._differentials == None) or (option != self._differentials[1]):
            # Computes the adjoint ideal using singular.
            k=QQ['z,w,u']
            singular.lib("paraplanecurves.lib")
            fnew = self.f(k.gen(0)/k.gen(2),k.gen(1)/k.gen(2)).numerator()
            L = singular(fnew).adjointIdeal(option).sage()
            adjointideal = L.gens()
            rightdeg = []
            # We need to take those generators which are of degree f.degree() - 3.
            # For those generators which are exactly of this degree, we do exactly that.
            # Otherwise, we need to multiply by monomials of appropriate degree.
            # This is what the below code does.
            belowdeg = filter(lambda x: x.degree() < self.f.degree() - 3, adjointideal)
            for x in belowdeg:
                # This bit of code computes all of the monomials of degree
                # less than or equal to tempdeg.
                # tempdeg is just the gap between degrees + 1.
                tempdeg = self.f.degree() - 2 - x.degree()
                monomials = []
                for i in range(tempdeg):
                    for j in range(i+1):
                        monomials += [self._R.gen(0)**(i-j)*self._R.gen(1)**j]
                rightdeg += [x(self._R.gen(0),self._R.gen(1),1)*m for m in monomials]
            # Then we add all of those generators which are precisely the right degree.
            rightdeg += [self._R(g(self._R.gen(0),self._R.gen(1),1)) 
                         for g in filter(lambda x: x.degree() == self.f.degree() - 3, adjointideal)]
            # The only linear dependency condition I check. If two generators are equal, I only need one.
            rightdegnew = set()
            for m in rightdeg:
                rightdegnew.add(m)
            rightdeg = list(rightdegnew)
            rightdeg.sort(key = lambda m: m.degree())
            # A sanity check. It's possible that I may have some linear dependencies of generators
            # leading to some phony results.. Currently, this is a lot of effort to resolve.
            if len(rightdeg) == self.genus:
                self._differentials = [rightdeg, option]
            else:
                raise ValueError('Length of differentials list is not equal to genus.')
        return self._differentials[0]

    def matrix_of_integral_values(self,differentials):
        """
        Compute the path integrals of the given differentials along the homology basis

        The returned answer has a row for each differential. If the Riemann surface is
        given by the equation f(z,w)=0, then the differentials are encoded by polynomials
        g, signifying the differential g(z,w)/(df/dw) * dz, where f(z,w)=0 is the equation
        defining the Riemann surface.

        INPUT:

        - ``differentials`` -- a list of polynomials.

        OUTPUT:

        A matrix, one row per differential, containing the values of the path integrals along
        the homology basis of the Riemann surface.

        EXAMPLES::
            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S=RiemannSurface(x^3+y^3+1)
            sage: B=S.cohomology_basis()
            sage: S.matrix_of_integral_values(B) #abs tol 1e-12
            [1.76663875028545 - 5.55111512312578e-17*I    0.883319375142723 + 1.52995403705719*I]


        """
        cycles = self.homology_basis()
        def normalize_pairs(L):
            """
            Returns a list of edges encoded by the path in L.
            The edges are normalized to be in the direction in which
            the homotopy continuation should have been computed along them.
            """
            R=[]
            for i in range(len(L)-1):
                if L[i][0]<L[i+1][0]:
                    R.append((L[i],L[i+1]))
                else:
                    R.append((L[i+1],L[i]))
            return R
        occurring_edges = set()
        occurring_edges.update(*[normalize_pairs(p[1]) for h in cycles for p in h])
        integral_dict=dict()
        for upstairs_edge in occurring_edges:
            integral_dict[upstairs_edge]=self.simple_vector_line_integral(upstairs_edge,differentials)
        rows=[]
        for cycle in cycles:
            V = VectorSpace(self._CC,self.genus).zero()
            for multiplicity,loop in cycle:
                for i in range(len(loop)-1):
                    if loop[i][0]<loop[i+1][0]:
                        direction=1
                        upstairs_edge=(loop[i],loop[i+1])
                    else:
                        direction=-1
                        upstairs_edge=(loop[i+1],loop[i])
                    V+=(multiplicity*direction)*integral_dict[upstairs_edge]
            rows.append(V)
        return Matrix(rows).transpose()

    def period_matrix(self,option=1):
        """
        Computes the period matrix of f given the homology basis and cohomology
        basis.

        INPUT:

        - ``option`` -- the computation options for singular's adjointIdeal
          computations. Default 1, options 1, 2, 3, 4. Changing this may make
          the computation faster.

        OUTPUT:

        A matrix of complex values.

        EXAMPLES:

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=200)
            sage: S.period_matrix() # abs tol 0.000001
            [-0.775197805903... - 0.61819962132...*I  0.96665585280... - 0.770882318822...*I    1.74185365871... - 1.38908194015...*I                  0.430202326362... + 0.893324335546...*I -0.238744279457... - 0.495757604603...*I    2.94725379097... + 1.11395722593...*I]
            [  2.17205598507... - 0.495757604603...*I  0.966655852808... + 0.220632890384...*I  -1.20540013226... - 0.275124714218...*I                1.74185365871... + 1.38908194015...*I    1.39685817917... + 1.11395722593...*I  -1.63560245862... - 0.618199621328...*I]
            [  0.536453526445... + 1.11395722593...*I   0.966655852808... - 2.00728156147...*I  0.430202326362... - 0.893324335546...*I               -1.20540013226... + 0.275124714218...*I   2.70850951152... - 0.618199621328...*I  -1.31165133234... - 0.495757604603...*I]
        """
        if (self._PM == None) or (self._differentials[1] != option):
            differentials = self.cohomology_basis(option=option)
            differentials = [fast_callable(omega,domain=self._CC) for omega in differentials]
            self._PM=self.matrix_of_integral_values(differentials)
        return self._PM

    def riemann_matrix(self, option=1):
        """
        Computes the Riemann matrix.

        INPUT:

        - ``option`` -- the computation options for singular's adjointIdeal
          computations. Default 1, options 1, 2, 3, 4. Changing this may make
          the computation faster.

        OUTPUT:

        A matrix of complex values.

        EXAMPLES:

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<z,w> = QQ[]
            sage: f = z^3*w + w^3 + z
            sage: S = RiemannSurface(f, prec=200)
            sage: S.riemann_matrix() #abs tol 0.0000001
            [-0.125000000000... + 0.992156741649...*I  0.375000000000... - 0.330718913883...*I -0.750000000000... + 0.661437827766...*I]
            [ 0.375000000000... - 0.330718913883...*I  0.874999999999... + 0.992156741649...*I -0.249999999999... - 0.661437827766...*I]
            [-0.750000000000... + 0.661437827766...*I -0.250000000000... - 0.661437827766...*I   0.499999999999... + 1.32287565553...*I]

        """
        # I won't cache this since the computations in it aren't so bad.
        PeriodMatrix = self.period_matrix(option=option)
        Am = PeriodMatrix[0:self.genus,0:self.genus]
        self._RM = (Am.inverse())*PeriodMatrix[0:self.genus,self.genus:2*self.genus]
        # Below are checks to see if an error was made.
        #posdef = matrix(self.genus,self.genus, [x.imag_part() for x in self._RM.list()]).eigenvalues()
        #sym = (self._RM - self._RM.transpose()).list()
        #for x in posdef:
        #    if x < 0:
        #        print 'Computation failed: resulting matrix is not positive definite.'
        #        break
        #for x in sym:
        #    if abs(x) > 10**(-5):
        #        print 'Computation failed: resulting matrix is not symmetric.'
        #        break
        return self._RM

    def plot_paths(self):
        """
        Make a graphical representation of the integration paths

        The plot returns a plot containing the branch points (in red)
        and the integration paths (obtained from the Voronoi cells of the
        branch points). The integration paths are plotted by plotting the points
        that have been computed for homotopy continuation, so the density
        gives an indication of where numerically sensitive features occur.

        OUTPUT:

        A 2d graphical object.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2-x^3-x)
            sage: S.plot_paths()
            Graphics object consisting of 2 graphics primitives

        """
        from sage.plot.point import point2d
        P=[]
        #trigger the computation of the homology basis, so that self._L is present
        _=self.homology_basis()
        for e in self._L.keys():
            z0=self._vertices[e[0]]
            z1=self._vertices[e[1]]
            def path(t):
                z=(1-t)*z0+t*z1
                return z
            T=self._L[e]
            P+=[path(t[0]) for t in T]
        plt=point2d(P,size=1)+point2d(self.branch_locus,color="red")
        return plt

    def plot_paths3d(self):
        """
        Return the homology basis as a graph in 3-space

        The homology basis of the surface is constructing by taking the Voronoi cells
        around the branch points and taking the inverse image of the edges on the
        Riemann surface. If the surface is given by the equation f(z,w), the returned
        object gives the image of this graph in 3-space with coordinates
        (Re(z),Im(z),Im(w)).

        OUTPUT:

        EXAMPLE::
            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S = RiemannSurface(y^2-x^3-x)
            sage: S.plot_paths3d()
            Graphics3d Object

        """
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes2 import point3d, line3d
        P = Graphics()
        #trigger the computation of the homology basis, so that self._L is present
        _=self.homology_basis()

        for e in self._L.keys():
            z0 = self._vertices[e[0]]
            z1 = self._vertices[e[1]]
            def path(t):
                z = (1-t)*z0+t*z1
                return (z.real_part(),z.imag_part())
            T = self._L[e]
            color = "blue"
            for i in range(self.degree):
                P += line3d([path(t[0])+(t[1][i].imag_part(),) for t in T],color=color,thickness=0.01)
            for z,ws in zip(self._vertices,self._wvalues):
                for w in ws:
                    P += point3d([z.real_part(),z.imag_part(),w.imag_part()],color="purple", size=20)
        return P

    def endomorphism_basis(self,b=None,r=None):
        """
        Numerical computation of a ZZ-basis for the endomorphism ring.

        Let (I|M) be the normalized period matrix (M is the gxg Riemann matrix).
        We consider the system of matrix equations
        MA+B=(MC+D)M
        where A,B,C,D are gxg integer matrices.
        We determine small integer (near) solutions using LLL reductions.
        These solutions are returned as 2g x 2g integer matrices obtained by
        stacking (D|B) on top of (C|A).

        INPUT:

        - ``b`` -- integer (default: (precision - 10)). The equation coefficients are
            scaled by 2^b before rounding to integers.

        - ``r`` -- integer (default: b/4). Solutions that are have all coefficients smaller than
            2^r in absolute value are reported as actual solutions.

        OUTPUT:

        A list of 2g x 2g integer matrices that, for large enough r and (b-r),
        generate the endomorphism ring.

        EXAMPLES::

            sage: from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
            sage: R.<x,y> = QQ[]
            sage: S=RiemannSurface(x^3+y^3+1)
            sage: S.endomorphism_basis()
            [
            [1 0]  [ 0 -1]
            [0 1], [ 1  1]
            ]

        """
        M = self.riemann_matrix()
        H = max(max( abs(m.real_part()) for m in M.list()), max( abs(m.imag_part()) for m in M.list()))
        if b is None:
            b = self._prec - 10
        if r is None:
            r = b//4
        S = 2**b
        if H*S > 2**(self._prec-5):
            raise ValueError("Insufficient precision for b=%s"%b)
        g = M.ncols()
        CC = M.base_ring()
        V = ["%s%s"%(n,i) for n in ["a","b","c","d"] for i in srange(1,1+g**2)]
        R = PolynomialRing(CC,V)
        A = Matrix(R,g,g,V[:g**2])
        B = Matrix(R,g,g,V[g**2:2*g**2])
        C = Matrix(R,g,g,V[2*g**2:3*g**2])
        D = Matrix(R,g,g,V[3*g**2:4*g**2])
        # Given the normalized period matrix ( I | M ) we multiply on the right by
        # the integer matric (D, B ; C, A) to get the result (D+MC | B+MA).
        # Bringing that matrix in normalized form gives (I | (D+MC)^(-1)(B+MA)).
        # Equating this to (I|M) and clearing denominators gives the equations
        # below.
        W = ((M*A+B) - (M*C+D)*M).list()
        vars = R.gens()
        #mt = Matrix(ZZ,[[1 if i==j else 0 for j in range(4*g**2)] +
        #  [((S*w.monomial_coefficient(vars[i]).real_part()).round() for w in W] +
        #  [S*w.monomial_coefficient(vars[i]).imag_part()).round() for w in W] for i in range(len(vars))])
        mt = Matrix(ZZ,[[1 if i==j else 0 for j in range(4*g**2)] +
          [(S*w.monomial_coefficient(vars[i]).real_part()).round() for w in W] +
          [(S*w.monomial_coefficient(vars[i]).imag_part()).round() for w in W] for i in range(len(vars))])
        #we compute an LLL-reduced basis of this lattice:
        mtL = mt.LLL()
        def vectomat(v,g):
            A = Matrix(g,g,v[:g**2].list())
            B = Matrix(g,g,v[g**2:2*g**2].list())
            C = Matrix(g,g,v[2*g**2:3*g**2].list())
            D = Matrix(g,g,v[3*g**2:4*g**2].list())
            return D.augment(B).stack(C.augment(A))
        c = 2**r
        return [vectomat(v,g) for v in mtL if all(abs(a) <= c for a in v[4*g**2:])]

r"""
Sets of Morphisms between Differentiable Manifolds

The class :class:`DifferentiableManifoldHomset` implements sets of morphisms between
two differentiable manifolds over the same topological field `K` (in most
applications, `K = \RR` or `K = \CC`), a morphism being a *differentiable map*
for the category of differentiable manifolds.

The subclass :class:`DifferentiableCurveSet` is devoted to the specific case of
differential curves, i.e. morphisms whose domain is an open interval of
`\RR`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [Lee2013]_
- [KN1963]_

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.manifolds.manifold_homset import TopologicalManifoldHomset
from sage.manifolds.differentiable.diff_map import DiffMap
from sage.manifolds.differentiable.curve import DifferentiableCurve
from sage.manifolds.differentiable.integrated_curve import IntegratedCurve
from sage.manifolds.differentiable.integrated_curve import IntegratedAutoparallelCurve
from sage.manifolds.differentiable.integrated_curve import IntegratedGeodesic
from sage.misc.cachefunc import cached_method

class DifferentiableManifoldHomset(TopologicalManifoldHomset):
    r"""
    Set of differentiable maps between two differentiable manifolds.

    Given two differentiable manifolds `M` and `N` over a topological field `K`,
    the class :class:`DifferentiableManifoldHomset` implements the set
    `\mathrm{Hom}(M,N)` of morphisms (i.e. differentiable maps)
    `M\rightarrow N`.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.manifolds.differentiable.diff_map.DiffMap`.

    INPUT:

    - ``domain`` -- differentiable manifold `M` (domain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
    - ``codomain`` -- differentiable manifold `N` (codomain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
    - ``name`` -- (default: ``None``) string; name given to the homset; if
      ``None``, Hom(M,N) will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      homset; if ``None``, `\mathrm{Hom}(M,N)` will be used

    EXAMPLES:

    Set of differentiable maps between a 2-dimensional differentiable manifold
    and a 3-dimensional one::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: N = Manifold(3, 'N')
        sage: Y.<u,v,w> = N.chart()
        sage: H = Hom(M, N) ; H
        Set of Morphisms from 2-dimensional differentiable manifold M to
         3-dimensional differentiable manifold N in Category of smooth
         manifolds over Real Field with 53 bits of precision
        sage: type(H)
        <class 'sage.manifolds.differentiable.manifold_homset.DifferentiableManifoldHomset_with_category'>
        sage: H.category()
        Category of homsets of topological spaces
        sage: latex(H)
        \mathrm{Hom}\left(M,N\right)
        sage: H.domain()
        2-dimensional differentiable manifold M
        sage: H.codomain()
        3-dimensional differentiable manifold N

    An element of ``H`` is a differentiable map from ``M`` to ``N``::

        sage: H.Element
        <class 'sage.manifolds.differentiable.diff_map.DiffMap'>
        sage: f = H.an_element() ; f
        Differentiable map from the 2-dimensional differentiable manifold M to the
         3-dimensional differentiable manifold N
        sage: f.display()
        M --> N
           (x, y) |--> (u, v, w) = (0, 0, 0)

    The test suite is passed::

        sage: TestSuite(H).run()

    When the codomain coincides with the domain, the homset is a set of
    *endomorphisms* in the category of differentiable manifolds::

        sage: E = Hom(M, M) ; E
        Set of Morphisms from 2-dimensional differentiable manifold M to
         2-dimensional differentiable manifold M in Category of smooth
         manifolds over Real Field with 53 bits of precision
        sage: E.category()
        Category of endsets of topological spaces
        sage: E.is_endomorphism_set()
        True
        sage: E is End(M)
        True

    In this case, the homset is a monoid for the law of morphism composition::

        sage: E in Monoids()
        True

    This was of course not the case for ``H = Hom(M, N)``::

        sage: H in Monoids()
        False

    The identity element of the monoid is of course the identity map of ``M``::

        sage: E.one()
        Identity map Id_M of the 2-dimensional differentiable manifold M
        sage: E.one() is M.identity_map()
        True
        sage: E.one().display()
        Id_M: M --> M
           (x, y) |--> (x, y)

    The test suite is passed by ``E``::

        sage: TestSuite(E).run()

    This test suite includes more tests than in the case of ``H``, since ``E``
    has some extra structure (monoid).

    """

    Element = DiffMap

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M, N) ; H
            Set of Morphisms from 2-dimensional differentiable manifold M to
             3-dimensional differentiable manifold N in Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: TestSuite(H).run()

        Test for an endomorphism set::

            sage: E = Hom(M, M) ; E
            Set of Morphisms from 2-dimensional differentiable manifold M to
             2-dimensional differentiable manifold M in Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: TestSuite(E).run()

        """
        from sage.manifolds.differentiable.manifold import \
                                                         DifferentiableManifold
        if not isinstance(domain, DifferentiableManifold):
            raise TypeError("domain = {} is not an ".format(domain) +
                            "instance of DifferentiableManifold")
        if not isinstance(codomain, DifferentiableManifold):
            raise TypeError("codomain = {} is not an ".format(codomain) +
                            "instance of DifferentiableManifold")
        TopologicalManifoldHomset.__init__(self, domain, codomain, name=name,
                                           latex_name=latex_name)

    #### Parent methods ####

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M,N)
            sage: H._coerce_map_from_(ZZ)
            False
            sage: H._coerce_map_from_(M)
            False
            sage: H._coerce_map_from_(N)
            False

        """
        #!# for the time being:
        return False

    #### End of parent methods ####


#******************************************************************************

class DifferentiableCurveSet(DifferentiableManifoldHomset):
    r"""
    Set of differentiable curves in a differentiable manifold.

    Given an open interval `I` of `\RR` (possibly `I = \RR`) and
    a differentiable manifold `M` over `\RR`, this is the set
    `\mathrm{Hom}(I,M)` of morphisms (i.e. differentiable curves) `I \to M`.

    INPUT:

    - ``domain`` --
      :class:`~sage.manifolds.differentiable.real_line.OpenInterval`
      if an open interval `I \subset \RR` (domain of the morphisms),
      or :class:`~sage.manifolds.differentiable.real_line.RealLine`
      if `I = \RR`
    - ``codomain`` --
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
      differentiable manifold `M` (codomain of the morphisms)
    - ``name`` -- (default: ``None``) string; name given to the set of
      curves; if ``None``, ``Hom(I, M)`` will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the set of curves; if ``None``, `\mathrm{Hom}(I,M)` will be used

    EXAMPLES:

    Set of curves `\RR \longrightarrow M`, where `M` is a 2-dimensional
    manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: R.<t> = RealLine() ; R
        Real number line R
        sage: H = Hom(R, M) ; H
        Set of Morphisms from Real number line R to 2-dimensional
         differentiable manifold M in Category of smooth manifolds over Real
         Field with 53 bits of precision
        sage: H.category()
        Category of homsets of topological spaces
        sage: latex(H)
        \mathrm{Hom}\left(\Bold{R},M\right)
        sage: H.domain()
        Real number line R
        sage: H.codomain()
        2-dimensional differentiable manifold M

    An element of ``H`` is a curve in ``M``::

        sage: c = H.an_element(); c
        Curve in the 2-dimensional differentiable manifold M
        sage: c.display()
        R --> M
           t |--> (x, y) = (1/(t^2 + 1) - 1/2, 0)

    The test suite is passed::

        sage: TestSuite(H).run()

    The set of curves `(0,1) \longrightarrow U`, where `U` is an open
    subset of `M`::

        sage: I = R.open_interval(0, 1) ; I
        Real interval (0, 1)
        sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1}) ; U
        Open subset U of the 2-dimensional differentiable manifold M
        sage: H = Hom(I, U) ; H
        Set of Morphisms from Real interval (0, 1) to Open subset U of the
         2-dimensional differentiable manifold M in Join of Category of
         subobjects of sets and Category of smooth manifolds over Real Field
         with 53 bits of precision

    An element of ``H`` is a curve in ``U``::

        sage: c = H.an_element() ; c
        Curve in the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: c.display()
        (0, 1) --> U
           t |--> (x, y) = (1/(t^2 + 1) - 1/2, 0)

    The set of curves `\RR \longrightarrow \RR` is a set of (manifold)
    endomorphisms::

        sage: E = Hom(R, R) ; E
        Set of Morphisms from Real number line R to Real number line R in
         Category of smooth manifolds over Real Field with 53 bits of precision
        sage: E.category()
        Category of endsets of topological spaces
        sage: E.is_endomorphism_set()
        True
        sage: E is End(R)
        True

    It is a monoid for the law of morphism composition::

        sage: E in Monoids()
        True

    The identity element of the monoid is the identity map of `\RR`::

        sage: E.one()
        Identity map Id_R of the Real number line R
        sage: E.one() is R.identity_map()
        True
        sage: E.one().display()
        Id_R: R --> R
           t |--> t

    A "typical" element of the monoid::

        sage: E.an_element().display()
        R --> R
           t |--> 1/(t^2 + 1) - 1/2

    The test suite is passed by ``E``::

        sage: TestSuite(E).run()

    Similarly, the set of curves `I \longrightarrow I` is a monoid, whose
    elements are (manifold) endomorphisms::

        sage: EI = Hom(I, I) ; EI
        Set of Morphisms from Real interval (0, 1) to Real interval (0, 1) in
         Join of Category of subobjects of sets and
             Category of smooth manifolds over Real Field with 53 bits of precision
        sage: EI.category()
        Category of endsets of subobjects of sets and topological spaces
        sage: EI is End(I)
        True
        sage: EI in Monoids()
        True

    The identity element and a "typical" element of this monoid::

        sage: EI.one()
        Identity map Id_(0, 1) of the Real interval (0, 1)
        sage: EI.one().display()
        Id_(0, 1): (0, 1) --> (0, 1)
           t |--> t
        sage: EI.an_element().display()
        (0, 1) --> (0, 1)
           t |--> 1/2/(t^2 + 1) + 1/4

    The test suite is passed by ``EI``::

        sage: TestSuite(EI).run()

    """
    Element = DifferentiableCurve

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: H = Hom(R, M); H
            Set of Morphisms from Real number line R to 3-dimensional
             differentiable manifold M in Category of smooth manifolds over
             Real Field with 53 bits of precision
            sage: TestSuite(H).run()
            sage: Hom(R, M) is Hom(R, M)
            True
            sage: H = Hom(R, R); H
            Set of Morphisms from Real number line R to Real number line R in
             Category of smooth manifolds over Real Field with 53 bits of
             precision
            sage: TestSuite(H).run()
            sage: I = R.open_interval(-1, 2)
            sage: H = Hom(I, M); H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Join of Category of subobjects of
             sets and Category of smooth manifolds over Real Field with 53 bits
             of precision
            sage: TestSuite(H).run()
            sage: H = Hom(I, I); H
            Set of Morphisms from Real interval (-1, 2) to Real interval (-1, 2)
             in Join of Category of subobjects of sets and Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: TestSuite(H).run()

        """
        from sage.manifolds.differentiable.real_line import OpenInterval
        if not isinstance(domain, OpenInterval):
            raise TypeError("{} is not an open real interval".format(domain))
        DifferentiableManifoldHomset.__init__(self, domain, codomain, name=name,
                                              latex_name=latex_name)

    #### Parent methods ####

    def _element_constructor_(self, coord_expression, name=None,
                              latex_name=None, is_isomorphism=False,
                              is_identity=False):
        r"""
        Construct an element of ``self``, i.e. a differentiable curve
        `I \to M`, where `I` is a real interval and `M` some
        differentiable manifold.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine() ; R
            Real number line R
            sage: H = Hom(R, M)
            sage: c = H({X: [sin(t), sin(2*t)/2]}, name='c') ; c
            Curve c in the 2-dimensional differentiable manifold M
            sage: c = Hom(R, R)({}, is_identity=True) ; c
            Identity map Id_R of the Real number line R

        """
        # Standard construction
        return self.element_class(self, coord_expression=coord_expression,
                                  name=name, latex_name=latex_name,
                                  is_isomorphism=is_isomorphism,
                                  is_identity=is_identity)

    def _an_element_(self):
        r"""
        Construct some element of ``self``.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = Hom(R,M)._an_element_() ; c
            Curve in the 3-dimensional differentiable manifold M
            sage: c.display()
            R --> M
               t |--> (x, y, z) = (1/(t^2 + 1) - 1/2, 0, 0)

        ::

            sage: I = R.open_interval(0, pi)
            sage: c = Hom(I,M)._an_element_() ; c
            Curve in the 3-dimensional differentiable manifold M
            sage: c.display()
            (0, pi) --> M
               t |--> (x, y, z) = (1/(t^2 + 1) - 1/2, 0, 0)

        ::

            sage: c = Hom(I,I)._an_element_() ; c
            Differentiable map from the Real interval (0, pi) to itself
            sage: c.display()
            (0, pi) --> (0, pi)
               t |--> 1/4*pi + 1/2*pi/(t^2 + 1)

        """
        from sage.rings.infinity import Infinity
        from sage.rings.rational_field import QQ
        dom = self.domain()
        codom = self.codomain()
        # A simple curve is constructed around a point of the codomain:
        chart2 = codom.default_chart()
        target_point = chart2.domain().an_element()
        target_coord = list(target_point.coord(chart2))
        bounds = chart2._bounds[0]  # bounds of first coordinate
        # Determination of an interval (x1, x2) arround target_point:
        xmin = bounds[0][0]
        xmax = bounds[1][0]
        one_half = QQ(1) / QQ(2)
        if xmin == -Infinity:
            if xmax == Infinity:
                x1 = - one_half
                x2 = one_half
            else:
                x1 = xmax - 3*one_half
                x2 = xmax - one_half
        else:
            if xmax == Infinity:
                x1 = xmin + one_half
                x2 = xmin + 3*one_half
            else:
                dx = (xmax - xmin) / 4
                x1 = xmin + dx
                x2 = xmax - dx
        # The coordinate function defining the curve:
        t = dom.canonical_coordinate()
        target_coord[0] = x1 + (x2-x1) / (1+t*t)
        coord_expression = {chart2: target_coord}
        return self.element_class(self, coord_expression)

#******************************************************************************

class IntegratedCurveSet(DifferentiableCurveSet):
    r"""
    Set of integrated curves in a differentiable manifold.

    INPUT:

    - ``domain`` --
      :class:`~sage.manifolds.differentiable.real_line.OpenInterval`
      open interval `I \subset \RR` with finite boundaries (domain of
      the morphisms)
    - ``codomain`` --
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
      differentiable manifold `M` (codomain of the morphisms)
    - ``name`` -- (default: ``None``) string; name given to the set of
      integrated curves; if ``None``, ``Hom_integrated(I, M)`` will be
      used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the set of integrated curves; if ``None``,
      `\mathrm{Hom_{integrated}}(I,M)` will be used

    EXAMPLE:

    This parent class needs to be imported::

        sage: from sage.manifolds.differentiable.manifold_homset import IntegratedCurveSet

    Integrated curves are only allowed to be defined on an interval with
    finite bounds.
    This forbids to define an instance of this parent class whose domain
    has infinite bounds::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: R.<t> = RealLine()
        sage: H = IntegratedCurveSet(R, M)
        Traceback (most recent call last):
        ...
        ValueError: Both boundaries of the interval defining the domain
         of a Homset of integrated curves need to be finite.

    An instance whose domain is an interval with finite bounds allows to
    build an integrated curve defined on the interval::

        sage: I = R.open_interval(-1, 2)
        sage: H = IntegratedCurveSet(I, M) ; H
        Set of Morphisms from Real interval (-1, 2) to 2-dimensional
         differentiable manifold M in Category of homsets of subobjects
         of sets and topological spaces which actually are integrated
         curves
        sage: eqns_rhs = [1,1]
        sage: vels = X.symbolic_velocities()
        sage: t = var('t')
        sage: p = M.point((3,4))
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((1,2))
        sage: c = H(eqns_rhs, vels, t, v, name='c') ; c
        Integrated curve c in the 2-dimensional differentiable
         manifold M

    A "typical" element of ``H`` is a curve in ``M``::

        sage: d = H.an_element(); d
        Integrated curve in the 2-dimensional differentiable manifold M
        sage: sys = d.system(verbose=True)
        Curve in the 2-dimensional differentiable manifold M integrated
         over the Real interval (-1, 2) as a solution to the following
         system, written w.r.t. Chart (M, (x, y)):
        <BLANKLINE>
        Initial point: Point on the 2-dimensional differentiable
         manifold M with coordinates [0, 0] w.r.t. Chart (M, (x, y))
        Initial tangent vector: Tangent vector at Point on the
         2-dimensional differentiable manifold M with components
         [1/4, 0] w.r.t. Chart (M, (x, y))
        <BLANKLINE>
        d(x)/dt = Dx
        d(y)/dt = Dy
        d(Dx)/dt = -1/4*sin(t + 1)
        d(Dy)/dt = 0
        <BLANKLINE>

    The test suite is passed::

        sage: TestSuite(H).run()

    More generally, an instance of this class may be defined with
    abstract bounds :MATH:`(a,b)`::

        sage: [a,b] = var('a b')
        sage: J = R.open_interval(a, b)
        sage: H = IntegratedCurveSet(J, M) ; H
        Set of Morphisms from Real interval (a, b) to 2-dimensional
         differentiable manifold M in Category of homsets of subobjects
         of sets and topological spaces which actually are integrated
         curves

    A "typical" element of ``H`` is a curve in ``M``::

        sage: f = H.an_element(); f
        Integrated curve in the 2-dimensional differentiable manifold M
        sage: sys = f.system(verbose=True)
        Curve in the 2-dimensional differentiable manifold M integrated
         over the Real interval (a, b) as a solution to the following
         system, written w.r.t. Chart (M, (x, y)):
        <BLANKLINE>
        Initial point: Point on the 2-dimensional differentiable
         manifold M with coordinates [0, 0] w.r.t. Chart (M, (x, y))
        Initial tangent vector: Tangent vector at Point on the
         2-dimensional differentiable manifold M with components
         [1/4, 0] w.r.t. Chart (M, (x, y))
        <BLANKLINE>
        d(x)/dt = Dx
        d(y)/dt = Dy
        d(Dx)/dt = -1/4*sin(-a + t)
        d(Dy)/dt = 0
        <BLANKLINE>

    Yet, even in the case of abstract bounds, considering any of them to
    be infinite is still prohibited since no numerical integration could
    be performed::

        sage: f.solve(parameters_values={a:-1, b:+oo})
        Traceback (most recent call last):
        ...
        ValueError: Both boundaries of the interval need to be finite.

    The set of integrated curves `J \longrightarrow J` is a set of
    numerical (manifold) endomorphisms::

        sage: H = IntegratedCurveSet(J, J); H
        Set of Morphisms from Real interval (a, b) to Real interval
         (a, b) in Category of endsets of subobjects of sets and
         topological spaces which actually are integrated curves
        sage: H.category()
        Category of endsets of subobjects of sets and topological spaces

    It is a monoid for the law of morphism composition::

        sage: H in Monoids()
        True

    The identity element of the monoid is a numerical version of the
    identity map of `J`::

        sage: e = H.one() ; e
        Integrated curve Id_(a, b) in the Real interval (a, b)

    A "typical" element of the monoid::

        sage: g = H.an_element() ; g
        Integrated curve in the Real interval (a, b)
        sage: sys = g.system(verbose=True)
        Curve in the Real interval (a, b) integrated over the Real interval (a, b) as a solution to the following system, written w.r.t. Chart ((a, b), (t,)):
        <BLANKLINE>
        Initial point: Point on the Real number line R with coordinates [0] w.r.t. Chart ((a, b), (t,))
        Initial tangent vector: Tangent vector at Point on the Real number line R with components [1/4] w.r.t. Chart ((a, b), (t,))
        <BLANKLINE>
        d(t)/ds = Dt
        d(Dt)/ds = -1/4*sin(-a + s)
        <BLANKLINE>

    The test suite is passed::

        sage: TestSuite(H).run()

    """

    Element = IntegratedCurve

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedCurveSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: H = IntegratedCurveSet(R, M)
            Traceback (most recent call last):
            ...
            ValueError: Both boundaries of the interval defining the
             domain of a Homset of integrated curves need to be finite.
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedCurveSet(I, M) ; H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Category of homsets of
             subobjects of sets and topological spaces which actually
             are integrated curves
            sage: TestSuite(H).run()
            sage: H = IntegratedCurveSet(I, I); H
            Set of Morphisms from Real interval (-1, 2) to Real interval
             (-1, 2) in Category of endsets of subobjects of sets and
             topological spaces which actually are integrated curves
            sage: TestSuite(H).run()

        """

        from sage.rings.infinity import Infinity

        DifferentiableCurveSet.__init__(self, domain, codomain,
                                       name=name, latex_name=latex_name)

        # checking argument 'domain': 't_min' and 't_max' are only
        # allowed to be either expressions of finite real values
        t_min = domain.lower_bound()
        t_max = domain.upper_bound()
        if t_min == -Infinity or t_max == +Infinity:
            raise ValueError("Both boundaries of the interval " +
                             "defining the domain of a Homset of " +
                             "integrated curves need to be finite.")

        if name is None:
            self._name = "Hom_integrated({},{})".format(domain._name,
                                                         codomain._name)
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = r"\mathrm{{Hom}_{integrated}}"
            self._latex_name += r"\left({},{}\right)".format(
                               domain._latex_name, codomain._latex_name)
        else:
            self._latex_name = latex_name

    #### Parent methods ####

    def _repr_(self):
        """
        TESTS::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedCurveSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedCurveSet(I, M) ; H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Category of homsets of
             subobjects of sets and topological spaces which actually
             are integrated curves

        """
        description = "Set of Morphisms "
        description += "from {} to {} in {} ".format(self._domain,
                                        self._codomain, self.category())
        description += "which actually are integrated curves"
        return description


    def _element_constructor_(self, equations_rhs, velocities,
                 curve_parameter, initial_tangent_vector, chart=None,
                 name=None, latex_name=None, is_isomorphism=False,
                 is_identity=False, verbose=False):
        r"""
        Construct an element of ``self``, i.e. an integrated curve
        `I \to M`, where `I` is a real interval and `M` some
        differentiable manifold.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedCurve`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedCurveSet
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedCurveSet(I, M)
            sage: eqns_rhs = [1,1]
            sage: vels = X.symbolic_velocities()
            sage: t = var('t')
            sage: p = M.point((3,4))
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,2))
            sage: c = H(eqns_rhs, vels, t, v, name='c') ; c
            Integrated curve c in the 2-dimensional differentiable
             manifold M

        """
        # Standard construction
        return self.element_class(self, equations_rhs, velocities,
                curve_parameter, initial_tangent_vector, chart=chart,
                name=name, latex_name=latex_name,
                is_isomorphism=is_isomorphism, is_identity=is_identity,
                verbose=verbose)

    def _an_element_(self):
        r"""
        Construct some element of ``self``.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedCurve`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedCurveSet
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedCurveSet(I, M)
            sage: c = H._an_element_() ; c
            Integrated curve in the 2-dimensional differentiable
             manifold M
            sage: sys = c.system(verbose=True)
            Curve in the 2-dimensional differentiable manifold M
             integrated over the Real interval (-1, 2) as a solution to
             the following system, written w.r.t. Chart (M, (x, y)):
            <BLANKLINE>
            Initial point: Point on the 2-dimensional differentiable
             manifold M with coordinates [0, 0] w.r.t. Chart (M, (x, y))
            Initial tangent vector: Tangent vector at Point on the
             2-dimensional differentiable manifold M with components
             [1/4, 0] w.r.t. Chart (M, (x, y))
            <BLANKLINE>
            d(x)/dt = Dx
            d(y)/dt = Dy
            d(Dx)/dt = -1/4*sin(t + 1)
            d(Dy)/dt = 0
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: c(1)
            [0.22732435599328793, 0.0]
            sage: H = IntegratedCurveSet(I, I)
            sage: c = H._an_element_() ; c
            Integrated curve in the Real interval (-1, 2)
            sage: sys = c.system(verbose=True)
            Curve in the Real interval (-1, 2) integrated over the Real
             interval (-1, 2) as a solution to the following system,
             written w.r.t. Chart ((-1, 2), (t,)):
            <BLANKLINE>
            Initial point: Point on the Real number line R with
             coordinates [1/2] w.r.t. Chart ((-1, 2), (t,))
            Initial tangent vector: Tangent vector at Point on the Real
             number line R with components [3/8] w.r.t.
             Chart ((-1, 2), (t,))
            <BLANKLINE>
            d(t)/ds = Dt
            d(Dt)/ds = -3/8*sin(s + 1)
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: c(1)
            [0.840986533989932]

        """

        from sage.categories.homset import Hom
        from sage.functions.trig import sin
        from sage.functions.trig import cos
        from sage.symbolic.ring import var
        from sage.symbolic.expression import Expression

        dom = self.domain()
        t = dom.canonical_coordinate()
        t_min = dom.lower_bound() # this is either an expression or a
        # finite value thanks to tests in '__init__'
        t_max = dom.upper_bound() # this is either an expression or a
        # finite value thanks to tests in '__init__'

        codom = self.codomain()
        dim = codom.dim()
        chart2 = codom.default_chart()
        # In case the codomain coincides with the domain,
        # it is important to distinguish between the canonical
        # coordinate, and the curve parameter since, in such a
        # situation, the coordinate should not be used to denote the
        # curve parameter, since it actually becomes a function of the
        # curve parameter, and such a function is an unknown of the
        # system defining the curve.
        # In other cases, it might still happen for a coordinate of the
        # codomain to be denoted the same as the canonical coordinate of
        # the domain (for instance, the codomain could be another
        # real interval, different from the domain, and yet with same
        # letter denoting its canonical coordinate).
        # In such case, an error is raised from method 'init'
        # of class IntegratedCurve; to solve it, the user is
        # free to change the name of the codomain coordinate in the
        # chart used on the codomain.
        if dom == codom:
            param = var('s')
            if t == param: # the canonical coordinate of the domain
            # might be the expression 's' even though it was affected
            # above to the variable 't'
                param = var('u')
        else:
            param = t

        # An analytical curve is used to find a region of the codomain
        # where a certain integrated curve may be defined:
        H = Hom(dom, codom)
        c = H.an_element()
        x0_A = c.expr()[0].substitute({t:1})
        x0_B = c.expr()[0].substitute({t:0}) # necessarily, x0_A < x0_B
        p_coords = [x0_A] + list(c.expr()[1:dim])
        p = codom.point(p_coords)

        # The initial tangent vector:
        v_comps = [(x0_B-x0_A)/2] + [0 for i in range(dim-1)]
        v = codom.tangent_space(p)(v_comps)

        # The equations defining the curve:
        eqns_rhs=[-(x0_B-x0_A)/2*sin(param-t_min)]+[0 for i in range(dim-1)]
        # combined with the initial components above, all velocities
        # vanish, except the first one, which is a cosine function.
        # This differential system results in a curve constant in all
        # its coordinates, except the first one, which oscillates around
        # the value 'x0_A' with an amplitude '(x0_B-x0_A)/2'

        # The symbolic expressions for the velocities:
        vels = chart2.symbolic_velocities()

        return self.element_class(self,eqns_rhs,vels,param,v)

    @cached_method
    def one(self):
        r"""
        Return the identity element of ``self`` considered as a monoid
        (case of a set of endomorphisms).

        This applies only when the codomain of the homset is equal to its
        domain, i.e. when the homset is of the type `\mathrm{Hom}(M,M)`.
        Indeed, `\mathrm{Hom}(M,M)` equipped with the law of morphisms
        composition is a monoid, whose identity element is nothing but the
        identity map of `M`.

        OUTPUT:

        - the identity map of `M`, as an instance of
          :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedCurve`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedCurveSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedCurveSet(I, I); H
            Set of Morphisms from Real interval (-1, 2) to Real interval
             (-1, 2) in Category of endsets of subobjects of sets and
             topological spaces which actually are integrated curves
            sage: e = H.one() ; e
            Integrated curve Id_(-1, 2) in the Real interval (-1, 2)

        """

        from sage.symbolic.ring import var

        if self.codomain() != self.domain():
            raise TypeError("{} is not a monoid".format(self))

        t = self.domain().canonical_coordinate()
        # It is important to distinguish between the canonical
        # coordinate, and the curve parameter since, in such a
        # situation, the coordinate should not be used to denote the
        # curve parameter, since it actually becomes a function of the
        # curve parameter, and such a function is an unknown of the
        # system defining the curve.
        param = var('s')
        if t == param: # the canonical coordinate of the domain
        # might be the expression 's' even though it was affected
        # above to the variable 't'
            param = var('u')

        vel = [var("D{}".format(t))]

        # Below, set argument to 'None' where the value does not have
        # any influence since argument 'is_identity' is set to 'True'
        return self.element_class(self, None, vel, param, None,
                                                       is_identity=True)

#******************************************************************************

class IntegratedAutoparallelCurveSet(IntegratedCurveSet):
    r"""
    Set of integrated autoparallel curves in a differentiable manifold.

    INPUT:

    - ``domain`` --
      :class:`~sage.manifolds.differentiable.real_line.OpenInterval`
      open interval `I \subset \RR` with finite boundaries (domain of
      the morphisms)
    - ``codomain`` --
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
      differentiable manifold `M` (codomain of the morphisms)
    - ``name`` -- (default: ``None``) string; name given to the set of
      integrated autoparallel curves; if ``None``,
      ``Hom_autoparallel(I, M)`` will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the set of integrated autoparallel curves; if ``None``,
      `\mathrm{Hom_{autoparallel}}(I,M)` will be used

    EXAMPLE:

    This parent class needs to be imported::

        sage: from sage.manifolds.differentiable.manifold_homset import IntegratedAutoparallelCurveSet

    Integrated autoparallel curves are only allowed to be defined on an
    interval with finite bounds.
    This forbids to define an instance of this parent class whose domain
    has infinite bounds::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: R.<t> = RealLine()
        sage: H = IntegratedAutoparallelCurveSet(R, M)
        Traceback (most recent call last):
        ...
        ValueError: Both boundaries of the interval defining the domain
         of a Homset of integrated autoparallel curves need to be finite.

    An instance whose domain is an interval with finite bounds allows to
    build a curve that is autoparallel w.r.t a connection defined on the
    codomain::

        sage: I = R.open_interval(-1, 2)
        sage: H = IntegratedAutoparallelCurveSet(I, M) ; H
        Set of Morphisms from Real interval (-1, 2) to 2-dimensional
         differentiable manifold M in Category of homsets of subobjects
         of sets and topological spaces which actually are integrated
         autoparallel curves w.r.t a certain affine connection
        sage: nab = M.affine_connection('nabla')
        sage: nab[0,1,0], nab[0,0,1] = 1,2
        sage: nab.torsion()[:]
        [[[0, -1], [1, 0]], [[0, 0], [0, 0]]]
        sage: t = var('t')
        sage: p = M.point((3,4))
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((1,2))
        sage: c = H(nab, t, v, name='c') ; c
        Integrated autoparallel curve c in the 2-dimensional
         differentiable manifold M

    A "typical" element of ``H`` is an autoparallel curve in ``M``::

        sage: d = H.an_element(); d
        Integrated autoparallel curve in the 2-dimensional
         differentiable manifold M
        sage: sys = d.system(verbose=True)
        Autoparallel curve in the 2-dimensional differentiable manifold
         M equipped with Affine connection nab on the Open subset U of
         the 2-dimensional differentiable manifold M, and integrated
         over the Real interval (-1, 2) as a solution to the following
         equations, written w.r.t. Chart (M, (x, y)):
        <BLANKLINE>
        Initial point: Point on the 2-dimensional differentiable
         manifold M with coordinates [0.0312500000000000, -1/2] w.r.t.
         Chart (M, (x, y))
        Initial tangent vector: Tangent vector at Point on the
         2-dimensional differentiable manifold M with components
         [0.156250000000000, 2.07698341289763] w.r.t. Chart (M, (x, y))
        <BLANKLINE>
        d(x)/dt = Dx
        d(y)/dt = Dy
        d(Dx)/dt = 0
        d(Dy)/dt = -Dx*Dy*cos(x)/sin(x)
        <BLANKLINE>

    The test suite is passed::

        sage: TestSuite(H).run()

    For any open interval `J` with finite bounds `(a,b)`, all curves are
    autoparallel w.r.t any connection.
    Therefore, the set of autoparallel curves `J \longrightarrow J` is a set of
    numerical (manifold) endomorphisms that is a monoid for the law of
    morphism composition::

        sage: [a,b] = var('a b')
        sage: J = R.open_interval(a, b)
        sage: H = IntegratedAutoparallelCurveSet(J, J); H
        Set of Morphisms from Real interval (a, b) to Real interval
         (a, b) in Category of endsets of subobjects of sets and
         topological spaces which actually are integrated autoparallel
         curves w.r.t a certain affine connection
        sage: H.category()
        Category of endsets of subobjects of sets and topological spaces
        sage: H in Monoids()
        True
        sage: e = H.one() ; e
        Integrated autoparallel curve Id_(a, b) in the Real
         interval (a, b)

    A "typical" element of the monoid::

        sage: g = H.an_element() ; g
        Integrated autoparallel curve in the Real interval (a, b)
        sage: sys = g.system(verbose=True)
        Autoparallel curve in the Real interval (a, b) equipped with Affine connection nab on the Open subset U of the Real number line R, and integrated over the Real interval (a, b) as a solution to the following equations, written w.r.t. Chart ((a, b), (t,)):
        <BLANKLINE>
        Initial point: Point on the Real number line R with coordinates [1/200] w.r.t. Chart ((a, b), (t,))
        Initial tangent vector: Tangent vector at Point on the Real number line R with components [-9999/400/(a - b)] w.r.t. Chart ((a, b), (t,))
        <BLANKLINE>
        d(t)/ds = Dt
        d(Dt)/ds = -Dt^2/t
        <BLANKLINE>

    The test suite is passed::

        sage: TestSuite(H).run()

    """

    Element = IntegratedAutoparallelCurve

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedAutoparallelCurveSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: H = IntegratedAutoparallelCurveSet(R, M)
            Traceback (most recent call last):
            ...
            ValueError: Both boundaries of the interval defining the
             domain of a Homset of integrated autoparallel curves need
             to be finite.
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedAutoparallelCurveSet(I, M) ; H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Category of homsets of
             subobjects of sets and topological spaces which actually are
             integrated autoparallel curves w.r.t a certain affine connection
            sage: TestSuite(H).run()
            sage: H = IntegratedAutoparallelCurveSet(I, I); H
            Set of Morphisms from Real interval (-1, 2) to Real interval
             (-1, 2) in Category of endsets of subobjects of sets and
             topological spaces which actually are integrated
             autoparallel curves w.r.t a certain affine connection
            sage: TestSuite(H).run()

        """

        from sage.rings.infinity import Infinity

        DifferentiableCurveSet.__init__(self, domain, codomain,
                                       name=name, latex_name=latex_name)

        # checking argument 'domain'
        t_min = domain.lower_bound()
        t_max = domain.upper_bound()
        if t_min == -Infinity or t_max == +Infinity:
            raise ValueError("Both boundaries of the interval " +
                             "defining the domain of a Homset of " +
                             "integrated autoparallel curves need to " +
                             "be finite.")

        if name is None:
            self._name = "Hom_autoparallel"
            self._name += "({},{})".format(domain._name, codomain._name)
        else:
            self._name = name
        if latex_name is None:
            self._latex_name=r"\mathrm{{Hom}_{autoparallel}}"
            self._latex_name+= r"\left({},{}\right)".format(
                               domain._latex_name, codomain._latex_name)
        else:
            self._latex_name = latex_name

    #### Parent methods ####

    def _repr_(self):
        """
        TESTS::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedAutoparallelCurveSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedAutoparallelCurveSet(I, M) ; H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Category of homsets of
             subobjects of sets and topological spaces which actually
             are integrated autoparallel curves w.r.t a certain affine
             connection

        """

        description = "Set of Morphisms "
        description += "from {} to {} in {} ".format(self._domain,
                                        self._codomain, self.category())
        description += "which actually are integrated autoparallel "
        description += "curves w.r.t a certain affine connection"
        return description


    def _element_constructor_(self, affine_connection, curve_parameter,
                    initial_tangent_vector, chart=None, name=None,
                    latex_name=None, is_isomorphism=False,
                    is_identity=False, verbose=False):
        r"""
        Construct an element of ``self``, i.e. an integrated
        autoparallel curve `I \to M`, where `I` is a real interval and
        `M` some differentiable manifold.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurve`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedAutoparallelCurveSet
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedAutoparallelCurveSet(I, M)
            sage: nab = M.affine_connection('nabla')
            sage: nab[0,1,0], nab[0,0,1] = 1,2
            sage: nab.torsion()[:]
            [[[0, -1], [1, 0]], [[0, 0], [0, 0]]]
            sage: t = var('t')
            sage: p = M.point((3,4))
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,2))
            sage: c = H(nab, t, v, name='c') ; c
            Integrated autoparallel curve c in the 2-dimensional
             differentiable manifold M

        """
        # Standard construction
        return self.element_class(self, affine_connection,
                 curve_parameter, initial_tangent_vector, chart=chart,
                 name=name,latex_name=latex_name,
                 is_isomorphism=is_isomorphism, is_identity=is_identity,
                 verbose=verbose)

    def _an_element_(self):
        r"""
        Construct some element of ``self``.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurve`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedAutoparallelCurveSet
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedAutoparallelCurveSet(I, M)
            sage: c = H._an_element_() ; c
            Integrated autoparallel curve in the 2-dimensional
             differentiable manifold M
            sage: sys = c.system(verbose=True)
            Autoparallel curve in the 2-dimensional differentiable
             manifold M equipped with Affine connection nab on the Open
             subset U of the 2-dimensional differentiable manifold M,
             and integrated over the Real interval (-1, 2) as a solution
             to the following equations, written w.r.t.
             Chart (M, (x, y)):
            <BLANKLINE>
            Initial point: Point on the 2-dimensional differentiable
             manifold M with coordinates [0.0312500000000000, -1/2]
             w.r.t. Chart (M, (x, y))
            Initial tangent vector: Tangent vector at Point on the
             2-dimensional differentiable manifold M with components
             [0.156250000000000, 2.07698341289763]
             w.r.t. Chart (M, (x, y))
            <BLANKLINE>
            d(x)/dt = Dx
            d(y)/dt = Dy
            d(Dx)/dt = 0
            d(Dy)/dt = -Dx*Dy*cos(x)/sin(x)
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: c(1)
            [0.34375000000000056, 0.49999813529956155]
            sage: H = IntegratedAutoparallelCurveSet(I, I)
            sage: c = H._an_element_() ; c
            Integrated autoparallel curve in the Real interval (-1, 2)
            sage: sys = c.system(verbose=True)
            Autoparallel curve in the Real interval (-1, 2) equipped
             with Affine connection nab on the Open subset U of the Real
             number line R, and integrated over the Real
             interval (-1, 2) as a solution to the following equations,
             written w.r.t. Chart ((-1, 2), (t,)):
            <BLANKLINE>
            Initial point: Point on the Real number line R with
             coordinates [1/2] w.r.t. Chart ((-1, 2), (t,))
            Initial tangent vector: Tangent vector at Point on the Real
             number line R with components [7/16] w.r.t.
             Chart ((-1, 2), (t,))
            <BLANKLINE>
            d(t)/ds = Dt
            d(Dt)/ds = -Dt^2/t
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: c(1)
            [1.0606601601128764]

        """

        from sage.categories.homset import Hom
        from sage.symbolic.ring import var
        from sage.functions.trig import sin
        from sage.functions.trig import tan
        from sage.functions.log import ln
        from sage.symbolic.constants import pi
        from sage.rings.infinity import Infinity
        from sage.rings.rational_field import QQ

        dom = self.domain()
        t = dom.canonical_coordinate()
        t_min = dom.lower_bound() # this is either an expression or a
        # finite value thanks to tests in '__init__'
        t_max = dom.upper_bound() # idem

        codom = self.codomain()
        dim = codom.dim()
        i0 = codom.start_index()
        chart2 = codom.default_chart()
        th = chart2[:][0]
        # In case the codomain coincides with the domain,
        # it is important to distinguish between the canonical
        # coordinate, and the curve parameter since, in such a
        # situation, the coordinate should not be used to denote the
        # curve parameter, since it actually becomes a function of the
        # curve parameter, and such a function is an unknown of the
        # system defining the curve.
        # In other cases, it might still happen for a coordinate of the
        # codomain to be denoted the same as the canonical coordinate of
        # the domain (for instance, the codomain could be another
        # real interval, different from the domain, and yet with same
        # letter denoting its canonical coordinate).
        # In such case, an error is raised from method 'init'
        # of class IntegratedCurve; to solve it, the user is
        # free to change the name of the codomain coordinate in the
        # chart used on the codomain.
        if dom == codom:
            param = var('s')
        else:
            param = t

        # An analytical curve is used to find a region of the codomain
        # where a certain integrated autoparallel curve may be defined:
        H = Hom(dom, codom)
        c = H.an_element()
        th_A = c.expr()[0].substitute({t:1})
        th_B = c.expr()[0].substitute({t:0}) # necessarily, th_A < th_B

        if dim == 1:
            # Redefine 'th_A' and 'th_B' so that they are of the same
            # sign and such that th_A < th_B still holds.
            # This will avoid division by zero in the equations
            # defining autoparallel curves w.r.t. the connection defined
            # below
            if th_A.numerical_approx().sign() != th_B.numerical_approx().sign():
            # knowing that th_A < th_B
                if th_B > 0:
                    th_A = th_B/100
                else: # th_A is necessarily strictly less than 0
                    th_B = th_A/100

            # The initial point:
            p = codom.point([th_A])

            # The initial tangent vector:
            th_dot_A = (th_B**2-th_A**2)/(2*th_A*(t_max-t_min))
            v = codom.tangent_space(p)([th_dot_A])

            U = codom.open_subset('U')
            chart2_U = chart2.restrict(U, [th > th_A, th < th_B]) # this
            # is to prevent from problems arising from a division by
            # zero in the equations defining autoparallel curves w.r.t.
            # the connection defined below
            nab = U.affine_connection('nab')
            nab.set_coef()[i0,i0,i0] = 1/th
            return self.element_class(self, nab, param, v) # the
            # autoparallel curve returned will correspond to the
            # following analytical solution:
            # th(t) = (th_A*(th_A + 2*th_dot_A*(t-t_min)))^(1/2)

        # else: (i.e. dim >= 2)

        # in dimension greater than 1, the idea is to consider the first
        # two coordinates of the chart to be the polar coordinates on the
        # unit 2-sphere, in order to reproduce autoparallel curves on
        # the sphere with respect to the Mercator connection (which is
        # described in the example introducing the class
        # IntegratedAutoparallelCurve)

        # Redefine 'th_A' and 'th_B' so that [th_A, th_B] is contained
        # within an interval (k*pi, (k+1)*pi), k being an integer.
        # This is needed to take the inverse of the sine on [th_A, th_B]
        # and to avoid problems from evaluating the tangent of 'th/2',
        # which appear in the equations:

        # th_A = k_A*pi + r_A, k_A integer, 0 <= r_A < pi
        # idem for th_B
        r_A = th_A.numerical_approx()%pi.numerical_approx()
        r_B = th_B.numerical_approx()%pi.numerical_approx()
        k_A = (th_A-r_A)/pi
        k_B = (th_B-r_B)/pi
        if k_A == k_B and r_A == 0:
            th_A = k_A*pi + (th_B-th_A)/16
        elif k_A != k_B:
            th_B = (k_A+1)*pi - 1/16*((k_A+1)*pi - th_A)
            if r_A == 0:
                th_A = k_A*pi + (th_B-th_A)/16

        bounds = chart2._bounds[1]  # bounds of second coordinate
        # Determination of an interval (x1, x2) arround target_point:
        ph_min = bounds[0][0]
        ph_max = bounds[1][0]
        one_half = QQ(1) / QQ(2)
        if ph_min == -Infinity:
            if ph_max == Infinity:
                ph_A = - one_half
                ph_B = one_half
            else:
                ph_A = ph_max - 3*one_half
                ph_B = ph_max - one_half
        else:
            if ph_max == Infinity:
                ph_A = ph_min + one_half
                ph_B = ph_min + 3*one_half
            else:
                dph = (ph_max - ph_min) / 4
                ph_A = ph_min + dph
                ph_B = ph_max - dph

        # The initial point:
        p_coords = [th_A] + [ph_A] + list(c.expr()[2:dim])
        p = codom.point(p_coords)

        # The initial tangent vector:
        th_dot_A = (th_B - th_A)/(t_max - t_min)
        num_tan_alpha = ln(tan(th_A/2))-ln(tan((th_dot_A*t_max+th_A)/2))
        denom_tan_alpha = ph_B - ph_A
        tan_alpha = num_tan_alpha/denom_tan_alpha
        ph_dot_A = - th_dot_A/(tan_alpha * sin(th_A))
        v_comps = [th_dot_A] + [ph_dot_A] + [0 for i in range(dim-2)]
        v = codom.tangent_space(p)(v_comps)

        # The Mercator connection is reproduced on the first two
        # coordinates
        U = codom.open_subset('U')
        chart2_U = chart2.restrict(U, [th > th_A, th < th_B]) # this is
        # to prevent from problems arising from defining the orhtonormal
        # basis 'epolar_ON' with the inverse of a sine
        epolar_U = chart2_U.frame()
        ch_basis = U.automorphism_field()
        ch_basis[i0,i0], ch_basis[i0+1,i0+1] = 1, 1/sin(th)
        if dim > 2:
            for i in range(2,dim):
                ch_basis[i0+i,i0+i] = 1
        epolar_ON_U = epolar_U.new_frame(ch_basis, 'epolar_ON')
        nab = U.affine_connection('nab') # the connection is only
        # defined on U since it was defined with respect to epolar_ON_U,
        # which was defined on U only
        nab.set_coef(epolar_ON_U)[:] # this sets all the coefficients to
        # zero w.r.t. 'epolar_ON'

        return self.element_class(self, nab, param, v)

    @cached_method
    def one(self):
        r"""
        Return the identity element of ``self`` considered as a monoid
        (case of a set of endomorphisms).

        This applies only when the codomain of the homset is equal to its
        domain, i.e. when the homset is of the type `\mathrm{Hom}(M,M)`.
        Indeed, `\mathrm{Hom}(M,M)` equipped with the law of morphisms
        composition is a monoid, whose identity element is nothing but the
        identity map of `M`.

        OUTPUT:

        - the identity map of `M`, as an instance of
          :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurve`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedAutoparallelCurveSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedAutoparallelCurveSet(I, I); H
            Set of Morphisms from Real interval (-1, 2) to Real interval
             (-1, 2) in Category of endsets of subobjects of sets and
             topological spaces which actually are integrated
             autoparallel curves w.r.t a certain affine connection
            sage: e = H.one() ; e
            Integrated autoparallel curve Id_(-1, 2) in the Real
             interval (-1, 2)

        """

        from sage.symbolic.ring import var

        if self.codomain() != self.domain():
            raise TypeError("{} is not a monoid".format(self))

        t = self.domain().canonical_coordinate()
        # It is important to distinguish between the canonical
        # coordinate, and the curve parameter since, in such a
        # situation, the coordinate should not be used to denote the
        # curve parameter, since it actually becomes a function of the
        # curve parameter, and such a function is an unknown of the
        # system defining the curve.
        param = var('s')
        if t == param: # the canonical coordinate of the domain
        # might be the expression 's' even though it was affected
        # above to the variable 't'
            param = var('u')

        vel = [var("D{}".format(t))]

        # Below, set argument to 'None' where the value does not have any
        # influence since argument 'is_identity' is set to 'True'
        # As a result, it does not need to be overriden in
        # IntegratedGeodesicSet since the signature of the 'init'
        # methods of the classes IntegratedAutoparallelCurve and
        # IntegratedGeodesic are the same, except argument
        # 'affine_connection' of the first one is replaced by 'metric'
        # in the second one.
        # Since both do not have any influence on the definition of the
        # curve when argument 'is_identity' is set to 'True', the class
        # IntegratedGeodesicSet may use this method 'one' with arguments
        # set to 'None' in the same places.
        return self.element_class(self, None, param, None,
                                                       is_identity=True)


#******************************************************************************

class IntegratedGeodesicSet(IntegratedAutoparallelCurveSet):
    r"""
    Set of integrated geodesic in a differentiable manifold.

    INPUT:

    - ``domain`` --
      :class:`~sage.manifolds.differentiable.real_line.OpenInterval`
      open interval `I \subset \RR` with finite boundaries (domain of
      the morphisms)
    - ``codomain`` --
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
      differentiable manifold `M` (codomain of the morphisms)
    - ``name`` -- (default: ``None``) string; name given to the set of
      integrated geodesics; if ``None``, ``Hom_geodesic(I, M)`` will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the set of integrated geodesics; if ``None``,
      `\mathrm{Hom_{geodesic}}(I,M)` will be used

    EXAMPLE:

    This parent class needs to be imported::

        sage: from sage.manifolds.differentiable.manifold_homset import IntegratedGeodesicSet

    Integrated geodesics are only allowed to be defined on an interval
    with finite bounds.
    This forbids to define an instance of this parent class whose domain
    has infinite bounds::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: R.<t> = RealLine()
        sage: H = IntegratedGeodesicSet(R, M)
        Traceback (most recent call last):
        ...
        ValueError: Both boundaries of the interval defining the domain
         of a Homset of integrated geodesics need to be finite.

    An instance whose domain is an interval with finite bounds allows to
    build a geodesic w.r.t. a metric defined on the codomain::

        sage: I = R.open_interval(-1, 2)
        sage: H = IntegratedGeodesicSet(I, M) ; H
        Set of Morphisms from Real interval (-1, 2) to 2-dimensional
         differentiable manifold M in Category of homsets of subobjects
         of sets and topological spaces which actually are integrated
         geodesics w.r.t a certain metric
        sage: g = M.metric('g')
        sage: g[0,0], g[1,1], g[0,1] = 1, 1, 2
        sage: t = var('t')
        sage: p = M.point((3,4))
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((1,2))
        sage: c = H(g, t, v, name='c') ; c
        Integrated geodesic c in the 2-dimensional differentiable
         manifold M

    A "typical" element of ``H`` is a geodesic in ``M``::

        sage: d = H.an_element(); d
        Integrated geodesic in the 2-dimensional differentiable
         manifold M
        sage: sys = d.system(verbose=True)
        Geodesic in the 2-dimensional differentiable manifold M equipped
         with Riemannian metric g on the 2-dimensional differentiable
         manifold M, and integrated over the Real interval (-1, 2) as a
         solution to the following geodesic equations, written w.r.t.
         Chart (M, (x, y)):
        <BLANKLINE>
        Initial point: Point on the 2-dimensional differentiable
         manifold M with coordinates [0.0312500000000000, 0] w.r.t.
         Chart (M, (x, y))
        Initial tangent vector: Tangent vector at Point on the
         2-dimensional differentiable manifold M with components
         [0.156250000000000, 0] w.r.t. Chart (M, (x, y))
        <BLANKLINE>
        d(x)/dt = Dx
        d(y)/dt = Dy
        d(Dx)/dt = Dy^2*cos(x)*sin(x)
        d(Dy)/dt = -2*Dx*Dy*cos(x)/sin(x)

    The test suite is passed::

        sage: TestSuite(H).run()

    For any open interval `J` with finite bounds `(a,b)`, all curves are
    geodesics w.r.t any metric.
    Therefore, the set of geodesics `J \longrightarrow J` is a set of
    numerical (manifold) endomorphisms that is a monoid for the law of
    morphism composition::

        sage: [a,b] = var('a b')
        sage: J = R.open_interval(a, b)
        sage: H = IntegratedGeodesicSet(J, J); H
        Set of Morphisms from Real interval (a, b) to Real interval
         (a, b) in Category of endsets of subobjects of sets and
         topological spaces which actually are integrated geodesics
         w.r.t a certain metric
        sage: H.category()
        Category of endsets of subobjects of sets and topological spaces
        sage: H in Monoids()
        True
        sage: e = H.one() ; e
        Integrated geodesic Id_(a, b) in the Real interval (a, b)

    A "typical" element of the monoid::

        sage: g = H.an_element() ; g
        Integrated geodesic in the Real interval (a, b)
        sage: sys = g.system(verbose=True)
        Geodesic in the Real interval (a, b) equipped with Riemannian
         metric g on the Open subset U of the Real number line R, and
         integrated over the Real interval (a, b) as a solution to the
         following geodesic equations, written w.r.t.
         Chart ((a, b), (t,)):
        <BLANKLINE>
        Initial point: Point on the Real number line R with coordinates
         [1/200] w.r.t. Chart ((a, b), (t,))
        Initial tangent vector: Tangent vector at Point on the Real
         number line R with components [-9999999999/1000/(a - b)] w.r.t.
         Chart ((a, b), (t,))
        <BLANKLINE>
        d(t)/ds = Dt
        d(Dt)/ds = -4*Dt^2/t
        <BLANKLINE>

    The test suite is passed::

        sage: TestSuite(H).run()

    """

    Element = IntegratedGeodesic

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedGeodesicSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: H = IntegratedGeodesicSet(R, M)
            Traceback (most recent call last):
            ...
            ValueError: Both boundaries of the interval defining the
             domain of a Homset of integrated geodesics need to be
             finite.
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedGeodesicSet(I, M) ; H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Category of homsets of
             subobjects of sets and topological spaces which actually
             are integrated geodesics w.r.t a certain metric
            sage: TestSuite(H).run()
            sage: H = IntegratedGeodesicSet(I, I); H
            Set of Morphisms from Real interval (-1, 2) to Real interval
             (-1, 2) in Category of endsets of subobjects of sets and
             topological spaces which actually are integrated geodesics
             w.r.t a certain metric
            sage: TestSuite(H).run()

        """

        from sage.rings.infinity import Infinity

        DifferentiableCurveSet.__init__(self, domain, codomain,
                                       name=name, latex_name=latex_name)

        # checking argument 'domain'
        t_min = domain.lower_bound()
        t_max = domain.upper_bound()
        if t_min == -Infinity or t_max == +Infinity:
            raise ValueError("Both boundaries of the interval " +
                             "defining the domain of a Homset of " +
                             "integrated geodesics need to be finite.")

        if name is None:
            self._name = "Hom_geodesic"
            self._name += "({},{})".format(domain._name, codomain._name)
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = r"\mathrm{{Hom}_{geodesic}}"
            self._latex_name+= r"\left({},{}\right)".format(
                               domain._latex_name, codomain._latex_name)
        else:
            self._latex_name = latex_name

    #### Parent methods ####

    def _repr_(self):
        """
        TESTS::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedGeodesicSet
            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedGeodesicSet(I, M) ; H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Category of homsets of
             subobjects of sets and topological spaces which actually
             are integrated geodesics w.r.t a certain metric

        """
        description = "Set of Morphisms "
        description += "from {} to {} in {} ".format(self._domain,
                                        self._codomain, self.category())
        description += "which actually are integrated geodesics "
        description += "w.r.t a certain metric"
        return description

    def _element_constructor_(self, metric, curve_parameter,
                    initial_tangent_vector, chart=None, name=None,
                    latex_name=None, is_isomorphism=False,
                    is_identity=False, verbose=False):
        r"""
        Construct an element of ``self``, i.e. an integrated geodesic
        `I \to M`, where `I` is a real interval and
        `M` some differentiable manifold.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedGeodesic`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedGeodesicSet
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedGeodesicSet(I, M)
            sage: g = M.metric('g')
            sage: g[0,0], g[1,1], g[0,1] = 1, 1, 2
            sage: t = var('t')
            sage: p = M.point((3,4))
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,2))
            sage: c = H(g, t, v, name='c') ; c
            Integrated geodesic c in the 2-dimensional differentiable
             manifold M

        """
        # Standard construction
        return self.element_class(self, metric, curve_parameter,
                 initial_tangent_vector, chart=chart, name=name,
                 latex_name=latex_name, is_isomorphism=is_isomorphism,
                 is_identity=is_identity, verbose=verbose)

    def _an_element_(self):
        r"""
        Construct some element of ``self``.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedGeodesic`

        EXAMPLE::

            sage: from sage.manifolds.differentiable.manifold_homset import IntegratedGeodesicSet
            sage: M = Manifold(4, 'M', start_index=1)
            sage: X.<w,x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(-1, 2)
            sage: H = IntegratedGeodesicSet(I, M)
            sage: c = H._an_element_() ; c
            Integrated geodesic in the 4-dimensional differentiable
             manifold M
            sage: sys = c.system(verbose=True)
            Geodesic in the 4-dimensional differentiable manifold M
             equipped with Riemannian metric g on the 4-dimensional
             differentiable manifold M, and integrated over the Real
             interval (-1, 2) as a solution to the following geodesic
             equations, written w.r.t. Chart (M, (w, x, y, z)):
            <BLANKLINE>
            Initial point: Point on the 4-dimensional differentiable
             manifold M with coordinates [0.0312500000000000, 0, 0, 0]
             w.r.t. Chart (M, (w, x, y, z))
            Initial tangent vector: Tangent vector at Point on the
             4-dimensional differentiable manifold M with components
             [0.156250000000000, 0, 0, 0] w.r.t. Chart (M, (w, x, y, z))
            <BLANKLINE>
            d(w)/dt = Dw
            d(x)/dt = Dx
            d(y)/dt = Dy
            d(z)/dt = Dz
            d(Dw)/dt = Dx^2*cos(w)*sin(w)
            d(Dx)/dt = -2*Dw*Dx*cos(w)/sin(w)
            d(Dy)/dt = 0
            d(Dz)/dt = 0
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: c(1)
            [0.34375000000000056, 0.0, 0.0, 0.0]
            sage: H = IntegratedGeodesicSet(I, I)
            sage: c = H._an_element_() ; c
            Integrated geodesic in the Real interval (-1, 2)
            sage: sys = c.system(verbose=True)
            Geodesic in the Real interval (-1, 2) equipped with
             Riemannian metric g on the Open subset U of the Real number
             line R, and integrated over the Real interval (-1, 2) as a
             solution to the following geodesic equations, written
             w.r.t. Chart ((-1, 2), (t,)):
            <BLANKLINE>
            Initial point: Point on the Real number line R with
             coordinates [1/2] w.r.t. Chart ((-1, 2), (t,))
            Initial tangent vector: Tangent vector at Point on the Real
             number line R with components [1031/320] w.r.t.
             Chart ((-1, 2), (t,))
            <BLANKLINE>
            d(t)/ds = Dt
            d(Dt)/ds = -4*Dt^2/t
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: c(1)
            [1.1495758053215723]

        """

        from sage.categories.homset import Hom
        from sage.symbolic.ring import var
        from sage.functions.trig import sin
        from sage.symbolic.constants import pi

        dom = self.domain()
        t = dom.canonical_coordinate()
        t_min = dom.lower_bound() # this is either an expression or a
        # finite value thanks to tests in '__init__'
        t_max = dom.upper_bound() # idem

        codom = self.codomain()
        dim = codom.dim()
        i0 = codom.start_index()
        chart2 = codom.default_chart()
        th = chart2[:][0]
        # In case the codomain coincides with the domain,
        # it is important to distinguish between the canonical
        # coordinate, and the curve parameter since, in such a
        # situation, the coordinate should not be used to denote the
        # curve parameter, since it actually becomes a function of the
        # curve parameter, and such a function is an unknown of the
        # system defining the curve.
        # In other cases, it might still happen for a coordinate of the
        # codomain to be denoted the same as the canonical coordinate of
        # the domain (for instance, the codomain could be another
        # real interval, different from the domain, and yet with same
        # letter denoting its canonical coordinate).
        # In such case, an error is raised from method 'init'
        # of class IntegratedCurve; to solve it, the user is
        # free to change the name of the codomain coordinate in the
        # chart used on the codomain.
        if dom == codom:
            param = var('s')
        else:
            param = t

        # An analytical curve is used to find a region of the codomain
        # where a certain integrated autoparallel curve may be defined:
        H = Hom(dom, codom)
        c = H.an_element()
        th_A = c.expr()[0].substitute({t:1})
        th_B = c.expr()[0].substitute({t:0}) # necessarily, th_A < th_B

        if dim == 1:
            # Redefine 'th_A' and 'th_B' so that they are of the same
            # sign and such that th_A < th_B still holds.
            # This will avoid division by zero in the equations
            # defining geodesics w.r.t. the metric defined below.
            if th_A.numerical_approx().sign() != th_B.numerical_approx().sign():
            # knowing that th_A < th_B
                if th_B > 0:
                    th_A = th_B/100
                else: # th_A is necessarily strictly less than 0
                    th_B = th_A/100

            # The initial point:
            p = codom.point([th_A])

            # The initial tangent vector:
            th_dot_A = (th_B**5-th_A**5)/(5*(th_A**4)*(t_max-t_min))
            v = codom.tangent_space(p)([th_dot_A])

            U = codom.open_subset('U')
            chart2_U = chart2.restrict(U, [th > th_A, th < th_B]) # this
            # is to prevent from problems arising from a division by
            # zero in the equations defining autoparallel curves w.r.t.
            # the connection defined below
            g = U.metric('g')
            g[i0,i0] = th**8
            return self.element_class(self, g, param, v) # the
            # autoparallel curve returned will correspond to the
            # following analytical solution:
            # th(t) = (th_A^4*(th_A + 5*th_dot_A*(t-t_min)))^(1/5)

        # else: (i.e. dim >= 2)

        # in dimension greater than 1, the idea is to consider the first
        # two coordinates of the chart to be the polar coordinates on
        # the unit 2-sphere, in order to reproduce a line of longitude,
        # which is a geodesic on the sphere w.r.t. the standard metric
        # induced by the 3D Euclidean metric on the sphere.

        # Redefine 'th_A' and 'th_B' so that [th_A, th_B] is contained
        # within an interval (k*pi, (k+1)*pi), k being an integer.
        # This is needed to take the inverse of the tangent function on
        # [th_A, th_B] in the geodesic equations:

        # th_A = k_A*pi + r_A, k_A integer, 0 <= r_A < pi
        # idem for th_B
        r_A = th_A.numerical_approx()%pi.numerical_approx()
        r_B = th_B.numerical_approx()%pi.numerical_approx()
        k_A = (th_A-r_A)/pi
        k_B = (th_B-r_B)/pi
        if k_A == k_B and r_A == 0:
            th_A = k_A*pi + (th_B-th_A)/16
        elif k_A != k_B:
            th_B = (k_A+1)*pi - 1/16*((k_A+1)*pi - th_A)
            if r_A == 0:
                th_A = k_A*pi + (th_B-th_A)/16

        # The initial point:
        p_coords = [th_A] + list(c.expr()[1:dim])
        p = codom.point(p_coords)

        # The initial tangent vector:
        th_dot_A = (th_B - th_A)/(t_max - t_min)
        v_comps = [th_dot_A] + [0 for i in range(dim-1)]
        v = codom.tangent_space(p)(v_comps)

        # The metric of the unit 2-sphere is reproduced on the first two
        # coordinates
        g = codom.metric('g')
        g[i0,i0] = 1
        g[i0+1,i0+1] = sin(th)**2
        if dim > 2:
            for i in range(2,dim):
                g[i0+i,i0+i] = 1

        return self.element_class(self, g, param, v)

r"""
Points of topological manifolds

The class :class:`TopManifoldPoint` implements points of a
topological manifold, in a coordinate independent manner:
a :class:`TopManifoldPoint` object can have coordinates in
various charts defined on the manifold. Two points are declared equal if they
have the same coordinates in the same chart.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

REFERENCES:

- J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed., Springer (New
  York) (2011)
- M. Berger & B. Gostiaux: *Geometrie differentielle, varietes, courbes et
  surfaces*, Presses Universitaires de France (Paris, 1987)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York,
  2013)

EXAMPLES:

Defining a point on `\RR^3` by its spherical coordinates::

    sage: M = TopManifold(3, 'R3', r'\mathcal{M}')
    sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
    sage: p = M.point((1, pi/2, 0), name='P') # coordinates in the manifold's default chart
    sage: p
    Point P on the 3-dimensional topological manifold R3
    sage: latex(p)
    P

Computing the coordinates of the point in a new chart::

    sage: c_cart.<x,y,z> = M.chart()
    sage: ch = c_spher.transition_map(c_cart, [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
    sage: p.coord(c_cart) # evaluate P's Cartesian coordinates
    (1, 0, 0)

Points can be compared::

    sage: p1 = M.point((1, pi/2, 0))
    sage: p == p1
    True
    sage: q = M.point((1,2,3), c_cart, name='Q') # point defined by its Cartesian coordinates
    sage: p == q
    False

Listing all the coordinates of a point in different charts::

    sage: p._coordinates # random (dictionary output)
    {Chart (R3, (r, th, ph)): (1, 1/2*pi, 0), Chart (R3, (x, y, z)): (1, 0, 0)}
    sage: q._coordinates
    {Chart (R3, (x, y, z)): (1, 2, 3)}

"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Element

class TopManifoldPoint(Element):
    r"""
    Point of a topological manifold.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.manifold.TopManifold`.

    INPUT:

    - ``subset`` -- the manifold subset to which the point belongs (can be
      the entire manifold)
    - ``coords`` -- (default: ``None``) the point coordinates (as a tuple or a list)
    - ``chart`` -- (default: ``None``) chart in which the coordinates are given;
      if none is provided, the coordinates are assumed
      to refer to the subset's default chart
    - ``name`` -- (default: ``None``) name given to the point
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the point; if
      none is provided, the LaTeX symbol is set to ``name``
    - ``check_coords`` -- (default: ``True``) determines whether ``coords`` are
      valid coordinates for the chart ``chart``; for symbolic coordinates, it
      is recommended to set ``check_coords`` to ``False``.

    EXAMPLES:

    A point on a 2-dimensional manifold::

        sage: M = TopManifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: (a, b) = var('a b') # generic coordinates for the point
        sage: p = M.point((a, b), name='P'); p
        Point P on the 2-dimensional topological manifold M
        sage: p.coord()  # coordinates of P in the subset's default chart
        (a, b)

    Since points are Sage *elements*, the (facade) *parent* of which being the
    subset on which they are defined, it is equivalent to write::

        sage: p = M((a, b), name='P'); p
        Point P on the 2-dimensional topological manifold M

    A point is an element of the manifold subset in which it has been defined::

        sage: p in M
        True
        sage: p.parent()
        2-dimensional topological manifold M
        sage: U = M.open_subset('U', coord_def={c_xy: x>0})
        sage: q = U.point((2,1), name='q')
        sage: q in U
        True
        sage: q in M
        True

    Note that the parent of a point is always the manifold, not the subset
    in which it has been defined (the latter being returned by the method
    :meth:`containing_set`)::

        sage: q.parent()
        2-dimensional topological manifold M
        sage: q.containing_set()
        Open subset U of the 2-dimensional topological manifold M

    By default, the LaTeX symbol of the point is deduced from its name::

        sage: latex(p)
        P

    But it can be set to any value::

        sage: p = M.point((a, b), name='P', latex_name=r'\mathcal{P}')
        sage: latex(p)
        \mathcal{P}

    Points can be drawn in 2D or 3D graphics thanks to the method :meth:`plot`.

    """
    def __init__(self, subset, coords=None, chart=None, name=None,
                 latex_name=None, check_coords=True):
        r"""
        Construct a manifold point.

        TESTS::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,3), name='p'); p
            Point p on the 2-dimensional topological manifold M
            sage: TestSuite(p).run()
            sage: U = M.open_subset('U', coord_def={X: x<0})
            sage: q = U((-1,2), name='q'); q
            Point q on the 2-dimensional topological manifold M
            sage: TestSuite(q).run()

        """
        Element.__init__(self, subset._manifold)
        self._manifold = subset._manifold
        self._subset = subset
        self._coordinates = {}
        if coords is not None:
            if len(coords) != self._manifold._dim:
                raise ValueError("The number of coordinates must be equal" +
                                 " to the manifold dimension.")
            if chart is None:
                chart = self._subset._def_chart
            else:
                if chart not in self._subset._atlas:
                    raise ValueError("The " + str(chart) +
                           " has not been defined on the " + str(self._subset))
            if check_coords:
                if not chart.valid_coordinates(*coords):
                    raise ValueError("The coordinates " + str(coords) +
                                     " are not valid on the " + str(chart))
            for schart in chart._supercharts:
                self._coordinates[schart] = tuple(coords)
            for schart in chart._subcharts:
                if schart != chart:
                    if schart.valid_coordinates(*coords):
                        self._coordinates[schart] = tuple(coords)
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

    def _repr_(self):
        r"""
        Return a string representation of the point.
        """
        description = "Point"
        if self._name is not None:
            description += " " + self._name
        description += " on the {}".format(self._manifold)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of the point.
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def containing_set(self):
        r"""
        Return a manifold subset that contains ``self``.

        A priori, this method returns the manifold subset (possibly the
        manifold itself) in which the point ``self`` has been defined.

        OUTPUT:

        - an instance of
          :class:`~sage.manifolds.subset.TopManifoldSubset`

        EXAMPLES:

        Points on a 2-dimensional manifold::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,3), name='p'); p
            Point p on the 2-dimensional topological manifold M
            sage: p.containing_set()
            2-dimensional topological manifold M
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: q = U.point((2,1), name='q'); q
            Point q on the 2-dimensional topological manifold M
            sage: q.containing_set()
            Open subset U of the 2-dimensional topological manifold M

        Note that in the present case, the containing set is tighter than the
        parent, which is always the manifold::

            sage: q.parent()
            2-dimensional topological manifold M
            sage: q.containing_set().is_subset(q.parent())
            True
            sage: q.containing_set() != q.parent()
            True

        """
        return self._subset

    def coord(self, chart=None, old_chart=None):
        r"""
        Return the point coordinates in the specified chart.

        If these coordinates are not already known, they are computed from
        known ones by means of change-of-chart formulas.

        INPUT:

        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the subset's default chart
        - ``old_chart`` -- (default: ``None``) chart from which the coordinates in
          ``chart`` are to be computed. If ``None``, a chart in which the point's
          coordinates are already known will be picked, priveleging the
          subset's default chart.

        EXAMPLES:

        Spherical coordinates of a point on `\RR^3`::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(3, 'R3', r'\mathcal{M}')
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi') # spherical coordinates
            sage: p = M.point((1, pi/2, 0))
            sage: p.coord()    # coordinates on the manifold's default chart
            (1, 1/2*pi, 0)
            sage: p.coord(c_spher) # with the chart c_spher specified (same result as above since this is the default chart)
            (1, 1/2*pi, 0)

        Computing the Cartesian coordinates from the spherical ones::

            sage: c_cart.<x,y,z> = M.chart()  # Cartesian coordinates
            sage: c_spher.transition_map(c_cart, [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
            Change of coordinates from Chart (R3, (r, th, ph)) to Chart (R3, (x, y, z))
            sage: p.coord(c_cart)  # the computation is performed by means of the above change of coordinates
            (1, 0, 0)

        Coordinates of a point on a 2-dimensional manifold::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: (a, b) = var('a b') # generic coordinates for the point
            sage: p = M.point((a, b), name='P')
            sage: p.coord()  # coordinates of P in the manifold's default chart
            (a, b)

        Coordinates of P in a new chart::

            sage: c_uv.<u,v> = M.chart()
            sage: ch_xy_uv = c_xy.transition_map(c_uv, [x-y, x+y])
            sage: p.coord(c_uv)
            (a - b, a + b)

        Coordinates of P in a third chart::

            sage: c_wz.<w,z> = M.chart()
            sage: ch_uv_wz = c_uv.transition_map(c_wz, [u^3, v^3])
            sage: p.coord(c_wz, old_chart=c_uv)
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)

        Actually, in the present case, it is not necessary to specify
        old_chart='uv'::

            sage: p.set_coord((a-b, a+b), c_uv) # erases all the coordinates except those in the chart c_uv
            sage: p._coordinates
            {Chart (M, (u, v)): (a - b, a + b)}
            sage: p.coord(c_wz)
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)
            sage: p._coordinates  # random (dictionary output)
            {Chart (M, (u, v)): (a - b, a + b),
             Chart (M, (w, z)): (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)}

        """
        if chart is None:
            dom = self._subset
            chart = dom._def_chart
            def_chart = chart
        else:
            dom = chart._domain
            def_chart = dom._def_chart
            if self not in dom:
                raise ValueError("The point does not belong to the domain " +
                                 "of " + str(chart))
        if chart not in self._coordinates:
            # Check whether chart corresponds to a superchart of a chart
            # in which the coordinates are known:
            for ochart in self._coordinates:
                if chart in ochart._supercharts or chart in ochart._subcharts:
                    self._coordinates[chart] = self._coordinates[ochart]
                    return self._coordinates[chart]
            # If this point is reached, some change of coordinates must be
            # performed
            if old_chart is not None:
                s_old_chart = old_chart
                s_chart = chart
            else:
                # A chart must be found as a starting point of the computation
                # The domain's default chart is privileged:
                if def_chart in self._coordinates \
                        and (def_chart, chart) in dom._coord_changes:
                    old_chart = def_chart
                    s_old_chart = def_chart
                    s_chart = chart
                else:
                    for ochart in self._coordinates:
                        for subchart in ochart._subcharts:
                            if (subchart, chart) in dom._coord_changes:
                                old_chart = ochart
                                s_old_chart = subchart
                                s_chart = chart
                                break
                        if old_chart is not None:
                            break
                if old_chart is None:
                    # Some search involving the subcharts of chart is
                    # performed:
                    for schart in chart._subcharts:
                        for ochart in self._coordinates:
                            for subchart in ochart._subcharts:
                                if (subchart, schart) in dom._coord_changes:
                                    old_chart = ochart
                                    s_old_chart = subchart
                                    s_chart = schart
                                    break
                            if old_chart is not None:
                                break
                        if old_chart is not None:
                            break
            if old_chart is None:
                raise ValueError("The coordinates of " + str(self) + \
                    " in the " + str(chart) + " cannot be computed" + \
                    " by means of known changes of charts.")
            else:
                chcoord = dom._coord_changes[(s_old_chart, s_chart)]
                self._coordinates[chart] = \
                                    chcoord(*self._coordinates[old_chart])
        return self._coordinates[chart]

    def set_coord(self, coords, chart=None):
        r"""
        Sets the point coordinates in the specified chart.

        Coordinates with respect to other charts are deleted, in order to
        avoid any inconsistency. To keep them, use the method :meth:`add_coord`
        instead.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the subset's default chart

        EXAMPLES:

        Setting coordinates to a point on `\RR^3`::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart.<x,y,z> = M.chart()
            sage: p = M.point()
            sage: p.set_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)

        A point defined in another coordinate system::

            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: q = M.point()
            sage: q.set_coord((1,2,3), c_spher)
            sage: cart_from_spher = c_spher.transition_map(c_cart, [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])

        If we set the coordinates of q in the chart c_cart, those in the chart c_spher
        are lost::

            sage: q.set_coord( cart_from_spher(*q.coord(c_spher)), c_cart)
            sage: q._coordinates
            {Chart (R3, (x, y, z)): (cos(3)*sin(2), sin(3)*sin(2), cos(2))}
            sage: p._coordinates
            {Chart (R3, (x, y, z)): (1, 2, 3)}

        """
        self._coordinates.clear()
        self.add_coord(coords, chart)

    def add_coord(self, coords, chart=None):
        r"""
        Adds some coordinates in the specified chart.

        The previous coordinates with respect to other charts are kept. To
        clear them, use :meth:`set_coord` instead.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the subset's default chart

        .. WARNING::

           If the point has already coordinates in other charts, it
           is the user's responsability to make sure that the coordinates
           to be added are consistent with them.

        EXAMPLES:

        Setting coordinates to a point on `\RR^3`::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart.<x,y,z> = M.chart()
            sage: p = M.point()
            sage: p.add_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)

        A point defined in another coordinate system::

            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: q = M.point()
            sage: q.add_coord((1,2,3), c_spher)
            sage: cart_from_spher = c_spher.transition_map(c_cart, [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
            sage: q.add_coord( cart_from_spher(*q.coord(c_spher)), c_cart)
            sage: q._coordinates  # random (dictionary output)
            {Chart (R3, (r, th, ph)): (1, 2, 3),
             Chart (R3, (x, y, z)): (cos(3)*sin(2), sin(3)*sin(2), cos(2))}
            sage: p._coordinates
            {Chart (R3, (x, y, z)): (1, 2, 3)}
            sage: p == q  # p and q should differ because the coordinates (1,2,3) are on different charts
            False

        Contrary to :meth:`set_coord`, the method :meth:`add_coord` does not
        the coordinates in other charts::

            sage: p = M.point((1,2,3), c_spher)
            sage: p._coordinates
            {Chart (R3, (r, th, ph)): (1, 2, 3)}
            sage: p.set_coord((4,5,6), c_cart)
            sage: p._coordinates
            {Chart (R3, (x, y, z)): (4, 5, 6)}
            sage: p.add_coord((7,8,9), c_spher)
            sage: p._coordinates # random (dictionary output)
            {Chart (R3, (x, y, z)): (4, 5, 6), chart (R3, (r, th, ph)): (7, 8, 9)}

        """
        if len(coords) != self._manifold._dim:
            raise ValueError("The number of coordinates must be equal " +
                             "to the manifold dimension.")
        if chart is None:
            chart = self._subset._def_chart
        else:
            if chart not in self._subset._atlas:
                raise ValueError("The " + str(chart) +
                    " has not been defined on the " + str(self._subset))
        self._coordinates[chart] = coords

    def __eq__(self, other):
        r"""
        Compares the current point with another one.

        EXAMPLES:

        Comparison with coordinates in the same chart::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,-3), chart=X)
            sage: q = M((2,-3), chart=X)
            sage: p == q
            True
            sage: q = M((-2,-3), chart=X)
            sage: p == q
            False

        Comparison with coordinates of other in a subchart::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: XU = X.restrict(U)
            sage: q = U((2,-3), chart=XU)
            sage: p == q and q == p
            True
            sage: q = U((1,-3), chart=XU)
            sage: p == q or q == p
            False

        Comparison requiring a change of chart::

            sage: Y.<u,v> = U.chart()
            sage: XU_to_Y = XU.transition_map(Y, (ln(x), x+y))
            sage: XU_to_Y.inverse()(u,v)
            (e^u, v - e^u)
            sage: q = U((ln(2),-1), chart=Y)
            sage: p == q and q == p
            True
            sage: q = U((ln(3),1), chart=Y)
            sage: p == q or q == p
            False

        """
        if not isinstance(other, TopManifoldPoint):
            return False
        if other._manifold != self._manifold:
            return False
        # Search for a common chart to compare the coordinates
        common_chart = None
        # the subset's default chart is privileged:
        def_chart = self._subset._def_chart
        if def_chart in self._coordinates and def_chart in other._coordinates:
            common_chart = def_chart
        else:
            for chart in self._coordinates:
                if chart in other._coordinates:
                    common_chart = chart
                    break
        if common_chart is None:
            # At this stage, a commont chart is searched via a coordinate
            # transformation:
            for chart in self._coordinates:
                try:
                    other.coord(chart)
                    common_chart = chart
                    break
                except ValueError:
                    pass
            else:
                # Attempt a coordinate transformation in the reverse way:
                for chart in other._coordinates:
                    try:
                        self.coord(chart)
                        common_chart = chart
                        break
                    except ValueError:
                        pass
        if common_chart is None:
            return False
            #!# Another option would be:
            # raise ValueError("No common chart has been found to compare " +
            #                 str(self) + " and " + str(other))
        return self._coordinates[common_chart] == \
                                              other._coordinates[common_chart]

    def __ne__(self, other):
        r"""
        Non-equality operator.
        """
        return not self.__eq__(other)

    def __cmp__(self, other):
        r"""
        Old-style (Python 2) comparison operator.

        This is provisory, until migration to Python 3 is achieved.

        """
        if self.__eq__(other):
            return 0
        else:
            return -1

    def __hash__(self):
        r"""
        This hash function is set to constant on a given manifold, to fulfill
        Python's credo:
        p == q  ==>  hash(p) == hash(q)
        This is necessary since p and q may be created in different coordinate
        systems and nevertheless be equal
        """
        return self._manifold.__hash__()

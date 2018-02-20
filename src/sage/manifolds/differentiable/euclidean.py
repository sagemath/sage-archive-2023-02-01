r"""
Euclidean spaces


AUTHORS:

- Eric Gourgoulhon (2018): initial version

"""

#*****************************************************************************
#       Copyright (C) 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.functions.trig import cos, sin, atan2
from sage.functions.other import sqrt
from sage.manifolds.differentiable.pseudo_riemannian import \
                                                       PseudoRiemannianManifold

###############################################################################

class EuclideanSpaceGeneric(PseudoRiemannianManifold):
    r"""
    Euclidean space.

    """
    def __init__(self, n, name=None, latex_name=None,
                 coordinates='Cartesian', symbols=None, metric_name='g',
                 metric_latex_name=None, start_index=0, ambient=None,
                 category=None, names=None, init_coordinates=None):
        r"""
        Construct an Euclidean space.
        """
        if name is None:
            name = 'E^{}'.format(n)
            if latex_name is None:
                latex_name = r'\mathbb{E}^{' + str(n) + '}'
        PseudoRiemannianManifold.__init__(self, n, name, metric_name=metric_name,
                                          signature=n, ambient=ambient,
                                          latex_name=latex_name,
                                          metric_latex_name=metric_latex_name,
                                          start_index=start_index,
                                          category=category)
        if names is not None and symbols is None:
            symbols = ''
            for x in names:
                symbols += x + ' '
            symbols = symbols[:-1]
        if symbols is None:
            if n == 1:
                if coordinates == 'Cartesian':
                    symbols = 'x'
                else:
                    raise TypeError("unkown coordinate type")
            elif n > 3:
                if coordinates == 'Cartesian':
                    symbols = ''
                    for i in self.irange():
                        symbols += 'x{}:x_{{}} '.format(i,i)
                    symbols = symbols[:-1]
                else:
                    raise TypeError("unkown coordinate type")
            else:
                raise NotImplementedError("dimension not implemented yet")
        self._named_charts = {}
        self._named_frames = {}
        if init_coordinates is None:
            self._init_coordinates = {'Cartesian':
                                      self._init_coordinates_cartesian}
        else:
            self._init_coordinates = init_coordinates
        self._init_coordinates[coordinates](symbols)

    def _first_ngens(self, n):
        r"""
        Return the list of coordinates of the default chart.

        This is useful only for the use of Sage preparser::

            sage: preparse("E.<x,y,z> = EuclideanSpace(3)")
            "E = EuclideanSpace(Integer(3), names=('x', 'y', 'z',)); (x, y, z,) = E._first_ngens(3)"

        """
        return self._def_chart[:]

    def _init_coordinates_cartesian(self, symbols):
        r"""
        Construct the chart of Cartesian coordinates and initialize the
        components of the metric tensor in it.
        """
        chart = self.chart(coordinates=symbols)
        self._named_charts['Cartesian'] = chart
        frame = chart.frame()
        self._named_frames['Cartesian'] = frame
        g = self.metric()
        gc = g.add_comp(frame)
        for i in self.irange():
            gc[i, i, chart] = 1
        g.connection().coef(frame)

    def cartesian_coordinates(self, symbols=None):
        r"""
        Return Cartesian coordinates.
        """
        return self.cartesian_chart(symbols=symbols)[:]

###############################################################################

class EuclideanPlane(EuclideanSpaceGeneric):
    r"""
    Euclidean plane.
    """
    def __init__(self, name=None, latex_name=None, coordinates='Cartesian',
                 symbols=None, metric_name='g', metric_latex_name=None,
                 start_index=0, ambient=None, category=None, names=None):
        r"""
        Construct an Euclidean plane.
        """
        if names is not None and symbols is None:
            symbols = ''
            for x in names:
                symbols += x + ' '
            symbols = symbols[:-1]
            if coordinates == 'polar':
                if names[1] in ['p', 'ph', 'phi']:
                    symbols += ':\\phi'
                elif names[1] in ['t', 'th', 'theta']:
                    symbols += ':\\theta'
        if symbols is None:
            if coordinates == 'Cartesian':
                symbols = 'x y'
            elif coordinates == 'polar':
                symbols = 'r ph:\\phi'
            else:
                raise TypeError("unkown coordinate type")
        init_coordinates = {'Cartesian': self._init_coordinates_cartesian,
                            'polar': self._init_coordinates_polar}
        EuclideanSpaceGeneric.__init__(self, 2, name=name,
                                       latex_name=latex_name,
                                       coordinates=coordinates,
                                       symbols=symbols,
                                       metric_name=metric_name,
                                       metric_latex_name=metric_latex_name,
                                       start_index=start_index,
                                       ambient=ambient, category=category,
                                       init_coordinates=init_coordinates)

    def _init_coordinates_polar(self, symbols):
        r"""
        Construct the chart of polar coordinates and initialize the
        components of the metric tensor in it.
        """
        coords = symbols.split()  # list of strings, one per coordinate
        # Adding the coordinate ranges:
        coordinates = coords[0] + ':(0,+oo) ' + coords[1] + ':(0,2*pi)'
        chart = self.chart(coordinates=coordinates)
        self._named_charts['polar'] = chart
        frame = chart.frame()
        self._named_frames['polar_coord'] = frame
        # Initialization of the metric components and the associated
        # Christoffel symbols
        g = self.metric()
        gc = g.add_comp(frame)
        i0 = self._sindex
        i1 = i0 + 1
        r = chart[i0]
        gc[i0, i0, chart] = 1
        gc[i1, i1, chart] = r**2
        nabla = g.connection()
        nabla.coef(frame)
        # Orthonormal frame associated with polar coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i0, i0, chart] = 1
        to_orthonormal[frame, i1, i1, chart] = 1/r
        oframe = frame.new_frame(to_orthonormal, 'e')
        self._named_frames['polar_ortho'] = oframe
        g.comp(oframe)
        nabla.coef(oframe)

    def _transition_polar_cartesian(self):
        r"""
        Transitions between polar and Cartesian coordinates.
        """
        # Transition maps polar chart <-> Cartesian chart
        chart_cart = self._named_charts['Cartesian']
        chart_pol = self._named_charts['polar']
        x, y = chart_cart[:]
        r, ph = chart_pol[:]
        pol_to_cart = chart_pol.transition_map(chart_cart,
                                               [r*cos(ph), r*sin(ph)])
        pol_to_cart.set_inverse(sqrt(x**2+y**2), atan2(y,x))
        # Automorphisms ortho polar frame <-> Cartesian frame
        oframe = self._named_frames['polar_ortho']
        cframe = chart_cart.frame()
        sframe = chart_pol.frame()
        changes = self._frame_changes
        cframe_to_oframe = changes[(sframe, oframe)] * \
                           changes[(cframe, sframe)]
        oframe_to_cframe = changes[(sframe, cframe)] * \
                           changes[(oframe, sframe)]
        changes[(cframe, oframe)] = cframe_to_oframe
        changes[(oframe, cframe)] = oframe_to_cframe
        vmodule = self.vector_field_module()
        vmodule._basis_changes[(cframe, oframe)] = cframe_to_oframe
        vmodule._basis_changes[(oframe, cframe)] = oframe_to_cframe

    def polar_chart(self, symbols=None):
        r"""
        Return the chart of polar coordinates.
        """
        if 'polar' not in self._named_charts:
            if symbols is None:
                symbols = 'r ph:\\phi'
            self._init_coordinates_polar(symbols)
            if 'Cartesian' in self._named_charts:
                self._transition_polar_cartesian()
        return self._named_charts['polar']

    def polar_coordinates(self, symbols=None):
        r"""
        Return polar coordinates.
        """
        return self.polar_chart(symbols=symbols)[:]

    def cartesian_chart(self, symbols=None):
        r"""
        Return the chart of Cartesian coordinates.
        """
        if 'Cartesian' not in self._named_charts:
            if symbols is None:
                symbols = 'x y'
            self._init_coordinates_cartesian(symbols)
            if 'polar' in self._named_charts:
                self._transition_polar_cartesian()
        return self._named_charts['Cartesian']


###############################################################################

def EuclideanSpace(n, name=None, latex_name=None, coordinates='Cartesian',
                   symbols=None, metric_name='g', metric_latex_name=None,
                   start_index=1, names=None):
    r"""
    Construct an Euclidean space.
    """
    if n == 2:
        return EuclideanPlane(name=name, latex_name=latex_name,
                              coordinates=coordinates, symbols=symbols,
                              metric_name=metric_name,
                              metric_latex_name=metric_latex_name,
                              start_index=start_index, names=names)
    return EuclideanSpaceGeneric(n, name=name, latex_name=latex_name,
                              coordinates=coordinates, symbols=symbols,
                              metric_name=metric_name,
                              metric_latex_name=metric_latex_name,
                              start_index=start_index, names=names)


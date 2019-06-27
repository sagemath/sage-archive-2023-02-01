# -*- coding: utf-8 -*-
"""
The Normaliz backend for polyhedral computations

.. NOTE::

    This backend requires `PyNormaliz <https://pypi.python.org/pypi/PyNormaliz/1.5>`_.
    To install PyNormaliz, type :code:`sage -i pynormaliz` in the terminal.

AUTHORS:

- Matthias Köppe (2016-12): initial version
- Jean-Philippe Labbé (2019-04): Expose normaliz features and added functionalities
"""

#*****************************************************************************
#  Copyright (C) 2016 Matthias Köppe <mkoeppe at math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function

from sage.structure.element import Element
from sage.misc.all import prod
from sage.features import PythonModule

from sage.rings.all import ZZ, QQ
from sage.arith.functions import LCM_list
from sage.misc.functional import denominator
from sage.matrix.constructor import vector

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ


def _format_function_call(fn_name, *v, **k):
    """
    Return a Python function call as a string.

    Keywords are sorted.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.backend_normaliz import _format_function_call
        sage: _format_function_call('foo', 17, hellooooo='goodbyeeee')
        "foo(17, hellooooo='goodbyeeee')"
    """
    args = [ repr(a) for a in v ] + [ "%s=%r" % (arg, val) for arg, val in sorted(k.items()) ]
    return "{}({})".format(fn_name, ", ".join(args))

#########################################################################
class Polyhedron_normaliz(Polyhedron_base):
    """
    Polyhedra with normaliz

    INPUT:

    - ``parent`` -- :class:`~sage.geometry.polyhedron.parent.Polyhedra`
      the parent

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``; the
      V-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the H-representation

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``; the
      H-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the V-representation

    - ``normaliz_cone`` -- a PyNormaliz wrapper of a normaliz cone

    Only one of ``Vrep``, ``Hrep``, or ``normaliz_cone`` can be different
    from ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)],   # optional - pynormaliz
        ....:                lines=[], backend='normaliz')
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - pynormaliz

    Two ways to get the full space::

        sage: Polyhedron(eqns=[[0, 0, 0]], backend='normaliz')             # optional - pynormaliz
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines
        sage: Polyhedron(ieqs=[[0, 0, 0]], backend='normaliz')             # optional - pynormaliz
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines

    A lower-dimensional affine cone; we test that there are no mysterious
    inequalities coming in from the homogenization::

        sage: P = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],             # optional - pynormaliz
        ....:                backend='normaliz')
        sage: P.n_inequalities()                                           # optional - pynormaliz
        1
        sage: P.equations()                                                # optional - pynormaliz
        (An equation (1, 0) x - 1 == 0,)

    The empty polyhedron::

        sage: P=Polyhedron(ieqs=[[-2, 1, 1], [-3, -1, -1], [-4, 1, -2]],   # optional - pynormaliz
        ....:              backend='normaliz')
        sage: P                                                            # optional - pynormaliz
        The empty polyhedron in QQ^2
        sage: P.Vrepresentation()                                          # optional - pynormaliz
        ()
        sage: P.Hrepresentation()                                          # optional - pynormaliz
        (An equation -1 == 0,)

    TESTS:

    Tests copied from various methods in :mod:`sage.geometry.polyhedron.base`::

        sage: p = Polyhedron(vertices = [[1,0,0], [0,1,0], [0,0,1]],       # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_equations()                                              # optional - pynormaliz
        1
        sage: p.n_inequalities()                                           # optional - pynormaliz
        3

        sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)],   # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_facets()                                                 # optional - pynormaliz
        8

        sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]], # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_vertices()                                               # optional - pynormaliz
        2

        sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]],       # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_rays()                                                   # optional - pynormaliz
        1

        sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]],      # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_lines()                                                  # optional - pynormaliz
        1

    """
    def __init__(self, parent, Vrep, Hrep, normaliz_cone=None, normaliz_data=None, **kwds):
        """
        Initializes the polyhedron.

        See :class:`Polyhedron_normaliz` for a description of the input
        data.

        TESTS:

        We skip the pickling test because pickling is currently
        not implemented::

            sage: p = Polyhedron(backend='normaliz')                 # optional - pynormaliz
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - pynormaliz
            sage: p = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],   # optional - pynormaliz
            ....:                backend='normaliz')
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - pynormaliz
            sage: p = Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],  # optional - pynormaliz
            ....:                backend='normaliz')
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - pynormaliz
        """
        if normaliz_cone:
            if Hrep is not None or Vrep is not None or normaliz_data is not None:
                raise ValueError("only one of Vrep, Hrep, normaliz_cone, or normaliz_data can be different from None")
            Element.__init__(self, parent=parent)
            self._init_from_normaliz_cone(normaliz_cone)
        elif normaliz_data:
            if Hrep is not None or Vrep is not None:
                raise ValueError("only one of Vrep, Hrep, normaliz_cone, or normaliz_data can be different from None")
            Element.__init__(self, parent=parent)
            self._init_from_normaliz_data(normaliz_data)
        else:
            Polyhedron_base.__init__(self, parent, Vrep, Hrep, **kwds)

    def _nmz_result(self, normaliz_cone, property):
        """
        Call PyNormaliz's NmzResult function.

        TESTS::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)],   # optional - pynormaliz
            ....:                lines=[], backend='normaliz')
            sage: p._nmz_result(p._normaliz_cone, 'EquivariantXyzzyModuleSeries')  # optional - pynormaliz
            Traceback (most recent call last):
            ...
            NormalizError: Some error in the normaliz input data detected: Unknown ConeProperty...
        """
        PythonModule("PyNormaliz", spkg="pynormaliz").require()
        import PyNormaliz
        return PyNormaliz.NmzResult(normaliz_cone, property)

    def _init_from_normaliz_cone(self, normaliz_cone):
        """
        Construct polyhedron from a PyNormaliz wrapper of a normaliz cone.

        TESTS::

            sage: p = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz   # optional - pynormaliz
            sage: Polyhedron_normaliz._init_from_Hrepresentation(p, [], [])  # indirect doctest  # optional - pynormaliz
        """
        if normaliz_cone and self._nmz_result(normaliz_cone, "AffineDim") < 0:
            # Empty polyhedron. Special case because Normaliz defines the
            # recession cone of an empty polyhedron given by an
            # H-representation as the cone defined by the homogenized system.
            self._init_empty_polyhedron()
        else:
            self._normaliz_cone = normaliz_cone
            self._init_Vrepresentation_from_normaliz()
            self._init_Hrepresentation_from_normaliz()

    def _init_from_normaliz_data(self, data, verbose=False):
        """
        Construct polyhedron from normaliz ``data`` (a dictionary).

        TESTS::

            sage: p = Polyhedron(backend='normaliz', ambient_dim=2)                             # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_QQ_normaliz  # optional - pynormaliz
            sage: data = {'inhom_inequalities': [[-1L, 2L, 0L], [0L, 0L, 1L], [2L, -1L, 0L]]}   # optional - pynormaliz
            sage: Polyhedron_QQ_normaliz._init_from_normaliz_data(p, data)                      # optional - pynormaliz
            sage: p.inequalities_list()                                                         # optional - pynormaliz
            [[0, -1, 2], [0, 2, -1]]

        """
        if verbose:
            import six
            if isinstance(verbose, six.string_types):
                print("# Wrote equivalent Normaliz input file to {}".format(verbose))
                self._normaliz_format(data, file_output=verbose)
            else:
                print("# ----8<---- Equivalent Normaliz input file ----8<----")
                print(self._normaliz_format(data), end='')
                print("# ----8<-------------------8<-------------------8<----")

        if verbose:
            print("# Calling {}".format(_format_function_call('PyNormaliz.NmzCone', **data)))
        PythonModule("PyNormaliz", spkg="pynormaliz").require()
        import PyNormaliz
        cone = PyNormaliz.NmzCone(**data)
        assert cone, "{} did not return a cone".format(_format_function_call('PyNormaliz.NmzCone', **data))

        self._init_from_normaliz_cone(cone)

    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)
            sage: p._is_zero(0)
            True
            sage: p._is_zero(1/100000)
            False
        """
        return x == 0

    def _is_nonneg(self, x):
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)
            sage: p._is_nonneg(1)
            True
            sage: p._is_nonneg(-1/100000)
            False
        """
        return x >= 0

    def _is_positive(self, x):
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)
            sage: p._is_positive(1)
            True
            sage: p._is_positive(0)
            False
        """
        return x > 0

    def _init_from_Vrepresentation(self, vertices, rays, lines, minimize=True, verbose=False):
        r"""
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of point; each point can be specified
           as any iterable container of
           :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``rays`` -- list of rays; each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``lines`` -- list of lines; each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz   # optional - pynormaliz
            sage: Polyhedron_normaliz._init_from_Vrepresentation(p, [], [], [])   # optional - pynormaliz
        """

        def vert_ray_line_QQ():
            nmz_vertices = []
            for v in vertices:
                d = LCM_list([denominator(v_i) for v_i in v])
                dv = [ d*v_i for v_i in v ]
                nmz_vertices.append(dv + [d])
            nmz_rays = []
            for r in rays:
                d = LCM_list([denominator(r_i) for r_i in r])
                dr = [ d*r_i for r_i in r ]
                nmz_rays.append(dr)
            nmz_lines = []
            for l in lines:
                d = LCM_list([denominator(l_i) for l_i in l])
                dl = [ d*l_i for l_i in l ]
                nmz_lines.append(dl)
            return nmz_vertices, nmz_rays, nmz_lines

        if vertices is None:
                vertices = []
        if rays is None:
                rays = []
        if lines is None:
                lines = []

        nmz_vertices, nmz_rays, nmz_lines = vert_ray_line_QQ()

        if not nmz_vertices and not nmz_rays and not nmz_lines:
            # Special case to avoid:
            #   error: Some error in the normaliz input data detected:
            #   All input matrices empty!
            self._init_empty_polyhedron()
        else:
            data = {"vertices": nmz_vertices,
                    "cone": nmz_rays,
                    "subspace": nmz_lines}

            self._init_from_normaliz_data(data, verbose=verbose)

    def _init_from_Hrepresentation(self, ieqs, eqns, minimize=True, verbose=False):
        r"""
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``eqns`` -- list of equalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``minimize`` -- boolean (default: ``True``); ignored

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz   # optional - pynormaliz
            sage: Polyhedron_normaliz._init_from_Hrepresentation(p, [], [])   # optional - pynormaliz
        """

        def nmz_ieqs_eqns_QQ():
            nmz_ieqs = []
            for ieq in ieqs:
                d = LCM_list([denominator(ieq_i) for ieq_i in ieq])
                dieq = [ ZZ(d*ieq_i) for ieq_i in ieq ]
                b = dieq[0]
                A = dieq[1:]
                nmz_ieqs.append(A + [b])
            nmz_eqns = []
            for eqn in eqns:
                d = LCM_list([denominator(eqn_i) for eqn_i in eqn])
                deqn = [ ZZ(d*eqn_i) for eqn_i in eqn ]
                b = deqn[0]
                A = deqn[1:]
                nmz_eqns.append(A + [b])
            return nmz_ieqs, nmz_eqns

        if ieqs is None:
            ieqs = []
        if eqns is None:
            eqns = []

        nmz_ieqs, nmz_eqns = nmz_ieqs_eqns_QQ()

        if not nmz_ieqs:
            # If normaliz gets an empty list of inequalities, it adds
            # nonnegativities. So let's add a tautological inequality to work
            # around this.
            nmz_ieqs.append([0]*self.ambient_dim() + [0])
        data = {"inhom_equations": nmz_eqns,
                "inhom_inequalities": nmz_ieqs}
        self._init_from_normaliz_data(data, verbose=verbose)

    def _init_Vrepresentation_from_normaliz(self):
        r"""
        Create the Vrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],  # indirect doctest # optional - pynormaliz
            ....:                backend='normaliz')
            sage: p.Hrepresentation()                               # optional - pynormaliz
            (An inequality (-5, 12) x + 10 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (1, 4) x - 2 >= 0)
            sage: p.Vrepresentation()                               # optional - pynormaliz
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
        """
        self._Vrepresentation = []
        parent = self.parent()
        base_ring = self.base_ring()
        cone = self._normaliz_cone
        for g in self._nmz_result(cone, "VerticesOfPolyhedron"):
            d = g[-1]
            if d == 1:
                parent._make_Vertex(self, g[:-1])
            else:
                parent._make_Vertex(self, [base_ring(x)/d for x in g[:-1]])
        for g in self._nmz_result(cone, "ExtremeRays"):
            parent._make_Ray(self, g[:-1])
        for g in self._nmz_result(cone, "MaximalSubspace"):
            parent._make_Line(self, g[:-1])
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_from_normaliz(self):
        r"""
        Create the Hrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest # optional - pynormaliz
            ....:                backend='normaliz')
            sage: p.Hrepresentation()                                 # optional - pynormaliz
            (An inequality (-5, 12) x + 10 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (1, 4) x - 2 >= 0)
            sage: p.Vrepresentation()                                 # optional - pynormaliz
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
        """
        self._Hrepresentation = []
        cone = self._normaliz_cone
        parent = self.parent()
        for g in self._nmz_result(cone, "SupportHyperplanes"):
            if all(x == 0 for x in g[:-1]):
                # Ignore vertical inequality
                pass
            else:
                parent._make_Inequality(self, (g[-1],) + tuple(g[:-1]))
        for g in self._nmz_result(cone, "Equations"):
            parent._make_Equation(self, (g[-1],) + tuple(g[:-1]))
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_empty_polyhedron(self):
        r"""
        Initializes an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='normaliz'); empty            # optional - pynormaliz
            The empty polyhedron in ZZ^0
            sage: empty.Vrepresentation()                                  # optional - pynormaliz
            ()
            sage: empty.Hrepresentation()                                  # optional - pynormaliz
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [], backend='normaliz')            # optional - pynormaliz
            The empty polyhedron in ZZ^0
            sage: Polyhedron(backend='normaliz')._init_empty_polyhedron()  # optional - pynormaliz
        """
        super(Polyhedron_normaliz, self)._init_empty_polyhedron()
        # Can't seem to set up an empty _normaliz_cone.
        # For example, PyNormaliz.NmzCone(vertices=[]) gives
        # error: Some error in the normaliz input data detected: All input matrices empty!
        self._normaliz_cone = None

    @classmethod
    def _from_normaliz_cone(cls, parent, normaliz_cone):
        r"""
        Initializes a polyhedron from a PyNormaliz wrapper of a normaliz cone.

        TESTS::

            sage: P=Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]],   # optional - pynormaliz
            ....:              backend='normaliz')
            sage: PI = P.integral_hull()                 # indirect doctest; optional - pynormaliz
        """
        return cls(parent, None, None, normaliz_cone=normaliz_cone)

    @staticmethod
    def _make_normaliz_cone(data, verbose=False):
        r"""
        Returns a normaliz cone from ``data``.

        INPUT:

        - ``data`` -- a dictionary

        - ``verbose`` -- a boolean (default: ``False``)

        TESTS::

            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz    # optional - pynormaliz
            sage: data = {'inhom_inequalities': [[-1L, 2L, 0L], [0L, 0L, 1L], [2L, -1L, 0L]]}  # optional - pynormaliz
            sage: nmz_cone = Polyhedron_normaliz._make_normaliz_cone(data,verbose=False)       # optional - pynormaliz
            sage: from PyNormaliz import NmzResult                                             # optional - pynormaliz
            sage: NmzResult(nmz_cone, "ExtremeRays")                                           # optional - pynormaliz
            [[1L, 2L, 0L], [2L, 1L, 0L]]
        """
        PythonModule("PyNormaliz", spkg="pynormaliz").require()
        import PyNormaliz
        if verbose:
            print("# Calling PyNormaliz.NmzCone(**{})".format(data))
        cone = PyNormaliz.NmzCone(**data)
        assert cone, "NmzCone(**{}) did not return a cone".format(data)
        return cone

    @staticmethod
    def _cone_generators(pynormaliz_cone):
        r"""
        Returns the generators of a pynormaliz cone.

        This is particularly useful to get the reordering of the vertices (or
        rays) that is internally used by normaliz.

        INPUT:

        - ``pynormaliz_cone`` -- a pynormaliz cone object.

        OUTPUT:

        - a tuple of generators for the cone.

        TESTS::

            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz     # optional - pynormaliz
            sage: data = {'inhom_inequalities': [[-1L, 2L, 0L], [0L, 0L, 1L], [2L, -1L, 0L]]}   # optional - pynormaliz
            sage: nmz_cone = Polyhedron_normaliz._make_normaliz_cone(data,verbose=False)        # optional - pynormaliz
            sage: Polyhedron_normaliz._cone_generators(nmz_cone)                                # optional - pynormaliz
            [[1L, 2L, 0L], [0L, 0L, 1L], [2L, 1L, 0L]]
        """
        PythonModule("PyNormaliz", spkg="pynormaliz").require()
        import PyNormaliz
        return PyNormaliz.NmzResult(pynormaliz_cone, "Generators")

    def _get_nmzcone_data(self):
        r"""
        Get the data necessary to reproduce the normaliz cone.

        OUTPUT:

        - ``data`` -- a dictionary.

        TESTS:

        The empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')                               # optional - pynormaliz
            sage: P._get_nmzcone_data()                                            # optional - pynormaliz
            {}

        Another simple example::

            sage: C = Polyhedron(backend='normaliz',rays=[[1,2],[2,1]])            # optional - pynormaliz
            sage: C._get_nmzcone_data()                                            # optional - pynormaliz
            {'cone': [[1L, 2L], [2L, 1L]],
             'inhom_equations': [],
             'inhom_inequalities': [[-1L, 2L, 0L], [0L, 0L, 1L], [2L, -1L, 0L]],
             'subspace': [],
             'vertices': [[0L, 0L, 1L]]}
        """
        if self.is_empty():
            return {}

        vertices = self._nmz_result(self._normaliz_cone, "VerticesOfPolyhedron")
        # get rid of the last 0 in rays:
        rays = [r[:-1] for r in self._nmz_result(self._normaliz_cone, "ExtremeRays")]
        lines = self._nmz_result(self._normaliz_cone, "MaximalSubspace")
        ineqs = self._nmz_result(self._normaliz_cone, "SupportHyperplanes")
        eqs = self._nmz_result(self._normaliz_cone, "Equations")

        data = {'vertices': vertices,
                'cone': rays,
                'subspace': lines,
                'inhom_equations': eqs,
                'inhom_inequalities': ineqs}

        return data

    def _normaliz_format(self, data, file_output=None):
        r"""
        Return a string containing normaliz format.

        INPUT:

        - ``data`` -- a dictionary of PyNormaliz cone input properties

        - ``file_output`` (string; optional) -- a filename to which the
          representation should be written. If set to ``None`` (default),
          representation is returned as a string.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0]], # indirect doctest; optional - pynormaliz
            ....:                backend='normaliz', verbose=True)
            # ----8<---- Equivalent Normaliz input file ----8<----
            amb_space 2
            cone 0
            subspace 0
            vertices 3
             0 0 1
             0 1 1
             1 0 1
            # ----8<-------------------8<-------------------8<----
            # Calling ...
        """
        def format_number(x):
            try:
                return '{}'.format(QQ(x))
            except (ValueError, TypeError):
                return '({})'.format(x.polynomial('a'))
        def format_field(key, value):
            if isinstance(value, list) or isinstance(value, tuple):
                s = '{} {}\n'.format(key, len(value))
                for e in value:
                    for x in e:
                        s += ' ' + format_number(x)
                    s += '\n'
                return s
            else:
                return '{} {}\n'.format(key, value)

        s = format_field('amb_space', self.ambient_dim())
        for key, value in sorted(data.items()):
            s += format_field(key, value)
        if file_output is not None:
            in_file = open(file_output, 'w')
            in_file.write(s)
            in_file.close()
        else:
            return s

    def integral_hull(self):
        r"""
        Return the integral hull in the polyhedron.

        This is a new polyhedron that is the convex hull of all integral
        points.

        EXAMPLES:

        Unbounded example from Normaliz manual, "a dull polyhedron"::

            sage: P = Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]], # optional - pynormaliz
            ....:              backend='normaliz')
            sage: PI = P.integral_hull()                                   # optional - pynormaliz
            sage: P.plot(color='yellow') + PI.plot(color='green')          # optional - pynormaliz
            Graphics object consisting of 10 graphics primitives
            sage: PI.Vrepresentation()                                     # optional - pynormaliz
            (A vertex at (-1, 0), A vertex at (0, 1), A ray in the direction (1, 0))

        Nonpointed case::

            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]],     # optional - pynormaliz
            ....:              lines=[[-1, 1]], backend='normaliz')
            sage: PI = P.integral_hull()                                   # optional - pynormaliz
            sage: PI.Vrepresentation()                                     # optional - pynormaliz
            (A vertex at (1, 0),
             A ray in the direction (1, 0),
             A line in the direction (1, -1))

        Empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: PI = P.integral_hull()                                   # optional - pynormaliz
            sage: PI.Vrepresentation()                                     # optional - pynormaliz
            ()
        """
        if self.is_empty():
            return self
        cone = self._nmz_result(self._normaliz_cone, "IntegerHull")
        return self.parent().element_class._from_normaliz_cone(parent=self.parent(),
                                                               normaliz_cone=cone)

    def ehrhart_series(self, variable='t'):
        r"""
        Return the Ehrhart series of a compact rational polyhedron.

        The Ehrhart series is the generating function where the coefficient of
        `t^k` is number of integer lattice points inside the `k`-th dilation of
        the polytope.

        INPUT:

        - ``variable`` -- string (default: ``'t'``)

        OUTPUT:

        A rational function.

        EXAMPLES::

            sage: S = Polyhedron(vertices=[[0,1],[1,0]], backend='normaliz')  # optional - pynormaliz
            sage: ES = S.ehrhart_series()                                     # optional - pynormaliz
            sage: ES.numerator()                                              # optional - pynormaliz
            1
            sage: ES.denominator().factor()                                   # optional - pynormaliz
            (t - 1)^2

            sage: C = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],backend='normaliz') # optional - pynormaliz
            sage: ES = C.ehrhart_series()            # optional - pynormaliz
            sage: ES.numerator()                     # optional - pynormaliz
            t^2 + 4*t + 1
            sage: ES.denominator().factor()          # optional - pynormaliz
            (t - 1)^4

        The following example is from the Normaliz manual contained in the file
        ``rational.in``::

            sage: rat_poly = Polyhedron(vertices=[[1/2,1/2],[-1/3,-1/3],[1/4,-1/2]],backend='normaliz') # optional - pynormaliz
            sage: ES = rat_poly.ehrhart_series()                                       # optional - pynormaliz
            sage: ES.numerator()                                                       # optional - pynormaliz
            2*t^6 + 3*t^5 + 4*t^4 + 3*t^3 + t^2 + t + 1
            sage: ES.denominator().factor()                                            # optional - pynormaliz
            (-1) * (t + 1)^2 * (t - 1)^3 * (t^2 + 1) * (t^2 + t + 1)

        The polyhedron should be compact::

            sage: C = Polyhedron(backend='normaliz',rays=[[1,2],[2,1]])  # optional - pynormaliz
            sage: C.ehrhart_series()                                     # optional - pynormaliz
            Traceback (most recent call last):
            ...
            NotImplementedError: Ehrhart series can only be computed for compact polyhedron

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.backend_normaliz.hilbert_series`
        """
        if self.is_empty():
            return 0

        if not self.is_compact():
            raise NotImplementedError("Ehrhart series can only be computed for compact polyhedron")

        cone = self._normaliz_cone
        e = self._nmz_result(cone, "EhrhartSeries")
        # The output format of PyNormaliz is a list with 3 things:
        # 1) the coefficients of the h^*-polynomial
        # 2) a list of the exponents e such that (1-t^e) appears as a factor in
        # the denominator
        # 3) a shifting of the generating function.

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.fraction_field import FractionField
        poly_ring = FractionField(PolynomialRing(ZZ, variable))
        t = poly_ring.gens()[0]
        es = sum([e[0][i]*t**i for i in range(len(e[0]))])
        for expo in range(len(e[1])):
            es = es / (1 - t**e[1][expo])

        # The shift:
        es = es * t**e[2]

        return es

    def ehrhart_quasipolynomial(self, variable='t'):
        r"""
        Return the Ehrhart quasi-polynomial of a compact rational polyhedron
        using Normaliz.

        INPUT:

        - ``variable`` -- string (default: ``'t'``)

        OUTPUT:

        If it is a polynomial, returns the polynomial. Otherwise, returns a
        tuple of rational polynomials whose length is the quasi-period of the
        quasi-polynomial and each rational polynomial describes a residue class.

        EXAMPLES::

            sage: C = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],backend='normaliz') # optional - pynormaliz
            sage: C.ehrhart_quasipolynomial()  # optional - pynormaliz
            t^3 + 3*t^2 + 3*t + 1

            sage: P = Polyhedron(vertices=[[0,0],[3/2,0],[0,3/2],[1,1]],backend='normaliz')  # optional - pynormaliz
            sage: P.ehrhart_quasipolynomial()  # optional - pynormaliz
            (3/2*t^2 + 2*t + 1, 3/2*t^2 + 2*t + 1/2)
            sage: P.ehrhart_quasipolynomial('x')  # optional - pynormaliz
            (3/2*x^2 + 2*x + 1, 3/2*x^2 + 2*x + 1/2)

        The quasi-polynomial evaluated at ``i`` counts the integral points
        in the ``i``-th dilate::

            sage: Q = Polyhedron(vertices = [[-1/3],[2/3]],backend='normaliz')  # optional - pynormaliz
            sage: p0,p1,p2 = Q.ehrhart_quasipolynomial()  # optional - pynormaliz
            sage: r0 = [p0(i) for i in range(15)]         # optional - pynormaliz
            sage: r1 = [p1(i) for i in range(15)]         # optional - pynormaliz
            sage: r2 = [p2(i) for i in range(15)]         # optional - pynormaliz
            sage: result = [None]*15                      # optional - pynormaliz
            sage: result[::3] = r0[::3]                   # optional - pynormaliz
            sage: result[1::3] = r1[1::3]                 # optional - pynormaliz
            sage: result[2::3] = r2[2::3]                 # optional - pynormaliz
            sage: result == [(i*Q).integral_points_count() for i in range(15)]  # optional - pynormaliz
            True

        The polyhedron should be compact::

            sage: C = Polyhedron(backend='normaliz',rays=[[1,2],[2,1]])  # optional - pynormaliz
            sage: C.ehrhart_quasipolynomial()                            # optional - pynormaliz
            Traceback (most recent call last):
            ...
            NotImplementedError: Ehrhart quasi-polynomial can only be computed for compact polyhedron

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.backend_normaliz.hilbert_series`,
            :meth:`~sage.geometry.polyhedron.backend_normaliz.ehrhart_series`
        """
        if self.is_empty():
            return 0

        if not self.is_compact():
            raise NotImplementedError("Ehrhart quasi-polynomial can only be computed for compact polyhedron")

        cone = self._normaliz_cone
        # Normaliz needs to compute the EhrhartSeries first
        PythonModule("PyNormaliz", spkg="pynormaliz").require()
        import PyNormaliz
        assert PyNormaliz.NmzCompute(cone, ["EhrhartSeries"])
        e = self._nmz_result(cone, "EhrhartQuasiPolynomial")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        poly_ring = PolynomialRing(QQ, variable)
        t = poly_ring.gens()[0]
        if len(e) == 2:
            # It is a polynomial
            es = sum([e[0][i]*t**i for i in range(len(e[0]))])
            return es / ZZ(e[1])
        else:
            # It is a quasi-polynomial
            polynomials = []
            for p in e[:-1]:
                es = sum([p[i]*t**i for i in range(len(p))]) / ZZ(e[-1])
                polynomials += [es]

        return tuple(polynomials)

    def hilbert_series(self, grading, variable='t'):
        r"""
        Return the Hilbert series of the polyhedron with respect to ``grading``.

        INPUT:

        - ``grading`` -- vector. The grading to use to form the Hilbert series

        - ``variable`` -- string (default: ``'t'``)

        OUTPUT:

        A rational function.

        EXAMPLES::

            sage: C = Polyhedron(backend='normaliz',rays=[[0,0,1],[0,1,1],[1,0,1],[1,1,1]]) # optional - pynormaliz
            sage: HS = C.hilbert_series([1,1,1]) # optional - pynormaliz
            sage: HS.numerator() # optional - pynormaliz
            t^2 + 1
            sage: HS.denominator().factor() # optional - pynormaliz
            (-1) * (t + 1) * (t - 1)^3 * (t^2 + t + 1)

        By changing the grading, you can get the Ehrhart series of the square
        lifted at height 1::

            sage: C.hilbert_series([0,0,1]) # optional - pynormaliz
            (t + 1)/(-t^3 + 3*t^2 - 3*t + 1)

        Here is an example ``2cone.in`` from the Normaliz manual::

            sage: C = Polyhedron(backend='normaliz',rays=[[1,3],[2,1]]) # optional - pynormaliz
            sage: HS = C.hilbert_series([1,1]) # optional - pynormaliz
            sage: HS.numerator() # optional - pynormaliz
            t^5 + t^4 + t^3 + t^2 + 1
            sage: HS.denominator().factor() # optional - pynormaliz
            (t + 1) * (t - 1)^2 * (t^2 + 1) * (t^2 + t + 1)

            sage: HS = C.hilbert_series([1,2]) # optional - pynormaliz
            sage: HS.numerator() # optional - pynormaliz
            t^8 + t^6 + t^5 + t^3 + 1
            sage: HS.denominator().factor() # optional - pynormaliz
            (t + 1) * (t - 1)^2 * (t^2 + 1) * (t^6 + t^5 + t^4 + t^3 + t^2 + t + 1)

        Here is the magic square example form the Normaliz manual::

            sage: eq = [[0,1,1,1,-1,-1,-1, 0, 0, 0],
            ....:       [0,1,1,1, 0, 0, 0,-1,-1,-1],
            ....:       [0,0,1,1,-1, 0, 0,-1, 0, 0],
            ....:       [0,1,0,1, 0,-1, 0, 0,-1, 0],
            ....:       [0,1,1,0, 0, 0,-1, 0, 0,-1],
            ....:       [0,0,1,1, 0,-1, 0, 0, 0,-1],
            ....:       [0,1,1,0, 0,-1, 0,-1, 0, 0]]
            sage: magic_square = Polyhedron(eqns=eq,backend='normaliz') & Polyhedron(rays=identity_matrix(9).rows()) # optional - pynormaliz
            sage: grading = [1,1,1,0,0,0,0,0,0]
            sage: magic_square.hilbert_series(grading) # optional - pynormaliz
            (t^6 + 2*t^3 + 1)/(-t^9 + 3*t^6 - 3*t^3 + 1)

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.backend_normaliz.ehrhart_series`
        """
        if self.is_empty():
            return 0

        data = self._get_nmzcone_data()
        data['grading'] = [grading]
        new_cone = self._make_normaliz_cone(data)
        h = self._nmz_result(new_cone, "HilbertSeries")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.fraction_field import FractionField
        poly_ring = FractionField(PolynomialRing(ZZ, variable))
        t = poly_ring.gens()[0]
        hs = sum([h[0][i]*t**i for i in range(len(h[0]))])
        for expo in range(len(h[1])):
            hs = hs / (1 - t**h[1][expo])

        # The shift:
        hs = hs * t**h[2]

        return hs

    def integral_points(self, threshold=10000):
        r"""
        Return the integral points in the polyhedron.

        Uses either the naive algorithm (iterate over a rectangular
        bounding box) or triangulation + Smith form.

        INPUT:

        - ``threshold`` -- integer (default: 10000); use the naïve
          algorithm as long as the bounding box is smaller than this

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],      # optional - pynormaliz
            ....:            backend='normaliz').integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = Polyhedron([(1,2,3), (2,3,7), (-2,-3,-11)],    # optional - pynormaliz
            ....:                      backend='normaliz')
            sage: simplex.integral_points()                                # optional - pynormaliz
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)],   # optional - pynormaliz
            ....:                      backend='normaliz')
            sage: simplex.integral_points()                                # optional - pynormaliz
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = Polyhedron([(2,3,7)],                            # optional - pynormaliz
            ....:                    backend='normaliz')
            sage: point.integral_points()                                  # optional - pynormaliz
            ((2, 3, 7),)

            sage: empty = Polyhedron(backend='normaliz')                   # optional - pynormaliz
            sage: empty.integral_points()                                  # optional - pynormaliz
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = Polyhedron(v, backend='normaliz'); simplex     # optional - pynormaliz
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
            sage: len(simplex.integral_points())                           # optional - pynormaliz
            49

        A rather thin polytope for which the bounding box method would
        be a very bad idea (note this is a rational (non-lattice)
        polytope, so the other backends use the bounding box method)::

            sage: P = Polyhedron(vertices=((0, 0), (178933,37121))) + 1/1000*polytopes.hypercube(2)
            sage: P = Polyhedron(vertices=P.vertices_list(),               # optional - pynormaliz
            ....:                backend='normaliz')
            sage: len(P.integral_points())                                 # optional - pynormaliz
            434

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ....:      (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v, backend='normaliz')                    # optional - pynormaliz
            sage: pts1 = P.integral_points()                               # optional - pynormaliz
            sage: all(P.contains(p) for p in pts1)                         # optional - pynormaliz
            True
            sage: pts2 = LatticePolytope(v).points()          # PALP
            sage: for p in pts1: p.set_immutable()                         # optional - pynormaliz
            sage: set(pts1) == set(pts2)                                   # optional - pynormaliz
            True

            sage: timeit('Polyhedron(v, backend='normaliz').integral_points()')   # not tested - random
            625 loops, best of 3: 1.41 ms per loop
            sage: timeit('LatticePolytope(v).points()')       # not tested - random
            25 loops, best of 3: 17.2 ms per loop

        TESTS:

        Test some trivial cases (see :trac:`17937`):

        Empty polyhedron in 1 dimension::

            sage: P = Polyhedron(ambient_dim=1, backend='normaliz')        # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Empty polyhedron in 0 dimensions::

            sage: P = Polyhedron(ambient_dim=0, backend='normaliz')        # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Single point in 1 dimension::

            sage: P = Polyhedron([[3]], backend='normaliz')                # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((3),)

        Single non-integral point in 1 dimension::

            sage: P = Polyhedron([[1/2]], backend='normaliz')              # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Single point in 0 dimensions::

            sage: P = Polyhedron([[]], backend='normaliz')                 # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((),)

        A polytope with no integral points (:trac:`22938`)::

            sage: ieqs = [[1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2],
            ....:         [0, -1, 0, 0], [0, 0, -1, 0],  [0, 0, 0, -1],
            ....:         [-1, -1, -1, -1], [1, 1, 0, 0], [1, 0, 1, 0],
            ....:         [1, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')            # optional - pynormaliz
            sage: P.bounding_box()                                         # optional - pynormaliz
            ((-3/4, -1/2, -1/4), (-1/2, -1/4, 0))
            sage: P.bounding_box(integral_hull=True)                       # optional - pynormaliz
            (None, None)
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Check the polytopes from :trac:`22984`::

            sage: base = [[0, 2, 0, -1, 0, 0, 0, 0, 0],
            ....:         [0, 0, 2, 0, -1, 0, 0, 0, 0],
            ....:         [1, -1, 0, 2, -1, 0, 0, 0, 0],
            ....:         [0, 0, -1, -1, 2, -1, 0, 0, 0],
            ....:         [0, 0, 0, 0, -1, 2, -1, 0, 0],
            ....:         [0, 0, 0, 0, 0, -1, 2, -1, 0],
            ....:         [1, 0, 0, 0, 0, 0, -1, 2, -1],
            ....:         [0, 0, 0, 0, 0, 0, 0, -1, 2],
            ....:         [0, -1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [0, 0, -1, 0, 0, 0, 0, 0, 0],
            ....:         [0, 0, 0, -1, 0, 0, 0, 0, 0],
            ....:         [0, 0, 0, 0, -1, 0, 0, 0, 0],
            ....:         [0, 0, 0, 0, 0, -1, 0, 0, 0],
            ....:         [0, 0, 0, 0, 0, 0, -1, 0, 0],
            ....:         [0, 0, 0, 0, 0, 0, 0, -1, 0],
            ....:         [0, 0, 0, 0, 0, 0, 0, 0, -1],
            ....:         [-1, -1, -1, -1, -1, -1, -1, -1, -1]]

            sage: ieqs = base + [
            ....:         [2, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:         [7, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:         [6, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:         [4, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:         [2, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:         [1, 0, 0, 0, 0, 0, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')            # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((-2, -2, -4, -5, -4, -3, -2, -1),
             (-2, -2, -4, -5, -4, -3, -2, 0),
             (-1, -2, -3, -4, -3, -2, -2, -1),
             (-1, -2, -3, -4, -3, -2, -1, 0),
             (-1, -1, -2, -2, -2, -2, -2, -1),
             (-1, -1, -2, -2, -1, -1, -1, 0),
             (-1, -1, -2, -2, -1, 0, 0, 0),
             (-1, 0, -2, -2, -2, -2, -2, -1),
             (0, -1, -1, -2, -2, -2, -2, -1),
             (0, 0, -1, -1, -1, -1, -1, 0))

            sage: ieqs = base + [
            ....:         [3, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:         [6, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:         [8, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:         [6, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:         [4, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:         [2, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:         [1, 0, 0, 0, 0, 0, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')            # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((-3, -4, -6, -8, -6, -4, -2, -1),
             (-3, -4, -6, -8, -6, -4, -2, 0),
             (-2, -2, -4, -5, -4, -3, -2, -1),
             (-2, -2, -4, -5, -4, -3, -2, 0),
             (-1, -2, -3, -4, -3, -2, -2, -1),
             (-1, -2, -3, -4, -3, -2, -1, 0),
             (-1, -1, -2, -2, -2, -2, -2, -1),
             (-1, -1, -2, -2, -1, -1, -1, 0),
             (-1, -1, -2, -2, -1, 0, 0, 0),
             (-1, 0, -2, -2, -2, -2, -2, -1),
             (0, -1, -1, -2, -2, -2, -2, -1),
             (0, 0, -1, -1, -1, -1, -1, 0))
        """
        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')
        # Trivial cases: polyhedron with 0 or 1 vertices
        if self.n_vertices() == 0:
            return ()
        if self.n_vertices() == 1:
            v = self.vertices_list()[0]
            try:
                return (vector(ZZ, v),)
            except TypeError:  # vertex not integral
                return ()
        # for small bounding boxes, it is faster to naively iterate over the points of the box
        if threshold > 1:
            box_min, box_max = self.bounding_box(integral_hull=True)
            if box_min is None:
                return ()
            box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))
            if box_points < threshold:
                from sage.geometry.integral_points import rectangular_box_points
                return rectangular_box_points(list(box_min), list(box_max), self)
        # Compute with normaliz
        points = []
        cone = self._normaliz_cone
        assert cone
        for g in self._nmz_result(cone, "ModuleGenerators"):
            assert g[-1] == 1
            points.append(vector(ZZ, g[:-1]))
        return tuple(points)

    def integral_points_generators(self):
        r"""
        Return the integral points generators of the polyhedron.

        Every integral point in the polyhedron can be written as a (unique)
        non-negative linear combination of integral points contained in the three
        defining parts of the polyhedron: the integral points (the compact
        part), the recession cone, and the lineality space.

        OUTPUT:

        A tuple consisting of the integral points, the Hilbert basis of the
        recession cone, and an integral basis for the lineality space.

        EXAMPLES:

        Normaliz gives a nonnegative integer basis of the lineality space::

            sage: P = Polyhedron(backend='normaliz',lines=[[2,2]])  # optional - pynormaliz
            sage: P.integral_points_generators()                    # optional - pynormaliz
            (((0, 0),), (), ((1, 1),))

        A recession cone generated by two rays::

            sage: C = Polyhedron(backend='normaliz',rays=[[1,2],[2,1]])  # optional - pynormaliz
            sage: C.integral_points_generators()                         # optional - pynormaliz
            (((0, 0),), ((1, 1), (1, 2), (2, 1)), ())

        Empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')  # optional - pynormaliz
            sage: P.integral_points_generators()      # optional - pynormaliz
            ((), (), ())
        """
        # Trivial cases: polyhedron with 0 vertices
        if self.n_vertices() == 0:
            return ((), (), ())
        # Compute with normaliz
        cone = self._normaliz_cone
        compact_part = []
        recession_cone_part = []
        lineality_part = []
        assert cone
        for g in self._nmz_result(cone, "ModuleGenerators"):
            assert g[-1] == 1
            compact_part.append(vector(ZZ, g[:-1]))

        for g in self._nmz_result(cone, "HilbertBasis"):
            assert g[-1] == 0
            recession_cone_part.append(vector(ZZ, g[:-1]))

        for g in self._nmz_result(cone, "MaximalSubspace"):
            assert g[-1] == 0
            lineality_part.append(vector(ZZ, g[:-1]))

        return tuple(compact_part), tuple(recession_cone_part), tuple(lineality_part)

    def _volume_normaliz(self, measure='euclidean'):
        r"""
        Computes the volume of a polytope using normaliz.

        INPUT:

        - ``measure`` -- (default: 'euclidean') the measure to take. 'euclidean'
          correspond to ``EuclideanVolume`` in normaliz and 'induced_lattice'
          correspond to ``Volume`` in normaliz.

        OUTPUT:

        A float value (when ``measure`` is 'euclidean') or a rational number
        (when ``measure`` is 'induced_lattice').

        .. NOTE::

            This function depends on Normaliz (i.e., the ``pynormaliz`` optional
            package). See the Normaliz documentation for further details.

        EXAMPLES:

        For normaliz, the default is the euclidean volume in the ambient
        space and the result is a float::

            sage: s = polytopes.simplex(3,backend='normaliz')  # optional - pynormaliz
            sage: s._volume_normaliz()                         # optional - pynormaliz
            0.3333333333333333

        The other possibility is to compute the scaled volume where a unimodual
        simplex has volume 1::

            sage: s._volume_normaliz(measure='induced_lattice')  # optional - pynormaliz
            1
            sage: v = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
            sage: cube = Polyhedron(vertices=v,backend='normaliz')  # optional - pynormaliz
            sage: cube._volume_normaliz()  # optional - pynormaliz
            1.0
            sage: cube._volume_normaliz(measure='induced_lattice')  # optional - pynormaliz
            6

        """
        cone = self._normaliz_cone
        assert cone
        if measure == 'euclidean':
            return self._nmz_result(cone, 'EuclideanVolume')
        elif measure == 'induced_lattice':
            n, d = self._nmz_result(cone, 'Volume')
            return ZZ(n) / ZZ(d)

    def _triangulate_normaliz(self):
        r"""
        Gives a triangulation of the polyhedron using normaliz

        OUTPUT:

        A tuple of pairs ``(simplex,simplex_volume)`` used in the
        triangulation.

        .. NOTE::

            This function depends on Normaliz (i.e. the ``pynormaliz`` optional
            package). See the Normaliz documentation for further details.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0,1],[1,0,1],[0,1,1],[1,1,1]],backend='normaliz')  #  optional - pynormaliz
            sage: P._triangulate_normaliz()  #  optional - pynormaliz
            [(0, 1, 2), (1, 2, 3)]
            sage: C1 = Polyhedron(rays=[[0,0,1],[1,0,1],[0,1,1],[1,1,1]],backend='normaliz')  #  optional - pynormaliz
            sage: C1._triangulate_normaliz()  #  optional - pynormaliz
            [(0, 1, 2), (1, 2, 3)]
            sage: C2 = Polyhedron(rays=[[1,0,1],[0,0,1],[0,1,1],[1,1,10/9]],backend='normaliz')  #  optional - pynormaliz
            sage: C2._triangulate_normaliz()  #  optional - pynormaliz
            [(0, 1, 2), (1, 2, 3)]
        """
        cone = self._normaliz_cone
        assert cone
        if self.lines():
            raise NotImplementedError("triangulation of non-compact not pointed polyhedron is not supported")
        if len(self.vertices_list()) >= 2 and self.rays_list():  # A mix of polytope and cone
            raise NotImplementedError("triangulation of non-compact polyhedra that are not cones is not supported")

        data = self._get_nmzcone_data()
        # Recreates a pointed cone. This is a hack and should be fixed once
        # Normaliz accepts compact polyhedron
        # For now, we lose the information about the volume?
        # if self.is_compact():
        #     data['cone'] = data['vertices']
        if not self.is_compact():
            data.pop('vertices', None)
        data.pop('inhom_equations', None)
        data.pop('inhom_inequalities', None)
        cone = self._make_normaliz_cone(data)

        nmz_triangulation = self._nmz_result(cone, "Triangulation")
        triang_indices = tuple(vector(ZZ, s[0]) for s in nmz_triangulation)

        # Get the Normaliz ordering of generators
        if self.is_compact():
            generators = [list(vector(ZZ, g)[:-1]) for g in self._cone_generators(cone)]
        else:
            generators = [list(vector(ZZ, g)) for g in self._cone_generators(cone)]

        # Get the Sage ordering of generators
        if self.is_compact():
            poly_gen = self.vertices_list()
        else:
            poly_gen = self.rays_list()

        # When triangulating, Normaliz uses the indexing of 'Generators' and
        # not necessarily the indexing of the V-representation. So we apply the
        # appropriate relabeling into the V-representation inside sage.
        triangulation = [tuple(sorted([poly_gen.index(generators[i]) for i in s])) for s in triang_indices]

        return triangulation

#########################################################################
class Polyhedron_QQ_normaliz(Polyhedron_normaliz, Polyhedron_QQ):
    r"""
    Polyhedra over `\QQ` with normaliz.

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - pynormaliz
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='normaliz', base_ring=QQ)
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - pynormaliz
    """
    pass

#########################################################################
class Polyhedron_ZZ_normaliz(Polyhedron_QQ_normaliz, Polyhedron_ZZ):
    r"""
    Polyhedra over `\ZZ` with normaliz.

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - pynormaliz
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='normaliz', base_ring=ZZ)
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - pynormaliz
    """
    pass

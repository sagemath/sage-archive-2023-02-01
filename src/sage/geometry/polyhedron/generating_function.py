r"""
Generating Function of Polyhedron's Integral Points

This module provides :func:`generating_function_of_polyhedron` which
computes the generating function of the integral points of a polyhedron.

The main function is accessible via
:meth:`sage.geometry.polyhedron.base.Polyhedron_base.generating_function_of_integral_points`
as well.

Various
=======

AUTHORS:

- Daniel Krenn (2016)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the Austrian Science Fund (FWF): P 24644-N26.


Functions
=========
"""

# *****************************************************************************
# Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import print_function
from __future__ import absolute_import
from six import iteritems, itervalues


def generating_function_of_integral_points(polyhedron, split=False,
                                      result_as_tuple=None, **kwds):
    r"""
    Return the multivariate generating function of the
    integral points of the polyhedron.

    To be precise, this returns

    .. MATH::

        \sum_{(r_0,\dots,r_{d-1}) \in \mathit{polyhedron}\cap \ZZ^d}
        y_0^{r_0} \dots y_{d-1}^{r_{d-1}}.

    INPUT:

    - ``polyhedron`` -- an instance of
      :class:`~sage.geometry.polyhedron.base.Polyhedron_base`
      (see also :mod:`sage.geometry.polyhedron.constructor`).

    - ``split`` -- (default: ``False``) ``False`` computes the generating
      function directly, whereas ``True`` splits the ``polyhedron``
      into several small disjoint polyhedra and adds the results.
      ``split`` may also be a list of disjoint polyhedra.

    - ``result_as_tuple`` -- (default: ``None``) a boolean or ``None``
      specifying whether the output is a (partial) factorization
      (``result_as_tuple=False``) or a sum of such (partial)
      factorizations (``result_as_tuple=True``). By default
      (``result_as_tuple=None``), this is automatically determined.
      If the output is a sum, it is represented as a tuple whose
      entries are the summands.

    - ``indices`` -- (default: ``None``) a list or tuple. If this
      is ``None``, this is automatically determined.

    - ``prefix_variable_name`` -- (default: ``'y'``) a string.
      The variable names of the laurent polynomial ring of the output
      are this string followed by an integer.

    - ``Factorization_sort`` (default: ``False``) and
      ``Factorization_simplify`` (default: ``True``) -- are passed on to
      :class:`sage.structure.factorization.Factorization` when creating
      the result.

    - ``sort_factors`` -- (default: ``False``) a boolean. If set, then
      the factors of the output are sorted such that the numerator is
      first and only then all factors of the denominator. It is ensured
      that the sorting is always the same; use this for doctesting.

    OUTPUT:

    The generating function as a (partial)
    :class:`~sage.structure.factorization.Factorization`
    of the result whose factors are laurent polynomials.
    The result might be a tuple of such factorizations
    (depending on the parameter ``result_as_tuple``) as well.

    .. NOTE::

        At the moment, only polyhedra with nonnegative coordinates
        (i.e. a polyhedon in the nonnegative orthant) are handeled.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.generating_function import generating_function_of_integral_points

    ::

        sage: P2 = (
        ....:   Polyhedron(ieqs=[(0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, -1)]),
        ....:   Polyhedron(ieqs=[(0, -1, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0)]))
        sage: generating_function_of_integral_points(P2[0], sort_factors=True)
        1 * (-y0 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1
        sage: generating_function_of_integral_points(P2[1], sort_factors=True)
        1 * (-y1 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y2 + 1)^-1
        sage: (P2[0] & P2[1]).Hrepresentation()
        (An equation (1, 0, -1) x + 0 == 0,
         An inequality (1, 0, 0) x + 0 >= 0,
         An inequality (0, 1, 0) x + 0 >= 0)
        sage: generating_function_of_integral_points(P2[0] & P2[1], sort_factors=True)
        1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1

    ::

        sage: P3 = (
        ....:   Polyhedron(
        ....:     ieqs=[(0, 0, 0, 0, 1), (0, 0, 0, 1, 0),
        ....:           (0, 0, 1, 0, -1), (-1, 1, 0, -1, -1)]),
        ....:   Polyhedron(
        ....:     ieqs=[(0, 0, -1, 0, 1), (0, 1, 0, 0, -1),
        ....:           (0, 0, 0, 1, 0), (0, 0, 1, 0, 0), (-1, 1, -1, -1, 0)]),
        ....:   Polyhedron(
        ....:     ieqs=[(1, -1, 0, 1, 1), (1, -1, 1, 1, 0),
        ....:           (0, 0, 0, 0, 1), (0, 0, 0, 1, 0), (0, 0, 1, 0, 0),
        ....:           (1, 0, 1, 1, -1), (0, 1, 0, 0, 0), (1, 1, 1, 0, -1)]),
        ....:   Polyhedron(
        ....:     ieqs=[(0, 1, 0, -1, 0), (0, -1, 0, 0, 1),
        ....:           (-1, 0, -1, -1, 1), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0)]),
        ....:   Polyhedron(
        ....:     ieqs=[(0, 1, 0, 0, 0), (0, 0, 1, 0, 0),
        ....:           (-1, -1, -1, 0, 1), (0, -1, 0, 1, 0)]))
        sage: def intersect(I):
        ....:     I = iter(I)
        ....:     result = next(I)
        ....:     for i in I:
        ....:         result &= i
        ....:     return result
        sage: for J in subsets(range(len(P3))):
        ....:     if not J:
        ....:         continue
        ....:     P = intersect([P3[j] for j in J])
        ....:     print('{}: {}'.format(J, P.Hrepresentation()))
        ....:     print(generating_function_of_integral_points(P, sort_factors=True))
        [0]: (An inequality (0, 0, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, -1) x + 0 >= 0,
              An inequality (1, 0, -1, -1) x - 1 >= 0)
        y0 * (-y0 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [1]: (An inequality (0, -1, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, -1, -1, 0) x - 1 >= 0,
              An inequality (1, 0, 0, -1) x + 0 >= 0)
        (-y0^2*y2*y3 - y0^2*y3 + y0*y3 + y0) *
        (-y0 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y3 + 1)^-1 *
        (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1]: (An equation (0, 1, 0, -1) x + 0 == 0,
                 An inequality (1, -1, -1, 0) x - 1 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0 * (-y0 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [2]: (An inequality (-1, 0, 1, 1) x + 1 >= 0,
             An inequality (-1, 1, 1, 0) x + 1 >= 0,
              An inequality (0, 0, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (0, 1, 1, -1) x + 1 >= 0,
              An inequality (1, 0, 0, 0) x + 0 >= 0,
              An inequality (1, 1, 0, -1) x + 1 >= 0)
        (y0^7*y1^6*y2^3*y3^5 + y0^7*y1^5*y2^4*y3^4 + y0^6*y1^7*y2^2*y3^5 -
         y0^7*y1^5*y2^3*y3^4 + y0^6*y1^6*y2^3*y3^4 - y0^6*y1^6*y2^2*y3^5 -
         2*y0^6*y1^6*y2^2*y3^4 - 3*y0^6*y1^5*y2^3*y3^4 - y0^6*y1^5*y2^2*y3^5 -
         y0^6*y1^5*y2^3*y3^3 - y0^6*y1^4*y2^4*y3^3 + y0^6*y1^5*y2^2*y3^4 -
         2*y0^5*y1^6*y2^2*y3^4 - 2*y0^6*y1^4*y2^3*y3^4 - y0^5*y1^6*y2*y3^5 +
         y0^6*y1^5*y2^2*y3^3 + y0^6*y1^4*y2^3*y3^3 - y0^5*y1^5*y2^3*y3^3 -
         y0^6*y1^3*y2^4*y3^3 + y0^6*y1^4*y2^2*y3^4 - y0^5*y1^5*y2^2*y3^4 +
         y0^5*y1^5*y2*y3^5 + 3*y0^5*y1^5*y2^2*y3^3 + y0^6*y1^3*y2^3*y3^3 +
         2*y0^5*y1^5*y2*y3^4 - y0^4*y1^6*y2*y3^4 + 4*y0^5*y1^4*y2^2*y3^4 +
         y0^5*y1^4*y2^3*y3^2 + 3*y0^5*y1^4*y2^2*y3^3 + 4*y0^5*y1^3*y2^3*y3^3 -
         y0^5*y1^4*y2*y3^4 + 3*y0^4*y1^5*y2*y3^4 + y0^5*y1^3*y2^2*y3^4 -
         y0^5*y1^4*y2^2*y3^2 + y0^5*y1^3*y2^3*y3^2 + y0^5*y1^2*y2^4*y3^2 -
         y0^5*y1^4*y2*y3^3 + 2*y0^4*y1^5*y2*y3^3 - 2*y0^5*y1^3*y2^2*y3^3 +
         5*y0^4*y1^4*y2^2*y3^3 + y0^5*y1^2*y2^3*y3^3 - y0^5*y1^3*y2^2*y3^2 -
         y0^5*y1^2*y2^3*y3^2 + 2*y0^4*y1^3*y2^3*y3^2 - 4*y0^4*y1^4*y2*y3^3 +
         2*y0^3*y1^5*y2*y3^3 - y0^5*y1^2*y2^2*y3^3 - y0^4*y1^3*y2^2*y3^3 +
         y0^3*y1^5*y3^4 - y0^4*y1^3*y2*y3^4 - y0^4*y1^4*y2*y3^2 -
         5*y0^4*y1^3*y2^2*y3^2 + y0^3*y1^4*y2^2*y3^2 - y0^4*y1^2*y2^3*y3^2 -
         2*y0^4*y1^3*y2*y3^3 - y0^3*y1^4*y2*y3^3 - 3*y0^4*y1^2*y2^2*y3^3 -
         y0^3*y1^4*y3^4 - y0^4*y1^2*y2^3*y3 + y0^4*y1^3*y2*y3^2 -
         3*y0^3*y1^4*y2*y3^2 - y0^4*y1^2*y2^2*y3^2 - 2*y0^3*y1^3*y2^2*y3^2 -
         y0^4*y1*y2^3*y3^2 - 2*y0^3*y1^4*y3^3 + y0^4*y1^2*y2*y3^3 -
         5*y0^3*y1^3*y2*y3^3 + y0^4*y1^2*y2^2*y3 - y0^3*y1^3*y2^2*y3 +
         y0^4*y1^2*y2*y3^2 - y0^3*y1^3*y2*y3^2 - y0^2*y1^4*y2*y3^2 +
         y0^4*y1*y2^2*y3^2 - 4*y0^3*y1^2*y2^2*y3^2 + y0^3*y1^3*y3^3 -
         2*y0^2*y1^4*y3^3 + y0^3*y1^2*y2*y3^3 + y0^3*y1^3*y2*y3 -
         y0^3*y1*y2^3*y3 + y0^3*y1^3*y3^2 + 5*y0^3*y1^2*y2*y3^2 -
         2*y0^2*y1^3*y2*y3^2 + y0^3*y1*y2^2*y3^2 + y0^2*y1^3*y3^3 +
         y0^3*y1^2*y2*y3 + y0^2*y1^3*y2*y3 + 2*y0^3*y1*y2^2*y3 -
         y0^2*y1^2*y2^2*y3 + 3*y0^2*y1^3*y3^2 + 4*y0^2*y1^2*y2*y3^2 +
         y0^2*y1^2*y3^3 - y0^3*y1*y2*y3 + 4*y0^2*y1^2*y2*y3 +
         2*y0^2*y1*y2^2*y3 + y0^2*y1^2*y3^2 + y0*y1^3*y3^2 +
         2*y0^2*y1*y2*y3^2 + y0^2*y1*y2^2 - y0^2*y1^2*y3 - y0^2*y1*y2*y3 +
         y0*y1^2*y2*y3 + y0^2*y2^2*y3 - y0^2*y1*y3^2 + y0*y1^2*y3^2 -
         y0^2*y1*y2 - y0^2*y1*y3 - y0*y1^2*y3 - y0^2*y2*y3 - 2*y0*y1*y3^2 -
         y0*y1*y2 - 3*y0*y1*y3 - 2*y0*y2*y3 -
         y0*y2 + y0*y3 - y1*y3 + y0 + y3 + 1) *
         (-y1 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y1*y3 + 1)^-1 *
         (-y0*y1*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 *
         (-y0*y2*y3 + 1)^-1 * (-y0*y1^2*y3 + 1)^-1 * (-y0^2*y1*y2*y3 + 1)^-1
        [0, 2]: (An equation (1, 0, -1, -1) x - 1 == 0,
                 An inequality (-1, 1, 1, 0) x + 1 >= 0,
                 An inequality (1, 0, -1, 0) x - 1 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [1, 2]: (An equation (1, -1, -1, 0) x - 1 == 0,
                 An inequality (0, -1, 0, 1) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (1, 0, 0, -1) x + 0 >= 0,
                 An inequality (1, -1, 0, 0) x - 1 >= 0)
        (-y0^2*y2*y3 + y0*y3 + y0) *
        (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 2]: (An equation (0, 1, 0, -1) x + 0 == 0,
                    An equation (1, -1, -1, 0) x - 1 == 0,
                    An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, -1, 0, 0) x - 1 >= 0)
        y0 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [3]: (An inequality (-1, 0, 0, 1) x + 0 >= 0,
              An inequality (0, -1, -1, 1) x - 1 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, 0, -1, 0) x + 0 >= 0)
        (-y0*y1*y3^2 - y0*y3^2 + y0*y3 + y3) *
        (-y3 + 1)^-1 * (-y0*y3 + 1)^-1 *
        (-y1*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 3]: (An equation -1 == 0,)
        0
        [1, 3]: (An equation (1, 0, 0, -1) x + 0 == 0,
                 An inequality (1, -1, -1, 0) x - 1 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0*y3 * (-y0*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 3]: (An equation -1 == 0,)
        0
        [2, 3]: (An equation (0, 1, 1, -1) x + 1 == 0,
                 An inequality (1, 0, -1, 0) x + 0 >= 0,
                 An inequality (-1, 1, 1, 0) x + 1 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0)
        (-y0*y1*y3^2 + y0*y3 + y3) *
        (-y1*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 3]: (An equation -1 == 0,)
        0
        [1, 2, 3]: (An equation (1, 0, 0, -1) x + 0 == 0,
                    An equation (1, -1, -1, 0) x - 1 == 0,
                    An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, -1, 0, 0) x - 1 >= 0)
        y0*y3 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 2, 3]: (An equation -1 == 0,)
        0
        [4]: (An inequality (-1, -1, 0, 1) x - 1 >= 0,
              An inequality (-1, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y2 + 1)^-1 * (-y3 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 4]: (An equation -1 == 0,)
        0
        [1, 4]: (An equation -1 == 0,)
        0
        [0, 1, 4]: (An equation -1 == 0,)
        0
        [2, 4]: (An equation (1, 1, 0, -1) x + 1 == 0,
                 An inequality (-1, 0, 1, 0) x + 0 >= 0,
                 An inequality (1, 0, 0, 0) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0)
        y3 * (-y2 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 4]: (An equation -1 == 0,)
        0
        [1, 2, 4]: (An equation -1 == 0,)
        0
        [0, 1, 2, 4]: (An equation -1 == 0,)
        0
        [3, 4]: (An equation (1, 0, -1, 0) x + 0 == 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (-1, -1, 0, 1) x - 1 >= 0,
                 An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y3 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 3, 4]: (An equation -1 == 0,)
        0
        [1, 3, 4]: (An equation -1 == 0,)
        0
        [0, 1, 3, 4]: (An equation -1 == 0,)
        0
        [2, 3, 4]: (An equation (1, 1, 0, -1) x + 1 == 0,
                    An equation (1, 0, -1, 0) x + 0 == 0,
                    An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 3, 4]: (An equation -1 == 0,)
        0
        [1, 2, 3, 4]: (An equation -1 == 0,)
        0
        [0, 1, 2, 3, 4]: (An equation -1 == 0,)
        0

    ::

        sage: P = Polyhedron(vertices=[[1], [5]])
        sage: P.generating_function_of_integral_points()
        y0^5 + y0^4 + y0^3 + y0^2 + y0

    .. SEEALSO::

        This function is accessible via
        :meth:`sage.geometry.polyhedron.base.Polyhedron_base.generating_function_of_integral_points`
        as well. More examples can be found there.

    TESTS::

        sage: generating_function_of_integral_points(
        ....:     Polyhedron(ieqs=[(0, 0, 1, 0, 0), (-1, 1, -1, 0, 0),
        ....:                      (0, 0, 0, 1, 0), (0, 0, 0, 0, 1)]),
        ....:     sort_factors=True)
        y0 * (-y0 + 1)^-1 * (-y2 + 1)^-1 * (-y3 + 1)^-1 * (-y0*y1 + 1)^-1
        sage: generating_function_of_integral_points(
        ....:     Polyhedron(ieqs=[(0, 0, -1, 0, 1), (0, 0, 1, 0, 0),
        ....:                      (0, 1, 0, 0, -1), (-1, 1, -1, 0, 0),
        ....:                      (0, 0, 0, 1, 0)]),
        ....:     sort_factors=True)
        (-y0^2*y3 + y0*y3 + y0) *
        (-y0 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1

        sage: generating_function_of_integral_points(
        ....:     Polyhedron(ieqs=[(0, 1, 0, -1, 0, 0), (0, 0, 0, 1, 0, 0)],
        ....:                eqns=[(0, 0, 0, 1, 0, -1), (0, 1, 0, 0, -1, 0),
        ....:                      (0, 1, -1, 0, 0, 0)]),
        ....:     sort_factors=True)
        1 * (-y0*y1*y3 + 1)^-1 * (-y0*y1*y2*y3*y4 + 1)^-1

    ::

        sage: G = generating_function_of_integral_points(P2[0], sort_factors=True)
        sage: S = generating_function_of_integral_points(P2[0], sort_factors=True,
        ....:                                       split=True)
        sage: sum(S) == G.value()
        True

        sage: G = generating_function_of_integral_points(P2[1], sort_factors=True)
        sage: S = generating_function_of_integral_points(P2[1], sort_factors=True,
        ....:                                       split=True)
        sage: sum(S) == G.value()
        True

        sage: G = generating_function_of_integral_points(P3[0], sort_factors=True)
        sage: S = generating_function_of_integral_points(P3[0], sort_factors=True,
        ....:                                       split=True)
        sage: sum(S) == G.value()
        True

    ::

        sage: generating_function_of_integral_points(
        ....:     Polyhedron(ieqs=[(0, 0, 1, 0, 0), (-1, 1, -1, 0, 0)]),
        ....:     sort_factors=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: Cannot compute the generating function of
        polyhedra with negative coordinates.
        sage: generating_function_of_integral_points(
        ....:     Polyhedron(ieqs=[(0, 0, -1, 0, 1), (0, 0, 1, 0, 0),
        ....:                      (0, 1, 0, 0, -1), (-1, 1, -1, 0, 0)]),
        ....:     sort_factors=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: Cannot compute the generating function of
        polyhedra with negative coordinates.

    ::

        sage: generating_function_of_integral_points(
        ....:     Polyhedron(ieqs=[(0, 0, 1, 0), (-1, 1, -1, 1),
        ....:                      (0, 0, 0, 1), (0, 1, 0, 0)]),
        ....:     prefix_variable_name='z',
        ....:     sort_factors=True)
        (-z0*z1*z2 - z0*z2 + z0 + z2) *
        (-z0 + 1)^-1 * (-z2 + 1)^-1 * (-z0*z1 + 1)^-1 * (-z1*z2 + 1)^-1
        sage: generating_function_of_integral_points(
        ....:     Polyhedron(ieqs=[(0, 0, 1, 0), (-1, 1, -1, 1),
        ....:                      (0, 0, 0, 1), (0, 1, 0, 0)]),
        ....:     prefix_variable_name='mu',
        ....:     sort_factors=True)
        (-mu0*mu1*mu2 - mu0*mu2 + mu0 + mu2) *
        (-mu0 + 1)^-1 * (-mu2 + 1)^-1 * (-mu0*mu1 + 1)^-1 * (-mu1*mu2 + 1)^-1

    ::

        sage: P = Polyhedron(rays=[(1, sqrt(2)), (0, 1)])
        sage: P.generating_function_of_integral_points()
        Traceback (most recent call last):
        ...
        TypeError: Base ring Symbolic Ring of the polyhedron is not ZZ or QQ.
    """
    import logging
    logger = logging.getLogger(__name__)

    from sage.combinat.permutation import Permutations
    from sage.geometry.polyhedron.constructor import Polyhedron
    from sage.rings.integer_ring import ZZ
    from sage.rings.rational_field import QQ

    if result_as_tuple is None:
        result_as_tuple = split

    if polyhedron.is_empty():
        from sage.structure.factorization import Factorization
        result = Factorization([], unit=0)
        if result_as_tuple:
            return (result,)
        else:
            return result

    if polyhedron.base_ring() not in (ZZ, QQ):
        raise TypeError('Base ring {} of the polyhedron is not '
                        'ZZ or QQ.'.format(polyhedron.base_ring()))

    d = polyhedron.ambient_dim()
    nonnegative_orthant = Polyhedron(ieqs=[dd*(0,) + (1,) + (d-dd)*(0,)
                                           for dd in range(1, d+1)])
    if polyhedron & nonnegative_orthant != polyhedron:
        raise NotImplementedError('Cannot compute the generating function of '
                                  'polyhedra with negative coordinates.')

    logger.info('%s', polyhedron)

    if split is False:
        result = _generating_function_of_integral_points_(polyhedron, **kwds)
        if result_as_tuple:
            return result
        else:
            if len(result) != 1:
                raise ValueError("Cannot unpack result. "
                                 "(Set 'result_as_tuple=True'.)")
            return result[0]

    if d <= 1:
        raise ValueError('Cannot do splitting with only '
                         'dimension {}.'.format(d))

    if split is True:
        split = iter(
            (Polyhedron(
                ieqs=[tuple(1 if i==b else (-1 if i==a or i==0 and a > b else 0)
                            for i in range(d+1))
                      for a, b in zip(pi[:-1], pi[1:])]),
             'b{}'.format(pi[0]-1) +
             ''.join((' <= ' if a < b else ' < ') +
                     'b{}'.format(b-1)
                     for a, b in zip(pi[:-1], pi[1:])))
            for pi in Permutations(d))
    else:
        split = iter((ph, ph.repr_pretty_Hrepresentation(prefix='b'))
                     for ph in split)

    result = []
    for split_polyhedron, pi_log in split:
        logger.info('split polyhedron by %s', pi_log)
        result.append(_generating_function_of_integral_points_(
            polyhedron & split_polyhedron, **kwds))
    if not result_as_tuple:
        raise ValueError("Cannot unpack result."
                         "(Unset 'result_as_tuple=False'.)")
    return sum(result, ())


def _generating_function_of_integral_points_(
        polyhedron, indices=None, **kwds):
    r"""
    Helper function for :func:`generating_function_of_integral_points` which
    does the mid-level stuff.

    TESTS::

        sage: from sage.geometry.polyhedron.generating_function import generating_function_of_integral_points

        sage: generating_function_of_integral_points(  # indirect doctest
        ....:     Polyhedron(ieqs=[(0, 1, 0, 0), (0, -1, 1, 0)],
        ....:                eqns=[(0, -1, -1, 2)]),
        ....:     result_as_tuple=True, sort_factors=True)
        ((-y0^2*y1^2*y2^2 + 1) * (-y1^2*y2 + 1)^-1 *
         (-y0^2*y1^2*y2^2 + 1)^-1 * (-y0^2*y1^2*y2^2 + 1)^-1,
         (-y0^3*y1^3*y2^3 + y0*y1*y2) * (-y1^2*y2 + 1)^-1 *
         (-y0^2*y1^2*y2^2 + 1)^-1 * (-y0^2*y1^2*y2^2 + 1)^-1)
    """
    import logging
    logger = logging.getLogger(__name__)

    logger.info('using polyhedron %s',
                polyhedron.repr_pretty_Hrepresentation(prefix='b'))

    if polyhedron.is_empty():
        from sage.structure.factorization import Factorization
        return (Factorization([], unit=0),)

    Hrepr = polyhedron.Hrepresentation()

    inequalities = tuple(tuple(entry)
                         for entry in Hrepr if entry.is_inequality())
    equations = tuple(tuple(entry)
                      for entry in Hrepr if entry.is_equation())
    if len(inequalities) + len(equations) != len(Hrepr):
        raise ValueError('Cannot handle {}.'.format(polyhedron))

    if not inequalities:
        raise NotImplementedError('no inequality given')

    if indices is None:
        indices = range(len(inequalities[0]) - 1)

    n = len(indices) + 1
    if any(len(e) != n for e in inequalities):
        raise ValueError('Not all coefficient vectors of the inequalities '
                         'have the same length.')
    if any(len(e) != n for e in equations):
        raise ValueError('Not all coefficient vectors of the equations '
                         'have the same length.')

    mods = _generate_mods_(equations)
    logger.debug('splitting by moduli %s', mods)

    return tuple(__generating_function_of_integral_points__(
        indices, inequalities, equations, mod, **kwds) for mod in mods)


def __generating_function_of_integral_points__(
        indices, inequalities, equations, mod,
        prefix_variable_name='y',
        Factorization_sort=False, Factorization_simplify=False,
        sort_factors=False):
    r"""
    Helper function for :func:`generating_function_of_integral_points` which
    does the actual computation of the generating function.

    TESTS::

        sage: from sage.geometry.polyhedron.generating_function import __generating_function_of_integral_points__

        sage: __generating_function_of_integral_points__(
        ....:     (0, 2), [(0, 1, 0)], [(1, -1, 2)],
        ....:     {0: (2, 1)}, sort_factors=True)
        y0 * (-y0^2*y2 + 1)^-1
        sage: __generating_function_of_integral_points__(
        ....:     srange(3), [(0, 1, 0, 0), (0, 0, 1, 0)], [(1, -1, 0, 2)],
        ....:     {0: (2, 1)}, sort_factors=True)
        y0 * (-y1 + 1)^-1 * (-y0^2*y2 + 1)^-1
        sage: __generating_function_of_integral_points__(
        ....:     srange(3), [(0, 1, 0, 0), (0, -1, 1, 0)], [(0, -1, -1, 2)],
        ....:     {0: (2, 1), 1: (2, 1)}, sort_factors=True)
        (-y0^3*y1^3*y2^3 + y0*y1*y2) *
        (-y1^2*y2 + 1)^-1 * (-y0^2*y1^2*y2^2 + 1)^-1 * (-y0^2*y1^2*y2^2 + 1)^-1
    """
    import logging
    logger = logging.getLogger(__name__)

    from .representation import repr_pretty
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.rings.polynomial.omega import _Omega_
    from sage.rings.polynomial.omega import partition
    from sage.structure.factorization import Factorization

    B = LaurentPolynomialRing(
        ZZ,
        tuple(prefix_variable_name + str(k) for k in indices),
        len(indices))

    logger.info('preprocessing %s inequalities and %s equations...',
                 len(inequalities), len(equations))

    extra_factor_mod, rules_mod, inequalities, equations = \
        _prepare_mod_(mod, B, inequalities, equations)

    extra_factor_equations, rules_equations, indices_equations = \
        _prepare_equations_(equations, B)

    inequalities, extra_factor_inequalities, rules_inequalities = \
        _prepare_inequalities_(inequalities, B)

    logger.info('%s inequalities left; using Omega...', len(inequalities))

    numerator = B(1)
    terms = B.gens()
    L = B
    mu = 'mu' if not prefix_variable_name.startswith('mu') else 'nu'
    for i, coeffs in enumerate(inequalities):
        L = LaurentPolynomialRing(L, mu + str(i), sparse=True)
        l = L.gen()
        logger.debug('mapping %s --> %s', l, repr_pretty(coeffs, 0))
        it_coeffs = iter(coeffs)
        numerator *= l**next(it_coeffs)
        assert numerator.parent() == L
        terms = tuple(l**c * t for c, t in zip(it_coeffs, terms))
    assert all(y == t for y, t in
               (zip(B.gens(), terms)[i] for i in indices_equations))
    terms = tuple(t for i, t in enumerate(terms)
                  if i not in indices_equations)

    logger.debug('terms denominator %s', terms)

    def decode_factor(factor):
        D = factor.dict()
        assert len(D) == 1
        exponent, coefficient = next(iteritems(D))
        return coefficient, exponent

    while repr(numerator.parent().gen()).startswith(mu):
        logger.info('applying Omega[%s]...', numerator.parent().gen())
        logger.debug('...on terms denominator %s', terms)
        logger.debug('...(numerator has %s)', numerator.number_of_terms())

        decoded_factors, other_factors = \
            partition((decode_factor(factor) for factor in terms),
                      lambda factor: factor[1] == 0)
        other_factors = tuple(factor[0] for factor in other_factors)
        numerator, factors_denominator = \
            _Omega_(numerator.dict(), tuple(decoded_factors))
        terms = other_factors + factors_denominator

    numerator = \
        (((numerator.subs(rules_inequalities) * extra_factor_inequalities
          ).subs(rules_equations) * extra_factor_equations
         ).subs(rules_mod) * extra_factor_mod)
    terms = tuple(
        t.subs(rules_inequalities).subs(rules_equations).subs(rules_mod)
        for t in terms)

    if sort_factors:
        def key(t):
            D = t.dict().popitem()[0]
            return (-sum(abs(d) for d in D), D)
        terms = sorted(terms, key=key, reverse=True)
    return Factorization([(numerator, 1)] +
                         list((1-t, -1) for t in terms),
                         sort=Factorization_sort,
                         simplify=Factorization_simplify)


def _prepare_inequalities_(inequalities, B):
    r"""
    Split off (simple) inequalities which can be handled better
    without passing them to Omega.

    INPUT:

    - ``inequalitites`` -- a list of tuples

    - ``B`` -- a Laurent polynomial ring

    OUTPUT:

    A triple ``(inequalities, factor, rules)`` with the following properties:

    - ``inequalities`` -- a list of tuples

      Determine the generating function of these inequalities instead
      of the input.

    - ``factor`` -- a Laurent polynomial

      The numerator of the generating function has to be multiplied
      with ``factor`` *after* substituting ``rules``.

    - ``rules`` -- a dictionary mapping Laurent polynomial variables to
      Laurent polynomials

      Substitute ``rules`` into the generating function.

    The generating function of the input ``inequalities`` is equal
    to the generating function of the output ``inequalities``
    in which ``rules`` were substitited and ``factor`` was multiplied.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.generating_function import _prepare_inequalities_

        sage: B = LaurentPolynomialRing(ZZ, 'y', 3)
        sage: _prepare_inequalities_([(0, -1, 1, 0), (2, -1, -1, 1)], B)
        ([(2, -2, -1, 1)], 1, {y2: y2, y1: y1, y0: y0*y1})
        sage: _prepare_inequalities_([(-1, -1, 1, 0), (2, -1, -1, 1)], B)
        ([(1, -2, -1, 1)], y1, {y2: y2, y1: y1, y0: y0*y1})
        sage: _prepare_inequalities_([(2, -1, 1, 0), (2, -1, -1, 1)], B)
        ([(4, -2, -1, 1)], y1^-2, {y2: y2, y1: y1, y0: y0*y1})

    TESTS::

        sage: B = LaurentPolynomialRing(ZZ, 'y', 3)
        sage: _prepare_inequalities_([(1, 1, -1, 0), (1, -1, 0, 1),
        ....:                       (2, -1, -1, 3)], B)
        ([(-3, 2, 1, 3)], y0^-1*y2^-2, {y2: y2, y1: y0*y1*y2, y0: y0*y2})

        sage: B = LaurentPolynomialRing(ZZ, 'y', 4)
        sage: _prepare_inequalities_([(-1, 1, -1, 0, 0)], B)
        ([], y0, {y3: y3, y2: y2, y1: y0*y1, y0: y0})
        sage: _prepare_inequalities_([(0, 0, -1, 0, 1), (0, 0, 1, 0, 0),
        ....:                       (0, 1, 0, 0, -1), (-1, 1, -1, 0, 0)], B)
        ([(1, 1, 0, 0, -1)], y0, {y3: y3, y2: y2, y1: y0*y1*y3, y0: y0})
        sage: _prepare_inequalities_([(-2, 1, -1, 0, 0)], B)
        ([], y0^2, {y3: y3, y2: y2, y1: y0*y1, y0: y0})

        sage: _prepare_inequalities_([(0, -1, 1, 0, 0), (-2, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-3, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: _prepare_inequalities_([(0, -1, 1, 0, 0), (-3, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-2, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
        sage: _prepare_inequalities_([(0, -1, 1, 0, 0), (-2, 0, -1, 0, 1),
        ....:                       (-3, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y2^3*y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: _prepare_inequalities_([(0, -1, 1, 0, 0), (-3, 0, -1, 0, 1),
        ....:                       (-2, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y2^2*y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
        sage: _prepare_inequalities_([(-2, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-3, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y1^2*y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: _prepare_inequalities_([(-3, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-2, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y1^3*y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
        sage: _prepare_inequalities_([(-2, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (-3, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y1^2*y2^3*y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: _prepare_inequalities_([(-3, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (-2, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y1^3*y2^2*y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
    """
    import logging
    logger = logging.getLogger(__name__)

    from itertools import takewhile
    from .representation import repr_pretty
    from sage.graphs.digraph import DiGraph
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
    from sage.rings.integer_ring import ZZ

    inequalities_filtered = []
    chain_links = {}
    for coeffs in inequalities:
        dim = len(coeffs)
        if all(c >= 0 for c in coeffs):
            logger.debug('skipping %s (all coefficients >= 0)',
                         repr_pretty(coeffs, 0))
            continue
        constant = coeffs[0]
        ones = tuple(i+1 for i, c in enumerate(coeffs[1:]) if c == 1)
        mones = tuple(i+1 for i, c in enumerate(coeffs[1:]) if c == -1)
        absgetwo = tuple(i+1 for i, c in enumerate(coeffs[1:]) if abs(c) >= 2)
        if len(ones) == 1 and not mones and not absgetwo:
            if constant < 0:
                # This case could be cleverly skipped...
                inequalities_filtered.append(coeffs)
        elif len(ones) == 1 and len(mones) == 1 and not absgetwo:
            logger.debug('handling %s',
                         repr_pretty(coeffs, 0))
            chain_links[(mones[0], ones[0])] = constant
        else:
            inequalities_filtered.append(coeffs)

    G = DiGraph(chain_links, format='list_of_edges')
    potential = {}
    paths = {}
    D = {}
    inequalities_extra = []
    for i in range(dim):
        D[(i, i)] = 1
    for v in G.topological_sort():
        NP = iter(sorted(((n, potential[n] + chain_links[(n, v)])
                          for n in G.neighbor_in_iterator(v)),
                         key=lambda k: k[1]))
        n, p = next(NP, (None, 0))
        potential[v] = p
        D[(0, v)] = -p
        paths[v] = paths.get(n, ()) + (v,)
        for u in paths[v]:
            D[(u, v)] = 1

        for n, p in NP:
            ell = len(tuple(takewhile(lambda u: u[0] == u[1],
                                      zip(paths[n], paths[v]))))
            coeffs = dim*[0]
            for u in paths[v][ell:]:
                coeffs[u] = 1
            for u in paths[n][ell:]:
                coeffs[u] = -1
            coeffs[0] = p - potential[v]
            inequalities_extra.append(tuple(coeffs))
    T = matrix(ZZ, dim, dim, D)

    inequalities = list(tuple(T*vector(ieq))
                        for ieq in inequalities_filtered) + \
                   inequalities_extra

    rules_pre = iter((y, B({tuple(row[1:]): 1}))
                     for y, row in zip((1,) + B.gens(), T.rows()))
    factor = next(rules_pre)[1]
    rules = dict(rules_pre)

    return inequalities, factor, rules


def _prepare_equations_transformation_(E):
    r"""
    Return a transformation matrix and indices which variables
    in the equation to "eliminate" and deal with later.

    INPUT:

    - ``E`` -- a matrix whose rows represent equations

    OUTPUT:

    A triple ``(TE, indices, indicesn)`` with the following properties:

    - ``TE`` -- a matrix

      This matrix arises from ``E`` by multiplying a transformation matrix
      on the left.

    - ``indices`` -- a sorted tuple of integers representing column indices

      The the sub-matrix of ``TE`` with columns ``indices``
      is the identity matrix.

    - ``indicesn`` -- a sorted tuple of integers representing column indices

      ``indicesn`` contains ``0`` and all indices of the columns of ``E``
      which are non-zero.

    TESTS::

        sage: from sage.geometry.polyhedron.generating_function import _prepare_equations_transformation_

        sage: _prepare_equations_transformation_(matrix([(0, 1, 0, -2)]))
        ([   0 -1/2    0    1], (3,), (0, 1))
        sage: _prepare_equations_transformation_(matrix([(0, 1, -2, 0), (0, 2, 0, -3)]))
        (
        [   0 -1/2    1    0]
        [   0 -2/3    0    1], (2, 3), (0, 1)
        )
    """
    indices_nonzero = tuple(i for i, col in enumerate(E.columns())
                            if i > 0 and not col.is_zero())
    indices = []
    r = 0
    for i in reversed(indices_nonzero):
        indices.append(i)
        r1 = E.matrix_from_columns(indices).rank()
        if r1 > r:
            r = r1
            if len(indices) >= E.nrows():
                break
        else:
            indices = indices[:-1]
    assert len(indices) == E.nrows()
    indices = tuple(reversed(indices))
    indicesn = (0,) + tuple(i for i in indices_nonzero if i not in indices)
    TE = E.matrix_from_columns(indices).inverse() * E
    return TE, indices, indicesn


def _prepare_equations_(equations, B):
    r"""
    Prepare the substitutions coming from "eliminated" variables
    in the given equations.

    INPUT:

    - ``equations`` -- a list of tuples

    - ``B`` -- a Laurent polynomial ring

    OUTPUT:

    A triple ``(factor, rules, indices)`` with the following properties:

    - ``factor`` -- a Laurent polynomial

      The numerator of the generating function has to be multiplied
      with ``factor`` *after* substituting ``rules``.

    - ``rules`` -- a dictionary mapping Laurent polynomial variables to
      Laurent polynomials

      Substitute ``rules`` into the generating function.

    - ``indices`` -- a sorted tuple of integers representing
      indices of Laurent polynomial ring variables

      These are exactly the "eliminated" variables which we take care of
      by ``factor`` and ``rules``.

    The generating function of some inequalities and ``equations`` is equal
    to the generating function of these inequalities
    in which ``rules`` were substitited and ``factor`` was multiplied.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.generating_function import _prepare_equations_

        sage: B = LaurentPolynomialRing(ZZ, 'y', 4)
        sage: _prepare_equations_([(1, 1, 1, -1, 0)], B)
        (y2, {y1: y1*y2, y0: y0*y2}, (2,))
        sage: _prepare_equations_([(0, 1, 0, -1, 0)], B)
        (1, {y0: y0*y2}, (2,))
        sage: _prepare_equations_([(-1, 0, 1, -1, -1), (1, 1, 0, 1, 2)], B)
        (y2^-1, {y1: y1*y2^2*y3^-1, y0: y0*y2*y3^-1}, (2, 3))

    TESTS::

        sage: B = LaurentPolynomialRing(ZZ, 'y', 4)
        sage: _prepare_equations_([(0, 0, 1, 0, -1), (-1, 1, -1, -1, 0)], B)
        (y2^-1, {y1: y1*y2^-1*y3, y0: y0*y2}, (2, 3))

        sage: B = LaurentPolynomialRing(ZZ, 'y', 5)
        sage: _prepare_equations_([(0, 0, 0, 1, 0, -1), (0, 1, 0, 0, -1, 0),
        ....:                    (0, 1, -1, 0, 0, 0)], B)
        (1, {y2: y2*y4, y0: y0*y1*y3}, (1, 3, 4))
    """
    from sage.matrix.constructor import matrix
    from sage.misc.misc_c import prod

    E = matrix(equations)
    if not E:
        return 1, {}, ()

    TE, indices, indicesn = _prepare_equations_transformation_(E)

    gens = (1,) + B.gens()
    z = tuple(gens[i] for i in indices)
    gens_cols = zip(gens, TE.columns())
    rules_pre = iter((y, y * prod(zz**(-c) for zz, c in zip(z, col)))
                     for y, col in (gens_cols[i] for i in indicesn))
    factor = next(rules_pre)[1]
    rules = dict(rules_pre)

    return factor, rules, tuple(i-1 for i in indices)


def _generate_mods_(equations):
    r"""
    Extract the moduli and residue classes implied
    by the equations.

    INPUT:

    - ``equations`` -- a list of tuples

    OUTPUT:

    A tuple where each entry represents one possible configuration.
    Each entry is a dictionary mapping ``i`` to ``(m, r)`` with the following 
    meaning: The ``i``th coordinate of each element of the polyhedron
    has to be congruent to ``r`` modulo ``m``.

    TESTS::

        sage: from sage.geometry.polyhedron.generating_function import _generate_mods_
        sage: _generate_mods_([(0, 1, 1, -2)])
        ({0: (2, 0), 1: (2, 0)}, {0: (2, 1), 1: (2, 1)})
    """
    from sage.arith.misc import lcm
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ
    from sage.rings.rational_field import QQ

    TE, TEi, TEin = _prepare_equations_transformation_(matrix(equations))
    TEin = TEin[1:]
    if TE.base_ring() == ZZ:
        mods = [{}]
    elif TE.base_ring() == QQ:
        m = lcm([e.denominator() for e in TE.list()])
        if m == 1:
            mods = [{}]
        else:
            cols = TE.columns()
            assert all(cols[j][i] == 1 for i, j in enumerate(TEi))
            pre_mods = compositions_mod((tuple(ZZ(cc*m) for cc in cols[i])
                                         for i in TEin),
                                        m, r=(-cc*m for cc in cols[0]),
                                        multidimensional=True)
            mods = tuple({i-1: (aa.modulus(), ZZ(aa))
                          for i, aa in zip(TEin, a) if aa.modulus() > 1}
                         for a in pre_mods)
    else:
        raise TypeError('Equations over ZZ or QQ expected, but got '
                        'equations over {}.'.format(TE.base_ring()))

    return mods


def _prepare_mod_(mod, B, *vecs):
    r"""
    Prepare the substitutions coming from the moduli.

    INPUT:

    - ``mod`` -- a dictionary mapping an index ``i`` to ``(m, r)``

      This is one entry of the output tuple of :func:`_generate_mods_`.

    - ``B`` -- a Laurent polynomial ring

    - ``*vecs`` -- each is a list of tuples

      Typically, two ``vec``-arguments are passed, namely
      ``inequalities`` and ``equations``.

    OUTPUT:

    A tuple ``(factor, rules, *vecs)`` with the following properties:

    - ``factor`` -- a Laurent polynomial

      The numerator of the generating function has to be multiplied
      with ``factor`` *after* substituting ``rules``.

    - ``rules`` -- a dictionary mapping Laurent polynomial variables to
      Laurent polynomials

      Substitute ``rules`` into the generating function.

    - ``*vecs`` -- a list of tuples

      Determine the generating function of these ``vecs`` (``inequalities``,
      ``equations``) instead of the input.

    The generating function of the input ``vecs`` is equal
    to the generating function of the output ``vecs``
    in which ``rules`` were substitited and ``factor`` was multiplied.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.generating_function import _prepare_mod_

        sage: B = LaurentPolynomialRing(ZZ, 'y', 3)
        sage: _prepare_mod_({0: (2, 1)}, B, [(1, -1, 0, 2)])
        (y0, {y2: y2, y1: y1, y0: y0^2}, ((0, -2, 0, 2),))
        sage: _prepare_mod_({0: (2, 1), 1: (2, 1)}, B,
        ....:             [(0, -1, -1, 2)], [(0, -1, 1, 0)])
        (y0*y1, {y2: y2, y1: y1^2, y0: y0^2},
         ((-2, -2, -2, 2),), ((0, -2, 2, 0),))
    """
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
    from sage.rings.integer_ring import ZZ

    if not mod:
        return (1, {}) + vecs

    n = len(B.gens()) + 1

    D = {(i, i): 1 for i in range(n)}
    for i, mr in iteritems(mod):
        D[(i+1, i+1)] = mr[0]
        D[(i+1, 0)] = mr[1]
    T = matrix(ZZ, n, n, D)

    rules_pre = iter((y, B({tuple(row[1:]): 1}))
                     for y, row in zip((1,) + B.gens(), T.columns()))
    factor = next(rules_pre)[1]
    rules = dict(rules_pre)

    vecs = tuple(tuple(tuple(vector(e)*T) for e in vec) for vec in vecs)

    return (factor, rules) + vecs


def compositions_mod(u, m, r=0, multidimensional=False):
    r"""
    Return an iterable of tuples `a` such that `a u^T \equiv r \mod m`.

    INPUT:

    - ``m`` -- the modulus as a positive integer.

    - ``multidimensional`` -- (default: ``False``) a boolean.

    If ``multidimensional=False``:

    - ``u`` -- the coefficients as a tuple.

    - ``r`` -- (default: `0`)
      the remainder as a nonnegative integer.

    If ``multidimensional=True``:

    - ``u`` -- the coefficients as a tuple of tuples (read column-wise).

    - ``r`` -- (default: the zero vector)
      the remainder as a tuple of nonnegative integers.

    OUTPUT:

    An iterable of tuples; all these tuples have the same size as ``u``.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.generating_function import compositions_mod

        sage: list(compositions_mod([1, 1], 2))
        [(0, 0), (1, 1)]
        sage: list(compositions_mod([1, 2, 3], 6))
        [(0, 0, 0), (1, 1, 1), (2, 2, 0), (3, 0, 1), (4, 1, 0), (5, 2, 1)]
        sage: list(compositions_mod([2, 2, 2], 6))
        [(0, 0, 0), (0, 1, 2), (0, 2, 1), (1, 0, 2),
         (1, 1, 1), (1, 2, 0), (2, 0, 1), (2, 1, 0), (2, 2, 2)]

    ::

        sage: list(compositions_mod([(1, 0), (0, 1)], 2,
        ....:                       multidimensional=True))
        [(0, 0)]
        sage: list(compositions_mod([(1, 2), (2, 2), (3, 2)], 6,
        ....:                       multidimensional=True))
        [(0, 0, 0), (1, 1, 1), (2, 2, 2)]

    TESTS::

        sage: list(compositions_mod([1, 0], 2))
        [(0, 0)]
    """
    from sage.modules.free_module_element import vector
    from sage.rings.finite_rings.integer_mod_ring import Zmod

    Z = Zmod(m)
    if not multidimensional:
        u = tuple(vector([Z(uu)]) for uu in u)
        r = vector([Z(r)])
    else:
        u = tuple(vector(Z(uuu) for uuu in uu) for uu in u)
        if r == 0:
            r = vector(Z(0) for _ in range(len(u[0])))
        else:
            r = vector(Z(rr) for rr in r)

    return _compositions_mod_(u, r)


def _compositions_mod_(u, r):
    r"""
    Helper function to :func:`compositions_mod`.

    TESTS::

        sage: from sage.geometry.polyhedron.generating_function import _compositions_mod_
        sage: Z = Zmod(2)
        sage: list(_compositions_mod_((vector([Z(1)]), vector([Z(1)])), vector([Z(0)])))
        [(0, 0), (1, 1)]
    """
    if not u:
        if all(rr == 0 for rr in r):
            yield ()
        return

    from itertools import product
    from sage.arith.srange import srange
    from sage.modules.free_module_element import vector
    from sage.rings.finite_rings.integer_mod_ring import Zmod

    v = u[0]
    m = max(vv.order() for vv in v)
    Z = Zmod(m)
    for j in srange(m):
        for a in _compositions_mod_(u[1:], r - j*v):
            yield (Z(j),) + a

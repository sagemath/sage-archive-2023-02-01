r"""
Matrix Backend for SDP solvers

It stores the SDP data in Sage matrices.  It allow users to specify a base ring
and can store exact SDPs with rational or algebraic data.

The class does not provide a solver.  It can be used as a base class for
user-defined classes implementing solvers.

"""

#*****************************************************************************
#       Copyright (C) 2020 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .generic_sdp_backend import GenericSDPBackend

class MatrixSDPBackend(GenericSDPBackend):

    def __init__(self, base_ring=None):

        super(MatrixSDPBackend, self).__init__()

        if base_ring is None:
            from sage.rings.all import QQ
            base_ring = QQ
        self._base_ring = base_ring

        self.set_sense(+1)

    def set_sense(self, sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.matrix_sdp_backend import MatrixSDPBackend
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver=MatrixSDPBackend)
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        self._sense = sense

    def is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.matrix_sdp_backend import MatrixSDPBackend
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver=MatrixSDPBackend)
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        return self._sense == +1

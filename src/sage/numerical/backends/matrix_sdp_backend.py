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

class MatrixSDPBackend:

    pass

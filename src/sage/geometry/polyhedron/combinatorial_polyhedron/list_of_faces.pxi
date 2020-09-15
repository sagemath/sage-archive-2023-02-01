"""
Cython methods for lists of faces.
"""
# ****************************************************************************
#       Copyright (C) 2020 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

include "sage/geometry/polyhedron/combinatorial_polyhedron/face.pxi"

from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces cimport *

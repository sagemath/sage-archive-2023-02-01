r"""
PathTableaux
"""
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

import sage.combinat.path_tableaux.catalog as path_tableaux

from .path_tableau import PathTableau, PathTableaux, CylindricalDiagram
from .dyck_path import DyckPath, DyckPaths

lazy_import('sage.combinat.path_tableaux.path_tableau', ['PathTableau', 'PathTableaux', 'CylindricalDiagram'])
lazy_import('sage.combinat.path_tableaux.dyck_path', ['DyckPath','DyckPaths'])

del absolute_import

r"""
Examples of Lie Conformal Algebras

We implement the following examples of Lie conformal algebras:

- :mod:`Affine Lie conformal algebra<.affine_lie_conformal_algebra>`
- :mod:`Neveu-Schwarz super Lie conformal
  algebra<.neveu_schwarz_lie_conformal_algebra>`
- :mod:`Virasoro Lie conformal algebra<.virasoro_lie_conformal_algebra>`

AUTHORS:

- Reimundo Heluani (2020-06-15): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .affine_lie_conformal_algebra import AffineLieConformalAlgebra as Affine
from .neveu_schwarz_lie_conformal_algebra import NeveuSchwarzLieConformalAlgebra as NeveuSchwarz
from .virasoro_lie_conformal_algebra import VirasoroLieConformalAlgebra as Virasoro

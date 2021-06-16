r"""
Invariant algebras
"""

# ****************************************************************************
#       Copyright (C) 2021 Trevor K. Karn <karnx018 at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
from sage.categories.algebras import Algebras

class FiniteDimensionalInvariantAlgebra(FiniteDimensionalInvariantModule):

	def __init__(self, R, *args, **kwargs):
		"""
        INPUTS::

        - ``R`` -- an instance of a ``Representation`` of a semigroup `S`
                   acting on the algebra `A`.

        OUTPUTS::

        - ``AS`` -- the invariant algebra of the semigroup action of `S` on `A`, or
                    equivalently, the isotypic component of the representation of
                    `S` carried by `A` corresponding to the trivial character.

        EXAMPLES::
        
            

        TESTS::

        """

        if R._module not in Algebras:

            raise ValueErrror(f'{R._module} is not an algebra')

        if R._module not in Algebras().FiniteDimensional().WithBasis():

            raise NotImplementedError(f'{R._module} must be finite-dimensional with a basis')

        super().__init__(R,*args,**kwargs)
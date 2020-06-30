"""
Lie Conformal Algebras With Basis

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.combinat.free_module import CombinatorialFreeModule

class LieConformalAlgebraWithBasis(CombinatorialFreeModule):
    def __init__(self,R, basis_keys=None, element_class=None, category=None,
                 prefix=None, **kwds):
        """
        Abstract base class for a Lie conformal algebra with a
        preferred basis.

        This class provides no functionality, it simply passes the
        arguments to :class:`CombinatorialFreeModule`.

        EXAMPLES::

            sage: R = VirasoroLieConformalAlgebra(QQbar);R
            The Virasoro Lie conformal algebra over Algebraic Field

        TESTS::

            sage: R = VirasoroLieConformalAlgebra(QQ)
            sage: R.0
            L
            sage: R._repr_generator(R.0)
            'L'
            sage: R = AffineLieConformalAlgebra(QQ, 'A1')
            sage: R.0
            B[alpha[1]]
            sage: R._repr_generator(R.0)
            'B[alpha[1]]'
            sage: R = AffineLieConformalAlgebra(QQ, 'A1', names = ('e', 'h','f'))
            sage: R.0
            e
            sage: R._repr_generator(R.0)
            'e'
        """
        if prefix is None:
            prefix = ''
            kwds['bracket'] = ''
            kwds['string_quotes'] = False

        category = LieConformalAlgebras(R).WithBasis().or_subcategory(category)
        super(LieConformalAlgebraWithBasis,self).__init__(R,
                basis_keys=basis_keys, element_class=element_class,
                category=category, prefix=prefix, names=None, **kwds)

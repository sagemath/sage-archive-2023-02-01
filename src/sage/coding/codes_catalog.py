r"""
Index of code constructions

The ``codes`` object may be used to access the codes that Sage can build.

{INDEX_OF_FUNCTIONS}

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.codes_catalog import *

"""
#*****************************************************************************
#       Copyright (C) 2009 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


# Implementation note:
#
# This module is imported as "codes" in all.py so that codes.<tab> is available
# in the global namespace.

from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import as _lazy_import

from .linear_code import LinearCode
from sage.coding.linear_code import LinearCode

_lazy_import('sage.coding.code_constructions',
        ['BCHCode', 'BinaryGolayCode', 'CyclicCodeFromGeneratingPolynomial',
         'CyclicCode', 'CyclicCodeFromCheckPolynomial', 'DuadicCodeEvenPair',
         'DuadicCodeOddPair', 'ExtendedBinaryGolayCode',
         'ExtendedQuadraticResidueCode', 'ExtendedTernaryGolayCode',
         'from_parity_check_matrix',
         'LinearCodeFromCheckMatrix', #deprecated
         'QuadraticResidueCode', 'QuadraticResidueCodeEvenPair',
         'QuadraticResidueCodeOddPair',
         'random_linear_code',
         'RandomLinearCode', #deprecated
         'ReedSolomonCode', 'TernaryGolayCode',
         'ToricCode', 'WalshCode'])

_lazy_import('sage.coding.bch', 'BCHCode')
_lazy_import('sage.coding.cyclic_code', 'CyclicCode')
_lazy_import('sage.coding.extended_code', 'ExtendedCode')
_lazy_import('sage.coding.golay_code', 'GolayCode')
_lazy_import('sage.coding.grs', 'GeneralizedReedSolomonCode')
_lazy_import('sage.coding.guava', ['QuasiQuadraticResidueCode',
                                    'RandomLinearCodeGuava'])
_lazy_import('sage.coding.hamming_code', 'HammingCode')
_lazy_import('sage.coding.parity_check_code', 'ParityCheckCode')
_lazy_import('sage.coding.punctured_code', 'PuncturedCode')
_lazy_import('sage.coding.reed_muller_code', ['BinaryReedMullerCode',
                                              'ReedMullerCode'])
_lazy_import('sage.coding.subfield_subcode', 'SubfieldSubcode')

from . import decoders_catalog as decoders
from . import encoders_catalog as encoders
from . import bounds_catalog as bounds

_lazy_import('sage.coding','databases')

del _lazy_import
from sage.misc.rest_index_of_methods import gen_rest_table_index as _gen_rest_table_index
import sys as _sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=_gen_rest_table_index(_sys.modules[__name__], only_local_functions=False))

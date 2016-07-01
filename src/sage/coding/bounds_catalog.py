r"""
Index of bounds

The ``codes.bounds`` object may be used to access the bounds that Sage can compute.

{INDEX_OF_FUNCTIONS}

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.bounds_catalog import *
"""
from sage.misc.lazy_import import lazy_import as _lazy_import
_lazy_import("sage.coding.code_bounds", ["codesize_upper_bound",
            "dimension_upper_bound",
            "volume_hamming",
            "gilbert_lower_bound",
            "plotkin_upper_bound",
            "griesmer_upper_bound",
            "elias_upper_bound",
            "hamming_upper_bound",
            "singleton_upper_bound",
            "gv_info_rate",
            "entropy",
            "gv_bound_asymp",
            "hamming_bound_asymp",
            "singleton_bound_asymp",
            "plotkin_bound_asymp",
            "elias_bound_asymp",
            "mrrw1_bound_asymp"])

from sage.misc.rest_index_of_methods import gen_rest_table_index as _gen_rest_table_index
import sys as _sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=_gen_rest_table_index(_sys.modules[__name__], only_local_functions=False))

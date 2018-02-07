from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import as _lazy_import

_lazy_import("sage.coding.code_constructions", ["permutation_action",
            "walsh_matrix"])

_lazy_import("sage.coding.linear_code", [
            "LinearCode",
            "LinearCodeFromVectorSpace",
            "self_orthogonal_binary_codes"])

# Functions removed from the global namespace

_lazy_import('sage.coding.delsarte_bounds','krawtchouk', "Krawtchouk",
    deprecation=(20908, "Krawtchouk will be removed from the global namespace. Please use codes.bounds.krawtchouk instead."))
_lazy_import('sage.coding.delsarte_bounds','krawtchouk', "Kravchuk",
    deprecation=(20908, "Kravchuk will be removed from the global namespace. Please use codes.bounds.krawtchouk instead."))

_lazy_import('sage.coding.delsarte_bounds',
    ["delsarte_bound_hamming_space", "delsarte_bound_additive_hamming_space"],
    deprecation=(20908, "This function will soon be removed from the global namespace. "
                "Please call it using codes.bounds.... instead"))

_lazy_import('sage.coding', 'codes_catalog', 'codes')
_lazy_import('sage.coding', 'channels_catalog', 'channels')

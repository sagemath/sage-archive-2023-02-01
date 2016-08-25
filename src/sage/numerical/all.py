from sage.misc.lazy_import import lazy_import
lazy_import("sage.numerical.optimize",
            ["find_fit", "find_local_maximum", "find_local_minimum",
             "find_root", "linear_program", "minimize", "minimize_constrained"])
lazy_import("sage.numerical.mip", ["MixedIntegerLinearProgram"])
lazy_import("sage.numerical.sdp", ["SemidefiniteProgram"])
lazy_import("sage.numerical.backends.generic_backend", ["default_mip_solver"])
lazy_import("sage.numerical.backends.generic_sdp_backend", ["default_sdp_solver"])

lazy_import("sage.numerical.interactive_simplex_method",
            ["InteractiveLPProblem", "InteractiveLPProblemStandardForm"])

from sage.misc.superseded import deprecated_callable_import

# Deprecation (#17867) -- two lines at the end of
# sage.numerical.interactive_simplex_method should be removed too.
deprecated_callable_import(17867,
                           "sage.numerical.interactive_simplex_method",
                           globals(),
                           locals(),
                           ["LPProblem"],
                           ("This class meant for **educational purposes only** has been renamed to InteractiveLPProblem"))

deprecated_callable_import(17867,
                           "sage.numerical.interactive_simplex_method",
                           globals(),
                           locals(),
                           ["LPProblemStandardForm"],
                           ("This class meant for **educational purposes only** has been renamed to InteractiveLPProblemStandardForm"))

import mpmath

# Patch mpmath to use Cythonized functions
import utils as _utils
#mpmath.libmp.normalize = mpmath.libmp.normalize1 = normalize = _utils.normalize
#mpmath.libmp.from_man_exp = from_man_exp = _utils.from_man_exp

# Also import internal functions
from mpmath.libmp import *

# Main namespace
from mpmath import *

# Utilities
from utils import call, mpmath_to_sage, sage_to_mpmath

# Use mpmath internal functions for constants, to avoid unnecessary overhead
_constants_funcs = {
  'glaisher' : glaisher_fixed,
  'khinchin' : khinchin_fixed,
  'twinprime' : twinprime_fixed,
  'mertens' : mertens_fixed
}

def eval_constant(name, ring):
    prec = ring.precision() + 20
    return ring(_constants_funcs[name](prec)) >> prec


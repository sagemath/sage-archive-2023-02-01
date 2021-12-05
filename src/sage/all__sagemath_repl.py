from .all__sagemath_objects import *

# FIXME: all.py should import from here and remove these imports.
from sage.doctest.all    import *
from sage.repl.all       import *

# For doctesting. These are overwritten later

Integer = int
RealNumber = float

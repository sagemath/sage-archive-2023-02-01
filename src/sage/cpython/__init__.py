# sage.cpython is an ordinary package, not a namespace package.

# This package is imported very early, which is why workarounds/monkey-patching
# are done in this file.

# Make sure that the correct zlib library is loaded. This is needed
# to prevent the system zlib to be loaded instead of the Sage one.
# See https://trac.sagemath.org/ticket/23122
import zlib as _zlib


# Monkey-patch ExtensionFileLoader to allow IPython to find the sources
# of Cython files. See https://trac.sagemath.org/ticket/24681
from importlib.machinery import ExtensionFileLoader as _ExtensionFileLoader
del _ExtensionFileLoader.get_source

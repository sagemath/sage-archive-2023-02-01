r"""
Exceptions

This module defines Sage-specific exceptions.
"""

class OptionalPackageNotFoundError(RuntimeError):
    """
    This class defines the exception that should be raised when a
    function, method, or class cannot detect an optional package that
    it depends on.  When an ``OptionalPackageNotFoundError`` is
    raised, this means one of the following:

    - The required optional package is not installed.

    - The required optional package is installed, but the relevant
      interface to that package is unable to detect the package.

    EXAMPLES::

        sage: from sage.misc.exceptions import OptionalPackageNotFoundError
        sage: def find_package(fav_package):
        ...       try:
        ...           raise OptionalPackageNotFoundError("Unable to detect optional package: %s" % fav_package)
        ...       except OptionalPackageNotFoundError:
        ...           raise
        ...
        sage: find_package("ham and spam")
        Traceback (most recent call last):
        ...
        OptionalPackageNotFoundError: Unable to detect optional package: ham and spam
    """
    pass

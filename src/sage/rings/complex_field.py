r"""
Deprecated import from :trac:`24483`.

TESTS::

    sage: from sage.rings.complex_field import ComplexField
    sage: ComplexField()
"""

from sage.misc.lazy_import import lazy_import

ComplexField = lazy_import('sage.rings.complex_mpfr', 'ComplexField', deprecation=24483)


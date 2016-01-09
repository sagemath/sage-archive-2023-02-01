r"""
Deprecated module (:trac:`16132`)

EXAMPLES::

    sage: from sage.misc.lazy_list import lazy_list
    sage: lazy_list([0,1,2])
    doctest:...: DeprecationWarning: Importing lazy_list from here is
    deprecated. If you need to use it, please import it directly from
    sage.data_structures.lazy_list
    See http://trac.sagemath.org/16132 for details.
    lazy list [0, 1, 2]
"""
from lazy_import import lazy_import

lazy_import('sage.data_structures.lazy_list', ['lazy_list'], deprecation=16132)

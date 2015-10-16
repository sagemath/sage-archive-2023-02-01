"""
Deprecated integer list module

TESTS::

    sage: from sage.combinat.integer_list import IntegerListsLex
    sage: IntegerListsLex(3)
    doctest:...: DeprecationWarning:
    Importing IntegerListsLex from here is deprecated. If you need to use it, please import it directly from sage.combinat.integer_lists
    See http://trac.sagemath.org/18109 for details.
    Integer lists of sum 3 satisfying certain constraints
"""
from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.integer_list_old', '*', deprecation=18109)
lazy_import('sage.combinat.integer_lists', '*', deprecation=18109)

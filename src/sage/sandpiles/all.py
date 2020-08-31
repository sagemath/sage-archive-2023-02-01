"""
Test for deprecations of import in global namespace from :trac:`30479`::

    sage: random_DAG
    doctest:warning...:
    DeprecationWarning:
    Importing random_DAG from here is deprecated. If you need to use it, please import it directly from sage.sandpiles.sandpile
    See https://trac.sagemath.org/30479 for details.
    ...
"""
from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

from .sandpile import (Sandpile,
                       SandpileDivisor,
                       SandpileConfig,
                       firing_graph,
                       parallel_firing_graph,
                       wilmes_algorithm,
                       triangle_sandpile)

lazy_import('sage.sandpiles.examples', 'sandpiles')

lazy_import('sage.sandpiles.sandpile', ['random_DAG'], deprecation=30479)

r"""
Catalog of Algebras

The ``algebras`` object may be used to access examples of various algebras
currently implemented in Sage. Using tab-completion on this object is an
easy way to discover and quickly create the algebras that are available
(as listed here).

Let ``<tab>`` indicate pressing the tab key.  So begin by typing
``groups.<tab>`` to the see primary divisions, followed by (for example)
``groups.matrix.<tab>`` to access various groups implemented as sets of matrices.

- Permutation Groups  (``groups.permutation.<tab>``)

  - :class:`groups.permutation.Symmetric <sage.groups.perm_gps.permgroup_named.SymmetricGroup>`
"""

from sage.algebras.free_algebra import FreeAlgebra as Free
from sage.algebras.iwahori_hecke_algebra import IwahoriHeckeAlgebra as IwahoriHecke
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra as Quaternion

from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.nil_coxeter_algebra', 'NilCoxeterAlgebra', 'NilCoxeter')
lazy_import('sage.algebras.hall_algebra', 'HallAlgebra', 'Hall')
lazy_import('sage.algebras.shuffle_algebra', 'ShuffleAlgebra', 'Shuffle')


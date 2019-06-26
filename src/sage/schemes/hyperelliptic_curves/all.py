from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

from .constructor import HyperellipticCurve
from .kummer_surface import KummerSurface
lazy_import('sage.schemes.hyperelliptic_curves.invariants',
            ['igusa_clebsch_invariants', 'absolute_igusa_invariants_kohel',
             'absolute_igusa_invariants_wamelen', 'clebsch_invariants'],
            deprecation=28064)
from .mestre import (Mestre_conic, HyperellipticCurve_from_invariants)
from . import monsky_washnitzer

del absolute_import

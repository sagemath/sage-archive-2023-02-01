from sage.misc.lazy_import import lazy_import

from pari_group import PariGroup

from matrix_gps.all import *
from abelian_gps.all import *

from perm_gps.all import *

from generic import *

lazy_import('sage.groups.class_function', 'ClassFunction')

from additive_abelian.all import *

lazy_import('sage.groups.conjugacy_classes', ['ConjugacyClass', 'ConjugacyClassGAP'])

lazy_import('sage.groups.free_group', 'FreeGroup')
lazy_import('sage.groups.braid', 'BraidGroup')

import groups_catalog as groups

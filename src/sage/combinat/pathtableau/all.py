from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.tableau',['Tableau','Tableaux'])

lazy_import('sage.combinat.tableau.pathtableaux', 'PathTableau')
lazy_import('sage.combinat.tableau.catalan', 'CatalanTableau')
lazy_import('sage.combinat.tableau.semistandard', 'DualSemistandardTableau')

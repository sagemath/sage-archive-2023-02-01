from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

#lazy_import('sage.combinat.tableau',['Tableau','Tableaux'])

lazy_import('sage.combinat.pathtableau.pathtableaux', ['PathTableau','PathTableaux'])
lazy_import('sage.combinat.pathtableau.catalan', ['CatalanTableau','CatalanTableaux'])
lazy_import('sage.combinat.pathtableau.semistandard', ['DualSemistandardTableau','DualSemistandardTableaux'])

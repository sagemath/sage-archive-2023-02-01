class Catalogue:
    pass

cat = Catalogue()

class Rings:
    pass

import sage.rings.all
rings = Rings()
rings.ZZ = sage.rings.all.ZZ
rings.QQ = sage.rings.all.QQ
rings.PolynomialRing = sage.rings.all.PolynomialRing
rings.LaurentSeriesRing = sage.rings.all.LaurentSeriesRing

cat.rings = rings


import sage.groups.all
class Groups:
    pass
groups = Groups()
groups.PermutationGroup = sage.groups.all.PermutationGroup

cat.groups = groups





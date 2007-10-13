"""
Tables of elliptic curves of given rank

AUTHOR:
    -- William Stein (2007-10-07): initial version
"""

import os

from ell_rational_field import (EllipticCurve_rational_field)

class EllipticCurves:
    def rank(self, rank, tors=0, n=10, labels=False):
        """
        Return a list of at most n non-isogenous curves with given
        rank and torsion order.

        If labels is true, return Cremona labels instead.

        EXAMPLES:
            sage: elliptic_curves.rank(n=5, rank=3, tors=2, labels=True)
            ['59450i1', '59450i2', '61376c1', '61376c2', '65481c1']
            sage: elliptic_curves.rank(n=5, rank=0, tors=5, labels=True)
            ['11a1', '11a3', '38b1', '50b1', '50b2']
            sage: elliptic_curves.rank(n=5, rank=1, tors=7, labels=True)
            ['574i1', '4730k1', '6378c1', '10766h1', '15918w1']
        """
        db = "%s/ellcurves/"%os.environ['SAGE_DATA']
        data = '%s/rank%s'%(db,rank)
        if not os.path.exists(data):
            return []
        v = []
        tors = int(tors)
        for w in open(data).readlines():
            N,iso,num,ainvs,r,t = w.split()
            if tors and tors != int(t):
                continue
            label = '%s%s%s'%(N,iso,num)
            if labels:
                v.append(label)
            else:
                E = EllipticCurve_rational_field(eval(ainvs))
                E._set_rank(r)
                E._set_torsion_order(t)
                E._set_conductor(N)
                E._set_cremona_label(label)
                v.append(E)
            if len(v) >= n:
                break
        return v

elliptic_curves = EllipticCurves()



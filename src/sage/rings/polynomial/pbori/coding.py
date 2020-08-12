from __future__ import print_function

from .PyPolyBoRi import Polynomial
from .nf import symmGB_F2_C
from .ll import ll_encode

try:
    from itertools import ifilter as filter
except ImportError:
    pass


class OccCounter(object):
    def __init__(self):
        self.impl = dict()

    def __getitem__(self, k):
        try:
            return self.impl[k]
        except KeyError:
            return 0

    def increase(self, k):
        try:
            self.impl[k] = self.impl[k] + 1
        except KeyError:
            self.impl[k] = 1

    def uniques(self):
        def filter_fun(k):
            return self.impl[k] == 1
        return filter(filter_fun, self.impl.keys())


def preprocess(I, prot=True):
    def min_gb(I):
        strat = symmGB_F2_C(I, opt_lazy=False, opt_exchange=False, prot=prot,
            selection_size=10000, opt_red_tail=True)
        return list(strat.minimalize_and_tail_reduce())
    I = [Polynomial(p) for p in I]
    lin = [p for p in I if p.deg() == 1]
  # lin_strat=symmGB_F2_C(lin, opt_lazy=False,opt_exchange=False,prot=prot,sele
  # ction_size=10000,opt_red_tail=True)
    lin = min_gb(lin)  # list(lin_strat.minimalize_and_tail_reduce())
    for m in sorted([p.lead() for p in lin]):
        print(m)
    lin_ll = ll_encode(lin)
    square = [p.lead() for p in I if p.deg() == 2 and len(p) == 1]
    assert(len(lin) + len(square) == len(I))
    res = list(lin)
    counter = OccCounter()

    def unique_index(s):
        for idx in s:
            if counter[idx] == 1:
                return idx
        never_come_here = False
        assert never_come_here

    for m in square:
        for idx in m:
            counter.increase(idx)
    systems = dict(((idx, []) for idx in counter.uniques()))
    for m in square:
        u_index = unique_index(m)
        systems[u_index].append(ll_red_nf(m / Variable(u_index), lin_ll))
    rewritings = dict()
    u_var = Variable(u)
    u_im = ll_red_nf(u_var, lin_ll)
    if not u_im.isConstant():
        if u_im != u_var:
            r = u_im.navigation().value()
        else:
            pass
    for u in counter.uniques():
        #print u
        u_var = Variable(u)
        u_im = ll_red_nf(u_var, lin_ll)
        if not u_im.isConstant():
            if u_im != u_var:
                r = u_im.navigation().value()
            else:
                pass
    for u in systems.keys():
        u_var = Variable(u)
        u_im = ll_red_nf(u_var, lin_ll)
        res.extend([u_im * p for p in min_gb(systems[u])])
      #print [u_im*p for p in min_gb(systems[u])]
    print("lin:", len(lin), "res:", len(res), "square:", len(square))
    res = [p for p in (Polynomial(p) for p in res) if not p.is_zero()]
    return res

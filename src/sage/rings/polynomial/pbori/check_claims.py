# encoding: utf-8
r"""
check_claims.py

Created by Michael Brickenstein on 2007-03-05.
Copyright (c) 2007 The PolyBoRi Team.
"""
from __future__ import print_function

import sys
from optparse import OptionParser
#if __name__ == "__main__":
#    import pathadjuster
from .PyPolyBoRi import Polynomial, Monomial, BooleConstant, BooleSet
from .PyPolyBoRi import recursively_insert
from .gbrefs import my_import, load_data, clean_data, load_file
from .blocks import IfThen
from copy import copy
from .ll import ll_encode, ll_red_nf_noredsb, ll_red_nf_redsb


def find_one(p, res=None):
    def zero_nav(n):
        return n.constant() and (not n.terminal_one())
    try:
        p = p.navigation()
    except AttributeError:
        pass
    if res is None:
        res = dict()
    if zero_nav(p):
        raise ValueError
    if p.terminal_one():
        return res
    else_branch = p.else_branch()
    if zero_nav(else_branch):
        res[Monomial(Variable(p.value()))] = 1
        find_one(p.then_branch(), res)
    else:
        res[Monomial(Variable(p.value()))] = 0
        find_one(else_branch, res)
    return res

parser = OptionParser()
NF3 = "nf3"
LINEAR_LEAD_NOREDSB = "ll"
parser.add_option("--method",
                  action="store", dest="method", type="choice",
                  choices=["nf3", "linear-lead-redsb", LINEAR_LEAD_NOREDSB],
                  default="linear-lead-redsb",
                  help="select method")


def my_red_nf(p, strat):
    if p.is_zero():
        return p
    hr = nf3(strat.reduction_strategy, p, p.lead())
    if hr.is_zero():
        return hr
    return red_tail(strat, hr)


def gen_strat(polys):
    polys = [Polynomial(p) for p in polys]
    polys = [p for p in polys if not p.is_zero()]
    assert len(set([p.lead() for p in polys])) == len(polys)
    assert polys
    strat = GroebnerStrategy(polys[0].ring())
    for p in polys:
        print("Adding")
        strat.add_generator(p)
    print("finished")
    return strat


def logicaland(l):
    res = BooleConstant(0)
    for p in l:
        res = 1 + (res + 1) * (p + 1)
    return res


def logicalor(l):
    res = BooleConstant(1)
    for p in l:
        res = res * p
    return res


def proof(ifthen, strat):
    ip = ifthen.ifpart
    it = ifthen.thenpart
    print("proofing:", ifthen)
    c = logicalor([1 + logicaland(ip), logicaland(it)])
    if c.is_zero():
        print("TRUE (trivial)")
        return
    else:
        c = nf3(strat.reduction_strategy, c, c.lead())
        if c.is_zero():
            print("TRUE")
            return
        else:
            print("FALSE")


def proofll(ifthen, reductors, redsb=True, prot=True):

    if prot and (not ifthen.supposedToBeValid):
        print("THIS THEOREM IS NOT SUPPOSED TO BE VALID")
    ip_pre = ifthen.ifpart
    ip = []

    for p in ip_pre:
        p = Polynomial(p)
        if p.is_zero():
            continue
        li = list(p.lead().variables())
        if len(li) == 1 and (not (li[0] in list(Polynomial(reductors).lead().
            variables()))):
            assert not Polynomial(reductors).is_zero()
            lead_index = li[0]
            if redsb:
                p = ll_red_nf_redsb(p, reductors)
                reductors = ll_red_nf_redsb(Polynomial(reductors), BooleSet(p.
                    set()))

            p_nav = p.navigation()
            reductors = recursively_insert(p_nav.else_branch(), p_nav.value(),
                reductors)
        else:
            ip.append(p)
    it = ifthen.thenpart
    if prot:
        print("proofing:", ifthen)
    ip = logicaland(ip)
    for c in it:
        if prot:
            print("proofing part:", c)
        c = logicalor([BooleConstant(1) + ip, c])

        if c.is_zero():
            if prot:
                print("TRUE (trivial)")
            return True
        else:
            c_orig = c
            if redsb:
                c = ll_red_nf_redsb(c, reductors)
            else:
                c = ll_red_nf_noredsb(c, reductors)
            if c.is_zero():
                if prot:
                    print("TRUE")
                return True
            else:
                if prot:
                    print("FAILED")
                    print("can construct COUNTER EXAMPLE with:", find_one(c))
                return False


def to_if_then(p):
    if isinstance(p, IfThen):
        return p
    else:
        return IfThen([], [p])


def main(argv=None):
    (opts, args) = parser.parse_args()
    mydata = load_file(args[0])
    claims = mydata.claims
    if opts.method == NF3:
        strat = gen_strat(mydata.ideal)
        for c in claims:
            proof(to_if_then(c), strat)
        del strat
        try:
            del c
        except NameError:
            pass
    else:
        if opts.method == LINEAR_LEAD_NOREDSB:
            reductors = ll_encode(mydata.ideal)
            for c in claims:
                proofll(to_if_then(c), reductors, redsb=False)
            del reductors
            try:
                del c
            except NameError:
                pass
        else:
            reductors = ll_encode(mydata.ideal, reduce=True)
            for c in claims:
                proofll(to_if_then(c), reductors)
            del reductors
            try:
                del c
            except NameError:
                pass
    return 0

if __name__ == "__main__":
    sys.exit(main())

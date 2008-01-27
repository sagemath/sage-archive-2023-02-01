r"""
Run difficult calculations that test the modular forms
functionality.

There's currently no good system for timing these doctests across
all platforms, so I'm turning these all into comments (so that they
aren't counted against are doctest coverage), noting that we should
use these when (if?) we one day have a good regression testing
system in place. All the functions below have just had the "def"
removed from the beginning of their declaration.
-Craig Citro

from sage.all import *

m=0; t=0; tw=0

pre():
    global m, t, tw
    m = get_memory_usage()
    t = cputime()
    tw = walltime()


post():
    global m,t
    print "total time: %s (wall: %.2f); memory usage diff: %sMB"%(cputime(t),
                 walltime(tw), get_memory_usage() - m)

test1():
    pre()
    for N in range(1,75):
        M = ModularForms(N,2)
        print M
        print M.basis()
    post()

test2():
    pre()
    for N in range(1,30):
        M = ModularForms(Gamma1(N),2)
        print M
        print M.basis()
    post()

test3():
    pre()
    for k in range(2,100):
        M = ModularForms(1,k)
        print M
        print M.basis()
    post()

test4():
    pre()
    for N in range(1,30):
        M = ModularForms(N,4, prec=20)
        print M
        print M.basis()
    post()

test5():
    pre()
    for N in range(1,50):
        M = ModularForms(N,2, prec=30)
        print M
        print M.basis()
    post()

test6():
    pre()
    for N in range(225,230):
        M = ModularForms(N,2,prec=40)
        print M
        print M.basis()
    post()

test7():
    pre()
    for k in range(2,30):
        M = ModularForms(2,k)
        print M
        print M.basis()
    post()

test8():
    pre()
    for k in range(2,20):
        M = ModularForms(Gamma1(3),k)
        print M
        print M.basis()
    post()

test9():
    pre()
    for k in range(2,11):
        M = ModularForms(Gamma1(8),k)
        M.set_precision(M.dimension()+2)
        print M
        print M.basis()
    post()

test10():
    pre()
    for k in range(2,11):
        M = ModularForms(k,k)
        M.set_precision(M.dimension()+2)
        print M
        print M.basis()
    post()

test11():
    pre()
    for i in range(100):
        M = ModularForms(randint(1,100),randint(2,6))
        print M
        print M.basis()
    post()

test12():
    S = CuspForms(23,2)
    print S
    print S.hecke_operator(2)
    print S.hecke_operator(2).matrix()

"""


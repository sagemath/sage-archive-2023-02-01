from misc import cputime

from sage.all import *

def benchmark(n=-1):
    """
    Run a well-chosen range of SAGE commands and record the time it
    takes for each to run.

    INPUT:
        n -- int (default: -1) the benchmark number; the default
             of -1 runs all the benchmarks.
    OUTPUT:
        list -- summary of timings for each benchmark.
        int -- if n == -1, also return the total time
    """

    if isinstance(n, list):
        t = cputime()
        v = [benchmark(m) for m in n]
        return v, cputime(t)

    if n != -1:
        print "Running benchmark %s"%n
        try:
            desc, t = eval("bench%s()"%n)
        except NameError:
            raise RuntimeError, "no benchmark %s"%n
        print desc
        print "Time: %s seconds"%t
        return (n, t, desc)

    t = cputime()
    m = 0
    v = []
    while True:
        try:
            v.append(benchmark(m))
            m += 1
        except RuntimeError:
            break
    return v, cputime(t)

def bench0():
    desc = """Benchmark 0: Factor the following polynomial over
    the rational numbers: (x^97+19*x+1)*(x^103-19*x^97+14)*(x^100-1)"""
    x = polygen(QQ,"x")
    f = (x**97+19*x+1)*(x**103-19*x**97+14)*(x**100-1)
    t = cputime()
    F = f.factor()
    return (desc, cputime(t))

def bench1():
    desc = """Find the Mordell-Weil group of the elliptic curve 5077A using mwrank"""
    E = mwrank_EllipticCurve([0, 0, 1, -7, 6])
    t = cputime()
    g = E.gens()
    return (desc, cputime(t))

def bench2():
    desc = """Some basic arithmetic with very large Integer numbers: '3^1000001 * 19^100001"""
    t = cputime()
    a = ZZ(3)**1000001 * ZZ(19)**100001
    return (desc, cputime(t))

def bench3():
    desc = """Some basic arithmetic with very large Rational numbers: '(2/3)^100001 * (17/19)^100001"""
    t = cputime()
    a = QQ('2/3')**100001 * QQ('17/19')**100001
    return (desc, cputime(t))

def bench4():
    desc = """Rational polynomial arithmetic using SAGE. Compute (x^29+17*x-5)^200."""
    t = cputime()
    f = x**29 + 17*x-5
    a = f**200
    return (desc, cputime(t))

def bench5():
    desc = """Rational polynomial arithmetic using SAGE. Compute (x^19 - 18*x + 1)^50 one hundred times."""
    t = cputime()
    f = x**19 - 18*x + 1
    w = [f**50 for _ in range(100)]
    return (desc, cputime(t))

def bench6():
    desc = """Compute the p-division polynomials of y^2 = x^3 + 37*x - 997 for primes p < 40."""
    E = EllipticCurve([0,0,0,37,-997])
    t = cputime()
    for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
        f = E.division_polynomial(p)
    return (desc, cputime(t))

def bench7():
    desc = """Compute the Mordell-Weil group of y^2 = x^3 + 37*x - 997."""
    E = EllipticCurve([0,0,0,37,-997])
    t = cputime()
    G = E.gens()
    return (desc, cputime(t))

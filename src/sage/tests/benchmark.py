"""
Benchmarks


COMMENTS:

Taken as a whole these benchmarks suggest that by far the fastest math
software is MAGMA, Mathematica, and SAGE (with appropriate tuning and
choices -- it can also be very slow !).  Maxima is very slow, at least
in the form that comes with SAGE (perhaps this is because of using
clisp, at least partly).  GP is slow at some thing and very fast at
others.

"""

from sage.all import *

def avg(X):
    s = sum(X,0)
    return s/float(len(X))


# list any CAS's you don't have here and they will
# never be tested by the run command below.
DO_NOT_HAVE=[]


class Benchmark:
    def run(self,
            systems=['sage', 'magma', 'maxima', 'macaulay2', 'gap', 'gp', 'pari', 'python'],
            timeout=60,
            trials=1,
            sort=False):
        if sort:
            systems.sort()
        print '\n\n\n' + str(self)
        #print "Timeout: %s seconds"%timeout
        print '  %-12s%-12s%-12s%-12s%-12s%15s'%('System', 'min',
                                         'avg', 'max', 'trials', 'cpu or wall')
        for S in systems:
            if S in DO_NOT_HAVE:
                continue
            try:
                X = []
                wall = False
                for i in range(trials):
                    alarm(timeout)
                    t = getattr(self, S)()
                    alarm(0)
                    if isinstance(t, tuple):
                        wall = True
                        t = t[1]
                    X.append(t)
                mn = str(min(X))[:8]
                mx = str(max(X))[:8]
                av = str(avg(X))[:8]
                s = '* %-12s%-12s%-12s%-12s%-12s'%(S, mn, av,
                                                 mx, trials)
                if wall:
                    s += '%15s'%'wall'
                print s
            except KeyboardInterrupt:
                print '%-12sinterrupted (timeout: %s seconds wall time)'%(
                    S, timeout)
            except AttributeError:
                pass
            except Exception, msg:
                print msg

    bench = run

    def runall(self, timeout=60, trials=1):
        self.run(timeout=timeout, systems = ['sage', 'magma', 'gap', 'gp', 'python', 'maple', 'mathematica', 'macaulay2', 'maxima'], trials=trials)

class Divpoly(Benchmark):
    def __init__(self, n):
        self.__n = n

    def __repr__(self):
        return "%s-Division polynomial"%self.__n

    def sage(self):
        n = self.__n
        t = cputime()
        E = EllipticCurve([1,2,3,4,5])
        f = E.division_polynomial(n)
        return cputime(t)

    def magma(self):
        n = self.__n
        t = magma.cputime()
        m = magma('DivisionPolynomial(EllipticCurve([1,2,3,4,5]), %s)'%n)
        return magma.cputime(t)

class PolySquare(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R

    def __repr__(self):
        return 'Square a polynomial of degree %s over %s'%(
            self.__n, self.__R)

    def sage(self):
        R = self.__R
        n = self.__n
        f = R['x'](range(1,n+1))
        t = cputime()
        g = f**2
        return cputime(t)

    def magma(self):
        R = magma(self.__R)
        f = magma('PolynomialRing(%s)![1..%s]'%(R.name(),self.__n))
        t = magma.cputime()
        g = f*f
        return magma.cputime(t)

    def maple(self):
        R = self.__R
        if not (R == ZZ or R == QQ):
            raise NotImplementedError
        n = self.__n
        f = maple(str(R['x'](range(1,n+1))))
        t = walltime()
        g = f*f
        return False, walltime(t)

class MPolynomialPower(Benchmark):
    def __init__(self, nvars=2, exp=10, base=QQ, allow_singular=True):
        self.nvars = nvars
        self.exp = exp
        self.base = base
        self.allow_singular = allow_singular

    def __repr__(self):
        s = 'Compute (x_0 + ... + x_%s)^%s over %s'%(
            self.nvars - 1, self.exp, self.base)
        if self.allow_singular:
            s += ' (use singular for SAGE mult.)'
        return s

    def sage(self):
        R = PolynomialRing(self.base, self.nvars)
        z = sum(R.gens())
        if self.allow_singular:
            z = singular(sum(R.gens()))
            t = walltime()
            w = z**self.exp
            return False, walltime(t)
        t = cputime()
        w = z**self.exp
        return cputime(t)

    def macaulay2(self):
        R = PolynomialRing(self.base, self.nvars, macaulay2=True)
        z = macaulay2(sum(R.gens()))
        t = walltime()
        w = z**self.exp
        return False, walltime(t)

    def maxima(self):
        R = PolynomialRing(self.base, self.nvars)
        z = maxima(str(sum(R.gens())))
        w = walltime()
        f = (z**self.exp).expand()
        return False, walltime(w)

    def maple(self):
        R = PolynomialRing(self.base, self.nvars)
        z = maple(str(sum(R.gens())))
        w = walltime()
        f = (z**self.exp).expand()
        return False, walltime(w)

    def mathematica(self):
        R = PolynomialRing(self.base, self.nvars)
        z = mathematica(str(sum(R.gens())))
        w = walltime()
        f = (z**self.exp).Expand()
        return False, walltime(w)

## this doesn't really expand out -- pari has no function to do so,
## as far as I know.
##     def gp(self):
##         R = PolynomialRing(self.base, self.nvars)
##         z = gp(str(sum(R.gens())))
##         gp.eval('gettime')
##         f = z**self.exp
##         return float(gp.eval('gettime/1000.0'))


    def magma(self):
        R = magma.PolynomialRing(self.base, self.nvars)
        z = R.gen(1)
        for i in range(2,self.nvars+1):
            z += R.gen(i)
        t = magma.cputime()
        w = z**magma(self.exp)
        return magma.cputime(t)



class MPolynomialMult(Benchmark):
    def __init__(self, nvars=2, base=QQ, allow_singular=True):
        if nvars%2:
            nvars += 1
        self.nvars = nvars
        self.base = base
        self.allow_singular = allow_singular

    def __repr__(self):
        s =  'Compute (x_0 + ... + x_%s) * (x_%s + ... + x_%s) over %s'%(
            self.nvars/2 - 1, self.nvars/2, self.nvars, self.base)
        if self.allow_singular:
            s += ' (use singular for SAGE mult.)'
        return s

    def maxima(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = maxima(str(sum(R.gens()[:k])))
        z1 = maxima(str(sum(R.gens()[k:])))
        w = walltime()
        f = (z0*z1).expand()
        return False, walltime(w)

    def maple(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = maple(str(sum(R.gens()[:k])))
        z1 = maple(str(sum(R.gens()[k:])))
        w = walltime()
        f = (z0*z1).expand()
        return False, walltime(w)

    def mathematica(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = mathematica(str(sum(R.gens()[:k])))
        z1 = mathematica(str(sum(R.gens()[k:])))
        w = walltime()
        f = (z0*z1).Expand()
        return False, walltime(w)

##     def gp(self):
##         R = PolynomialRing(self.base, self.nvars)
##         k = int(self.nvars/2)
##         z0 = gp(str(sum(R.gens()[:k])))
##         z1 = gp(str(sum(R.gens()[k:])))
##         gp.eval('gettime')
##         f = z0*z1
##         return float(gp.eval('gettime/1000.0'))

    def sage(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = sum(R.gens()[:k])
        z1 = sum(R.gens()[k:])
        if self.allow_singular:
            z0 = singular(z0)
            z1 = singular(z1)
            t = walltime()
            w = z0*z1
            return False, walltime(t)
        else:
            t = cputime()
            w = z0 * z1
            return cputime(t)

    def macaulay2(self):
        R = PolynomialRing(self.base, self.nvars, macaulay2=True)
        k = int(self.nvars/2)
        z0 = macaulay2(sum(R.gens()[:k]))
        z1 = macaulay2(sum(R.gens()[k:]))
        t = walltime()
        w = z0*z1
        return False, walltime(t)

    def magma(self):
        R = magma.PolynomialRing(self.base, self.nvars)
        z0 = R.gen(1)
        k = int(self.nvars/2)
        for i in range(2,k+1):
            z0 += R.gen(i)
        z1 = R.gen(k + 1)
        for i in range(k+1, self.nvars + 1):
            z1 += R.gen(i)
        t = magma.cputime()
        w = z0 * z1
        return magma.cputime(t)

class MPolynomialMult2(Benchmark):
    def __init__(self, nvars=2, base=QQ, allow_singular=True):
        if nvars%2:
            nvars += 1
        self.nvars = nvars
        self.base = base
        self.allow_singular = allow_singular

    def __repr__(self):
        s =  'Compute (x_1 + 2*x_2 + 3*x_3 + ... + %s*x_%s) * (%s * x_%s + ... + %s*x_%s) over %s'%(
            self.nvars/2, self.nvars/2, self.nvars/2+1, self.nvars/2+1,
            self.nvars+1, self.nvars+1, self.base)
        if self.allow_singular:
            s += ' (use singular for SAGE mult.)'
        return s

##     def gp(self):
##         R = PolynomialRing(self.base, self.nvars)
##         k = int(self.nvars/2)
##         z0 = R(0)
##         z1 = R(0)
##         for i in range(k):
##             z0 += (i+1)*R.gen(i)
##         for i in range(k,self.nvars):
##             z1 += (i+1)*R.gen(i)
##         z0 = gp(str(z0))
##         z1 = gp(str(z1))
##         gp.eval('gettime')
##         f = z0*z1
##         print f
##         return float(gp.eval('gettime/1000.0'))

    def maxima(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i+1)*R.gen(i)
        for i in range(k,self.nvars):
            z1 += (i+1)*R.gen(i)
        z0 = maxima(str(z0))
        z1 = maxima(str(z1))
        w = walltime()
        f = (z0*z1).expand()
        return False, walltime(w)

    def macaulay2(self):
        R = PolynomialRing(self.base, self.nvars, macaulay2=True)
        k = int(self.nvars/2)
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i+1)*R.gen(i)
        for i in range(k,self.nvars):
            z1 += (i+1)*R.gen(i)
        z0 = macaulay2(z0)
        z1 = macaulay2(z1)
        t = walltime()
        w = z0*z1
        return False, walltime(t)

    def maple(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i+1)*R.gen(i)
        for i in range(k,self.nvars):
            z1 += (i+1)*R.gen(i)
        z0 = maple(str(z0))
        z1 = maple(str(z1))
        w = walltime()
        f = (z0*z1).expand()
        return False, walltime(w)

    def mathematica(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i+1)*R.gen(i)
        for i in range(k,self.nvars):
            z1 += (i+1)*R.gen(i)
        z0 = mathematica(str(z0))
        z1 = mathematica(str(z1))
        w = walltime()
        f = (z0*z1).Expand()
        return False, walltime(w)

    def sage(self):
        R = PolynomialRing(self.base, self.nvars)
        k = int(self.nvars/2)
        z0 = R(0)
        z1 = R(0)
        for i in range(k):
            z0 += (i+1)*R.gen(i)
        for i in range(k,self.nvars):
            z1 += (i+1)*R.gen(i)
        if self.allow_singular:
            z0 = singular(z0)
            z1 = singular(z1)
            t = walltime()
            w = z0*z1
            return False, walltime(t)
        else:
            t = cputime()
            w = z0 * z1
            return cputime(t)

    def magma(self):
        R = magma.PolynomialRing(self.base, self.nvars)
        z0 = R.gen(1)
        k = int(self.nvars/2)
        for i in range(2,k+1):
            z0 += magma(i)*R.gen(i)
        z1 = R.gen(k + 1)
        for i in range(k+1, self.nvars + 1):
            z1 += magma(i)*R.gen(i)
        t = magma.cputime()
        w = z0 * z1
        return magma.cputime(t)



class CharPolyTp(Benchmark):
    def __init__(self, N=37,k=2,p=2,sign=1):
        self.N = N
        self.k = k
        self.p = p
        self.sign = sign

    def matrix(self):
        try:
            return self._matrix
        except AttributeError:
            self._matrix = ModularSymbols(group=self.N, weight=self.k, sign=self.sign).T(self.p).matrix()
        return self._matrix

    def __repr__(self):
        return "Compute the charpoly (given the matrix) of T_%s on S_%s(Gamma_0(%s)) with sign %s."%(
            self.p, self.k, self.N, self.sign)

    def sage(self):
        m = self.matrix()
        t = cputime()
        f = m.charpoly()
        return cputime(t)

    def gp(self):
        m = gp(self.matrix())
        gp.eval('gettime')
        f = m.charpoly()
        return float(gp.eval('gettime/1000.0'))

    def pari(self):
        m = pari(self.matrix())
        t = cputime()
        f = m.charpoly()
        return cputime(t)

    def magma(self):
        m = magma(self.matrix())
        t = magma.cputime()
        f = m.CharacteristicPolynomial()
        return magma.cputime(t)




class PolyFactor(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R

    def __repr__(self):
        return "Factor a product of 2 polynomials of degree %s over %s."%(
            self.__n, self.__R)

    def sage(self):
        R = PolynomialRing(self.__R)
        f = R(range(1,self.__n+1))
        g = R(range(self.__n+1,2*(self.__n+1)))
        h = f*g
        t = cputime()
        h.factor()
        return cputime(t)

    def magma(self):
        R = magma(self.__R)
        f = magma('PolynomialRing(%s)![1..%s]'%(R.name(),self.__n))
        g = magma('PolynomialRing(%s)![%s+1..2*(%s+1)]'%(
            R.name(),self.__n,self.__n))
        h = f*g
        t = magma.cputime()
        h.Factorization()
        return magma.cputime(t)

    def gp(self):
        R = PolynomialRing(self.__R)
        f = R(range(1,self.__n+1))
        g = R(range(self.__n+1,2*(self.__n+1)))
        h = f*g
        f = gp(h)
        gp.eval('gettime')
        f.factor()
        return float(gp.eval('gettime/1000.0'))



class SquareInts(Benchmark):
    def __init__(self, base=10, ndigits=10**5):
        self.__ndigits = ndigits
        self.base =base

    def __repr__(self):
        return "Square the integer %s^%s"%(self.base, self.__ndigits)

    def sage(self):
        n = Integer(self.base)**self.__ndigits
        t = cputime()
        m = n**2
        return cputime(t)

    def gp(self):
        n = gp('%s^%s'%(self.base,self.__ndigits))
        gp.eval('gettime')
        m = n**2
        return float(gp.eval('gettime/1000.0'))

    def maxima(self):
        n = maxima('%s^%s'%(self.base,self.__ndigits))
        t = walltime()
        m = n**2
        return False, walltime(t)

    def magma(self):
        n = magma('%s^%s'%(self.base,self.__ndigits))
        t = magma.cputime()
        m = n**2
        return magma.cputime(t)

    def python(self):
        n = self.base**self.__ndigits
        t = cputime()
        m = n**2
        return cputime(t)

    def maple(self):
        n = maple('%s^%s'%(self.base,self.__ndigits))
        t = walltime()
        m = n**2
        return False, walltime(t)

    def gap(self):
        n = gap('%s^%s'%(self.base,self.__ndigits))
        t = walltime()
        m = n**2
        return False, walltime(t)

    def mathematica(self):
        n = mathematica('%s^%s'%(self.base,self.__ndigits))
        t = walltime()
        m = n**2
        return False, walltime(t)


class MatrixSquare(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R

    def __repr__(self):
        return 'Square a matrix of degree %s over %s'%(
            self.__n, self.__R)

    def sage(self):
        R = self.__R
        n = self.__n
        f = MatrixSpace(R,n)(range(n*n))
        t = cputime()
        g = f**2
        return cputime(t)

    def magma(self):
        R = magma(self.__R)
        f = magma('MatrixAlgebra(%s, %s)![0..%s^2-1]'%(
            R.name(),self.__n, self.__n))
        t = magma.cputime()
        g = f*f
        return magma.cputime(t)

    def gp(self):
        n = self.__n
        m = gp('matrix(%s,%s,m,n,%s*(m-1)+(n-1))'%(n,n,n))
        gp('gettime')
        n = m*m
        return float(gp.eval('gettime/1000.0'))

    def gap(self):
        n = self.__n
        m = gap(str([range(n*k,n*(k+1)) for k in range(n)]))
        t = walltime()
        j = m*m
        return False, walltime(t)


class Factorial(Benchmark):
    def __init__(self, n):
        self.__n = n

    def __repr__(self):
        return "Compute the factorial of %s"%self.__n

    def sage(self):
        t = cputime()
        n = factorial(self.__n)
        return cputime(t)

    def magma(self):
        t = magma.cputime()
        n = magma('&*[1..%s]'%self.__n)  # &* is way better than Factorial!!
        return magma.cputime(t)

    def maple(self):
        n = maple(self.__n)
        t = walltime()
        m = n.factorial()
        return False, walltime(t)

    def gp(self):
        gp.eval('gettime')
        n = gp('%s!'%self.__n)
        return float(gp.eval('gettime/1000.0'))

class Fibonacci(Benchmark):
    def __init__(self, n):
        self.__n = n

    def __repr__(self):
        return "Compute the %s-th Fibonacci number"%self.__n

    def sage(self):
        t = cputime()
        n = fibonacci(self.__n)
        return cputime(t)

    def magma(self):
        t = magma.cputime()
        n = magma('Fibonacci(%s)'%self.__n)
        return magma.cputime(t)

    def gap(self):
        n = gap(self.__n)
        t = walltime()
        m = n.Fibonacci()
        return False, walltime(t)

    def mathematica(self):
        n = mathematica(self.__n)
        t = walltime()
        m = n.Fibonacci()
        return False, walltime(t)

    def gp(self):
        gp.eval('gettime')
        n = gp('fibonacci(%s)'%self.__n)
        return float(gp.eval('gettime/1000.0'))



class SEA(Benchmark):
    def __init__(self, p):
        self.__p = p

    def __repr__(self):
        return "Do SEA on an elliptic curve over GF(%s)"%self.__p

    def sage(self):
        E = EllipticCurve([1,2,3,4,5])
        t = walltime()
        n = E.SEA(self.__p)
        return False, walltime(t)

    def magma(self):
        magma(0)
        t = magma.cputime()
        m = magma('#EllipticCurve([GF(%s)|1,2,3,4,5])'%(self.__p))
        return magma.cputime(t)

class MatrixKernel(Benchmark):
    def __init__(self, n, R):
        self.__n = n
        self.__R = R

    def __repr__(self):
        return 'Kernel of a matrix of degree %s over %s'%(
            self.__n, self.__R)

    def sage(self):
        R = self.__R
        n = self.__n
        f = MatrixSpace(R,n,2*n)(range(n*(2*n)))
        t = cputime()
        g = f.kernel()
        return cputime(t)

    def magma(self):
        R = magma(self.__R)
        f = magma('RMatrixSpace(%s, %s, %s)![0..(%s*2*%s)-1]'%(
            R.name(),self.__n, 2*self.__n, self.__n, self.__n))
        t = magma.cputime()
        g = f.Kernel()
        return magma.cputime(t)

    def gp(self):
        n = self.__n
        m = gp('matrix(%s,%s,m,n,%s*(m-1)+(n-1))'%(n,2*n,n))
        gp('gettime')
        n = m.matker()
        return float(gp.eval('gettime/1000.0'))

class ComplexMultiply(Benchmark):
    def __init__(self, bits_prec, times):
        self.__bits_prec = bits_prec
        self.__times = times

    def __repr__(self):
        return "List of multiplies of two complex numbers with %s bits of precision %s times"%(self.__bits_prec, self.__times)

    def sage(self):
        CC = ComplexField(self.__bits_prec)
        s = CC(2).sqrt() + (CC.gen()*2).sqrt()
        t = cputime()
        v = [s*s for _ in range(self.__times)]
        return cputime(t)

    def magma(self):
        # decimal digits (despite magma docs that say bits!!)
        n = int(self.__bits_prec/log(10,2)) + 1
        CC = magma.ComplexField(n)
        s = CC(2).Sqrt() + CC.gen(1).Sqrt()
        t = magma.cputime()
        magma.eval('s := %s;'%s.name())
        v = magma('[s*s : i in [1..%s]]'%self.__times)
        return magma.cputime(t)

    def gp(self):
        n = int(self.__bits_prec/log(10,2)) + 1
        gp.set_real_precision(n)
        gp.eval('s = sqrt(2) + sqrt(2*I);')
        gp.eval('gettime;')
        v = gp('vector(%s,i,s*s)'%self.__times)
        return float(gp.eval('gettime/1000.0'))

class ModularSymbols1(Benchmark):
    def __init__(self, N, k=2):
        self.__N = N
        self.__k = k

    def __repr__(self):
        return 'Presentation for modular symbols on Gamma_0(%s) of weight %s'%(self.__N, self.__k)

    def sage(self):
        t = cputime()
        M = ModularSymbols(self.__N, self.__k)
        return cputime(t)

    def magma(self):
        magma = Magma() # new instance since otherwise modsyms are cached, and cache can't be cleared
        t = magma.cputime()
        M = magma('ModularSymbols(%s, %s)'%(self.__N, self.__k))
        return magma.cputime(t)

class ModularSymbolsDecomp1(Benchmark):
    def __init__(self, N, k=2, sign=1, bnd=10):
        self.N = N
        self.k = k
        self.sign = sign
        self.bnd = bnd

    def __repr__(self):
        return 'Decomposition of modular symbols on Gamma_0(%s) of weight %s and sign %s'%(self.N, self.k, self.sign)

    def sage(self):
        t = cputime()
        M = ModularSymbols(self.N, self.k, sign=self.sign, use_cache=False)
        D = M.decomposition(self.bnd)
        return cputime(t)

    def magma(self):
        m = Magma() # new instance since otherwise modsyms are cached, and cache can't be cleared
        t = m.cputime()
        D = m.eval('Decomposition(ModularSymbols(%s, %s, %s),%s);'%(
            self.N, self.k, self.sign, self.bnd))
        return m.cputime(t)

class EllipticCurveTraces(Benchmark):
    def __init__(self, B):
        self.B = B

    def __repr__(self):
        return "Compute all a_p for the elliptic curve [1,2,3,4,5], for p < %s"%self.B

    def sage(self):
        E = EllipticCurve([1,2,3,4,5])
        t = cputime()
        v = E.anlist(self.B, pari_ints=True)
        return cputime(t)

    def magma(self):
        E = magma.EllipticCurve([1,2,3,4,5])
        t = magma.cputime()
        v = E.TracesOfFrobenius(self.B)
        return magma.cputime(t)

class EllipticCurvePointMul(Benchmark):
    def __init__(self, n):
        self.n = n

    def __repr__(self):
        return "Compute %s*(0,0) on the elliptic curve [0, 0, 1, -1, 0] over QQ"%self.n

    def sage(self):
        E = EllipticCurve([0, 0, 1, -1, 0])
        P = E([0,0])
        t = cputime()
        Q = self.n * P
        return cputime(t)

    def magma(self):
        E = magma.EllipticCurve('[0, 0, 1, -1, 0]')
        P = E('[0,0]')
        t = magma.cputime()
        Q = magma(self.n) * P
        return magma.cputime(t)

    def gp(self):
        E = gp.ellinit('[0, 0, 1, -1, 0]')
        gp.eval('gettime')
        P = gp([0,0])
        Q = E.ellpow(P, self.n)
        return float(gp.eval('gettime/1000.0'))

    def pari(self):
        E = pari('ellinit([0, 0, 1, -1, 0])')
        pari('gettime')
        P = pari([0,0])
        Q = E.ellpow(P, self.n)
        return float(pari('gettime/1000.0'))

class EllipticCurveMW(Benchmark):
    def __init__(self, ainvs):
        self.ainvs = ainvs

    def __repr__(self):
        return "Compute generators for the Mordell-Weil group of the elliptic curve %s over QQ"%self.ainvs

    def sage(self):
        E = EllipticCurve(self.ainvs)
        t = walltime()
        G = E.gens()
        return False, walltime(t)

    def magma(self):
        E = magma.EllipticCurve(str(self.ainvs))
        t = magma.cputime()
        G = E.Generators()
        return magma.cputime(t)



"""
TODO:
   * multiply reals
   * modular degree
   * anlist
   * MW group
   * symbolic det
   * poly factor
   * multivariate poly factor

"""




def suite1():
    PolySquare(10000,QQ).run()
    PolySquare(20000,ZZ).run()
    PolySquare(50000,GF(5)).run()
    PolySquare(20000,Integers(8)).run()

    SquareInts(10,2000000).run()

    MatrixSquare(200,QQ).run()
    MatrixSquare(50,ZZ).run()

    SquareInts(10,150000).run()

    Factorial(2*10**6).run(systems = ['sage', 'magma'])
    Fibonacci(10**6).run()
    Fibonacci(2*10^7).run(systems=["sage", "magma", "mathematica"])

    MatrixKernel(150,QQ).run()

    ComplexMultiply(100000,1000)
    ComplexMultiply(100,100000)
    ComplexMultiply(53,100000)

    PolyFactor(300,ZZ)
    PolyFactor(300,GF(19))
    PolyFactor(700,GF(19))

    PolyFactor(500,GF(49))
    PolyFactor(100,GF(10007^3))

    CharPolyTp(54,4).run()
    CharPolyTp(389,2).run()
    CharPolyTp(389,2,sign=0,p=3).run()
    CharPolyTp(1000,2,sign=1,p=2).run(systems=['sage','magma'])
    CharPolyTp(1,100,sign=1,p=5).run(systems=['sage','magma'])   # SAGE's multimodular really sucks here! (GP is way better, even)
    CharPolyTp(512,sign=1,p=3).run(systems=['sage','magma','gp'])
    CharPolyTp(512,sign=0,p=3).run(systems=['sage','magma','gp'])
    CharPolyTp(1024,sign=1,p=3).run(systems=['sage','magma','gp'])
    CharPolyTp(2006,sign=1,p=2).run(systems=['sage','magma','gp'])
    CharPolyTp(2006,sign=1,p=2).run(systems=['sage','magma'])    # gp takes > 1 minute.

def mpoly():
    # This includes a maxima benchmark.  Note that
    # maxima is *shockingly* slow in comparison to Singular or MAGMA.
    # It is so slow as to be useless, basically, i.e., factor
    # of 5000 slower than Singular on this example!
    MPolynomialPower(nvars=6,exp=10).run()

    main = ['sage', 'magma']  # just the main competitors
    MPolynomialPower(nvars=2,exp=200, allow_singular=False).run(main)
    MPolynomialPower(nvars=5,exp=10, allow_singular=False).run(main)
    MPolynomialPower(nvars=5,exp=30, allow_singular=True).run(main)
    MPolynomialPower(nvars=2,exp=1000, allow_singular=True).run(main)
    MPolynomialPower(nvars=10,exp=10, allow_singular=True).run(main)
    MPolynomialPower(nvars=4,exp=350, base=GF(7), allow_singular=True).run(main)
    MPolynomialMult(200, allow_singular=False).run(main)
    MPolynomialMult(400, allow_singular=True).run(main)
    MPolynomialMult(800, allow_singular=True).run(main)
    MPolynomialMult2(500, allow_singular=True).run(main)


def mpoly_all(include_maple=False):
    """
    Runs benchmarks for multipoly arithmetic on all systems (except
    Maxima, since it is very very slow).  You must have mathematica,
    maple, and magma.

    NOTES:
       * maple is depressingly slow on these benchmarks.
       * Singular (i.e., SAGE) does shockingly well.
       * mathematica is sometimes amazing.
       * macaulay2 is also quite bad (though not as bad as maple).
    """
    systems = ['sage', 'magma', 'mathematica', 'macaulay2']
    if include_maple:
        systems.append('maple')
    MPolynomialMult(200).run(systems=systems)
    MPolynomialMult(400).run(systems=systems)
    MPolynomialMult2(256).run(systems=systems)
    MPolynomialMult2(512).run(systems=systems)
    MPolynomialPower(nvars=4,exp=50).run(systems=systems)   # mathematica wins
    MPolynomialPower(nvars=10,exp=10).run(systems=systems)

def modsym_present():
    ModularSymbols1(2006,2)
    ModularSymbols1(1,50)
    ModularSymbols1(1,100)
    ModularSymbols1(1,150)
    ModularSymbols1(30,8)
    ModularSymbols1(225,4)
    ModularSymbols1(2,50)
    ModularSymbols1(2,100)

def modsym_decomp():
    ModularSymbolsDecomp1(1,24).run()
    ModularSymbolsDecomp1(125,2).run()
    ModularSymbolsDecomp1(389,2).run()
    ModularSymbolsDecomp1(1,100).run()
    ModularSymbolsDecomp1(54,4).run()

def elliptic_curve():
    EllipticCurveTraces(100000).run()
    EllipticCurveTraces(500000).run()
    Divpoly(59).run()
    EllipticCurvePointMul(1000).run()
    EllipticCurvePointMul(2000).run()
    EllipticCurvePointMul(2500).run()      # sage is clearly using the wrong algorithm -- maybe need a balanced rep!?

    # NOTE -- SAGE can also do these using Simon's program, which is
    # *way* *way* faster than MAGMA...
    EllipticCurveMW([5,6,7,8,9]).run()
    EllipticCurveMW([50,6,7,8,9]).run()
    EllipticCurveMW([1, -1, 0, -79, 289]).run(trials=1)   # rank 4
    EllipticCurveMW([0, 0, 1, -79, 342]).run(trials=1)    # rank 5  (SAGE wins)

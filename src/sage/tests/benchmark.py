from sage.all import *

class Benchmark:
    def run(self,
            systems=['sage', 'magma', 'gap', 'gp', 'python'],
            timeout=60):
        print self
        print "(timeout: %s seconds)"%timeout
        for S in systems:
            try:
                alarm(timeout)
                t = getattr(self, S)()
                if isinstance(t, tuple):
                    t = '%s (walltime)'%t[1]
                print '%-12s%s'%(S, t)
            except KeyboardInterrupt:
                print '%-12sinterrupted'%S
            except AttributeError:
                pass
            except Exception, msg:
                print msg
            alarm(0)

    bench = run

    def runall(self, timeout=60):
        self.run(timeout=timeout, systems = ['sage', 'magma', 'gap', 'gp', 'python', 'maple', 'mathematica'])

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
    def __init__(self, ndigits):
        self.__ndigits = ndigits

    def __repr__(self):
        return "Square an integer with about %s digits"%self.__ndigits

    def sage(self):
        n = Integer(10)**self.__ndigits
        t = cputime()
        m = n**2
        return cputime(t)

    def gp(self):
        n = gp('10^%s'%self.__ndigits)
        gp.eval('gettime')
        m = n**2
        return float(gp.eval('gettime/1000.0'))

    def magma(self):
        n = magma('10^%s'%self.__ndigits)
        t = magma.cputime()
        m = n**2
        return magma.cputime(t)

    def python(self):
        n = 10**self.__ndigits
        t = cputime()
        m = n**2
        return cputime(t)

    def maple(self):
        n = maple('10^%s'%self.__ndigits)
        t = walltime()
        m = n**2
        return False, walltime(t)

    def gap(self):
        n = gap('10^%s'%self.__ndigits)
        t = walltime()
        m = n**2
        return False, walltime(t)

    #def mathematica(self):
    #    n = mathematica('10^%s'%self.__ndigits)
    #    t = walltime()
    #    m = n**2
    #    return False, walltime(t)


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
        t = magma.cputime()
        M = magma('ModularSymbols(%s, %s)'%(self.__N, self.__k))
        return magma.cputime(t)

"""
TODO:
   * multiply reals
   * multiply complexes
   * modular degree
   * anlist
   * MW group
   * symbolic det
   * poly factor
   * multivariate poly factor

"""




def suite1():
    Divpoly(59).run()

    PolySquare(10000,QQ).run()
    PolySquare(20000,ZZ).run()
    PolySquare(50000,GF(5)).run()
    PolySquare(20000,Integers(8)).run()

    SquareInts(2000000).run()

    MatrixSquare(200,QQ).run()
    MatrixSquare(50,ZZ).run()

    SquareInts(150000).run()

    Factorial(2*10**6).run(systems = ['sage', 'magma'])

    MatrixKernel(150,QQ).run()

    ComplexMultiply(100000,1000)
    ComplexMultiply(100,100000)
    ComplexMultiply(53,100000)

    PolyFactor(300,ZZ)
    PolyFactor(300,GF(19))
    PolyFactor(700,GF(19))

    PolyFactor(500,GF(49))
    PolyFactor(100,GF(10007^3))


def modsym():
    ModularSymbols1(2006,2)
    ModularSymbols1(1,50)
    ModularSymbols1(1,100)
    ModularSymbols1(1,150)
    ModularSymbols1(30,8)
    ModularSymbols1(225,4)
    ModularSymbols1(2,50)
    ModularSymbols1(2,100)

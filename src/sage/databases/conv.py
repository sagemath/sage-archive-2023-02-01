import os

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.arith import GCD
import sage.misc.db as db
from sage.misc.misc import SAGE_ROOT

# TODO:
# I messed up and now self.d is i and self.i is d,
# i.e., the degree and number are swapped.


PATH = SAGE_ROOT + "/src/tables/modular/gamma0/db/"
if not os.path.exists(PATH):
    os.makedirs(PATH)

class ModularForm:
    def __init__(self, N, d, i, Wq, r, charpolys, disc):
        self.N = N
        self.d = d
        self.i = i
        self.Wq = Wq
        self.r = r
        self.charpolys = charpolys
        self.disc = disc

    def __repr__(self):
        s = "Modular form: level %s, degree %s, number %s"%(
            self.N,self.d,self.i) + \
            ", Wq: %s, an rank bnd: %s"%(
            self.Wq, self.r)
        return s

    def congruence_multiple(self, f, anemic=True):
        """
        Return an integer C such that any prime of congruence between
        this modular form and f is a divisor of C.

        If anemic=True (the default), include only coefficients a_p
        with p coprime to the levels of self and f.

        This function returns the gcd of the resultants of each of the
        charpolys of coefficients of self and f that are both known.
        """
        C = 0
        N = self.N * f.N
        for p in self.charpolys.keys():
            if p in f.charpolys.keys():
                if (not anemic) or (anemic and N%p!=0):
                    ap = self.charpolys[p]
                    bp = f.charpolys[p]
                    C = GCD(C, ap.resultant(bp))
        return C

    def torsion_multiple(self):
        """
        Returns a multiple of the order of the torsion subgroup of the
        corresponding abelian variety.
        """
        T = 0
        for p in self.charpolys.keys():
            if self.N % p != 0:
                f = self.charpolys[p]
                T = GCD(T, long(f(p+1)))
        return T


def conv(N):
    file = "data/%s"%N
    if not os.path.exists(file):
        raise RuntimeError, "Data for level %s does not exist."%N

    F = open(file).read()
    i = F.find(":=")
    if i == -1:
        raise RuntimeError, "Syntax error in file for level %s."%N
    F = F[i+2:]
    TRANS = [("[*", "["), ("*]", "]"), ("<","["), (">","]"), \
             (";",""), ("\n",""), ("^","**")]
    for z,w in TRANS:
        F = F.replace(z,w)
    X = []
    # Define x so the eval below works.
    R = PolynomialRing(RationalField())
    x = R.gen()
    print "starting eval."
    #print "F = ", F
    for f in eval(F):
        print "creating object from f=",f[:4]
        cp = {}
        disc = 0
        for z in f[5]:
            g = R(z[1])
            disc = GCD(disc,g.discriminant())
            cp[z[0]] = g
        X.append(ModularForm(f[0],f[1],f[2],f[3],f[4],cp,disc))
    return X


def newforms(N, recompute=False):
    if N <= 10:
        return []
    p = "%s/%s"%(PATH,N)
    if os.path.exists(p + ".bz2") and not recompute:
        return db.load(p,bzip2=True)
    c = conv(N)
    db.save(c,p, bzip2=True)
    return c


def withdisc(v, degree, D):
    ans = []
    for N in v:
        for f in newforms(N):
            if f.d == degree and f.disc % D == 0:
                print f
                ans.append(f)
        if N%100==0:
            print N
    return ans


# 1 - x^2 - x^3 - x^4 + x^6  --> x^3 - 4*x - 1,  disc 122
# 1 - x^3 - x^4 - x^5 + x^8  --> x^4 - 4*x^2 - x + 1  disc 1957
# 1 + x - x^3 - x^4 - x^5 - x^6 - x^7 + x^9 + x^10 --> x^5 + x^4 - 5*x^3 - 5*x^2 + 4*x + 3 disc 36497


def convert(Nstart, Nstop, recompute=False):
    for N in range(Nstart, Nstop):
        print N
        try:
            newforms(N, recompute=recompute)
        except:
            print "Error for N=",N


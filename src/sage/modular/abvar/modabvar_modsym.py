import modabvar
from   sage.modular.congroup import is_Gamma0, Gamma1, Gamma0, GammaH
from   sage.rings.all import QQ, ZZ, gcd, primes
from   sage.modular.modsym.space import is_ModularSymbolsSpace
from   sage.modular.modsym.modsym import ModularSymbols

def J1(N):
    M = ModularSymbols(Gamma1(N), weight=2, sign=1).cuspidal_submodule()
    return ModAbVar_modsym_plus(M,check=False,
                        repr = 'J_1(%s)'%N)

def J0(N):
    M = ModularSymbols(Gamma0(N), weight=2, sign=1).cuspidal_submodule()
    return ModAbVar_modsym_plus(M,check=False,
                        repr = 'J_0(%s)'%N)

def JH(N, H):
    M = ModularSymbols(GammaH(N, H), weight=2, sign=1).cuspidal_submodule()
    return ModAbVar_modsym_plus(M,check=False,
                        repr = 'J_H(%s,%s)'%(N,H))

class ModAbVar_modsym_plus(modabvar.ModularAbelianVariety):
    def __init__(self, M, check=True, repr=''):
        if check:
            if not is_ModularSymbolsSpace(M):
                raise TypeError, "M must be a modular symbols space"
            if not M.sign() == 1:
                raise ValueError, "M must have sign 1"
            if not M.is_cuspidal():
                raise ValueError, "M must be cuspidal."
            if not M.base_ring() == QQ:
                raise ValueError, "M must be defined over QQ"
            if not M.weight() == 2:
                raise ValueError, "M must have weight 2."
            # etc.
        modabvar.ModularAbelianVariety.__init__(self, QQ, M.level())
        self._M = M
        self._repr = repr

    def _repr_(self):
        return "Modular abelian variety%s%s of dimension %s and level %s over %s"%(
            ' ' if self._repr != '' else '', self._repr,
            self.dimension(), self.level(), self.base_ring())

    def dimension(self):
        return self._M.dimension()

    def hecke_polynomial(self, n, var='x'):
        n = ZZ(n)
        key = (n, str(var))
        try:
            return self._hecke_polynomial[key]
        except AttributeError:
            self._hecke_polynomial = {}
        except KeyError:
            pass
        f = self._M.hecke_polynomial(n, var)
        self._hecke_polynomial[key] = f
        return f

    def torsion_multiple(self, stop=20):
        M = self._M
        if M.weight() != 2 or not is_Gamma0(M.group()):
            raise NotImplementedError
        N = self.level()
        P = [p for p in primes(3,stop) if N %p != 0]
        if len(P) ==0:
            return ZZ(0)
        else:
            return gcd([(self.hecke_polynomial(p))(p+1) for p in P])



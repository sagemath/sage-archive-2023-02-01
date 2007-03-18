## def ModularAbelianVariety(X):
##     try:
##         return X.modular_abelian_variety()
##     except AttributeError:
##         raise ValueError, "No known way to associate a modular abelian variety to %s"%X


def J0(N):
    from sage.modular.congroup import Gamma0
    return Gamma0(N).modular_abelian_variety()

def J1(N):
    from sage.modular.congroup import Gamma1
    return Gamma1(N).modular_abelian_variety()

def JH(N, H):
    from sage.modular.congroup import GammaH
    return GammaH(N, H).modular_abelian_variety()

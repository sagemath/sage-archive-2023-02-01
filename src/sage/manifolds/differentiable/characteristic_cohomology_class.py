from sage.manifolds.differentiable.de_rham_cohomology import DeRhamCohomologyClass

def CharacteristicCohomologyClass():
    pass

class CharacteristicCohomologyClass_base(DeRhamCohomologyClass):
    pass

class CharacteristicCohomologyClass_complex(CharacteristicCohomologyClass_base):
    pass

class CharacteristicCohomologyClass_real_even(CharacteristicCohomologyClass_base):
    pass

class CharacteristicCohomologyClass_real_odd(CharacteristicCohomologyClass_base):
    pass

#*****************************************************************************
# ALGORITHMS
#*****************************************************************************

from sage.misc.fast_methods import Singleton
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.manifolds.differentiable.affine_connection import AffineConnection
from sage.manifolds.differentiable.bundle_connection import BundleConnection

class Algorithm_generic(SageObject):
    r"""
    Algorithm class to generate characteristic forms.
    """
    @cached_method
    def get(self, nabla):
        r"""
        Return the global characteristic form w.r.t. a given connection.
        """
        if isinstance(nabla, AffineConnection):
            vbundle = nabla._domain.tangent_bundle()
        elif isinstance(nabla, BundleConnection):
            vbundle = nabla._vbundle
        else:
            raise TypeError(f'{nabla} must be a connection')
        dom = nabla._domain
        res = dom.mixed_form()
        for frame in dom._get_min_covering(nabla._coefficients):
            cmatrix = [[nabla.curvature_form(i, j, frame)
                        for j in vbundle.irange()]
                       for i in vbundle.irange()]
            res_loc = self.get_local(cmatrix)
            res.set_restriction(res_loc)
        return res

    @abstract_method
    def get_local(self, cmat):
        pass

class ChernAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Chern forms.
    """
    def get_local(self, cmat):
        r"""
        Return the local Chern forms for a given curvature matrix.
        """
        from sage.symbolic.constants import pi
        from sage.libs.pynac.pynac import I

        rk = len(cmat)
        dim = cmat[0][0]._domain._dim
        ran = min(rk, dim//2)
        fac = I / (2*pi)
        res = []
        m = cmat
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / k
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= I / (2*pi)
            m = [[sum(cmat[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        res -= fac * sum(m[i][i] for i in range(rk)) / ran
        return res


class PontryaginAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Pontryagin forms.
    """
    def get_local(self, cmat):
        r"""
        Return the local Pontryagin forms for a given curvature matrix.
        """
        from sage.symbolic.constants import pi

        rk = len(cmat)
        dim = cmat[0][0]._domain._dim
        ran = min(rk//2, dim//4)
        fac = -1 / (2*pi)**2
        res = []
        m = cmat2 = [[sum(cmat[i][l].wedge(cmat[l][j])
                          for l in range(rk))
                      for j in range(rk)] for i in range(rk)]
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / (2*k)
            res.append(fac * c)
            for i in range(rk):
                m[i][i] += c
            fac *= -1 / (2*pi)**2
            m = [[sum(cmat2[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        res -= fac * sum(m[i][i] for i in range(rk)) / (2*ran)
        return res

class EulerAlgorithm(Singleton, Algorithm_generic):
    r"""
    Algorithm class to generate Euler forms.
    """
    def get_local(self, cmat):
        r"""
        Return the local Euler form for a given curvature matrix.
        """
        from sage.symbolic.constants import pi

        rk = len(cmat)
        ran = rk // 2
        m = a = [cmat[i].copy() for i in range(rk)]
        for i in range(0, rk, 2):
            m[i], m[i+1] = m[i+1], m[i]  # swap entries
            for k in range(rk):
                m[k][i+1] = -m[k][i+1]
        for k in range(1, ran):
            c = -sum(m[i][i] for i in range(rk)) / (2*k)
            for i in range(rk):
                m[i][i] += c
            m = [[sum(a[i][l].wedge(m[l][j]) for l in range(rk))
                  for j in range(rk)] for i in range(rk)]
        c = -sum(m[i][i] for i in range(rk)) / (2*rk)
        return (-1/(2*pi))**rk * c

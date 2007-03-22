"""
Cuspidal subgroups of modular abelian varieties

TESTS:
    sage: C = J0(54).cuspidal_subgroup()
    sage: loads(dumps(C)) == C
    True
    sage: D = J0(54).rational_cuspidal_subgroup()
    sage: loads(dumps(D)) == D
    True
"""

from finite_subgroup         import FiniteSubgroup
from sage.rings.all          import infinity, QQ, gcd
from sage.matrix.all         import matrix
from sage.modular.congroup   import is_Gamma0
from sage.modular.cusps      import Cusp

class CuspidalSubgroup_generic(FiniteSubgroup):
    def _compute_generators(self, rational_only=False):
        """
        Return a list of vectors that define elements of the rational
        homology that generate this finite subgroup.

        EXAMPLES:
            sage: J = J0(37)
            sage: C = J.cuspidal_subgroup()
            sage: C._generators()
            [(0, 0, 0, 1/3)]
            sage: J = J0(43)
            sage: C = J.cuspidal_subgroup()
            sage: C._generators()
            [(0, -1/7, 0, 1/7, 0, 2/7)]
            sage: J = J0(22)
            sage: C = J.cuspidal_subgroup()
            sage: C._generators()
            [(0, 0, 0, -1/5), (-1/5, -1/5, 1/5, 2/5), (-1/5, -1/5, 1/5, -2/5)]
            sage: J = J1(13)
            sage: C = J.cuspidal_subgroup()
            sage: len(C._generators())
            11
            sage: C._generators()[:3]
            [(0, 1/19, 1/19, -1/19), (6/19, -2/19, -2/19, -1/19), (4/19, -2/19, -2/19, 0)]
        """
        A = self.abelian_variety()
        M = A.modular_symbols()
        I = M.integral_period_mapping().matrix()
        Amb = M.ambient_module()
        C = Amb.cusps()
        N = Amb.level()
        if rational_only:
            if not is_Gamma0(A.group()):
                raise NotImplementedError, 'computation of rational cusps only implemented in Gamma0 case.'
            if not N.is_squarefree():
                data = [n for n in range(2,N) if gcd(n,N) == 1]
                C = [c for c in C if is_rational_cusp_gamma0(c, N, data)]


        G = [Amb([infinity, alpha]).element() for alpha in C]
        J = matrix(QQ, len(G), Amb.dimension(), G)
        R = (J * I).rows()
        return [x for x in R if x.denominator() != 1]

class CuspidalSubgroup(CuspidalSubgroup_generic):
    def _repr_(self):
        return "Cuspidal subgroup of %s"%self.abelian_variety()

    def _generators(self):
        return self._compute_generators(rational_only = False)



class RationalCuspidalSubgroup(CuspidalSubgroup_generic):
    def _repr_(self):
        return "Rational cuspidal subgroup of %s"%self.abelian_variety()

    def _generators(self):
        return self._compute_generators(rational_only = True)


def is_rational_cusp_gamma0(c, N, data):
    num = c.numerator()
    den = c.denominator()
    for d in data:
        if not c.is_gamma0_equiv(Cusp(num,d*den), N):
            return False
    return True

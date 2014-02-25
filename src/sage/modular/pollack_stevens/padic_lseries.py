r"""
P-adic L-series attached to overconvergent eigensymbols
"""
#*****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.padics.all import pAdicField
from sage.rings.all import ZZ, QQ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.big_oh import O
from sage.rings.arith import binomial, gcd, kronecker

from sage.structure.sage_object import SageObject
from sigma0 import Sigma0
from fund_domain import M2Z

class pAdicLseries(SageObject):
    r"""
    The `p`-adic `L`-series associated to an overconvergent eigensymbol.

    INPUT:
    
    - ``symb`` -- overconvergent eigensymbol
    - ``gamma`` -- topological generator of `1 + pZ_p`
    - ``quadratic_twist`` -- conductor of quadratic twist `\chi`, default 1
    - ``precision`` -- if None is specified, the correct precision bound is
      computed and the answer is returned modulo that accuracy

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
        sage: E = EllipticCurve('37a')
        sage: p = 5
        sage: prec = 4
        sage: phi = ps_modsym_from_elliptic_curve(E)
        sage: phi_stabilized = phi.p_stabilize(p,20)
        sage: Phi = phi_stabilized.lift(p,prec,algorithm='stevens',eigensymbol=True)
        sage: L = pAdicLseries(Phi)
        sage: L[1]
        2 + 3*5 + O(5^3)
        sage: L[0]
        O(5^3)

    Using the existing algorithm in Sage, it seems we're off by a factor of 2:

        sage: L = E.padic_lseries(5)
        sage: L.series(4)[1]
        1 + 4*5 + 2*5^2 + O(5^3)

    But here, we're correct without the factor of 2:

        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
        sage: E = EllipticCurve('57a')
        sage: p = 5
        sage: prec = 4
        sage: phi = ps_modsym_from_elliptic_curve(E)
        sage: phi_stabilized = phi.p_stabilize(p,M = prec+3)
        sage: Phi = phi_stabilized.lift(p=5, M=prec, alpha=None, algorithm='stevens', eigensymbol=True)
        sage: L = pAdicLseries(Phi)
        sage: L[1]
        3*5 + 5^2 + O(5^3)

        sage: L1 = E.padic_lseries(5)
        sage: L1.series(4)[1]
        3*5 + 5^2 + O(5^3)

    An example of a `p`-adic `L`-series associated to a modular abelian surface:
        
        sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
        sage: A = ModularSymbols(103,2,1).cuspidal_submodule().new_subspace().decomposition()[0]
        sage: p = 19
        sage: prec = 4
        sage: phi = ps_modsym_from_simple_modsym_space(A)
        sage: ap = phi.Tq_eigenvalue(p,prec)
        sage: c1,c2 = phi.completions(p,prec)
        sage: phi1,psi1 = c1
        sage: phi2,psi2 = c2
        sage: phi1p = phi1.p_stabilize_and_lift(p,ap = psi1(ap), M = prec, algorithm='stevens') # long time
        sage: L1 = pAdicLseries(phi1p) # long time
        sage: phi2p = phi2.p_stabilize_and_lift(p,ap = psi2(ap), M = prec, algorithm='stevens') # long time
        sage: L2  = pAdicLseries(phi2p) # long time
        sage: L1[1]*L2[1] # long time
        13 + 9*19 + O(19^2)
    """
    def __init__(self, symb, gamma=None, quadratic_twist=1, precision=None):
        r"""

        EXAMPLE::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('37a')
            sage: p = 37
            sage: prec = 3
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: Phi = phi.lift(p,prec,algorithm='stevens',eigensymbol=True)
            sage: L = pAdicLseries(Phi)
            sage: L[1]
            4 + 37 + O(37^2)

            sage: TestSuite(L).run()
        """
        self._coefficients = {}
        
        if symb.parent().prime() == None:
            raise ValueError ("Not a p-adic overconvergent modular symbol.")
        
        self._symb = symb

        if gamma == None:
            gamma = 1 + self._symb.parent().prime()

        self._gamma = gamma
        self._quadratic_twist = quadratic_twist
        self._precision = precision

    def __getitem__(self, n):
        """
        Returns the `n`-th coefficient of the `p`-adic `L`-series
        
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('57a')
            sage: p = 5
            sage: prec = 4
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec+3)
            sage: Phi = phi_stabilized.lift(p=p,M=prec,alpha=None,algorithm='stevens',eigensymbol=True)
            sage: L = pAdicLseries(Phi)
            sage: L[1]
            3*5 + 5^2 + O(5^3)

            sage: L1 = E.padic_lseries(5)
            sage: L1.series(4)[1]
            3*5 + 5^2 + O(5^3)
            
        """
        if self._coefficients.has_key(n):
            return self._coefficients[n]
        else:
            p = self.prime()
            symb = self.symb()
            ap = symb.Tq_eigenvalue(p)
            gamma = self._gamma
            precision = self._precision
            
            S = QQ[['z']]
            z = S.gen()
            M = symb.precision_absolute()
            K = pAdicField(p, M)
            dn = 0
            if n == 0:
                precision = M
                lb = [1] + [0 for a in range(M-1)]
            else:
                lb = log_gamma_binomial(p, gamma, z, n, 2*M)
                if precision == None:
                    precision = min([j + lb[j].valuation(p) for j in range(M, len(lb))])
                lb = [lb[a] for a in range(M)]

            for j in range(len(lb)):
                cjn = lb[j]
                temp = sum((ZZ(K.teichmuller(a))**(-j)) * self._basic_integral(a, j) for a in range(1, p))
                dn = dn + cjn*temp
            self._coefficients[n] = dn + O(p**precision)
            return self._coefficients[n]

    def __cmp__(self, other):
        r"""
        Compare self and other.

        EXAMPLE::

            sage: E = EllipticCurve('11a')
            sage: S = sage.modular.pollack_stevens.space.ps_modsym_from_elliptic_curve(E)
            sage: SS = S.lift(11, M=10, algorithm='stevens')
            sage: L = pAdicLseries(SS)
            sage: L == loads(dumps(L)) # indirect doctest
            True
        """
        return cmp(type(self), type(other)) \
            or cmp(self._symb, other._symb) \
            or cmp(self._quadratic_twist, other._quadratic_twist) \
            or cmp(self._gamma, other._gamma) \
            or cmp(self._precision, other._precision)

    def symb(self):
        r"""
        Returns the overconvergent modular symbol

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('37a')
            sage: p = 5
            sage: prec = 6
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec)
            sage: Phi = phi_stabilized.lift(p=p,M=prec,alpha=None,algorithm='stevens',eigensymbol=True)
            sage: L = pAdicLseries(Phi)
            sage: L.symb() is Phi
            True
        """
        return self._symb

    def prime(self):
        r"""
        Returns the prime associatd to the OMS

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('37a')
            sage: p = 5
            sage: prec = 6
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec)
            sage: Phi = phi_stabilized.lift(p,prec,None,algorithm='stevens',eigensymbol=True)
            sage: L = pAdicLseries(Phi)
            sage: L.prime()
            5
        """
        return self._symb.parent().prime()
    
    def quadratic_twist(self):
        r"""
        Returns the discriminant of the quadratic twist

        EXAMPLES::
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('37a')
            sage: p = 5
            sage: prec = 6
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec)
            sage: Phi = phi_stabilized.lift(p,prec,None,algorithm='stevens',eigensymbol=True)
            sage: L = pAdicLseries(Phi)
            sage: L.quadratic_twist()
            1
        """
        return self._quadratic_twist

    def _repr_(self):
        r"""
        Returns the print representation

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('37a')
            sage: p = 5
            sage: prec = 6
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec)
            sage: Phi = phi_stabilized.lift(p,prec,None,algorithm='stevens',eigensymbol=True)
            sage: L = pAdicLseries(Phi)
            sage: L._repr_()
            '5-adic L-series of Modular symbol of level 37 with values in Space of 5-adic distributions with k=0 action and precision cap 7'
        """
        s = "%s-adic L-series of %s"%(self.prime(), self.symb())
        return s

    def series(self, n, prec):
        r"""
        Returns the `n`-th approximation to the `p`-adic `L`-series
        associated to self, as a power series in `T` (corresponding to
        `\gamma-1` with `\gamma= 1 + p` as a generator of `1+p\ZZ_p`).

        EXAMPLES::
        
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('57a')
            sage: p = 5
            sage: prec = 4
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec+3)
            sage: Phi = phi_stabilized.lift(p,prec,None,algorithm='stevens',eigensymbol=True)
            sage: L = pAdicLseries(Phi)
            sage: L.series(3,4)
            O(5^3) + (3*5 + 5^2 + O(5^3))*T + (5 + O(5^2))*T^2
            
            sage: L1 = E.padic_lseries(5)
            sage: L1.series(4)
            O(5^6) + (3*5 + 5^2 + O(5^3))*T + (5 + 4*5^2 + O(5^3))*T^2 + (4*5^2 + O(5^3))*T^3 + (2*5 + 4*5^2 + O(5^3))*T^4 + O(T^5)
        
        """
        p = self.prime()
        M = self.symb().precision_absolute()
        K = pAdicField(p, M)
        R = PowerSeriesRing(K, names = 'T')
        T = R.gens()[0]
        R.set_default_prec(prec)
        return sum(self[i] * T**i for i in range(n))

    def interpolation_factor(self, ap,chip=1, psi = None):
        r"""
        Returns the interpolation factor associated to self

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('57a')
            sage: p = 5
            sage: prec = 4
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec)
            sage: Phi = phi_stabilized.lift(p,prec,None,algorithm='stevens')
            sage: L = pAdicLseries(Phi)
            sage: ap = phi.Tq_eigenvalue(p)
            sage: L.interpolation_factor(ap)
            4 + 2*5 + 4*5^3 + O(5^4)

        Comparing against a different implementation:

            sage: L = E.padic_lseries(5)
            sage: (1-1/L.alpha(prec=4))^2
            4 + 2*5 + 4*5^3 + O(5^4)

        """
        M = self.symb().precision_absolute()
        p = self.prime()
        if p == 2:
            R = pAdicField(2, M + 1)
        else:
            R = pAdicField(p, M)
        if psi != None:
            ap = psi(ap)
        ap = ap*chip
        sdisc = R(ap**2 - 4*p).sqrt()
        v0 = (R(ap) + sdisc) / 2
        v1 = (R(ap) - sdisc) / 2
        if v0.valuation() > 0:
            v0, v1 = v1, v0
        alpha = v0
        return (1 - 1/alpha)**2
    
    def eval_twisted_symbol_on_Da(self, a): # rename! should this be in modsym?
        """
        Returns `\Phi_{\chi}(\{a/p}-{\infty})` where `Phi` is the OMS and
        `\chi` is a the quadratic character corresponding to self

        INPUT:
            - ``a`` -- integer in range(p)

        OUTPUT:

        The distribution `\Phi_{\chi}(\{a/p\}-\{\infty\})`.

        EXAMPLES:
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('57a')
            sage: p = 5
            sage: prec = 4
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: ap = phi.Tq_eigenvalue(p,prec)
            sage: Phi = phi.p_stabilize_and_lift(p,ap = ap, M = prec, algorithm='stevens')
            sage: L = pAdicLseries(Phi)
            sage: L.eval_twisted_symbol_on_Da(1)
            (2 + 2*5 + 2*5^2 + 2*5^3 + O(5^4), 2 + 3*5 + 2*5^2 + O(5^3), 4*5 + O(5^2), 3 + O(5))

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('40a4')
            sage: p = 7
            sage: prec = 4
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: ap = phi.Tq_eigenvalue(p,prec)
            sage: Phi = phi.p_stabilize_and_lift(p,ap = ap, M = prec, algorithm='stevens')
            sage: L = pAdicLseries(Phi)
            sage: L.eval_twisted_symbol_on_Da(1)
            (4 + 6*7 + 3*7^2 + O(7^4), 2 + 7 + O(7^3), 4 + 6*7 + O(7^2), 6 + O(7))

        """
        symb = self.symb()
        p = symb.parent().prime()
        S0p = Sigma0(p)
        Dists = symb.parent().coefficient_module()
        M = Dists.precision_cap()
        p = Dists.prime()
        twisted_dist = Dists.zero_element()
        m_map = symb._map
        D = self._quadratic_twist
        for b in range(1, abs(D) + 1):
            if gcd(b, D) == 1:
                M1 = S0p([1, (b / abs(D)) % p**M, 0, 1])
                new_dist = m_map(M1 * M2Z([a, 1, p, 0]))*M1
                new_dist = new_dist.scale(kronecker(D, b)).normalize()
                twisted_dist = twisted_dist + new_dist
                #ans = ans + self.eval(M1 * M2Z[a, 1, p, 0])._right_action(M1)._lmul_(kronecker(D, b)).normalize()
        return twisted_dist.normalize()

    def _basic_integral(self, a, j):
        r"""
        Returns `\int_{a+pZ_p} (z-{a})^j d\Phi(0-infty)`
        -- see formula [Pollack-Stevens, sec 9.2]

        INPUT:

        - ``a`` -- integer in range(p)
        - ``j`` -- integer in range(self.symb().precision_absolute())

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_elliptic_curve
            sage: E = EllipticCurve('57a')
            sage: p = 5
            sage: prec = 4
            sage: phi = ps_modsym_from_elliptic_curve(E)
            sage: phi_stabilized = phi.p_stabilize(p,M = prec+3)
            sage: Phi = phi_stabilized.lift(p,prec,None,algorithm = 'stevens',eigensymbol = True)
            sage: L = pAdicLseries(Phi)
            sage: L.eval_twisted_symbol_on_Da(1)
            (2 + 2*5 + 2*5^2 + 2*5^3 + O(5^4), 2 + 3*5 + 2*5^2 + O(5^3), 4*5 + O(5^2), 3 + O(5))
            sage: L._basic_integral(1,2)
            2*5^3 + O(5^4)
        
        """
        symb = self.symb()
        M = symb.precision_absolute()
        if j > M:
            raise PrecisionError ("Too many moments requested")
        p = self.prime()
        ap = symb.Tq_eigenvalue(p)
        D = self._quadratic_twist
        ap = ap * kronecker(D, p)
        K = pAdicField(p, M)
        symb_twisted = self.eval_twisted_symbol_on_Da(a)
        return sum(binomial(j, r) * ((a - ZZ(K.teichmuller(a)))**(j - r)) *
                (p**r) * symb_twisted.moment(r) for r in range(j + 1)) / ap

def log_gamma_binomial(p,gamma,z,n,M):
    r"""
    Returns the list of coefficients in the power series
    expansion (up to precision `M`) of `{\log_p(z)/\log_p(\gamma) \choose n}`

    INPUT:

        - ``p`` --  prime
        - ``gamma`` -- topological generator e.g., `1+p`
        - ``z`` -- variable
        - ``n`` -- nonnegative integer
        - ``M`` -- precision

    OUTPUT:

    The list of coefficients in the power series expansion of
    `{\log_p(z)/\log_p(\gamma) \choose n}`

    EXAMPLES:

        sage: R.<z> = QQ['z']
        sage: from sage.modular.pollack_stevens.padic_lseries import log_gamma_binomial
        sage: log_gamma_binomial(5,1+5,z,2,4)
        [0, -3/205, 651/84050, -223/42025]
        sage: log_gamma_binomial(5,1+5,z,3,4)
        [0, 2/205, -223/42025, 95228/25845375]
    """
    L = sum([ZZ(-1)**j / j*z**j for j in range (1,M)]) #log_p(1+z)
    loggam = L / (L(gamma - 1))                  #log_{gamma}(1+z)= log_p(1+z)/log_p(gamma)
    return z.parent()(binomial(loggam,n)).truncate(M).list()

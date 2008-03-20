"""
Elliptic curves over number fields

EXAMPLES:
    sage: k.<i> = NumberField(x^2+1)
    sage: E = EllipticCurve([i,2])
    sage: E.j_invariant()
    -23328/365*i + 864/365
    sage: E.simon_two_descent()
    (1, 1, [(2*i : -2*i + 2 : 1)])
    sage: P = E([2*i,-2*i+2])
    sage: P+P
    (15/32*i + 3/4 : 139/256*i + 339/256 : 1)

"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                          William Stein   <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from ell_field import EllipticCurve_field
import ell_point
from sage.rings.ring import Ring
from sage.rings.arith import lcm
from sage.misc.misc import prod
import sage.databases.cremona

from gp_simon import simon_two_descent
from constructor import EllipticCurve
from sage.rings.all import PolynomialRing, QQ, ZZ, is_Ideal, is_NumberFieldElement, is_NumberFieldFractionalIdeal
from sage.misc.misc import verbose, forall
from sage.misc.functional import ideal
from kodaira_symbol import KodairaSymbol
from sage.rings.integer import Integer
from sage.structure.element import RingElement
from sage.rings.infinity import Infinity # just for verbose output

class EllipticCurve_number_field(EllipticCurve_field):
    """
    Elliptic curve over a number field.

    EXAMPLES:
        sage: K.<i>=NumberField(x^2+1)
        sage: EllipticCurve([i, i - 1, i + 1, 24*i + 15, 14*i + 35])
        Elliptic Curve defined by y^2 + i*x*y + (i+1)*y = x^3 + (i-1)*x^2 + (24*i+15)*x + (14*i+35) over Number Field in i with defining polynomial x^2 + 1
    """
    def __init__(self, x, y=None):
        """
        Allow some ways to create an ellitpic curve over a number
        field in addition to the generic ones:

        EXAMPLES:

        A curve from the database of curves over Q, but over a larger field:

            sage: K.<i>=NumberField(x^2+1)
            sage: EllipticCurve(K,'389a1')
            Elliptic Curve defined by y^2 + y = x^3 + x^2 + (-2)*x over Number Field in i with defining polynomial x^2 + 1

        Making the field of definition explicitly larger:

            sage: EllipticCurve(K,[0,-1,1,0,0])
            Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 over Number Field in i with defining polynomial x^2 + 1

        """
        if y is None:
            if isinstance(x, list):
                ainvs = x
                field = ainvs[0].parent()
        else:
            if isinstance(y, str):
                field = x
                X = sage.databases.cremona.CremonaDatabase()[y]
                ainvs = X.a_invariants()
            else:
                field = x
                ainvs = y
        if not (isinstance(field, Ring) and isinstance(ainvs,list)):
            raise TypeError

        EllipticCurve_field.__init__(self, [field(x) for x in ainvs])
        self._point_class = ell_point.EllipticCurvePoint_field

    def simon_two_descent(self, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Computes (probably) the rank of the Mordell-Weil group,
        with certainty the rank of the 2-Selmer group, and a list
        of independent points on the Weierstrass model of self.

        If the curve has 2-torsion, only the probable rank is returned.

        INPUT:
            verbose -- integer, 0,1,2,3; (default: 0), the verbosity level
            lim1    -- (default: 5) limite des points triviaux sur les quartiques
            lim3    -- (default: 50) limite des points sur les quartiques ELS
            limtriv -- (default: 10) limite des points triviaux sur la
                                     courbe elliptique
            maxprob -- (default: 20)
            limbigprime -- (default: 30)  to distinguish between small and large prime
                                          numbers. Use probabilistic tests for large
                                          primes. If 0, don't use probabilistic tests.

        OUTPUT:
            integer -- "probably" the rank of self
            integer -- the 2-rank of the Selmer group
            list    -- list of independent points on the Weierstrass model

        NOTE: For non-quadratic number fields, this code does return, but it takes a long time.

        IMPLEMENTATION: Uses {\bf Denis Simon's} GP/PARI scripts from
                         \url{http://www.math.unicaen.fr/~simon/}

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: E = EllipticCurve(K, '37')
            sage: E.simon_two_descent()
            (2, 2, [(-1 : 0 : 1), (1/2*a - 5/2 : -1/2*a - 13/2 : 1)])

            sage: K.<a> = NumberField(x^2 + 7, 'a')
            sage: E = EllipticCurve(K, [0,0,0,1,a]); E
            Elliptic Curve defined by y^2  = x^3 + x + a over Number Field in a with defining polynomial x^2 + 7
            sage: v = E.simon_two_descent(verbose=1); v
            courbe elliptique : Y^2 = x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7)
            A = 0
            B = Mod(1, y^2 + 7)
            C = Mod(y, y^2 + 7)
            LS2gen = [Mod(Mod(-5, y^2 + 7)*x^2 + Mod(-3*y, y^2 + 7)*x + Mod(8, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7)), Mod(Mod(1, y^2 + 7)*x^2 + Mod(1/2*y + 1/2, y^2 + 7)*x - 1, x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))]
            #LS2gen = 2
             Recherche de points triviaux sur la courbe
            points triviaux sur la courbe = [[1, 1, 0], [Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            zc = Mod(Mod(-5, y^2 + 7)*x^2 + Mod(-3*y, y^2 + 7)*x + Mod(8, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
             symbole de Hilbert (Mod(2, y^2 + 7),Mod(-5, y^2 + 7)) = -1
            zc = Mod(Mod(1, y^2 + 7)*x^2 + Mod(1/2*y + 1/2, y^2 + 7)*x + Mod(-1, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
             vient du point trivial [Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]
            m1 = 1
            m2 = 1
            #S(E/K)[2]    = 2
            #E(K)/2E(K)   = 2
            #III(E/K)[2]  = 1
            rang(E/K)     = 1
            listpointsmwr = [[Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            (1, 1, [(1/2*a + 3/2 : -a - 2 : 1)])

        A curve with 2-torsion
            sage: K.<a> = NumberField(x^2 + 7, 'a')
            sage: E = EllipticCurve(K, '15a')
            sage: v = E.simon_two_descent(); v  # long time (about 10 seconds), points can vary
            (1, 3, [...])
        """
        x = PolynomialRing(self.base_ring(), 'x').gen(0)
        t = simon_two_descent(self,
                              verbose=verbose, lim1=lim1, lim3=lim3, limtriv=limtriv,
                              maxprob=maxprob, limbigprime=limbigprime)
        prob_rank = Integer(t[0])
        two_selmer_rank = Integer(t[1])
        prob_gens = [self(P) for P in t[2]]
        return prob_rank, two_selmer_rank, prob_gens

    def is_local_integral_model(self,*P):
        r"""
        Tests if self is integral at the prime ideal $P$, or at all the
        primes if P is a list or tuple

        EXAMPLES:
            sage: K.<i>=NumberField(x^2+1)
            sage: P1,P2 =(K.factor_integer(5)[j][0] for j in (0,1))
            sage: E=EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: E.is_local_integral_model(P1,P2)
            False
            sage: Emin=E.local_integral_model(P1,P2)
            sage: Emin.is_local_integral_model(P1,P2)
            True
        """
        if len(P)==1: P=P[0]
        if isinstance(P,(tuple,list)):
            return forall(P, lambda x : self.is_local_integral_model(x))[0]
        return forall(self.ainvs(), lambda x : x.valuation(P) >= 0)[0]

    def local_integral_model(self,*P):
        r"""
        Return a model of self which is integral at the prime ideal $P$
        NB Does not affect integrality at other primes even if P non-principal

        EXAMPLES:
            sage: K.<i>=NumberField(x^2+1)
            sage: P1,P2 =(K.factor_integer(5)[j][0] for j in (0,1))
            sage: E=EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: E.local_integral_model((P1,P2))
            Elliptic Curve defined by y^2 + (-i)*x*y + (-25*i)*y = x^3 + 5*i*x^2 + 125*i*x + 3125*i over Number Field in i with defining polynomial x^2 + 1
        """
        if len(P)==1: P=P[0]
        if isinstance(P,(tuple,list)):
            E=self
            for Pi in P: E=E.local_integral_model(Pi)
            return E
        ai = self.a_invariants()
        e  = min([(ai[i].valuation(P)/[1,2,3,4,6][i]) for i in range(5)]).floor()
        pi = self.base_field().uniformizer(P, 'negative')
        return EllipticCurve([ai[i]/pi**(e*[1,2,3,4,6][i]) for i in range(5)])

    def is_global_integral_model(self):
        r"""
        Return true iff self is integral at all primes

        EXAMPLES:
            sage: K.<i>=NumberField(x^2+1)
            sage: E=EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1,P2 = (K.factor_integer(5)[j][0] for j in (0,1))
            sage: Emin=E.global_integral_model()
            sage: Emin.is_global_integral_model()
            True
        """
        return forall(self.a_invariants(), lambda x : x.is_integral())[0]

    def global_integral_model(self):
        r"""
        Return a model of self which is integral at all primes

        EXAMPLES:
            sage: K.<i>=NumberField(x^2+1)
            sage: E=EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1,P2 = (K.factor_integer(5)[j][0] for j in (0,1))
            sage: E.global_integral_model()
            Elliptic Curve defined by y^2 + (-i)*x*y + (-25*i)*y = x^3 + 5*i*x^2 + 125*i*x + 3125*i over Number Field in i with defining polynomial x^2 + 1

        """
        K = self.base_field()
        ai = self.a_invariants()
        for a in ai:
            if not a.is_integral():
               for P, _ in K.ideal(a.denominator()).factor():
                   pi=K.uniformizer(P,'negative')
                   e  = min([(ai[i].valuation(P)/[1,2,3,4,6][i]) for i in range(5)]).floor()
                   ai = [ai[i]/pi**(e*[1,2,3,4,6][i]) for i in range(5)]
        for z in ai:
            assert z.denominator() == 1, "bug in global_integral_model: %s" % ai
        return EllipticCurve(ai)

    def integral_model(self):
        r"""
        Return a weierstrass model, $F$, of self with integral coefficients,
        along with a morphism $\phi$ of points on self to points on $F$.

        EXAMPLES:
            sage: E = EllipticCurve([1/2,0,0,5,1/3])
            sage: F, phi = E.integral_model()
            sage: F
            Elliptic Curve defined by y^2 + 3*x*y  = x^3 + 6480*x + 15552 over Rational Field
            sage: phi
            Generic morphism:
              From: Abelian group of points on Elliptic Curve defined by y^2 + 1/2*x*y  = x^3 + 5*x + 1/3 over Rational Field
              To:   Abelian group of points on Elliptic Curve defined by y^2 + 3*x*y  = x^3 + 6480*x + 15552 over Rational Field
              Via:  (u,r,s,t) = (1/6, 0, 0, 0)
            sage: P = E([4/9,41/27])
            sage: phi(P)
            (16 : 328 : 1)
            sage: phi(P) in F
            True
        """
        F = self.global_integral_model()
        return F, self.isomorphism_to(F)

    def _tidy_model(self):
        """
        Transforms the elliptic curve to a model in which a1, a2, a3
        are reduced modulo 2, 3, 2 respectively.

        This only works on integral models, ie it requires that a1, a2
        and a3 lie in the ring of integers of the base field.

        EXAMPLES:
            sage: K.<a>=NumberField(x^2-38)

            sage: E=EllipticCurve([a, -5*a + 19, -39*a + 237, 368258520200522046806318224*a - 2270097978636731786720858047, 8456608930180227786550494643437985949781*a - 52130038506835491453281450568107193773505])
            sage: E.ainvs()
            [a,
            -5*a + 19,
            -39*a + 237,
            368258520200522046806318224*a - 2270097978636731786720858047,
            8456608930180227786550494643437985949781*a - 52130038506835491453281450568107193773505]

            sage: E._tidy_model().ainvs()
            [-a,
            a + 1,
            -11*a + 151,
            368258520200522046806318520*a - 2270097978636731786720859535,
            8456608930173478039472018047583706317255*a - 52130038506793883217874390501829588398139]
        """
        ZK = self.base_ring().maximal_order()
        (a1, a2, a3, a4, a6) = [ZK(a) for a in self.a_invariants()]
        # N.B. Must define s, r, t in the right order.
        s = -ZK([(a/2).round('away') for a in a1.list()])
        r = -ZK([(a/3).round('away') for a in (a2 + s*a1 +s*s).list()])
        t = -ZK([(a/2).round('away') for a in (a3 - r*a1).list()])

        return self.rst_transform(r, s, t)

    def local_information(self, P=None, proof = None):
        """
        Tate's algorithm for an elliptic curve over a number field.

        If a prime P of the base field is specified, computes local
        reduction data at the prime ideal P and a local minimal model.
        If no P is specified, computes local information at all bad primes.

        The model is not required to be integral on input.
        If P is principal, uses a generator as uniformizer, so it will not affect
        integrality or minimality at other primes.
        If P is not principal, the minimal model returned will preserve integrality
        at other primes, but not minimality.

        INPUT:
            self -- an elliptic curve over a number field.
            P    -- either None or a prime ideal of the base field of self.
            proof -- whether to only use provably correct methods (default controled by
                     global proof module).  Note that the proof module is number_field,
                     not elliptic_curves, since the functions that actually need the flag
                     are in number fields.

        OUTPUT:
            If P specified, returns a 6-tuple with the following data:
              Emin -- a model (integral and) minimal at P
              p    -- the residue characteristic
              vpd  -- the valuation of the local minimal discriminant
              fp   -- valuation of the conductor
              KS   -- Kodaira symbol
              cp   -- Tamagawa number
            Otherwise, for all primes dividing the discriminant, returns a pair with the first
            member of the pair being that prime P, and the second being a tuple with the above
            data for that P.

        EXAMPLES
            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([1 + i  ,0  ,1  ,0  ,0  ])
            sage: E.local_information()
            [(Fractional ideal (-3*i - 2),
            (Elliptic Curve defined by y^2 + (i+1)*x*y + (4*i+7)*y = x^3 + 12*x^2 + (-i+47)*x + (-4*i+58) over Number Field in i with defining polynomial x^2 + 1, 13, 2, 1, I2, 2)),
            (Fractional ideal (2*i + 1),
            (Elliptic Curve defined by y^2 + (i+1)*x*y + (4*i+7)*y = x^3 + 12*x^2 + (-i+47)*x + (-4*i+58) over Number Field in i with defining polynomial x^2 + 1, 5, 1, 1, I1, 1))]
            sage: E.local_information(K.ideal(3))
            (Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1,
            3, 0, 0, I0, 1)
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")
        if P is None:
            primes = [f[0] for f in self.base_ring().ideal(self.discriminant()).factor()]
            return [(pr, self._tate(pr, proof)) for pr in primes]
        if not (is_NumberFieldFractionalIdeal(P) and P.is_prime() # and P.order() == self.base_ring().integers()
                or is_NumberFieldElement(P) and P.parent().ideal(P).is_prime()
                or self.base_ring() is QQ and (isinstance(P, Integer) and P.is_prime()
                                               or is_Ideal(P) and P.base_ring() is ZZ and P.is_prime())):
            raise TypeError, "second argument must be a prime ideal"
        if isinstance(P, RingElement):
            P = self.base_ring().ideal(P)
        return self.integral_model()[0]._tate(P, proof)

    def local_minimal_model(self, P, proof = None):
        """
        Returns a model which is integral at all primes and minimal at P

        The model is not required to be integral on input.
        If P is principal, uses a generator as uniformizer, so it will not affect
        integrality or minimality at other primes.
        If P is not principal, the minimal model returned will preserve integrality
        at other primes, but not minimality.

        INPUT:
            self -- an elliptic curve over a number field.
            P    -- a prime ideal of the base field of self.
            proof -- whether to only use provably correct methods (default controled by
                     global proof module).  Note that the proof module is number_field,
                     not elliptic_curves, since the functions that actually need the flag
                     are in number fields.

        OUTPUT:
            Emin -- a model (integral and) minimal at P

        EXAMPLES:
            sage: K.<a>=NumberField(x^2-5)
            sage: E=EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
            sage: P=K.ideal(a)
            sage: E.local_minimal_model(P).ainvs()
            [20, -87, 30, a - 277, 2*a - 213]
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")
        if not (is_NumberFieldFractionalIdeal(P) and P.is_prime() # and P.order() == self.base_ring().integers()
                or is_NumberFieldElement(P) and P.parent().ideal(P).is_prime()
                or self.base_ring() is QQ and (isinstance(P, Integer) and P.is_prime()
                                               or is_Ideal(P) and P.base_ring() is ZZ and P.is_prime())):
            raise TypeError, "second argument must be a prime ideal"
        if isinstance(P, RingElement):
            P = self.base_ring().ideal(P)
        return self.local_information(P, proof)[0]

    def conductor(self):
        """
        Returns the conductor of this elliptic curve as a fractional ideal of the base field.

        OUTPUT:
            a fractional ideal

        EXAMPLES:
            sage: K.<i>=NumberField(x^2+1)
            sage: EllipticCurve([i, i - 1, i + 1, 24*i + 15, 14*i + 35]).conductor()
            Fractional ideal (21*i - 3)
            sage: K.<a>=NumberField(x^2-x+3)
            sage: EllipticCurve([1 + a  ,-1 + a  ,1 + a  ,-11 + a  ,5 -9*a  ]).conductor()
            Fractional ideal (-6*a)

            A not so well known curve with everywhere good reduction:
            sage: K.<a>=NumberField(x^2-38)
            sage: E=EllipticCurve([0,0,0, 21796814856932765568243810*a - 134364590724198567128296995, 121774567239345229314269094644186997594*a - 750668847495706904791115375024037711300])
            sage: E.conductor()
            Fractional ideal (1)
        """
        ## Ported from John Cremona's code implementing Tate's algorithm.
        primes = [f[0] for f in self.base_ring().ideal(self.discriminant()).factor()]
        ans = self.base_ring().ideal(1)
        for P in primes:
            ans *= P**(self._tate(P)[3])
        return ans

    def global_minimal_model(self, proof = None):
        """
        Returns a model of self that is minimal at all primes, and the conductor of self.

        Note that this only works for class number 1.

        INPUT:
            self -- an elliptic curve over a number field of class number
            proof -- whether to only use provably correct methods (default controled by
                     global proof module).  Note that the proof module is number_field,
                     not elliptic_curves, since the functions that actually need the flag
                     are in number fields.

        OUTPUT:
            A 2-tuple consisting of a global minimal model, and
            the conductor of self as a fractional ideal of the base
            field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2-38)
            sage: E = EllipticCurve([0,0,0, 21796814856932765568243810*a - 134364590724198567128296995, 121774567239345229314269094644186997594*a - 750668847495706904791115375024037711300])
            sage: E.global_minimal_model()
            (Elliptic Curve defined by y^2 + a*x*y + (-39*a+237)*y = x^3 + (-5*a+19)*x^2 + (368258520200522046806318224*a-2270097978636731786720858047)*x + (8456608930180227786550494643437985949781*a-52130038506835491453281450568107193773505) over Number Field in a with defining polynomial x^2 - 38, Fractional ideal (1))
        """
        ## Ported from John Cremona's code implementing Tate's algorithm.
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")
        K = self.base_ring()
        if K.class_number() != 1:
            raise ValueError, "global minimal models only exist in general for class number 1"
        primes = [f[0] for f in self.base_ring().ideal(self.discriminant()).factor()]
        N = K.ideal(1)
        E = self
        for P in primes:
            local_info = E._tate(P, proof)
            N *= P**local_info[3]
            E = local_info[0]
        return (E._tidy_model(), N)

    def _tate(self, P, proof = None):
        """
        Tate's algorithm for an elliptic curve over a number field:
        computes local reduction data at the prime ideal P and a local
        minimal model.

        The model is not required to be integral on input.  If P is
        principal, uses a generator as uniformizer, so it will not
        affect integrality or minimality at other primes.  If P is not
        principal, the minimal model returned will preserve
        integrality at other primes, but not minimality.

        INPUT:
            self -- an elliptic curve over a number field.
            P    -- a prime ideal of the base field of self.

        OUTPUT:
            Emin -- a model (integral and) minimal at P
            p    -- the residue characteristic
            vpd  -- the valuation of the local minimal discriminant
            fp   -- valuation of the conductor
            KS   -- Kodaira symbol
            cp   -- Tamagawa number

        EXAMPLES: see self.local_information()
        """
        ## Ported from John Cremona's code implementing Tate's algorithm.
        K = self.base_ring()
        OK = K.maximal_order()
        t = verbose("Running Tate's algorithm with P = %s"%P, level=1)
        F = OK.residue_field(P)
        p = F.characteristic()
        if P.is_principal():
            pi = P.gens_reduced()[0]
            verbose("P is principal, generator pi = %s"%pi, t, 1)
        else:
            pi = K.uniformizer(P, 'negative')
            verbose("P is not principal, uniformizer pi = %s"%pi, t, 1)

        def _pval(x):
            """
            Local function returning the valuation of x at P
            """
            if x==0: return Infinity
            return K.ideal(x).valuation(P)
        def _pdiv(x):
            """
            Local function returning True iff P divides x
            """
            return x==0 or _pval(x) > 0
        def _pinv(x):
            """
            Local function returning an inverse of x mod P
            """
            return F.lift(~F(x))
        def _proot(x, e):
            """
            Local function returning an e'th root of x mod P
            """
            L = F(x).nth_root(e, extend = False, all = True)
            assert len(L) > 0, "no e'th root exists mod P"
            return F.lift(L[0])
        def _preduce(x):
            """
            Local function returning x reduced modulo P
            """
            return F.lift(F(x))
        def _pquadroots(a, b, c):
            r"""
            Local function returning True iff $ax^2 + bx + c$ has roots modulo P
            """
            (a, b, c) = (F(a), F(b), F(c))
            if a == 0:
                return (b != 0) or (c == 0)
            elif p == 2:
                return len(PolynomialRing(F, "x")([c,b,a]).roots()) > 0
            else:
                return (b**2 - 4*a*c).is_square()
        def _pcubicroots(b, c, d):
            r"""
            Local function returning the number of roots of $x^3 +
            b*x^2 + c*x + d$ modulo P, counting multiplicities
            """
            return sum([rr[1] for rr in PolynomialRing(F, 'x')([d, c, b, 1]).roots()],0)

        if p == 2:
            halfmodp = OK(Integer(0))
        else:
            halfmodp = _pinv(Integer(2))

        A = self.a_invariants()
        A = [0, A[0], A[1], A[2], A[3], 0, A[4]]
        indices = [1,2,3,4,6]
        if min([_pval(a) for a in A if a != 0]) < 0:
            verbose("Non-integral model at P: valuations are %s; making integral"%([_pval(a) for a in A if a != 0]), t, 1)
            e = 0
            for i in range(7):
                if A[i] != 0:
                    e = max(e, (-_pval(A[i])/i).ceil())
            pie = pi**e
            for i in range(7):
                if A[i] != 0:
                    A[i] *= pie**i
            verbose("P-integral model is %s, with valuations %s"%([A[i] for i in indices], [_pval(A[i]) for i in indices]), t, 1)

        (a1, a2, a3, a4, a6) = (A[1], A[2], A[3], A[4], A[6])
        while True:
            C = EllipticCurve([a1, a2, a3, a4, a6]);
            (b2, b4, b6, b8) = C.b_invariants()
            (c4, c6) = C.c_invariants()
            delta = C.discriminant()
            vpd = _pval(delta)

            if vpd == 0:
                ## Good reduction already
                cp = 1
                fp = 0
                KS = KodairaSymbol("I0")
                break #return

            # Otherwise, we change coordinates so that p | a3, a4, a6
            if p == 2:
                if _pdiv(b2):
                    r = _proot(a4, 2)
                    t = _proot(((r + a2)*r + a4)*r + a6, 2)
                else:
                    temp = _pinv(a1)
                    r = temp * a3
                    t = temp * (a4 + r*r)
            elif p == 3:
                if _pdiv(b2):
                    r = _proot(-b6, 3)
                else:
                    r = -_pinv(b2) * b4
                t = a1 * r + a3
            else:
                if _pdiv(c4):
                    r = -_pinv(12) * b2
                else:
                    r = -_pinv(12*c4) * (c6 + b2 * c4)
                t = -halfmodp * (a1 * r + a3)
            r = _preduce(r)
            t = _preduce(t)
            # print "Before first tranform C = %s"%C
            # print "[a1,a2,a3,a4,a6] = %s"%([a1, a2, a3, a4, a6])
            C = C.rst_transform(r, 0, t)
            (a1, a2, a3, a4, a6) = C.a_invariants()
            (b2, b4, b6, b8) = C.b_invariants()
            if min([_pval(a) for a in (a1, a2, a3, a4, a6) if a != 0]) < 0:
                raise RuntimeError, "Non-integral model after first transform!"
            verbose("After first transform %s\n, [a1,a2,a3,a4,a6] = %s\n, valuations = %s"%([r, 0, t], [a1, a2, a3, a4, a6], [_pval(a1), _pval(a2), _pval(a3), _pval(a4), _pval(a6)]), t, 2)
            if _pval(a3) == 0:
                raise RuntimeError, "p does not divide a3 after first transform!"
            if _pval(a4) == 0:
                raise RuntimeError, "p does not divide a4 after first transform!"
            if _pval(a6) == 0:
                raise RuntimeError, "p does not divide a6 after first transform!"

            # Now we test for Types In, II, III, IV
            # Do we not have to update the c invariants?
            if not _pdiv(c4):
                ## Type In (n = vpd)
                if _pquadroots(1, a1, -a2):
                    cp = vpd
                elif Integer(2).divides(vpd):
                    cp = 2
                else:
                    cp = 1
                KS = KodairaSymbol("I%s"%vpd)
                fp = 1
                break #return
            if _pval(a6) < 2:
                ## Type II
                KS = KodairaSymbol("II")
                fp = vpd
                cp = 1
                break #return
            if _pval(b8) < 3:
                ## Type III
                KS = KodairaSymbol("III")
                fp = vpd - 1
                cp = 2
                break #return
            if _pval(b6) < 3:
                ## Type IV
                if _pquadroots(1, a3 / pi, -a6/(pi*pi)):
                    cp = 3
                else:
                    cp = 1
                KS = KodairaSymbol("IV")
                fp = vpd - 2
                break #return

            # If our curve is none of these types, we change types so that p | a1, a2 and p^2 | a3, a4 and p^3 | a6
            if p == 2:
                s = _proot(a2, 2)
                t = pi*_proot(a6/(pi*pi), 2)
            elif p == 3:
                s = a1
                t = a3
            else:
                s = -a1*halfmodp
                t = -a3*halfmodp
            C = C.rst_transform(0, s, t)
            (a1, a2, a3, a4, a6) = C.a_invariants()
            (b2, b4, b6, b8) = C.b_invariants()
            verbose("After second transform %s\n[a1, a2, a3, a4, a6] = %s\nValuations: %s"%([0, s, t], [a1,a2,a3,a4,a6],[_pval(a1),_pval(a2),_pval(a3),_pval(a4),_pval(a6)]), t, 2)
            if _pval(a1) == 0:
                raise RuntimeError, "p does not divide a1 after second transform!"
            if _pval(a2) == 0:
                raise RuntimeError, "p does not divide a2 after second transform!"
            if _pval(a3) < 2:
                raise RuntimeError, "p^2 does not divide a3 after second transform!"
            if _pval(a4) < 2:
                raise RuntimeError, "p^2 does not divide a4 after second transform!"
            if _pval(a6) < 3:
                raise RuntimeError, "p^3 does not divide a6 after second transform!"
            if min(_pval(a1), _pval(a2), _pval(a3), _pval(a4), _pval(a6)) < 0:
                raise RuntimeError, "Non-integral model after second transform!"

            # Analyze roots of the cubic T^3 + bT^2 + cT + d = 0, where b = a2/p, c = a4/p^2, d = a6/p^3
            b = a2/pi
            c = a4/(pi*pi)
            d = a6/(pi**3)
            bb = b*b
            cc = c*c
            bc = b*c
            w = 27*d*d - bb*cc + 4*b*bb*d - 18*bc*d + 4*c*cc
            x = 3*c - bb
            if _pdiv(w):
                if _pdiv(x):
                    sw = 3
                else:
                    sw = 2
            else:
                sw = 1
            verbose("Analyzing roots of cubic T^3 + %s*T^2 + %s*T + %s, case %s"%(b, c, d, sw), t, 1)
            if sw == 1:
                ## Three distinct roots - Type I*0
                verbose("Distinct roots", t, 1)
                KS = KodairaSymbol("I0*")
                cp = 1 + _pcubicroots(b, c, d)
                fp = vpd - 4
                break #return
            elif sw == 2:
                ## One double root - Type I*m for some m
                verbose("One double root", t, 1)
                ## Change coords so that the double root is T = 0 mod p
                if p == 2:
                    r = _proot(c, 2)
                elif p == 3:
                    r = c * _pinv(b)
                else:
                    r = (bc - 9*d)*_pinv(2*x)
                r = pi * _preduce(r)
                C = C.rst_transform(r, 0, 0)
                (a1, a2, a3, a4, a6) = C.a_invariants()
                (b2, b4, b6, b8) = C.b_invariants()
                ix = 3; iy = 3; mx = pi*pi; my = pi*pi
                while True:
                    a2t = a2 / pi
                    a3t = a3 / my
                    a4t = a4 / (pi*mx)
                    a6t = a6 / (mx*my)
                    if _pdiv(a3t*a3t + 4*a6t):
                        if p == 2:
                            t = my*_proot(a6t, 2)
                        else:
                            t = my*_preduce(-a3t*halfmodp)
                        C = C.rst_transform(0, 0, t)
                        (a1, a2, a3, a4, a6) = C.a_invariants()
                        (b2, b4, b6, b8) = C.b_invariants()
                        my = my*pi
                        iy += 1
                        a2t = a2/pi
                        a3t = a3/my
                        a4t = a4/(pi*mx)
                        a6t = a6/(mx*my)
                        if _pdiv(a4t*a4t - 4*a6t*a2t):
                            if p == 2:
                                r = mx*_proot(a6t*_pinv(a2t), 2)
                            else:
                                r = mx*_preduce(-a4t*_pinv(2*a2t))
                            C = C.rst_transform(r, 0, 0)
                            (a1, a2, a3, a4, a6) = C.a_invariants()
                            (b2, b4, b6, b8) = C.b_invariants()
                            mx = mx*pi
                            ix += 1 # and stay in loop
                        else:
                            if _pquadroots(a2t, a4t, a6t):
                                cp = 4
                            else:
                                cp = 2
                            break # exit loop
                    else:
                        if _pquadroots(1, a3t, -a6t):
                            cp = 4
                        else:
                            cp = 2
                        break
                KS = KodairaSymbol("I%s*"%(ix+iy-5))
                fp = vpd - ix - iy + 1
                break #return
            else: # sw == 3
                ## The cubic has a triple root
                verbose("Triple root", t, 1)
                ## First we change coordinates so that T = 0 mod p
                if p == 2:
                    r = b
                elif p == 3:
                    r = _proot(-d, 3)
                else:
                    r = -b * _pinv(3)
                r = pi*_preduce(r)
                C = C.rst_transform(r, 0, 0)
                (a1, a2, a3, a4, a6) = C.a_invariants()
                (b2, b4, b6, b8) = C.b_invariants()
                verbose("After third transform %s\n[a1,a2,a3,a4,a6] = %s\nValuations: %s"%([r,0,0],[a1,a2,a3,a4,a6],[_pval(ai) for ai in [a1,a2,a3,a4,a6]]), t, 2)
                if min(_pval(ai) for ai in [a1,a2,a3,a4,a6]) < 0:
                    raise RuntimeError, "Non-integral model after third transform!"
                if _pval(a2) < 2 or _pval(a4) < 3 or _pval(a6) < 4:
                    raise RuntimeError, "Cubic after transform does not have a triple root at 0"
                a3t = a3/(pi*pi)
                a6t = a6/(pi**4)
                # We test for Type IV*
                if not _pdiv(a3t*a3t + 4*a6t):
                    cp = 3 if _pquadroots(1, a3t, -a6t) else 1
                    KS = KodairaSymbol("IV*")
                    fp = vpd - 6
                    break #return
                # Now change coordinates so that p^3|a3, p^5|a6
                t =        -pi*pi*_proot(a6t, 2) if p==2 \
                      else  pi*pi*_preduce(-a3t*halfmodp)
                C = C.rst_transform(0, 0, t)
                (a1, a2, a3, a4, a6) = C.a_invariants()
                (b2, b4, b6, b8) = C.b_invariants()
                # We test for types III* and II*
                if _pval(a4) < 4:
                    ## Type III*
                    KS = KodairaSymbol("III*")
                    fp = vpd - 7
                    cp = 2
                    break #return
                if _pval(a6) < 6:
                    ## Type II*
                    KS = KodairaSymbol("II*")
                    fp = vpd - 8
                    cp = 1
                    break #return
                a1 /= pi
                a2 /= pi**2
                a3 /= pi**3
                a4 /= pi**4
                a6 /= pi**6
                verbose("Non-minimal equation, dividing out...\nNew model is %s"%([a1, a2, a3, a4, a6]), t, 1)
        return (C, p, vpd, fp, KS, cp)


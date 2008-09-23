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
from sage.rings.arith import lcm, gcd
from sage.misc.misc import prod
import sage.databases.cremona
import ell_torsion

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
        Allow some ways to create an elliptic curve over a number
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
        self._point_class = ell_point.EllipticCurvePoint_number_field

    def simon_two_descent(self, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Computes (probably) the rank of the Mordell-Weil group,
        with certainty the rank of the 2-Selmer group, and a list
        of independent points on the Weierstrass model of self.

        If the curve has 2-torsion, only the probable rank is returned.

        INPUT:
            verbose -- 0, 1, 2, or 3 (default: 0), the verbosity level
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
            LS2gen = [Mod(Mod(-5, y^2 + 7)*x^2 + Mod(-3*y, y^2 + 7)*x + Mod(8, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7)), Mod(Mod(1, y^2 + 7)*x^2 + Mod(1/2*y - 1/2, y^2 + 7)*x - 1, x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))]
            #LS2gen = 2
             Recherche de points triviaux sur la courbe
            points triviaux sur la courbe = [[1, 1, 0], [Mod(1/2*y + 3/2, y^2 + 7), Mod(-y - 2, y^2 + 7), 1]]
            zc = Mod(Mod(-5, y^2 + 7)*x^2 + Mod(-3*y, y^2 + 7)*x + Mod(8, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
             symbole de Hilbert (Mod(2, y^2 + 7),Mod(-5, y^2 + 7)) = -1
             zc = Mod(Mod(1, y^2 + 7)*x^2 + Mod(1/2*y - 1/2, y^2 + 7)*x + Mod(-1, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
             symbole de Hilbert (Mod(-2*y + 2, y^2 + 7),Mod(1, y^2 + 7)) = 0
             sol de Legendre = [1, 0, 1]~
             zc*z1^2 = Mod(Mod(2*y - 2, y^2 + 7)*x + Mod(2*y + 10, y^2 + 7), x^3 + Mod(1, y^2 + 7)*x + Mod(y, y^2 + 7))
             quartique : (-1/2*y + 1/2)*Y^2 = x^4 + (-3*y - 15)*x^2 + (-8*y - 16)*x + (-11/2*y - 15/2)
             reduite: Y^2 = (-1/2*y + 1/2)*x^4 - 4*x^3 + (-3*y + 3)*x^2 + (2*y - 2)*x + (1/2*y + 3/2)
             non ELS en [2, [0, 1]~, 1, 1, [1, 1]~]
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
            sage: K.<i> = NumberField(x^2+1)
            sage: P1,P2 = K.primes_above(5)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: E.is_local_integral_model(P1,P2)
            False
            sage: Emin = E.local_integral_model(P1,P2)
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
            sage: K.<i> = NumberField(x^2+1)
            sage: P1,P2 = K.primes_above(5)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
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
            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1,P2 = K.primes_above(5)
            sage: Emin = E.global_integral_model()
            sage: Emin.is_global_integral_model()
            True
        """
        return forall(self.a_invariants(), lambda x : x.is_integral())[0]

    def global_integral_model(self):
        r"""
        Return a model of self which is integral at all primes

        EXAMPLES:
            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([i/5,i/5,i/5,i/5,i/5])
            sage: P1,P2 = K.primes_above(5)
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

    integral_model = global_integral_model

    def _tidy_model(self):
        """
        Transforms the elliptic curve to a model in which a1, a2, a3
        are reduced modulo 2, 3, 2 respectively.

        This only works on integral models, i.e. it requires that a1, a2
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
            [a,
            a + 1,
            a + 1,
            368258520200522046806318444*a - 2270097978636731786720859345,
            8456608930173478039472018047583706316424*a - 52130038506793883217874390501829588391299]
            sage: EllipticCurve([101,202,303,404,505])._tidy_model().ainvs()
            [1, 1, 0, -2509254, 1528863051]
            sage: EllipticCurve([-101,-202,-303,-404,-505])._tidy_model().ainvs()
            [1, -1, 0, -1823195, 947995262]
        """
        ZK = self.base_ring().maximal_order()
        try:
            (a1, a2, a3, a4, a6) = [ZK(a) for a in self.a_invariants()]
        except TypeError:
            raise TypeError, "_tidy_model() requires an integral model."
        # N.B. Must define s, r, t in the right order.
        if ZK.degree() == 1:
            s = ((-a1)/2).round('up')
            r = ((-a2 + s*a1 +s*s)/3).round()
            t = ((-a3 - r*a1)/2).round('up')
        else:
            s = ZK([(a/2).round('up') for a in (-a1).list()])
            r = ZK([(a/3).round() for a in (-a2 + s*a1 +s*s).list()])
            t = ZK([(a/2).round('up') for a in (-a3 - r*a1).list()])

        return self.rst_transform(r, s, t)

    def local_information(self, P=None, proof=None):
        r"""
        \code{local_information} has been renamed \code{local_data}
        and is being deprecated.
        """
        raise DeprecationWarning, "local_information is deprecated; use local_data instead"
        return self.local_data(P,proof)

    def local_data(self, P=None, proof = None):
        r"""
        Local data for this elliptic curve at a prime.

        If a prime $P$ of the base field is specified, computes local
        reduction data at the prime ideal $P$ and a local minimal model.
        If no $P$ is specified, computes local information at all primes
        in the support of the discriminant of this model.

        The model is not required to be integral on input.

        If $P$ is principal, uses a generator as uniformizer, so it
        will not affect integrality or minimality at other primes.  If
        $P$ is not principal, the minimal model returned will preserve
        integrality at other primes, but not minimality.

        INPUT:
            self -- an elliptic curve over a number field.
            $P$    -- either None or a prime ideal of the base field of self.
            proof -- whether to only use provably correct methods
                     (default controled by global proof module).  Note
                     that the proof module is number_field, not
                     elliptic_curves, since the functions that
                     actually need the flag are in number fields.

        OUTPUT:

            If $P$ specified, returns the EllipticCurveLocalData object
            associated to the prime P for this curve.  Otherwise,
            returns a list of such objects, one for each prime $P$
            in the support of the discriminant.

        EXAMPLES
            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve([1 + i  ,0  ,1  ,0  ,0  ])
            sage: E.local_data()
            [Local data at Fractional ideal (-3*i - 2) of Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1:
            Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1
            Minimal discriminant valuation: 2
            Conductor exponent: 1
            Kodaira Symbol: I2
            Tamagawa Number: 2, Local data at Fractional ideal (2*i + 1) of Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1:
            Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1
            Minimal discriminant valuation: 1
            Conductor exponent: 1
            Kodaira Symbol: I1
            Tamagawa Number: 1]
            sage: E.local_data(K.ideal(3))
            Local data at Fractional ideal (3) of Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1:
            Local minimal model: Elliptic Curve defined by y^2 + (i+1)*x*y + y = x^3 over Number Field in i with defining polynomial x^2 + 1
            Minimal discriminant valuation: 0
            Conductor exponent: 0
            Kodaira Symbol: I0
            Tamagawa Number: 1

        An example raised in \#3897:
            sage: E = EllipticCurve([1,1])
            sage: E.local_data(3)
            Local data at Principal ideal (3) of Integer Ring of Elliptic Curve defined by y^2  = x^3 + x +1 over Rational Field:
            Local minimal model: Elliptic Curve defined by y^2  = x^3 + x +1 over Rational Field
            Minimal discriminant valuation: 0
            Conductor exponent: 0
            Kodaira Symbol: I0
            Tamagawa Number: 1
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        if P is None:
            primes = self.base_ring()(self.discriminant()).support()
            return [self._get_local_data(pr, proof) for pr in primes]

        from sage.schemes.elliptic_curves.ell_local_data import check_prime
        P = check_prime(self.base_field(),P)

        return self._get_local_data(P,proof)

    def _get_local_data(self, P, proof):
        """
        Internal function to create data for this elliptic curve at a prime.

        This function handles the caching of local data.  It is called
        by local_data() which is the user interface and which parses
        the input parameters P and proof.
        """
        try:
            return self._local_data[P]
        except AttributeError:
            self._local_data = {}
        except KeyError:
            pass
        from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
        self._local_data[P] = EllipticCurveLocalData(self, P, proof)
        return self._local_data[P]

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
            [0, 1, 0, a - 33, -2*a + 64]
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        return self.local_data(P, proof).minimal_model()

    def tamagawa_number(self, P, proof = None):
        """
        Returns the Tamagawa number of this elliptic curve at the prime P.

        INPUT:
            self -- an elliptic curve over a number field.
            P -- a prime ideal of the base field of self, or a field
                 element generating such an ideal.
            proof -- whether to only use provably correct methods
                     (default controled by global proof module).  Note
                     that the proof module is number_field, not
                     elliptic_curves, since the functions that
                     actually need the flag are in number fields.

        OUTPUT:
            c -- (positive integer) the Tamagawa number of the curve at P.

        EXAMPLES:
            sage: K.<a>=NumberField(x^2-5)
            sage: E=EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
            sage: [E.tamagawa_number(P) for P in E.discriminant().support()]
            [1, 1, 1, 1]
            sage: K.<a> = QuadraticField(-11)
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: [E.tamagawa_number(P) for P in K(11).support()]
            [10]
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        return self.local_data(P, proof).tamagawa_number()

    def kodaira_symbol(self, P, proof = None):
        """
        Returns the Kodaira Symbol of this elliptic curve at the prime P.

        INPUT:
            self -- an elliptic curve over a number field.
            P -- a prime ideal of the base field of self, or a field
                 element generating such an ideal.
            proof -- whether to only use provably correct methods
                     (default controled by global proof module).  Note
                     that the proof module is number_field, not
                     elliptic_curves, since the functions that
                     actually need the flag are in number fields.

        OUTPUT:
            K -- (string) the Kodaira Symbol of the curve at P.

        EXAMPLES:
            sage: K.<a>=NumberField(x^2-5)
            sage: E=EllipticCurve([20, 225, 750, 625*a + 6875, 31250*a + 46875])
            sage: bad_primes = E.discriminant().support(); bad_primes
            [Fractional ideal (7/2*a - 81/2),
            Fractional ideal (a + 52),
            Fractional ideal (-a),
            Fractional ideal (2)]
            sage: [E.kodaira_symbol(P) for P in bad_primes]
            [I1, I1, I0, II]
            sage: K.<a> = QuadraticField(-11)
            sage: E = EllipticCurve('11a1').change_ring(K)
            sage: [E.kodaira_symbol(P) for P in K(11).support()]
            [I10]
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")

        return self.local_data(P, proof).kodaira_symbol()


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
        try:
            return self._conductor
        except AttributeError:
            pass

        # Note: for number fields other than QQ we could initialize
        # N=K.ideal(1) or N=OK.ideal(1), which are the same, but for
        # K==QQ it has to be ZZ.ideal(1).
        OK = self.base_ring().ring_of_integers()
        self._conductor = prod([d.prime()**(d.conductor_valuation()) \
                                for d in self.local_data()],\
                               OK.ideal(1))
        return self._conductor

    def global_minimal_model(self, proof = None):
        """
        Returns a model of self that is integral, minimal at all primes.

        Note that this is only implemented for class number 1; in
        general such a model may or may not exist.

        INPUT:
            self -- an elliptic curve over a number field of class number
            proof -- whether to only use provably correct methods
                     (default controled by global proof module).  Note
                     that the proof module is number_field, not
                     elliptic_curves, since the functions that
                     actually need the flag are in number fields.

        OUTPUT:
            A global integral and minimal model.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2-38)
            sage: E = EllipticCurve([0,0,0, 21796814856932765568243810*a - 134364590724198567128296995, 121774567239345229314269094644186997594*a - 750668847495706904791115375024037711300])
            sage: E.global_minimal_model()
            Elliptic Curve defined by y^2 + a*x*y + (a+1)*y = x^3 + (a+1)*x^2 + (368258520200522046806318444*a-2270097978636731786720859345)*x + (8456608930173478039472018047583706316424*a-52130038506793883217874390501829588391299) over Number Field in a with defining polynomial x^2 - 38
        """
        if proof is None:
            import sage.structure.proof.proof
            # We use the "number_field" flag because the actual proof dependence is in Pari's number field functions.
            proof = sage.structure.proof.proof.get_flag(None, "number_field")
        K = self.base_ring()
        if K.class_number() != 1:
            raise ValueError, "global minimal models only exist in general for class number 1"

        E = self.global_integral_model()
        primes = self.base_ring()(self.discriminant()).support()
        for P in primes:
            E = E.local_data(P,proof).minimal_model()
        return E._tidy_model()

    def reduction(self,place):
       """
       Return the reduction of the elliptic curve at a place of good reduction

       INPUT:
            place -- a prime ideal in the base field of the curve

       OUTPUT:
            an elliptic curve over a finite field

       EXAMPLES:
           sage: K.<i> = QuadraticField(-1)
           sage: EK = EllipticCurve([0,0,0,i,i+3])
           sage: v = K.fractional_ideal(2*i+3)
           sage: EK.reduction(v)
           Elliptic Curve defined by y^2  = x^3 + 5*x + 8 over Residue field of Fractional ideal (2*i + 3)
           sage: EK.reduction(K.ideal(1+i))
           Traceback (most recent call last):
           ...
           AttributeError: The curve must have good reduction at the place.
           sage: EK.reduction(K.ideal(2))
           Traceback (most recent call last):
           ...
           AttributeError: The ideal must be prime.
       """
       K = self.base_field()
       OK = K.ring_of_integers()
       try:
           place = K.ideal(place)
       except TypeError:
           raise TypeError, "The parameter must be an ideal of the base field of the elliptic curve"
       if not place.is_prime():
           raise AttributeError, "The ideal must be prime."
       disc = self.discriminant()
       if not K.ideal(disc).valuation(place) == 0:
           raise AttributeError, "The curve must have good reduction at the place."
       Fv = OK.residue_field(place)
       return self.change_ring(Fv)

    def _torsion_bound(self,number_of_places = 20):
        r"""
        Computes an upper bound on the order of the torsion group of
        the elliptic curve by counting points modulo several primes of
        good reduction.  Note that the upper bound returned by this
        function is a multiple of the order of the torsion group.

        INPUT:
            number_of_places (default = 20) -- the number of places
                                that will be used to find the bound

        OUTPUT:
            integer -- the upper bound

        EXAMPLES:
            sage: CDB=CremonaDatabase()
            sage: [E._torsion_bound() for E in CDB.iter([14])]
            [6, 6, 6, 6, 6, 6]
            sage: [E.torsion_order() for E in CDB.iter([14])]
            [6, 6, 2, 6, 2, 6]
        """
        E = self
        bound = 0
        k = 0
        K = E.base_field()
        OK = K.ring_of_integers()
        disc = E.discriminant()
        p = Integer(1)
        # runs through primes, decomposes them into prime ideals
        while k < number_of_places :
            p = p.next_prime()
            f = K.primes_above(p)
            # runs through prime ideals above p
            for qq in f:
                fqq = qq.residue_class_degree()
                charqq = qq.smallest_integer()
                # take only places with small residue field (so that the
                # number of points will be small)
                if fqq == 1 or charqq**fqq < 3*number_of_places:
                    # check if the model is integral at the place
                    if min([K.ideal(a).valuation(qq) for a in E.a_invariants()]) >= 0:
                        eqq = qq.ramification_index()
                        # check if the formal group at the place is torsion-free
                        # if so the torsion injects into the reduction
                        if eqq < charqq - 1 and disc.valuation(qq) == 0:
                            Etilda = E.reduction(qq)
                            Npp = Etilda.cardinality()
                            bound = gcd(bound,Npp)
                            if bound == 1:
                                return bound
                            k += 1
        return bound

    def torsion_subgroup(self):
        """
        Returns the torsion subgroup of this elliptic curve.

        OUTPUT:
            The EllipticCurveTorsionSubgroup instance associated to this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve('11a1')
            sage: K.<t>=NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK=E.base_extend(K)
            sage: tor = EK.torsion_subgroup()
            sage: tor
            Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C5 x C5 associated to the Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 + (-10)*x + (-20) over Number Field in t with defining polynomial x^4 + x^3 + 11*x^2 + 41*x + 101
            sage: tor.gens()
            ((16 : 60 : 1), (t : 1/11*t^3 + 6/11*t^2 + 19/11*t + 48/11 : 1))

            sage: E = EllipticCurve('15a1')
            sage: K.<t>=NumberField(x^2 + 2*x + 10)
            sage: EK=E.base_extend(K)
            sage: EK.torsion_subgroup()
            Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C4 x C4 associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + (-10)*x + (-10) over Number Field in t with defining polynomial x^2 + 2*x + 10

            sage: E = EllipticCurve('19a1')
            sage: K.<t>=NumberField(x^9-3*x^8-4*x^7+16*x^6-3*x^5-21*x^4+5*x^3+7*x^2-7*x+1)
            sage: EK=E.base_extend(K)
            sage: EK.torsion_subgroup()
            Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C9 associated to the Elliptic Curve defined by y^2 + y = x^3 + x^2 + (-9)*x + (-15) over Number Field in t with defining polynomial x^9 - 3*x^8 - 4*x^7 + 16*x^6 - 3*x^5 - 21*x^4 + 5*x^3 + 7*x^2 - 7*x + 1

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve([0,0,0,i,i+3])
            sage: EK.torsion_subgroup ()
            Torsion Subgroup isomorphic to Trivial Abelian Group associated to the Elliptic Curve defined by y^2  = x^3 + i*x + (i+3) over Number Field in i with defining polynomial x^2 + 1
        """
        try:
            return self.__torsion_subgroup
        except AttributeError:
            self.__torsion_subgroup = ell_torsion.EllipticCurveTorsionSubgroup(self)
            return self.__torsion_subgroup

    def torsion_order(self):
        """
        Returns the order of the torsion subgroup of this elliptic curve.

        OUTPUT:
            An integer.

        EXAMPLES:
            sage: E = EllipticCurve('11a1')
            sage: K.<t> = NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            25

            sage: E = EllipticCurve('15a1')
            sage: K.<t> = NumberField(x^2 + 2*x + 10)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            16

            sage: E = EllipticCurve('19a1')
            sage: K.<t> = NumberField(x^9-3*x^8-4*x^7+16*x^6-3*x^5-21*x^4+5*x^3+7*x^2-7*x+1)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_order()
            9

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve([0,0,0,i,i+3])
            sage: EK.torsion_order()
            1
         """
        try:
            return self.__torsion_order
        except AttributeError:
            self.__torsion_order = self.torsion_subgroup().order()
            return self.__torsion_order

    def torsion_points(self):
        """
        Returns a list of the torsion points of this elliptic curve.

        OUTPUT:
            A sorted list of points.

        EXAMPLES:
            sage: E = EllipticCurve('11a1')
            sage: E.torsion_points()
            [(0 : 1 : 0), (5 : -6 : 1), (5 : 5 : 1), (16 : -61 : 1), (16 : 60 : 1)]
            sage: K.<t> = NumberField(x^4 + x^3 + 11*x^2 + 41*x + 101)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_points()
            [(t : 1/11*t^3 + 6/11*t^2 + 19/11*t + 48/11 : 1),
            (1/11*t^3 - 5/11*t^2 + 19/11*t - 40/11 : -6/11*t^3 - 3/11*t^2 - 26/11*t - 321/11 : 1),
            (1/11*t^3 - 5/11*t^2 + 19/11*t - 40/11 : 6/11*t^3 + 3/11*t^2 + 26/11*t + 310/11 : 1),
            (t : -1/11*t^3 - 6/11*t^2 - 19/11*t - 59/11 : 1),
            (16 : 60 : 1),
            (-3/55*t^3 - 7/55*t^2 - 2/55*t - 133/55 : 6/55*t^3 + 3/55*t^2 + 25/11*t + 156/55 : 1),
            (14/121*t^3 - 15/121*t^2 + 90/121*t + 232/121 : 16/121*t^3 - 69/121*t^2 + 293/121*t - 46/121 : 1),
            (-26/121*t^3 + 20/121*t^2 - 219/121*t - 995/121 : -15/121*t^3 - 156/121*t^2 + 232/121*t - 2887/121 : 1),
            (10/121*t^3 + 49/121*t^2 + 168/121*t + 73/121 : -32/121*t^3 - 60/121*t^2 + 261/121*t + 686/121 : 1),
            (5 : 5 : 1),
            (-9/121*t^3 - 21/121*t^2 - 127/121*t - 377/121 : -7/121*t^3 + 24/121*t^2 + 197/121*t + 16/121 : 1),
            (3/55*t^3 + 7/55*t^2 + 2/55*t + 78/55 : 7/55*t^3 - 24/55*t^2 + 9/11*t + 17/55 : 1),
            (-5/121*t^3 + 36/121*t^2 - 84/121*t + 24/121 : -34/121*t^3 + 27/121*t^2 - 305/121*t - 829/121 : 1),
            (5/121*t^3 - 14/121*t^2 - 158/121*t - 453/121 : 49/121*t^3 + 129/121*t^2 + 315/121*t + 86/121 : 1),
            (5 : -6 : 1),
            (5/121*t^3 - 14/121*t^2 - 158/121*t - 453/121 : -49/121*t^3 - 129/121*t^2 - 315/121*t - 207/121 : 1),
            (-5/121*t^3 + 36/121*t^2 - 84/121*t + 24/121 : 34/121*t^3 - 27/121*t^2 + 305/121*t + 708/121 : 1),
            (3/55*t^3 + 7/55*t^2 + 2/55*t + 78/55 : -7/55*t^3 + 24/55*t^2 - 9/11*t - 72/55 : 1),
            (-9/121*t^3 - 21/121*t^2 - 127/121*t - 377/121 : 7/121*t^3 - 24/121*t^2 - 197/121*t - 137/121 : 1),
            (16 : -61 : 1),
            (10/121*t^3 + 49/121*t^2 + 168/121*t + 73/121 : 32/121*t^3 + 60/121*t^2 - 261/121*t - 807/121 : 1),
            (-26/121*t^3 + 20/121*t^2 - 219/121*t - 995/121 : 15/121*t^3 + 156/121*t^2 - 232/121*t + 2766/121 : 1),
            (14/121*t^3 - 15/121*t^2 + 90/121*t + 232/121 : -16/121*t^3 + 69/121*t^2 - 293/121*t - 75/121 : 1),
            (-3/55*t^3 - 7/55*t^2 - 2/55*t - 133/55 : -6/55*t^3 - 3/55*t^2 - 25/11*t - 211/55 : 1),
            (0 : 1 : 0)]

            sage: E = EllipticCurve('15a1')
            sage: K.<t> = NumberField(x^2 + 2*x + 10)
            sage: EK = E.base_extend(K)
            sage: EK.torsion_points()
            [(t : t - 5 : 1),
            (-1 : 0 : 1),
            (t : -2*t + 4 : 1),
            (8 : 18 : 1),
            (1/2 : 5/4*t + 1/2 : 1),
            (-2 : 3 : 1),
            (-7 : 5*t + 8 : 1),
            (3 : -2 : 1),
            (-t - 2 : 2*t + 8 : 1),
            (-13/4 : 9/8 : 1),
            (-t - 2 : -t - 7 : 1),
            (8 : -27 : 1),
            (-7 : -5*t - 2 : 1),
            (-2 : -2 : 1),
            (1/2 : -5/4*t - 2 : 1),
            (0 : 1 : 0)]

            sage: K.<i> = QuadraticField(-1)
            sage: EK = EllipticCurve(K,[0,0,0,0,-1])
            sage: EK.torsion_points ()
            [(-2 : -3*i : 1),
            (0 : -i : 1),
            (1 : 0 : 1),
            (0 : i : 1),
            (-2 : 3*i : 1),
            (0 : 1 : 0)]
         """
        T = self.torsion_subgroup() # make sure it is cached
        return T.points()           # these are also cached in T

    def period_lattice(self, real_embedding):
        r"""
        Returns the period lattice of the elliptic curve for the given real embedding of its base field.

        NOTE:

            Period lattices have not yet been implemented for non-real
            complex embeddings.

        EXAMPLES:
        Define a field with two real embeddings:
            sage: K.<a> = NumberField(x^2-2)
            sage: E=EllipticCurve([0,0,0,a,2])
            sage: embs=K.embeddings(RR); len(embs)
            2

        For each embedding we have a different period lattice:
            sage: E.period_lattice(embs[0])
            Period lattice associated to Elliptic Curve defined by y^2  = x^3 + a*x + 2 over Number Field in a with defining polynomial x^2 - 2 with respect to the real embedding Ring morphism:
            From: Number Field in a with defining polynomial x^2 - 2
            To:   Real Field with 53 bits of precision
            Defn: a |--> -1.414213562373...
            sage: E.period_lattice(embs[1])
            Period lattice associated to Elliptic Curve defined by y^2  = x^3 + a*x + 2 over Number Field in a with defining polynomial x^2 - 2 with respect to the real embedding Ring morphism:
            From: Number Field in a with defining polynomial x^2 - 2
            To:   Real Field with 53 bits of precision
            Defn: a |--> 1.414213562373...

        Although the original embeddings have only the default
        precision, we can obtain the basis with higher precision
        later:
            sage: L=E.period_lattice(embs[0])
            sage: L.basis()
            (4.13107185270501..., -2.06553592635250... + 0.988630424469107...*I)
            sage: L.basis(prec=75)
            (4.131071852705016774309696...,
            -2.065535926352508387154848... + 0.9886304244691077723690104...*I)
            sage: L.basis(prec=100)
            (4.1310718527050167743096955262475367...,
            -2.0655359263525083871548477631237683... + 0.98863042446910777236901069433960633...*I)

        Once increased, the precision does not decrease
            sage: L.basis(prec=10)
            (4.13107185, -2.06553592 + 0.988630424*I) # 32-bit
            (4.13107185270501677, -2.06553592635250838 + 0.988630424469107772*I) # 64-bit
            sage: [CC(w) for w in L.basis()]
            [4.13107185270..., -2.06553592635... + 0.988630424469...*I]
        """
        from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
        return PeriodLattice_ell(self,real_embedding)



r"""
Local data for elliptic curves over number fields (including Q) at primes.

AUTHORS:
    -- John Cremona: First version 2008-09-21 (refactoring code from
       ell_number_field.py and ell_rational_field.py)
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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


from sage.structure.sage_object import SageObject
from sage.misc.misc import verbose, forall

from sage.rings.all import PolynomialRing, QQ, ZZ, Integer, is_Ideal, is_NumberFieldElement, is_NumberFieldFractionalIdeal, is_NumberField
from sage.structure.element import RingElement
from constructor import EllipticCurve
from kodaira_symbol import KodairaSymbol

class EllipticCurveLocalData(SageObject):
    """
    The class for the local reduction data of an elliptic curve.

    Currently supported are elliptic curves defined over Q, and
    elliptic curves defined over a number field, at an arbitrary prime
    or prime ideal.
    """

    def __init__(self, E, P, proof=None, algorithm="pari"):
        """
        Initializes the reduction data for the elliptic curve E at the prime P.

        INPUT:
            E -- an elliptic curve defined over a number field (or QQ)
            P -- a prime ideal of the field
                 (or a prime integer if the field is QQ)
            proof -- whether to only use provably correct methods
                     (default controlled by global proof module).  Note
                     that the proof module is number_field, not
                     elliptic_curves, since the functions that
                     actually need the flag are in number fields.
            algorithm -- str, (default: "pari")
                   (ignored unless E.base_field() is QQ)
                   "pari"   -- use the PARI C-library ellglobalred
                               implementation of Tate's algorithm over QQ.
                   "generic" -- use the general number field implementation.

        EXAMPLES:
            This function is not normally called directly by the user.

            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve('14a1')
            sage: EllipticCurveLocalData(E,2)
            Local data at Principal ideal (2) of Integer Ring:
            Reduction type: bad non-split multiplicative
            Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field
            Minimal discriminant valuation: 6
            Conductor exponent: 1
            Kodaira Symbol: I6
            Tamagawa Number: 2
            sage: EllipticCurveLocalData(E,3)
            Local data at Principal ideal (3) of Integer Ring:
            Reduction type: good
            Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field
            Minimal discriminant valuation: 0
            Conductor exponent: 0
            Kodaira Symbol: I0
            Tamagawa Number: 1
            sage: EllipticCurveLocalData(E,7)
            Local data at Principal ideal (7) of Integer Ring:
            Reduction type: bad split multiplicative
            Local minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field
            Minimal discriminant valuation: 3
            Conductor exponent: 1
            Kodaira Symbol: I3
            Tamagawa Number: 3
        """
        self._curve = E
        K = E.base_field()
        self._prime = check_prime(K,P) # error handling done in that function
        self._reduction_type = None

        if algorithm=="pari" and K is QQ:
            p = self._prime.gen()
            Eint = E.integral_model()
            data = Eint.pari_curve().elllocalred(p)
            self._fp = data[0].python()
            self._KS = KodairaSymbol(data[1].python())
            self._cp = data[3].python()
            # We use a global minimal model since we can:
            self._Emin = Eint.minimal_model()
            self._val_disc = self._Emin.discriminant().valuation(p)
            if self._fp>0:
                self._reduction_type = Eint.ap(p) # = 0,-1 or +1
        else:
            p = self._prime
            self._Emin, ch, self._val_disc, self._fp, self._KS, self._cp, self._split = self._tate(proof)
            if self._fp>0:
                if self._Emin.c4().valuation(p)>0:
                    self._reduction_type = 0
                elif self._split:
                    self._reduction_type = +1
                else:
                    self._reduction_type = -1

    def __repr__(self):
        """
        Returns the string representation of this reduction data.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve('14a1')
            sage: EllipticCurveLocalData(E,2).__repr__()
            'Local data at Principal ideal (2) of Integer Ring:\nReduction type: bad non-split multiplicative\nLocal minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field\nMinimal discriminant valuation: 6\nConductor exponent: 1\nKodaira Symbol: I6\nTamagawa Number: 2'
            sage: EllipticCurveLocalData(E,3).__repr__()
            'Local data at Principal ideal (3) of Integer Ring:\nReduction type: good\nLocal minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field\nMinimal discriminant valuation: 0\nConductor exponent: 0\nKodaira Symbol: I0\nTamagawa Number: 1'
            sage: EllipticCurveLocalData(E,7).__repr__()
            'Local data at Principal ideal (7) of Integer Ring:\nReduction type: bad split multiplicative\nLocal minimal model: Elliptic Curve defined by y^2 + x*y + y = x^3 + 4*x - 6 over Rational Field\nMinimal discriminant valuation: 3\nConductor exponent: 1\nKodaira Symbol: I3\nTamagawa Number: 3'
        """
        red_type = "good"
        if not self._reduction_type is None:
            red_type = ["bad non-split multiplicative","bad additive","bad split multiplicative"][1+self._reduction_type]
        return "Local data at %s:\nReduction type: %s\nLocal minimal model: %s\nMinimal discriminant valuation: %s\nConductor exponent: %s\nKodaira Symbol: %s\nTamagawa Number: %s"%(self._prime,red_type,self._Emin,self._val_disc,self._fp,self._KS,self._cp)

    def minimal_model(self):
        """
        Return the (local) minimal model from this local reduction data.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2  = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.minimal_model()
            Elliptic Curve defined by y^2  = x^3 +1 over Rational Field
            sage: data.minimal_model() == E.local_minimal_model(2)
            True
        """
        return self._Emin

    def prime(self):
        """
        Return the prime ideal associated with this local reduction data.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2  = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.prime()
            Principal ideal (2) of Integer Ring
        """
        return self._prime

    def conductor_valuation(self):
        """
        Return the valuation of the conductor from this local reduction data.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2  = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.conductor_valuation()
            2
        """
        return self._fp

    def kodaira_symbol(self):
        """
        Return the Kodaira symbol from this local reduction data.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2  = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.kodaira_symbol()
            IV
        """
        return self._KS

    def tamagawa_number(self):
        r"""
        Return the Tamagawa number from this local reduction data.

        This is the index $[E(K_v):E^0(K_v)]$.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve([0,0,0,0,64]); E
            Elliptic Curve defined by y^2  = x^3 + 64 over Rational Field
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.tamagawa_number()
            3
        """
        return self._cp

    def tamagawa_exponent(self):
        r"""
        Return the Tamagawa index from this local reduction data.

        This is the exponent of $E(K_v)/E^0(K_v)$; in most cases it is
        the same as the Tamagawa index.

        EXAMPLES:
            sage: from sage.schemes.elliptic_curves.ell_local_data import EllipticCurveLocalData
            sage: E = EllipticCurve('816a1')
            sage: data = EllipticCurveLocalData(E,2)
            sage: data.kodaira_symbol()
            I2*
            sage: data.tamagawa_number()
            4
            sage: data.tamagawa_exponent()
            2

            sage: E = EllipticCurve('200c4')
            sage: data = EllipticCurveLocalData(E,5)
            sage: data.kodaira_symbol()
            I4*
            sage: data.tamagawa_number()
            4
            sage: data.tamagawa_exponent()
            2
        """
        cp = self._cp
        if not cp==4:
            return cp
        ks = self._KS
        if ks._roman==1 and ks._n%2==0 and ks._starred:
            return 2
        return 4

    def bad_reduction_type(self):
        r"""
        Return the type of bad reduction.

        OUTPUT:
            +1 for split multiplicative reduction
            -1 for non-split multiplicative reduction
            0  for additive reduction
            None for good reduction

        EXAMPLES:
            sage: E=EllipticCurve('14a1')
            sage: [(p,E.local_data(p).bad_reduction_type()) for p in prime_range(15)]
            [(2, -1), (3, None), (5, None), (7, 1), (11, None), (13, None)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.local_data(p).bad_reduction_type()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), None), (Fractional ideal (2*a + 1), 0)]
       """
        return self._reduction_type

    def has_good_reduction(self):
        r"""
        Return True if there is good reduction.

        EXAMPLES:
            sage: E=EllipticCurve('14a1')
            sage: [(p,E.local_data(p).has_good_reduction()) for p in prime_range(15)]
            [(2, False), (3, True), (5, True), (7, False), (11, True), (13, True)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.local_data(p).has_good_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), True),
            (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type is None

    def has_bad_reduction(self):
        r"""
        Return True if there is bad reduction.

        EXAMPLES:
            sage: E=EllipticCurve('14a1')
            sage: [(p,E.local_data(p).has_bad_reduction()) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.local_data(p).has_bad_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
            (Fractional ideal (2*a + 1), True)]
        """
        return not self._reduction_type is None

    def has_multiplicative_reduction(self):
        r"""
        Return True if there is multiplicative reduction.

        See also has_split_multiplicative_reduction() and
                 has_nonsplit_multiplicative_reduction().

        EXAMPLES:
            sage: E=EllipticCurve('14a1')
            sage: [(p,E.local_data(p).has_multiplicative_reduction()) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.local_data(p).has_multiplicative_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False), (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type in (-1,+1)

    def has_split_multiplicative_reduction(self):
        r"""
        Return True if there is split multiplicative  reduction.

        EXAMPLES:
            sage: E=EllipticCurve('14a1')
            sage: [(p,E.local_data(p).has_split_multiplicative_reduction()) for p in prime_range(15)]
            [(2, False), (3, False), (5, False), (7, True), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.local_data(p).has_split_multiplicative_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
            (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type == +1

    def has_nonsplit_multiplicative_reduction(self):
        r"""
        Return True if there is non-split multiplicative  reduction.

        EXAMPLES:
            sage: E=EllipticCurve('14a1')
            sage: [(p,E.local_data(p).has_nonsplit_multiplicative_reduction()) for p in prime_range(15)]
            [(2, True), (3, False), (5, False), (7, False), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.local_data(p).has_nonsplit_multiplicative_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False), (Fractional ideal (2*a + 1), False)]
        """
        return self._reduction_type == -1

    def has_additive_reduction(self):
        r"""
        Return True if there is additive reduction.

        EXAMPLES:
            sage: E=EllipticCurve('27a1')
            sage: [(p,E.local_data(p).has_additive_reduction()) for p in prime_range(15)]
            [(2, False), (3, True), (5, False), (7, False), (11, False), (13, False)]

            sage: K.<a>=NumberField(x^3-2)
            sage: P17a, P17b = [P for P,e in K.factor(17)]
            sage: E = EllipticCurve([0,0,0,0,2*a+1])
            sage: [(p,E.local_data(p).has_additive_reduction()) for p in [P17a,P17b]]
            [(Fractional ideal (4*a^2 - 2*a + 1), False),
            (Fractional ideal (2*a + 1), True)]
        """
        return self._reduction_type == 0

    def _tate(self, proof = None):
        """
        Tate's algorithm for an elliptic curve over a number field:
        computes local reduction data at a prime ideal and a local
        minimal model.

        The model is not required to be integral on input.  If P is
        principal, uses a generator as uniformizer, so it will not
        affect integrality or minimality at other primes.  If P is not
        principal, the minimal model returned will preserve
        integrality at other primes, but not minimality.

        INPUT:
            Called only by the EllipticCurveLocalData.__init__()

        OUTPUT:
            Emin -- a model (integral and) minimal at P
            p    -- the residue characteristic
            val_disc  -- the valuation of the local minimal discriminant
            fp   -- valuation of the conductor
            KS   -- Kodaira symbol
            cp   -- Tamagawa number

        """
        E = self._curve
        P = self._prime
        K = E.base_ring()
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
        prime = pi if K is QQ else P

        pval = lambda x: x.valuation(prime)
        pdiv = lambda x: x.is_zero() or pval(x) > 0
        pinv = lambda x: F.lift(~F(x))
        proot = lambda x,e: F.lift(F(x).nth_root(e, extend = False, all = True)[0])
        preduce = lambda x: F.lift(F(x))

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
            halfmodp = pinv(Integer(2))

        A = E.a_invariants()
        A = [0, A[0], A[1], A[2], A[3], 0, A[4]]
        indices = [1,2,3,4,6]
        if min([pval(a) for a in A if a != 0]) < 0:
            verbose("Non-integral model at P: valuations are %s; making integral"%([pval(a) for a in A if a != 0]), t, 1)
            e = 0
            for i in range(7):
                if A[i] != 0:
                    e = max(e, (-pval(A[i])/i).ceil())
            pie = pi**e
            for i in range(7):
                if A[i] != 0:
                    A[i] *= pie**i
            verbose("P-integral model is %s, with valuations %s"%([A[i] for i in indices], [pval(A[i]) for i in indices]), t, 1)

        split = None # only relevant for multiplicative reduction

        (a1, a2, a3, a4, a6) = (A[1], A[2], A[3], A[4], A[6])
        while True:
            C = EllipticCurve([a1, a2, a3, a4, a6]);
            (b2, b4, b6, b8) = C.b_invariants()
            (c4, c6) = C.c_invariants()
            delta = C.discriminant()
            val_disc = pval(delta)

            if val_disc == 0:
                ## Good reduction already
                cp = 1
                fp = 0
                KS = KodairaSymbol("I0")
                break #return

            # Otherwise, we change coordinates so that p | a3, a4, a6
            if p == 2:
                if pdiv(b2):
                    r = proot(a4, 2)
                    t = proot(((r + a2)*r + a4)*r + a6, 2)
                else:
                    temp = pinv(a1)
                    r = temp * a3
                    t = temp * (a4 + r*r)
            elif p == 3:
                if pdiv(b2):
                    r = proot(-b6, 3)
                else:
                    r = -pinv(b2) * b4
                t = a1 * r + a3
            else:
                if pdiv(c4):
                    r = -pinv(12) * b2
                else:
                    r = -pinv(12*c4) * (c6 + b2 * c4)
                t = -halfmodp * (a1 * r + a3)
            r = preduce(r)
            t = preduce(t)
            # print "Before first tranform C = %s"%C
            # print "[a1,a2,a3,a4,a6] = %s"%([a1, a2, a3, a4, a6])
            C = C.rst_transform(r, 0, t)
            (a1, a2, a3, a4, a6) = C.a_invariants()
            (b2, b4, b6, b8) = C.b_invariants()
            if min([pval(a) for a in (a1, a2, a3, a4, a6) if a != 0]) < 0:
                raise RuntimeError, "Non-integral model after first transform!"
            verbose("After first transform %s\n, [a1,a2,a3,a4,a6] = %s\n, valuations = %s"%([r, 0, t], [a1, a2, a3, a4, a6], [pval(a1), pval(a2), pval(a3), pval(a4), pval(a6)]), t, 2)
            if pval(a3) == 0:
                raise RuntimeError, "p does not divide a3 after first transform!"
            if pval(a4) == 0:
                raise RuntimeError, "p does not divide a4 after first transform!"
            if pval(a6) == 0:
                raise RuntimeError, "p does not divide a6 after first transform!"

            # Now we test for Types In, II, III, IV
            # NB the c invariants never change.

            if not pdiv(c4):
                # Multiplicative reduction: Type In (n = val_disc)
                split = False
                if _pquadroots(1, a1, -a2):
                    cp = val_disc
                    split = True
                elif Integer(2).divides(val_disc):
                    cp = 2
                else:
                    cp = 1
                KS = KodairaSymbol("I%s"%val_disc)
                fp = 1
                break #return

            # Additive reduction

            if pval(a6) < 2:
                ## Type II
                KS = KodairaSymbol("II")
                fp = val_disc
                cp = 1
                break #return
            if pval(b8) < 3:
                ## Type III
                KS = KodairaSymbol("III")
                fp = val_disc - 1
                cp = 2
                break #return
            if pval(b6) < 3:
                ## Type IV
                if _pquadroots(1, a3 / pi, -a6/(pi*pi)):
                    cp = 3
                else:
                    cp = 1
                KS = KodairaSymbol("IV")
                fp = val_disc - 2
                break #return

            # If our curve is none of these types, we change types so that p | a1, a2 and p^2 | a3, a4 and p^3 | a6
            if p == 2:
                s = proot(a2, 2)
                t = pi*proot(a6/(pi*pi), 2)
            elif p == 3:
                s = a1
                t = a3
            else:
                s = -a1*halfmodp
                t = -a3*halfmodp
            C = C.rst_transform(0, s, t)
            (a1, a2, a3, a4, a6) = C.a_invariants()
            (b2, b4, b6, b8) = C.b_invariants()
            verbose("After second transform %s\n[a1, a2, a3, a4, a6] = %s\nValuations: %s"%([0, s, t], [a1,a2,a3,a4,a6],[pval(a1),pval(a2),pval(a3),pval(a4),pval(a6)]), t, 2)
            if pval(a1) == 0:
                raise RuntimeError, "p does not divide a1 after second transform!"
            if pval(a2) == 0:
                raise RuntimeError, "p does not divide a2 after second transform!"
            if pval(a3) < 2:
                raise RuntimeError, "p^2 does not divide a3 after second transform!"
            if pval(a4) < 2:
                raise RuntimeError, "p^2 does not divide a4 after second transform!"
            if pval(a6) < 3:
                raise RuntimeError, "p^3 does not divide a6 after second transform!"
            if min(pval(a1), pval(a2), pval(a3), pval(a4), pval(a6)) < 0:
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
            if pdiv(w):
                if pdiv(x):
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
                fp = val_disc - 4
                break #return
            elif sw == 2:
                ## One double root - Type I*m for some m
                verbose("One double root", t, 1)
                ## Change coords so that the double root is T = 0 mod p
                if p == 2:
                    r = proot(c, 2)
                elif p == 3:
                    r = c * pinv(b)
                else:
                    r = (bc - 9*d)*pinv(2*x)
                r = pi * preduce(r)
                C = C.rst_transform(r, 0, 0)
                (a1, a2, a3, a4, a6) = C.a_invariants()
                (b2, b4, b6, b8) = C.b_invariants()
                ix = 3; iy = 3; mx = pi*pi; my = pi*pi
                while True:
                    a2t = a2 / pi
                    a3t = a3 / my
                    a4t = a4 / (pi*mx)
                    a6t = a6 / (mx*my)
                    if pdiv(a3t*a3t + 4*a6t):
                        if p == 2:
                            t = my*proot(a6t, 2)
                        else:
                            t = my*preduce(-a3t*halfmodp)
                        C = C.rst_transform(0, 0, t)
                        (a1, a2, a3, a4, a6) = C.a_invariants()
                        (b2, b4, b6, b8) = C.b_invariants()
                        my = my*pi
                        iy += 1
                        a2t = a2/pi
                        a3t = a3/my
                        a4t = a4/(pi*mx)
                        a6t = a6/(mx*my)
                        if pdiv(a4t*a4t - 4*a6t*a2t):
                            if p == 2:
                                r = mx*proot(a6t*pinv(a2t), 2)
                            else:
                                r = mx*preduce(-a4t*pinv(2*a2t))
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
                fp = val_disc - ix - iy + 1
                break #return
            else: # sw == 3
                ## The cubic has a triple root
                verbose("Triple root", t, 1)
                ## First we change coordinates so that T = 0 mod p
                if p == 2:
                    r = b
                elif p == 3:
                    r = proot(-d, 3)
                else:
                    r = -b * pinv(3)
                r = pi*preduce(r)
                C = C.rst_transform(r, 0, 0)
                (a1, a2, a3, a4, a6) = C.a_invariants()
                (b2, b4, b6, b8) = C.b_invariants()
                verbose("After third transform %s\n[a1,a2,a3,a4,a6] = %s\nValuations: %s"%([r,0,0],[a1,a2,a3,a4,a6],[pval(ai) for ai in [a1,a2,a3,a4,a6]]), t, 2)
                if min(pval(ai) for ai in [a1,a2,a3,a4,a6]) < 0:
                    raise RuntimeError, "Non-integral model after third transform!"
                if pval(a2) < 2 or pval(a4) < 3 or pval(a6) < 4:
                    raise RuntimeError, "Cubic after transform does not have a triple root at 0"
                a3t = a3/(pi*pi)
                a6t = a6/(pi**4)
                # We test for Type IV*
                if not pdiv(a3t*a3t + 4*a6t):
                    cp = 3 if _pquadroots(1, a3t, -a6t) else 1
                    KS = KodairaSymbol("IV*")
                    fp = val_disc - 6
                    break #return
                # Now change coordinates so that p^3|a3, p^5|a6
                t =        -pi*pi*proot(a6t, 2) if p==2 \
                      else  pi*pi*preduce(-a3t*halfmodp)
                C = C.rst_transform(0, 0, t)
                (a1, a2, a3, a4, a6) = C.a_invariants()
                (b2, b4, b6, b8) = C.b_invariants()
                # We test for types III* and II*
                if pval(a4) < 4:
                    ## Type III*
                    KS = KodairaSymbol("III*")
                    fp = val_disc - 7
                    cp = 2
                    break #return
                if pval(a6) < 6:
                    ## Type II*
                    KS = KodairaSymbol("II*")
                    fp = val_disc - 8
                    cp = 1
                    break #return
                a1 /= pi
                a2 /= pi**2
                a3 /= pi**3
                a4 /= pi**4
                a6 /= pi**6
                verbose("Non-minimal equation, dividing out...\nNew model is %s"%([a1, a2, a3, a4, a6]), t, 1)
        C = C._tidy_model()
        return (C, p, val_disc, fp, KS, cp, split)


def check_prime(K,P):
    """
    Function to check that P determines a prime of K, and return that ideal.

    INPUT:
        K -- a number field (including QQ)
        P -- an element of K or a (fractional) ideal of K

    OUTPUT:
        A prime ideal equal to or generated by P, or raise an error if invalid.

    EXAMPLES:
        sage: from sage.schemes.elliptic_curves.ell_local_data import check_prime
        sage: check_prime(QQ,3)
        Principal ideal (3) of Integer Ring
        sage: check_prime(QQ,ZZ.ideal(31))
        Principal ideal (31) of Integer Ring
        sage: K.<a>=NumberField(x^2-5)
        sage: check_prime(K,a)
        Fractional ideal (a)
        sage: check_prime(K,a+1)
        Fractional ideal (a + 1)
        sage: [check_prime(K,P) for P in K.primes_above(31)]
        [Fractional ideal (-5/2*a - 1/2), Fractional ideal (-5/2*a + 1/2)]
    """
    if K is QQ:
        if isinstance(P, (int,long,Integer)):
            P = ZZ.ideal(Integer(P))
        if is_Ideal(P) and P.base_ring() is ZZ and P.is_prime():
            return P
        raise TypeError, "%s is not a prime ideal of %s"%(P,ZZ)

    if not is_NumberField(K):
        raise TypeError, "%s is not a number field"%K

    if is_NumberFieldFractionalIdeal(P):
        if P.is_prime():
            return P
        else:
            raise TypeError, "%s is not a prime ideal of %s"%(P,K)

    if is_NumberFieldElement(P):
        if P in K:
            P = K.ideal(P)
        else:
            raise TypeError, "%s is not an element of %s"%(P,K)
        if P.is_prime():
            return P
        else:
            raise TypeError, "%s is not a prime ideal of %s"%(P,K)

    raise TypeError, "%s is not a valid prime of %s"%(P,K)

"""
Elliptic curves over the rational numbers

AUTHORS:
   -- William Stein (2005): first version
   -- William Stein (2006-02-26): fixed Lseries_extended which didn't work
            because of changes elsewhere in SAGE.
   -- David Harvey (2006-09): Added padic_E2, padic_sigma, padic_height,
            padic_regulator methods.
   -- David Harvey (2007-02): reworked padic-height related code
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

import ell_point
import formal_group
import rational_torsion
from ell_field import EllipticCurve_field

import sage.groups.all
import sage.rings.arith as arith
import sage.rings.all as rings
import sage.rings.number_field.number_field as number_field
import sage.misc.misc as misc
import sage.functions.constants as constants
import sage.modular.modform.constructor
import sage.modular.modform.element
from sage.misc.functional import log
from sage.rings.padics.zp import Zp
from sage.rings.padics.qp import Qp

# Use some interval arithmetic to guarantee correctness.  We assume
# that alpha is computed to the precision of a float.
IR = rings.RIF
#from sage.rings.interval import IntervalRing; IR = IntervalRing()

import sage.matrix.all as matrix
import sage.databases.cremona
from   sage.libs.pari.all import pari
import sage.functions.transcendental as transcendental
import math
import sage.libs.mwrank.all as mwrank
import constructor
from sage.interfaces.all import gp

import ell_modular_symbols
import padic_lseries

import mod5family

from sage.rings.all import (
    PowerSeriesRing, LaurentSeriesRing, O,
    infinity as oo,
    Integer,
    IntegerRing, RealField,
    ComplexField, RationalField)

import gp_cremona
import padic_height
import monsky_washnitzer
import sea

from gp_simon import simon_two_descent

factor = arith.factor
sqrt = math.sqrt
exp = math.exp
mul = misc.mul
next_prime = arith.next_prime

Q = RationalField()
Z = IntegerRing()
C = ComplexField()
R = RealField()


_MAX_HEIGHT=21

class EllipticCurve_rational_field(EllipticCurve_field):
    """
    Elliptic curve over the Rational Field.
    """
    def __init__(self, ainvs, extra=None):
        if extra != None:   # possibility of two arguments (the first would be the field)
            ainvs = extra
        if isinstance(ainvs, str):
            label = ainvs
            X = sage.databases.cremona.CremonaDatabase()[label]
            EllipticCurve_field.__init__(self, [Q(a) for a in X.a_invariants()])
            for attr in ['rank', 'torsion_order', 'cremona_label', 'conductor',
                         'modular_degree', 'gens', 'regulator']:
                s = "_EllipticCurve_rational_field__"+attr
                if hasattr(X,s):
                    setattr(self, s, getattr(X, s))
            return
        EllipticCurve_field.__init__(self, [Q(x) for x in ainvs])
        self.__np = {}
        if self.base_ring() != Q:
            raise TypeError, "Base field (=%s) must be the Rational Field."%self.base_ring()

    def _set_rank(self, r):
        self.__rank = r
    def _set_torsion_order(self, t):
        self.__torsion_order = t
    def _set_cremona_label(self, L):
        self.__cremona_label = L
    def _set_conductor(self, N):
        self.__conductor_pari = Z(N)
    def _set_modular_degree(self, deg):
        self.__modular_degree = deg

    def _set_gens(self, gens):
        self.__gens = [self.point(x, check=True) for x in gens]
        self.__gens_certain = True
        self.__gens.sort()

    def is_integral(self):
        try:
            return self.__is_integral
        except AttributeError:
            one = Z(1)
            self.__is_integral = bool(misc.mul([x.denominator() == 1 for x in self.ainvs()]))
            return self.__is_integral


    def mwrank(self, options=''):
        """
        Run Cremona's mwrank program on this elliptic curve and
        return the result as a string.

        INPUT:
            options -- string; passed when starting mwrank.  The format is
        q p<precision> v<verbosity> b<hlim_q> x<naux>  c<hlim_c> l t o s d>]

        OUTPUT:
            string -- output of mwrank on this curve
        """
        if options == "":
            from sage.interfaces.all import mwrank
        else:
            from sage.interfaces.all import Mwrank
            mwrank = Mwrank(options=options)
        return mwrank(self.a_invariants())

    def conductor(self, algorithm="pari"):
        """
        Returns the conductor of the elliptic curve.

        INPUT:
            algorithm -- str, (default: "pari")
                   "pari"   -- use the PARI C-library ellglobalred
                               implementation of Tate's algorithm
                   "mwrank" -- use Cremona's mwrank implementation of
                               Tate's algorithm; can be faster if the
                               curve has integer coefficients (TODO:
                               limited to small conductor until mwrank
                               gets integer factorization)
                   "gp" -- use the GP interpreter.
                   "all" -- use both implementations, verify that the
                            results are the same (or raise an error),
                            and output the common value.

        EXAMPLE:
            sage: E = EllipticCurve([1, -1, 1, -29372, -1932937])
            sage: E.conductor(algorithm="pari")
            3006
            sage: E.conductor(algorithm="mwrank")
            3006
            sage: E.conductor(algorithm="gp")
            3006
            sage: E.conductor(algorithm="all")
            3006

        NOTE: The conductor computed using each algorithm is cached separately.
        Thus calling E.conductor("pari"), then E.conductor("mwrank") and
        getting the same result checks that both systems compute the same answer.
        """

        if algorithm == "pari":
            try:
                return self.__conductor_pari
            except AttributeError:
                self.__conductor_pari = Z(self.pari_mincurve().ellglobalred()[0])
            return self.__conductor_pari

        elif algorithm == "gp":
            try:
                return self.__conductor_gp
            except AttributeError:
                self.__conductor_gp = Z(gp.eval('ellglobalred(ellinit(%s,0))[1]'%self.a_invariants()))
                return self.__conductor_gp

        elif algorithm == "mwrank":
            try:
                return self.__conductor_mwrank
            except AttributeError:
                if self.is_integral():
                    self.__conductor_mwrank = Z(self.mwrank_curve().conductor())
                else:
                    self.__conductor_mwrank = Z(self.minimal_model().mwrank_curve().conductor())
            return self.__conductor_mwrank

        elif algorithm == "all":
            N1 = self.conductor("pari")
            N2 = self.conductor("mwrank")
            N3 = self.conductor("gp")
            if N1 != N2 or N2 != N3:
                raise ArithmeticError, "Pari, mwrank and gp compute different conductors (%s,%s,%s) for %s"%(
                    N1, N2, N3, self)
            return N1
        else:
            raise RuntimeError, "algorithm '%s' is not known."%algorithm

    ####################################################################
    #  Access to PARI curves related to this curve.
    ####################################################################

    def __pari_double_prec(self):
        """
        Double the precision of computations with this curve.
        """
        self.pari_curve()
        prec = 2 * self.__pari_curve[1]
        self.__pari_curve = (pari(self.ainvs()).ellinit(precision=prec), prec)
        try:
            del self.__pari_mincurve
        except AttributeError:
            pass

    def pari_curve(self):
        """
        Return the PARI curve corresponding to this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve([0, 0,1,-1,0])
            sage: e = E.pari_curve()
            sage: type(e)
            <type 'sage.libs.pari.gen.gen'>
            sage: e.type()
            't_VEC'
            sage: e.ellan(10)
            [1, -2, -3, 2, -2, 6, -1, 0, 6, 4]

            sage: E = EllipticCurve(RationalField(), ['1/3', '2/3'])
            sage: e = E.pari_curve()
            sage: e.type()
            't_VEC'
            sage: e[:5]
            [0, 0, 0, 1/3, 2/3]
        """
        try:
            return self.__pari_curve[0]
        except AttributeError:
            self.__pari_curve = (pari(self.a_invariants()).ellinit(precision=10), 10)
        return self.__pari_curve[0]

    def pari_mincurve(self):
        """
        Return the PARI curve corresponding to a minimal model
        for this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve(RationalField(), ['1/3', '2/3'])
            sage: e = E.pari_mincurve()
            sage: e[:5]
            [0, 0, 0, 27, 486]
            sage: E.conductor()
            47232
            sage: e.ellglobalred()
            [47232, [1, 0, 0, 0], 2]
        """
        try:
            return self.__pari_mincurve
        except AttributeError:
            e = self.pari_curve()
            mc, change = e.ellminimalmodel()
            self.__pari_mincurve = mc
            # self.__min_transform = change
        return self.__pari_mincurve

    def database_curve(self):
        """
        Return the curve in the elliptic curve database isomorphic to
        this curve, if possible.  Otherwise raise a RuntimeError
        exception.

        EXAMPLES:
            sage: E = EllipticCurve([0,1,2,3,4])
            sage: E.database_curve()
            Elliptic Curve defined by y^2  = x^3 + x^2 + 3*x + 5 over Rational Field

        NOTES: The model of the curve in the database can be different
               than the Weierstrass model for this curve, e.g.,
               database models are always minimal.
        """
        try:
            return self.__database_curve
        except AttributeError:
            misc.verbose("Looking up %s in the database."%self)
            D = sage.databases.cremona.CremonaDatabase()
            ainvs = self.minimal_model().ainvs()
            try:
                self.__database_curve = D.elliptic_curve_from_ainvs(self.conductor(), ainvs)
            except RuntimeError:
                raise RuntimeError, "Elliptic curve %s not in the database."%self
            return self.__database_curve


    def Np(self, p):
        """
        The number of points on E modulo p, where p is a prime, not
        necessarily of good reduction.  (When p is a bad prime, also
        counts the singular point.)

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.Np(2)
            5
            sage: E.Np(3)
            5
            sage: E.conductor()
            11
            sage: E.Np(11)
            11
        """
        if self.conductor() % p == 0:
            return p + 1 - self.ap(p)
        #raise ArithmeticError, "p (=%s) must be a prime of good reduction"%p
        if p < 1125899906842624:   # TODO: choose more wisely?
            return p+1 - self.ap(p)
        else:
            return self.sea(p)

    def sea(self, p, early_abort=False):
        r"""
        Return the number of points on $E$ over $\F_p$ computed using
        the SEA algorithm, as implemented in PARI by Christophe Doche
        and Sylvain Duquesne.

        INPUT:
            p -- a prime number
            early_abort -- bool (default: Falst); if True an early abort technique
                       is used and the computation is interrupted as soon
                       as a small divisor of the order is detected.

        \note{As of 2006-02-02 this function does not work on
        Microsoft Windows under Cygwin (though it works under
        vmware of course).}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.sea(next_prime(10^30))
            1000000000000001426441464441649
        """
        p = rings.Integer(p)
        return sea.ellsea(self.minimal_model().a_invariants(), p, early_abort=early_abort)

    def __pari_double_prec(self):
        EllipticCurve_field._EllipticCurve__pari_double_prec(self)
        try:
            del self.__pari_mincurve
        except AttributeError:
            pass

    ####################################################################
    #  Access to mwrank
    ####################################################################
    def mwrank_curve(self, verbose=False):
        try:
            return self.__mwrank_curve
        except AttributeError:
            pass
        self.__mwrank_curve = mwrank.mwrank_EllipticCurve(
            self.ainvs(), verbose=verbose)
        return self.__mwrank_curve

    def two_descent(self, verbose=True,
                    selmer_only = False,
                    first_limit = 20,
                    second_limit = 8,
                    n_aux = -1,
                    second_descent = 1):
        """
        Compute 2-descent data for this curve.

        INPUT:
            verbose     -- (default: True) print what mwrank is doing
            selmer_only -- (default: False) selmer_only switch
            first_limit -- (default: 20) firstlim is bound on |x|+|z|
            second_limit-- (default: 8)  secondlim is bound on log max {|x|,|z| },
                                         i.e. logarithmic
            n_aux       -- (default: -1) n_aux only relevant for general
                           2-descent when 2-torsion trivial; n_aux=-1 causes default
                           to be used (depends on method)
            second_descent -- (default: True) second_descent only relevant for
                           descent via 2-isogeny
        OUTPUT:
            Nothing -- nothing is returned (though much is printed)
        """
        self.mwrank_curve().two_descent(verbose, selmer_only,
                                        first_limit, second_limit,
                                        n_aux, second_descent)


    ####################################################################
    #  Etc.
    ####################################################################

    def aplist(self, pmax):
        """
        Return list of pairs (p, a_p(E)) for p up to pmax.
        """
        v = arith.prime_range(pmax)
        return [(p,self.ap(p)) for p in v]

    def anlist(self, n, pari_ints=False):
        """
        The Fourier coefficients up to and including $a_n$ of the
        modular form attached to this elliptic curve.  The ith element
        of the return list is a[i].

        INPUT:
            n -- integer
            pari_ints -- bool (default: False); if True return a list of
                      PARI ints instead of SAGE integers; this can
                      be much faster for large n.

        OUTPUT:
            -- list of integers

        If pari_ints is False, the result is cached.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.anlist(3)
            [0, 1, -2, -1]

            sage: E = EllipticCurve([0,1])
            sage: E.anlist(20)
            [0, 1, 0, 0, 0, 0, 0, -4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 8, 0]
        """
        n = int(n)
        if not pari_ints:
            try:
                if len(self.__anlist) > n:
                    return self.__anlist[:n+1]
            except AttributeError:
                pass
        E = self.pari_mincurve()
        if n >= 2147483648:
            raise RuntimeError, "anlist: n (=%s) must be < 2147483648."%n

        if not pari_ints:
            ZZ = rings.Integer
            v = [ZZ(0)] + [ZZ(x) for x in E.ellan(n)]
        else:
            v = E.ellan(n)
        if not pari_ints:
            self.__anlist = v
        return v


        # There is some overheard associated with coercing the PARI
        # list back to Python, but it's not bad.  It's better to do it
        # this way instead of trying to eval the whole list, since the
        # int conversion is done very sensibly.  NOTE: This would fail
        # if a_n won't fit in a C int, i.e., is bigger than
        # 2147483648; however, we wouldn't realistically compute
        # anlist for n that large anyways.
        #
        # Some relevant timings:
        #
        # E <--> [0, 1, 1, -2, 0]   389A
        #  E = EllipticCurve([0, 1, 1, -2, 0]);   // SAGE or MAGMA
        #  e = E.pari_mincurve()
        #  f = ellinit([0,1,1,-2,0]);
        #
        #  Computation                                              Time (1.6Ghz Pentium-4m laptop)
        #  time v:=TracesOfFrobenius(E,10000);  // MAGMA            0.120
        #  gettime;v=ellan(f,10000);gettime/1000                    0.046
        #  time v=e.ellan (10000)                                   0.04
        #  time v=E.anlist(10000)                                   0.07

        #  time v:=TracesOfFrobenius(E,100000);  // MAGMA           1.620
        #  gettime;v=ellan(f,100000);gettime/1000                   0.676
        #  time v=e.ellan (100000)                                  0.7
        #  time v=E.anlist(100000)                                  0.83

        #  time v:=TracesOfFrobenius(E,1000000);  // MAGMA          20.850
        #  gettime;v=ellan(f,1000000);gettime/1000                  9.238
        #  time v=e.ellan (1000000)                                 9.61
        #  time v=E.anlist(1000000)                                 10.95  (13.171 in cygwin vmware)

        #  time v:=TracesOfFrobenius(E,10000000);  //MAGMA          257.850
        #  gettime;v=ellan(f,10000000);gettime/1000      FAILS no matter how many allocatemem()'s!!
        #  time v=e.ellan (10000000)                                139.37
        #  time v=E.anlist(10000000)                                136.32
        #
        #  The last SAGE comp retries with stack size 40MB,
        #  80MB, 160MB, and succeeds last time.  It's very interesting that this
        #  last computation is *not* possible in GP, but works in py_pari!
        #

    def q_expansion(self, prec):
        """
        Return the q-expansion to precision prec of the newform attached to this
        elliptic curve.

        INPUT:
            prec -- an integer
        """
        return PowerSeriesRing(Q, 'q')(self.anlist(prec), prec, check=True)

    def modular_form(self):
        r"""
        Return the cuspidal modular form associated to this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: f = E.modular_form()
            sage: f
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)

        NOTE: If you just want the $q$-expansion, use
        \code{self.q_expansion(prec)}.
        """
        try:
            return self.__modular_form
        except AttributeError:
            M = sage.modular.modform.constructor.ModularForms(self.conductor(),weight=2)
            f = sage.modular.modform.element.ModularFormElement_elliptic_curve(M, self, None)
            self.__modular_form = f
            return f

    def modular_symbol_space(self, sign=1, base_ring=Q):
        r"""
        Return the space of cuspidal modular symbols associated to
        this elliptic curve, with given sign and base ring.

        INPUT:
            sign -- 0, -1, or 1
            base_ring -- a ring

        EXAMPLES:

        NOTE: If you just want the $q$-expansion, use
        \code{self.q_expansion(prec)}.
        """
        typ = (sign, base_ring)
        try:
            return self.__modular_symbol_space[typ]
        except AttributeError:
            self.__modular_symbol_space = {}
        except KeyError:
            pass
        M = ell_modular_symbols.modular_symbol_space(self, sign, base_ring)
        self.__modular_symbol_space[typ] = M
        return M

    def modular_symbol(self, sign=1, base_ring=Q):
        r"""
        Return the modular symbol associated to this elliptic curve,
        with given sign and base ring.  This is the map that sends r/s
        to a fixed multiple of 2*pi*I*f(z)dz from oo to r/s,
        normalized so that all values of this map take values in QQ.

        NOTE: Currently there is no guarantee about how this map is
        normalized.  This will be added.

        INPUT:
            sign -- -1, or 1
            base_ring -- a ring

        NOTE: If you just want the $q$-expansion, use
        \code{self.q_expansion(prec)}.
        """
        typ = (sign, base_ring)
        try:
            return self.__modular_symbol[typ]
        except AttributeError:
            self.__modular_symbol = {}
        except KeyError:
            pass
        M = ell_modular_symbols.ModularSymbol(self, sign, base_ring)
        self.__modular_symbol[typ] = M
        return M

    def padic_lseries(self, p, prec=20):
        """
        Return the p-adic Lseries of self at p with given p-adic precision.

        INPUT:
            p -- prime
            prec -- precision of p-adic computations

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: L = E.padic_lseries(5); L
            5-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: type(L)
            <class 'sage.schemes.elliptic_curves.padic_lseries.pAdicLseriesOrdinary'>
        """
        key = (p,prec)
        try:
            return self._padic_lseries[key]
        except AttributeError:
            self._padic_lseries = {}
        except KeyError:
            pass
        if self.ap(p) % p != 0:
            Lp = padic_lseries.pAdicLseriesOrdinary(self, p, prec)
        else:
            Lp = padic_lseries.pAdicLseriesSupersingular(self, p, prec)
        self._padic_lseries[key] = Lp
        return Lp

    def newform(self):
        """
        Same as \code{self.modular_form()}.
        """
        return self.modular_form()

    def q_eigenform(self, prec):
        """
        Synonym for self.q_expansion(prec).
        """
        return self.q_expansion(prec)

    def analytic_rank(self, algorithm="cremona"):
        r"""
        Return an integer that is \emph{probably} the analytic rank of
        this elliptic curve.

        INPUT:
            algorithm:
                -- 'cremona' (default) --  Use the Buhler-Gross algorithm
                    as implemented in GP by Tom Womack and John Cremona,
                    who note that their implementation is practical for
                    any rank and conductor $\leq 10^{10}$ in 10 minutes.

                -- 'ec' -- use Watkins's program ec (this has bugs if more
                    than a million terms of the L-series are required, i.e.,
                    only use this for conductor up to about $10^11$).

                -- 'sympow' --use Watkins's program sympow

                -- 'rubinstein' -- use Rubinstein's L-function C++ program lcalc.

                -- 'magma' -- use MAGMA

                -- 'all' -- compute with all other free algorithms, check that the
                            answers agree, and return the common answer.

        \note{If the curve is loaded from the large Cremona database,
        then the modular degree is taken from the database.}

        Of the three above, probably Rubinstein's is the most
        efficient (in some limited testing I've done).

        \note{It is an open problem to \emph{prove} that \emph{any}
        particular elliptic curve has analytic rank $\geq 4$.}

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage.: E.analytic_rank(algorithm='ec')
            2
            sage: E.analytic_rank(algorithm='cremona')
            2
            sage: E.analytic_rank(algorithm='rubinstein')
            2
            sage: E.analytic_rank(algorithm='sympow')
            2
            sage: E.analytic_rank(algorithm='magma')    # optional
            2
            sage.: E.analytic_rank(algorithm='all')
            2
        """
        if algorithm == 'ec' and misc.is_64_bit:
            algorithm = 'sympow'

        if algorithm == 'cremona':
            return rings.Integer(gp_cremona.ellanalyticrank(self.minimal_model().a_invariants()))
        elif algorithm == 'ec':
            return rings.Integer(self.watkins_data()['analytic rank'])
        elif algorithm == 'rubinstein':
            from sage.lfunctions.lcalc import lcalc
            return lcalc.analytic_rank(L=self)
        elif algorithm == 'sympow':
            from sage.lfunctions.sympow import sympow
            return sympow.analytic_rank(self)[0]
        elif algorithm == 'magma':
            return rings.Integer(self._magma_().AnalyticRank())
        elif algorithm == 'all':
            S = list(set([self.analytic_rank('cremona'), self.analytic_rank('ec'),
                     self.analytic_rank('rubinstein'), self.analytic_rank('sympow')]))
            if len(S) != 1:
                raise RuntimeError, "Bug in analytic rank; algorithms don't agree! (E=%s)"%E
            return S[0]
        else:
            raise ValueError, "algorithm %s not defined"%algorithm

    def p_isogenous_curves(self, p=None):
        r"""
        Return a list of pairs $(p, L)$ where $p$ is a prime and $L$
        is a list of the elliptic curves over $\Q$ that are
        $p$-isogenous to this elliptic curve.

        INPUT:
            p -- prime or None (default: None); if a prime, returns
                 a list of the p-isogenous curves.  Otherwise, returns
                 a list of all prime-degree isogenous curves sorted
                 by isogeny degree.

        This is implemented using Cremona's GP script \code{allisog.gp}.

        EXAMPLES:
            sage: E = EllipticCurve([0,-1,0,-24649,1355209])
            sage: E.p_isogenous_curves()
            [(2, [Elliptic Curve defined by y^2  = x^3 - x^2 - 91809*x - 9215775 over Rational Field, Elliptic Curve defined by y^2  = x^3 - x^2 - 383809*x + 91648033 over Rational Field, Elliptic Curve defined by y^2  = x^3 - x^2 + 1996*x + 102894 over Rational Field])]

        The isogeny class of the curve 11a2 has three curves in it.
        But \code{p_isogenous_curves} only returns one curves, since
        there is only one curve $5$-isogenous to 11a2.
            sage: E = EllipticCurve('11a2')
            sage: E.p_isogenous_curves()
            [(5, [Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field])]
            sage: E.p_isogenous_curves(5)
            [Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field]
            sage: E.p_isogenous_curves(3)
            []

        In contrast, the curve 11a1 admits two $5$-isogenies:
            sage: E = EllipticCurve('11a1')
            sage: E.p_isogenous_curves(5)
            [Elliptic Curve defined by y^2 + y = x^3 - x^2 - 7820*x - 263580 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 - x^2 over Rational Field]
        """
        if p is None:
            X = eval(gp_cremona.allisog(self.minimal_model().a_invariants()))
            Y = [(p, [constructor.EllipticCurve(ainvs) for ainvs in L]) for p, L in X]
            Y.sort()
            return Y
        else:
            X = eval(gp_cremona.p_isog(self.minimal_model().a_invariants(), p))
            Y = [constructor.EllipticCurve(ainvs) for ainvs in X]
            Y.sort()
            return Y

    def simon_two_descent(self, verbose=0, lim1=5, lim3=50, limtriv=10, maxprob=20, limbigprime=30):
        r"""
        Given a curve with no 2-torsion, computes (probably) the rank
        of the Mordell-Weil group, with certainty the rank of the
        2-Selmer group, and a list of independent points on the
        Weierstrass model of self.

        \note{The points are not translated back to self only because
        I haven't written code to do this yet.}

        INPUT:
            verbose -- integer, 0,1,2,3; (default: 0), the verbosity level
            lim1    -- (default: 5) limite des points triviaux sur les quartiques
            lim3    -- (default: 50) limite des points sur les quartiques ELS
            limtriv -- (default: 10) limite des points triviaux sur la
                                     courbe elliptique
            maxprob -- (default: 20)
            limbigprime -- (default: 30)  pour distinguer un petit 1nombre premier
                                     d'un grand utilise un test probabiliste pour
                                     les grands si LIMBIGPRIME = 0, n'utilise
                                     aucun test probabiliste

        OUTPUT:
            integer -- "probably" the rank of self
            integer -- the 2-rank of the Selmer group
            list    -- list of independent points on the Weierstrass model

        IMPLEMENTATION: Uses {\bf Denis Simon's} GP/PARI scripts from
                         \url{http://www.math.unicaen.fr/~simon/}

        EXAMPLES:
        We compute the ranks of the curves of lowest known conductor up to rank $8$.
        Amazingly, each of these computations finishes almost instantly!

            sage: E = EllipticCurve('11a1')
            sage: E.simon_two_descent()
            (0, 0, [])
            sage: E = EllipticCurve('37a1')
            sage: E.simon_two_descent()
            (1, 1, [(0 : 108 : 1)])
            sage: E = EllipticCurve('389a1')
            sage: E.simon_two_descent()
            (2, 2, [(57/4 : 621/8 : 1), (57 : 243 : 1)])
            sage: E = EllipticCurve('5077a1')
            sage: E.simon_two_descent()
            (3, 3, [(9 : 459 : 1), (153/4 : 189/8 : 1), (100 : 620 : 1)])


        In this example Simon's program does not find any points, though
        it does correctly compute the rank of the 2-Selmer group.
            sage: E = EllipticCurve([1, -1, 0, -751055859, -7922219731979])     # long (0.6 seconds)
            sage: E.simon_two_descent ()
            (1, 1, [])

        The rest of these entries were taken from Tom Womack's page
        \url{http://tom.womack.net/maths/conductors.htm}

            sage: E = EllipticCurve([1, -1, 0, -79, 289])
            sage: E.simon_two_descent()
            (4, 4, [(8415/49 : 10800/343 : 1), (-9 : 3672 : 1), (207 : 432 : 1), (-369 : 432 : 1)])
            sage: E = EllipticCurve([0, 0, 1, -79, 342])
            sage: E.simon_two_descent()        # random output
            (5, 5, [(0 : 3996 : 1), (-380 : 44 : 1), (52 : 3284 : 1), (110628/289 : 28166508/4913 : 1), (23364/25 : 3392388/125 : 1)])
            sage: E = EllipticCurve([1, 1, 0, -2582, 48720])
            sage: r, s, G = E.simon_two_descent(); r,s
            (6, 6)
            sage: E = EllipticCurve([0, 0, 0, -10012, 346900])
            sage: r, s, G = E.simon_two_descent(); r,s
            (7, 7)
            sage: E = EllipticCurve([0, 0, 1, -23737, 960366])
            sage: r, s, G = E.simon_two_descent(); r,s            # long time (1 second)
            (8, 8)
        """
        if self.torsion_order() % 2 == 0:
            raise ArithmeticError, "curve must not have rational 2-torsion\nThe *only* reason for this is that I haven't finished implementing the wrapper\nin this case.  It wouldn't be too difficult.\nPerhaps you could do it?!  Email me (wstein@gmail.com)."
        F = self.integral_weierstrass_model()
        a1,a2,a3,a4,a6 = F.a_invariants()
        t = simon_two_descent(a2,a4,a6, verbose=verbose, lim1=lim1, lim3=lim3, limtriv=limtriv,
                              maxprob=maxprob, limbigprime=limbigprime)
        prob_rank = rings.Integer(t[0])
        two_selmer_rank = rings.Integer(t[1])
        prob_gens = [F(P) for P in t[2]]
        return prob_rank, two_selmer_rank, prob_gens

    two_descent_simon = simon_two_descent

    def three_selmer_rank(self, bound=0, method=2):
        r"""
        Return the 3-selmer rank of this elliptic curve, computed
        using Magma.

        This is not implemented for all curves; a NotImplementedError
        exception is raised when this function is called on curves
        for which 3-descent isn't implemented.

        \note{Use a slightly modified version of Michael Stoll's MAGMA
        file \code{3descent.m}.  You must have Magma to use this
        function.}

        EXAMPLES:
            sage: EllipticCurve('37a').three_selmer_rank()  # optional & long -- Magma
            1

            sage: EllipticCurve('14a1').three_selmer_rank()      # optional
            Traceback (most recent call last):
            ...
            NotImplementedError:  Currently, only the case with irreducible phi3 is implemented.
        """
        import magma_3descent
        try:
            return magma_3descent.three_selmer_rank(self, bound, method)
        except RuntimeError, msg:
            msg = str(msg)
            i = msg.rfind(':')
            raise NotImplementedError, msg[i+1:]


    def rank(self, use_database=False, verbose=False,
                   only_use_mwrank=True,
                   algorithm='mwrank_shell'):
        """
        Return the rank of this elliptic curve, assuming no conjectures.

        If we fail to provably compute the rank, raises a RuntimeError
        exception.

        INPUT:
            use_database (bool) -- (default: False), if True, try to
                  look up the regulator in the Cremona database.
            verbose -- (default: None), if specified changes the
                       verbosity of mwrank computations.
            algorithm -- 'mwrank_shell' -- call mwrank shell command
                      -- 'mwrank_lib' -- call mwrank c library
            only_use_mwrank -- (default: True) if False try using
                       analytic rank methods first.

        OUTPUT:
            rank (int) -- the rank of the elliptic curve.

        IMPLEMENTATION: Uses L-functions and mwrank.
        """
        try:
            return self.__rank
        except AttributeError:
            if use_database:
                try:
                    self.__rank = self.database_curve().rank()
                    return self.__rank
                except (AttributeError, RuntimeError):
                    pass
            if not only_use_mwrank:
                N = self.conductor()
                prec = int(4*float(sqrt(N))) + 10
                if self.root_number() == 1:
                    L, err = self.Lseries_at1(prec)
                    if abs(L) > err + R(0.0001):  # definitely doesn't vanish
                        misc.verbose("rank 0 because L(E,1)=%s"%L)
                        self.__rank = 0
                        return self.__rank
                else:
                    Lprime, err = self.Lseries_deriv_at1(prec)
                    if abs(Lprime) > err + R(0.0001):  # definitely doesn't vanish
                        misc.verbose("rank 1 because L'(E,1)=%s"%Lprime)
                        self.__rank = 1
                        return self.__rank

            if algorithm == 'mwrank_lib':
                misc.verbose("using mwrank lib")
                C = self.mwrank_curve()
                C.set_verbose(verbose)
                r = C.rank()
                if not C.certain():
                    del self.__mwrank_curve
                    raise RuntimeError, "Unable to compute the rank with certainty (lower bound=%s).  This could be because Sha(E/Q)[2] is nontrivial."%C.rank() + "\nTrying calling something like two_descent(second_limit=13) on the curve then trying this command again.  You could also try rank with only_use_mwrank=False."
                self.__rank = r
            elif algorithm == 'mwrank_shell':
                misc.verbose("using mwrank shell")
                X = self.mwrank()
                if not 'The rank and full Mordell-Weil basis have been determined unconditionally' in X:
                    raise RuntimeError, '%s\nRank not provably correct (maybe try rank with only_use_mwrank=False).'%X
                i = X.find('Rank = ')
                assert i != -1
                j = i + X[i:].find('\n')
                self.__rank = Integer(X[i+7:j])
        return self.__rank

    def gens(self, verbose=False, rank1_search=10,
             algorithm='mwrank_shell',
             only_use_mwrank=True,
             proof = True):
        """
        Compute and return generators for the Mordell-Weil group E(Q)
        *modulo* torsion.

        HINT: If you would like to control the height bounds used
        in the 2-descent, first call the two_descent function with
        those height bounds.

        TODO: Right now this function assumes that the input curve is
        in minimal Weierstrass form.  This restriction will be removed
        in the near future.  This function raises a
        NotImplementedError if a non-minimal curve is given as input.

        WARNING: If the program fails to give a provably correct
        result, it prints a warning message, but does not raise an
        exception.  Use the gens_certain command to find out if
        this warning message was printed.

        INPUT:
            verbose -- (default: None), if specified changes the
                       verbosity of mwrank computations.
            rank1_search -- (default: 16), if the curve has analytic
                       rank 1, try to find a generator by a direct
                       search up to this logarithmic height.  If this
                       fails the usual mwrank procedure is called.
            algorithm -- 'mwrank_shell' (default) -- call mwrank shell command
                      -- 'mwrank_lib' -- call mwrank c library
        OUTPUT:
            generators -- List of generators for the Mordell-Weil group.

        IMPLEMENTATION: Uses Cremona's mwrank C library.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: E.gens()                 # random output
            [(-1 : 1 : 1), (0 : 0 : 1)]
        """
        try:
            return list(self.__gens)  # return copy so not changed
        except AttributeError:
            pass
        if self.conductor() > 10**7:
            only_use_mwrank = True

        if not only_use_mwrank:
            try:
                misc.verbose("Trying to compute rank.")
                r = self.rank()
                misc.verbose("Got r = %s."%r)
                if r == 0:
                    misc.verbose("Rank = 0, so done.")
                    self.__gens = []
                    self.__regulator = R(1)
                    self.__gens_certain = True
                    return self.__gens
                if r == 1 and rank1_search:
                    misc.verbose("Rank = 1, so using direct search.")
                    h = 6
                    while h <= rank1_search:
                        misc.verbose("Trying direct search up to height %s"%h)
                        G = self.point_search(h, verbose)
                        G = [P for P in G if P.order() == oo]
                        if len(G) > 0:
                            misc.verbose("Direct search succeeded.")
                            G, _, reg = self.saturation(G, verbose=verbose)
                            misc.verbose("Computed saturation.")
                            self.__gens = G
                            self.__gens_certain = True
                            self.__regulator = reg
                            return self.__gens
                        h += 2
                    misc.verbose("Direct search FAILED.")
            except RuntimeError:
                pass
        # end if (not_use_mwrank)
        if not self.is_integral():
            raise NotImplementedError, "gens via mwrank only implemented for curves with integer coefficients."
        if algorithm == "mwrank_lib":
            misc.verbose("Calling mwrank C++ library.")
            C = self.mwrank_curve(verbose)
            if not (verbose is None):
                C.set_verbose(verbose)
            G = C.gens()
            self.__gens_certain = C.certain()
            if not self.__gens_certain:
                del self.__mwrank_curve
                raise RuntimeError, "Unable to compute the rank, hence generators, with certainty (lower bound=%s).  This could be because Sha(E/Q)[2] is nontrivial."%C.rank() + \
                      "\nTrying calling something like two_descent(second_limit=13) on the curve then trying this command again."
        else:
            X = self.mwrank()
            misc.verbose("Calling mwrank shell.")
            if not 'The rank and full Mordell-Weil basis have been determined unconditionally' in X:
                msg = 'Generators not provably computed.'
                if proof:
                    raise RuntimeError, '%s\n%s'%(X,msg)
                else:
                    misc.verbose("Warning -- %s"%msg, level=0)
            G = []
            i = X.find('Generator ')
            while i != -1:
                j = i + X[i:].find(';')
                k = i + X[i:].find('[')
                G.append(eval(X[k:j].replace(':',',')))
                X = X[j:]
                i = X.find('Generator ')
            i = X.find('Regulator = ')
            j = i + X[i:].find('\n')
            self.__regulator = R(X[i+len('Regulator = '):j])
        ####
        self.__gens = [self.point(x, check=True) for x in G]
        self.__gens.sort()
        self.__rank = len(self.__gens)
        return self.__gens

    def gens_certain(self):
        """
        Return True if the generators have been proven correct.
        """
        try:
            return self.__gens_certain
        except AttributeError:
            self.gens()
            return self.__gens_certain

    def ngens(self):
        return len(self.gens())

    def regulator(self, use_database=True, verbose=None):
        """
        Returns the regulator of this curve, which must be defined
        over Q.

        INPUT:
            use_database -- bool (default: False), if True, try to
                  look up the regulator in the Cremona database.
            verbose -- (default: None), if specified changes the
                  verbosity of mwrank computations.

        EXAMPLES:
            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: E.regulator()                           # long time (1 second)
            0.051111408239968799
        """
        try:
            return self.__regulator
        except AttributeError:
            if use_database:
                try:
                    self.__regulator = R(self.database_curve().db_extra[3])
                    return self.__regulator
                except (AttributeError, RuntimeError):
                    pass
            G = self.gens()
            try:  # in some cases self.gens() efficiently computes regulator.
                return self.__regulator
            except AttributeError:
                pass
            C = self.mwrank_curve()
            reg = R(C.regulator())
            if not C.certain():
                raise RuntimeError, "Unable to compute the rank, hence regulator, with certainty (lower bound=%s)."%C.rank()
            self.__regulator = reg
        return self.__regulator

    def saturation(self, points, verbose=False, max_prime=0, odd_primes_only=False):
        """
        Given a list of rational points on E, compute the saturation
        in E(Q) of the subgroup they generate.

        INPUT:
            points (list) -- list of points on E
            verbose  (bool) -- (default: False), if True, give verbose output
            max_prime (int) -- (default: 0), saturation is performed
                               for all primes up to max_prime.  If max_prime==0,
                               perform saturation at *all* primes, i.e., compute
                               the true saturation.
            odd_primes_only (bool) -- only do saturation at odd primes

        OUTPUT:
            saturation (list) -- points that form a basis for the saturation
            index (int) -- the index of the group generated by points in their saturation
            regulator (float) -- regulator of saturated points.

        IMPLEMENTATION: Uses Cremona's mwrank package.  With max_prime=0, we call
            mwrank with successively larger prime bounds until the full saturation is
            provably found.  The results of saturation at the previous primes is stored
            in each case, so this should be reasonably fast.
        """
        if not isinstance(points, list):
            raise TypeError, "points (=%s) must be a list."%points

        v = []
        for P in points:
            if not isinstance(P, ell_point.EllipticCurvePoint_field):
                P = self(P)
            elif P.curve() != self:
                raise ArithmeticError, "point (=%s) must be %s."%(P,self)
            x, y = P.xy()
            d = x.denominator().lcm(y.denominator())
            v.append((x*d, y*d, d))
        c = self.mwrank_curve()
        mw = mwrank.mwrank_MordellWeil(c, verbose)
        mw.process(v)
        if max_prime == 0:
            repeat_until_saturated = True
            max_prime = 97
        while True:
            ok, index, unsat = mw.saturate(max_prime=max_prime, odd_primes_only = odd_primes_only)
            reg = mw.regulator()
            if not ok and repeat_until_saturated:
                max_prime = arith.next_prime(max_prime + 100)
                ok, index, unsat = mw.saturate(max_prime=max_prime, odd_primes_only = odd_primes_only)
                reg = mw.regulator()
            else:
                break
        sat = mw.points()
        sat = [self(P) for P in sat]
        return sat, index, R(reg)

    def CPS_height_bound(self):
        """
        Return the Cremona-Prickett-Siksek height bound.  This is a
        floating point number B such that if P is a point on the curve,
        then the naive logarithmetic height of P is off from the
        canonical height by at most B.

        EXAMPLES:
            sage: E = EllipticCurve("11a")
            sage: E.CPS_height_bound()
            2.8774743273580445
            sage: E = EllipticCurve("5077a")
            sage: E.CPS_height_bound()
            0.0
            sage: E = EllipticCurve([1,2,3,4,1])
            sage: E.CPS_height_bound()
            Traceback (most recent call last):
            ...
            RuntimeError: curve must be minimal.
            sage: F = E.quadratic_twist(-19)
            sage: F
            Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 + 1376*x - 130 over Rational Field
            sage: F.CPS_height_bound()
            0.65551583769728516

        IMPLEMENTATION:
            Call the corresponding mwrank C++ library function.
        """
        if not self.is_minimal():
            raise RuntimeError, "curve must be minimal."
        return self.mwrank_curve().CPS_height_bound()


    def silverman_height_bound(self):
        """
        Return the Silverman height bound.  This is a floating point
        number B such that if P is a point on the curve, then the
        naive logarithmetic height of P is off from the canonical
        height by at most B.

        Note that the CPS_height_bound is typically much better than
        the Silverman bound.
        """
        return self.mwrank_curve().silverman_bound()


    def point_search(self, height_limit, verbose=True):
        """
        Search for points on a curve up to an input bound on the naive logarithmic height.

        INPUT:
            height_limit (float) -- bound on naive height (at most 21, or mwrank overflows)
            verbose  (bool) -- (default: True)

        OUTPUT:
            points (list) -- list of points found

        IMPLEMENTATION: Uses Cremona's mwrank package.
        """
        height_limit = float(height_limit)
        if height_limit > _MAX_HEIGHT:
            raise OverflowError, "height_limit (=%s) must be at most %s."%(height_limit,_MAX_HEIGHT)
        c = self.mwrank_curve()
        mw = mwrank.mwrank_MordellWeil(c, verbose)
        mw.search(height_limit, verbose=verbose)
        v = mw.points()
        return [self(P) for P in v]

    def two_torsion_rank(self):
        r"""
        Return the dimension of the 2-torsion subgroup of $E(\Q)$.

        EXAMPLES:
        """
        A = self.torsion_subgroup().invariants()
        if len(A) == 2:
            return rings.Integer(2)
        elif len(A) == 1 and A[0] % 2 == 0:
            return rings.Integer(1)
        else:
            return rings.Integer(0)

    def selmer_rank_bound(self):
        """
        Bound on the rank of the curve, computed using the
        2-selmer group.  This is the rank of the curve
        minus the rank of the 2-torsion, minus a number
        determined by whatever mwrank was able to determine
        related to Sha[2].  Thus in many cases, this is
        the actual rank of the curve.

        EXAMPLE:
        The following is the curve 960D1, which has rank 0,
        but Sha of order 4.
            sage: E = EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.selmer_rank_bound()
            0

        It gives 0 instead of 2, because it knows Sha is nontrivial.
        In contrast, for the curve 571A, also with rank 0 and Sha
        of order 4, we get a worse bound:
            sage: E = EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.selmer_rank_bound()
            2
            sage: E.rank(only_use_mwrank=False)   # uses L-function
            0
        """
        try:
            return self.__selmer_rank_bound
        except AttributeError:
            C = self.mwrank_curve()
            self.__selmer_rank_bound = C.selmer_rank_bound()
            return self.__selmer_rank_bound


    def an(self, n):
        """
        The n-th Fourier coefficient of the modular form corresponding
        to this elliptic curve, where n is a positive integer.
        """
        return Integer(self.pari_mincurve().ellak(n))

    def ap(self, p):
        """
        The p-th Fourier coefficient of the modular form corresponding
        to this elliptic curve, where p is prime.
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        return Integer(self.pari_mincurve().ellap(p))

    def quadratic_twist(self, D):
        return EllipticCurve_field.quadratic_twist(self, D).minimal_model()

    def minimal_model(self):
        r"""
        Return the unique minimal Weierstrass equation for this
        elliptic curve.  This is the model with minimal discriminant
        and $a_1,a_2,a_3 \in \{0,\pm 1\}$.
        """
        try:
            return self.__minimal_model
        except AttributeError:
            F = self.pari_mincurve()
            self.__minimal_model = EllipticCurve_rational_field([Q(F[i]) for i in range(5)])
            return self.__minimal_model

    def is_minimal(self):
        return self.ainvs() == self.minimal_model().ainvs()

    def is_integral(self):
        for n in self.ainvs():
            if n.denominator() != 1:
                return False
        return True

    def is_isomorphic(self, E):
        if not isinstance(E, EllipticCurve_rational_field):
            raise TypeError, "E (=%s) must be an elliptic curve over the rational numbers"%E
        return E.minimal_model() == self.minimal_model()

    def kodaira_type(self, p):
        """
        Local Kodaira type of the elliptic curve at $p$.

        1 means good reduction (type $I_0$), 2, 3 and 4 mean types II,
        III and IV, respectively, $4 + \\nu$ with $\\nu > 0$ means
        type $I_{\\nu}$; finally the opposite values -1, -2,
        etc. refer to the starred types $I_0^*$, $II^*$, etc.

        EXAMPLES:
            sage: E = EllipticCurve('124a')
            sage: E.kodaira_type(2)
            '4'
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        try:
            self.__kodaira_type
        except AttributeError:
            self.__kodaira_type = {}
            self.__tamagawa_number = {}
        if not self.__kodaira_type.has_key(p):
            v = self.pari_mincurve().elllocalred(p)
            self.__kodaira_type[p] = str(v[1])
            self.__tamagawa_number[p] = int(v[3])
        return self.__kodaira_type[p]

    def tamagawa_number(self, p):
        """
        The Tamagawa number of the elliptic curve at $p$.

        EXAMPLES:
            sage: E = EllipticCurve('11a')
            sage: E.tamagawa_number(11)
            5
            sage: E = EllipticCurve('37b')
            sage: E.tamagawa_number(37)
            3
        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        try:
            return self.__tamagawa_number[p]
        except (AttributeError, KeyError):
            self.kodaira_type(p)
            return self.__tamagawa_number[p]

    def tamagawa_numbers(self):
        return [self.tamagawa_number(p) for p in arith.prime_divisors(self.conductor())]

    def tamagawa_product(self):
        """
        Returns the product of the Tamagawa numbers.

        EXAMPLES:
            sage: E = EllipticCurve('54a')
            sage: E.tamagawa_product ()
            3
        """
        try:
            return self.__tamagawa_product
        except AttributeError:
            self.__tamagawa_product = self.pari_mincurve().ellglobalred()[2].python()
            return self.__tamagawa_product

    def real_components(self):
        """
        Returns 1 if there is 1 real component and 2 if there are 2.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.real_components ()
            2
            sage: E = EllipticCurve('37b')
            sage: E.real_components ()
            2
            sage: E = EllipticCurve('11a')
            sage: E.real_components ()
            1
        """
        invs = self.weierstrass_model().ainvs()
        x = rings.polygen(self.base_ring())
        f = x**3 + invs[3]*x + invs[4]
        if f.discriminant() > 0:
            return 2
        else:
            return 1

    def period_lattice(self):
        r"""
        Return a basis for the period lattice of the elliptic curve
        over $\Q$ as a 2-tuple.

        The basis has the form $[\omega_1, \omega_2]$, where
        $\Im(\omega_1/\omega_2) > 0$ and $\omega_1$ is real.

        TODO: The precision is determined by the state of the PARI C
        library, which is not good.

        INPUT:
            -- an elliptic curve
        OUTPUT:
            omega_1 -- complex number
            omega_2 -- complex number

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.period_lattice ()
            (2.993458646231959629832009979452508177797583791370132985340523378563250356987, 2.451389381986790060854224831866525225349617289144796614656471406129152899999*I)    # 32-bit
            (2.99345864623195962983200997945250817779758379137013298534052337856325035698668290412039406739705147343584052710494728819414438723737202525437537667109326, 2.45138938198679006085422483186652522534961728914479661465647140612915289999925689289113212802918108871268421886966184797547519986661675580167893816478303*I)   # 64-bit
        """
        return tuple(self.pari_curve().omega().python())

    def omega(self):
        """
        Returns the real period.

        If self is given by a \emph{minimal Weierstrass equation} then
        this is the correct period in the BSD conjecture, i.e., it is
        the least real period * 2 when the period lattice is
        rectangular.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.omega()
            5.986917292463919259664019958905016355595167582740265970681046757126500713973     # 32-bit
            5.98691729246391925966401995890501635559516758274026597068104675712650071397336580824078813479410294687168105420989457638828877447474405050875075334218652        # 64-bit


        This is not a minimal model.
            sage: E = EllipticCurve([0,-432*6^2])
            sage: E.omega()
            0.4861093857100564298972304561738255425526098601923921971195408561181781048715    # 32-bit
            0.486109385710056429897230456173825542552609860192392197119540856118178104871498709353307487730842084963451161261340032305532890657753313985258848453458110       # 64-bit

        If you were to plug the above omega into the BSD conjecture, you
        would get nonsense.   The following works though:
            sage: F = E.minimal_model()
            sage: F.omega()
            0.9722187714201128597944609123476510851052197203847843942390817122363562097430    # 32-bit
             0.972218771420112859794460912347651085105219720384784394239081712236356209742997418706614975461684169926902322522680064611065781315506627970517696906916220      # 64-bit
        """
        return self.period_lattice()[0] * self.real_components()

    def complex_area(self):
        """
        Return the area of a fundamental domain for the period lattice
        of the elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.complex_area()
            7.338132740789576739070721003332305588006176586123733717543180156079096606979     # 32-bit
            7.33813274078957673907072100333230558800617658612373371754318015607909660697945809438214607592923817142863798604406024044503049908233884534256274529672707        # 64-bit
        """
        w1,w2 = self.period_lattice()
        return (w1*w2.imag()).real()

    def Lseries_dokchitser(self, prec=53,
                           max_imaginary_part=0,
                           max_asymp_coeffs=40,
                           algorithm='gp'):
        r"""
        Return interface to Tim Dokchitser's program for computing
        with the L-series of this elliptic curve; this provides a way
        to compute Taylor expansions and higher derivatives of
        $L$-series.

        INPUT:
            prec -- integer (bits precision)
            max_imaginary_part -- real number
            max_asymp_coeffs -- integer
            algorithm -- string: 'gp' or 'magma'

        \note{If algorithm='magma', then the precision is in digits rather
        than bits and the object returned is a Magma L-series, which has
        different functionality from the SAGE L-series.}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: L = E.Lseries_dokchitser()
            sage: L(2)
            0.381575408260711
            sage.: L = E.Lseries_dokchitser(algorithm='magma')         # optional  (not auto tested)
            sage.: L.Evaluate(2)                                       # optional  (not auto tested)
            0.38157540826071121129371040958008663667709753398892116
        """
        if algorithm == 'magma':
            from sage.interfaces.all import magma
            return magma(self).LSeries(Precision = prec)

        from sage.lfunctions.all import Dokchitser
        L = Dokchitser(conductor = self.conductor(),
                       gammaV = [0,1],
                       weight = 2,
                       eps = self.root_number(),
                       poles = [],
                       prec = prec)
        gp = L.gp()
        s = 'e = ellinit(%s);'%self.minimal_model().a_invariants()
        s += 'a(k) = ellak(e, k);'
        L.init_coeffs('a(k)', 1, pari_precode = s,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        L.rename('Dokchitser L-function associated to %s'%self)
        return L

    def Lseries_sympow(self, n, prec):
        r"""
        Return $L(\Sym^{(n)}(E, \text{edge}))$ to prec digits
        of precision.

        INPUT:
            n -- integer
            prec -- integer

        OUTPUT:
            string -- real number to prec digits of precision as a string.

        \note{Before using this function for the first time for
        a given $n$, you may have to type \code{sympow('-new_data <n>')},
        where \code{<n>} is replaced by your value of $n$.  This
        command takes a long time to run.}

        EXAMPLES:
            sage.: E = EllipticCurve('37a')
            sage.: a = E.Lseries_sympow(2,16); a
            '2.492262044273650E+00'
            sage.: RR(a)
            2.4922620442736498
        """
        from sage.lfunctions.sympow import sympow
        return sympow.L(self, n, prec)

    def Lseries_sympow_derivs(self, n, prec, d):
        r"""
        Return $0$th to $d$th derivatives of $L(\Sym^{(n)}(E,
        \text{edge}))$ to prec digits of precision.

        INPUT:
            n -- integer
            prec -- integer
            d -- integer

        OUTPUT:
            a string, exactly as output by sympow

        \note{To use this function you may have to run a few commands
        like \code{sympow('-new_data 1d2')}, each which takes a few
        minutes.  If this function fails it will indicate what
        commands have to be run.}

        EXAMPLES:
            sage.: E = EllipticCurve('37a')
            sage.: E.Lseries_sympow_derivs(1,16,2)
            ...
            1n0: 3.837774351482055E-01
            1w0: 3.777214305638848E-01
            1n1: 3.059997738340522E-01
            1w1: 3.059997738340524E-01
            1n2: 1.519054910249753E-01
            1w2: 1.545605024269432E-01
        """
        from sage.lfunctions.sympow import sympow
        return sympow.Lderivs(self, n, prec, d)

    def Lseries_zeros(self, n):
        """
        Return the imaginary parts of the first $n$ nontrivial zeros
        on the critical line of the L-function in the upper half
        plane, as 32-bit reals.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.Lseries_zeros(2)
            [0.000000000, 5.00317001]

            sage: a = E.Lseries_zeros(20)      # long time
            sage: point([(1,x) for x in a])    # graph  (long time)

        AUTHOR:
            -- Uses Rubinstein's L-functions calculator.
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.zeros(n, L=self)

    def Lseries_zeros_in_interval(self, x, y, stepsize):
        r"""
        Return the imaginary parts of (most of) the nontrivial zeros
        on the critical line $\Re(s)=1$ with positive imaginary part
        between $x$ and $y$, along with a technical quantity for each.

        INPUT:
            x, y, stepsize -- positive floating point numbers

        OUTPUT:
            list of pairs (zero, S(T)).

        Rubinstein writes: The first column outputs the imaginary part
        of the zero, the second column a quantity related to S(T) (it
        increases roughly by 2 whenever a sign change, i.e. pair of
        zeros, is missed). Higher up the critical strip you should use
        a smaller stepsize so as not to miss zeros.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.Lseries_zeros_in_interval(6, 10, 0.1)      # long
            [(6.87039122, 0.248922780), (8.01433081, -0.140168533), (9.93309835, -0.129943029)]
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.zeros_in_interval(x, y, stepsize, L=self)

    def Lseries_values_along_line(self, s0, s1, number_samples):
        """
        Return values of $L(E, s)$ at \code{number_samples}
        equally-spaced sample points along the line from $s_0$ to
        $s_1$ in the complex plane.

        \note{The L-series is normalized so that the center of the
        critical strip is 1.}

        INPUT:
            s0, s1 -- complex numbers
            number_samples -- integer

        OUTPUT:
            list -- list of pairs (s, zeta(s)), where the s are
                    equally spaced sampled points on the line from
                    s0 to s1.

        EXAMPLES:
            sage: I = CC.0
            sage: E = EllipticCurve('37a')
            sage: E.Lseries_values_along_line(1, 0.5+20*I, 5)     # long time
            [(0.50000000000, 0), (0.40000000002 + 4.0000000000*I, 3.3192024464 - 2.6002805391*I), (0.30000000005 + 8.0000000000*I, -0.88634118531 - 0.42264033738*I), (0.20000000001 + 12.000000000*I, -3.5055893594 - 0.10853169035*I), (0.10000000001 + 16.000000000*I, -3.8704328826 - 1.8804941061*I)]
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.values_along_line(s0-RationalField()('1/2'),
                                       s1-RationalField()('1/2'),
                                       number_samples, L=self)

    def Lseries_twist_values(self, s, dmin, dmax):
        r"""
        Return values of $L(E, s, \chi_d)$ for each quadratic
        character $\chi_d$ for $d_{\min} \leq d \leq d_{\max}$.

        \note{The L-series is normalized so that the center of the
        critical strip is 1.}

        INPUT:
            s -- complex numbers
            dmin -- integer
            dmax -- integer

        OUTPUT:
            list -- list of pairs (d, L(E, s,chi_d))

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.Lseries_twist_values(1, -12, -4)    # slightly random output depending on architecture
            [(-11, 1.4782434171), (-8, 0), (-7, 1.8530761916), (-4, 2.4513893817)]
            sage: F = E.quadratic_twist(-8)
            sage: F.rank()
            1
            sage: F = E.quadratic_twist(-7)
            sage: F.rank()
            0
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.twist_values(s - RationalField()('1/2'), dmin, dmax, L=self)

    def Lseries_twist_zeros(self, n, dmin, dmax):
        r"""
        Return first $n$ real parts of nontrivial zeros of
        $L(E,s,\chi_d)$ for each quadratic character $\chi_d$ with
        $d_{\min} \leq d \leq d_{\max}$.

        \note{The L-series is normalized so that the center of the
        critical strip is 1.}

        INPUT:
            n -- integer
            dmin -- integer
            dmax -- integer

        OUTPUT:
            dict -- keys are the discriminants $d$, and
                    values are list of corresponding zeros.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.Lseries_twist_zeros(3, -4, -3)         # long
            {-4: [1.60813783, 2.96144840, 3.89751747], -3: [2.06170900, 3.48216881, 4.45853219]}
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.twist_zeros(n, dmin, dmax, L=self)

    def Lseries_at1(self, k=0):
        r"""
        Compute $L(E,1)$ using $k$ terms of the series for $L(E,1)$ as
        explained on page 406 of Henri Cohen's book"A Course in Computational
        Algebraic Number Theory".  If the argument $k$ is not specified,
        then it defaults to $\sqrt(N)$, where $N$ is the conductor.

        The real precision used in each step of the computation is the
        precision of machine floats.

        INPUT:
            k -- (optional) an integer, defaults to sqrt(N).

        OUTPUT:
            float -- L(E,1)
            float -- a bound on the error in the approximation; this
                     is a proveably correct upper bound on the sum
                     of the tail end of the series used to compute L(E,1).

        This function is disjoint from the PARI \code{elllseries}
        command, which is for a similar purpose.  To use that command
        (via the PARI C library), simply type
                \code{E.pari_mincurve().elllseries(1)}

        ALGORITHM:
        \begin{enumerate}
            \item Compute the root number eps.  If it is -1, return 0.

            \item Compute the Fourier coefficients a_n, for n up to and
               including k.

            \item Compute the sum
            $$
                 2 * sum_{n=1}^{k} (a_n / n) * exp(-2*pi*n/Sqrt(N)),
            $$
               where N is the conductor of E.

            \item Compute a bound on the tail end of the series, which is
            $$
                 2 * e^(-2 * pi * (k+1) / sqrt(N)) / (1 - e^(-2*pi/sqrt(N))).
            $$
               For a proof see [Grigov-Jorza-Patrascu-Patrikis-Stein].
        \end{enumerate}

        EXAMPLES:
            sage: E = EllipticCurve('37b')
            sage: E.Lseries_at1(100)
            (0.725681061936153, 1.52437502288743e-45)
        """
        if self.root_number() == -1:
            return 0
        sqrtN = float(self.conductor().sqrt())
        k = int(k)
        if k == 0: k = int(math.ceil(sqrtN))
        an = self.anlist(k)           # list of SAGE ints
        # Compute z = e^(-2pi/sqrt(N))
        pi = 3.14159265358979323846
        z = exp(-2*pi/sqrtN)
        zpow = z
        s = 0.0
        for n in xrange(1,k+1):
            s += (zpow * float(an[n]))/n
            zpow *= z

        error = 2*zpow / (1 - z)

        return R(2*s), R(error)

    def Lseries_deriv_at1(self, k=0):
        r"""
        Compute $L'(E,1)$ using$ k$ terms of the series for $L'(E,1)$.

        The algorithm used is from page 406 of Henri Cohen's book ``A
        Course in Computational Algebraic Number Theory.''

        The real precision of the computation is the precision of
        Python floats.

        INPUT:
            k -- int; number of terms of the series

        OUTPUT:
            real number -- an approximation for L'(E,1)
            real number -- a bound on the error in the approximation

        ALGORITHM:
        \begin{enumerate}
            \item Compute the root number eps.  If it is 1, return 0.

            \item Compute the Fourier coefficients $a_n$, for $n$ up to and
                  including $k$.

            \item Compute the sum
               $$
                 2 * \sum_{n=1}^{k} (a_n / n) * E_1(2 \pi n/\sqrt{N}),
               $$
               where $N$ is the conductor of $E$, and $E_1$ is the
               exponential integral function.

            \item Compute a bound on the tail end of the series, which is
               $$
                 2 * e^{-2 \pi (k+1) / \sqrt{N}} / (1 - e^{-2 \ pi/\sqrt{N}}).
               $$
               For a proof see [Grigorov-Jorza-Patrascu-Patrikis-Stein].  This
               is exactly the same as the bound for the approximation to
               $L(E,1)$ produced by \code{Lseries_at1}.
        \end{enumerate}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.Lseries_deriv_at1(100)
            (0.305999773834879, 1.52437502288740e-45)
        """
        if self.root_number() == 1: return 0
        k = int(k)
        sqrtN = float(self.conductor().sqrt())
        if k == 0: k = int(math.ceil(sqrtN))
        an = self.anlist(k)           # list of C ints
        # Compute z = e^(-2pi/sqrt(N))
        pi = 3.14159265358979323846
        v = transcendental.exponential_integral_1(2*pi/sqrtN, k)
        L = 2*float(sum([ (v[n-1] * an[n])/n for n in xrange(1,k+1)]))
        error = 2*exp(-2*pi*(k+1)/sqrtN)/(1-exp(-2*pi/sqrtN))
        return R(L), R(error)

    def Lseries(self, s):
        r"""
        Returns the value of the L-series of the elliptic curve E at s, where s
        must be a real number.

        Use self.Lseries_extended for s complex.

        \note{If the conductor of the curve is large, say $>10^{12}$,
        then this function will take a very long time, since it uses
        an $O(\sqrt{N})$ algorithm.}

        EXAMPLES:
            sage: E = EllipticCurve([1,2,3,4,5])
            sage: E.Lseries(1)
            0.000000000000000
            sage: E.Lseries('1.1')       # long time (!)
            0.28549100767814833

        TODO: Planned massive improvement -- use Micheal Rubinstein's
        L-functions package and/or Tim Dokchitser's.  Both are already
        available via other function calls.  Note that Dokchitser's
        program is vastly faster than PARI, e.g., at computing
        E.Lseries(1.1) above, even with all the startup overhead, etc,
        e.g., 10 seconds versus 0.25 seconds.
        """
        try:
            s = R(s)
        except TypeError:
            raise TypeError, "for non-real input, use self.Lseries_extended instead."
        if s <= 0 and s.frac() == 0:
            # The L-series vanishes at negative integers, but PARI
            # is broken for this.
            return R(0)
        return R(self.pari_mincurve().elllseries(s))

    def Lambda(self, s, prec):
        r"""
        Returns the value of the Lambda-series of the elliptic curve E
        at s, where s can be any complex number.

        IMPLEMENTATION: Fairly *slow* computation using the definitions
        and implemented in Python.

        Uses prec terms of the power series.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: E.Lambda(1.4+0.5*I, 50)
            -0.354172680517671 + 0.874518681720170*I
        """
        s = C(s)
        N = self.conductor()
        pi = R(constants.pi)
        Gamma = transcendental.gamma
        Gamma_inc = transcendental.gamma_inc
        a = self.anlist(prec)
        eps = self.root_number()
        sqrtN = float(N.sqrt())
        def F(n, t):
            return Gamma_inc(t+1, 2*pi*n/sqrtN) * C(sqrtN/(2*pi*n))**(t+1)
        return sum([a[n]*(F(n,s-1) + eps*F(n,1-s)) for n in xrange(1,prec+1)])

    def Lseries_extended(self, s, prec):
        r"""
        Returns the value of the L-series of the elliptic curve E at s
        can be any complex number using prec terms of the power series
        expansion.


        WARNING: This may be slow.  Consider using \code{Lseries_dokchitser()}
        instead.

        INPUT:
            s -- complex number
            prec -- integer

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: E.Lseries_extended(1 + I, 50)
            -0.638409959099589 + 0.715495262192901*I
            sage: E.Lseries_extended(1 + 0.1*I, 50)
            -0.00761216538818315 + 0.000434885704670107*I

        NOTE: You might also want to use Tim Dokchitser's
        L-function calculator, which is available by typing
        L = E.Lseries_dokchitser(), then evaluating L.  It
        gives the same information but is sometimes much faster.

        """
        try:
            s = C(s)
        except TypeError:
            raise TypeError, "Input argument %s must be coercible to a complex number"%s
        prec = int(prec)
        if abs(s.imag()) < R(0.0000000000001):
            return self.Lseries(s.real())
        N = self.conductor()
        pi = R(constants.pi)
        Gamma = transcendental.gamma
        Gamma_inc = transcendental.gamma_inc
        a = self.anlist(prec)
        eps = self.root_number()
        sqrtN = float(N.sqrt())
        def F(n, t):
            return Gamma_inc(t+1, 2*pi*n/sqrtN) * C(sqrtN/(2*pi*n))**(t+1)
        return C(N)**(-s/2) * C(2*pi)**s * Gamma(s)**(-1)\
               * sum([a[n]*(F(n,s-1) + eps*F(n,1-s)) for n in xrange(1,prec+1)])

    def sigma(self, z, flag=0):
        """
        Returns the value of the Weierstrass sigma function of the lattice
        associated to this elliptic curve E.

        INPUT:
            z -- a complex number
            flag -- 0 - default ???
                    1 - computes an arbitrary determination of log(sigma(z))
                    2, 3 - same using the product expansion instead of theta series.
                           ???
        OUTPUT:
            a complex number

        NOTE: The reason for the ???'s above, is that the PARI documentation for
              ellsigma is very vague.
        """
        return self.pari_curve().ellsigma(z, flag)

    def weierstrass_model(self):
        r"""
        Return a model of the form $y^2 = x^3 + a*x + b$ for this curve.

        More precisely, we have $a = c_4 / (2^4 \cdot 3)$ and
        $b = -c_6 / (2^5\cdot 3^3)$, where $c_4, c_6$ are the $c$-invariants
        for a minimal Weierstrass equation for $E$.

        Use \code{self.integral_weierstrass_model()} for a model with
        $a,b\in\ZZ$.
        """
        F = self.minimal_model()
        return EllipticCurve_field.weierstrass_model(F)

    def integral_weierstrass_model(self):
        """
        Return a model of the form $y^2 = x^3 + a*x + b$ for this curve with $a,b\in\Z$.
        """
        F = self.minimal_model()
        a0, a1, a2, a3, a4 = F.ainvs()
        return constructor.EllipticCurve([-27*a0**4 - 216*a0**2*a1 + 648*a0*a2 - 432*a1**2 + 1296*a3, \
                                          54*a0**6 + 648*a0**4*a1 - 1944*a0**3*a2 + 2592*a0**2*a1**2 -\
                                          3888*a0**2*a3 - 7776*a0*a1*a2 + 3456*a1**3 - \
                                          15552*a1*a3 + 11664*a2**2 + 46656*a4])

    def watkins_data(self):
        """
        Return a dict of the data computed by Mark Watkins's ec
        program applied to this elliptic curve.
        """
        try:
            return self.__watins_data
        except AttributeError:
            try:
                import sage.libs.ec.all
            except ImportErrror:
                raise NotImplementedError
            self.__watkins_data = sage.libs.ec.all.ec(self.ainvs())
            return self.__watkins_data

    def modular_degree(self, algorithm='sympow'):
        r"""
        Return the modular degree of this elliptic curve.

        The result is cached.  Subsequence calls, even with a
        different algorithm, just returned the cached result.

        INPUT:
           algorithm -- string:
              'sympow' -- (default) use Mark Watkin's (newer) C program sympow
              'ec' -- use Mark Watkins's C program ec
              'magma' -- requires that MAGMA be installed (also implemented
                         by Mark Watkins)

        \note{On 64-bit computers ec does not work, so \sage uses
        sympow even if ec is selected on a 64-bit computer.}

        The correctness of this function when called with algorithm "ec"
        is subject to the following three hypothesis:

        \begin{itemize}

            \item Manin's conjecture: the Manin constant is 1

            \item Steven's conjecture: the $X_1(N)$-optimal quotient
            is the curve with minimal Faltings height.  (This is proved
            in most cases.)

            \item The modular degree fits in a machine double, so it
            better be less than about 50-some bits.  (If you use sympow
            this contraint does not apply.)

        \end{itemize}

        Moreover for all algorithms, computing a certain value of an
        $L$-function ``uses a heuristic method that discerns when the
        real-number approximation to the modular degree is within epsilon
        [=0.01 for algorithm="sympow"] of the same integer for 3
        consecutive trials (which occur maybe every 25000 coefficients
        or so). Probably it could just round at some point. For
        rigour, you would need to bound the tail by assuming
        (essentially) that all the $a_n$ are as large as possible, but
        in practise they exhibit significant (square root)
        cancellation. One difficulty is that it doesn't do the sum in
        1-2-3-4 order; it uses 1-2-4-8---3-6-12-24--9-18-- (Euler
        product style) instead, and so you have to guess ahead of time
        at what point to curtail this expansion.''  (Quote from an
        email of Mark Watkins.)

        \note{If the curve is loaded from the large Cremona database,
        then the modular degree is taken from the database.}

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: E.modular_degree()
            1
            sage: E = EllipticCurve('5077a')
            sage: E.modular_degree()
            1984
            sage: factor(1984)
            2^6 * 31

            sage: EllipticCurve([0, 0, 1, -7, 6]).modular_degree()
            1984
            sage: EllipticCurve([0, 0, 1, -7, 6]).modular_degree(algorithm='sympow')
            1984
            sage: EllipticCurve([0, 0, 1, -7, 6]).modular_degree(algorithm='magma')  # optional
            1984

        We compute the modular degree of the curve with rank four having smallest
        (known) conductor:

            sage: E = EllipticCurve([1, -1, 0, -79, 289])
            sage: factor(E.conductor())
            2 * 117223
            sage: factor(E.modular_degree())
            2^7 * 2617
        """
        try:
            return self.__modular_degree

        except AttributeError:

            if misc.is_64_bit and algorithm == 'ec':
                misc.verbose('64-bit computer -- switching to algorithm sympow')
                algorithm = 'sympow'

            if algorithm == 'ec':
                v = self.watkins_data()
                m = rings.Integer(v["Modular degree"])
            elif algorithm == 'sympow':
                from sage.lfunctions.all import sympow
                m = sympow.modular_degree(self)
            elif algorithm == 'magma':
                m = rings.Integer(self._magma_().ModularDegree())
            else:
                raise ValueError, "unknown algorithm %s"%algorithm
            self.__modular_degree = m
            return m

    def modular_parametrization(self):
        """
        Computes and returns ...
        """
        return self.pari_mincurve().elltaniyama()

    def cremona_label(self, space=False):
        """
        Return the Cremona label associated to (the minimal model) of this curve,
        if it is known.  If not, raise a RuntimeError exception.
        """
        try:
            if not space:
                return self.__cremona_label.replace(' ','')
            return self.__cremona_label
        except AttributeError:
            try:
                X = self.database_curve()
            except RuntimeError:
                raise RuntimeError, "Cremona label not known for %s."%self
            self.__cremona_label = X.__cremona_label
            return self.cremona_label(space)

    def label(self):
        return self.cremona_label()

    def torsion_order(self):
        """
        Return the order of the torsion subgroup.
        """
        try:
            return self.__torsion_order
        except AttributeError:
            self.__torsion_order = self.torsion_subgroup().order()
            return self.__torsion_order

    def torsion_subgroup(self, flag=0):
        """
        Returns the torsion subgroup of this elliptic curve.

        INPUT:
            flag -- (default: 0)  chooses PARI algorithm:
              flag = 0: uses Doud algorithm
              flag = 1: uses Lutz-Nagell algorithm

        OUTPUT:
            The EllipticCurveTorsionSubgroup instance associated to this elliptic curve.

        EXAMPLES:
            sage: EllipticCurve('11a').torsion_subgroup()
            Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C5 associated to the Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: EllipticCurve('37b').torsion_subgroup()
            Torsion Subgroup isomorphic to Multiplicative Abelian Group isomorphic to C3 associated to the Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field
        """
        try:
            return self.__torsion_subgroup
        except AttributeError:
            self.__torsion_subgroup = rational_torsion.EllipticCurveTorsionSubgroup(self, flag)
            return self.__torsion_subgroup

    ## def newform_eval(self, z, prec):
##         """
##         The value of the newform attached to this elliptic curve at
##         the point z in the complex upper half plane, computed using
##         prec terms of the power series expansion.  Note that the power
##         series need not converge well near the real axis.
##         """
##         raise NotImplementedError

    def root_number(self):
        """
        Returns the root number of this elliptic curve.

        This is 1 if the order of vanishing of the L-function L(E,s)
        at 1 is even, and -1 if it is odd.
        """

        try:
            return self.__root_number
        except AttributeError:
            self.__root_number = int(self.pari_mincurve().ellrootno())
        return self.__root_number


    def has_cm(self):
        return self.j_invariant() in [0, 54000, -12288000, 1728, \
                                      287496, -3375, 16581375, 8000, \
                                      -32768,  -884736, -884736000,\
                                      -147197952000, -262537412640768000]

    def quadratic_twist(self, D):
        """
        Return the global minimal model of the quadratic twist of this curve by D.
        """
        return EllipticCurve_field.quadratic_twist(self, D).minimal_model()


    ##########################################################
    # Isogeny class (currently just uses Cremona database.)
    ##########################################################
    def isogeny_class(self, algorithm="mwrank", verbose=False):
        r"""
        Return all curves over $\Q$ in the isogeny class of this
        elliptic curve.

        INPUT:
            algorithm -- string:
                 "mwrank"   -- (default) use the mwrank C++ library
                 "database" -- use the Cremona database (only works if
                               curve is isomorphic to a curve in the database)

        OUTPUT:
            Returns the sorted list of the curves isogenous to self.
            If algorithm is "mwrank", also returns the isogeny matrix (otherwise
            returns None as second return value).

        \note{The result is \emph{not} provably correct, in the sense
            that when the numbers are huge isogenies could be missed
            because of precision issues.}

        \note{The ordering depends on which algorithm is used.}

        EXAMPLES:
            sage: I, A = EllipticCurve('37b').isogeny_class('mwrank')
            sage: I   # randomly ordered
            [Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 3*x +1 over Rational Field]
            sage: A
            [0 3 3]
            [3 0 0]
            [3 0 0]

            sage: I, _ = EllipticCurve('37b').isogeny_class('database'); I
            [Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field,
             Elliptic Curve defined by y^2 + y = x^3 + x^2 - 3*x +1 over Rational Field]

        This is an example of a curve with a $37$-isogeny:
            sage: E = EllipticCurve([1,1,1,-8,6])
            sage: E.isogeny_class ()
            ([Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 8*x + 6 over Rational Field,
              Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 208083*x - 36621194 over Rational Field],
             [ 0 37]
             [37  0])

        This curve had numerous $2$-isogenies:
        sage: e=EllipticCurve([1,0,0,-39,90])
            sage: e.isogeny_class ()
            ([Elliptic Curve defined by y^2 + x*y  = x^3 - 39*x + 90 over Rational Field,
              Elliptic Curve defined by y^2 + x*y  = x^3 - 4*x -1 over Rational Field,
              Elliptic Curve defined by y^2 + x*y  = x^3 + x over Rational Field,
              Elliptic Curve defined by y^2 + x*y  = x^3 - 49*x - 136 over Rational Field,
              Elliptic Curve defined by y^2 + x*y  = x^3 - 34*x - 217 over Rational Field,
              Elliptic Curve defined by y^2 + x*y  = x^3 - 784*x - 8515 over Rational Field],
             [0 2 0 0 0 0]
             [2 0 2 2 0 0]
             [0 2 0 0 0 0]
             [0 2 0 0 2 2]
             [0 0 0 2 0 0]
             [0 0 0 2 0 0])

        See \url{http://modular.ucsd.edu/Tables/nature/} for more interesting
        examples of isogeny structures.
        """
        #if algorithm == "gp":

       #     return sum([L for _, L in self.isogenous_curves(algorithm="gp")], [self])

        if algorithm == "mwrank":
            try:
                E = self.mwrank_curve()
            except ValueError:
                E = self.minimal_model().mwrank_curve()
            I, A = E.isogeny_class(verbose=verbose)
            mat = matrix.MatrixSpace(rings.IntegerRing(), len(A))(A)
            I = [constructor.EllipticCurve(ainvs) for ainvs in I]
            return I, mat

        elif algorithm == "database":

            try:
                label = self.cremona_label(space=False)
            except RuntimeError:
                raise RuntimeError, "unable to to find %s in the database"%self
            db = sage.databases.cremona.CremonaDatabase()
            I = db.isogeny_class(label)
            I.sort()
            return I, None

        else:

            raise ValueError, "unknown algorithm '%s'%"%algorithm


    ##########################################################
    # Galois Representations
    ##########################################################

    def is_reducible(self, p):
        """
        Return True if the mod-p representation attached
        to E is reducible.

        EXAMPLES:
            sage: E = EllipticCurve('121a'); E
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 30*x - 76 over Rational Field
            sage: E.is_reducible(7)
            False
            sage: E.is_reducible(11)
            True
            sage: EllipticCurve('11a').is_reducible(5)
            True
            sage: e = EllipticCurve('11a2')
            sage: e.is_reducible(5)
            True
            sage: e.torsion_order()
            1
        """
        # we do is_surjective first, since this is
        # much easier than computing isogeny_class
        t, why = self.is_surjective(p)
        if t == True:
            return False  # definitely not reducible
        isogeny_matrix = self.isogeny_class()[ 1 ]
        v = isogeny_matrix[0]  # first row
        for a in v:
            if a != 0 and a % p == 0:
                return True
        return False

    def is_irreducible(self, p):
        """
        Return True if the mod p represenation is irreducible.
        """
        return not self.is_reducible()

    def is_surjective(self, p, A=1000):
        """
        Return True if the mod-p representation attached to E
        is surjective, False if it is not, or None if we were
        unable to determine whether it is or not.

        INPUT:
            p -- int (a prime number)
            A -- int (a bound on the number of a_p to use)

        OUTPUT:
            a 2-tuple:
            -- surjective or (probably) not
            -- information about what it is if not surjective

        EXAMPLES:

        REMARKS:

            1.  If p >= 5 then the mod-p representation is surjective
                if and only if the p-adic representation is
                surjective.  When p = 2, 3 there are counterexamples.
                See a very recent paper of Elkies for more details
                when p=3.

            2.  When p <= 3 this function always gives the correct
                result irregardless of A, since it explicitly
                determines the p-division polynomial.

        """
        if not arith.is_prime(p):
            raise TypeError, "p (=%s) must be prime."%p

        T = self.torsion_subgroup().order()
        if T % p == 0:
            return False, "%s-torsion"%p

        if p == 2:
            invs = self.weierstrass_model().ainvs()
            R = rings.PolynomialRing(self.base_ring(), 'x')
            x = R.gen()
            f = x**3 + invs[3]*x + invs[4]
            if not f.is_irreducible():
                return False, '2-torsion'
            if arith.is_square(f.discriminant()):
                return False, "A3"
            return True, None

        if p == 3:
            # Algorithm: Let f be the 3-division polynomial, which is
            # a polynomial of degree 4.  Then I claim that this
            # polynomial has Galois group S_4 if and only if the
            # representation rhobar_{E,3} is surjective.  If the group
            # is S_4, then S_4 is a quotient of the image of
            # rhobar_{E,3}.  Since S_4 has order 24 and GL_2(F_3)
            # has order 48, the only possibility we have to consider
            # is that the image of rhobar is isomorphic to S_4.
            # But this is not the case because S_4 is not a subgroup
            # of GL_2(F_3).    If it were, it would be normal, since
            # it would have index 2.  But there is a *unique* normal
            # subgroup of GL_2(F_3) of index 2, namely SL_2(F_3),
            # and SL_2(F_3) is not isomorphic to S_4 (S_4 has a normal
            # subgroup of index 2 and SL_2(F_3) does not.)
            # (What's a simple way to see that SL_2(F_3) is the
            # unique index-2 normal subgroup?  I didn't see an obvious
            # reason, so just used the NormalSubgroups command in MAGMA
            # and it output exactly one of index 2.)

            # Here's Noam Elkies proof for the other direction:

            #> Let E be an elliptic curve over Q.  Is the mod-3
            #> representation E[3]  surjective if and only if the
            #> (degree 4) division polynomial has Galois group S_4?  I
            #> can see why the group being S_4 implies the
            #> representation is surjective, but the converse is not
            #> clear to me.
            # I would have thought that this is the easier part: to
            # say that E[3] is surjective is to say the 3-torsion
            # field Q(E[3]) has Galois group GL_2(Z/3) over Q.  Let
            # E[3]+ be the subfield fixed by the element -1 of
            # GL_2(Z/3).  Then E[3] has Galois group PGL_2(Z/3), which
            # is identified with S_4 by its action on the four
            # 3-element subgroups of E[3].  Each such subgroup is in
            # turn determined by the x-coordinate shared by its two
            # nonzero points.  So, if E[3] is surjective then any
            # permutation of those x-coordinates is realized by some
            # element of Gal(E[3]+/Q).  Thus the Galois group of the
            # division polynomial (whose roots are those
            # x-coordinates) maps surjectively to S_4, which means it
            # equals S_4.


            f = self.division_polynomial(3)
            if not f.is_irreducible():
                return False, "reducible_3-divpoly"
            n = pari(f).polgalois()[0]
            if n == 24:
                return True, None
            else:
                return False, "3-divpoly_galgroup_order_%s"%n

        if self.has_cm():
            return False, "CM"
        an = self.anlist(A)
        ell = 0
        Np = self.conductor() * p
        signs = []
        while True:
            ell = arith.next_prime(ell)
            if ell >= A: break
            if Np % ell != 0:
                a_ell = an[int(ell)]
                if a_ell % p != 0:
                    s = arith.kronecker(a_ell**2 - 4*ell, p)
                    #print ell, s
                    if s == 0: continue
                    if not (s in signs):
                        signs.append(s)
                        if len(signs) == 2:
                            return True, None

        # could do something further here...
        return False, signs

    def is_semistable(self):
        if self.base_ring() != Q:
            raise NotImplementedError, "is_semistable only implemented for curves over the rational numbers."
        return self.conductor().is_squarefree()

    def reducible_primes(self):
        r"""
        Returns a list of the primes $p$ such that the mod $p$
        representation $\rho_{E,p}$ is reducible.  For all other
        primes the representation is irreducible.

        NOTE -- this is \emph{not} provably correct in general.
        See the documentation for \code{self.isogeny_class}.

        EXAMPLES:
            sage: E = EllipticCurve('225a')
            sage: E.reducible_primes()
            [3]
        """
        try:
            return self.__reducible_primes
        except AttributeError:
            pass
        C, I = self.isogeny_class(algorithm='mwrank')
        X = set(I.list())
        R = [p for p in X if arith.is_prime(p)]
        self.__reducible_primes = R
        return R

    def non_surjective(self, A=1000):
        r"""
        Returns a list of primes p such that the mod-p representation
        $\rho_{E,p}$ *might* not be surjective (this list usually
        contains 2, because of shortcomings of the algorithm).  If p
        is not in the returned list, then rho_{E,p} is provably
        surjective (see A. Cojocaru's paper).  If the curve has CM
        then infinitely many representations are not surjective, so we
        simply return the sequence [(0,"cm")] and do no further computation.

        INPUT:
            A -- an integer
        OUTPUT:
            list -- if curve has CM, returns [(0,"cm")].  Otherwise, returns a
                    list of primes where mod-p representation very likely
                    not surjective.   At any prime not in this list,
                    the representation is definitely surjective.
        EXAMPLES:
            sage: E = EllipticCurve([0, 0, 1, -38, 90])  # 361A
            sage: E.non_surjective()   # CM curve
            [(0, 'cm')]

            sage: E = EllipticCurve([0, -1, 1, 0, 0]) # X_1(11)
            sage: E.non_surjective()
            [(5, '5-torsion')]

            sage: E = EllipticCurve([0, 0, 1, -1, 0]) # 37A
            sage: E.non_surjective()
            []

            sage: E = EllipticCurve([0,-1,1,-2,-1])   # 141C
            sage: E.non_surjective()
            [(13, [1])]

        ALGORITHM:
            When p<=3 use division polynomials.  For 5 <= p <= B,
            where B is Cojocaru's bound, use the results in Section 2
            of Serre's inventiones paper"Sur Les Representations Modulaires Deg
            Degre 2 de Galqbar Over Q."
        """
        if self.has_cm():
            misc.verbose("cm curve")
            return [(0,"cm")]
        N = self.conductor()
        if self.is_semistable():
            C = 11
            misc.verbose("semistable -- so bound is 11")
        else:
            C = 1 + 4*sqrt(6)*int(N)/3 * sqrt(mul([1+1.0/int(p) for p,_ in factor(N)]))
            misc.verbose("conductor = %s, and bound is %s"%(N,C))
        C = 1 + 4*sqrt(6)*int(N)/3 * sqrt(mul([1+1.0/int(p) for p,_ in factor(N)]))
        misc.verbose("conductor = %s, and bound is %s"%(N,C))
        B = []
        p = 2
        while p <= C:
            t, v = self.is_surjective(p, A=A)
            misc.verbose("(%s,%s,%s)"%(p,t,v))
            if not t:
                B.append((p,v))
            p = next_prime(p)
        return B

    def is_ordinary(self, p, ell=None):
        """
        Return True precisely when the mod-p representation attached
        to this elliptic curve is ordinary at ell.

        INPUT:
            p -- a prime
            ell - a prime (default: p)

        OUTPUT:
            bool
        """
        if ell is None:
            ell = p
        return self.ap(ell) % p != 0

    def is_good(self, p, check=True):
        """
        Return True if $p$ is a prime of good reduction for $E$.

        INPUT:
            p -- a prime

        OUTPUT:
            bool

        EXAMPLES:
            sage: e = EllipticCurve('11a')
            sage: e.is_good(-8)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime
            sage: e.is_good(-8, check=False)
            True
        """
        if check:
            if not arith.is_prime(p):
                raise ValueError, "p must be prime"
        return self.conductor() % p != 0


    def is_supersingular(self, p, ell=None):
        """
        Return True precisely when the mod-p representation attached
        to this elliptic curve is supersingular at ell.

        INPUT:
            p -- a prime
            ell - a prime (default: p)

        OUTPUT:
            bool
        """
        if ell is None:
            ell = p
        return not self.is_ordinary(p, ell)

    def eval_modular_form(self, points, prec):
        if not isinstance(points, (list,xrange)):
            try:
                points = list(points)
            except TypeError:
                return self.eval_modular_form([points],prec)
        an = self.pari_mincurve().ellan(prec)
        s = 0
        I = pari("I")
        pi = pari("Pi")
        c = pari(2)*pi*I
        ans = []
        for z in points:
            s = pari(0)
            r0 = (c*z).exp()
            r = r0
            for n in xrange(1,prec):
                s += an[n-1]*r
                r *= r0
            ans.append(s.python())
        return ans

    ########################################################################
    # Functions related to the BSD conjecture.
    ########################################################################
    def sha_an(self, use_database=False):
        """
        Returns the Birch and Swinnerton-Dyer conjectural order of
        Sha, unless the analytic rank is > 1, in which case this
        function returns 0.

        This result is proved correct if the order of vanishing is 0
        and the Manin constant is <= 2.

        If the optional parameter use_database is True (default:
        False), this function returns the analytic order of Sha as
        listed in Cremona's tables, if this curve appears in Cremona's
        tables.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.sha_an()
            1
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.sha_an()
            1

        The smallest conductor curve with nontrivial Sha:
            sage: E = EllipticCurve([1,1,1,-352,-2689])     # 66B3
            sage: E.sha_an()
            4

        The four optimal quotients with nontrivial Sha and conductor <= 1000:
            sage: E = EllipticCurve([0, -1, 1, -929, -10595])       # 571A
            sage: E.sha_an()
            4
            sage: E = EllipticCurve([1, 1, 0, -1154, -15345])       # 681B
            sage: E.sha_an()
            9
            sage: E = EllipticCurve([0, -1, 0, -900, -10098])       # 960D
            sage: E.sha_an()
            4
            sage: E = EllipticCurve([0, 1, 0, -20, -42])            # 960N
            sage: E.sha_an()
            4

        The smallest conductor curve of rank > 1:
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.sha_an()
            0

        The following are examples that require computation of the Mordell-Weil
        group and regulator:

            sage: E = EllipticCurve([0, 0, 1, -1, 0])                     # 37A  (rank 1)
            sage: E.sha_an()
            1

            sage: E = EllipticCurve("1610f3")
            sage: E.sha_an()
            4

        In this case the input curve is not minimal, and if this function didn't
        transform it to be minimal, it would give nonsense:
            sage: E = EllipticCurve([0,-432*6^2])
            sage: E.sha_an()
            1
        """
#            sage: e = EllipticCurve([1, 0, 0, -19491080, -33122512122])   # 15834T2
#            sage: e.sha_an()                          # takes a long time (way too long!!)
#            25
        if hasattr(self, '__sha_an'):
            return self.__sha_an
        if use_database:
            try:
                self.__sha_an = int(round(float(self.database_curve().db_extra[4])))
                return self.__sha_an
            except RuntimeError, AttributeError:
                pass

        # it's critical to switch to the minimal model.
        E = self.minimal_model()
        eps = E.root_number()
        if eps == 1:
            L1_over_omega = E.L_ratio()
            if L1_over_omega == 0:
                return 0
            T = E.torsion_subgroup().order()
            Sha = (L1_over_omega * T * T) / Q(E.tamagawa_product())
            try:
                Sha = Z(Sha)
            except ValueError:
                raise RuntimeError, \
                      "There is a bug in sha_an, since the computed conjectural order of Sha is %s, which is not an integer."%Sha
            if not arith.is_square(Sha):
                raise RuntimeError, \
                      "There is a bug in sha_an, since the computed conjectural order of Sha is %s, which is not a square."%Sha
            E.__sha_an = Sha
            self.__sha_an = Sha
            return Sha

        else:  # rank > 0  (Not provably correct)
            L1, error_bound = E.Lseries_deriv_at1(10*sqrt(E.conductor()) + 10)
            if abs(L1) < error_bound:
                E.__sha_an = 0
                self.__sha_an = 0
                return 0   # vanishes to order > 1, to computed precision
            regulator = E.regulator()   # this could take a *long* time; and could fail...?
            T = E.torsion_subgroup().order()
            omega = E.omega()
            Sha = int(round ( (L1 * T * T) / (E.tamagawa_product() * regulator * omega) ))
            try:
                Sha = Z(Sha)
            except ValueError:
                raise RuntimeError, \
                      "There is a bug in sha_an, since the computed conjectural order of Sha is %s, which is not an integer."%Sha
            if not arith.is_square(Sha):
                raise RuntimeError, \
                      "There is a bug in sha_an, since the computed conjectural order of Sha is %s, which is not a square."%Sha
            E.__sha_an = Sha
            self.__sha_an = Sha
            return Sha

    def sha_an_padic(self, p):
        """
        Return the power of p that divides Sha(E)(p), according to the
        p-adic analogue of the BSD conjecture.

        INPUT:
            p -- a prime

        OUTPUT:
            integer -- power of p that conjecturally equals #Sha(E)(p)

        Note that in many cases this conjecture has been proved.
        """
        try:
            return self.__sha_an_padic[p]
        except AttributeError:
            self.__sha_an_padic = {}
        except KeyError:
            pass

        if self.is_ordinary(p) and self.is_good(p):
            S = self._sha_an_padic_good_ordinary(p)
        else:
            raise NotImplementedError, "only the good ordinary case is implemented."
        self.__sha_an_padic[p] = S
        return S

    def _sha_an_padic_good_ordinary(self, p):
        """
        Return the power of p that divides Sha(E)(p), according to the
        p-adic analogue of the BSD conjecture.

        INPUT:
            p -- a prime of good ordinary reduction for E

        OUTPUT:
            integer -- power of p that conjecturally equals #Sha(E)(p)
        """



    def L_ratio(self):
        r"""
        Returns the ratio $L(E,1)/\Omega$ as an exact rational
        number. The result is \emph{provably} correct if the Manin
        constant of the associated optimal quotient is $\leq 2$.  This
        hypothesis on the Manin constant is true for all semistable
        curves (i.e., squarefree conductor), by a theorem of Mazur
        from his \emph{Rational Isogenies of Prime Degree} paper.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.L_ratio()
            1/5
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.L_ratio()
            1/25
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.L_ratio()
            0
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.L_ratio()
            0
            sage: E = EllipticCurve([0, 0, 1, -38, 90])     # 361A (CM curve))
            sage: E.L_ratio()
            0
            sage: E = EllipticCurve([0,-1,1,-2,-1])         # 141C (13-isogeny)
            sage: E.L_ratio()
            1
            sage: E = EllipticCurve(RationalField(), [1, 0, 0, 1/24624, 1/886464])
            sage: E.L_ratio()
            2

        WARNING: It's conceivable that machine floats are not large
        enough precision for the computation; if this could be the
        case a RuntimeError is raised.  The curve's real period would
        have to be very small for this to occur.

        ALGORITHM: Compute the root number.  If it is -1 then L(E,s)
        vanishes to odd order at 1, hence vanishes.  If it is +1, use
        a result about modular symbols and Mazur's "Rational Isogenies"
        paper to determine a provably correct bound (assuming Manin
        constant is <= 2) so that we can determine whether L(E,1) = 0.

        AUTHOR: William Stein, 2005-04-20.
        """
        try:
            return self.__lratio
        except AttributeError:
            pass

        if not self.is_minimal():
            self.__lratio = self.minimal_model().L_ratio()
            return self.__lratio

        if self.root_number() == -1:
            self.__lratio = Q(0)
            return self.__lratio

        # Even root number.  Decide if L(E,1) = 0.  If E is a modular
        # *optimal* quotient of J_0(N) elliptic curve, we know that T *
        # L(E,1)/omega is an integer n, where T is the order of the
        # image of the rational torsion point (0)-(oo) in E(Q), and
        # omega is the least real Neron period.  (This is proved in my
        # Ph.D. thesis, but is probably well known.)  We can easily
        # compute omega to very high precision using AGM.  So to prove
        # that L(E,1) = 0 we compute T/omega * L(E,1) to sufficient
        # precision to determine it as an integer.  If eps is the
        # error in computation of L(E,1), then the error in computing
        # the product is (2T/Omega_E) * eps, and we need this to be
        # less than 0.5, i.e.,
        #          (2T/Omega_E) * eps < 0.5,
        # so
        #          eps < 0.5 * Omega_E / (2T) = Omega_E / (4*T).
        #
        # Since in general E need not be optimal, we have to choose
        # eps = Omega_E/(8*t*B), where t is the exponent of E(Q)_tor,
        # and B is a bound on the degree of any isogeny.   A liberal
        # bound on the degrees of cyclic N-isogenies is 200, by Mazur's
        # "Rational Isogenies of Prime Degree" paper, so we take B=200.
        #
        # NOTES: We *do* have to worry about the Manin constant, since
        # we are using the Neron model to compute omega, not the
        # newform.  My theorem replaces the omega above by omega/c,
        # where c is the Manin constant, and the bound must be
        # correspondingly smaller.  If the level is square free, then
        # the Manin constant is 1 or 2, so there's no problem (since
        # we took 8 instead of 4 in the denominator).  If the level
        # is divisible by a square, then the Manin constant could
        # be a divisible by an arbitrary power of that prime, except
        # that Edixhoven claims the primes that appear are <= 7.

        t = self.torsion_subgroup().exponent()
        omega = self.period_lattice()[0]
        C = 8*200*t
        eps = omega / C
        #   coercion of 10**(-15) to our real field is needed to make unambiguous comparison
        if eps < R(10**(-15)):  # liberal bound on precision of float
            raise RuntimeError, "Insufficient machine precision (=%s) for computation."%eps
        sqrtN = 2*int(sqrt(self.conductor()))
        k = sqrtN + 10
        while True:
            L1, error_bound = self.Lseries_at1(k)
            if error_bound < eps:
                n = int(round(L1*C/omega))
                quo = Q(n) / Q(C)
                self.__lratio = quo / self.real_components()
                return self.__lratio
            k += sqrtN
            misc.verbose("Increasing precision to %s terms."%k)


    def L1_vanishes(self):
        """
        Returns whether or not L(E,1) = 0. The result is provably
        correct if the Manin constant of the associated optimal
        quotient is <= 2.  This hypothesis on the Manin constant
        is true for all curves of conductor <= 40000 (by Cremona) and
        all semistable curves (i.e., squarefree conductor).

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.L1_vanishes()
            False
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.L1_vanishes()
            False
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.L1_vanishes()
            True
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.L1_vanishes()
            True
            sage: E = EllipticCurve([0, 0, 1, -38, 90])     # 361A (CM curve))
            sage: E.L1_vanishes()
            True
            sage: E = EllipticCurve([0,-1,1,-2,-1])         # 141C (13-isogeny)
            sage: E.L1_vanishes()
            False

        WARNING: It's conceivable that machine floats are not large
        enough precision for the computation; if this could be the
        case a RuntimeError is raised.  The curve's real period would
        have to be very small for this to occur.

        ALGORITHM: Compute the root number.  If it is -1 then L(E,s)
        vanishes to odd order at 1, hence vanishes.  If it is +1, use
        a result about modular symbols and Mazur's "Rational Isogenies"
        paper to determine a provably correct bound (assuming Manin
        constant is <= 2) so that we can determine whether L(E,1) = 0.

        AUTHOR: William Stein, 2005-04-20.
        """
        return self.L_ratio() == 0

    ########################################################################
    # Functions related to bounding the order of Sha (provably correctly!)
    # Heegner points and Kolyvagin's theorem
    ########################################################################
    def two_selmer_shabound(self):
        """
        Returns a bound on the dimension of Sha(E)[2], computed using
        a 2-descent.
        """
        S = self.selmer_rank_bound()
        r = self.rank()
        t = self.two_torsion_rank()
        b = S - r - t
        if b % 2 != 0:
            raise ArithmeticError, "There is a bug in two_selmer_shabound since it's %s, but it must be even."%b
        return b

    def satisfies_heegner_hypothesis(self, D):
        """
        Returns True precisely when D is a fundamental discriminant
        that satisfies the Heegner hypothesis for this elliptic curve.
        """
        if not number_field.is_fundamental_discriminant(D):
            return False
        if arith.GCD(D, self.conductor()) != 1:
            return False
        K = number_field.QuadraticField(D, 'a')
        for p, _ in factor(self.conductor()):
            if len(K.factor_integer(p)) != 2:
                return False
        return True

    def heegner_discriminants(self, bound):
        return [-D for D in xrange(1,bound) if self.satisfies_heegner_hypothesis(-D)]

    def heegner_discriminants_list(self, n):
        """
        List of the first n Heegner discriminants for self.
        """
        v = []
        D = -5
        while len(v) < n:
            while not self.satisfies_heegner_hypothesis(D):
                D -= 1
            v.append(D)
            D -= 1
        return v

    def heegner_point_height(self, D, prec=2):
        """
        Use the Gross-Zagier formula to compute the Neron-Tate
        canonical height over K of the Heegner point corresponding to
        D, as an Interval (since it's computed to some precision using
        L-functions).

        INPUT:
            D (int) -- fundamental discriminant (=/= -3, -4)
            prec (int) -- (default: 2), use prec*sqrt(N) + 20 terms
                          of L-series in computations, where N is the
                          conductor.

        OUTPUT:
            Interval that contains the height of the Heegner point.

        EXAMPLE:
            sage: E = EllipticCurve('11a')
            sage: E.heegner_point_height(-7)
            [0.22226977 ... 0.22227479]
        """

        if not self.satisfies_heegner_hypothesis(D):
            raise ArithmeticError, "Discriminant (=%s) must be a fundamental discriminant that satisfies the Heegner hypothesis."%D
        if D == -3 or D == -4:
            raise ArithmeticError, "Discriminant (=%s) must not be -3 or -4."%D
        eps = self.root_number()
        L1_vanishes = self.L1_vanishes()
        if eps == 1 and L1_vanishes:
            return IR(0) # rank even hence >= 2, so Heegner point is torsion.
        alpha = R(sqrt(abs(D)))/(2*self.complex_area())
        F = self.quadratic_twist(D)
        E = self
        k_E = prec*sqrt(E.conductor()) + 20
        k_F = prec*sqrt(F.conductor()) + 20

        IR = rings.RealIntervalField(20)
        MIN_ERR = R('1e-6')   # we assume that regulator and
                            # discriminant, etc., computed to this accuracy.
                            # this should be made more intelligent / rigorous relative
                             # to the rest of the system.
        if eps == 1:   # E has even rank
            LF1, err_F = F.Lseries_deriv_at1(k_F)
            LE1, err_E = E.Lseries_at1(k_E)
            err_F = max(err_F, MIN_ERR)
            err_E = max(err_E, MIN_ERR)
            return IR(alpha-MIN_ERR,alpha+MIN_ERR) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F)

        else:          # E has odd rank
            LE1, err_E = E.Lseries_deriv_at1(k_E)
            LF1, err_F = F.Lseries_at1(k_F)
            err_F = max(err_F, MIN_ERR)
            err_E = max(err_E, MIN_ERR)
            return IR(alpha-MIN_ERR,alpha+MIN_ERR) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F)


    def heegner_index(self, D,  min_p=3, prec=5, verbose=False):
        """
        Return an interval that contains the SQUARE of the index of
        the Heegner point in the group of K-rational points *modulo
        torsion* on the twist of the elliptic curve by D, computed
        using the Gross-Zagier formula and/or a point search.

        WARNING: This function uses the Gross-Zagier formula.
        When E is 681b and D=-8 for some reason the returned index
        is 9/4 which is off by a factor of 4.   Evidently the
        GZ formula must be modified when D=-8.

        If 0 is in the interval of the height of the Heegner point
        computed to the given prec, then this function returns 0.

        INPUT:
            D (int) -- Heegner discriminant
            min_p (int) -- (default: 3) only rule out primes >= min_p
                           dividing the index.
            verbose (bool) -- (default: False); print lots of mwrank search status
                                                information when computing regulator
            prec (int) -- (default: 5), use prec*sqrt(N) + 20 terms
                          of L-series in computations, where N is the
                          conductor.

        OUTPUT:
            an interval that contains the index

        EXAMPLES:
            sage: E = EllipticCurve('11a')
            sage: E.heegner_discriminants(50)
            [-7, -8, -19, -24, -35, -39, -40, -43]
            sage: E.heegner_index(-7)
            [0.99998760 ... 1.0000134]

            sage: E = EllipticCurve('37b')
            sage: E.heegner_discriminants(100)
            [-3, -4, -7, -11, -40, -47, -67, -71, -83, -84, -95]
            sage: E.heegner_index(-95)          # long time (1 second)
            [3.9999771 ... 4.0000229]

        Current discriminants -3 and -4 are not supported:
            sage: E.heegner_index(-3)
            Traceback (most recent call last):
            ...
            ArithmeticError: Discriminant (=-3) must not be -3 or -4.
        """
        # First compute upper bound on height of Heegner point.
        tm = misc.verbose("computing heegner point height...")
        h0 = self.heegner_point_height(D, prec=prec)

        # We divide by 2 to get the height **over Q** of the
        # Heegner point on the twist.

        ht = h0/2
        misc.verbose('Height of heegner point = %s'%ht, tm)

        if self.root_number() == 1:
            F = self.quadratic_twist(D)
        else:
            F = self
        h  = ht.upper()
        misc.verbose("Heegner height bound = %s"%h)
        B = F.CPS_height_bound()
        misc.verbose("CPS bound = %s"%B)
        c = h/(min_p**2) + B
        misc.verbose("Search would have to be up to height = %s"%c)

        if c > _MAX_HEIGHT or F is self:
            misc.verbose("Doing direct computation of MW group.")
            reg = F.regulator(verbose=verbose)
            return ht/IR(reg)

        # Do naive search to eliminate possibility that Heegner point
        # is divisible by p<min_p, without finding Heegner point.
        misc.verbose("doing point search")
        P = F.point_search(c, verbose=verbose)
        misc.verbose("done with point search")
        P = [x for x in P if x.order() == oo]
        if len(P) == 0:
            return IR(1)
        misc.verbose("saturating")
        S, I, reg = F.saturation(P, verbose=verbose)
        misc.verbose("done saturating")
        return ht/IR(reg)


    def heegner_index_bound(self, D=0,  prec=5, verbose=True, max_height=_MAX_HEIGHT):
        """
        Assume self has rank 0.

        Return a list v of primes such that if an odd prime p divides
        the index of the the Heegner point in the group of rational
        points *modulo torsion*, then p is in v.

        If 0 is in the interval of the height of the Heegner point
        computed to the given prec, then this function returns v = 0.
        This does not mean that the Heegner point is torsion, just
        that it is very likely torsion.

        If we obtain no information from a search up to max_height, e.g.,
        if the Siksek et al. bound is bigger than max_height, then
        we return v = -1.

        INPUT:
            D (int) -- (deault: 0) Heegner discriminant; if 0, use the
                       first discriminant < -4 that satisfies the Heegner hypothesis
            verbose (bool) -- (default: True)
            prec (int) -- (default: 5), use prec*sqrt(N) + 20 terms
                          of L-series in computations, where N is the
                          conductor.
            max_height (float) -- should be <= 21; bound on logarithmic naive height
                                  used in point searches.  Make smaller to make this
                                  function faster, at the expense of possibly obtaining
                                  a worse answer.  A good range is between 13 and 21.

        OUTPUT:
            v -- list or int (bad primes or 0 or -1)
            D -- the discriminant that was used (this is useful if D was
                 automatically selected).
        """
        max_height = min(float(max_height), _MAX_HEIGHT)
        if self.root_number() != 1:
            raise RuntimeError, "The rank must be 0."

        if D == 0:
            D = -5
            while not self.satisfies_heegner_hypothesis(D):
                D -= 1

        # First compute upper bound on Height of Heegner point.
        ht = self.heegner_point_height(D, prec=prec)
        if 0 in ht:
            return 0, D
        F = self.quadratic_twist(D)
        h  = ht.upper()
        misc.verbose("Heegner height bound = %s"%h)
        B = F.CPS_height_bound()
        misc.verbose("CPS bound = %s"%B)
        H = h
        p = 3
        while True:
            c = h/(2*p**2) + B
            if c < max_height:
                break
            if p > 100:
                break
            p = next_prime(p)
        misc.verbose("Using p = %s"%p)

        if c > max_height:
            misc.verbose("No information by searching only up to max_height (=%s)."%c)
            return -1, D

        misc.verbose("Searching up to height = %s"%c)
        eps = 10e-5

        def _bound(P):
            """
            We will use this function below in two places.  It
            bounds the index using a nontrivial point.
            """
            assert len(P) == 1

            S, I, reg = F.saturation(P, verbose=verbose)
            h = IR(reg-eps,reg+eps)
            ind2 = ht/(h/2)
            misc.verbose("index squared = %s"%ind2)
            ind = ind2.sqrt()
            misc.verbose("index = %s"%ind)
            # Compute upper bound on square root of index.
            if ind.absolute_diameter() < 1:
                t, i = ind.is_int()
                if t:   # unique integer in interval, so we've found exact index squared.
                    return arith.prime_divisors(i), D
            raise RuntimeError, "Unable to compute bound for e=%s, D=%s (try increasing precision)"%(self,D)

        # First try a quick search, in case we get lucky and find
        # a generator.
        P = F.point_search(13, verbose=verbose)
        P = [x for x in P if x.order() == oo]
        if len(P) > 0:
            return _bound(P)

        # Do search to eliminate possibility that Heegner point is
        # divisible by primes up to p, without finding Heegner point.
        P = F.point_search(c, verbose=verbose)
        P = [x for x in P if x.order() == oo]
        if len(P) == 0:
            # We've eliminated the possibility of a divisor up to p.
            return arith.prime_range(3,p), D
        else:
            return _bound(P)


    def shabound_kolyvagin(self, D=0, regulator=None,
                           ignore_nonsurj_hypothesis=False):
        """
        Given a fundamental discriminant D (!= -3,-4) that satisfies the
        Heegner hypothesis, return a list of primes so that
        Kolyvagin's theorem (as in Gross's paper) implies that any
        prime divisor of $\#Sha$ is in this list.

        INPUT:
            D -- (optional) a fundamental discriminant < -4 that satisfies the
                 Heegner hypothesis for E; if not given, use the first such D

            regulator -- (optional) regulator of E(K); if not given, will
                         be computed (which could take a long time)


            ignore_nonsurj_hypothesis (optional: default False) --
                      If True, then gives the bound coming from Heegner point
                      index, but without any hypothesis on surjectivity
                      of the mod-p representation.


        OUTPUT:
            bound and index

        More precisely:

                0 -- if E/K has complex multiplication or analytic rank >= 2
            or
                B -- list of primes such that if p divides Sha(E/K), then p
                     is in B.

            and

                I -- the odd part of the index of the Heegner point in the full
                     group of K-rational points on E.  (If E has CM, returns 0.)

        REMARKS:
            (1) We do not have to assume that the Manin constant is 1
                (or a power of 2).  If the Manin constant were
                divisible by a prime, that prime would get included in
                the list of bad primes.

            (2) We assume the Gross-Zagier theorem is True under the
                hypothesis that gcd(N,D) = 1, instead of the stronger
                hypothesis gcd(2*N,D)=1 that is in the original
                Gross-Zagier paper.  That Gross-Zagier is true when
                gcd(N,D)=1 is"well-known" to the experts, but doesn't
                seem to written up well in the literature.

            (3) Correctness of the computation is guaranteed using
                interval arithmetic, under the assumption that the
                regulator, square root, and period lattice are
                computed to precision at least $10^{-10}$, i.e., they are
                correct up to addition or a real number with absolute
                value less than $10^{-10}$.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: E.shabound_kolyvagin()
            ([2], 1)
            sage: E = EllipticCurve('141a')
            sage: E.sha_an()
            1
            sage: E.shabound_kolyvagin()
            ([2, 7], 49)

        We get no information the curve has rank $2$.
            sage: E = EllipticCurve('389a')
            sage: E.shabound_kolyvagin()
            (0, 0)
            sage: E = EllipticCurve('681b')
            sage: E.sha_an()
            9
            sage: E.shabound_kolyvagin()
            ([2, 3], 9)

        """
        if self.has_cm():
            return 0, 0

        if D == 0:
            D = -5
            while not self.satisfies_heegner_hypothesis(D):
                D -= 1

        if not self.satisfies_heegner_hypothesis(D):
            raise ArithmeticError, "Discriminant (=%s) must be a fundamental discriminant that satisfies the Heegner hypothesis."%D
        if D == -3 or D == -4:
            raise ArithmeticError, "Discriminant (=%s) must not be -3 or -4."%D
        eps = self.root_number()
        L1_vanishes = self.L1_vanishes()
        if eps == 1 and L1_vanishes:
            return 0, 0        # rank even hence >= 2, so Kolyvagin gives nothing.
        alpha = sqrt(abs(D))/(2*self.complex_area())
        F = self.quadratic_twist(D)
        E = self
        k_E = 2*sqrt(E.conductor()) + 10
        k_F = 2*sqrt(F.conductor()) + 10
        #k_E = 2
        #k_F = 2

        MIN_ERR = 1e-10   # we assume that regulator and
                          # discriminant, etc., computed to this accuracy.
        tries = 0
        while True:
            tries += 1
            if tries >= 6:
                raise RuntimeError, "Too many precision increases in shabound_kolyvagin"
            if eps == 1:   # E has even rank
                misc.verbose("Conductor of twist = %s"%F.conductor())
                LF1, err_F = F.Lseries_deriv_at1(k_F)
                LE1, err_E = E.Lseries_at1(k_E)
                err_F = max(err_F, MIN_ERR)
                err_E = max(err_E, MIN_ERR)
                if regulator != None:
                    hZ = regulator/2
                else:
                    hZ = F.regulator(use_database=True)/2
                #print  alpha * LE1 * LF1 / hZ
                I = IR(alpha) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F) / hZ
                #print I

            else:          # E has odd rank

                if regulator != None:
                    hZ = regulator/2
                else:
                    hZ = self.regulator(use_database=True)/2
                LE1, err_E = E.Lseries_deriv_at1(k_E)
                LF1, err_F = F.Lseries_at1(k_F)
                err_F = max(err_F, MIN_ERR)
                err_E = max(err_E, MIN_ERR)
                #I = alpha * LE1 * LF1 / hZ

                I = IR(alpha) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F) / hZ

            misc.verbose('interval = %s'%I)
            t, n = I.is_int()
            if t:
                break
            elif I.absolute_diameter() < 1:
                raise RuntimeError, "Problem in shabound_kolyvagin; square of index is not an integer -- D=%s, I=%s."%(D,I)
            misc.verbose("Doubling bounds")
            k_E *= 2
            k_F *= 2
        # end while

        # We include 2 since Kolyvagin (in Gross) says nothing there
        if n == 0:  return 0, 0  # no bound
        F = factor(n)
        B = [2]
        for p, e in factor(n):
            if p > 2:
                if e%2 != 0:
                    raise RuntimeError, "Problem in shabound_kolyvagin; square of index is not a perfect square!  D=%s, I=%s, n=%s, e=%s."%(D,I,n,e)
                B.append(p)
            else:
                n /= 2**e  # replace n by its odd part
        if not ignore_nonsurj_hypothesis:
            for p, _ in self.non_surjective():
                B.append(p)
        B = list(set([int(x) for x in B]))
        B.sort()
        return B, n


    def shabound_kato(self):
        """
        Returns a list p of primes such that the theorems of Kato's
        and others (e.g., as explained in a paper/thesis of Grigor
        Grigorov) imply that if p divides $\\#Sha(E)$ then $p$ is in
        the list.

        If L(E,1) = 0, then Kato's theorem gives no information, so
        this function returns False.

        THEOREM (Kato): Suppose p >= 5 is a prime so the p-adic
        representation rho_{E,p} is surjective.  Then $ord_p(\\#Sha(E))$
        divides $ord_p(L(E,1)/Omega_E)$.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.shabound_kato()
            [2, 3, 5]
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.shabound_kato()
            [2, 3, 5]
            sage: E = EllipticCurve([1,1,1,-352,-2689])     # 66B3
            sage: E.shabound_kato()
            [2, 3]

        For the following curve one really has 25 | $\\#Sha$ (by Grigorov-Stein paper):
            sage: E = EllipticCurve([1, -1, 0, -332311, -73733731])   # 1058D1
            sage: E.shabound_kato()                 # long time (about 1 second)
            [2, 3, 5]
            sage: E.non_surjective()                # long time (about 1 second)
            []

        For this one, Sha is divisible by 7.
            sage: E = EllipticCurve([0, 0, 0, -4062871, -3152083138])   # 3364C1
            sage: E.shabound_kato()                 # long time (< 10 seconds)
            [2, 3, 7]

        No information about curves of rank > 0:
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.shabound_kato()
            False
        """
        if self.has_cm():
            return False
        if self.L1_vanishes():
            return False
        B = [2,3]
        for p, _ in self.non_surjective():   # for p >= 5, mod-p surj => p-adic surj
            if p > 3:
                B.append(p)

        # The only other p that might divide B are those that divide
        # the integer 2*#E(Q)_tor^2 * L(E,1)/omega.  So we compute
        # that to sufficient precision to determine it.  Note that
        # we have to assume the Manin constant is <=2 in order to provably
        # compute L(E,1)/omega.
        for p, n in factor(self.sha_an()):
            if n >= 2:    # use parity of Sha
                B.append(int(p))
        B = list(set(B))
        B.sort()
        return B

    def shabound(self):
        """
        Compute a provably correct bound on the order of the Shafarevich-Tate
        group of this curve. The bound is a either False (no bound) or a list
        B of primes such that any divisor of Sha is in this list.
        """
        if self.L1_vanishes():
            B = self.shabound_kolyvagin()
        else:
            B = self.shabound_kato()
        return B

    def __check_padic_hypotheses(self, p):
        p = rings.Integer(p)
        if not p.is_prime():
            raise ValueError, "p = (%s) must be prime"%p
        if p == 2:
            raise ValueError, "p must be odd"
        if self.conductor() % p == 0 or self.ap(p) % p == 0:
            raise ArithmeticError, "p must be a good ordinary prime"
        return p


    # This is the old version of padic_height that requires MAGMA:
    def padic_height_magma(self, p, P, prec=20):
        """
        Return the cyclotomic $p$-adic height of $P$, in the sense
        of Mazur and Tate.

        \note{This function requires that Magma to be installed on your
        computer.}

        INPUT:
            p -- prime
            P -- point
            prec -- integer (default: 20) affects the precision; the
                    precision is *not* guaranteed to be this high!
        OUTPUT:
            p-adic number
        """
        p = self.__check_padic_hypotheses(p)
        if not P in self:
            raise ArithmeticError, "P = (%s) must be a point on this curve"%P
        return padic_height.padic_height(self.a_invariants(), p, P, prec)


    # This is the old version of padic_regulator that requires MAGMA:
    def padic_regulator_magma(self, p, prec=20):
        """
        Return the cyclotomic $p$-adic regulator of $P$, in the sense
        of Mazur and Tate.

        \note{This function requires that Magma to be installed on your
        computer.}

        INPUT:
            p -- prime
            prec -- integer (default: 20) affects the precision; the
                    precision is *not* guaranteed to be this high!
        OUTPUT:
            p-adic number
        """
        p = self.__check_padic_hypotheses(p)
        return padic_height.padic_regulator(self.a_invariants(),
                                            p,
                                            self.gens(),
                                            prec)


    def padic_regulator(self, p, prec=20, height=None, check_hypotheses=True):
        r"""
        Computes the cyclotomic p-adic regulator of this curve.

        This curve must be in minimal weierstrass form.

        INPUT:
            p -- prime >= 5
            prec -- answer will be returned modulo p^prec (unless p is
                    anomalous; see below)
            height -- precomputed height function. If not supplied, this
                 function will call padic_height to compute it.
            check_hypotheses -- boolean, whether to check that this is a
                 curve for which the p-adic height makes sense

        OUTPUT:
            The p-adic cyclotomic regulator of this curve, to the requested
            precision. HOWEVER, if $p$ is anomalous for this curve, then
            there may be some precision loss (at most $2r$ p-adic digits,
            where $r$ is the rank of the curve). This is caused by the fact
            that the matrix of height pairings may contain some denominators
            (even though each entry is computed to the correct precision).

            If the rank is 0, we output 1.

        TODO:
            -- remove restriction that curve must be in minimal weierstrass
            form. This is currently required for E.gens().

        AUTHORS:
            -- Liang Xiao, original implementation at the 2006 MSRI graduate
            workshop on modular forms
            -- David Harvey (2006-09-13), cleaned up and integrated into SAGE,
            removed some redundant height computations

        EXAMPLES:
            sage: E = EllipticCurve("37a")
            sage: E.padic_regulator(5, 10)
            1 + 5 + 5^2 + 3*5^5 + 4*5^6 + 5^8 + 5^9 + O(5^10)

        A rank zero example:
            sage: EllipticCurve('11a').padic_regulator(3)
            1 + O(3^20)

        An anomalous case:
            sage: E.padic_regulator(53, 10)
            26*53^-2 + 30*53^-1 + 20 + 47*53 + 10*53^2 + 32*53^3 + 9*53^4 + 22*53^5 + 35*53^6 + 30*53^7 + O(53^8)

        An anomalous case where the precision drops some:
            sage: E = EllipticCurve("5077a")
            sage: E.padic_regulator(5, 10)
            5^-2 + 5^-1 + 4 + 2*5 + 2*5^2 + 2*5^3 + 4*5^4 + 2*5^5 + 5^6 + O(5^7)

        Check that answers agree over a range of precisions:
            sage: max_prec = 30    # make sure we get past p^2    # long time
            sage: full = E.padic_regulator(5, max_prec)           # long time
            sage: for prec in range(1, max_prec):                 # long time
            ...       assert E.padic_regulator(5, prec) == full   # long time

        """
        d = self.padic_height_pairing_matrix(p=p, prec=prec,
                              height=height, check_hypotheses=check_hypotheses)
        if d.nrows() == 0:
            return d.base_ring()(1)
        return d.determinant()

    def padic_height_pairing_matrix(self, p, prec=20, height=None, check_hypotheses=True):
        r"""
        Computes the cyclotomic $p$-adic height pairing matrix of this
        curve with respect to the basis self.gens() for the
        Mordell-Weil group.

        This curve must be in minimal weierstrass form.

        INPUT:
            p -- prime >= 5
            prec -- answer will be returned modulo p^prec (unless p is
                    anomalous; see below)
            height -- precomputed height function. If not supplied, this
                 function will call padic_height to compute it.
            check_hypotheses -- boolean, whether to check that this is a
                 curve for which the p-adic height makes sense

        OUTPUT:
            The p-adic cyclotomic height pairing matrix of this curve
            to the given precision. HOWEVER, if $p$ is anomalous for
            this curve, then there may be some precision loss (at most
            $2r$ p-adic digits, where $r$ is the rank of the
            curve).

        TODO:
            -- remove restriction that curve must be in minimal weierstrass
            form. This is currently required for E.gens().

        AUTHORS:
            -- David Harvey, Liang Xiao, Robert Bradshaw, Jennifer
            Balakrishnan, original implementation at the 2006 MSRI
            graduate workshop on modular forms
            -- David Harvey (2006-09-13), cleaned up and integrated into SAGE,
            removed some redundant height computations

        EXAMPLES:
            sage: E = EllipticCurve("37a")
            sage: E.padic_height_pairing_matrix(5, 10)
            [1 + 5 + 5^2 + 3*5^5 + 4*5^6 + 5^8 + 5^9 + O(5^10)]

        A rank two example:
            E = sage: EllipticCurve('389a').padic_height_pairing_matrix(5,4)
            [    3 + 2*5 + 5^3 + O(5^4) 3 + 5 + 5^2 + 5^3 + O(5^4)]
            [3 + 5 + 5^2 + 5^3 + O(5^4)     2*5^2 + 3*5^3 + O(5^4)]

        An anomalous rank 3 example:
            sage: E = EllipticCurve("5077a")
            sage: EllipticCurve('5077a').padic_height_pairing_matrix(5,2)   # somewhat random precision depending on architecture.
            [4*5^-1 + 2 + 4*5 + O(5^2)              4*5 + O(5^2)       4*5^-1 + 1 + O(5^2)]
            [             4*5 + O(5^2) 4*5^-1 + 3 + 4*5 + O(5^2)       2*5^-1 + 5 + O(5^2)]
            [      4*5^-1 + 1 + O(5^2)       2*5^-1 + 5 + O(5^2)              2*5 + O(5^2)]
        """
        if check_hypotheses:
            p = self.__check_padic_hypotheses(p)

        K = Qp(p, prec=prec)

        rank = self.rank()
        M = matrix.matrix(K, rank, rank, 0)
        if rank == 0:
            return M

        basis = self.gens()

        if height is None:
            height = self.padic_height(p, prec, check_hypotheses=False)

        # Use <P, Q> = h(P) + h(Q) - h(P + Q)

        point_height = [height(P) for P in basis]
        for i in range(rank):
            for j in range(i, rank):
                M[i, j] = M[j, i] = point_height[i] + point_height[j] \
                                    - height(basis[i] + basis[j])
        return M



    class _DivPolyContext:
        r"""
        This class implements the algorithm in Section 3 of "Efficient
        Computation of p-adic Heights" (David Harvey, still in draft form).

        The constructor takes as input an elliptic curve $E/\QQ$ with integer
        coefficients, a rational point $P$ that reduces to a non-singular
        point at all primes, and a ring $R$ in which $2$ is invertible.
        Typically $R$ will be $\ZZ/R\ZZ$ for some positive integer $R$.

        One then calls triple(m) for an integer $m \geq 1$; it returns a
        triple $(a', b', d')$ such that if the point $mQ$ has coordinates
        $(a/d^2, b/d^3)$, then we have $a' \equiv a$, $b' \equiv \pm b$,
        $d' \equiv \pm d$ all modulo $R$.

        Note the ambiguity of signs for $b'$ and $d'$. There's not much one
        can do about this, but at least one can say that the sign for $b'$
        will match the sign for $d'$.

        Complexity is soft $O(\log R \log^2 m)$.

        AUTHOR:
            -- David Harvey (2007-02)

        EXAMPLES:

        37a has trivial tamagawa numbers so all points have nonsingular
        reduction at all primes:
            sage: E = EllipticCurve("37a")
            sage: P = E.gens()[0]; P
             (0 : 0 : 1)
            sage: 19*P
             (-59997896/67387681 : -641260644409/553185473329 : 1)
            sage: R = Integers(625)
            sage: C = E._DivPolyContext(E, R, P)
            sage: C.triple(19)
             (229, 34, 541)
            sage: -59997896 % 625
             229
            sage: 67387681.sqrt()
             8209
            sage: -8209 % 625          # note sign is flipped
             541
            sage: 641260644409 % 625   # sign flipped here too
             34

        Test over a range of $n$ for a single curve with fairly random
        coefficients:
            sage: R = Integers(625)
            sage: E = EllipticCurve([4, -11, 17, -8, -10])
            sage: P = E.gens()[0] * LCM(E.tamagawa_numbers())
            sage: C = E._DivPolyContext(E, R, P)
            sage: Q = E(0)
            sage: for n in range(1, 25):
            ...      Q = Q + P
            ...      naive = R(Q[0].numerator()),  \
            ...              R(Q[1].numerator()),  \
            ...              R(Q[0].denominator().sqrt())
            ...      triple = C.triple(n)
            ...      assert (triple[0] == naive[0]) and ( \
            ...        (triple[1] == naive[1] and triple[2] == naive[2]) or \
            ...        (triple[1] == -naive[1] and triple[2] == -naive[2])), \
            ...           "_DivPolyContext.triple() gave an incorrect answer"

        """
        def __init__(self, E, R, P):
            alpha = R(P[0].numerator())
            beta = R(P[1].numerator())
            d = R(P[0].denominator().sqrt())

            a1 = R(E.a1()) * d
            a3 = R(E.a3()) * d**3

            b2 = R(E.b2()) * d**2
            b4 = R(E.b4()) * d**4
            b6 = R(E.b6()) * d**6
            b8 = R(E.b8()) * d**8

            B4 = 6*alpha**2 + b2*alpha + b4
            B6 = 4*alpha**3 + b2*alpha**2 + 2*b4*alpha + b6
            B8 = 3*alpha**4 + b2*alpha**3 + 3*b4*alpha**2 + 3*b6*alpha + b8

            self.E = E
            self.R = R
            self.alpha = alpha
            self.beta = beta
            self.d = d
            self.a1 = a1
            self.a3 = a3
            self.b2 = b2
            self.b4 = b4
            self.b6 = b6
            self.b8 = b8
            self.B4 = B4
            self.B6 = B6
            self.B6_sqr = B6*B6
            self.B8 = B8
            self.T = 2*beta + a1*alpha + a3

            self.g_cache = \
                  {0 : R(0), 1 : R(1), 2 : R(-1), 3 : B8, 4 : B6**2 - B4*B8}
            self.psi_cache = {}
            self.theta_cache = {}
            self.omega_cache = {1 : beta}


        def g(self, n):
            r"""
            Returns $\hat g_n(P)$ as defined in the paper mentioned above.
            """
            # try to get cached value
            try:
                return self.g_cache[n]
            except KeyError:
                pass

            # didn't work, have to compute it
            g = self.g
            m = n // 2
            if n & 1:
                prod1 = g(m+2) * g(m)**3
                prod2 = g(m-1) * g(m+1)**3

                if m & 1:
                    X = prod1 - self.B6_sqr * prod2
                else:
                    X = self.B6_sqr * prod1 - prod2
            else:
                X = g(m) * (g(m-2) * g(m+1)**2 - g(m+2) * g(m-1)**2)

            self.g_cache[n] = X
            return X


        def psi(self, n):
            r"""
            Returns $\hat \psi_n(P)$ as defined in the paper mentioned above.
            """
            # try to get cached value
            try:
                return self.psi_cache[n]
            except KeyError:
                assert n >= 0

            # didn't work, have to compute it
            X = self.g(n)
            if n & 1 == 0:
                X = X * self.T

            self.psi_cache[n] = X
            return X


        def theta(self, n):
            r"""
            Returns $\hat \theta g_n(P)$ as defined in the paper mentioned above.
            """
            # try to get cached value
            try:
                return self.theta_cache[n]
            except KeyError:
                assert n >= 1

            # didn't work, have to compute it
            X = self.alpha * self.psi(n)**2 - self.psi(n-1) * self.psi(n+1)
            self.theta_cache[n] = X
            return X


        def omega(self, n):
            r"""
            Returns $\hat \omega g_n(P)$ as defined in the paper mentioned above.
            """
            # try to get cached value
            try:
                return self.omega_cache[n]
            except KeyError:
                assert n >= 2

            # didn't work, have to compute it
            X = self.g(n-2) * self.g(n+1)**2 - self.g(n+2) * self.g(n-1)**2
            if n & 1 == 1:
                X = X * self.T
            X = (X + (self.a1 * self.theta(n) + self.a3 * self.psi(n)**2) \
                                                          * self.psi(n)) / -2

            self.omega_cache[n] = X
            return X


        def triple(self, n):
            assert n >= 1
            return self.theta(n), self.omega(n), self.psi(n) * self.d



    def padic_height(self, p, prec=20, sigma=None, check_hypotheses=True):
        r"""
        Computes the cyclotomic p-adic height, as defined by Mazur and Tate.
        The height is normalised to take values in $\Z_p$, unless $p$ is
        anomalous in which case it takes values in $(1/p^2)\Z_p$.

        The equation of the curve must be minimal at $p$.

        INPUT:
            p -- prime >= 5 for which the curve has good ordinary reduction
            prec -- integer >= 1, desired precision of result
            sigma -- precomputed value of sigma. If not supplied, this function
                 will call padic_sigma to compute it.
            check_hypotheses -- boolean, whether to check that this is a
                 curve for which the p-adic height makes sense

        OUTPUT:
            A function that accepts two parameters:
              * a Q-rational point on the curve whose height should be computed
              * optional boolean flag "check": if False, it skips some
                input checking,
            and returns the p-adic height of that point to the desired
            precision.

        AUTHORS:
            -- Jennifer Balakrishnan: original code developed at the 2006
            MSRI graduate workshop on modular forms
            -- David Harvey (2006-09-13): integrated into SAGE, optimised
            to speed up repeated evaluations of the returned height function,
            addressed some thorny precision questions
            -- David Harvey (2006-09-30): rewrote to use division polynomials
            for computing denominator of $nP$.
            -- David Harvey (2007-02): cleaned up according to algorithms
            in "Efficient Computation of p-adic Heights"

        TODO:
            -- Probably this code is broken when P is a torsion point.

        EXAMPLES:
            sage: E = EllipticCurve("37a")
            sage: P = E.gens()[0]
            sage: h = E.padic_height(5, 10)
            sage: h(P)
             2 + 4*5 + 5^2 + 2*5^3 + 2*5^4 + 3*5^5 + 2*5^6 + 4*5^7 + 5^8 + 4*5^9 + O(5^10)

          An anomalous case:
            sage: h = E.padic_height(53, 10)
            sage: h(P)
             40*53^-2 + 37*53^-1 + 42 + 2*53 + 21*53^2 + 10*53^3 + 48*53^4 + 41*53^5 + 8*53^6 + 11*53^7 + O(53^8)

          Boundary case:
            sage: E.padic_height(5, 1)(P)
             2 + O(5)

          A case that works the division polynomial code a little harder:
            sage: E.padic_height(5, 10)(5*P)
            2*5^2 + 4*5^3 + 5^4 + 2*5^5 + 2*5^6 + 3*5^7 + O(5^8)

          Check that answers agree over a range of precisions:
            sage: max_prec = 30    # make sure we get past p^2    # long time
            sage: full = E.padic_height(5, max_prec)(P)           # long time
            sage: for prec in range(1, max_prec):                 # long time
            ...       assert E.padic_height(5, prec)(P) == full   # long time

        """
        if check_hypotheses:
            p = self.__check_padic_hypotheses(p)

        prec = int(prec)
        if prec < 1:
            raise ValueError, "prec (=%s) must be at least 1" % prec

        # For notation and definitions, see "Efficient Computation of
        # p-adic Heights", David Harvey (unpublished)

        n1 = self.change_ring(rings.GF(p)).cardinality()
        n2 = arith.LCM(self.tamagawa_numbers())
        n = arith.LCM(n1, n2)
        m = int(n / n2)

        adjusted_prec = prec + 2 * arith.valuation(n, p) + 1
        R = rings.Integers(p ** adjusted_prec)

        if sigma is None:
            sigma = self.padic_sigma(p, adjusted_prec, check_hypotheses=False)

        # K is the field for the final result
        K = Qp(p, prec=prec)
        E = self


        def height(P, check=True):
            if check:
                assert P.curve() == E, "the point P must lie on the curve " \
                       "from which the height function was created"

            Q = n2 * P
            C = E._DivPolyContext(E, R, Q)

            alpha, beta, d = C.triple(m)

            assert beta.lift() % p != 0, "beta should be a unit!"
            assert d.lift() % p == 0, "d should not be a unit!"

            t = -d * alpha / beta

            total = R(1)
            t_power = t
            for k in range(2, adjusted_prec + 1):
                # yuck... should just be able to multiply without the lift here
                total = total + t_power * R(sigma[k].lift())
                t_power = t_power * t

            L = Qp(p, prec=adjusted_prec)
            total = (-alpha / beta) * total
            total = L(total.lift(), adjusted_prec)   # yuck... get rid of this lift!
            answer = total.log() / n**2 / p

            if check:
                assert answer.precision_absolute() >= prec, "we should have got an " \
                       "answer with precision at least prec, but we didn't."
            return K(answer.lift(), prec - answer.valuation())


        # (man... I love python's local function definitions...)
        return height


    def padic_sigma(self, p, N=20, E2=None, check=False, check_hypotheses=True):
        r"""
        Computes the p-adic sigma function with respect to the standard
        invariant differential $dx/(2y + a_1 x + a_3)$, as defined by
        Mazur and Tate, as a power series in the usual uniformiser $t$ at the
        origin.

        The equation of the curve must be minimal at $p$.

        INPUT:
            p -- prime >= 5 for which the curve has good ordinary reduction
            N -- integer >= 1, indicates precision of result; see OUTPUT
                 section for description
            E2 -- precomputed value of E2. If not supplied, this function will
                 call padic_E2 to compute it. The value supplied must be
                 correct mod $p^{N-2}$.
            check -- boolean, whether to perform a consistency check (i.e.
                 verify that the computed sigma satisfies the defining
                 differential equation -- note that this does NOT guarantee
                 correctness of all the returned digits, but it comes pretty
                 close :-))
            check_hypotheses -- boolean, whether to check that this is a
                 curve for which the p-adic sigma function makes sense

        OUTPUT:
            A power series $t + \cdots$ with coefficients in $\Z_p$.

            The output series will be truncated at $O(t^{N+1})$, and the
            coefficient of $t^n$ for $n \geq 1$ will be correct to precision
            $O(p^{N-n+1})$.

            In practice this means the following. If $t_0 = p^k u$, where
            $u$ is a $p$-adic unit with at least $N$ digits of precision,
            and $k \geq 1$, then the returned series may be used to compute
            $\sigma(t_0)$ correctly modulo $p^{N+k}$ (i.e. with $N$ correct
            $p$-adic digits).

        ALGORITHM:
           Described in ``Efficient Computation of p-adic Heights''
           (David Harvey), which is basically an optimised version of the
           algorithm from ``p-adic Heights and Log Convergence'' (Mazur,
           Stein, Tate).

           Running time is soft-$O(N^2 \log p)$, plus whatever time is
           necessary to compute $E_2$.

        AUTHOR:
            -- David Harvey (2006-09-12)
            -- David Harvey (2007-02): rewrote

        EXAMPLES:
            sage: EllipticCurve([-1, 1/4]).padic_sigma(5, 10)
            (1 + O(5^20))*t + (3 + 2*5^2 + 3*5^3 + 3*5^6 + 4*5^7 + O(5^8))*t^3 + (2 + 4*5^2 + 4*5^3 + 5^4 + 5^5 + O(5^6))*t^5 + (2 + 2*5 + 5^2 + 4*5^3 + O(5^4))*t^7 + (1 + 2*5 + O(5^2))*t^9 + O(t^11)

          Run it with a consistency check:
            sage: EllipticCurve("37a").padic_sigma(5, 10, check=True)
            (1 + O(5^20))*t + (3 + 2*5^2 + 3*5^3 + 3*5^6 + 4*5^7 + O(5^8))*t^3 + (3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + O(5^7))*t^4 + (2 + 4*5^2 + 4*5^3 + 5^4 + 5^5 + O(5^6))*t^5 + (2 + 3*5 + 5^4 + O(5^5))*t^6 + (4 + 3*5 + 2*5^2 + O(5^4))*t^7 + (2 + 3*5 + 2*5^2 + O(5^3))*t^8 + (4*5 + O(5^2))*t^9 + (1 + O(5))*t^10 + O(t^11)

          Boundary cases:
            sage: EllipticCurve([1, 1, 1, 1, 1]).padic_sigma(5, 1)
             (1 + O(5^20))*t + O(t^2)
            sage: EllipticCurve([1, 1, 1, 1, 1]).padic_sigma(5, 2)
             (1 + O(5^20))*t + (3 + O(5))*t^2 + O(t^3)

          Check that sigma is ``weight 1''. [This test is disabled until
          trac \#254 is addressed. The lines f(2*t)/2 and g should return
          exactly the same answer. Currently there is some precision loss.]
            sage.: f = EllipticCurve([-1, 3]).padic_sigma(5, 10)
            sage.: g = EllipticCurve([-1*(2**4), 3*(2**6)]).padic_sigma(5, 10)
            sage.: t = f.parent().gen()
            sage.: f(2*t)/2
            sage.: g
             (1 + O(5^20))*t + (4 + 3*5 + 3*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 3*5^6 + 5^7 + 2*5^8 + O(5^9))*t^3 + (3 + 3*5^2 + 5^4 + 2*5^5 + 3*5^6 + 2*5^7 + 3*5^8 + O(5^9))*t^5 + (4 + 5 + 3*5^3 + 2*5^4 + 2*5^5 + 5^6 + O(5^7))*t^7 + (4 + 2*5 + 4*5^2 + 2*5^4 + 4*5^5 + 5^6 + O(5^7))*t^9 + O(t^11)

          Test that it returns consistent results over a range of precision:
          [NOTE: this test currently FAILS due to trac \#255. It seems to be
          essentially correct though. Should be revisited when that bug is
          resolved.]
            sage: max_N = 30   # get up to at least p^2         # long time
            sage: E = EllipticCurve([1, 1, 1, 1, 1])            # long time
            sage: p = 5                                         # long time
            sage: E2 = E.padic_E2(5, max_N)                     # long time
            sage: max_sigma = E.padic_sigma(p, max_N, E2=E2)    # long time
            sage: for N in range(3, max_N):                     # long time
            ...      sigma = E.padic_sigma(p, N, E2=E2)         # long time
            ...      for n in range(sigma.prec()):              # long time
            ...          assert sigma[n] == max_sigma[n]        # long time

        """
        if check_hypotheses:
            p = self.__check_padic_hypotheses(p)

        # todo: implement the p == 3 case
        # NOTE: If we ever implement p == 3, it's necessary to check over
        # the precision loss estimates (below) vey carefully; I think it
        # may become necessary to compute E2 to an even higher precision.
        if p < 5:
            raise NotImplementedError, "p (=%s) must be at least 5" % p

        N = int(N)
        if N <= 2:
            # a few special cases for small N
            if N < 1:
                raise ValueError, "N (=%s) must be at least 1" % prec

            K = Qp(p)

            if N == 1:
                # return simply t + O(t^2)
                return PowerSeriesRing(K, "t")([K(0), K(1)], prec=2)

            if N == 2:
                # return t + a_1/2 t^2 + O(t^3)
                return PowerSeriesRing(K, "t")([K(0), K(1),
                                                K(self.a1()/2, 1)], prec=3)

        if self.discriminant().valuation(p) != 0:
            raise NotImplementedError, "equation of curve must be minimal at p"

        if E2 is None:
            E2 = self.padic_E2(p, N-2, check_hypotheses=False)
        elif E2.big_oh() < N-2:
            raise ValueError, "supplied E2 has insufficient precision"

        QQt = LaurentSeriesRing(RationalField(), "x")

        R = rings.Integers(p**(N-2))
        X = self.change_ring(R)
        c = (X.a1()**2 + 4*X.a2() - R(E2)) / 12

        f = X.formal_group().differential(N+2)   # f = 1 + ... + O(t^{N+2})
        x = X.formal_group().x(N)                # x = t^{-2} + ... + O(t^N)

        Rt = x.parent()

        A  = (x + c) * f
        # do integral over QQ, to avoid divisions by p
        A = Rt(QQt(A).integral())
        A = (-X.a1()/2 - A) * f

        # Convert to a power series and remove the -1/x term.
        # Also we artifically bump up the accuracy from N-2 to to N-1 digits;
        # the constant term needs to be known to N-1 digits, so we compute
        # it directly
        assert A.valuation() == -1 and A[-1] == 1
        A = A - A.parent().gen() ** (-1)
        A = A.power_series().list()
        R = rings.Integers(p**(N-1))
        A = [R(u) for u in A]
        A[0] = self.change_ring(R).a1()/2     # fix constant term
        A = PowerSeriesRing(R, "x")(A, len(A))

        theta = _brent(A, p, N)
        sigma = theta * theta.parent().gen()

        # Convert the answer to power series over p-adics; drop the precision
        # of the $t^k$ coefficient to $p^(N-k+1)$.
        # [Note: there are actually more digits available, but it's a bit
        # tricky to figure out exactly how many, and we only need $p^(N-k+1)$
        # for p-adic height purposes anyway]
        K = rings.pAdicField(p)
        sigma = sigma.padded_list(N+1)
        sigma[0] = K(0)
        sigma[1] = K(1)
        for n in range(2, N+1):
            sigma[n] = K(sigma[n].lift(), N - n + 1)

        S = rings.PowerSeriesRing(K, "t", N+1)
        sigma = S(sigma, N+1)

        # if requested, check that sigma satisfies the appropriate
        # differential equation
        if check:
            R = rings.Integers(p**N)
            X = self.change_ring(R)
            x = X.formal_group().x(N+5)       # few extra terms for safety
            f = X.formal_group().differential(N+5)
            c = (X.a1()**2 + 4*X.a2() - R(E2)) / 12

            # convert sigma to be over Z/p^N
            s = f.parent()(sigma)

            # apply differential equation
            temp = (s.derivative() / s / f).derivative() / f + c + x

            # coefficient of t^k in the result should be zero mod p^(N-k-2)
            for k in range(N-2):
                assert temp[k].lift().valuation(p) >= N - k - 2, \
                            "sigma correctness check failed!"

        return sigma



    def padic_E2(self, p, prec=20, check=False, check_hypotheses=True,
                 algorithm="auto"):
        r"""
        Returns the value of the $p$-adic modular form $E2$ for $(E, \omega)$
        where $\omega$ is the usual invariant differential
        $dx/(2y + a_1 x + a_3)$.

        INPUT:
            p -- prime (>= 5) for which $E$ is good and ordinary
            prec -- (relative) p-adic precision (>= 1) for result
            check -- boolean, whether to perform a consistency check.
                 This will slow down the computation by a constant factor < 2.
                 (The consistency check is to compute the whole matrix of
                 frobenius on Monsky-Washnitzer cohomology, and verify that
                 its trace is correct to the specified precision. Otherwise,
                 the trace is used to compute one column from the other one
                 (possibly after a change of basis).)
            check_hypotheses -- boolean, whether to check that this is a
                 curve for which the p-adic sigma function makes sense
            algorithm -- one of "standard", "sqrtp", or "auto". This
                 selects which version of Kedlaya's algorithm is used. The
                 "standard" one is the one described in Kedlaya's paper.
                 The "sqrtp" one has better performance for large $p$;
                 however it is not as well tested yet, and might not work
                 properly when $p$ is too small. The "auto" option selects
                 "sqrtp" when $p$ is at least 3000, and "standard" otherwise.
                 Note that if the "sqrtp" algorithm is used, a consistency
                 check will automatically be applied, regardless of the
                 setting of the "check" flag.

        OUTPUT:
            p-adic number to precision prec

        NOTES:
            -- If the discriminant of the curve has nonzero valuation at p,
               then the result will not be returned mod $p^\text{prec}$, but
               it still *will* have prec *digits* of precision.

        AUTHOR:
            -- David Harvey (2006-09-01): partly based on code written by
               Robert Bradshaw at the MSRI 2006 modular forms workshop

        ACKNOWLEDGMENT:
           -- discussion with Eyal Goren that led to the trace trick.

        EXAMPLES:
        Here is the example discussed in the paper ``Computation of p-adic
        Heights and Log Convergence'' (Mazur, Stein, Tate):
            sage: EllipticCurve([-1, 1/4]).padic_E2(5)
            2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + O(5^20)

        Let's try to higher precision (this is the same answer the MAGMA
        implementation gives):
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 100)
            2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + 4*5^20 + 5^21 + 4*5^22 + 2*5^23 + 3*5^24 + 3*5^26 + 2*5^27 + 3*5^28 + 2*5^30 + 5^31 + 4*5^33 + 3*5^34 + 4*5^35 + 5^36 + 4*5^37 + 4*5^38 + 3*5^39 + 4*5^41 + 2*5^42 + 3*5^43 + 2*5^44 + 2*5^48 + 3*5^49 + 4*5^50 + 2*5^51 + 5^52 + 4*5^53 + 4*5^54 + 3*5^55 + 2*5^56 + 3*5^57 + 4*5^58 + 4*5^59 + 5^60 + 3*5^61 + 5^62 + 4*5^63 + 5^65 + 3*5^66 + 2*5^67 + 5^69 + 2*5^70 + 3*5^71 + 3*5^72 + 5^74 + 5^75 + 5^76 + 3*5^77 + 4*5^78 + 4*5^79 + 2*5^80 + 3*5^81 + 5^82 + 5^83 + 4*5^84 + 3*5^85 + 2*5^86 + 3*5^87 + 5^88 + 2*5^89 + 4*5^90 + 4*5^92 + 3*5^93 + 4*5^94 + 3*5^95 + 2*5^96 + 4*5^97 + 4*5^98 + 2*5^99 + O(5^100)

        Check it works at low precision too:
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 1)
            2 + O(5)
            sage: EllipticCurve([1, 1, 1, 1, 1]).padic_E2(5, 1)
            5 + O(5^2)
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 2)
            2 + 4*5 + O(5^2)
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 3)
            2 + 4*5 + O(5^3)

          Check it works for different models of the same curve (37a),
          even when the discriminant changes by a power of p (note that
          E2 depends on the differential too, which is why it gets scaled
          in some of the examples below):

            sage: X1 = EllipticCurve([-1, 1/4])
            sage: X1.j_invariant(), X1.discriminant()
             (110592/37, 37)
            sage: X1.padic_E2(5, 10)
             2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

            sage: X2 = EllipticCurve([0, 0, 1, -1, 0])
            sage: X2.j_invariant(), X2.discriminant()
             (110592/37, 37)
            sage: X2.padic_E2(5, 10)
             2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

            sage: X3 = EllipticCurve([-1*(2**4), 1/4*(2**6)])
            sage: X3.j_invariant(), X3.discriminant() / 2**12
             (110592/37, 37)
            sage: 2**(-2) * X3.padic_E2(5, 10)
             2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

            sage: X4 = EllipticCurve([-1*(7**4), 1/4*(7**6)])
            sage: X4.j_invariant(), X4.discriminant() / 7**12
             (110592/37, 37)
            sage: 7**(-2) * X4.padic_E2(5, 10)
             2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

            sage: X5 = EllipticCurve([-1*(5**4), 1/4*(5**6)])
            sage: X5.j_invariant(), X5.discriminant() / 5**12
             (110592/37, 37)
            sage: 5**(-2) * X5.padic_E2(5, 10)
             2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

            sage: X6 = EllipticCurve([-1/(5**4), 1/4/(5**6)])
            sage: X6.j_invariant(), X6.discriminant() * 5**12
             (110592/37, 37)
            sage: 5**2 * X6.padic_E2(5, 10)
             2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + O(5^10)

        Test check=True vs check=False:
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 1, check=False)
            2 + O(5)
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 1, check=True)
            2 + O(5)
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 30, check=False)
            2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + 4*5^20 + 5^21 + 4*5^22 + 2*5^23 + 3*5^24 + 3*5^26 + 2*5^27 + 3*5^28 + O(5^30)
            sage: EllipticCurve([-1, 1/4]).padic_E2(5, 30, check=True)
            2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11 + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + 4*5^20 + 5^21 + 4*5^22 + 2*5^23 + 3*5^24 + 3*5^26 + 2*5^27 + 3*5^28 + O(5^30)

        Here's one using the $p^{1/2}$ algorithm:
            sage: EllipticCurve([-1, 1/4]).padic_E2(3001, 3, algorithm="sqrtp")  # long time (< 10s)
            1907 + 2819*3001 + 1124*3001^2 + O(3001^3)

        """
        if check_hypotheses:
            p = self.__check_padic_hypotheses(p)

        # The cutoff p = 3000 is pretty arbitrary. It seems to work well
        # on sage.math for low precision problems. In reality the optimal
        # crossover will depend on the precision, and also on the architecture.
        if algorithm == "auto":
            algorithm = "standard" if p < 3000 else "sqrtp"

        if algorithm not in ["standard", "sqrtp"]:
            raise ValueError, "unknown algorithm '%s'" % algorithm

        # todo: maybe it would be good if default prec was None, and then
        # it selects an appropriate precision based on how large the prime
        # is

        # todo: implement the p == 3 case
        if p < 5:
            raise NotImplementedError, "p (=%s) must be at least 5" % p

        prec = int(prec)
        if prec < 1:
            raise ValueError, "prec (=%s) must be at least 1" % prec

        # To run matrix_of_frobenius(), we need to have the equation in the
        # form y^2 = x^3 + ax + b, whose discriminant is invertible mod p.
        # When we change coordinates like this, we might scale the invariant
        # differential, so we need to account for this. We do this by
        # comparing discriminants: if the discrimimants differ by u^12,
        # then the differentials differ by u. There's a sign ambiguity here,
        # but it doesn't matter because E2 changes by u^2 :-)

        # todo: the weierstrass_model() function is overkill here, because
        # it finds a *minimal* model. I imagine it has to do some factoring
        # to do this. We only need a model that's minimal at p.

        # todo: In fact, there should be available a function that returns
        # exactly *which* coordinate change was used. If I had that I could
        # avoid the discriminant circus at the end.

        # todo: The following strategy won't work at all for p = 2, 3.

        X = self.weierstrass_model()

        assert X.discriminant().valuation(p) == 0, "Something's gone wrong. " \
               "The discriminant of the weierstrass model should be a unit " \
               " at p."

        if algorithm == "standard":
            # Need to increase precision a little to compensate for precision
            # losses during the computation. (See monsky_washnitzer.py
            # for more details.)
            adjusted_prec = monsky_washnitzer.adjusted_prec(p, prec)

            if check:
                trace = None
            else:
                trace = self.ap(p)

            base_ring = rings.Integers(p**adjusted_prec)
            output_ring = rings.pAdicField(p, prec)

            R, x = rings.PolynomialRing(base_ring, 'x').objgen()
            Q = x**3 + base_ring(X.a4()) * x + base_ring(X.a6())
            frob_p = monsky_washnitzer.matrix_of_frobenius(
                             Q, p, adjusted_prec, trace)

        else:   # algorithm == "sqrtp"
            base_ring = rings.Integers(p**prec)
            output_ring = rings.pAdicField(p, prec)
            frob_p = monsky_washnitzer.matrix_of_frobenius_alternate(
                          X.a4(), X.a6(), p, prec)

            # let's force a trace-check since this algorithm is fairly new
            # and we don't want to get caught with our pants down...
            trace = self.ap(p)
            check = True


        if check:
            trace_of_frobenius = frob_p.trace().lift() % p**prec
            correct_trace = self.ap(p) % p**prec
            assert trace_of_frobenius == correct_trace, \
                    "Consistency check failed! (correct = %s, actual = %s)" % \
                    (correct_trace, trace_of_frobenius)

        frob_p_n = frob_p**prec

        # todo: write a coercion operator from integers mod p^n to the
        # p-adic field (it doesn't seem to currently exist)
        # see trac #4

        # todo: think about the sign of this. Is it correct?

        E2_of_X = output_ring( (-12 * frob_p_n[0,1] / frob_p_n[1,1]).lift() ) \
                  + O(p**prec)

        # Take into account the coordinate change.
        fudge_factor = (X.discriminant() / self.discriminant()).nth_root(6)

        # todo: here I should be able to write:
        #  return E2_of_X / fudge_factor
        # However, there is a bug in SAGE (#51 on trac) which makes this
        # crash sometimes when prec == 1. For example,
        #    EllipticCurve([1, 1, 1, 1, 1]).padic_E2(5, 1)
        # makes it crash. I haven't figured out exactly what the bug
        # is yet, but for now I use the following workaround:
        fudge_factor_inverse = Qp(p, prec=(E2_of_X.precision_absolute() + 1))(1 / fudge_factor)
        return output_ring(E2_of_X * fudge_factor_inverse)


    # This is the old version of padic_E2 that requires MAGMA:
    def padic_E2_magma(self, p, prec=20):
        """
        Return the value of the $p$-adic.
        """
        p = self.__check_padic_hypotheses(p)
        c4, c6 = self.c_invariants()
        return padic_height.padic_E2_of_c4c6(c4, c6, p, prec)


    # 	def	weierstrass_p(self):
    #         # TODO: add allowing negative valuations for power series
    #         return 1/t**2 + a1/t + rings.frac(1,12)*(a1-8*a2) -a3*t \
    #                - (a4+a1*a3)*t**2  + O(t**3)


    def mod5family(self):
        """
        Return the family of all elliptic curves with the same mod-5
        representation as self.
        """
        E = self.weierstrass_model()
        a = E.a4()
        b = E.a6()
        return mod5family.mod5family(a,b)


def cremona_curves(conductors):
    """
    Return iterator over all known curves (in database) with conductor
    in the list of conductors.
    """
    if isinstance(conductors, (int,long, rings.RingElement)):
        conductors = [conductors]
    return sage.databases.cremona.CremonaDatabase().iter(conductors)

def cremona_optimal_curves(conductors):
    """
    Return iterator over all known optimal curves (in database) with
    conductor in the list of conductors.
    """
    if isinstance(conductors, (int,long,rings.RingElement)):
        conductors = [conductors]
    return sage.databases.cremona.CremonaDatabase().iter_optimal(conductors)

def _brent(F, p, N):
    r"""
    This is an internal function; it is used by padic_sigma().

    $F$ is a assumed to be a power series over $R = \Z/p^{N-1}\Z$.

    It solves the differential equation $G'(t)/G(t) = F(t)$ using Brent's
    algorithm, with initial condition $G(0) = 1$. It is assumed that the
    solution $G$ has $p$-integral coefficients.

    More precisely, suppose that $f(t)$ is a power series with genuine
    $p$-adic coefficients, and suppose that $g(t)$ is an exact solution to
    $g'(t)/g(t) = f(t)$. Let $I$ be the ideal $(p^N, p^{N-1} t, \ldots,
    p t^{N-1}, t^N)$. The input $F(t)$ should be a finite-precision
    approximation to $f(t)$, in the sense that $\int (F - f) dt$ should lie
    in $I$. Then the function returns a series $G(t)$ such that $(G - g)(t)$
    lies in $I$.

    Complexity should be about $O(N^2 \log^2 N \log p)$, plus some log-log
    factors.

    For more information, and a proof of the precision guarantees,
    see Lemma 4 in ``Efficient Computation of p-adic Heights'' (David
    Harvey).

    AUTHOR:
        -- David Harvey (2007-02)

    EXAMPLES:
    Carefully test the precision guarantees:
        sage: brent = sage.schemes.elliptic_curves.ell_rational_field._brent
        sage: for p in [2, 3, 5]:
        ...     for N in [2, 3, 10, 50]:
        ...       R = Integers(p**(N-1))
        ...       Rx, x = PowerSeriesRing(R, "x").objgen()
        ...       for _ in range(5):
        ...         g = [R.random_element() for i in range(N)]
        ...         g[0] = R(1)
        ...         g = Rx(g, len(g))
        ...         f = g.derivative() / g
        ...         # perturb f by something whose integral is in I
        ...         err = [R.random_element() * p**(N-i) for i in range(N+1)]
        ...         err = Rx(err, len(err))
        ...         err = err.derivative()
        ...         F = f + err
        ...         # run brent() and compare output modulo I
        ...         G = brent(F, p, N)
        ...         assert G.prec() >= N, "not enough output terms"
        ...         err = (G - g).list()
        ...         for i in range(len(err)):
        ...           assert err[i].lift().valuation(p) >= (N - i), \
        ...                  "incorrect precision output"

    """
    Rx = F.parent()           # Rx = power series ring over Z/p^{N-1} Z
    R = Rx.base_ring()        # R = Z/p^{N-1} Z
    Qx = PowerSeriesRing(RationalField(), "x")

    # initial approximation:
    G = Rx(1)

    # loop over an appropriate increasing sequence of lengths s
    for s in misc.newton_method_sizes(N):
        # zero-extend to s terms
        # todo: there has to be a better way in SAGE to do this...
        G = Rx(G.list(), s)

        # extend current approximation to be correct to s terms
        H = G.derivative() / G - F
        # Do the integral of H over QQ[x] to avoid division by p problems
        H = Rx(Qx(H).integral())
        G = G * (1 - H)

    return G

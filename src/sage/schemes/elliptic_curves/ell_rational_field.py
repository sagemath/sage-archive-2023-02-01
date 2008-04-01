"""
Elliptic curves over the rational numbers

AUTHORS:
   -- William Stein (2005): first version
   -- William Stein (2006-02-26): fixed Lseries_extended which didn't work
            because of changes elsewhere in SAGE.
   -- David Harvey (2006-09): Added padic_E2, padic_sigma, padic_height,
            padic_regulator methods.
   -- David Harvey (2007-02): reworked padic-height related code
   -- Christian Wuthrich (2007): added padic sha computation
   -- David Roe (2007-9): moved sha, l-series and p-adic functionality to separate files.
   -- John Cremona (2008-01)
"""

#*****************************************************************************
#       Copyright (C) 2005,2006,2007 William Stein <wstein@gmail.com>
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
from ell_generic import EllipticCurve_generic
from ell_number_field import EllipticCurve_number_field

import sage.groups.all
import sage.rings.arith as arith
import sage.rings.all as rings
import sage.rings.number_field.number_field as number_field
import sage.misc.misc as misc
from sage.misc.all import verbose
import sage.functions.constants as constants
import sage.modular.modform.constructor
import sage.modular.modform.element
from sage.misc.functional import log
from sage.rings.padics.factory import Zp, Qp

# Use some interval arithmetic to guarantee correctness.  We assume
# that alpha is computed to the precision of a float.
# IR = rings.RIF
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
import padics

from lseries_ell import Lseries_ell

import mod5family

from sage.rings.all import (
    PowerSeriesRing, LaurentSeriesRing, O,
    infinity as oo,
    Integer,
    Integers,
    IntegerRing, RealField,
    ComplexField, RationalField)

import gp_cremona
import sea

from gp_simon import simon_two_descent

import ell_tate_curve

factor = arith.factor
sqrt = math.sqrt
exp = math.exp
mul = misc.mul
next_prime = arith.next_prime
kronecker_symbol = arith.kronecker_symbol

Q = RationalField()
C = ComplexField()
R = RealField()
IR = rings.RealIntervalField(20)

_MAX_HEIGHT=21

# CMJ is a dict of pairs (j,D) where j is a rational CM j-invariant
# and D is the corresponding quadratic discriminant

CMJ={ 0: -3, 54000: -12, -12288000: -27, 1728: -4, 287496: -16,
      -3375: -7, 16581375: -28, 8000: -8, -32768: -11, -884736: -19,
      -884736000: -43, -147197952000: -67, -262537412640768000: -163}



class EllipticCurve_rational_field(EllipticCurve_number_field):
    """
    Elliptic curve over the Rational Field.
    """
    def __init__(self, ainvs, extra=None):
        if extra != None:   # possibility of two arguments (the first would be the field)
            ainvs = extra
        if isinstance(ainvs, str):
            label = ainvs
            X = sage.databases.cremona.CremonaDatabase()[label]
            EllipticCurve_number_field.__init__(self, Q, X.a_invariants())
            for attr in ['rank', 'torsion_order', 'cremona_label', 'conductor',
                         'modular_degree', 'gens', 'regulator']:
                s = "_EllipticCurve_rational_field__"+attr
                if hasattr(X,s):
                    setattr(self, s, getattr(X, s))
            return
        EllipticCurve_number_field.__init__(self, Q, ainvs)
        self.__np = {}
        self.__gens = {}
        self.__rank = {}
        self.__regulator = {}
        if self.base_ring() != Q:
            raise TypeError, "Base field (=%s) must be the Rational Field."%self.base_ring()

    def _set_rank(self, r):
        """
        Internal function to set the cached rank of this elliptic curve to r.

        WARNING:  No checking is done!  Not intended for use by users.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E._set_rank(99)  # bogus value -- not checked
            sage: E.rank()         # returns bogus cached value
            99
            sage: E.gens()         # causes actual rank to be computed
            [(0 : -1 : 1)]
            sage: E.rank()         # the correct rank
            1
        """
        self.__rank = {}
        self.__rank[True] = Integer(r)

    def _set_torsion_order(self, t):
        """
        Internal function to set the cached torsion order of this elliptic curve to t.

        WARNING:  No checking is done!  Not intended for use by users.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E._set_torsion_order(99)  # bogus value -- not checked
            sage: E.torsion_order()         # returns bogus cached value
            99
            sage: T = E.torsion_subgroup()  # causes actual torsion to be computed
            sage: E.torsion_order()         # the correct value
            1
        """
        self.__torsion_order = Integer(t)

    def _set_cremona_label(self, L):
        """
        Internal function to set the cached label of this elliptic curve to L.

        WARNING:  No checking is done!  Not intended for use by users.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E._set_cremona_label('bogus')
            sage: E.label()
            'bogus'
            sage: E.database_curve().label()
            '37a1'
            sage: E.label() # no change
            'bogus'
            sage: E._set_cremona_label(E.database_curve().label())
            sage: E.label() # now it is correct
            '37a1'
        """
        self.__cremona_label = L

    def _set_conductor(self, N):
        """
        Internal function to set the cached conductor of this elliptic curve to N.

        WARNING: No checking is done!  Not intended for use by users.
        Setting to the wrong value will cause strange problems (see
        examples).

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E._set_conductor(99)      # bogus value -- not checked
            sage: E.conductor()             # returns bogus cached value
            99

        This will not work since the conductor is used when searching
        the database:
            sage: E._set_conductor(E.database_curve().conductor())
            Traceback (most recent call last):
            ...
            RuntimeError: Elliptic curve ... not in the database.
            sage: E._set_conductor(EllipticCurve(E.a_invariants()).database_curve().conductor())
            sage: E.conductor()             # returns correct value
            37
        """
        self.__conductor_pari = Integer(N)

    def _set_modular_degree(self, deg):
        """
        Internal function to set the cached modular degree of this elliptic curve to deg.

        WARNING:  No checking is done!

        EXAMPLES:
            sage: E=EllipticCurve('5077a1')
            sage: E.modular_degree()
            1984
            sage: E._set_modular_degree(123456789)
            sage: E.modular_degree()
            123456789
            sage: E._set_modular_degree(E.database_curve().modular_degree())
            sage: E.modular_degree()
            1984
        """
        self.__modular_degree = Integer(deg)

    def _set_gens(self, gens):
        """
        Internal function to set the cached generators of this elliptic curve to gens.

        WARNING:  No checking is done!

        EXAMPLES:
            sage: E=EllipticCurve('5077a1')
            sage: E.rank()
            3
            sage: E.gens() # random
            [(-2 : 3 : 1), (-7/4 : 25/8 : 1), (1 : -1 : 1)]
            sage: E._set_gens([]) # bogus list
            sage: E.rank()        # unchanged
            3
            sage: E._set_gens(E.database_curve().gens())
            sage: E.gens()
            [(-2 : 3 : 1), (-7/4 : 25/8 : 1), (1 : -1 : 1)]
        """
        self.__gens = {}
        self.__gens[True] = [self.point(x, check=True) for x in gens]
        self.__gens[True].sort()

    def is_integral(self):
        """
        Returns True if this elliptic curve has integral coefficients (in Z)

        EXAMPLES:
            sage: E=EllipticCurve(QQ,[1,1]); E
            Elliptic Curve defined by y^2  = x^3 + x +1 over Rational Field
            sage: E.is_integral()
            True
            sage: E2=E.change_weierstrass_model(2,0,0,0); E2
            Elliptic Curve defined by y^2  = x^3 + 1/16*x + 1/64 over Rational Field
            sage: E2.is_integral()
            False
         """
        try:
            return self.__is_integral
        except AttributeError:
            one = Integer(1)
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

        NOTE: The output is a raw string and completely illegible
            using automatic display, so it is recommended to use print
            for legible outout

        EXAMPLES:
            sage: E = EllipticCurve('37a1')
            sage: E.mwrank() #random
            ...
            sage: print E.mwrank()
            Curve [0,0,1,-1,0] :        Basic pair: I=48, J=-432
            disc=255744
            ...
            Generator 1 is [0:-1:1]; height 0.05111...
            <BLANKLINE>
            Regulator = 0.05111...
            <BLANKLINE>
            The rank and full Mordell-Weil basis have been determined unconditionally.
            ...

        Options to mwrank can be passed:
            sage: E = EllipticCurve([0,0,0,877,0])

        Run mwrank with 'verbose' flag set to 0 but list generators if found
            sage: print E.mwrank('-v0 -l')
            Curve [0,0,0,877,0] :   0 <= rank <= 1
            Regulator = 1

        Run mwrank again, this time with a higher bound for point
        searching on homogeneous spaces:

            sage: print E.mwrank('-v0 -l -b11')
            Curve [0,0,0,877,0] :   Rank = 1
            Generator 1 is [29604565304828237474403861024284371796799791624792913256602210:-256256267988926809388776834045513089648669153204356603464786949:490078023219787588959802933995928925096061616470779979261000]; height 95.980371987964
            Regulator = 95.980371987964
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
                self.__conductor_pari = Integer(self.pari_mincurve().ellglobalred()[0])
            return self.__conductor_pari

        elif algorithm == "gp":
            try:
                return self.__conductor_gp
            except AttributeError:
                self.__conductor_gp = Integer(gp.eval('ellglobalred(ellinit(%s,0))[1]'%self.a_invariants()))
                return self.__conductor_gp

        elif algorithm == "mwrank":
            try:
                return self.__conductor_mwrank
            except AttributeError:
                if self.is_integral():
                    self.__conductor_mwrank = Integer(self.mwrank_curve().conductor())
                else:
                    self.__conductor_mwrank = Integer(self.minimal_model().mwrank_curve().conductor())
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

    def pari_curve(self, prec = None, factor = 1):
        """
        Return the PARI curve corresponding to this elliptic curve.

        INPUT:
            prec -- The precision of quantities calculated for the
                    returned curve (in decimal digits).  if None, defaults
                    to factor * the precision of the largest cached curve
                    (or 10 if none yet computed)
            factor -- the factor to increase the precision over the
                      maximum previously computed precision.  Only used if
                      prec (which gives an explicit precision) is None.

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
        if prec is None:
            try:
                L = self._pari_curve.keys()
                L.sort()
                if factor == 1:
                    return self._pari_curve[L[-1]]
                else:
                    prec = int(factor * L[-1])
                    self._pari_curve[prec] = pari(self.a_invariants()).ellinit(precision=prec)
                    return self._pari_curve[prec]
            except AttributeError:
                pass
        try:
            return self._pari_curve[prec]
        except AttributeError:
            prec = 10
            self._pari_curve = {}
        except KeyError:
            pass
        self._pari_curve[prec] = pari(self.a_invariants()).ellinit(precision=prec)
        return self._pari_curve[prec]

    def pari_mincurve(self, prec = None):
        """
        Return the PARI curve corresponding to a minimal model
        for this elliptic curve.

        INPUT:
        prec -- The precision of quantities calculated for the
                returned curve (in decimal digits).  if None, defaults
                to the precision of the largest cached curve (or 10 if
                none yet computed)

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
        if prec is None:
            try:
                L = self._pari_mincurve.keys()
                L.sort()
                return self._pari_mincurve[L[len(L) - 1]]
            except AttributeError:
                pass
        try:
            return self._pari_mincurve[prec]
        except AttributeError:
            prec = 10
            self._pari_mincurve = {}
        except KeyError:
            pass
        e = self.pari_curve(prec)
        mc, change = e.ellminimalmodel()
        self._pari_mincurve[prec] = mc
        # self.__min_transform = change
        return mc

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

    #def __pari_double_prec(self):
    #    EllipticCurve_number_field._EllipticCurve__pari_double_prec(self)
    #    try:
    #        del self._pari_mincurve
    #    except AttributeError:
    #        pass

    ####################################################################
    #  Access to mwrank
    ####################################################################
    def mwrank_curve(self, verbose=False):
        """
        Construct an mwrank_EllipticCurve from this elliptic curve

        The resulting mwrank_EllipticCurve has available methods from
        John Cremona's eclib library.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: EE=E.mwrank_curve()
            sage: EE
            y^2+ y = x^3 - x^2 - 10*x - 20
            sage: type(EE)
            <class 'sage.libs.mwrank.interface.mwrank_EllipticCurve'>
            sage: EE.isogeny_class()
            ([[0, -1, 1, -10, -20], [0, -1, 1, -7820, -263580], [0, -1, 1, 0, 0]],
            [[0, 5, 5], [5, 0, 0], [5, 0, 0]])
        """
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
                           If False, *no output* is
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

            Nothing -- nothing is returned (though much is printed
                                            unless verbose=False)

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.two_descent(verbose=False) # no output
        """
        self.mwrank_curve().two_descent(verbose, selmer_only,
                                        first_limit, second_limit,
                                        n_aux, second_descent)


    ####################################################################
    #  Etc.
    ####################################################################

    def aplist(self, n, python_ints=False):
        r"""
        The Fourier coefficients $a_p$ of the modular form attached to
        this elliptic curve, for all primes $p\leq n$.

        INPUT:
            n -- integer
            python_ints -- bool (default: False); if True return a list of
                      Python ints instead of SAGE integers.

        OUTPUT:
            -- list of integers

        EXAMPLES:
            sage: e = EllipticCurve('37a')
            sage: e.aplist(1)
            []
            sage: e.aplist(2)
            [-2]
            sage: e.aplist(10)
            [-2, -3, -2, -1]
            sage: v = e.aplist(13); v
            [-2, -3, -2, -1, -5, -2]
            sage: type(v[0])
            <type 'sage.rings.integer.Integer'>
            sage: type(e.aplist(13, python_ints=True)[0])
            <type 'int'>
        """
        # How is this result dependant on the real precision in pari?  At all?
        e = self.pari_mincurve()
        v = e.ellaplist(n, python_ints=True)
        if python_ints:
            return v
        else:
            return [Integer(a) for a in v]



    def anlist(self, n, python_ints=False):
        """
        The Fourier coefficients up to and including $a_n$ of the
        modular form attached to this elliptic curve.  The ith element
        of the return list is a[i].

        INPUT:
            n -- integer
            python_ints -- bool (default: False); if True return a list of
                      Python ints instead of SAGE integers.

        OUTPUT:
            -- list of integers

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])
            sage: E.anlist(3)
            [0, 1, -2, -1]

            sage: E = EllipticCurve([0,1])
            sage: E.anlist(20)
            [0, 1, 0, 0, 0, 0, 0, -4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 8, 0]
        """
        n = int(n)
        # How is this result dependent on the real precision in Pari?  At all?
        e = self.pari_mincurve()
        if n >= 2147483648:
            raise RuntimeError, "anlist: n (=%s) must be < 2147483648."%n

        v = [0] + e.ellan(n, python_ints=True)
        if not python_ints:
            v = [Integer(x) for x in v]
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
        r"""
        Return the $q$-expansion to precision prec of the newform
        attached to this elliptic curve.

        INPUT:
            prec -- an integer

        OUTPUT:
            a power series (in th evariable 'q')

        NOTE: If you want the output to be a modular form and not just
        a $q$-expansion, use \code{self.modular_form()}.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.q_expansion(20)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + 4*q^10 - 5*q^11 - 6*q^12 - 2*q^13 + 2*q^14 + 6*q^15 - 4*q^16 - 12*q^18 + O(q^20)

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

            If you need to see more terms in the $q$-expansion:
            sage: f.q_expansion(20)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + 4*q^10 - 5*q^11 - 6*q^12 - 2*q^13 + 2*q^14 + 6*q^15 - 4*q^16 - 12*q^18 + O(q^20)


        NOTE: If you just want the $q$-expansion, use
        \code{self.q_expansion(prec)}.
        """
        try:
            return self.__modular_form
        except AttributeError:
            M = sage.modular.modform.constructor.ModularForms(self.conductor(),weight=2)
            f = sage.modular.modform.element.ModularFormElement_elliptic_curve(M, self)
            self.__modular_form = f
            return f

    def modular_symbol_space(self, sign=1, base_ring=Q, bound=None):
        r"""
        Return the space of cuspidal modular symbols associated to
        this elliptic curve, with given sign and base ring.

        INPUT:
            sign -- 0, -1, or 1
            base_ring -- a ring

        EXAMPLES:
            sage: f = EllipticCurve('37b')
            sage: f.modular_symbol_space()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(37) of weight 2 with sign 1 over Rational Field
            sage: f.modular_symbol_space(-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(37) of weight 2 with sign -1 over Rational Field
            sage: f.modular_symbol_space(0, bound=3)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field

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
        M = ell_modular_symbols.modular_symbol_space(self, sign, base_ring, bound=bound)
        self.__modular_symbol_space[typ] = M
        return M

    def modular_symbol(self, sign=1, normalize=True):
        r"""
        Return the modular symbol associated to this elliptic curve,
        with given sign and base ring.  This is the map that sends r/s
        to a fixed multiple of 2*pi*I*f(z)dz from oo to r/s,
        normalized so that all values of this map take values in QQ.

        If sign=1, the normalization is such that the p-adic
        L-function associated to this modular symbol is correct.
        I.e., the normalization is the same as for the integral period
        mapping divided by 2.

        INPUT:
            sign -- -1, or 1
            base_ring -- a ring
            normalize -- (default: True); if True, the modular symbol
                is correctly normalized (up to possibly a factor of
                -1 or 2).  If False, the modular symbol is almost certainly
                not correctly normalized, i.e., all values will be a
                fixed scalar multiple of what they should be.  But
                the initial computation of the modular symbol is
                much faster, though evaluation of it after computing
                it won't be any faster.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: M=E.modular_symbol(); M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: M(1/2)
            0
            sage: M(1/5)
            1/2
        """
        typ = (sign, normalize)
        try:
            return self.__modular_symbol[typ]
        except AttributeError:
            self.__modular_symbol = {}
        except KeyError:
            pass
        M = ell_modular_symbols.ModularSymbol(self, sign, normalize)
        self.__modular_symbol[typ] = M
        return M

    padic_lseries = padics.padic_lseries

    def newform(self):
        r"""
        Same as \code{self.modular_form()}.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.newform()
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
            sage: E.newform() == E.modular_form()
            True
        """
        return self.modular_form()

    def q_eigenform(self, prec):
        r"""
        Synonym for \code{self.q_expansion(prec)}.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.q_eigenform(10)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + O(q^10)
            sage: E.q_eigenform(10) == E.q_expansion(10)
            True
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
            sage: E.analytic_rank(algorithm='cremona')
            2
            sage: E.analytic_rank(algorithm='rubinstein')
            2
            sage: E.analytic_rank(algorithm='sympow')
            2
            sage: E.analytic_rank(algorithm='magma')    # optional
            2
            sage: E.analytic_rank(algorithm='all')
            2
        """
        if algorithm == 'cremona':
            return rings.Integer(gp_cremona.ellanalyticrank(self.minimal_model().a_invariants()))
        elif algorithm == 'rubinstein':
            from sage.lfunctions.lcalc import lcalc
            return lcalc.analytic_rank(L=self)
        elif algorithm == 'sympow':
            from sage.lfunctions.sympow import sympow
            return sympow.analytic_rank(self)[0]
        elif algorithm == 'magma':
            return rings.Integer(self._magma_().AnalyticRank())
        elif algorithm == 'all':
            S = list(set([self.analytic_rank('cremona'),
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
        2-Selmer group, and a list of independent points on the curve.

        INPUT:
            verbose -- integer, 0,1,2,3; (default: 0), the verbosity level
            lim1    -- (default: 5) limite des points triviaux sur les quartiques
            lim3    -- (default: 50) limite des points sur les quartiques ELS
            limtriv -- (default: 10) limite des points triviaux sur la
                                     courbe elliptique
            maxprob -- (default: 20)
            limbigprime -- (default: 30)  to distinguish between small and large prime
                                          numbers. Use probabilistic tests for large
                                          primes. If 0, don't any probabilistic tests.

        OUTPUT:
            integer -- "probably" the rank of self
            integer -- the 2-rank of the Selmer group
            list    -- list of independent points on the curve.

        IMPLEMENTATION: Uses {\bf Denis Simon's} GP/PARI scripts from
                         \url{http://www.math.unicaen.fr/~simon/}

        EXAMPLES:
        These computations use pseudo-random numbers, so we set the
        seed for reproducible testing.

        We compute the ranks of the curves of lowest known conductor up to rank $8$.
        Amazingly, each of these computations finishes almost instantly!

            sage: E = EllipticCurve('11a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (0, 0, [])
            sage: E = EllipticCurve('37a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (1, 1, [(0 : 0 : 1)])
            sage: E = EllipticCurve('389a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (2, 2, [(1 : 0 : 1), (-11/9 : -55/27 : 1)])
            sage: E = EllipticCurve('5077a1')
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (3, 3, [(1 : 0 : 1), (2 : -1 : 1), (0 : 2 : 1)])


        In this example Simon's program does not find any points, though
        it does correctly compute the rank of the 2-Selmer group.
            sage: E = EllipticCurve([1, -1, 0, -751055859, -7922219731979])     # long (0.6 seconds)
            sage: set_random_seed(0)
            sage: E.simon_two_descent ()
            (1, 1, [])

        The rest of these entries were taken from Tom Womack's page
        \url{http://tom.womack.net/maths/conductors.htm}

            sage: E = EllipticCurve([1, -1, 0, -79, 289])
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (4, 4, [(6 : -5 : 1), (4 : 3 : 1), (5 : -3 : 1), (8 : -15 : 1)])
            sage: E = EllipticCurve([0, 0, 1, -79, 342])
            sage: set_random_seed(0)
            sage: E.simon_two_descent()
            (5, 5, [(5 : 8 : 1), (10 : 23 : 1), (3 : 11 : 1), (4 : -10 : 1), (0 : 18 : 1)])
            sage: E = EllipticCurve([1, 1, 0, -2582, 48720])
            sage: set_random_seed(0)
            sage: r, s, G = E.simon_two_descent(); r,s
            (6, 6)
            sage: E = EllipticCurve([0, 0, 0, -10012, 346900])
            sage: set_random_seed(0)
            sage: r, s, G = E.simon_two_descent(); r,s
            (7, 7)
            sage: E = EllipticCurve([0, 0, 1, -23737, 960366])
            sage: set_random_seed(0)
            sage: r, s, G = E.simon_two_descent(); r,s
            (8, 8)
        """
        t = simon_two_descent(self, verbose=verbose, lim1=lim1, lim3=lim3, limtriv=limtriv,
                              maxprob=maxprob, limbigprime=limbigprime)
        prob_rank = rings.Integer(t[0])
        two_selmer_rank = rings.Integer(t[1])
        prob_gens = [self(P) for P in t[2]]
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
                   algorithm='mwrank_shell',
                   proof=None):
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
            proof -- bool or None (default: None, see proof.elliptic_curve or
                       sage.structure.proof).  Note that results
                       obtained from databases are considered proof = True

        OUTPUT:
            rank (int) -- the rank of the elliptic curve.

        IMPLEMENTATION: Uses L-functions, mwrank, and databases.

        EXAMPLES:
            sage: EllipticCurve('11a').rank()
            0
            sage: EllipticCurve('37a').rank()
            1
            sage: EllipticCurve('389a').rank()
            2
            sage: EllipticCurve('5077a').rank()
            3
            sage: EllipticCurve([1, -1, 0, -79, 289]).rank()   # long time.  This will use the default proof behavior of True.
            4
            sage: EllipticCurve([0, 0, 1, -79, 342]).rank(proof=False)  # long time -- but under a minute
            5
            sage: EllipticCurve([0, 0, 1, -79, 342]).simon_two_descent()[0]  # much faster -- almost instant.
            5

        Examples with denominators in defining equations:
            sage: E = EllipticCurve( [0, 0, 0, 0, -675/4])
            sage: E.rank()
            0
            sage: E = EllipticCurve( [0, 0, 1/2, 0, -1/5])
            sage: E.rank()
            1
            sage: E.minimal_model().rank()
            1
        """
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "elliptic_curve")
        else:
            proof = bool(proof)
        try:
            return self.__rank[proof]
        except KeyError:
            if proof is False and self.__rank.has_key(True):
                return self.__rank[True]
        if use_database:
            try:
                self.__rank[True] = self.database_curve().rank()
                return self.__rank[True]
            except (AttributeError, RuntimeError):
                pass
        if not only_use_mwrank:
            N = self.conductor()
            prec = int(4*float(sqrt(N))) + 10
            if self.root_number() == 1:
                L, err = self.lseries().at1(prec)
                if abs(L) > err + R(0.0001):  # definitely doesn't vanish
                    misc.verbose("rank 0 because L(E,1)=%s"%L)
                    self.__rank[proof] = 0
                    return self.__rank[proof]
            else:
                Lprime, err = self.lseries().deriv_at1(prec)
                if abs(Lprime) > err + R(0.0001):  # definitely doesn't vanish
                    misc.verbose("rank 1 because L'(E,1)=%s"%Lprime)
                    self.__rank[proof] = 1
                    return self.__rank[proof]

        if algorithm == 'mwrank_lib':
            misc.verbose("using mwrank lib")
            C = self.mwrank_curve()
            C.set_verbose(verbose)
            r = C.rank()
            if not C.certain():
                del self.__mwrank_curve
                raise RuntimeError, "Unable to compute the rank with certainty (lower bound=%s).  This could be because Sha(E/Q)[2] is nontrivial."%C.rank() + "\nTrying calling something like two_descent(second_limit=13) on the curve then trying this command again.  You could also try rank with only_use_mwrank=False."
            self.__rank[proof] = r
        elif algorithm == 'mwrank_shell':
            misc.verbose("using mwrank shell")
            X = self.mwrank()
            if not 'The rank and full Mordell-Weil basis have been determined unconditionally' in X:
                if proof:
                    raise RuntimeError, '%s\nRank not provably correct.'%X
                else:
                    misc.verbose("Warning -- rank not provably correct", level=1)
            elif proof is False:
                proof = True #since we actually provably found the rank
            i = X.find('Rank = ')
            assert i != -1
            j = i + X[i:].find('\n')
            self.__rank[proof] = Integer(X[i+7:j])
        return self.__rank[proof]

    def gens(self, verbose=False, rank1_search=10,
             algorithm='mwrank_shell',
             only_use_mwrank=True,
             proof = None):
        """
        Compute and return generators for the Mordell-Weil group E(Q)
        *modulo* torsion.

        HINT: If you would like to control the height bounds used in
        the 2-descent, first call the two_descent function with those
        height bounds. However that function, while it displays a lot
        of output, returns no values.

        TODO: (1) Right now this function assumes that the input curve
        is in minimal Weierstrass form.  This restriction will be
        removed in the future.  This function raises a
        NotImplementedError if a non-minimal curve is given as input.

              (2) Allow passing of command-line parameters to mwrank.

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
            only_use_mwrank -- bool (default True) if false, attempts to
                       first use more naive, natively implemented methods.
            proof -- bool or None (default None, see proof.elliptic_curve or
                       sage.structure.proof).
        OUTPUT:
            generators -- List of generators for the Mordell-Weil group.

        IMPLEMENTATION: Uses Cremona's mwrank C library.

        EXAMPLES:
            sage: E = EllipticCurve('389a')
            sage: E.gens()                 # random output
            [(-1 : 1 : 1), (0 : 0 : 1)]
        """
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "elliptic_curve")
        else:
            proof = bool(proof)
        try:
            return list(self.__gens[proof])  # return copy so not changed
        except KeyError:
            if proof is False and self.__gens.has_key(True):
                return self.__gens[True]
        if self.conductor() > 10**7:
            only_use_mwrank = True

        if not only_use_mwrank:
            try:
                misc.verbose("Trying to compute rank.")
                r = self.rank(only_use_mwrank = False)
                misc.verbose("Got r = %s."%r)
                if r == 0:
                    misc.verbose("Rank = 0, so done.")
                    self.__gens[True] = []
                    self.__regulator[True] = R(1)
                    return self.__gens[True]
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
                            self.__gens[True] = G
                            self.__regulator[True] = reg
                            return self.__gens[True]
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
            if proof is True and C.certain() is False:
                del self.__mwrank_curve
                raise RuntimeError, "Unable to compute the rank, hence generators, with certainty (lower bound=%s).  This could be because Sha(E/Q)[2] is nontrivial."%C.rank() + \
                      "\nTrying calling something like two_descent(second_limit=13) on the curve then trying this command again."
            else:
                proof = C.certain()
        else:
            X = self.mwrank()
            misc.verbose("Calling mwrank shell.")
            if not 'The rank and full Mordell-Weil basis have been determined unconditionally' in X:
                msg = 'Generators not provably computed.'
                if proof:
                    raise RuntimeError, '%s\n%s'%(X,msg)
                else:
                    misc.verbose("Warning -- %s"%msg, level=1)
            elif proof is False:
                proof = True
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
            self.__regulator[proof] = R(X[i+len('Regulator = '):j])
        ####
        self.__gens[proof] = [self.point(x, check=True) for x in G]
        self.__gens[proof].sort()
        self.__rank[proof] = len(self.__gens[proof])
        return self.__gens[proof]

    def gens_certain(self):
        """
        Return True if the generators have been proven correct.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.gens()
            [(0 : -1 : 1)]
            sage: E.gens_certain()
            True
       """
        return self.__gens.has_key(True)

    def ngens(self, proof = None):
        """
        Return the number of generators of this elliptic curve.

        NOTE: See gens() for further documentation.  The function
        ngens() calls gens() if not already done, but only with
        default parameters.  Better results may be obtained by calling
        mwrank() with carefully chosen parameters.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.ngens()
            1

        TO DO: This example should not cause a run-time error.
            sage: E=EllipticCurve([0,0,0,877,0])
            sage: # E.ngens()  ######## causes run-time error

            sage: print E.mwrank('-v0 -b12 -l')
            Curve [0,0,0,877,0] :   Rank = 1
            Generator 1 is [29604565304828237474403861024284371796799791624792913256602210:-256256267988926809388776834045513089648669153204356603464786949:490078023219787588959802933995928925096061616470779979261000]; height 95.980371987964
            Regulator = 95.980...


        """
        return len(self.gens(proof = proof))

    def regulator(self, use_database=True, verbose=None, proof=None):
        """
        Returns the regulator of this curve, which must be defined
        over Q.

        INPUT:
            use_database -- bool (default: False), if True, try to
                  look up the regulator in the Cremona database.
            verbose -- (default: None), if specified changes the
                  verbosity of mwrank computations.
            proof -- bool or None (default: None, see proof.[tab] or
                       sage.structure.proof).  Note that results from
                       databases are considered proof = True

        EXAMPLES:
            sage: E = EllipticCurve([0, 0, 1, -1, 0])
            sage: E.regulator()              # long time (1 second)
            0.0511114082399688
            sage: EllipticCurve('11a').regulator()
            1.00000000000000
            sage: EllipticCurve('37a').regulator()
            0.0511114082399688
            sage: EllipticCurve('389a').regulator()
            0.152460177943144
            sage: EllipticCurve('5077a').regulator()    # random low order bit
            0.417143558758385
            sage: EllipticCurve([1, -1, 0, -79, 289]).regulator()  # long time (seconds)
            1.50434488827528
            sage: EllipticCurve([0, 0, 1, -79, 342]).regulator(proof=False)  # long time (seconds)
            14.7905275701310
        """
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "elliptic_curve")
        else:
            proof = bool(proof)
        try:
            return self.__regulator[proof]
        except KeyError:
            if proof is False and self.__regulator.has_key(True):
                return self.__regulator[True]
        if use_database:
            try:
                self.__regulator[True] = R(self.database_curve().db_extra[3])
                return self.__regulator[True]
            except (AttributeError, RuntimeError):
                pass
        G = self.gens(proof=proof)
        try:  # in some cases self.gens() efficiently computes regulator.
            return self.__regulator[proof]
        except KeyError:
            if proof is False and self.__regulator.has_key(True):
                return self.__regulator[True]
        C = self.mwrank_curve()
        reg = R(C.regulator())
        if proof is True and not C.certain():
            raise RuntimeError, "Unable to compute the rank, hence regulator, with certainty (lower bound=%s)."%C.rank()
        proof = C.certain()
        self.__regulator[proof] = reg
        return self.__regulator[proof]

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
            index (int) -- the index of the group generated by points
                           in their saturation
            regulator (float) -- regulator of saturated points.

        IMPLEMENTATION:
            Uses Cremona's mwrank package.  With max_prime=0, we call
            mwrank with successively larger prime bounds until the
            full saturation is provably found.  The results of
            saturation at the previous primes is stored in each case,
            so this should be reasonably fast.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: P=E.gens()[0]
            sage: Q=5*P; Q
            (1/4 : -3/8 : 1)
            sage: E.saturation([Q])
            ([(0 : -1 : 1)], '5', 0.0511114075779915)
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
        r""" Return the Silverman height bound.  This is a positive real
        (floating point) number B such that for all rational points
        $P$ on the curve,
        $$
             h(P) \le \hat{h}(P) + B
        $$
        where h(P) is the logarithmic height of $P$ and $\hat{h}(P)$
        is the canonical height.

        Note that the CPS_height_bound is often better (i.e. smaller)
        than the Silverman bound.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.silverman_height_bound()
            4.8254007581809182
            sage: E.CPS_height_bound()
            0.16397076103046915
        """
        return self.mwrank_curve().silverman_bound()


    def point_search(self, height_limit, verbose=True):
        """
        Search for points on a curve up to an input bound on the naive
        logarithmic height.

        INPUT:

            height_limit (float) -- bound on naive height (at most 21,
                                    or mwrank overflows -- see below)

            verbose (bool) -- (default: True)

                              If True, report on each point as found
                              together with linear relations between
                              the points found and the saturation
                              process.

                              If False, just return the result.

        OUTPUT: points (list) -- list of independent points which
            generate the subgroup of the Mordell-Weil group generated
            by the points found and then p-saturated for p<20.

        WARNING:
            height_limit is logarithmic, so increasing by 1 will cause
            the running time to increase by a factor of approximately
            4.5 (=exp(1.5)).  The limit of 21 is to prevent overflow,
            but in any case using height_limit=20 takes rather a long
            time!

        IMPLEMENTATION: Uses Cremona's mwrank package. At the heart of
        this function is Cremona's port of Stoll's ratpoints program
        (version 1.4).

        EXAMPLES:
            sage: E=EllipticCurve('389a1')
            sage: E.point_search(5, verbose=False)
            [(0 : -1 : 1), (-1 : 1 : 1)]

        Increasing the height_limit takes longer, but finds no more points:
            sage: E.point_search(10, verbose=False)
            [(0 : -1 : 1), (-1 : 1 : 1)]

        In fact this curve has rank 2 so no more than 2 points will
        ever be output, but we are not using this fact.

            sage: E.saturation(_)
            ([(0 : -1 : 1), (-1 : 1 : 1)], '1', 0.152460172772408)

        What this shows is that if the rank is 2 then the points
        listed do generate the Mordell-Weil group (mod torsion).
        Finally,

            sage: E.rank()
            2

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

        This will be 0, 1 or 2.

        NOTE: As s side-effect of calling this function, the full
        torsion subgroup of the curve is computed (if not already
        cached).  A simpler implementation of this function would be
        possible (by counting the roots of the 2-division polynomial),
        but the full torsion subgroup computation is not expensive.

        EXAMPLES:
            sage: EllipticCurve('11a1').two_torsion_rank()
            0
            sage: EllipticCurve('14a1').two_torsion_rank()
            1
            sage: EllipticCurve('15a1').two_torsion_rank()
            2
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

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: [E.an(n) for n in range(20) if n>0]
            [1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0]
        """
        return Integer(self.pari_mincurve().ellak(n))

    def ap(self, p):
        """
        The p-th Fourier coefficient of the modular form corresponding
        to this elliptic curve, where p is prime.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: [E.ap(p) for p in prime_range(50)]
            [-2, -3, -2, -1, -5, -2, 0, 0, 2, 6, -4, -1, -9, 2, -9]

        """
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        return Integer(self.pari_mincurve().ellap(p))

    def quadratic_twist(self, D):
        """
        Return the quadratic twist of this elliptic curve by D.

        D must be a nonzero rational number.

        NOTE: This function overrides the generic quadratic_twist()
            function for elliptic curves, returning a minimal model

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E2=E.quadratic_twist(2); E2
            Elliptic Curve defined by y^2  = x^3 - 4*x + 2 over Rational Field
            sage: E2.conductor()
            2368
            sage: E2.quadratic_twist(2) == E
            True
       """
        return EllipticCurve_number_field.quadratic_twist(self, D).minimal_model()

    def minimal_model(self):
        r"""
        Return the unique minimal Weierstrass equation for this
        elliptic curve.  This is the model with minimal discriminant
        and $a_1,a_2,a_3 \in \{0,\pm 1\}$.

        EXAMPLES:
            sage: E=EllipticCurve([10,100,1000,10000,1000000])
            sage: E.minimal_model()
            Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + x +1 over Rational Field
        """
        try:
            return self.__minimal_model
        except AttributeError:
            F = self.pari_mincurve()
            self.__minimal_model = EllipticCurve_rational_field([Q(F[i]) for i in range(5)])
            return self.__minimal_model

    def is_minimal(self):
        r"""
        Return True iff this elliptic curve is a reduced minimal model.

        the unique minimal Weierstrass equation for this
        elliptic curve.  This is the model with minimal discriminant
        and $a_1,a_2,a_3 \in \{0,\pm 1\}$.

        TO DO: This is not very efficient since it just computes the
        minimal model and compares.  A better implementation using the
        Kraus conditions would be preferable.

        EXAMPLES:
            sage: E=EllipticCurve([10,100,1000,10000,1000000])
            sage: E.is_minimal()
            False
            sage: E=E.minimal_model()
            sage: E.is_minimal()
            True
        """
        return self.ainvs() == self.minimal_model().ainvs()

#    Duplicate!
#
#    def is_integral(self):
#        for n in self.ainvs():
#            if n.denominator() != 1:
#                return False
#        return True

    def kodaira_type(self, p):
        """
        Local Kodaira type of the elliptic curve at $p$.

        INPUT:
           -- p, an integral prime
        OUTPUT:
           -- the kodaira type of this elliptic curve at p, as a KodairaSymbol.

        EXAMPLES:
            sage: E = EllipticCurve('124a')
            sage: E.kodaira_type(2)
            IV
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
            from kodaira_symbol import KodairaSymbol
            self.__kodaira_type[p] = KodairaSymbol(v[1])
            self.__tamagawa_number[p] = Integer(v[3])
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
        """
        Return a list of all Tamagawa numbers for all prime divisors of
        the conductor (in order).

        EXAMPLES:
            sage: e = EllipticCurve('30a1')
            sage: e.tamagawa_numbers()
            [2, 3, 1]
            sage: vector(e.tamagawa_numbers())
            (2, 3, 1)
        """
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
        invs = self.short_weierstrass_model().ainvs()
        x = rings.polygen(self.base_ring())
        f = x**3 + invs[3]*x + invs[4]
        if f.discriminant() > 0:
            return 2
        else:
            return 1

    def period_lattice(self):
        r"""
        Returns the period lattice of the elliptic curve.

        EXAMPLES:
        sage: E = EllipticCurve('37a')
        sage: E.period_lattice()
        Period lattice associated to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        from sage.schemes.elliptic_curves.period_lattice import PeriodLattice_ell
        return PeriodLattice_ell(self)

    def lseries(self):
        """
        Returns the L-series of this elliptic curve.

        Further documentation is available for the functions which
        apply to the L-series.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.lseries()
            Complex L-series of the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        try:
            return self.__lseries
        except AttributeError:
            from lseries_ell import Lseries_ell
            self.__lseries = Lseries_ell(self)
            return self.__lseries

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
            -0.354172680517... + 0.874518681720...*I
        """
        s = C(s)
        N = self.conductor()
        pi = R(constants.pi)
        Gamma = transcendental.gamma
        Gamma_inc = transcendental.gamma_inc
        a = self.anlist(prec)
        eps = self.root_number()
        sqrtN = float(N.sqrt())
        def _F(n, t):
            return Gamma_inc(t+1, 2*pi*n/sqrtN) * C(sqrtN/(2*pi*n))**(t+1)
        return sum([a[n]*(_F(n,s-1) + eps*_F(n,1-s)) for n in xrange(1,prec+1)])

    def is_local_integral_model(self,*p):
        r"""
        Tests if self is integral at the prime $p$, or at all the
        primes if $p$ is a list or tuple of primes

        EXAMPLES:
            sage: E=EllipticCurve([1/2,1/5,1/5,1/5,1/5])
            sage: [E.is_local_integral_model(p) for p in (2,3,5)]
            [False, True, False]
            sage: E.is_local_integral_model(2,3,5)
            False
            sage: Eint2=E.local_integral_model(2)
            sage: Eint2.is_local_integral_model(2)
            True
        """
        if len(p)==1: p=p[0]
        if isinstance(p,(tuple,list)):
            return misc.forall(p, lambda x : self.is_local_integral_model(x))[0]
        assert p.is_prime(), "p must be prime in is_local_integral_model()"
        return misc.forall(self.ainvs(), lambda x : x.valuation(p) >= 0)[0]

    def local_integral_model(self,p):
        r"""
        Return a model of self which is integral at the prime $p$

        EXAMPLES:
            sage: E=EllipticCurve([0, 0, 1/216, -7/1296, 1/7776])
            sage: E.local_integral_model(2)
	     Elliptic Curve defined by y^2 + 1/27*y = x^3 - 7/81*x + 2/243 over Rational Field
            sage: E.local_integral_model(3)
             Elliptic Curve defined by y^2 + 1/8*y = x^3 - 7/16*x + 3/32 over Rational Field
            sage: E.local_integral_model(2).local_integral_model(3) == EllipticCurve('5077a1')
            True
        """
        assert p.is_prime(), "p must be prime in local_integral_model()"
        ai = self.a_invariants()
        e  = min([(ai[i].valuation(p)/[1,2,3,4,6][i]) for i in range(5)]).floor()
        return constructor.EllipticCurve([ai[i]/p**(e*[1,2,3,4,6][i]) for i in range(5)])

    def is_global_integral_model(self):
        r"""
        Return true iff self is integral at all primes

        EXAMPLES:
            sage: E=EllipticCurve([1/2,1/5,1/5,1/5,1/5])
            sage: E.is_global_integral_model()
            False
            sage: Emin=E.global_integral_model()
            sage: Emin.is_global_integral_model()
            True
        """
        return self.is_integral()

    def global_integral_model(self):
        r"""
        Return a model of self which is integral at all primes

        EXAMPLES:
            sage: E = EllipticCurve([0, 0, 1/216, -7/1296, 1/7776])
            sage: F = E.global_integral_model(); F
            Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
            sage: F == EllipticCurve('5077a1')
            True
        """
        ai = self.a_invariants()
        for a in ai:
            if not a.is_integral():
               for p, _ in a.denom().factor():
                  e  = min([(ai[i].valuation(p)/[1,2,3,4,6][i]) for i in range(5)]).floor()
                  ai = [ai[i]/p**(e*[1,2,3,4,6][i]) for i in range(5)]
        for z in ai:
            assert z.denominator() == 1, "bug in global_integral_model: %s" % ai
        return constructor.EllipticCurve(ai)

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

    def integral_weierstrass_model(self):
        r"""
        Return a model of the form $y^2 = x^3 + a*x + b$ for this curve with $a,b\in\Z$.

        EXAMPLES:
            sage: E = EllipticCurve('17a1')
            sage: E.integral_weierstrass_model()
            Elliptic Curve defined by y^2  = x^3 - 11*x - 890 over Rational Field
        """
        F = self.minimal_model().short_weierstrass_model()
        _,_,_,A,B = F.ainvs()
        for p in [2,3]:
            e=min(A.valuation(p)/4,B.valuation(p)/6).floor()
            A /= Integer(p**(4*e))
            B /= Integer(p**(6*e))
        return constructor.EllipticCurve([A,B])

    def modular_degree(self, algorithm='sympow'):
        r"""
        Return the modular degree of this elliptic curve.

        The result is cached.  Subsequence calls, even with a
        different algorithm, just returned the cached result.

        INPUT:
           algorithm -- string:
              'sympow' -- (default) use Mark Watkin's (newer) C program sympow
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

        We compute the modular degree of the curve with rank 4 having smallest
        (known) conductor:

            sage: E = EllipticCurve([1, -1, 0, -79, 289])
            sage: factor(E.conductor())  # conductor is 234446
            2 * 117223
            sage: factor(E.modular_degree())
            2^7 * 2617
        """
        try:
            return self.__modular_degree

        except AttributeError:
            if algorithm == 'sympow':
                from sage.lfunctions.all import sympow
                m = sympow.modular_degree(self)
            elif algorithm == 'magma':
                m = rings.Integer(self._magma_().ModularDegree())
            else:
                raise ValueError, "unknown algorithm %s"%algorithm
            self.__modular_degree = m
            return m

    def modular_parametrization(self):
        r"""
        Computes and returns the modular parametrization of this elliptic curve.

        The curve is converted to a minimal model.

        OUTPUT: A list of two Larent series [X(x),Y(x)] of degrees -2,
        -3 respectively, which satisfy the equation of the (minimal
        model of the) elliptic curve.  The are modular functions on
        $\Gamma_0(N)$ where $N$ is the conductor.

        X.deriv()/(2*Y+a1*X+a3) should equal f(q)dq/q where f is
        self.q_expansion().

        EXAMPLES:
            sage: E=EllipticCurve('389a1')
            sage: X,Y=E.modular_parametrization()
            sage: X
            q^-2 + 2*q^-1 + 4 + 7*q + 13*q^2 + 18*q^3 + 31*q^4 + 49*q^5 + 74*q^6 + 111*q^7 + 173*q^8 + 251*q^9 + 379*q^10 + 560*q^11 + 824*q^12 + 1199*q^13 + 1773*q^14 + 2365*q^15 + 3463*q^16 + 4508*q^17 + O(q^18)
            sage: Y
            -q^-3 - 3*q^-2 - 8*q^-1 - 17 - 33*q - 61*q^2 - 110*q^3 - 186*q^4 - 320*q^5 - 528*q^6 - 861*q^7 - 1383*q^8 - 2218*q^9 - 3472*q^10 - 5451*q^11 - 8447*q^12 - 13020*q^13 - 20083*q^14 - 29512*q^15 - 39682*q^16 + O(q^17)

            The following should give 0, but only approximately:
            sage: q = X.parent().gen()
            sage: E.defining_polynomial()(X,Y,1) + O(q^11) == 0
            True

            Note that below we have to change variable from x to q
            sage: a1,_,a3,_,_=E.a_invariants()
            sage: f=E.q_expansion(17)
            sage: q=f.parent().gen()
            sage: f/q == (X.derivative()/(2*Y+a1*X+a3))
            True

        """
        R = LaurentSeriesRing(RationalField(),'q')
        XY = self.pari_mincurve().elltaniyama()
        return [1/R(1/XY[0]),1/R(1/XY[1])]

    def cremona_label(self, space=False):
        """
        Return the Cremona label associated to (the minimal model) of
        this curve, if it is known.  If not, raise a RuntimeError
        exception.

        EXAMPLES:
            sage: E=EllipticCurve('389a1')
            sage: E.cremona_label()
            '389a1'

        The default database only contains conductors up to 10000, so
        any curve with conductor greater than that will cause an error
        to be raised:

            sage: E=EllipticCurve([1,2,3,4,5]);
            sage: E.conductor()
            10351
            sage: E.cremona_label()
            Traceback (most recent call last):
            ...
            RuntimeError: Cremona label not known for Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field.
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

    label = cremona_label

    def torsion_order(self):
        """
        Return the order of the torsion subgroup.

        EXAMPLES:
            sage: e = EllipticCurve('11a')
            sage: e.torsion_order()
            5
            sage: type(e.torsion_order())
            <type 'sage.rings.integer.Integer'>
            sage: e = EllipticCurve([1,2,3,4,5])
            sage: e.torsion_order()
            1
            sage: type(e.torsion_order())
            <type 'sage.rings.integer.Integer'>
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

            sage: e = EllipticCurve([-1386747,368636886]);e
            Elliptic Curve defined by y^2  = x^3 - 1386747*x + 368636886 over Rational Field
            sage: G = e.torsion_subgroup(); G
            Torsion Subgroup isomorphic to Multiplicative Abelian
            Group isomorphic to C8 x C2 associated to the Elliptic
            Curve defined by y^2 = x^3 - 1386747*x + 368636886 over
            Rational Field
            sage: G.0
            (1227 : 22680 : 1)
            sage: G.1
            (282 : 0 : 1)
            sage: list(G)
            [1, P1, P0, P0*P1, P0^2, P0^2*P1, P0^3, P0^3*P1, P0^4, P0^4*P1, P0^5, P0^5*P1, P0^6, P0^6*P1, P0^7, P0^7*P1]
        """
        try:
            return self.__torsion_subgroup
        except AttributeError:
            self.__torsion_subgroup = rational_torsion.EllipticCurveTorsionSubgroup(self, flag)
            self.__torsion_order = self.__torsion_subgroup.order()
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

        EXAMPLES:
            sage: EllipticCurve('11a1').root_number()
            1
            sage: EllipticCurve('37a1').root_number()
            -1
            sage: EllipticCurve('389a1').root_number()
            1
        """
        try:
            return self.__root_number
        except AttributeError:
            self.__root_number = int(self.pari_mincurve().ellrootno())
        return self.__root_number


    def has_cm(self):
        """
        Returns True iff this elliptic curve has Complex Multiplication

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.has_cm()
            False
            sage: E=EllipticCurve('32a1')
            sage: E.has_cm()
            True
            sage: E.j_invariant()
            1728
        """

        return CMJ.has_key(self.j_invariant())

    def cm_discriminant(self):
        """
        Returns the associated quadratic discriminant if this elliptic
        curve has Complex Multiplication.

        A ValueError is raised if the curve does not have CM (see the
        function has_cm())

        EXAMPLES:
            sage: E=EllipticCurve('32a1')
            sage: E.cm_discriminant()
            -4
            sage: E=EllipticCurve('121b1')
            sage: E.cm_discriminant()
            -11
            sage: E=EllipticCurve('37a1')
           sage: E.cm_discriminant()
           Traceback (most recent call last):
           ...
           ValueError: Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field does not have CM
        """

        try:
            return CMJ[self.j_invariant()]
        except KeyError:
            raise ValueError, "%s does not have CM"%self


    def quadratic_twist(self, D):
        """
        Return the global minimal model of the quadratic twist of this
        curve by D.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E7=E.quadratic_twist(7); E7
            Elliptic Curve defined by y^2  = x^3 - 784*x + 5488 over Rational Field
            sage: E7.conductor()
            29008
            sage: E7.quadratic_twist(7) == E
            True
        """
        return EllipticCurve_number_field.quadratic_twist(self, D).minimal_model()


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

    def isogeny_graph(self):
        r"""
        Returns a graph representing the isogeny class of this elliptic curve,
        where the vertices are isogenous curves over $\Q$ and the edges are
        prime degree isogenies labeled by their degree.

        EXAMPLES:
            sage: LL = []
            sage: for e in cremona_optimal_curves(range(1, 38)):
            ...    G = e.isogeny_graph()
            ...    already = False
            ...    for H in LL:
            ...        if G.is_isomorphic(H):
            ...            already = True
            ...            break
            ...    if not already:
            ...        LL.append(G)
            ...
            sage: graphs_list.show_graphs(LL)

            sage: E = EllipticCurve('195a')
            sage: G = E.isogeny_graph()
            sage: for v in G: print v, G.get_vertex(v)
            ...
            0 Elliptic Curve defined by y^2 + x*y  = x^3 - 110*x + 435 over Rational Field
            1 Elliptic Curve defined by y^2 + x*y  = x^3 - 115*x + 392 over Rational Field
            2 Elliptic Curve defined by y^2 + x*y  = x^3 + 210*x + 2277 over Rational Field
            3 Elliptic Curve defined by y^2 + x*y  = x^3 - 520*x - 4225 over Rational Field
            4 Elliptic Curve defined by y^2 + x*y  = x^3 + 605*x - 19750 over Rational Field
            5 Elliptic Curve defined by y^2 + x*y  = x^3 - 8125*x - 282568 over Rational Field
            6 Elliptic Curve defined by y^2 + x*y  = x^3 - 7930*x - 296725 over Rational Field
            7 Elliptic Curve defined by y^2 + x*y  = x^3 - 130000*x - 18051943 over Rational Field
            sage: G.plot(edge_labels=True).save('isogeny_graph.png')
        """
        from sage.graphs.graph import Graph
        L, M = self.isogeny_class()
        G = Graph(M, format='weighted_adjacency_matrix')
        d = {}
        for v in G.vertices():
            d[v] = L[v]
        G.set_vertices(d)
        return G

    ##########################################################
    # Galois Representations
    ##########################################################

    def is_reducible(self, p):
        """
        Return True if the mod-p representation attached
        to E is reducible.

        INPUT:
            p -- a prime number

        NOTE: The answer is cached.

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
        try:
            return self.__is_reducible[p]
        except AttributeError:
            self.__is_reducible = {}
        except KeyError:
            pass

        if not arith.is_prime(p):
            raise ValueError, 'p (=%s) must be prime'%p
        # we do is_surjective first, since this is
        # much easier than computing isogeny_class
        t, why = self.is_surjective(p)
        if t == True:
            self.__is_reducible[p] = False
            return False  # definitely not reducible
        isogeny_matrix = self.isogeny_class()[ 1 ]
        v = isogeny_matrix.row(0) # first row
        for a in v:
            if a != 0 and a % p == 0:
                self.__is_reducible[p] = True
                return True
        self.__is_reducible[p] = False
        return False

    def is_irreducible(self, p):
        """
        Return True if the mod p represenation is irreducible.

        EXAMPLES:
            sage: e = EllipticCurve('37b')
            sage: e.is_irreducible(2)
            True
            sage: e.is_irreducible(3)
            False
            sage: e.is_reducible(2)
            False
            sage: e.is_reducible(3)
            True
        """
        return not self.is_reducible(p)

    def is_surjective(self, p, A=1000):
        """
        Return True if the mod-p representation attached to E
        is surjective, False if it is not, or None if we were
        unable to determine whether it is or not.

        NOTE: The answer is cached.

        INPUT:
            p -- int (a prime number)
            A -- int (a bound on the number of a_p to use)

        OUTPUT:
            a 2-tuple:
            -- surjective or (probably) not
            -- information about what it is if not surjective

        EXAMPLES:
            sage: e = EllipticCurve('37b')
            sage: e.is_surjective(2)
            (True, None)
            sage: e.is_surjective(3)
            (False, '3-torsion')


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
        A = int(A)
        key = (p, A)
        try:
            return self.__is_surjective[key]
        except KeyError:
            pass
        except AttributeError:
            self.__is_surjective = {}

        ans = self._is_surjective(p, A)
        self.__is_surjective[key] = ans
        return ans

    def _is_surjective(self, p, A):
        T = self.torsion_subgroup().order()
        if T % p == 0:
            return False, "%s-torsion"%p

        if p == 2:
            # E is isomorphic to  [0,b2,0,8*b4,16*b6]
            b2,b4,b6,b8=self.b_invariants()
            R = rings.PolynomialRing(self.base_ring(), 'x')
            x = R.gen()
            f = x**3 + b2*x**2 + 8*b4*x + 16*b6
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
        """
        Return True iff this elliptic curve is semi-stable at all primes.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.is_semistable()
            True
            sage: E=EllipticCurve('90a1')
            sage: E.is_semistable()
            False
        """
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

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.is_ordinary(37)
            True
            sage: E=EllipticCurve('32a1')
            sage: E.is_ordinary(2)
            False
            sage: [p for p in prime_range(50) if E.is_ordinary(p)]
            [5, 13, 17, 29, 37, 41]
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
        Return True precisely when p is a prime of good reduction and
        the mod-p representation attached to this elliptic curve is
        supersingular at ell.

        INPUT:
            p -- a prime
            ell - a prime (default: p)

        OUTPUT:
            bool

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.is_supersingular(37)
            False
            sage: E=EllipticCurve('32a1')
            sage: E.is_supersingular(2)
            False
            sage: E.is_supersingular(7)
            True
            sage: [p for p in prime_range(50) if E.is_supersingular(p)]
            [3, 7, 11, 19, 23, 31, 43, 47]
        """
        if ell is None:
            ell = p
        return self.is_good(p) and not self.is_ordinary(p, ell)

    def supersingular_primes(self, B):
        """
        Return a list of all supersingular primes for this elliptic curve
        up to and possibly including B.

        EXAMPLES:
            sage: e = EllipticCurve('11a')
            sage: e.aplist(20)
            [-2, -1, 1, -2, 1, 4, -2, 0]
            sage: e.supersingular_primes(1000)
            [2, 19, 29, 199, 569, 809]

            sage: e = EllipticCurve('27a')
            sage: e.aplist(20)
            [0, 0, 0, -1, 0, 5, 0, -7]
            sage: e.supersingular_primes(97)
            [2, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89]
            sage: e.ordinary_primes(97)
            [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97]
            sage: e.supersingular_primes(3)
            [2]
            sage: e.supersingular_primes(2)
            [2]
            sage: e.supersingular_primes(1)
            []
        """
        v = self.aplist(max(B, 3))
        P = arith.prime_range(max(B,3)+1)
        N = self.conductor()
        return [P[i] for i in [0,1] if P[i] <= B and v[i]%P[i]==0 and N%P[i] != 0] + \
                      [P[i] for i in range(2,len(v)) if v[i] == 0 and N%P[i] != 0]

    def ordinary_primes(self, B):
        """
        Return a list of all ordinary primes for this elliptic curve
        up to and possibly including B.

        EXAMPLES:
            sage: e = EllipticCurve('11a')
            sage: e.aplist(20)
            [-2, -1, 1, -2, 1, 4, -2, 0]
            sage: e.ordinary_primes(97)
            [3, 5, 7, 11, 13, 17, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
            sage: e = EllipticCurve('49a')
            sage: e.aplist(20)
            [1, 0, 0, 0, 4, 0, 0, 0]
            sage: e.supersingular_primes(97)
            [3, 5, 13, 17, 19, 31, 41, 47, 59, 61, 73, 83, 89, 97]
            sage: e.ordinary_primes(97)
            [2, 11, 23, 29, 37, 43, 53, 67, 71, 79]
            sage: e.ordinary_primes(3)
            [2]
            sage: e.ordinary_primes(2)
            [2]
            sage: e.ordinary_primes(1)
            []
        """
        v = self.aplist(max(B, 3) )
        P = arith.prime_range(max(B,3) +1)
        return [P[i] for i in [0,1] if P[i] <= B and v[i]%P[i]!=0] +\
               [P[i] for i in range(2,len(v)) if v[i] != 0]

    def eval_modular_form(self, points, prec):
        """
        Evaluate the L-series of this elliptic curve at points in CC

        INPUT:
              points--  a list of points in the half-plane of convergence
              prec--    precision

        OUTPUT:
              A list of values L(E,s) for s in points

        NOTE:  Better examples are welcome.

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: E.eval_modular_form([1.5+I,2.0+I,2.5+I],0.000001)
            [0, 0, 0]
        """
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

    def _multiple_of_degree_of_isogeny_to_optimal_curve(self):
        """
        Internal function returning an integer m such that the degree of
        the isogeny between this curve and the optimal curve in its
        isogeny class is a divisor of m.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: E._multiple_of_degree_of_isogeny_to_optimal_curve()
            25
            sage: E=EllipticCurve('11a2')
            sage: E._multiple_of_degree_of_isogeny_to_optimal_curve()
            5
            sage: E=EllipticCurve('11a3')
            sage: E._multiple_of_degree_of_isogeny_to_optimal_curve()
            5
        """
        M = self.isogeny_class()[1]
        return Integer(misc.prod([x for x in M.row(0) if x], 1))

    ########################################################################
    # Functions related to bounding the order of Sha (provably correctly!)
    # Heegner points and Kolyvagin's theorem
    ########################################################################

    def sha(self):
        """
        Return an object of class
        'sage.schemes.elliptic_curves.sha.Sha' attached to this
        elliptic curve.

        This can be used in functions related to bounding the order of
        Sha (The Tate-Shafarevich group of the curve).

        EXAMPLES:
            sage: E=EllipticCurve('37a1')
            sage: S=E.sha()
            sage: S
            <class 'sage.schemes.elliptic_curves.sha.Sha'>
            sage: S.bound_kolyvagin()
            ([2], 1)
        """
        try:
            return self.__sha
        except AttributeError:
            from sha import Sha
            self.__sha = Sha(self)
            return self.__sha

    def satisfies_heegner_hypothesis(self, D):
        """
        Returns True precisely when D is a fundamental discriminant
        that satisfies the Heegner hypothesis for this elliptic curve.

        EXAMPLES:
            sage: E = EllipticCurve('11a1')
            sage: E.satisfies_heegner_hypothesis(-7)
            True
            sage: E.satisfies_heegner_hypothesis(-11)
            False
        """
        if not number_field.is_fundamental_discriminant(D):
            return False
        if arith.GCD(D, self.conductor()) != 1:
            return False
        for p, _ in factor(self.conductor()):
            if kronecker_symbol(D,p) != 1:
                return False
        return True

    def heegner_discriminants(self, bound):
        """
        Return the list of self's Heegner discriminants between -1 and -bound.

        INPUT:
            bound (int) -- upper bound for -discriminant

        OUTPUT:
            The list of Heegner discriminants between -1 and -bound for the given elliptic curve.

        EXAMPLE:
            sage: E=EllipticCurve('11a')
            sage: E.heegner_discriminants(30)
            [-7, -8, -19, -24]
        """
        return [-D for D in xrange(1,bound) if self.satisfies_heegner_hypothesis(-D)]

    def heegner_discriminants_list(self, n):
        """
        Return the list of self's first n Heegner discriminants smaller than -5.

        INPUT:
            n (int) -- the number of discriminants to compute

        OUTPUT:
            The list of the first n Heegner discriminants smaller than -5 for the given elliptic curve.

        EXAMPLE:
            sage: E=EllipticCurve('11a')
            sage: E.heegner_discriminants_list(4)
            [-7, -8, -19, -24]
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
        L1_vanishes = self.lseries().L1_vanishes()
        if eps == 1 and L1_vanishes:
            return IR(0) # rank even hence >= 2, so Heegner point is torsion.
        alpha = R(sqrt(abs(D)))/(2*self.period_lattice().complex_area())
        F = self.quadratic_twist(D)
        E = self
        k_E = prec*sqrt(E.conductor()) + 20
        k_F = prec*sqrt(F.conductor()) + 20

        MIN_ERR = R('1e-6')  # we assume that regulator and
                             # discriminant, etc., computed to this accuracy (which is easily the case).
                             # this should be made more intelligent / rigorous relative
                             # to the rest of the system.
        if eps == 1:   # E has even rank
            LF1, err_F = F.lseries().deriv_at1(k_F)
            LE1, err_E = E.lseries().at1(k_E)
            err_F = max(err_F, MIN_ERR)
            err_E = max(err_E, MIN_ERR)
            return IR(alpha-MIN_ERR,alpha+MIN_ERR) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F)

        else:          # E has odd rank
            LE1, err_E = E.lseries().deriv_at1(k_E)
            LF1, err_F = F.lseries().at1(k_F)
            err_F = max(err_F, MIN_ERR)
            err_E = max(err_E, MIN_ERR)
            return IR(alpha-MIN_ERR,alpha+MIN_ERR) * IR(LE1-err_E,LE1+err_E) * IR(LF1-err_F,LF1+err_F)


    def heegner_index(self, D,  min_p=2, prec=5, verbose=False):
        r"""
        Return an interval that contains the index of the Heegner
        point $y_K$ in the group of K-rational points modulo torsion
        on this elliptic curve, computed using the Gross-Zagier
        formula and/or a point search, or the index divided by $2$.

        NOTES: If \code{min_p} is bigger than 2 then the index can be
        off by any prime less than \code{min_p}.   This function
        returns the index divided by $2$ exactly when $E(\Q)_{/tor}$
        has index $2$ in $E(K)_{/tor}$.

        INPUT:
            D (int) -- Heegner discriminant
            min_p (int) -- (default: 2) only rule out primes >= min_p
                           dividing the index.
            verbose (bool) -- (default: False); print lots of mwrank
                              search status information when computing
                              regulator
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
            [0.99999332 .. 1.0000077]

            sage: E = EllipticCurve('37b')
            sage: E.heegner_discriminants(100)
            [-3, -4, -7, -11, -40, -47, -67, -71, -83, -84, -95]
            sage: E.heegner_index(-95)          # long time (1 second)
            [1.9999923 .. 2.0000077]

        Current discriminants -3 and -4 are not supported:
            sage: E.heegner_index(-3)
            Traceback (most recent call last):
            ...
            ArithmeticError: Discriminant (=-3) must not be -3 or -4.

        The curve 681b returns an interval that contains $3/2$.
        This is because $E(\Q)$ is not saturated in $E(K)$.  The
        true index is $3$:
            sage: E = EllipticCurve('681b')
            sage: I = E.heegner_index(-8); I
            [1.4999942 .. 1.5000058]
            sage: 2*I
            [2.9999885 .. 3.0000115]

        In fact, whenever the returned index has a denominator
        of $2$, the true index is got by multiplying the returned
        index by $2$.  Unfortunately, this is not an if and only
        if condition, i.e., sometimes the index must be multiplied
        by $2$ even though the denominator is not $2$.
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
            return self.__adjust_heegner_index(ht/IR(reg))

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
        return self.__adjust_heegner_index(ht/IR(reg))

    def __adjust_heegner_index(self, a):
        r"""
        Take the square root of the interval that contains the Heegner
        index.

        EXAMPLES:
            sage: E = EllipticCurve('11a1')
            sage: a = RIF(sqrt(2))-1.4142135623730951
            sage: E._EllipticCurve_rational_field__adjust_heegner_index(a)
            [0.0000000... .. 1.490116...e-8]
        """
        if a.lower() < 0:
            a = IR((0, a.upper()))
        return a.sqrt()


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
                       first discriminant < -4 that satisfies the Heegner
                       hypothesis
            verbose (bool) -- (default: True)
            prec (int) -- (default: 5), use prec*sqrt(N) + 20 terms
                          of L-series in computations, where N is the
                          conductor.
            max_height (float) -- should be <= 21; bound on logarithmic
                                  naive height used in point searches.
                                  Make smaller to make this function
                                  faster, at the expense of possibly
                                  obtaining a worse answer.  A good
                                  range is between 13 and 21.

        OUTPUT:
            v -- list or int (bad primes or 0 or -1)
            D -- the discriminant that was used (this is useful if D was
                 automatically selected).

        EXAMPLES:
            sage: E = EllipticCurve('11a1')
            sage: E.heegner_index_bound(verbose=False)
            ([2], -7)
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

    padic_regulator = padics.padic_regulator

    padic_height_pairing_matrix = padics.padic_height_pairing_matrix

    padic_height = padics.padic_height
    padic_height_via_multiply = padics.padic_height_via_multiply

    padic_sigma = padics.padic_sigma
    padic_sigma_truncated = padics.padic_sigma_truncated

    padic_E2 = padics.padic_E2

    matrix_of_frobenius = padics.matrix_of_frobenius

    # def weierstrass_p(self):
    #         # TODO: add allowing negative valuations for power series
    #         return 1/t**2 + a1/t + rings.frac(1,12)*(a1-8*a2) -a3*t \
    #                - (a4+a1*a3)*t**2  + O(t**3)


    def mod5family(self):
        """
        Return the family of all elliptic curves with the same mod-5
        representation as self.

        EXAMPLES:
        sage: E=EllipticCurve('32a1')
        sage: E.mod5family()
        Elliptic Curve defined by y^2  = x^3 + 4*x over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        E = self.short_weierstrass_model()
        a = E.a4()
        b = E.a6()
        return mod5family.mod5family(a,b)

    def tate_curve(self, p):
        r"""
        Creates the Tate Curve over the $p$-adics associated to this elliptic curves.

        This Tate curve a $p$-adic curve with split multiplicative
        reduction of the form $y^2+xy=x^3+s_4 x+s_6$ which is
        isomorphic to the given curve over the algebraic closure of
        $\QQ_p$.  Its points over $\QQ_p$ are isomorphic to
        $\QQ_p^{\times}/q^{\Z}$ for a certain parameter $q\in\Z_p$.

        INPUT:

            p -- a prime where the curve has multiplicative reduction.


        EXAMPLES:
            sage: e = EllipticCurve('130a1')
            sage: e.tate_curve(2)
            2-adic Tate curve associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 - 33*x + 68 over Rational Field

        The input curve must have multiplicative reduction at the prime.
            sage: e.tate_curve(3)
            Traceback (most recent call last):
            ...
            ValueError: The elliptic curve must have multiplicative reduction at 3

        We compute with $p=5$:
            sage: T = e.tate_curve(5); T
            5-adic Tate curve associated to the Elliptic Curve defined by y^2 + x*y + y = x^3 - 33*x + 68 over Rational Field

        We find the Tate parameter $q$:
            sage: T.parameter(prec=5)
            3*5^3 + 3*5^4 + 2*5^5 + 2*5^6 + 3*5^7 + O(5^8)

        We compute the $L$-invariant of the curve:
            sage: T.L_invariant(prec=10)
            5^3 + 4*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 3*5^8 + 5^9 + O(5^10)
        """
        try:
            return self._tate_curve[p]
        except AttributeError:
            self._tate_curve = {}
        except KeyError:
            pass

        Eq = ell_tate_curve.TateCurve(self,p)
        self._tate_curve[p] = Eq
        return Eq


def cremona_curves(conductors):
    """
    Return iterator over all known curves (in database) with conductor
    in the list of conductors.

    EXAMPLES:
        sage: [(E.label(), E.rank()) for E in cremona_curves(srange(35,40))]
        [('35a1', 0),
        ('35a2', 0),
        ('35a3', 0),
        ('36a1', 0),
        ('36a2', 0),
        ('36a3', 0),
        ('36a4', 0),
        ('37a1', 1),
        ('37b1', 0),
        ('37b2', 0),
        ('37b3', 0),
        ('38a1', 0),
        ('38a2', 0),
        ('38a3', 0),
        ('38b1', 0),
        ('38b2', 0),
        ('39a1', 0),
        ('39a2', 0),
        ('39a3', 0),
        ('39a4', 0)]
    """
    if isinstance(conductors, (int,long, rings.RingElement)):
        conductors = [conductors]
    return sage.databases.cremona.CremonaDatabase().iter(conductors)

def cremona_optimal_curves(conductors):
    """
    Return iterator over all known optimal curves (in database) with
    conductor in the list of conductors.

    EXAMPLES:
        sage: [(E.label(), E.rank()) for E in cremona_optimal_curves(srange(35,40))]
        [('35a1', 0),
        ('36a1', 0),
        ('37a1', 1),
        ('37b1', 0),
        ('38a1', 0),
        ('38b1', 0),
        ('39a1', 0)]
    """
    if isinstance(conductors, (int,long,rings.RingElement)):
        conductors = [conductors]
    return sage.databases.cremona.CremonaDatabase().iter_optimal(conductors)


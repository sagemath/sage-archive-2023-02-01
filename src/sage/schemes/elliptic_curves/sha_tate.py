from sage.structure.sage_object import SageObject
from sage.rings.all import (
    Integer,
    RealField,
    RationalField,
    RIF)
from sage.misc.functional import log
from math import sqrt
from sage.misc.all import verbose
import sage.rings.arith as arith

factor = arith.factor
Q = RationalField()

class Sha(SageObject):
    """
    The Shafarevich-Tate group associated to an elliptic curve.

    If `E` is an elliptic curve over a field `K`, the Shafarevich-Tate group
    is the subgroup of `H^1(K,E)` which maps to zero under every global-to-local
    restriction map `H^1(K,E) \rightarrow H^1(K_v,E)`, one for each place `v`
    of `K`.

    EXAMPLES::
        sage: E = EllipticCurve('389a')
        sage: E.sha()
        Shafarevich-Tate group for the Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    """
    def __init__(self, E):
        self.E = E

    def __repr__(self):
        return "Shafarevich-Tate group for the " + repr(self.E)

    ########################################################################
    # Functions related to the BSD conjecture.
    ########################################################################

    def an_numerical(self, prec = None,
                         use_database=True, proof=None):
        """
        Return the numerical analytic order of Sha, which is
        a floating point number in all cases.

        INPUT:
            prec -- integer (default: 53) bits precision -- used
                    for the L-series computation, period,  regulator, etc.
            use_database -- whether the rank and generators should
                    be looked up in the database if possible.
            proof -- bool or None (default: None, see proof.[tab] or
                           sage.structure.proof) proof option passed
                    onto regulator and rank computation.

        NOTE: See also the an() command, which will return a
        provably correct integer when the rank is 0 or 1.

        WARNING: If the curve's generators are not known, computing
        them may be very time-consuming.  Also, computation of the
        L-series derivative will be time-consuming for large rank and
        large conductor, and the computation time for this may
        increase substantially at greater precision.  However, use of
        very low precision less than about 10 can cause the underlying
        pari library functions to fail.

        EXAMPLES:
            sage: EllipticCurve('11a').sha().an_numerical()
            1.00000000000000
            sage: EllipticCurve('37a').sha().an_numerical() # long time
            1.00000000000000
            sage: EllipticCurve('389a').sha().an_numerical() # long time
            1.00000000000000
            sage: EllipticCurve('66b3').sha().an_numerical()
            4.00000000000000
            sage: EllipticCurve('5077a').sha().an_numerical() # long time
            1.00000000000000

        A rank 4 curve:
            sage: EllipticCurve([1, -1, 0, -79, 289]).sha().an_numerical()   # long time
            1.00000000000000

        A rank 5 curve:
            sage: EllipticCurve([0, 0, 1, -79, 342]).sha().an_numerical(prec=10, proof=False) # long time -- about 30 seconds.
            1.0

            # See trac #1115
            sage: sha=EllipticCurve('37a1').sha()
            sage: [sha.an_numerical(prec) for prec in xrange(30,100,10)] # long time
            [1.0000000,
            1.0000000000,
            1.0000000000000,
            1.0000000000000000,
            1.0000000000000000000,
            1.0000000000000000000000,
            1.0000000000000000000000000]
        """
        if prec is None:
            prec = RealField().precision()
        RR = RealField(prec)
        prec2 = prec+2
        RR2 = RealField(prec2)
        try:
            an = self.__an_numerical
            if an.parent().precision() >= prec:
                return RR(an)
            else: # cached precision too low
                pass
        except AttributeError:
            pass
        r = Integer(self.E.rank(use_database=use_database, proof=proof))
        L = self.E.lseries().dokchitser(prec=prec2)
        Lr= RR2(L.derivative(1,r))  # L.derivative() returns a Complex
        Om = RR2(self.E.period_lattice().omega(prec2))
        Reg = self.E.regulator(use_database=use_database, proof=proof, precision=prec2)
        T = self.E.torsion_order()
        cp = self.E.tamagawa_product()
        Sha = RR((Lr*T*T)/(r.factorial()*Om*cp*Reg))
        self.__an_numerical = Sha
        return Sha

    def an(self, use_database=False):
        """
        Returns the Birch and Swinnerton-Dyer conjectural order of Sha
        as a provably corret integer, unless the analytic rank is > 1,
        in which case this function returns a numerical value.

        INPUT:
            use_database -- bool (default: False); if True, try to use any
            databases installed to lookup the analytic order of Sha, if
            possible.  The order of Sha is computed if it can't be looked up.

        This result is proved correct if the order of vanishing is 0
        and the Manin constant is <= 2.

        If the optional parameter use_database is True (default:
        False), this function returns the analytic order of Sha as
        listed in Cremona's tables, if this curve appears in Cremona's
        tables.

        EXAMPLES:
            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.sha().an()
            1
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.sha().an()
            1

            sage: EllipticCurve('14a4').sha().an()
            1
            sage: EllipticCurve('14a4').sha().an(use_database=True)   # will be faster if you have large Cremona database installed
            1

        The smallest conductor curve with nontrivial Sha:
            sage: E = EllipticCurve([1,1,1,-352,-2689])     # 66b3
            sage: E.sha().an()
            4

        The four optimal quotients with nontrivial Sha and conductor <= 1000:
            sage: E = EllipticCurve([0, -1, 1, -929, -10595])       # 571A
            sage: E.sha().an()
            4
            sage: E = EllipticCurve([1, 1, 0, -1154, -15345])       # 681B
            sage: E.sha().an()
            9
            sage: E = EllipticCurve([0, -1, 0, -900, -10098])       # 960D
            sage: E.sha().an()
            4
            sage: E = EllipticCurve([0, 1, 0, -20, -42])            # 960N
            sage: E.sha().an()
            4

        The smallest conductor curve of rank > 1:
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.sha().an()
            1.00000000000000

        The following are examples that require computation of the Mordell-Weil
        group and regulator:

            sage: E = EllipticCurve([0, 0, 1, -1, 0])                     # 37A  (rank 1)
            sage: E.sha().an()
            1

            sage: E = EllipticCurve("1610f3")
            sage: E.sha().an()
            4

        In this case the input curve is not minimal, and if this function didn't
        transform it to be minimal, it would give nonsense:
            sage: E = EllipticCurve([0,-432*6^2])
            sage: E.sha().an()
            1
        """
        if hasattr(self, '__an'):
            return self.__an
        if use_database:
            d = self.E.database_curve()
            if hasattr(d, 'db_extra'):
                self.__an = int(round(float(d.db_extra[4])))
                return self.__an

        # it's critical to switch to the minimal model.
        E = self.E.minimal_model()
        eps = E.root_number()
        if eps == 1:
            L1_over_omega = E.lseries().L_ratio()
            if L1_over_omega == 0:
                return self.an_numerical(use_database=use_database)
            T = E.torsion_subgroup().order()
            Sha = (L1_over_omega * T * T) / Q(E.tamagawa_product())
            try:
                Sha = Integer(Sha)
            except ValueError:
                raise RuntimeError, \
                      "There is a bug in an, since the computed conjectural order of Sha is %s, which is not an integer."%Sha
            if not arith.is_square(Sha):
                raise RuntimeError, \
                      "There is a bug in an, since the computed conjectural order of Sha is %s, which is not a square."%Sha
            E.__an = Sha
            self.__an = Sha
            return Sha

        else:  # rank > 0  (Not provably correct)
            L1, error_bound = E.lseries().deriv_at1(10*sqrt(E.conductor()) + 10)
            if abs(L1) < error_bound:
                s = self.an_numerical()
                E.__an = s
                self.__an = s
                return s

            regulator = E.regulator(use_database=use_database)   # this could take a *long* time; and could fail...?
            T = E.torsion_subgroup().order()
            omega = E.period_lattice().omega()
            Sha = int(round ( (L1 * T * T) / (E.tamagawa_product() * regulator * omega) ))
            try:
                Sha = Integer(Sha)
            except ValueError:
                raise RuntimeError, \
                      "There is a bug in an, since the computed conjectural order of Sha is %s, which is not an integer."%Sha
            if not arith.is_square(Sha):
                raise RuntimeError, \
                      "There is a bug in an, since the computed conjectural order of Sha is %s, which is not a square."%Sha
            E.__an = Sha
            self.__an = Sha
            return Sha

    def an_padic(self, p, prec=0):
        r"""
        Returns the conjectural order -- up to sign -- of Sha(E),
        according to the $p$-adic analogue of the BSD conjecture.

        INPUT:
            p -- a prime > 3
            prec (optional) -- the precision used in the computation of the
            p-adic L-Series

        OUTPUT:
            p-adic number -- that conjecturally equals $\#Sha(E)(p)$ or $-\#Sha(E)(p)$.

        NOTE:
            If prec is set to zero (default) then the precision is set so that
            at least the first p-adic digit of conjectural $\#Sha(E)(p)$ is
            determined.

        BUG:
            Currently for supersingular primes for curves of rank > 0, only the
            first digit will be correct. More is hard to compute anyway.

        EXAMPLES:

        Good ordinary examples:

            sage: EllipticCurve('11a1').sha().an_padic(5)  #rank 0
            1 + O(5^2)
            sage: EllipticCurve('43a1').sha().an_padic(5)  #rank 1
            1 + O(5)
            sage: EllipticCurve('389a1').sha().an_padic(5,4) #rank 2   (long time)
            1 + O(5^3)
            sage: EllipticCurve('858k2').sha().an_padic(7)  #rank 0, non trivial sha  (long time)
            7^2 + O(7^3)

        Exceptional cases:

            sage: EllipticCurve('11a1').sha().an_padic(11) #rank 0
            1 + O(11)

        The output maybe be only up to sign, as the following two examples illustrate:

            sage: EllipticCurve('123a1').sha().an_padic(41) #rank 1    (long time) -- random output (can be either 1 + O(41) or 40 + O(41)).
            40 + O(41)
            sage: EllipticCurve('817a1').sha().an_padic(43) #rank 2    (long time)
            42 + O(43)

        Supersingular cases:

            sage: EllipticCurve('34a1').sha().an_padic(5) # rank 0     (long time)
            1 + O(5^3)
            sage: EllipticCurve('43a1').sha().an_padic(7) # rank 1    (very long time -- nearly a minute)
            1 + O(7)
            sage: EllipticCurve('1483a1').sha().an_padic(5) # rank 2   (long time)
            1 + O(5)
        """
        try:
            return self.__an_padic[(p,prec)]
        except AttributeError:
            self.__an_padic = {}
        except KeyError:
            pass

        if self.E.has_cm():
            raise NotImplementedError, "I don't know about curves with complex multiplication."

        tam = self.E.tamagawa_product()
        tors = self.E.torsion_order()**2
        reg = self.E.padic_regulator(p)
        r = self.E.rank()

        #         shortcut to introduce later ::
        #   ( also we could avoid the computation of the modular symbols altogether )
        # if r == 0:
        #    lstar = lp.modular_symbol(0) * lp._quotient_of_periods
        #    sha = lstar/bsdp
        #    sha = lstar/bsdp[0] in the supersingular case

        lp = self.E.padic_lseries(p)

        if self.E.is_ordinary(p):
            K = reg.parent()
            lg = log(K(1+p))

            if (self.E.is_good(p) or self.E.ap(p) == -1):
                eps = (1-1/lp.alpha())**2
                # according to the p-adic BSD this should be equal to the leading term of the p-adic L-series:
                bsdp = tam * reg * eps/tors/lg**r
            else:
                r += 1   # exceptional zero
                eq = self.E.tate_curve(p)
                Li = eq.L_invariant()

                # according to the p-adic BSD (Mazur-Tate-Teitelbaum)
                # this should be equal to the leading term of the p-adic L-series:
                bsdp = tam * reg * Li/tors/lg**r


            v = bsdp.valuation()
            if v > 0:
                verbose("the prime is irregular.")

            # determine how much prec we need to prove at least the triviality of
            # the p-primary part od Sha

            if prec == 0:
                n = max(v,2)
                bounds = lp._prec_bounds(n,r+1)
                while bounds[r] <= v:
                    n += 1
                    bounds = lp._prec_bounds(n,r+1)
                verbose("set precision to %s"%n)
            else:
                n = max(2,prec)

            not_yet_enough_prec = True
            while not_yet_enough_prec:
                lps = lp.series(n,prec=r+1)
                lstar = lps[r]
                if (lstar != 0) or (prec != 0):
                    not_yet_enough_prec = False
                else:
                    n += 1
                    verbose("increased precision to %s"%n)

            shan = lstar/bsdp

        elif self.E.is_supersingular(p):
            K = reg[0].parent()
            lg = log(K(1+p))


            # according to the p-adic BSD this should be equal to the leading term of the D_p - valued
            # L-series :
            bsdp = tam /tors/lg**r * reg
            # note this is an element in Q_p^2

            verbose("the algebraic leading terms : %s"%bsdp)

            v = [bsdp[0].valuation(),bsdp[1].valuation()]

            if prec == 0:
                n = max(min(v)+2,3)
            else:
                n = max(3,prec)

            verbose("...computing the p-adic L-series")
            not_yet_enough_prec = True
            while not_yet_enough_prec:
                lps = lp.Dp_valued_series(n,prec=r+1)
                lstar = [lps[0][r],lps[1][r]]
                verbose("the leading terms : %s"%lstar)
                if (lstar[0] != 0 or lstar[1] != 0) or ( prec != 0):
                    not_yet_enough_prec = False
                else:
                    n += 1
                    verbose("increased precision to %s"%n)

            verbose("...putting things together")
            if bsdp[0] != 0:
                shan0 = lstar[0]/bsdp[0]
            else:
                shan0 = 0   # this should actully never happen
            if bsdp[1] != 0:
                shan1 = lstar[1]/bsdp[1]
            else:
                shan1 = 0   # this should conjecturally only happen when the rank is 0
            verbose("the two values for Sha : %s"%[shan0,shan1])

            # check consistency (the first two are only here to avoid a bug in the p-adic L-series
            # (namely the coefficients of zero-relative precision are treated as zero)
            if shan0 != 0 and shan1 != 0 and (shan0 - shan1 != 0 and shan0 + shan1 != 0):
                raise RuntimeError, "There must be a bug in the supersingular routines for the p-adic BSD."

            #take the better
            if shan1 == 0 or shan0.precision_relative() > shan1.precision_relative():
                shan = shan0
            else:
                shan = shan1

        else:
            raise ValueError, "The curve has to have semi-stable reduction at p."

        self.__an_padic[(p,prec)] = shan
        return shan


    def p_primary_bound(self, p):
        r"""
        Returns an upper bound of $\#Sha(E)(p)$.

        INPUT:
            p -- a prime > 3

        OUTPUT:
            integer -- power of p that bounds #Sha(E)(p) from above

        NOTE:
            The result is a proven upper bound on the order of $Sha(E)(p)$.
            So in particular it proves it finiteness even if the rank of
            the curve is larger than 1. Note also that this bound is sharp
            if one assumes the main conjecture of Iwasawa theory of
            elliptic curves (and this is known in certain cases).

        EXAMPLES:
            sage: e = EllipticCurve('858k2')
            sage: e.sha().p_primary_bound(3)           # long time
            0
            sage: e.sha().p_primary_bound(7)           # long time
            2
        """

        if self.E.is_ordinary(p) or self.E.is_good(p):
            shan = self.an_padic(p,prec = 0)
            if shan == 0:
                raise RuntimeError, "There is a bug in an_padic."
            S = shan.valuation()
            if not self.E.is_surjective(p) and not self.E.is_reducible(p):
                raise ValueError, "The mod-p Galois representation is neither surjective nor contained in a Borel group. Current knowledge about Euler systems does not provide an upper bound in this case. Try an_padic for a conjectural bound."
        else:
            raise ValueError, "The curve has to have semi-stable reduction at p."

        return S

    def two_selmer_bound(self):
        """
        Returns a bound on the dimension of Sha(E)[2], computed using
        a 2-descent.
        """
        S = self.E.selmer_rank_bound()
        r = self.E.rank()
        t = self.E.two_torsion_rank()
        b = S - r - t
        if b % 2 != 0:
            raise ArithmeticError, "There is a bug in two_selmer_bound since it's %s, but it must be even."%b
        return b

    def bound_kolyvagin(self, D=0, regulator=None,
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
            sage: E.sha().bound_kolyvagin()
            ([2], 1)
            sage: E = EllipticCurve('141a')
            sage: E.sha().an()
            1
            sage: E.sha().bound_kolyvagin()
            ([2, 7], 49)

        We get no information the curve has rank $2$.
            sage: E = EllipticCurve('389a')
            sage: E.sha().bound_kolyvagin()
            (0, 0)
            sage: E = EllipticCurve('681b')
            sage: E.sha().an()
            9
            sage: E.sha().bound_kolyvagin()
            ([2, 3], 9)

        """
        E = self.E
        if E.has_cm():
            return 0, 0

        if D == 0:
            D = -5
            while not E.satisfies_heegner_hypothesis(D):
                D -= 1

        if not E.satisfies_heegner_hypothesis(D):
            raise ArithmeticError, "Discriminant (=%s) must be a fundamental discriminant that satisfies the Heegner hypothesis."%D
        if D == -3 or D == -4:
            raise ArithmeticError, "Discriminant (=%s) must not be -3 or -4."%D
        eps = E.root_number()
        L1_vanishes = E.lseries().L1_vanishes()
        if eps == 1 and L1_vanishes:
            return 0, 0        # rank even hence >= 2, so Kolyvagin gives nothing.
        alpha = sqrt(abs(D))/(2*E.period_lattice().complex_area())
        F = E.quadratic_twist(D)
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
                raise RuntimeError, "Too many precision increases in bound_kolyvagin"
            if eps == 1:   # E has even rank
                verbose("Conductor of twist = %s"%F.conductor())
                LF1, err_F = F.lseries().deriv_at1(k_F)
                LE1, err_E = E.lseries().at1(k_E)
                err_F = max(err_F, MIN_ERR)
                err_E = max(err_E, MIN_ERR)
                if regulator != None:
                    hZ = regulator/2
                else:
                    hZ = F.regulator(use_database=True)/2
                #print  alpha * LE1 * LF1 / hZ
                I = RIF(alpha) * RIF(LE1-err_E,LE1+err_E) * RIF(LF1-err_F,LF1+err_F) / hZ
                #print I

            else:          # E has odd rank

                if regulator != None:
                    hZ = regulator/2
                else:
                    hZ = E.regulator(use_database=True)/2
                LE1, err_E = E.lseries().deriv_at1(k_E)
                LF1, err_F = F.lseries().at1(k_F)
                err_F = max(err_F, MIN_ERR)
                err_E = max(err_E, MIN_ERR)
                #I = alpha * LE1 * LF1 / hZ

                I = RIF(alpha) * RIF(LE1-err_E,LE1+err_E) * RIF(LF1-err_F,LF1+err_F) / hZ

            verbose('interval = %s'%I)
            t, n = I.is_int()
            if t:
                break
            elif I.absolute_diameter() < 1:
                raise RuntimeError, "Problem in bound_kolyvagin; square of index is not an integer -- D=%s, I=%s."%(D,I)
            verbose("Doubling bounds")
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
                    raise RuntimeError, "Problem in bound_kolyvagin; square of index is not a perfect square!  D=%s, I=%s, n=%s, e=%s."%(D,I,n,e)
                B.append(p)
            else:
                n /= 2**e  # replace n by its odd part
        if not ignore_nonsurj_hypothesis:
            for p, _ in E.non_surjective():
                B.append(p)
        B = list(set([int(x) for x in B]))
        B.sort()
        return B, n


    def bound_kato(self):
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
            sage: E.sha().bound_kato()
            [2, 3, 5]
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.sha().bound_kato()
            [2, 3, 5]
            sage: E = EllipticCurve([1,1,1,-352,-2689])     # 66B3
            sage: E.sha().bound_kato()
            [2, 3]

        For the following curve one really has 25 | $\\#Sha$ (by Grigorov-Stein paper):
            sage: E = EllipticCurve([1, -1, 0, -332311, -73733731])   # 1058D1
            sage: E.sha().bound_kato()                 # long time (about 1 second)
            [2, 3, 5]
            sage: E.non_surjective()                # long time (about 1 second)
            []

        For this one, Sha is divisible by 7.
            sage: E = EllipticCurve([0, 0, 0, -4062871, -3152083138])   # 3364C1
            sage: E.sha().bound_kato()                 # long time (< 10 seconds)
            [2, 3, 7]

        No information about curves of rank > 0:
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.sha().bound_kato()
            False
        """
        if self.E.has_cm():
            return False
        if self.E.lseries().L1_vanishes():
            return False
        B = [2,3]
        for p, _ in self.E.non_surjective():   # for p >= 5, mod-p surj => p-adic surj
            if p > 3:
                B.append(p)

        # The only other p that might divide B are those that divide
        # the integer 2*#E(Q)_tor^2 * L(E,1)/omega.  So we compute
        # that to sufficient precision to determine it.  Note that
        # we have to assume the Manin constant is <=2 in order to provably
        # compute L(E,1)/omega.
        for p, n in factor(self.an()):
            if n >= 2:    # use parity of Sha
                B.append(int(p))
        B = list(set(B))
        B.sort()
        return B

    def bound(self):
        """
        Compute a provably correct bound on the order of the Shafarevich-Tate
        group of this curve. The bound is a either False (no bound) or a list
        B of primes such that any divisor of Sha is in this list.

        EXAMPLES:
            sage: EllipticCurve('37a').sha().bound()
            ([2], 1)
        """
        if self.E.lseries().L1_vanishes():
            B = self.bound_kolyvagin()
        else:
            B = self.bound_kato()
        return B


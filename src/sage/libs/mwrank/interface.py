r"""
Cremona's MWRANK Library

This is the SAGE interface to John Cremona's MWRANK C++ library for
arithmetic on elliptic curves.  The classes defined in this module
give SAGE interpreter-level to access much of the functionality of
MWRANK.  For most purposes it is not necessary to directly use these
classes; instead, one can create an \class{EllipticCurve} and call
methods that are implemented using this module.

\note{This interface is a direct library-level interface to mwrank.}
"""

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import IntegerRing

def set_precision(n):
    """
    Set the global NTL real number precision.  This has a massive
    effect on the speed of mwrank calculations.  The default is n=15,
    but it might have to be increased if a computation fails.  In this
    case, one must recreate the mwrank curve from scratch after
    resetting this precision.

    NOTE -- this change is global and affects all of SAGE.

    INPUT:
        n -- long

    EXAMPLES:
        sage: mwrank_set_precision(20)
    """
    # don't want to load mwrank every time SAGE starts up, so we do
    # the import here.
    from sage.libs.mwrank.mwrank import set_precision
    set_precision(n)

class mwrank_EllipticCurve(SageObject):
    r"""
    The \class{mwrank_EllipticCurve} class represents an MWRANK
    elliptic curve.
    """
    def __init__(self, ainvs, verbose=False):
        r"""
        Create the mwrank elliptic curve with invariants
        \code{a_invs}, which is a list of $\leq 5$ \emph{integers}  $a_1$,
        $a_2$, $a_3$, $a_4$, and $a_6$.

        If strictly less than 5 invariants are given, then the first
        ones are set to 0, so, e.g., \code{[3,4]} means $a_1=a_2=a_3=0$ and
        $a_4=3$, $a_6=4$.

        INPUT:
            ainvs  --  a list of <= 5 integers
            verbose -- if True, then all Selmer group computations
                       will be verbose. (default: False).

        EXAMPLES:
        We create the elliptic curve $y^2 + y = x^3 + x^2 - 2x$.

            sage: e = mwrank_EllipticCurve([0, 1, 1, -2, 0])
            sage: e.ainvs()
            [0, 1, 1, -2, 0]

        This example illustrates that omitted $a$-invariants default to $0$:
            sage: e = mwrank_EllipticCurve([3, -4])
            sage: e
            y^2 = x^3 + 3*x - 4
            sage: e.ainvs()
            [0, 0, 0, 3, -4]

        The entries of the input list are coerced to \class{int}.
        If this is impossible then an error is raised:
            sage: e = mwrank_EllipticCurve([3, -4.8]); e
            Traceback (most recent call last):
            ...
            TypeError: ainvs must be a list of integers.

        When you enter a singular model you get an exception:
            sage: e = mwrank_EllipticCurve([0, 0])
            Traceback (most recent call last):
            ...
            ArithmeticError: Invariants (= [0, 0, 0, 0, 0]) do not describe an elliptic curve.
        """
        # import here to save time during startup (mwrank takes a while to init)

        from sage.libs.mwrank.mwrank import _Curvedata

        if not isinstance(ainvs, list) and len(ainvs) <= 5:
            raise TypeError, "ainvs must be a list of length at most 5."

        # Pad ainvs on the beginning by 0's, so e.g.
        # [a4,a6] works.
        ainvs = [0]*(5-len(ainvs)) + ainvs

        # Convert each entry to an int
        try:
            a_int = [IntegerRing()(x) for x in ainvs]
        except (TypeError, ValueError):
            raise TypeError, "ainvs must be a list of integers."
        self.__ainvs = a_int
        self.__curve = _Curvedata(a_int[0], a_int[1], a_int[2],
                                  a_int[3], a_int[4])

        if verbose:
            self.__verbose = True
        else:
            self.__verbose = False

        # place holders
        self.__saturate = -2  # not yet saturated

    def __reduce__(self):
        return mwrank_EllipticCurve, (self.__ainvs, self.__verbose)


    def set_verbose(self, verbose):
        """
        Set the verbosity of printing of output by the 2-descent and
        other functions.
        INPUT:
            verbosity -- bool; if True print lots of output
        """
        self.__verbose = verbose


    def _curve_data(self):
        return self.__curve

    def ainvs(self):
        return self.__ainvs

    def isogeny_class(self, verbose=False):
        return self.__curve.isogeny_class(verbose)

    def __repr__(self):
        a = self.ainvs()
        s = "y^2"
        if a[0] == -1:
            s += "- x*y "
        elif a[0] == 1:
            s += "+ x*y "
        elif a[0] != 0:
            s += "+ %s*x*y "%a[0]
        if a[2] == -1:
            s += " - y"
        elif a[2] == 1:
            s += "+ y"
        elif a[2] != 0:
            s += "+ %s*y"%a[2]
        s += " = x^3 "
        if a[1] == -1:
            s += "- x^2 "
        elif a[1] == 1:
            s += "+ x^2 "
        elif a[1] != 0:
            s += "+ %s*x^2 "%a[1]
        if a[3] == -1:
            s += "- x "
        elif a[3] == 1:
            s += "+ x "
        elif a[3] != 0:
            s += "+ %s*x "%a[3]
        if a[4] == -1:
            s += "-1"
        elif a[4] == 1:
            s += "+1"
        elif a[4] != 0:
            s += "+ %s"%a[4]
        s = s.replace("+ -","- ")
        return s


    def two_descent(self,
                    verbose = True,
                    selmer_only = False,
                    first_limit = 20,
                    second_limit = 8,
                    n_aux = -1,
                    second_descent = True):
        """
        Compute 2-descent data for this curve.

        INPUT:
            verbose     -- (default: True) print what mwrank is doing
            selmer_only -- (default: False) selmer_only switch
            first_limit -- (default: 20) firstlim is bound on |x|+|z|
            second_limit-- (default: 5)  secondlim is bound on log max {|x|,|z| },
                                         i.e. logarithmic
            n_aux       -- (default: -1) n_aux only relevant for general
                           2-descent when 2-torsion trivial; n_aux=-1 causes default
                           to be used (depends on method)
            second_descent -- (default: True) second_descent only relevant for
                           descent via 2-isogeny
        OUTPUT:
            Nothing -- nothing is returned
        """
        from sage.libs.mwrank.mwrank import _two_descent # import here to save time
        first_limit = int(first_limit)
        second_limit = int(second_limit)
        n_aux = int(n_aux)
        second_descent = int(second_descent)
        # TODO:  Don't allow limits above some value...???
        #   (since otherwise mwrank just sets limit tiny)
        self.__descent = _two_descent()
        self.__descent.do_descent(self.__curve,
                                      verbose,
                                      selmer_only,
                                      first_limit,
                                      second_limit,
                                      n_aux,
                                      second_descent)
        if not self.__two_descent_data().ok():
            raise RuntimeError, "A 2-descent didn't not complete successfully."

    def __two_descent_data(self):
        try:
            return self.__descent
        except AttributeError:
            self.two_descent(self.__verbose)
            return self.__descent

    def conductor(self):
        """
        Return the conductor of this curve, computed using Cremona's
        implementation of Tate's algorithm.

        NOTE: This is independent of PARI's.
        """
        return self.__curve.conductor()

    def rank(self):
        """
        Returns the rank of this curve, computed using 2-descent.
        """
        return self.__two_descent_data().getrank()

    def selmer_rank_bound(self):
        r"""
        Bound on the rank of the curve, computed using the 2-selmer
        group.  This is the rank of the curve minus the rank of the
        2-torsion, minus a number determined by whatever mwrank was
        able to determine related to $\Sha(E)[2]$ (e.g., using a
        second descent or if there is a rational $2$-torsion point,
        then there may be an isogeny to a curve with trivial
        $\Sha(E)$).  In many cases, this is the actual rank of the
        curve, but in general it is just $\geq$ the true rank.

        EXAMPLES:
        The following is the curve 960D1, which has rank 0,
        but Sha of order 4.

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.selmer_rank_bound()
            0

        In this case this was resolved using a second descent.

            sage: E = mwrank_EllipticCurve([0, -1, 0, -900, -10098])
            sage: E.two_descent(second_descent = False, verbose=False)
            sage: E.selmer_rank_bound()
            2

        Above, the \method{selmer_rank_bound} gives 0 instead of 2,
        because it knows Sha is nontrivial.  In contrast, for the
        curve 571A, also with rank 0 and $\Sha$ of order 4, we obtain
        a worse bound:

            sage: E = mwrank_EllipticCurve([0, -1, 1, -929, -10595])
            sage: E.selmer_rank_bound()
            2
            sage: E.rank()
            0
        """
        return self.__two_descent_data().getselmer()

    def regulator(self):
        """
        Return the regulator of the saturated Mordell-Weil group.

            sage: E = mwrank_EllipticCurve([0, 0, 1, -1, 0])
            sage: E.regulator()
            0.05111140823996884

        """
        self.saturate()
        if not self.certain():
            raise RuntimeError, "Unable to saturate Mordell-Weil group."
        R = self.__two_descent_data().regulator()
        return float(R)

    def saturate(self, bound=-1):
        """
        Compute the saturation of the Mordell-Weil group at all
        primes up to bound.

        INPUT:
            bound -- int (default: -1)   -1 saturate at *all* primes,
                                          0 -- do not saturate
                                          n -- saturate at least at all primes <= n.
        """
        bound = int(bound)
        if self.__saturate < bound:
            self.__two_descent_data().saturate(bound)
            self.__saturate = bound

    def gens(self):
        """
        Return a list of the generators for the Mordell-Weil group.
        """
        self.saturate()
        return eval(self.__two_descent_data().getbasis().replace(":",","))

    def certain(self):
        r"""
        True if the last \method{two_descent} call provably correctly
        computed the rank.  If \method{two_descent} hasn't been
        called, then it is first called by \method{certain}
        using the default parameters.

        EXAMPLES:
        A $2$-descent does not determine $E(\Q)$ with certainty
        for the curve $y^2 + y = x^3 - x^2 - 120x - 2183$.

            sage: E = mwrank_EllipticCurve([0, -1, 1, -120, -2183])
            sage: E.two_descent(False)
            ...
            sage: E.certain()
            False
            sage: E.rank()
            0

        The rank of $E$is actually 0 (as one could see by computing
        the L-function), but $\Sha$ has order 4 and the $2$-torsion is
        trivial, so mwrank does not conclusively determine the rank.

            sage: E.selmer_rank_bound()
            2
        """
        return bool(self.__two_descent_data().getcertain())

    #def fullmw(self):
    #    return self.__two_descent_data().getfullmw()

    def CPS_height_bound(self):
        r"""
        Return the Cremona-Prickett-Siksek height bound.  This is a
        floating point number $B$ such that if $P$ is a point on the curve,
        then the naive logarithmetic height of $P$ is off from the
        canonical height by at most $B$.

        \begin{notice}
        We assume the model is minimal!
        \end{notice}
        """
        return self.__curve.cps_bound()

    def silverman_bound(self):
        r"""
        Return the Silverman height bound.  This is a number $B$ such
        that if $P$ is a point on the curve, then the naive
        logarithmetic height of $P$ is off from the canonical height by
        at most $B$.

        \begin{notice}
        We assume the model is minimal!
        \end{notice}
        """
        return self.__curve.silverman_bound()


class mwrank_MordellWeil(SageObject):
    r"""
    The \class{mwrank_MordellWeil} class represents a subgroup of a
    Mordell-Weil group.  Use this class to saturate a specified list
    of points on an \class{mwrank_EllipticCurve}, or to search for
    points up to some bound.
    """
    def __init__(self, curve, verbose=True, pp=1, maxr=999):
        r"""
        Create a \class{mwrank_MordellWeil} instance.

        INPUT:
            curve -- \class{mwrank_EllipticCurve} instance
            verbose -- bool
            pp -- int
            maxr -- int
        """
        if not isinstance(curve, mwrank_EllipticCurve):
            raise TypeError, "curve (=%s) must be an mwrank_EllipticCurve"%curve
        self.__curve = curve
        self.__verbose = verbose
        self.__pp = pp
        self.__maxr = maxr
        if verbose:
            verb = 1
        else:
            verb = 0
        from sage.libs.mwrank.mwrank import _mw # import here to save time
        self.__mw = _mw(curve._curve_data(), verb, pp, maxr)

    def __reduce__(self):
        return mwrank_MordellWeil, (self.__curve, self.__verbose, self.__pp, self.__maxr)

    def __repr__(self):
        return "Subgroup of Mordell Weil group: %s"%self.__mw

    def process(self, v, sat=0):
        """
        This function allows one to add points to a mwrank_MordellWeil object.

        Process points in the list v, with saturation at primes up to
        sat.  If sat = 0 (the default), then saturate at all primes.

        INPUT:

            v -- a point (3-tuple of ints), or a
                 list of 3-tuples of integers, which define points on the curve.

            sat -- int, saturate at primes up to sat, or at all primes if sat=0.

        """
        if not isinstance(v, list):
            raise TypeError, "v (=%s) must be a list"%v
        sat = int(sat)
        for P in v:
            if not isinstance(P, (list,tuple)) or len(P) != 3:
                raise TypeError, "v (=%s) must be a list of 3-tuples of ints"%v
            self.__mw.process(P, sat)

    def regulator(self):
        """
        Return the regulator of the points in this subgroup of
        the Mordell-Weil group.
        """
        return self.__mw.regulator()

    def rank(self):
        """
        Return the rank of this subgroup of the Mordell-Weil group.
        """
        return self.__mw.getrank()

    def saturate(self, max_prime=-1, odd_primes_only=False):
        r"""
        Saturate this subgroup of the Mordell-Weil group.

        INPUT:
            max_prime (int) -- (default: 97), saturation is performed
                               for all primes up to max_prime

            odd_primes_only (bool) -- only do saturation at odd primes

        OUTPUT:
            ok (bool) -- True if and only if the saturation
                         is provably correct at \emph{all} primes.
            index (int) -- The index of the group generated by
                           points in their saturation
            saturation (list) -- list of points that form
                                 a basis for the saturation

        \begin{notice}
        We emphasize that if this function returns True as the first
        return argument, then the points it found are saturated at
        \emph{all} primes, i.e., saturating at the primes up to
        \code{max_prime} are sufficient to saturate at all primes.
        Note that the function might not have needed to saturate at
        all primes up to \code{max_prime}.
        It has worked out what prime you need to saturate up to,
        and that prime is $\leq $ \code{max_prime}.

        \end{notice}

        \begin{notice}
        Currently (July 2005), this does not remember the result of
        calling search.  So calling search up to height 20 then
        calling saturate results in another search up to height 18.
        \end{notice}


        """
        ok, index, unsat = self.__mw.saturate(int(max_prime), odd_primes_only)
        return bool(ok), str(index), unsat

    def search(self, height_limit=18, verbose=False):
        r"""
        Search for new points, and add them to this subgroup of the
        Mordell-Weil group.

        INPUT:
            height_limit -- float (default: 18) search up to
                            this logarithmetic height.
                   On 32-bit machines, h_lim MUST be < 21.48 else
                   $\exp(h_lim)>2^31$ and overflows.
        """
        # On 64-bit machines, h_lim must be < 43.668, but
        #           it's unlikely (in 2005) that you would ever
        #          search higher than 21.
        height_limit = float(height_limit)
        if height_limit >= 21.4:
            raise ValueError, "The height limit must be < 21.4."

        moduli_option = 0  # Use Stoll's sieving program... see strategies in ratpoints-1.4.c

        ##            moduli_option -- int (default: 0); if > 0; a flag used to determine
        ##                    the moduli that are used in sieving
        ##                       1 -- first 10 odd primes; the first one you
        ##                            would think of.
        ##                       2 -- three composites; $2^6\cdot 3^4$, ... TODO
        ##                             (from German mathematician J. Gebel;
        ##                               personal conversation about SIMATH)
        ##                       3 -- nine prime powers; $2^5, \ldots$;
        ##                            like 1 but includes powers of small primes
        ##                    TODO: Extract the meaning from mwprocs.cc; line 776 etc.

        verbose == bool(verbose)
        self.__mw.search(height_limit, moduli_option, verbose)

    def points(self):
        """
        Return a list of the generating points in this Mordell-Weil
        group.
        """
        return eval(self.__mw.getbasis().replace(":",","))



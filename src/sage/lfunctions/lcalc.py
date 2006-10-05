r"""
Rubinstein's $L$-function Calculator

This a standard part of \sage.  This interface provides complete
access to Rubinstein's lcalc calculator with extra PARI functionality
compiled in.

\note{Each call to \code{lcalc} runs a complete \code{lcalc} process.
On a typical Linux system, this entails about 0.3 seconds overhead.}

AUTHOR:
    -- Michael Rubinstein (2005): released under GPL the C++ program lcalc
    -- William Stein (2006-03-05): wrote SAGE interface to lcalc
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, weakref
from sage.structure.sage_object import SageObject
from sage.misc.all import pager
import sage.rings.all

prec = 32

class LCalc(SageObject):
    r"""
    Rubinstein's $L$-functions Calculator

    Type \code{lcalc.[tab]} for a list of useful commands that are
    implemented using the command line interface, but return objects
    that make sense in SAGE.  For each command the possible inputs
    for the L-function are:
    \begin{itemize}
       \item   \code{ ''} -- (default) the Riemann zeta function
       \item   \code{ 'tau'} -- the L function of the Ramanujan delta function
       \item   elliptic curve E -- where E is an elliptic curve over $\QQ$; defines $L(E,s)$
    \end{itemize}

    You can also use the complete command-line interface of
    Rubinstein's $L$-functions calculations program via this class.
    Type \code{lcalc.help()} for a list of commands and how to call
    them.
    """
    def _repr_(self):
        return "Rubinsteins L-function Calculator"

    def __call__(self, args):
        cmd = 'lcalc %s'%args
        return os.popen(cmd).read().strip()

    def _compute_L(self, L):
        if isinstance(L, str):
            if L == 'tau':
                return '--tau'
            return L
        import sage.schemes.all
        if sage.schemes.all.is_EllipticCurve(L):
            if L.base_ring() == sage.rings.all.RationalField():
                L = L.minimal_model()
                return '-e --a1 %s --a2 %s --a3 %s --a4 %s --a6 %s'%tuple(L.a_invariants())
        raise TypeError, "$L$-function of %s not known"%L

    def help(self):
        try:
            h = self.__help
        except AttributeError:
            h = "-"*70 + '\n'
            h += "   Call lcalc with one argument, e.g., \n"
            h += "      sage: lcalc('--tau -z 1000')\n"
            h += "   is translated into the command line\n"
            h += "      $ lcalc --tau -z 1000\n"
            h += "-"*70 + '\n'
            h += '\n' + self('--help')
            self.__help = h
        pager()(h)

    def zeros(self, n, L=''):
        """
        Return the imaginary parts of the first $n$ nontrivial zeros
        of the $L$-function in the upper half plane, as 32-bit reals.

        INPUT:
            n -- integer
            L -- defines $L$-function (default: Riemann zeta function)

        This function also checks the Riemann Hypothesis and makes
        sure no zeros are missed.  This means it looks for several
        dozen zeros to make sure none have been missed before
        outputting any zeros at all, so takes longer than
        \code{self.zeros_of_zeta_in_interval(...)}.

        EXAMPLES:
            sage: lcalc.zeros(4)                           # long
            [14.134725142, 21.022039637, 25.010857582, 30.424876124]
            sage: lcalc.zeros(5, L='--tau')                # long
            [9.2223794013, 13.907549862, 17.442776978, 19.656513140, 22.336103640]
            sage: lcalc.zeros(3, EllipticCurve('37a'))     # long
            [0.00000000000, 5.0031700134, 6.8703912161]
        """
        L = self._compute_L(L)
        RR = sage.rings.all.RealField(prec)
        X = self('-z %s %s'%(int(n), L))
        return [RR(z) for z in X.split()]

    def zeros_in_interval(self, x, y, stepsize, L=''):
        r"""
        Return the imaginary parts of (most of) the nontrivial zeros
        of the $L$-function on the line $\Re(s)=1/2$ with
        positive imaginary part between $x$ and $y$, along with a
        technical quantity for each.

        INPUT:
            x, y, stepsize -- positive floating point numbers
            L -- defines $L$-function (default: Riemann zeta function)

        OUTPUT:
            list of pairs (zero, S(T)).

        Rubinstein writes: The first column outputs the imaginary part
        of the zero, the second column a quantity related to $S(T)$ (it
        increases roughly by 2 whenever a sign change, i.e. pair of
        zeros, is missed). Higher up the critical strip you should use
        a smaller stepsize so as not to miss zeros.

        EXAMPLES:
            sage: lcalc.zeros_in_interval(10, 30, 0.1)
            [(14.134725142, 0.18467291567), (21.022039637, -0.067789328954), (25.010857582, -0.055587278126)]
        """
        L = self._compute_L(L)
        RR = sage.rings.all.RealField(prec)
        X = self('--zeros-interval -x %s -y %s --stepsize=%s %s'%(
            float(x), float(y), float(stepsize), L))
        return [tuple([RR(z) for z in t.split()]) for t in X.split('\n')]

    def value(self, s, L=''):
        r"""
        Return $L(s)$ for $s$ a complex number.

        INPUT:
            s -- complex number
            L -- defines $L$-function (default: Riemann zeta function)

        EXAMPLES:
            sage: I = CC.0
            sage: lcalc.value(0.5 + 100*I)
            2.6926198853 - 0.020386029602*I

        Note, SAGE can also compute zeta at complex numbers (using
        the PARI C library):
            sage: (0.5 + 100*I).zeta()
            2.6926198856813239 - 0.020386029602598162*I
        """
        L = self._compute_L(L)
        CC = sage.rings.all.ComplexField(prec)
        s = CC(s)
        x, y = self('-v -x %s -y %s %s'%(s.real(), s.imag(), L)).split()
        return CC((float(x), float(y)))

    def values_along_line(self, s0, s1, number_samples, L=''):
        r"""
        Return values of $L(s)$ at \code{number_samples}
        equally-spaced sample points along the line from $s_0$ to
        $s_1$ in the complex plane.

        INPUT:
            s0, s1 -- complex numbers
            number_samples -- integer
            L -- defines $L$-function (default: Riemann zeta function)

        OUTPUT:
            list -- list of pairs (s, zeta(s)), where the s are
                    equally spaced sampled points on the line from
                    s0 to s1.

        EXAMPLES:
            sage: I = CC.0
            sage: lcalc.values_along_line(0.5, 0.5+20*I, 5)
            [(0.50000000000, -1.4603545088), (0.50000000000 + 4.0000000000*I, 0.60678376444 + 0.091112139984*I), (0.50000000000 + 8.0000000000*I, 1.2416151054 + 0.36004758836*I), (0.50000000000 + 12.000000000*I, 1.0159366508 - 0.74511247221*I), (0.50000000000 + 16.000000000*I, 0.93854540843 + 1.2165878159*I)]
        """
        L = self._compute_L(L)
        CC = sage.rings.all.ComplexField(prec)
        s0 = CC(s0)
        s1 = CC(s1)
        v = self('--value-line-segment -x %s -y %s -X %s -Y %s --number-samples %s %s'%(
            (s0.real(), s0.imag(), s1.real(), s1.imag(), int(number_samples), L)))
        w = []
        for a in v.split('\n'):
            x0,y0,x1,y1 = a.split()
            w.append((CC(x0,y0), CC(x1,y1)))
        return w

    def twist_values(self, s, dmin, dmax, L=''):
        r"""
        Return values of $L(s, \chi_d)$ for each quadratic
        character $\chi_d$ for $d_{\min} \leq d \leq d_{\max}$.

        INPUT:
            s -- complex numbers
            dmin -- integer
            dmax -- integer
            L -- defines $L$-function (default: Riemann zeta function)

        OUTPUT:
            list -- list of pairs (d, L(s,chi_d))

        EXAMPLES:
            sage: lcalc.twist_values(0.5, -10, 10)
            [(-8, 1.1004214096), (-7, 1.1465856670), (-4, 0.66769145709), (-3, 0.48086755769), (5, 0.23175094748), (8, 0.37369171286)]
        """
        L = self._compute_L(L)
        CC = sage.rings.all.ComplexField(prec)
        Z = sage.rings.all.Integer
        s = CC(s)
        typ = '--twist-quadratic'
        dmin = int(dmin)
        dmax = int(dmax)
        v = self('-v -x %s -y %s %s --start %s --finish %s %s'%(
            (s.real(), s.imag(), typ, dmin, dmax, L)))
        w = []
        if len(v) == 0:
            return w
        if len(v) == 0:
            return w
        for a in v.split('\n'):
            d,x,y = a.split()
            w.append((Z(d), CC(x,y)))
        return w

    def twist_zeros(self, n, dmin, dmax, L=''):
        r"""
        Return first $n$ real parts of nontrivial zeros for each
        quadratic character $\chi_d$ for $d_{\min} \leq d \leq
        d_{\max}$.

        INPUT:
            n -- integer
            dmin -- integer
            dmax -- integer
            L -- defines $L$-function (default: Riemann zeta function)

        OUTPUT:
            dict -- keys are the discriminants $d$, and
                    values are list of corresponding zeros.

        EXAMPLES:
            sage: lcalc.twist_zeros(3, -3, 6)
            {-3: [8.0397371575, 11.249206208, 15.704619177], 5: [6.6484533455, 9.8314444311, 11.958845627]}
        """
        L = self._compute_L(L)
        RR = sage.rings.all.RealField(prec)
        Z = sage.rings.all.Integer
        typ = '--twist-quadratic'
        n = int(n)
        v = self('-z %s %s --start %s --finish %s %s'%(
            (n, typ, dmin, dmax, L)))
        w = {}
        if len(v) == 0:
            return w
        for a in v.split('\n'):
            d, x = a.split()
            x = RR(x)
            d = Z(d)
            if w.has_key(d):
                w[d].append(x)
            else:
                w[d] = [x]
        return w

    def analytic_rank(self, L=''):
        r"""
        Return the analytic rank of the $L$-function at the central
        critical point.

        INPUT:
            L -- defines $L$-function (default: Riemann zeta function)

        OUTPUT:
            integer

        \note{Of course this is not provably correct in general, since
        it is an open problem to compute analytic ranks provably
        correctly in general.}

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: lcalc.analytic_rank(E)
            1
        """
        L = self._compute_L(L)
        Z = sage.rings.all.Integer
        s = self('--rank-compute %s'%L)
        i = s.find('equals')
        return Z(s[i+6:])



# An instance
lcalc = LCalc()



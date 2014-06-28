r"""
Rubinstein's `L`-function Calculator

This interface provides complete
access to Rubinstein's lcalc calculator with extra PARI
functionality compiled in
and is a standard part of Sage.

.. note::

   Each call to ``lcalc`` runs a complete
   ``lcalc`` process. On a typical Linux system, this
   entails about 0.3 seconds overhead.

AUTHORS:

- Michael Rubinstein (2005): released under GPL the C++ program lcalc

- William Stein (2006-03-05): wrote Sage interface to lcalc
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import os
from sage.structure.sage_object import SageObject
from sage.misc.all import pager
import sage.rings.all
import sage.schemes.elliptic_curves.ell_generic

prec = 32

class LCalc(SageObject):
    r"""
    Rubinstein's `L`-functions Calculator

    Type ``lcalc.[tab]`` for a list of useful commands that
    are implemented using the command line interface, but return
    objects that make sense in Sage. For each command the possible
    inputs for the L-function are:


    -  ``"`` - (default) the Riemann zeta function

    -  ``'tau'`` - the L function of the Ramanujan delta
       function

    -  elliptic curve E - where E is an elliptic curve over
       `\mathbb{Q}`; defines `L(E,s)`


    You can also use the complete command-line interface of
    Rubinstein's `L`-functions calculations program via this
    class. Type ``lcalc.help()`` for a list of commands and
    how to call them.
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
        if sage.schemes.elliptic_curves.ell_generic.is_EllipticCurve(L):
            if L.base_ring() == sage.rings.all.RationalField():
                L = L.minimal_model()
                return '-e --a1 %s --a2 %s --a3 %s --a4 %s --a6 %s'%tuple(L.a_invariants())
        raise TypeError("$L$-function of %s not known"%L)

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
        Return the imaginary parts of the first `n` nontrivial
        zeros of the `L`-function in the upper half plane, as
        32-bit reals.

        INPUT:


        -  ``n`` - integer

        -  ``L`` - defines `L`-function (default:
           Riemann zeta function)


        This function also checks the Riemann Hypothesis and makes sure no
        zeros are missed. This means it looks for several dozen zeros to
        make sure none have been missed before outputting any zeros at all,
        so takes longer than
        ``self.zeros_of_zeta_in_interval(...)``.

        EXAMPLES::

            sage: lcalc.zeros(4)                           # long time
            [14.1347251, 21.0220396, 25.0108576, 30.4248761]
            sage: lcalc.zeros(5, L='--tau')                # long time
            [9.22237940, 13.9075499, 17.4427770, 19.6565131, 22.3361036]
            sage: lcalc.zeros(3, EllipticCurve('37a'))     # long time
            [0.000000000, 5.00317001, 6.87039122]
        """
        L = self._compute_L(L)
        RR = sage.rings.all.RealField(prec)
        X = self('-z %s %s'%(int(n), L))
        return [RR(z) for z in X.split()]

    def zeros_in_interval(self, x, y, stepsize, L=''):
        r"""
        Return the imaginary parts of (most of) the nontrivial zeros of the
        `L`-function on the line `\Re(s)=1/2` with positive
        imaginary part between `x` and `y`, along with a
        technical quantity for each.

        INPUT:


        -  ``x, y, stepsize`` - positive floating point
           numbers

        -  ``L`` - defines `L`-function (default:
           Riemann zeta function)


        OUTPUT: list of pairs (zero, S(T)).

        Rubinstein writes: The first column outputs the imaginary part of
        the zero, the second column a quantity related to `S(T)`
        (it increases roughly by 2 whenever a sign change, i.e. pair of
        zeros, is missed). Higher up the critical strip you should use a
        smaller stepsize so as not to miss zeros.

        EXAMPLES::

            sage: lcalc.zeros_in_interval(10, 30, 0.1)
            [(14.1347251, 0.184672916), (21.0220396, -0.0677893290), (25.0108576, -0.0555872781)]
        """
        L = self._compute_L(L)
        RR = sage.rings.all.RealField(prec)
        X = self('--zeros-interval -x %s -y %s --stepsize=%s %s'%(
            float(x), float(y), float(stepsize), L))
        return [tuple([RR(z) for z in t.split()]) for t in X.split('\n')]

    def value(self, s, L=''):
        r"""
        Return `L(s)` for `s` a complex number.

        INPUT:


        -  ``s`` - complex number

        -  ``L`` - defines `L`-function (default:
           Riemann zeta function)


        EXAMPLES::

            sage: I = CC.0
            sage: lcalc.value(0.5 + 100*I)
            2.69261989 - 0.0203860296*I

        Note, Sage can also compute zeta at complex numbers (using the PARI
        C library)::

            sage: (0.5 + 100*I).zeta()
            2.69261988568132 - 0.0203860296025982*I
        """
        L = self._compute_L(L)
        CC = sage.rings.all.ComplexField(prec)
        s = CC(s)
        x, y = self('-v -x %s -y %s %s'%(s.real(), s.imag(), L)).split()
        return CC((float(x), float(y)))

    def values_along_line(self, s0, s1, number_samples, L=''):
        r"""
        Return values of `L(s)` at ``number_samples``
        equally-spaced sample points along the line from `s_0` to
        `s_1` in the complex plane.

        INPUT:


        -  ``s0, s1`` - complex numbers

        -  ``number_samples`` - integer

        -  ``L`` - defines `L`-function (default:
           Riemann zeta function)


        OUTPUT:


        -  ``list`` - list of pairs (s, zeta(s)), where the s
           are equally spaced sampled points on the line from s0 to s1.


        EXAMPLES::

            sage: I = CC.0
            sage: lcalc.values_along_line(0.5, 0.5+20*I, 5)
            [(0.500000000, -1.46035451), (0.500000000 + 4.00000000*I, 0.606783764 + 0.0911121400*I), (0.500000000 + 8.00000000*I, 1.24161511 + 0.360047588*I), (0.500000000 + 12.0000000*I, 1.01593665 - 0.745112472*I), (0.500000000 + 16.0000000*I, 0.938545408 + 1.21658782*I)]

        Sometimes warnings are printed (by lcalc) when this command is
        run::

            sage: E = EllipticCurve('389a')
            sage: E.lseries().values_along_line(0.5, 3, 5)
            [(0.000000000, 0.209951303),
             (0.500000000, -...e-16),
             (1.00000000, 0.133768433),
             (1.50000000, 0.360092864),
             (2.00000000, 0.552975867)]
        """
        L = self._compute_L(L)
        CC = sage.rings.all.ComplexField(prec)
        s0 = CC(s0)
        s1 = CC(s1)
        v = self('--value-line-segment -x %s -y %s -X %s -Y %s --number-samples %s %s'%(
            (s0.real(), s0.imag(), s1.real(), s1.imag(), int(number_samples), L)))
        w = []
        for a in v.split('\n'):
            try:
                x0,y0,x1,y1 = a.split()
                w.append((CC(x0,y0), CC(x1,y1)))
            except ValueError:
                print 'lcalc: ', a
        return w

    def twist_values(self, s, dmin, dmax, L=''):
        r"""
        Return values of `L(s, \chi_k)` for each quadratic
        character `\chi_k` whose discriminant `d` satisfies
        `d_{\min} \leq d \leq d_{\max}`.

        INPUT:


        -  ``s`` - complex numbers

        -  ``dmin`` - integer

        -  ``dmax`` - integer

        -  ``L`` - defines `L`-function (default:
           Riemann zeta function)


        OUTPUT:


        -  ``list`` - list of pairs (d, L(s,chi_d))


        EXAMPLES::

            sage: lcalc.twist_values(0.5, -10, 10)
            [(-8, 1.10042141), (-7, 1.14658567), (-4, 0.667691457), (-3, 0.480867558), (5, 0.231750947), (8, 0.373691713)]
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
        Return first `n` real parts of nontrivial zeros for each
        quadratic character `\chi_k` whose discriminant `d` satisfies
        `d_{\min} \leq d \leq d_{\max}`.

        INPUT:


        -  ``n`` - integer

        -  ``dmin`` - integer

        -  ``dmax`` - integer

        -  ``L`` - defines `L`-function (default:
           Riemann zeta function)


        OUTPUT:


        -  ``dict`` - keys are the discriminants `d`,
           and values are list of corresponding zeros.


        EXAMPLES::

            sage: lcalc.twist_zeros(3, -3, 6)
            {-3: [8.03973716, 11.2492062, 15.7046192], 5: [6.64845335, 9.83144443, 11.9588456]}
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
            if d in w:
                w[d].append(x)
            else:
                w[d] = [x]
        return w

    def analytic_rank(self, L=''):
        r"""
        Return the analytic rank of the `L`-function at the central
        critical point.

        INPUT:


        -  ``L`` - defines `L`-function (default:
           Riemann zeta function)


        OUTPUT: integer

        .. note::

           Of course this is not provably correct in general, since it
           is an open problem to compute analytic ranks provably
           correctly in general.

        EXAMPLES::

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



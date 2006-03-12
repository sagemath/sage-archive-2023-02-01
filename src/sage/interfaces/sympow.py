r"""
Watkins Symmetric Power $L$-function Calculator

SYMPOW is a package to compute special values of symmetric power
elliptic curve L-functions. It can compute up to about 64 digits of
precision.  This interface provides complete access to sympow, which
is a standard part of \sage (and includes the extra data files).

\note{Each call to \code{sympow} runs a complete \code{sympow}
process.}

AUTHOR:
    -- Mark Watkins  (2005-2006): wrote and released sympow
    -- William Stein (2006-03-05): wrote SAGE interface

ACKNOWLEDGEMENT (from sympow readme):
\begin{itemize}
\item The quad-double package was modified from David Bailey's package:
\url{http://crd.lbl.gov/~dhbailey/mpdist/}

\item The \code{squfof} implementation was modified from Allan Steel's
version of Arjen Lenstra's original LIP-based code.

\item The \code{ec_ap} code was originally written for the kernel of
MAGMA, but was modified to use small integers when possible.

\item SYMPOW was originally developed using PARI, but due to licensing
difficulties, this was eliminated. SYMPOW also does not use the
standard math libraries unless Configure is run with the -lm option.
SYMPOW still uses GP to compute the meshes of inverse Mellin
transforms (this is done when a new symmetric power is added to
datafiles).

\end{itemize}
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, weakref

from sage.ext.sage_object import SageObject
from sage.misc.all import pager, verbose
import sage.rings.all

class Sympow(SageObject):
    r"""
    Watkins Symmetric Power $L$-function Calculator

    Type \code{sympow.[tab]} for a list of useful commands that are
    implemented using the command line interface, but return objects
    that make sense in SAGE.

    You can also use the complete command-line interface of sympow via
    this class.  Type \code{sympow.help()} for a list of commands and
    how to call them.
    """
    def _repr_(self):
        return "Watkins Symmetric Power L-function Calculator"

    def __call__(self, args):
        cmd = 'sympow %s'%args
        v = os.popen(cmd).read().strip()
        verbose(v, level=2)
        return v

    def _fix_err(self, err):
        w = err
        j = w.rfind('./sympow')
        if j != -1:
            w = w[:j-1] + "sympow('" + w[j+9:] + ')'
        return w

    def _curve_str(self, E):
        return '-curve "%s"'%(str(E.minimal_model().a_invariants()).replace(' ',''))

    def L(self, E, n, prec):
        r"""
        Return $L(\Sym^{(n)}(E, \text{edge}))$ to prec digits
        of precision.

        INPUT:
            E -- elliptic curve
            n -- integer
            prec -- integer

        OUTPUT:
            string -- real number to prec digits of precision as a string.

        \note{Before using this function for the first time for
        a given $n$, you may have to type \code{sympow('-new_data <n>')},
        where \code{<n>} is replaced by your value of $n$.}

        If you would like to see the extensive output sympow prints
        when running this function, just type \code{set_verbose(2)}.

        EXAMPLES:
            sage.: a = sympow.L(EllipticCurve('11a'), 2, 16); a
            '1.057599244590958E+00'
            sage.: RR(a)
            1.0575992445909579
        """
        if prec > 64:
            raise ValueError, "prec (=%s) must be at most 64"%prec
        if prec < 1:
            raise ValueError, "prec (=%s) must be at least 1"%prec
        v = self('-sp %sp%s %s'%(n, prec, self._curve_str(E)))
        i = v.rfind(': ')
        if i == -1:
            print self._fix_err(v)
            raise RuntimeError, "failed to compute symmetric power"
        x = v[i+2:]
        return x

    def Lderivs(self, E, n, prec, d):
        r"""
        Return $0$th to $d$th derivatives of
        $L(\Sym^{(n)}(E,\text{edge}))$ to prec digits
        of precision.

        INPUT:
            E -- elliptic curve
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
            sage.: print sympow.Lderivs(EllipticCurve('11a'), 1, 16, 2)
            ...
             1n0: 2.538418608559107E-01
             1w0: 2.538418608559108E-01
             1n1: 1.032321840884568E-01
             1w1: 1.059251499158892E-01
             1n2: 3.238743180659171E-02
             1w2: 3.414818600982502E-02
        """
        if prec > 64:
            raise ValueError, "prec (=%s) must be at most 64"%prec
        if prec < 1:
            raise ValueError, "prec (=%s) must be at least 1"%prec
        v = self('-sp %sp%sd%s %s'%(n, prec, d, self._curve_str(E)))
        return self._fix_err(v)

    def modular_degree(self, E):
        """
        Return the modular degree of the elliptic curve E, assuming
        the Stevens conjecture.

        INPUT:
            E -- elliptic curve over Q

        OUTPUT:
            integer -- modular degree

        EXAMPLES:
        We compute the modular degrees of the lowest known conductor
        curves of the first few ranks:

            sage: sympow.modular_degree(EC('11a'))
            1
            sage: sympow.modular_degree(EC('37a'))
            2
            sage: sympow.modular_degree(EC('389a'))
            40
            sage: sympow.modular_degree(EC('5077a'))
            1984
            sage: sympow.modular_degree(EllipticCurve([1, -1, 0, -79, 289]))
            334976
        """
        v = self('%s -moddeg'%self._curve_str(E))
        s = 'Modular Degree is '
        i = v.find(s)
        if i == -1:
            print self._fix_err(v)
            raise RuntimeError, "failed to compute modular degree"
        return sage.rings.all.Integer(v[i+len(s):])

    def analytic_rank(self, E):
        r"""
        Return the analytic rank and leading $L$-value of the elliptic
        curve $E$.

        INPUT:
            E -- elliptic curve over Q

        OUTPUT:
            integer -- analytic rank
            string -- leading coefficient (as string)

        \note{The analytic rank is \emph{not} computed provably
        correctly in general.}

        EXAMPLES:
        We compute the analytic ranks of the lowest known conductor
        curves of the first few ranks:

            sage: sympow.analytic_rank(EllipticCurve('11a'))
            (0, '2.53842e-01')
            sage: sympow.analytic_rank(EllipticCurve('37a'))
            (1, '3.06000e-01')
            sage: sympow.analytic_rank(EllipticCurve('389a'))
            (2, '7.59317e-01')
            sage: sympow.analytic_rank(EllipticCurve('5077a'))
            (3, '1.73185e+00')
            sage: sympow.analytic_rank(EllipticCurve([1, -1, 0, -79, 289]))
            (4, '8.94385e+00')
            sage: sympow.analytic_rank(EllipticCurve([0, 0, 1, -79, 342]))  # long
            (5, '3.02857e+01')
            sage: sympow.analytic_rank(EllipticCurve([1, 1, 0, -2582, 48720]))  # long
            (6, '3.20781e+02')
            sage: sympow.analytic_rank(EllipticCurve([0, 0, 0, -10012, 346900]))  # long
            (7, '1.32517e+03')
        """
        v = self('%s -analrank'%self._curve_str(E))
        s = 'Analytic Rank is '
        i = v.rfind(s)
        if i == -1:
            print self._fix_err(v)
            raise RuntimeError, "failed to compute analytic rank"
        j = v.rfind(':')
        r = sage.rings.all.Integer(v[i+len(s):j])
        i = v.rfind(' ')
        L = v[i+1:]
        return r, L

    def new_data(self, n):
        """
        Pre-compute data files needed for computation of
        n-th symmetric powers.
        """
        print self('-new_data %s'%n)

    def help(self):
        h = """
     sympow('-sp 2p16 -curve "[1,2,3,4,5]"')
will compute L(Sym^2 E,edge) for E=[1,2,3,4,5] to 16 digits of precision
The result
 2n0: 8.370510845377639E-01
 2w0: 8.370510845377586E-01
consists of two calculations using different parameters to test the
functional equation (to see if these are sufficiently close).

     sympow('-sp 3p12d2,4p8 -curve "[1,2,3,4,5]"')
will compute the 0th-2nd derivatives of L(Sym^3 E,center) to 12 digits
and L(Sym^4 E,edge) to 8 digits

Special cases: When a curve has CM, Hecke power can be used instead

     sympow('-sp 7p12d1 -hecke -curve "[0,0,1,-16758,835805]"')

will compute the 0th-1st derivatives of L(Sym^7 psi,special) to 12 digits

Bloch-Kato numbers can be obtained for powers not congruent to 0 mod 4:

     sympow('-sp 2bp16 -curve "[1,2,3,4,5]"')

should return
 2n0: 4.640000000000006E+02
 2w0: 4.639999999999976E+02
which can be seen to be very close to the integer 464.

Modular degrees can be computed with the -moddeg option.

     sympow('-curve "[1,2,3,4,5]" -moddeg')

should return
Modular Degree is 464

Analytic ranks can be computed with the -analrank option.

     sympow('-curve "[1,2,3,4,5]" -analrank')

should return
Analytic Rank is 1 : L'-value 3.51873e+00

and (if the mesh file for the fifth derivative is present)

     sympow('-curve "[0,0,1,-79,342]" -analrank')

should return
Analytic Rank is 5 : leading L-term 3.02857e+01

========================================================================

Adding new symmetric powers:

SYMPOW keeps data for symmetric powers in the directory datafiles
If a pre-computed mesh of inverse Mellin transform values is not
available for a given symmetric power, SYMPOW will fail. The command

     sympow('-new_data 2')

will add the data the 2nd symmetric power, while

     sympow('-new_data 3d2')

will add the data for the 2nd derivative and 3rd symmetric power,

     sympow('-new_data 6d0h')

will add the data for the 0th derivative of the 6th Hecke power, and

     sympow('-new_data 4c')

will add data for the 4th symmetric power for curves with CM
(these need to be done separately for powers divisible by 4).

The mesh files are stored in binary form, and thus endian-ness is a
worry when moving from one platform to another.

To enable modular degree computations, the 2nd symmetric power must
be extant, and analytic rank requires derivatives of the 1st power.
"""
        pager(h)




# An instance
sympow = Sympow()



r"""
Conductor and Reduction Types for Genus 2 Curves

AUTHORS:

- Qing Liu and Henri Cohen (1994-1998): wrote genus2reduction C
  program

- William Stein (2006-03-05): wrote Sage interface to genus2reduction

ACKNOWLEDGMENT: (From Liu's website:) Many thanks to Henri Cohen
who started writing this program. After this program is available,
many people pointed out to me (mathematical as well as programming)
bugs : B. Poonen, E. Schaefer, C. Stahlke, M. Stoll, F. Villegas.
So thanks to all of them. Thanks also go to Ph. Depouilly who help
me to compile the program.

Also Liu has given me explicit permission to include
genus2reduction with Sage and for people to modify the C source
code however they want.
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


import os
from expect import Expect
from sage.structure.sage_object import SageObject

class Genus2reduction_expect(Expect):
    def __init__(self, server=None, server_tmpdir=None, logfile=None):
        Expect.__init__(self,
                        name = 'genus2reduction',
                        prompt = 'enter',
                        command = 'genus2reduction',
                        server = server,
                        server_tmpdir = server_tmpdir,
                        maxread = 10000,
                        restart_on_ctrlc = True,
                        logfile = logfile,
                        verbose_start = False)

    def __getattr__(self, attrname):
        raise AttributeError

class ReductionData(SageObject):
    r"""
    Reduction data for a genus 2 curve.

    How to read ``local_data`` attribute, i.e., if this
    class is R, then the following is the meaning of
    ``R.local_data[p]``.

    For each prime number `p` dividing the discriminant of
    `y^2+Q(x)y=P(x)`, there are two lines.

    The first line contains information about the stable reduction
    after field extension. Here are the meanings of the symbols of
    stable reduction :

    (I) The stable reduction is smooth (i.e. the curve has potentially
    good reduction).

    (II) The stable reduction is an elliptic curve `E` with an
    ordinary double point. `j` mod `p` is the modular
    invariant of `E`.

    (III) The stable reduction is a projective line with two ordinary
    double points.

    (IV) The stable reduction is two projective lines crossing
    transversally at three points.

    (V) The stable reduction is the union of two elliptic curves
    `E_1` and `E_2` intersecting transversally at one
    point. Let `j_1`, `j_2` be their modular
    invariants, then `j_1+j_2` and `j_1 j_2` are
    computed (they are numbers mod `p`).

    (VI) The stable reduction is the union of an elliptic curve
    `E` and a projective line which has an ordinary double
    point. These two components intersect transversally at one point.
    `j` mod `p` is the modular invariant of
    `E`.

    (VII) The stable reduction is as above, but the two components are
    both singular.

    In the cases (I) and (V), the Jacobian `J(C)` has
    potentially good reduction. In the cases (III), (IV) and (VII),
    `J(C)` has potentially multiplicative reduction. In the two
    remaining cases, the (potential) semi-abelian reduction of
    `J(C)` is extension of an elliptic curve (with modular
    invariant `j` mod `p`) by a torus.

    The second line contains three data concerning the reduction at
    `p` without any field extension.


    #. The first symbol describes the REDUCTION AT `p` of
       `C`. We use the symbols of Namikawa-Ueno for the type of
       the reduction (Namikawa, Ueno:"The complete classification of
       fibers in pencils of curves of genus two", Manuscripta Math., vol.
       9, (1973), pages 143-186.) The reduction symbol is followed by the
       corresponding page number (or just an indiction) in the above
       article. The lower index is printed by , for instance, [I2-II-5]
       means [I_2-II-5]. Note that if `K` and `K'` are
       Kodaira symbols for singular fibers of elliptic curves, [K-K'-m]
       and [K'-K-m] are the same type. Finally, [K-K'-1] (not the same as
       [K-K'-1]) is [K'-K-alpha] in the notation of Namikawa-Ueno. The
       figure [2I_0-m] in Namikawa-Ueno, page 159 must be denoted by
       [2I_0-(m+1)].

    #. The second datum is the GROUP OF CONNECTED COMPONENTS (over an
       ALGEBRAIC CLOSURE (!) of `\GF{p}`) of the Neron
       model of J(C). The symbol (n) means the cyclic group with n
       elements. When n=0, (0) is the trivial group (1).
       ``Hn`` is isomorphic to (2)x(2) if n is even and to (4)
       otherwise.

       Note - The set of rational points of `\Phi` can be computed
       using Theorem 1.17 in S. Bosch and Q. Liu "Rational points of the
       group of components of a Neron model", Manuscripta Math. 98 (1999),
       275-293.

    #. Finally, `f` is the exponent of the conductor of
       `J(C)` at `p`.


    .. warning::

       Be careful regarding the formula:

       .. math::

          \text{valuation of the naive minimal discriminant} = f + n - 1 + 11c(X).

       (Q. Liu : "Conducteur et discriminant minimal de courbes de genre
       2", Compositio Math. 94 (1994) 51-79, Theoreme 2) is valid only if
       the residual field is algebraically closed as stated in the paper.
       So this equality does not hold in general over
       `\QQ_p`. The fact is that the minimal discriminant
       may change after unramified extension. One can show however that,
       at worst, the change will stabilize after a quadratic unramified
       extension (Q. Liu : "Modeles entiers de courbes hyperelliptiques
       sur un corps de valuation discrete", Trans. AMS 348 (1996),
       4577-4610, Section 7.2, Proposition 4).
    """
    def __init__(self, raw, P, Q, minimal_equation, minimal_disc,
                 local_data, conductor, prime_to_2_conductor_only):
        self.raw = raw
        self.P = P
        self.Q = Q
        self.minimal_equation = minimal_equation
        self.minimal_disc = minimal_disc
        self.local_data = local_data
        self.conductor = conductor
        self.prime_to_2_conductor_only = prime_to_2_conductor_only

    def _repr_(self):
        if self.prime_to_2_conductor_only:
            ex = ' (away from 2)'
        else:
            ex = ''
        if self.Q == 0:
            yterm = ''
        else:
            yterm = '+ (%s)*y '%self.Q
        s =  'Reduction data about this proper smooth genus 2 curve:\n'
        s += '\ty^2 %s= %s\n'%(yterm, self.P)
        s += 'A Minimal Equation (away from 2):\n\ty^2 = %s\n'%self.minimal_equation
        s += 'Minimal Discriminant (away from 2):  %s\n'%self.minimal_disc
        s += 'Conductor%s: %s\n'%(ex, self.conductor)
        s += 'Local Data:\n%s'%self._local_data_str()
        return s

    def _local_data_str(self):
        s = ''
        D = self.local_data
        K = sorted(D.keys())
        for p in K:
            s += 'p=%s\n%s\n'%(p, D[p])
        s = '\t' + '\n\t'.join(s.split('\n'))
        return s

class Genus2reduction(SageObject):
    r"""
    Conductor and Reduction Types for Genus 2 Curves.

    Use ``R = genus2reduction(Q, P)`` to obtain reduction
    information about the Jacobian of the projective smooth curve
    defined by `y^2 + Q(x)y = P(x)`. Type ``R?``
    for further documentation and a description of how to interpret the
    local reduction data.

    EXAMPLES::

        sage: x = QQ['x'].0
        sage: R = genus2reduction(x^3 - 2*x^2 - 2*x + 1, -5*x^5)
        sage: R.conductor
        1416875
        sage: factor(R.conductor)
        5^4 * 2267

    This means that only the odd part of the conductor is known.

    ::

        sage: R.prime_to_2_conductor_only
        True

    The discriminant is always minimal away from 2, but possibly not at
    2.

    ::

        sage: factor(R.minimal_disc)
        2^3 * 5^5 * 2267

    Printing R summarizes all the information computed about the curve

    ::

        sage: R
        Reduction data about this proper smooth genus 2 curve:
            y^2 + (x^3 - 2*x^2 - 2*x + 1)*y = -5*x^5
        A Minimal Equation (away from 2):
            y^2 = x^6 - 240*x^4 - 2550*x^3 - 11400*x^2 - 24100*x - 19855
        Minimal Discriminant (away from 2):  56675000
        Conductor (away from 2): 1416875
        Local Data:
            p=2
            (potential) stable reduction:  (II), j=1
            p=5
            (potential) stable reduction:  (I)
            reduction at p: [V] page 156, (3), f=4
            p=2267
            (potential) stable reduction:  (II), j=432
            reduction at p: [I{1-0-0}] page 170, (1), f=1

    Here are some examples of curves with modular Jacobians::

        sage: R = genus2reduction(x^3 + x + 1, -2*x^5 - 3*x^2 + 2*x - 2)
        sage: factor(R.conductor)
        23^2
        sage: factor(genus2reduction(x^3 + 1, -x^5 - 3*x^4 + 2*x^2 + 2*x - 2).conductor)
        29^2
        sage: factor(genus2reduction(x^3 + x + 1, x^5 + 2*x^4 + 2*x^3 + x^2 - x - 1).conductor)
        5^6

    EXAMPLE::

        sage: genus2reduction(0, x^6 + 3*x^3 + 63)
        Reduction data about this proper smooth genus 2 curve:
                y^2 = x^6 + 3*x^3 + 63
        A Minimal Equation (away from 2):
                y^2 = x^6 + 3*x^3 + 63
        Minimal Discriminant (away from 2):  10628388316852992
        Conductor (away from 2): 2893401
        Local Data:
                p=2
                (potential) stable reduction:  (V), j1+j2=0, j1*j2=0
                p=3
                (potential) stable reduction:  (I)
                reduction at p: [III{9}] page 184, (3)^2, f=10
                p=7
                (potential) stable reduction:  (V), j1+j2=0, j1*j2=0
                reduction at p: [I{0}-II-0] page 159, (1), f=2

    In the above example, Liu remarks that in fact at `p=2`,
    the reduction is [II-II-0] page 163, (1), `f=8`. So the
    conductor of J(C) is actually `2 \cdot 2893401=5786802`.

    A MODULAR CURVE:

    Consider the modular curve `X_1(13)` defined by an
    equation

    .. math::

                   y^2 + (x^3-x^2-1)y = x^2 - x.



    We have::

        sage: genus2reduction(x^3-x^2-1, x^2 - x)
        Reduction data about this proper smooth genus 2 curve:
                y^2 + (x^3 - x^2 - 1)*y = x^2 - x
        A Minimal Equation (away from 2):
                y^2 = x^6 + 58*x^5 + 1401*x^4 + 18038*x^3 + 130546*x^2 + 503516*x + 808561
        Minimal Discriminant (away from 2):  169
        Conductor: 169
        Local Data:
                p=13
                (potential) stable reduction:  (V), j1+j2=0, j1*j2=0
                reduction at p: [I{0}-II-0] page 159, (1), f=2

    So the curve has good reduction at 2. At `p=13`, the stable
    reduction is union of two elliptic curves, and both of them have 0
    as modular invariant. The reduction at 13 is of type [I_0-II-0]
    (see Namikawa-Ueno, page 159). It is an elliptic curve with a cusp.
    The group of connected components of the Neron model of
    `J(C)` is trivial, and the exponent of the conductor of
    `J(C)` at `13` is `f=2`. The conductor of
    `J(C)` is `13^2`. (Note: It is a theorem of
    Conrad-Edixhoven-Stein that the component group of
    `J(X_1(p))` is trivial for all primes `p`.)
    """
    def __init__(self):
        self.__expect = Genus2reduction_expect()

    def _repr_(self):
        return "Genus 2 reduction program"

    def console(self):
        genus2reduction_console()

    def raw(self, Q, P):
        r"""
        Return the raw output of running the
        ``genus2reduction`` program on the hyperelliptic curve
        `y^2 + Q(x)y = P(x)` as a string.

        INPUT:


        -  ``Q`` - something coercible to a univariate
           polynomial over Q.

        -  ``P`` - something coercible to a univariate
           polynomial over Q.


        OUTPUT:


        -  ``string`` - raw output

        -  ``Q`` - what Q was actually input to auxiliary
           genus2reduction program

        -  ``P`` - what P was actually input to auxiliary
           genus2reduction program


        EXAMPLES::

            sage: x = QQ['x'].0
            sage: print genus2reduction.raw(x^3 - 2*x^2 - 2*x + 1, -5*x^5)[0]
            a minimal equation over Z[1/2] is :
            y^2 = x^6-240*x^4-2550*x^3-11400*x^2-24100*x-19855
            <BLANKLINE>
            factorization of the minimal (away from 2) discriminant :
            [2,3;5,5;2267,1]
            <BLANKLINE>
            p=2
            (potential) stable reduction :  (II), j=1
            p=5
            (potential) stable reduction :  (I)
            reduction at p : [V] page 156, (3), f=4
            p=2267
            (potential) stable reduction :  (II), j=432
            reduction at p : [I{1-0-0}] page 170, (1), f=1
            <BLANKLINE>
            the prime to 2 part of the conductor is 1416875
            in factorized form : [2,0;5,4;2267,1]

        Verify that we fix trac 5573::

            sage: genus2reduction(x^3 + x^2 + x,-2*x^5 + 3*x^4 - x^3 - x^2 - 6*x - 2)
            Reduction data about this proper smooth genus 2 curve:
            y^2 + (x^3 + x^2 + x)*y = -2*x^5 + 3*x^4 - x^3 - x^2 - 6*x - 2
            ...
        """
        from sage.rings.all import QQ
        R = QQ['x']
        P = R(P)
        Q = R(Q)
        if P.degree() > 6:
            raise ValueError, "P (=%s) must have degree at most 6"%P
        if Q.degree() >=4:
            raise ValueError, "Q (=%s) must have degree at most 3"%Q

        E = self.__expect
        try:
            E.eval(str(Q).replace(' ',''))
            s = E.eval(str(P).replace(' ',''))
        except RuntimeError:
            # If something goes wrong genus2reduction often goes into
            # a bad state, and quitting it fixes things, since next
            # time it is used, it is started cleanly.  See trac 5573.
            E.quit()
            raise ValueError, "error in input; possibly singular curve? (Q=%s, P=%s)"%(Q,P)
        i = s.find('a minimal')
        j = s.rfind(']')
        return s[i:j+2], Q, P

    def __call__(self, Q, P):
        from sage.rings.all import ZZ, QQ
        from sage.misc.all import sage_eval

        s, Q, P = self.raw(Q, P)
        raw = s

        if 'the prime to 2 part of the conductor' in s:
            prime_to_2_conductor_only = True
        else:
            prime_to_2_conductor_only = False
        x = QQ['x'].gen(0)
        i = s.find('y^2 = ') + len('y^2 = ')
        j = i + s[i:].find('\n')
        minimal_equation = sage_eval(s[i:j], locals={'x':x})


        s = s[j+1:]
        i = s.find('[')
        j = s.find(']')
        minimal_disc = ZZ(eval(s[i+1:j].replace(',','**').replace(';','*')))

        phrase = 'the conductor is '
        j = s.find(phrase)
        assert j != -1
        k = s[j:].find('\n')
        prime_to_2_conductor = ZZ(s[j+len(phrase):j+k])

        local_data = {}
        while True:
            i = s.find('p=')
            if i == -1:
                break
            j = s[i+2:].find('p=')
            if j == -1:
                j = s.find('\n  \n')
                assert j != -1
            else:
                j = j + (i+2)
            k = i + s[i:].find('\n')
            p = ZZ(s[i+2:k])
            data = s[k+1:j].strip().replace(' : ',': ')
            local_data[p] = data
            s = s[j:]

        return ReductionData(raw, P, Q, minimal_equation, minimal_disc, local_data,
                             prime_to_2_conductor, prime_to_2_conductor_only)


    def __reduce__(self):
        return _reduce_load_Genus2reduction, tuple([])

# An instance
genus2reduction = Genus2reduction()

def _reduce_load_genus2reduction():
    return genus2reduction

import os
def genus2reduction_console():
    os.system('genus2reduction')


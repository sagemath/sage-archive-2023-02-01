r"""
Mathematical constants

The following standard mathematical constants are defined in \sage,
along with support for coercing them into GAP, GP/PARI, KASH, Maxima,
Mathematica, Maple, Octave, and Singular:

    sage: pi
    Pi
    sage: e             # base of the natural logarithm
    e
    sage: NaN           # Not a number
    NaN
    sage: golden_ratio
    golden_ratio
    sage: log2          # natural logarithm of the real number 2
    log2
    sage: euler_gamma   # Euler's gamma constant
    euler_gamma
    sage: catalan       # the Catalon constant
    K
    sage: khinchin      # Khinchin's constant
    khinchin

Suppose for coercion into the various systems means that if, e.g.,
you want to create $\pi$ in Maxima and Singular, you don't have
to figure out the special notation for each system.  You just
type the following:

    sage: pi.str()
    '3.1415926535897931'
    sage: maxima(pi)
    %pi
    sage: singular(pi)
    3.14159265358979323
    sage: gap(pi)
    "3.1415926535897931"
    sage: gp(pi)
    3.141592653589793238462643383     # 32-bit
    3.1415926535897932384626433832795028842   # 64-bit
    sage: pari(pi)
    3.141592653589793238462643383     # 32-bit
    3.1415926535897932384626433832795028842   # 64-bit
    sage: kash(pi)                    # optional
    3.14159265358979323846264338328
    sage: mathematica(pi)             # optional
    Pi
    sage: maple(pi)                   # optional
    Pi
    sage: octave(pi)                  # optional
    3.14159

Arithmetic operations with constants also yield constants, which
can be coerced into other systems or evaluated.
    sage: a = pi + e*4/5; a
    (Pi + ((e*4)/5))
    sage: maxima(a)
    %pi + 4*%e/5
    sage: a.str(15)      # 15 *bits* of precision
    '5.31616'
    sage: gp(a)
    5.316218116357029426750873360            # 32-bit
    5.3162181163570294267508733603616328824  # 64-bit
    sage: mathematica(a)                     # optional
             4 E
            --- + Pi
              5

EXAMPLES: Decimal expansions of constants

We can obtain floating point approximations to each of these constants
by coercing into the real field with given precision.  For example, to
200 decimal places we have the following:

    sage: R = RealField(200); R
    Real Field with 200 bits of precision

    sage: R(pi)
    3.1415926535897932384626433832795028841971693993751058209749445

    sage: R(e)
    2.7182818284590452353602874713526624977572470936999595749669679

    sage: R(NaN)
    NaN

    sage: R(golden_ratio)
    1.6180339887498948482045868343656381177203091798057628621354484

    sage: R(log2)
    0.69314718055994530941723212145817656807550013436025525412067998

    sage: R(euler_gamma)
    0.57721566490153286060651209008240243104215933593992359880576723

    sage: R(catalan)
    0.91596559417721901505460351493238411077414937428167213426649835

    sage: R(khinchin)
    2.6854520010653064453097148354817956938203822939944629530511514


EXAMPLES: Arithmetic with constants

    sage: pp = pi+pi; pp
    (Pi + Pi)
    sage: R(pp)
    6.2831853071795864769252867665590057683943387987502116419498890

    sage: s = (1 + e^pi);s
    (1 + (e^Pi))
    sage: R(s)
    24.140692632779269005729086367948547380266106242600211993445043
    sage: R(s-1)
    23.140692632779269005729086367948547380266106242600211993445043

    sage: l = (1-log2)/(1+log2);l
    ((1 - log2)/(1 + log2))
    sage: R(l)
    0.18123221829928249948761381864650311423330609774776013488055837

    sage: pim = maxima(pi)
    sage: maxima.eval('fpprec : 100')
    '100'
    sage: pim.bfloat()
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068B0



AUTHORS:
    -- Alex Clemesha  <aclemesh@ucsd.edu>, 2006-01-15
    -- William Stein
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#       Copyright (C) 2006 Alex Clemesha <aclemesh@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# TODO: in all.py only import what is really needed from this file, not *.
# (this is for me, William, to do.)

# TODO: See
#
#   http://numbers.computation.free.fr/Constants/Miscellaneous/Records.html
#
# for a bunch of constants, many related to number theory.  Often they
# aren't given on that site, but maybe you can track them down and
# hardcode them up to 1000 digits, say, and for more precision give an
# error message? I've hardcode khinchin's below, as an example.  For
# *many* purposes knowing up to 1000 digits is more than enough.  In
# some cases computing much more is *very* difficult.  I did some
# timing experiments and I think Khinchin is hardcoded in mathematica
# to 5000 digits, since anything above that takes forever, but it is
# instant.  So I'm not against hardcoding either.

import math

from sage.structure.all import RingElement, Element_cmp_
from sage.rings.all import CommutativeRing, RealField
import sage.interfaces.all
import sage.rings.all
import operator
from sage.libs.pari.all import pari

# Default real field used for coercion.

RR = sage.rings.all.RealField(53)

class ConstantRing_class(CommutativeRing):
    def __init__(self):
        self._default_precision = 53  # default bits of precision

    def _repr_(self):
        return "Ring of Mathematical Constants"

    def set_precision(self, prec):
        """
        Change the precision of the default real field used for coercing constants.

        EXAMPLES:
            sage: ConstantRing.set_precision(200)
            sage: pi.str()
            '3.1415926535897932384626433832795028841971693993751058209749445'
            sage: gap(pi)
            "3.1415926535897932384626433832795028841971693993751058209749445"
            sage: ConstantRing.set_precision(53)
            sage: pi.str()
            '3.1415926535897931'

        Note that there seems to be no real numbers in GAP, which is why that
        coercion returns a GAP string.

        The precision is only used when there is no way to describe
        the constant to higher precision (or exactly) to the
        interface:
             sage: maxima(pi)
             %pi

        If we coerce into a GP/PARI session, then the resulting number will
        have precision the precision of that session irregardless of the
        default precision of ConstantRing.  Note that GP/PARI precision is set in
        decimal digits, whereas in SAGE precisions are always in binary.

            sage: gp = Gp()  # new session so as not to affect overwrite global gp session
            sage: gp(pi)
            3.141592653589793238462643383             # 32-bit
            3.1415926535897932384626433832795028842   # 64-bit
            sage: gp.eval('\\p 100')
            '   realprecision = 105 significant digits (100 digits displayed)'      # 32-bit
            '   realprecision = 115 significant digits (100 digits displayed)'      # 64-bit
            sage: gp(pi)
            3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
            sage: ConstantRing.set_precision(5)   # has no affect on gp(pi)
            sage: gp(pi)
            3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
        """
        prec = int(prec)
        self._default_precision = prec
        global RR
        RR = sage.rings.all.RealField(prec)

    def default_precision(self):
        """
        Get the precision of the default real field used for coercing constants.
        """
        return self._default_precision

    def __call__(self, x):
        if isinstance(x, Constant):
            return x
        return Constant_gen(x)

    def characteristic(self):
        return sage.rings.all.Integer(0)

ConstantRing = ConstantRing_class()

class Constant(Element_cmp_, RingElement):
    def __init__(self, conversions={}):
        self.__conversions = conversions
        RingElement.__init__(self, ConstantRing)

    def str(self, bits=None):
        """
        Return a string representation of self as a decimal number
        to at least the given bits of precision.
        """
        if not bits is None:
            R = sage.rings.all.RealField(bits)
        else:
            R = RR
        return str(self._mpfr_(R))

    def _gap_(self, gap):
        try:
            return gap(self.__conversions['gap'])
        except KeyError:
            return gap('"%s"'%self.str())  # no reals in GAP!

    def _gp_(self, gp):
        try:
            return gp(self.__conversions['gp'])
        except KeyError:
            return gp(self.str())

    def _kash_(self, kash):
        try:
            return kash(self.__conversions['kash'])
        except KeyError:
            return kash(self.str())

    def _maxima_(self, maxima):
        try:
            return maxima(self.__conversions['maxima'])
        except KeyError:
            return maxima(self.str())

    def _mathematica_(self, mathematica):
        try:
            return mathematica(self.__conversions['mathematica'])
        except KeyError:
            return mathematica(self.str())

    def _maple_(self, maple):
        try:
            return maple(self.__conversions['maple'])
        except KeyError:
            return maple(self.str())

    def _octave_(self, octave):
        try:
            return octave(self.__conversions['octave'])
        except KeyError:
            return octave(self.str())

    def _pari_(self):
        try:
            return pari(self.__conversions['pari'])
        except KeyError:
            return pari(self.str())

    def _singular_(self, singular):
        try:
	    singular.lib('general')   # Needed for number_e, number_pi
            return singular(self.__conversions['singular'])
        except KeyError:
            return singular(self.str())

    def _mpfr_(self, R):
        raise NotImplementedError, "computation of this constant not implemented yet."

    # The following adds formal arithmetic support for constants
    def _add_(self, right):
        return Constant_arith(self, right, operator.add)

    def _sub_(self, right):
        return Constant_arith(self, right, operator.sub)

    def _mul_(self, right):
        return Constant_arith(self, right, operator.mul)

    def _div_(self, right):
        return Constant_arith(self, right, operator.div)

    def __pow__(self, right):
        right = self.parent()(right)
        return Constant_arith(self, right, operator.pow)

    def _cmp_(self, right):
        """
        EXAMPLES:
            sage: s = e + pi
            sage: s == 0
            False
            sage: t = e^2  +pi
            sage: s == t
            False
            sage: s == s
            True
            sage: t == t
            True
            sage: s < t
            True
            sage: t < s
            False
        """
        if self is right:
            return 0
        R = RealField()
        c = cmp(R(self), R(right))
        if c: return c
        raise NotImplementedError, "equality testing for general mathematical constants not yet implemented (can't prove equality)."

class Pi(Constant):
    """
    The ratio of a circle's circumference to its diameter.

    EXAMPLES:
        sage: pi
        Pi
        sage: float(pi)
        3.1415926535897931
        sage: gp(pi)
        3.141592653589793238462643383            # 32-bit
        3.1415926535897932384626433832795028842  # 64-bit
        sage: R(pi)
        3.1415926535897931
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(pi)
        3.1415926535897932384626433832795028841971693993751058209749445
        sage: pp = pi+pi; pp
        (Pi + Pi)
        sage: R(pp)
        6.2831853071795864769252867665590057683943387987502116419498890
        sage: maxima(pi)
        %pi
        sage: maxima(pi).float()
        3.141592653589793
    """
    def __init__(self):
        Constant.__init__(self,
            {'maxima':'%pi','gp':'Pi','kash':'PI','mathematica':'Pi',
             'matlab':'pi','maple':'Pi','octave':'pi','pari':'Pi'})

    def _repr_(self):
        return "Pi"

    def _latex_(self):
        return "\\pi"

    def __float__(self):
        return math.pi

    def _mpfr_(self, R):
        return R.pi()

    def _singular_(self, singular):
        singular.lib('general')
        return singular('number_pi(%s)'%int((ConstantRing._default_precision/3)+1))

pi = Pi()


class E(Constant):
    """
    The base of the natural logarithm.

    EXAMPLES:
        sage: R(e)
        2.7182818284590451
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(e)
        2.7182818284590452353602874713526624977572470936999595749669679
        sage: em = 1 + e^(1-e); em
        (1 + (e^(1 - e)))
        sage: R(em)
        1.1793740787340171819619895873183164984596816017589156131573701
        sage: maxima(e).float()
        2.718281828459045
        sage: t = mathematica(e)               # optional
        sage: t                                # optional
        E
        sage: float(t)                         # optional
        2.7182818284590451
    """
    def __init__(self):
        Constant.__init__(self,
            {'maxima':'%e','gp':'exp(1)','kash':'E','pari':'exp(1)',
             'mathematica':'E','maple':'exp(1)','octave':'e',
	     'singular':'number_e(20)'}) #singular precision?

    def _repr_(self):
        return 'e'

    def __float__(self):
        return math.e

    def _mpfr_(self, R):
        return R(1).exp()

e = E()
ee = e



class NotANumber(Constant):
    """
    Not a Number
    """
    def __init__(self):
        Constant.__init__(self,
	    {'matlab':'NaN'})

    def _repr_(self):
        return 'NaN'

    def _mpfr_(self,R):
        return R('NaN') #??? nan in mpfr: void mpfr_set_nan (mpfr_t x)

NaN = NotANumber()

class GoldenRatio(Constant):
    """
    The number (1+sqrt(5))/2

    EXAMPLES:
        sage: gr = golden_ratio
        sage: R(gr)
        1.6180339887498949
        sage: R = RealField(200)
        sage: R(gr)
        1.6180339887498948482045868343656381177203091798057628621354484
        sage: grm = maxima(golden_ratio);grm
        (sqrt(5) + 1)/2
        sage: grm + grm
        sqrt(5) + 1
        sage: float(grm + grm)
        3.2360679774997898
    """
    def __init__(self):
        # TODO: Here you should put a string that *symbolically* evaluates
        # to the golden ratio in each system, if possible.  E.g., in
        # mathematica I think Sqrt[5] means sqrt(5) algebraically.
	#
	#no sqrt function in singular?
        Constant.__init__(self,{'mathematica':'N[(1+Sqrt[5])/2]','gp':'(1+sqrt(5))/2',
				'maple':'(1+sqrt(5))/2','maxima':'(1+sqrt(5))/2',
				'pari':'(1+sqrt(5))/2','octave':'(1+sqrt(5))/2',
				'kash':'(1+Sqrt(5))/2'})
    def _repr_(self):
        return 'golden_ratio'

    def _latex_(self): #Should we do this
        return '\\phi'

    def __float__(self):
        return float(0.5)*(float(1.0)+math.sqrt(float(5.0)))

    def _mpfr_(self,R):  #is this OK for _mpfr_ ?
	return (R(1)+R(5).sqrt())*R(0.5)

golden_ratio = GoldenRatio()

class Log2(Constant):
    """
    The natural logarithm of the real number 2.

    EXAMPLES:
        sage: log2
        log2
        sage: float(log2)
        0.69314718055994529
        sage: R(log2)
        0.69314718055994529
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(log2)
        0.69314718055994530941723212145817656807550013436025525412067998
        sage: l = (1-log2)/(1+log2);l
        ((1 - log2)/(1 + log2))
        sage: R(l)
        0.18123221829928249948761381864650311423330609774776013488055837
        sage: maxima(log2)
        log(2)
        sage: maxima(log2).float()
        .6931471805599453
        sage: gp(log2)
        0.6931471805599453094172321215             # 32-bit
        0.69314718055994530941723212145817656808   # 64-bit
    """
    def __init__(self):
        # TODO: Here you should put a string that *symbolically* evaluates
        # to sqrt(2) in each system, if possible.  E.g., in
        # mathematica I think Sqrt[2] means sqrt(2) algebraically.
	#
	#no log function in singular?
         Constant.__init__(self,{'mathematica':'N[Log[2]]','kash':'Log(2)',
				'maple':'log(2)','maxima':'log(2)','gp':'log(2)',
				'pari':'log(2)','octave':'log(2)'})

    def _repr_(self):
        return 'log2'

    def __float__(self):
        return math.log(2) #is this OK?

    def _mpfr_(self,R):
        return R.log2()

log2 = Log2()

class EulerGamma(Constant):
    """
    The limiting difference between the harmonic series and the natural logarithm.

    EXAMPLES:
        sage: R = RealField()
        sage: R(euler_gamma)
	0.57721566490153287
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(euler_gamma)
        0.57721566490153286060651209008240243104215933593992359880576723
        sage: eg = euler_gamma + euler_gamma;eg
        (euler_gamma + euler_gamma)
        sage: R(eg)
        1.1544313298030657212130241801648048620843186718798471976115345
    """
    def __init__(self):
        Constant.__init__(self,
	    {'kash':'EulerGamma(R)','maple':'gamma',
             'mathematica':'EulerGamma','pari':'Euler'})

    def _repr_(self):
        return 'euler_gamma'

    def _latex_(self):
        return '\\gamma'

    def _mpfr_(self,R):
        return R.euler_constant()

euler_gamma = EulerGamma()

class Catalan(Constant):
    """
    A number appaering in combinatorics defined as the Dirichlet beta
    function evaluated at the number 2.

    EXAMPLES:
    """
    def __init__(self):
        Constant.__init__(self,
             {'mathematica':'Catalan','kash':'Catalan(R)', #kash: R is default prec
	      'maple':'Catalan'})

    def _repr_(self):
        return 'K'

    def _mpfr_(self, R):
        return R.catalan_constant()


catalan = Catalan()

class Khinchin(Constant):
    """
    The geometric mean of the continued fraction expansion
    of any (almost any) real number.

    EXAMPLES:
        sage: float(khinchin)
        2.6854520010653062
        sage: khinchin.str(100)
        '2.6854520010653064453097148354799'
        sage: m = mathematica(khinchin); m             # optional
        Khinchin
        sage: m.N(200)                                 # optional
        2.6854520010653064453097148354817956938203822939944629530511523455572188595371520028011411749318476979951534659052880900828976777164109630517925334832596683818523154213321194996260393285220448194096181
    """
    def __init__(self):
        Constant.__init__(self,
             {'mathematica':'Khinchin'}) #Khinchin is only implemented in Mathematica

        # digits come from http://pi.lacim.uqam.ca/piDATA/khintchine.txt
        self.__value = "2.6854520010653064453097148354817956938203822939944629530511523455572188595371520028011411749318476979951534659052880900828976777164109630517925334832596683818523154213321194996260393285220448194096180686641664289308477880620360737053501033672633577289049904270702723451702625237023545810686318501032374655803775026442524852869468234189949157306618987207994137235500057935736698933950879021244642075289741459147693018449050601793499385225470404203377985639831015709022233910000220772509651332460444439191691460859682348212832462282927101269069741823484776754573489862542033926623518620867781366509696583146995271837448054012195366666049648269890827548115254721177330319675947383719393578106059230401890711349624673706841221794681074060891827669566711716683740590473936880953450489997047176390451343232377151032196515038246988883248709353994696082647818120566349467125784366645797409778483662049777748682765697087163192938512899314199518611673792654620563505951385713761697126872299805327673278710513763"
        self.__bits = len(self.__value)*3-1   # underestimate

    def _repr_(self):
        return 'khinchin'

    def _mpfr_(self, R):
        if R.precision() < self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Khinchin's constant only available up to %s bits"%self.__bits

    def __float__(self):
        return 2.685452001065306445309714835481795693820

khinchin  = Khinchin()


#################################################################
#
# Support for arithmetic with constants.
#
#################################################################

symbols = {operator.add:' + ', operator.sub:' - ',
           operator.mul:'*', operator.div:'/',
           operator.pow:'^'}

class Constant_arith(Constant):
    """
    Generic arithmetic with constants is supported.

    EXAMPLES:
        sage: s = (pi + pi) * e + e
        sage: s
        (((Pi + Pi)*e) + e)
        sage: R(s)
        19.797750273806177
        sage: maxima(s)
        2*%e*%pi + %e

        sage: t = e^2 + pi + 2/3; t
        (((e^2) + Pi) + 2/3)
        sage: R(t)
        11.197315419187108
        sage: maxima(t)
        %pi + %e^2 + 2/3
        sage: t^e
        ((((e^2) + Pi) + 2/3)^e)
        sage: R(t^e)
        710.86524768885772
    """
    def __init__(self, x, y, op):
        if not isinstance(x, Constant) or not isinstance(y, Constant):
            raise TypeError
        Constant.__init__(self, {})
        self.__x = x
        self.__y = y
        self.__op = op  # a binary operator, e.g., operator.sub

    def _repr_(self):
        return '(%s%s%s)'%(self.__x, symbols[self.__op], self.__y)

    def _gap_(self, gap):
        return self.__op(self.__x._gap_(gap), self.__y._gap_(gap))

    def _gp_(self, gp):
        return self.__op(self.__x._gp_(gp), self.__y._gp_(gp))

    def _kash_(self, kash):
        return self.__op(self.__x._kash_(kash), self.__y._kash_(kash))

    def _maple_(self, maple):
        return self.__op(self.__x._maple_(maple), self.__y._maple_(maple))

    def _mathematica_(self, mathematica):
        return self.__op(self.__x._mathematica_(mathematica), self.__y._mathematica_(mathematica))

    def _maxima_(self, maxima):
        return self.__op(self.__x._maxima_(maxima), self.__y._maxima_(maxima))

    def _octave_(self, octave):
        return self.__op(self.__x._octave_(octave), self.__y._octave_(octave))

    def _pari_(self):
        return self.__op(self.__x._pari_(), self.__y._pari_())

    def _singular_(self, singular):
        return self.__op(self.__x._singular_(singular), self.__y._singular_(singular))

    def _mpfr_(self, R):
        return self.__op(self.__x._mpfr_(R), self.__y._mpfr_(R))


class Constant_gen(Constant):
    """
    Constant defined by a generic SAGE object.  This makes possible
    symbolic expressions like the following:

    EXAMPLES:
        sage: a = pi/2 + e
        sage: a
        ((Pi/2) + e)
        sage: maxima(a)
        %pi/2 + %e
        sage: RR(a)
        4.2890781552539412
        sage: RealField(200)(a)
        4.2890781552539418545916091629924139398558317933875124854544427

        sage: b = e + 5/7
        sage: maxima(b)
        %e + 5/7
        sage: RR(b)
        3.4325675427447595
    """
    def __init__(self, x):
        Constant.__init__(self)
        self.__x = x

    def _repr_(self):
        return str(self.__x)

    def _mpfr_(self, R):
        return R(self.__x)

    def str(self, bits=None):
        if bits is None:
            return str(self.__x)
        else:
            R = sage.rings.all.RealField(53)
            return str(R(self))

    def _gap_(self, gap):
        return gap(self.__x)

    def _gp_(self, gp):
        return gp(self.__x)

    def _kash_(self, kash):
        return kash(self.__x)

    def _maxima_(self, maxima):
        return maxima(self.__x)

    def _mathematica_(self, mathematica):
        return mathematica(self.__x)

    def _maple_(self, maple):
        return maple(self.__x)

    def _octave_(self, octave):
        return octave(self.__x)

    def _pari_(self):
        return pari(self.__x)

    def _singular_(self, singular):
        return singular(self.__x)


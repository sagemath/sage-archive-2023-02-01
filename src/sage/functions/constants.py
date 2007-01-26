r"""
Mathematical constants

The following standard mathematical constants are defined in \sage,
along with support for coercing them into GAP, GP/PARI, KASH, Maxima,
Mathematica, Maple, Octave, and Singular:

    sage: pi
    pi
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

Support for coercion into the various systems means that if, e.g.,
you want to create $\pi$ in Maxima and Singular, you don't have
to figure out the special notation for each system.  You just
type the following:

    sage: pi.str()
    '3.14159265358979'
    sage: maxima(pi)
    %pi
    sage: singular(pi)
    3.14159265358979
    sage: gap(pi)
    "3.14159265358979"
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
    (pi + ((e*4)/5))
    sage: maxima(a)
    %pi + 4*%e/5
    sage: a.str(15)      # 15 *bits* of precision
    '5.316'
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
    3.1415926535897932384626433832795028841971693993751058209749

    sage: R(e)
    2.7182818284590452353602874713526624977572470936999595749669

    sage: R(NaN)
    NaN

    sage: R(golden_ratio)
    1.6180339887498948482045868343656381177203091798057628621354

    sage: R(log2)
    0.69314718055994530941723212145817656807550013436025525412067

    sage: R(euler_gamma)
    0.57721566490153286060651209008240243104215933593992359880576

    sage: R(catalan)
    0.91596559417721901505460351493238411077414937428167213426649

    sage: R(khinchin)
    2.6854520010653064453097148354817956938203822939944629530511


EXAMPLES: Arithmetic with constants

    sage: pp = pi+pi; pp
    (pi + pi)
    sage: R(pp)
    6.2831853071795864769252867665590057683943387987502116419498

    sage: s = (1 + e^pi);s
    (1 + (e^pi))
    sage: R(s)
    24.140692632779269005729086367948547380266106242600211993445
    sage: R(s-1)
    23.140692632779269005729086367948547380266106242600211993445

    sage: l = (1-log2)/(1+log2);l
    ((1 - log2)/(1 + log2))
    sage: R(l)
    0.18123221829928249948761381864650311423330609774776013488055

    sage: pim = maxima(pi)
    sage: maxima.eval('fpprec : 100')
    '100'
    sage: pim.bfloat()
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068b0



AUTHORS:
    -- Alex Clemesha  <aclemesh@ucsd.edu>, 2006-01-15
    -- William Stein
    -- Alex Clemesha \& William Stein (2006-02-20): added new constants; removed todos
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 Alex Clemesha <aclemesh@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import math, operator

from sage.rings.all import CommutativeRing, RealField, Integer, RingElement
import sage.interfaces.all
import sage.rings.all
from sage.libs.pari.all import pari
from sage.misc.latex import latex


# Default real field used for coercion.

from functions import Function_gen, Function_arith, Function, FunctionRing_class

######################
# Ring of Constants
######################

class ConstantRing_class(FunctionRing_class):
    def _repr_(self):
        return "Ring of Real Mathematical Constants"

    def __cmp__(self, right):
        if isinstance(right, ConstantRing_class):
            return 0
        return -1

    def __call__(self, x):
        try:
            return self._coerce_(x)
        except TypeError:
            return Constant_gen(x)

    def _coerce_impl(self, x):
        if isinstance(x, (sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            return Constant_gen(x)
        raise TypeError, 'no canonical coercion of element into self.'

ConstantRing = ConstantRing_class()


######################
# Constant functions
######################

class Constant(Function):
    def __init__(self, conversions={}):
        self._conversions = conversions
        RingElement.__init__(self, ConstantRing)

    def _neg_(self):
        return -Integer(1)*self

    def __call__(self, x):
        return self

    def floor(self):
        return Integer(int(float(self)))

    # The following adds formal arithmetic support for generic constant
    def _add_(self, right):
        return Constant_arith(self, right, operator.add)

    def _sub_(self, right):
        return Constant_arith(self, right, operator.sub)

    def _mul_(self, right):
        return Constant_arith(self, right, operator.mul)

    def _div_(self, right):
        return Constant_arith(self, right, operator.div)

    def __pow__(self, right):
        try:
            right = self.parent()._coerce_(right)
        except TypeError:
            raise TypeError, "computation of %s^%s not defined"%(self, right)
        return Constant_arith(self, right, operator.pow)

    def _interface_is_cached_(self):
        """
        Return False, since coercion of functions to interfaces
        is not cached.

        We do not cache coercions of functions to interfaces, since
        the precision of the interface may change.

        EXAMPLES:
            sage: gp(pi)
            3.141592653589793238462643383              # 32-bit
            3.1415926535897932384626433832795028842    # 64-bit
            sage: old_prec = gp.set_precision(100)
            sage: gp(pi)
            3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
            sage: _ = gp.set_precision(old_prec)
            sage: gp(pi)
            3.141592653589793238462643383              # 32-bit
            3.1415926535897932384626433832795028842    # 64-bit
        """
        return False


class Constant_gen(Constant, Function_gen):
    def __init__(self, x):
        Function_gen.__init__(self, x)
        Constant.__init__(self)

    def __call__(self, x):
        return self.obj()

    def __repr__(self):
        return Function_gen._repr_(self)

class Constant_arith(Constant, Function_arith):
    def __init__(self, x, y, op):
        Function_arith.__init__(self, x, y, op)
        Constant.__init__(self)


class Pi(Constant):
    """
    The ratio of a circle's circumference to its diameter.

    EXAMPLES:
        sage: pi
        pi
        sage: float(pi)
        3.1415926535897931
        sage: gp(pi)
        3.141592653589793238462643383            # 32-bit
        3.1415926535897932384626433832795028842  # 64-bit
        sage: RR(pi)
        3.14159265358979
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(pi)
        3.1415926535897932384626433832795028841971693993751058209749
        sage: pp = pi+pi; pp
        (pi + pi)
        sage: R(pp)
        6.2831853071795864769252867665590057683943387987502116419498
        sage: maxima(pi)
        %pi
        sage: maxima(pi).float()
        3.141592653589793
    """
    def __init__(self):
        Constant.__init__(self,
            {'axiom':'%pi',
             'maxima':'%pi','gp':'Pi','kash':'PI','mathematica':'Pi',
             'matlab':'pi','maple':'Pi','octave':'pi','pari':'Pi'})

    def _repr_(self):
        return "pi"

    def _latex_(self):
        return "\\pi"

    def _mathml_(self):
        return "<mi>&pi;</mi>"

    def __float__(self):
        return math.pi

    def _mpfr_(self, R):
        return R.pi()

    def _real_double_(self):
        return sage.rings.all.RD.pi()

    def __abs__(self):
        if self.str()[0] != '-':
            return self
        return -self

    def floor(self):
        return Integer(3)

    # This just gives a string in singular anyways, and it's
    # *REALLY* slow!
    #def _singular_(self, singular):
    #    singular.lib('general')
    #    return singular('number_pi(%s)'%int((ConstantRing._default_precision/3)+1))

pi = Pi()


class E(Constant):
    """
    The base of the natural logarithm.

    EXAMPLES:
        sage: RR(e)
        2.71828182845904
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(e)
        2.7182818284590452353602874713526624977572470936999595749669
        sage: em = 1 + e^(1-e); em
        (1 + (e^(1 - e)))
        sage: R(em)
        1.1793740787340171819619895873183164984596816017589156131573
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
            {'axiom':'%e',
             'maxima':'%e',
             'gp':'exp(1)',
             'kash':'E',
             'pari':'exp(1)',
             'mathematica':'E',
             'maple':'exp(1)',
             'octave':'e'})

    def _repr_(self):
        return 'e'

    def _latex_(self):
        return 'e'

    def __float__(self):
        return math.e

    def _mpfr_(self, R):
        return R(1).exp()

    def floor(self):
        return Integer(2)

    def _real_double_(self):
        return sage.rings.all.RD.e()

    # This just gives a string in singular anyways, and it's
    # *REALLY* slow!
    #def _singular_(self, singular):
    #    singular.lib('general')
    #    return singular('number_e(%s)'%int((ConstantRing._default_precision/3)+1))

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

    def _real_double_(self):
        return sage.rings.all.RD.nan()

NaN = NotANumber()

class GoldenRatio(Constant):
    """
    The number (1+sqrt(5))/2

    EXAMPLES:
        sage: gr = golden_ratio
        sage: RR(gr)
        1.61803398874989
        sage: R = RealField(200)
        sage: R(gr)
        1.6180339887498948482045868343656381177203091798057628621354
        sage: grm = maxima(golden_ratio);grm
        (sqrt(5) + 1)/2
        sage: grm + grm
        sqrt(5) + 1
        sage: float(grm + grm)
        3.2360679774997898
    """
    def __init__(self):
        Constant.__init__(self,{'mathematica':'N[(1+Sqrt[5])/2]','gp':'(1+sqrt(5))/2',
				'maple':'(1+sqrt(5))/2','maxima':'(1+sqrt(5))/2',
				'pari':'(1+sqrt(5))/2','octave':'(1+sqrt(5))/2',
				'kash':'(1+Sqrt(5))/2'})
    def _repr_(self):
        return 'golden_ratio'

    def _latex_(self):
        return '\\phi'

    def __float__(self):
        return float(0.5)*(float(1.0)+math.sqrt(float(5.0)))

    def _real_double_(self):
        """
        EXAMPLES:
            sage: RDF(golden_ratio)
            1.61803398875
        """
        return sage.rings.all.RDF(1.61803398874989484820458)

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
        sage: RR(log2)
        0.693147180559945
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(log2)
        0.69314718055994530941723212145817656807550013436025525412067
        sage: l = (1-log2)/(1+log2);l
        ((1 - log2)/(1 + log2))
        sage: R(l)
        0.18123221829928249948761381864650311423330609774776013488055
        sage: maxima(log2)
        log(2)
        sage: maxima(log2).float()
        .6931471805599453
        sage: gp(log2)
        0.6931471805599453094172321215             # 32-bit
        0.69314718055994530941723212145817656807   # 64-bit
    """
    def __init__(self):
        # TODO: Here you should put a string that *symbolically* evaluates
        # to sqrt(2) in each system, if possible.  E.g., in
        # mathematica I think Sqrt[2] means sqrt(2) algebraically.
        Constant.__init__(self,{'mathematica':'N[Log[2]]','kash':'Log(2)',
				'maple':'log(2)','maxima':'log(2)','gp':'log(2)',
				'pari':'log(2)','octave':'log(2)'})

    def _repr_(self):
        return 'log2'

    def _latex_(self):
        return '\\log(2)'

    def __float__(self):
        return math.log(2)

    def _real_double_(self):
        """
        EXAMPLES:
            sage: RDF(log2)
            0.69314718056
        """
        return sage.rings.all.RD.log2()

    def _mpfr_(self,R):
        return R.log2()

    def floor(self):
        return Integer(0)

log2 = Log2()

class EulerGamma(Constant):
    """
    The limiting difference between the harmonic series and the natural logarithm.

    EXAMPLES:
        sage: R = RealField()
        sage: R(euler_gamma)
        0.577215664901532
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(euler_gamma)
        0.57721566490153286060651209008240243104215933593992359880576
        sage: eg = euler_gamma + euler_gamma;eg
        (euler_gamma + euler_gamma)
        sage: R(eg)
        1.1544313298030657212130241801648048620843186718798471976115
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

    def _real_double_(self):
        """
        EXAMPLES:
            sage: RDF(euler_gamma)
            0.577215664902
        """
        return sage.rings.all.RD.euler()

    def floor(self):
        return Integer(0)

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

    def _latex_(self):
        return 'K'

    def _mpfr_(self, R):
        return R.catalan_constant()

    def _real_double_(self):
        """
        EXAMPLES:
        We coerce to the real double field:
            sage: RDF(catalan)
            0.915965594177
        """
        return sage.rings.all.RDF(0.91596559417721901505460351493252)

    def __float__(self):
        """
        EXAMPLES:
            sage: float(catalan)
            0.91596559417721901
        """
        return 0.91596559417721901505460351493252

    def floor(self):
        return Integer(0)

catalan = Catalan()

class Khinchin(Constant):
    """
    The geometric mean of the continued fraction expansion
    of any (almost any) real number.

    EXAMPLES:
        sage: float(khinchin)
        2.6854520010653062
        sage: khinchin.str(100)
        '2.6854520010653064453097148354'
        sage: m = mathematica(khinchin); m             # optional
        Khinchin
        sage: m.N(200)                                 # optional
             2.685452001065306445309714835481795693820382293994462953051152345557     # 32-bit
        >    218859537152002801141174931847697995153465905288090082897677716410963051 # 32-bit
        >    7925334832596683818523154213321194996260393285220448194096181            # 32-bit
             2.685452001065306445309714835481795693820382293994462953051152345557     # 64-bit
        >    218859537152002801141174931847697995153465905288090082897677716410963051 # 64-bit
        >    7925334832596683818523154213321194996260393285220448194096181            # 64-bit
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
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Khinchin's constant only available up to %s bits"%self.__bits

    def _real_double_(self):
        """
        EXAMPLES:
            sage: RDF(khinchin)
            2.68545200107
        """
	return sage.rings.all.RDF(2.685452001065306445309714835481795693820)


    def __float__(self):
        return 2.685452001065306445309714835481795693820

    def floor(self):
        return Integer(2)

khinchin  = Khinchin()


class TwinPrime(Constant):
    r"""
    The Twin Primes constant is defined as $\prod 1 - 1/(p-1)^2$
    for primes $p > 2$.

    EXAMPLES:
	sage: float(twinprime)
        0.66016181584686962
        sage: R=RealField(200);R
        Real Field with 200 bits of precision
        sage: R(twinprime)
        0.66016181584686957392781211001455577843262336028473341331944
    """
    def __init__(self):
        Constant.__init__(self,{}) #Twin prime is not implemented in any other algebra systems.

        #digits come from http://www.gn-50uma.de/alula/essays/Moree/Moree-details.en.shtml

        self.__value = "0.660161815846869573927812110014555778432623360284733413319448423335405642304495277143760031413839867911779005226693304002965847755123366227747165713213986968741097620630214153735434853131596097803669932135255299767199302474590593101082978291553834469297505205916657133653611991532464281301172462306379341060056466676584434063501649322723528968010934966475600478812357962789459842433655749375581854814173628678098705969498703841243363386589311969079150040573717814371081810615401233104810577794415613125444598860988997585328984038108718035525261719887112136382808782349722374224097142697441764455225265548994829771790977784043757891956590649994567062907828608828395990394287082529070521554595671723599449769037800675978761690802426600295711092099633708272559284672129858001148697941855401824639887493941711828528382365997050328725708087980662201068630474305201992394282014311102297265141514194258422242375342296879836738796224286600285358098482833679152235700192585875285961205994728621007171131607980572"

        self.__bits = len(self.__value)*3-1   # underestimate

    def floor(self):
        return Integer(0)

    def _repr_(self):
        return 'twinprime'

    def _mpfr_(self, R):
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Twin Prime constant only available up to %s bits"%self.__bits

    def _real_double_(self):
        """
        EXAMPLES:
            sage: RDF(twinprime)
            0.660161815847
        """
	return sage.rings.all.RDF(0.660161815846869573927812110014555778432)

    def __float__(self):
        """
        EXAMPLES:
            sage: float(twinprime)
            0.66016181584686962
        """
	return 0.660161815846869573927812110014555778432

twinprime = TwinPrime()


class Merten(Constant):
    """
    The Merten constant is related to the Twin Primes constant
    and appears in Merten's second theorem.

    EXAMPLES:
        sage: float(merten)
        0.26149721284764277
        sage: R=RealField(200);R
        Real Field with 200 bits of precision
        sage: R(merten)
        0.26149721284764278375542683860869585905156664826119920619206
    """
    def __init__(self):
        Constant.__init__(self,{}) #Merten's constant is not implemented in any other algebra systems.

        # digits come from Sloane's tables at http://www.research.att.com/~njas/sequences/table?a=77761&fmt=0

        self.__value = "0.261497212847642783755426838608695859051566648261199206192064213924924510897368209714142631434246651051617"

        self.__bits = len(self.__value)*3-1   # underestimate

    def floor(self):
        return Integer(0)

    def _repr_(self):
        return 'merten'

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(1000)(merten)
            Traceback (most recent call last):
            ...
            NotImplementedError: Merten's constant only available up to 320 bits
            sage: RealField(320)(merten)
            0.261497212847642783755426838608695859051566648261199206192064213924924510897368209714142631434246
        """
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Merten's constant only available up to %s bits"%self.__bits

    def _real_double_(self):
        """
        EXAMPLES:
            sage: RDF(merten)
            0.261497212848
        """
        return sage.rings.all.RDF(0.261497212847642783755426838608695859051)

    def __float__(self):
        """
        EXAMPLES:
            sage: float(merten)
            0.26149721284764277
        """
	return 0.261497212847642783755426838608695859051

merten = Merten()

class Brun(Constant):
    """
    Brun's constant is the sum of reciprocals of odd twin primes.

    It is not known to very high precision; calculating the number
    using twin primes up to $10^16$ (Sebah 2002) gives the number
    $1.9021605831040$.

    EXAMPLES:
	sage: float(brun)
        1.902160583104
        sage: R = RealField(41); R
        Real Field with 41 bits of precision
        sage: R(brun)
        1.90216058310
    """
    def __init__(self):
        Constant.__init__(self,{}) #Brun's constant is not implemented in any other algebra systems.

        # digits come from Sloane's tables at http://www.research.att.com/~njas/sequences/table?a=65421&fmt=0

        self.__value = "1.902160583104"

        self.__bits = len(self.__value)*3-1 # bits  -- todo: make more intelligent in a general function!!!

    def floor(self):
        return Integer(1)

    def _repr_(self):
        return 'brun'

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(53)(brun)
            Traceback (most recent call last):
            ...
            NotImplementedError: Brun's constant only available up to 41 bits
            sage: RealField(41)(brun)
            1.90216058310
        """
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Brun's constant only available up to %s bits"%self.__bits

    def _real_double_(self):
        """
        EXAMPLES:
            sage: RDF(brun)
            1.9021605831
        """
        return sage.rings.all.RDF(1.9021605831040)

    def __float__(self):
        """
        EXAMPLES:
            sage: float(brun)
            1.902160583104
        """
	return 1.9021605831040

brun=Brun()

class UniversalPolynomialElement(Constant):
    """
    A universal indeterminate.

    EXAMPLES:
        sage: x
        x
        sage: x.parent()
        Univariate Polynomial Ring in x over Rational Field
    """
    def __init__(self, name):
        self._name = name
        Constant.__init__(self,
            {'axiom':name,
             'maxima':name,
             'gp':name,'kash':name,
             'mathematica':name,
             'matlab':name,
             'maple':name,
             'octave':name,
             'pari':name})

    def _repr_(self):
        return self._name

    def _latex_(self):
        return self._name

    def _mathml_(self):
        return "<mi>%s</mi>"%self._name

    def __float__(self):
        raise TypeError

    def _mpfr_(self, R):
        raise TypeError

    def _real_double_(self):
        raise TypeError

    def __abs__(self):
        raise TypeError

    def floor(self):
        raise TypeError

x = UniversalPolynomialElement('x')

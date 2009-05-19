"""
Mathematical constants

The following standard mathematical constants are defined in Sage,
along with support for coercing them into GAP, GP/PARI, KASH,
Maxima, Mathematica, Maple, Octave, and Singular::

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
    catalan
    sage: khinchin      # Khinchin's constant
    khinchin
    sage: twinprime
    twinprime
    sage: merten
    merten
    sage: brun
    brun

Support for coercion into the various systems means that if, e.g.,
you want to create `\pi` in Maxima and Singular, you don't
have to figure out the special notation for each system. You just
type the following::

    sage: maxima(pi)
    %pi
    sage: singular(pi)
    pi
    sage: gap(pi)
    pi
    sage: gp(pi)
    3.141592653589793238462643383     # 32-bit
    3.1415926535897932384626433832795028842   # 64-bit
    sage: pari(pi)
    3.14159265358979
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

::

    sage: a = pi + e*4/5; a
    pi + 4/5*e
    sage: maxima(a)
    %pi+4*%e/5
    sage: RealField(15)(a)           # 15 *bits* of precision
    5.316
    sage: gp(a)
    5.316218116357029426750873360              # 32-bit
    5.3162181163570294267508733603616328824    # 64-bit
    sage: print mathematica(a)                     # optional
     4 E
     --- + Pi
      5

EXAMPLES: Decimal expansions of constants

We can obtain floating point approximations to each of these
constants by coercing into the real field with given precision. For
example, to 200 decimal places we have the following::

    sage: R = RealField(200); R
    Real Field with 200 bits of precision

::

    sage: R(pi)
    3.1415926535897932384626433832795028841971693993751058209749

::

    sage: R(e)
    2.7182818284590452353602874713526624977572470936999595749670

::

    sage: R(NaN)
    NaN

::

    sage: R(golden_ratio)
    1.6180339887498948482045868343656381177203091798057628621354

::

    sage: R(log2)
    0.69314718055994530941723212145817656807550013436025525412068

::

    sage: R(euler_gamma)
    0.57721566490153286060651209008240243104215933593992359880577

::

    sage: R(catalan)
    0.91596559417721901505460351493238411077414937428167213426650

::

    sage: R(khinchin)
    2.6854520010653064453097148354817956938203822939944629530512

EXAMPLES: Arithmetic with constants

::

    sage: f = I*(e+1); f
    I*e + I
    sage: f^2
    (I*e + I)^2
    sage: _.expand()
    -2*e - e^2 - 1

::

    sage: pp = pi+pi; pp
    2*pi
    sage: R(pp)
    6.2831853071795864769252867665590057683943387987502116419499

::

    sage: s = (1 + e^pi); s
    e^pi + 1
    sage: R(s)
    24.140692632779269005729086367948547380266106242600211993445
    sage: R(s-1)
    23.140692632779269005729086367948547380266106242600211993445

::

    sage: l = (1-log2)/(1+log2); l
    -(log2 - 1)/(log2 + 1)
    sage: R(l)
    0.18123221829928249948761381864650311423330609774776013488056

::

    sage: pim = maxima(pi)
    sage: maxima.eval('fpprec : 100')
    '100'
    sage: pim.bfloat()
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068b0

AUTHORS:

- Alex Clemesha (2006-01-15)

- William Stein

- Alex Clemesha, William Stein (2006-02-20): added new constants;
  removed todos

- Didier Deshommes (2007-03-27): added constants from RQDF (deprecated)

TESTS:

Coercing the sum of a bunch of the constants to many different
floating point rings::

    sage: a = pi + e + golden_ratio + log2 + euler_gamma + catalan + khinchin + twinprime + merten; a
    pi + euler_gamma + catalan + golden_ratio + log2 + khinchin + twinprime + merten + e
    sage: parent(a)
    Symbolic Ring
    sage: RR(a)
    13.2713479401972
    sage: RealField(212)(a)
    13.2713479401972493100988191995758139408711068200030748178329712
    sage: RealField(230)(a)
    13.271347940197249310098819199575813940871106820003074817832971189555
    sage: CC(a)
    13.2713479401972
    sage: CDF(a)
    13.2713479402
    sage: ComplexField(230)(a)
    13.271347940197249310098819199575813940871106820003074817832971189555
    sage: RDF(a)
    13.2713479402
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#                     2009 Mike Hansen <mhansen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################
import math
from functools import partial
from sage.rings.infinity import (infinity, minus_infinity,
                                 unsigned_infinity)

constants_table = {}
constants_name_table = {}
constants_name_table[repr(infinity)] = infinity
constants_name_table[repr(unsigned_infinity)] = unsigned_infinity
constants_name_table[repr(minus_infinity)] = minus_infinity

def unpickle_Constant(class_name, name, conversions, latex, mathml, domain):
    """
    EXAMPLES::

        sage: from sage.symbolic.constants import unpickle_Constant
        sage: a = unpickle_Constant('Constant', 'a', {}, 'aa', '', 'positive')
        sage: a.domain()
        'positive'
        sage: latex(a)
        aa

    Note that if the name already appears in the
    ``constants_name_table``, then that will be returned instead of
    constructing a new object::

        sage: pi = unpickle_Constant('Pi', 'pi', None, None, None, None)
        sage: pi._maxima_init_()
        '%pi'
    """
    if name in constants_name_table:
        return constants_name_table[name]
    if class_name == "Constant":
        return Constant(name, conversions=conversions, latex=latex,
                        mathml=mathml, domain=domain)
    else:
        cls = globals()[class_name]
        return cls(name=name)

class Constant(object):
    def __init__(self, name, conversions=None, latex=None, mathml="",
                 domain='complex'):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: loads(dumps(p))
            p
        """
        self._conversions = conversions if conversions is not None else {}
        self._latex = latex if latex is not None else name
        self._mathml = mathml
        self._name = name
        self._domain = domain

        for system, value in self._conversions.items():
            setattr(self, "_%s_"%system, partial(self._generic_interface, value))
            setattr(self, "_%s_init_"%system, partial(self._generic_interface_init, value))

        from sage.symbolic.constants_c import PynacConstant
        self._pynac = PynacConstant(self._name, self._latex, self._domain)
        self._serial = self._pynac.serial()
        constants_table[self._serial] = self
        constants_name_table[self._name] = self

        from sage.symbolic.pynac import register_symbol
        register_symbol(self.expression(), self._conversions)

    def __reduce__(self):
        """
        Adds support for pickling constants.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: p.__reduce__()
            (<function unpickle_Constant at 0x...>,
             ('Constant', 'p', {}, 'p', '', 'complex'))
            sage: loads(dumps(p))
            p

            sage: pi.pyobject().__reduce__()
            (<function unpickle_Constant at 0x...>,
             ('Pi',
              'pi',
              ...,
              '\\pi',
              '<mi>&pi;</mi>',
              'positive'))
            sage: loads(dumps(pi.pyobject()))
            pi

        """
        return (unpickle_Constant, (self.__class__.__name__, self._name,
                                    self._conversions, self._latex,
                                    self._mathml, self._domain))


    def domain(self):
        """
        Returns the domain of this constant.  This is either positive,
        real, or complex, and is used by Pynac to make inferences
        about expressions containing this constant.

        EXAMPLES::

            sage: p = pi.pyobject(); p
            pi
            sage: type(_)
            <class 'sage.symbolic.constants.Pi'>
            sage: p.domain()
            'positive'
        """
        return self._domain

    def expression(self):
        """
        Returns an expression for this constant.

        EXAMPLES::

            sage: a = pi.pyobject()
            sage: pi2 = a.expression()
            sage: pi2
            pi
            sage: pi2 + 2
            pi + 2
            sage: pi - pi2
            0
        """
        return self._pynac.expression()

    def name(self):
        """
        Returns the name of this constant.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c')
            sage: c.name()
            'c'
        """
        return self._name

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c')
            sage: c
            c
        """
        return self._name

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c', latex=r'\xi')
            sage: latex(c)
            \xi
        """
        return self._latex

    def _mathml_(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: c = Constant('c', mathml=r'<mi>c</mi>')
            sage: mathml(c)
            <mi>c</mi>
        """
        return self._mathml

    def _generic_interface(self, value, I):
        """
        This is a helper method used in defining the ``_X_`` methods
        where ``X`` is the name of some interface.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p', conversions=dict(maxima='%pi'))
            sage: p._maxima_(maxima)
            %pi

        The above ``_maxima_`` is constructed like ``m`` below::

            sage: from functools import partial
            sage: m = partial(p._generic_interface, '%pi')
            sage: m(maxima)
            %pi

        """
        return I(value)

    def _generic_interface_init(self, value):
        """
        This is a helper method used in defining the ``_X_init_`` methods
        where ``X`` is the name of some interface.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p', conversions=dict(maxima='%pi'))
            sage: p._maxima_init_()
            '%pi'

        The above ``_maxima_init_`` is constructed like ``mi`` below::

            sage: from functools import partial
            sage: mi = partial(p._generic_interface_init, '%pi')
            sage: mi()
            '%pi'

        """
        return value

    def _interface_(self, I):
        """
        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p', conversions=dict(maxima='%pi'))
            sage: p._interface_(maxima)
            %pi
        """
        try:
            s = self._conversions[I.name()]
            return I(s)
        except KeyError:
            pass

        try:
            return getattr(self, "_%s_"%(I.name()))(I)
        except AttributeError:
            pass

        raise NotImplementedError

    def _gap_(self, gap):
        """
        Returns the constant as a string in GAP. Since GAP does not
        have floating point numbers, we simply return the constant as
        a string.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: gap(p)
            p
        """
        return gap('"%s"'%self)

    def _singular_(self, singular):
        """
        Returns the constant as a string in Singular. Since Singular
        does not have floating point numbers, we simply return the
        constant as a string.

        EXAMPLES::

            sage: from sage.symbolic.constants import Constant
            sage: p = Constant('p')
            sage: singular(p)
            p
        """
        return singular('"%s"'%self)


class Pi(Constant):
    def __init__(self, name="pi"):
        """
        TESTS::

            sage: pi._latex_()
            '\\pi'
            sage: latex(pi)
            \pi
            sage: mathml(pi)
            <mi>&pi;</mi>

        """
        conversions = dict(axiom='%pi', maxima='%pi', gp='Pi', kash='PI',
                           mathematica='Pi', matlab='pi', maple='pi',
                           octave='pi', pari='Pi', pynac='Pi')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r"\pi", mathml="<mi>&pi;</mi>",
                          domain='positive')

    def __float__(self):
        """
        EXAMPLES::

            sage: float(pi)
            3.1415926535897931
        """
        return math.pi

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: pi._mpfr_(RealField(100))
            3.1415926535897932384626433833
        """
        return R.pi()

    def _real_double_(self, R):
        """
         EXAMPLES::

             sage: pi._real_double_(RDF)
             3.14159265359
         """
        return R.pi()

    def _sympy_(self):
        """
        Converts pi to sympy pi.

        EXAMPLES::

            sage: import sympy
            sage: sympy.pi == pi # indirect doctest
            True
        """
        import sympy
        return sympy.pi

pi = Pi().expression()

class I_class(Constant):
    def __init__(self, name="I"):
        """
        The formal square root of -1.

        .. warning::

           Note that calling :meth:`pyobject` on ``I`` from within
           Sage does not return an instance of this class.  Instead,
           it returns a wrapper around a number field element.

        EXAMPLES::

            sage: I
            I
            sage: I^2
            -1

        Note that conversions to real fields will give TypeErrors::

            sage: float(I)
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex to float; use abs(z)
            sage: gp(I)
            I
            sage: RR(I)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert x (='1.00000000000000*I') to real number.

        We can convert to complex fields::

            sage: C = ComplexField(200); C
            Complex Field with 200 bits of precision
            sage: C(I)
            1.0000000000000000000000000000000000000000000000000000000000*I
            sage: I._complex_mpfr_field_(ComplexField(53))
            1.00000000000000*I

            sage: I._complex_double_(CDF)
            1.0*I
            sage: CDF(I)
            1.0*I

            sage: z = I + I; z
            2*I
            sage: C(z)
            2.0000000000000000000000000000000000000000000000000000000000*I
            sage: 1e8*I
            1.00000000000000e8*I

            sage: complex(I)
            1j

            sage: QQbar(I)
            1*I

            sage: abs(I)
            1

            sage: I.minpoly()
            x^2 + 1
            sage: maxima(2*I)
            2*%i

        TESTS::

            sage: repr(I)
            'I'
            sage: latex(I)
            I
        """
        conversions = dict(axiom='%i', maxima='%i', gp='I',
                           mathematica='I', matlab='i',maple='I',
                           octave='i', pari='I')
        Constant.__init__(self, name, conversions=conversions,
                          latex='i', mathml="<mi>&i;</mi>")


    def expression(self, constant=False):
        """
        Returns an Expression for I.  If *constnat* is True, then it
        returns a wrapper around a Pynac constant.  If *constant* is
        False, then it returns a wrapper around a NumberFieldElement.

        EXAMPLES::

            sage: from sage.symbolic.constants import I_class
            sage: a = I_class()
            sage: I_constant = a.expression(constant=True)
            sage: type(I_constant.pyobject())
            <class 'sage.symbolic.constants.I_class'>
            sage: I_nf = a.expression()
            sage: type(I_nf.pyobject())
            <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
        """
        if constant:
            return Constant.expression(self)
        else:
            from sage.symbolic.pynac import I
            return I

a = I_class()
I_constant = a.expression(constant=True)
i = I = a.expression(constant=False)

class E(Constant):
    def __init__(self, name="e"):
        """
        The base of the natural logarithm.

        EXAMPLES::

            sage: RR(e)
            2.71828182845905
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(e)
            2.7182818284590452353602874713526624977572470936999595749670
            sage: em = 1 + e^(1-e); em
            e^(-e + 1) + 1
            sage: R(em)
            1.1793740787340171819619895873183164984596816017589156131574
            sage: maxima(e).float()
            2.718281828459045
            sage: t = mathematica(e)               # optional
            sage: t                                # optional
            E
            sage: float(t)                         # optional
            2.7182818284590451

            sage: loads(dumps(e))
            e
        """
        conversions = dict(axiom='%e', maxima='%e', gp='exp(1)',
                           kash='E', pari='exp(1)', mathematica='E',
                           maple='exp(1)', octave='e')
        Constant.__init__(self, name, conversions=conversions,
                          latex='e', domain='real')

    def expression(self):
        """
        .. note::

           For e, we don't return a wrapper around a Pynac constant.
           Instead, we return exp(1) so that Pynac can perform appropiate.

        EXAMPLES::

            sage: e + 2
            e + 2
            sage: e.operator()
            exp
            sage: e.operands()
            [1]
        """
        from sage.symbolic.ring import SR
        return SR(1).exp()

    def __float__(self):
        """
        EXAMPLES::

            sage: float(e)
            2.7182818284590451
            sage: e.__float__()
            2.7182818284590451
        """
        return math.e

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: e._mpfr_(RealField(100))
            2.7182818284590452353602874714
        """
        return R(1).exp()

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: e._real_double_(RDF)
            2.71828182846
        """
        return R(1).exp()

    def _sympy_(self):
        """
        Converts e to sympy E.

        EXAMPLES::

            sage: import sympy
            sage: sympy.E == e # indirect doctest
            True
        """
        import sympy
        return sympy.E

e = E().expression()

class NotANumber(Constant):
    """
    Not a Number
    """
    def __init__(self, name="NaN"):
        """
        EXAMPLES::

            sage: loads(dumps(NaN))
            NaN
        """
        conversions=dict(matlab='NaN')
        Constant.__init__(self, name, conversions=conversions)

    def _mpfr_(self,R):
        """
        EXAMPLES::

            sage: NaN._mpfr_(RealField(53))
            NaN
            sage: type(_)
            <type 'sage.rings.real_mpfr.RealNumber'>
        """
        return R('NaN') #??? nan in mpfr: void mpfr_set_nan (mpfr_t x)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(NaN)
            NaN
        """
        return R.NaN()

    def _sympy_(self):
        """
        Converts NaN to SymPy NaN.

        EXAMPLES::

            sage: import sympy
            sage: sympy.nan == NaN # indirect doctest
            True
        """
        import sympy
        return sympy.nan

NaN = NotANumber().expression()

class GoldenRatio(Constant):
    """
    The number (1+sqrt(5))/2

    EXAMPLES::

        sage: gr = golden_ratio
        sage: RR(gr)
        1.61803398874989
        sage: R = RealField(200)
        sage: R(gr)
        1.6180339887498948482045868343656381177203091798057628621354
        sage: grm = maxima(golden_ratio);grm
        (sqrt(5)+1)/2
        sage: grm + grm
        sqrt(5)+1
        sage: float(grm + grm)
        3.2360679774997898
    """
    def __init__(self, name='golden_ratio'):
        """
        EXAMPLES::

            sage: loads(dumps(golden_ratio))
            golden_ratio
        """
        conversions = dict(mathematica='N[(1+Sqrt[5])/2]', gp='(1+sqrt(5))/2',
                           maple='(1+sqrt(5))/2', maxima='(1+sqrt(5))/2',
                           pari='(1+sqrt(5))/2', octave='(1+sqrt(5))/2',
                           kash='(1+Sqrt(5))/2')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r'\phi', domain='positive')

    def minpoly(self, bits=None, degree=None, epsilon=0):
        """
        EXAMPLES::

            sage: golden_ratio.minpoly()
            x^2 - x - 1
        """
        from sage.rings.all import QQ
        x = QQ['x'].gen(0)
        return x**2 - x - 1

    def __float__(self):
        """
        EXAMPLES::

            sage: float(golden_ratio)
            1.6180339887498949
            sage: golden_ratio.__float__()
            1.6180339887498949
        """
        return float(0.5)*(float(1.0)+math.sqrt(float(5.0)))

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(golden_ratio)
            1.61803398875
        """
        return R('1.61803398874989484820458')


    def _mpfr_(self,R):
        """
        EXAMPLES::

            sage: golden_ratio._mpfr_(RealField(100))
            1.6180339887498948482045868344
            sage: RealField(100)(golden_ratio)
            1.6180339887498948482045868344
        """
        return (R(1)+R(5).sqrt())/R(2)

    def _algebraic_(self, field):
        """
        EXAMPLES::

            sage: golden_ratio._algebraic_(QQbar)
            1.618033988749895?
            sage: QQbar(golden_ratio)
            1.618033988749895?
        """
        import sage.rings.qqbar
        return field(sage.rings.qqbar.get_AA_golden_ratio())

    def _sympy_(self):
        """
        Converts golden_ratio to SymPy GoldenRation.

        EXAMPLES::

            sage: import sympy
            sage: sympy.GoldenRatio == golden_ratio # indirect doctest
            True
        """
        import sympy
        return sympy.GoldenRatio

golden_ratio = GoldenRatio().expression()

class Log2(Constant):
    """
    The natural logarithm of the real number 2.

    EXAMPLES::

        sage: log2
        log2
        sage: float(log2)
        0.69314718055994529
        sage: RR(log2)
        0.693147180559945
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(log2)
        0.69314718055994530941723212145817656807550013436025525412068
        sage: l = (1-log2)/(1+log2); l
        -(log2 - 1)/(log2 + 1)
        sage: R(l)
        0.18123221829928249948761381864650311423330609774776013488056
        sage: maxima(log2)
        log(2)
        sage: maxima(log2).float()
        .6931471805599453
        sage: gp(log2)
        0.6931471805599453094172321215             # 32-bit
        0.69314718055994530941723212145817656807   # 64-bit
    """
    def __init__(self, name='log2'):
        """
        EXAMPLES::

            sage: loads(dumps(log2))
            log2
        """
        conversions = dict(mathematica='N[Log[2]]', kash='Log(2)',
                           maple='log(2)', maxima='log(2)', gp='log(2)',
                           pari='log(2)', octave='log(2)')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r'\log(2)', domain='positive')

    def __float__(self):
        """
        EXAMPLES::

            sage: float(log2)
            0.69314718055994529
            sage: log2.__float__()
            0.69314718055994529
        """
        return math.log(2)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(log2)
            0.69314718056
        """
        return R.log2()


    def _mpfr_(self,R):
        """
        EXAMPLES::

            sage: RealField(100)(log2)
            0.69314718055994530941723212146
            sage: log2._mpfr_(RealField(100))
            0.69314718055994530941723212146
        """
        return R.log2()

log2 = Log2().expression()

class EulerGamma(Constant):
    """
    The limiting difference between the harmonic series and the natural
    logarithm.

    EXAMPLES::

        sage: R = RealField()
        sage: R(euler_gamma)
        0.577215664901533
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(euler_gamma)
        0.57721566490153286060651209008240243104215933593992359880577
        sage: eg = euler_gamma + euler_gamma; eg
        2*euler_gamma
        sage: R(eg)
        1.1544313298030657212130241801648048620843186718798471976115
    """
    def __init__(self, name='euler_gamma'):
        """
        EXAMPLES::

            sage: loads(dumps(euler_gamma))
            euler_gamma
        """
        conversions = dict(kash='EulerGamma(R)', maple='gamma',
                           mathematica='EulerGamma', pari='Euler',
                           maxima='%gamma', pynac='Euler')
        Constant.__init__(self, name, conversions=conversions,
                          latex=r'\gamma', domain='positive')

    def _mpfr_(self,R):
        """
        EXAMPLES::

            sage: RealField(100)(euler_gamma)
            0.57721566490153286060651209008
            sage: euler_gamma._mpfr_(RealField(100))
            0.57721566490153286060651209008
        """
        return R.euler_constant()

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(euler_gamma)
            0.577215664902
        """
        return R.euler_constant()

    def _sympy_(self):
        """
        Converts euler_gamma to SymPy EulerGamma.

        EXAMPLES::

            sage: import sympy
            sage: sympy.EulerGamma == euler_gamma # indirect doctest
            True
        """
        import sympy
        return sympy.EulerGamma

euler_gamma = EulerGamma().expression()

class Catalan(Constant):
    """
    A number appaering in combinatorics defined as the Dirichlet beta
    function evaluated at the number 2.

    EXAMPLES::

        sage: catalan^2 + merten
        merten + catalan^2
    """
    def __init__(self, name='catalan'):
        """
        EXAMPLES::

            sage: loads(dumps(catalan))
            catalan
        """
        #kash: R is default prec
        conversions = dict(mathematica='Catalan', kash='Catalan(R)',
                           maple='Catalan', maxima='catalan',
                           pynac='Catalan')
        Constant.__init__(self, name, conversions=conversions,
                          domain='positive')


    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: RealField(100)(catalan)
            0.91596559417721901505460351493
            sage: catalan._mpfr_(RealField(100))
            0.91596559417721901505460351493
        """
        return R.catalan_constant()

    def _real_double_(self, R):
        """
        EXAMPLES: We coerce to the real double field::

            sage: RDF(catalan)
            0.915965594177
        """
        return R('0.91596559417721901505460351493252')


    def __float__(self):
        """
        EXAMPLES::

            sage: float(catalan)
            0.91596559417721901
        """
        return 0.91596559417721901505460351493252

    def _sympy_(self):
        """
        Converts catalan to SymPy Catalan.

        EXAMPLES::

            sage: import sympy
            sage: sympy.Catalan == catalan # indirect doctest
            True
        """
        import sympy
        return sympy.Catalan

catalan = Catalan().expression()

###############################
# Limited precision constants #
###############################

# NOTE: LimitedPrecisionConstant is an *abstract base class*.  It is
# not instantiated.  It doesn't make sense to pickle and unpickle
# this, so there is no loads(dumps(...)) doctest below.

class LimitedPrecisionConstant(Constant):
    def __init__(self, name, value, **kwds):
        """
        A class for constants that are only known to a limited
        precision.

        EXAMPLES::

            sage: from sage.symbolic.constants import LimitedPrecisionConstant
            sage: a = LimitedPrecisionConstant('a', '1.234567891011121314').expression(); a
            a
            sage: RDF(a)
            1.23456789101
            sage: RealField(200)(a)
            Traceback (most recent call last):
            ...
            NotImplementedError: a is only available up to 59 bits
        """
        self._value = value
        self._bits = len(self._value)*3-1 #FIXME
        Constant.__init__(self, name, **kwds)

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: RealField(100)(khinchin)
            2.6854520010653064453097148355
            sage: RealField(20000)(khinchin)
            Traceback (most recent call last):
            ...
            NotImplementedError: khinchin is only available up to 3005 bits
        """
        if R.precision() <= self._bits:
            return R(self._value)
        raise NotImplementedError, "%s is only available up to %s bits"%(self.name(), self._bits)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(khinchin)
            2.68545200107
        """
        return R(self._value)

    def __float__(self):
        """
        EXAMPLES::

            sage: float(khinchin)
            2.6854520010653062
        """
        return float(self._value)

class Khinchin(LimitedPrecisionConstant):
    """
    The geometric mean of the continued fraction expansion of any
    (almost any) real number.

    EXAMPLES::

        sage: float(khinchin)
        2.6854520010653062
        sage: m = mathematica(khinchin); m             # optional
        Khinchin
        sage: m.N(200)                                 # optional
        2.68545200106530644530971483548179569382038229399446295305115234555721885953715200280114117493184769799515346590528809008289767771641096305179253348325966838185231542133211949962603932852204481940961807          # 32-bit
        2.6854520010653064453097148354817956938203822939944629530511523455572188595371520028011411749318476979951534659052880900828976777164109630517925334832596683818523154213321194996260393285220448194096181                # 64-bit
    """
    def __init__(self, name='khinchin'):
        """
        EXAMPLES::

            sage: loads(dumps(khinchin))
            khinchin
        """
        conversions = dict(maxima='khinchin', mathematica='Khinchin')
        # digits come from http://pi.lacim.uqam.ca/piDATA/khintchine.txt
        value = "2.6854520010653064453097148354817956938203822939944629530511523455572188595371520028011411749318476979951534659052880900828976777164109630517925334832596683818523154213321194996260393285220448194096180686641664289308477880620360737053501033672633577289049904270702723451702625237023545810686318501032374655803775026442524852869468234189949157306618987207994137235500057935736698933950879021244642075289741459147693018449050601793499385225470404203377985639831015709022233910000220772509651332460444439191691460859682348212832462282927101269069741823484776754573489862542033926623518620867781366509696583146995271837448054012195366666049648269890827548115254721177330319675947383719393578106059230401890711349624673706841221794681074060891827669566711716683740590473936880953450489997047176390451343232377151032196515038246988883248709353994696082647818120566349467125784366645797409778483662049777748682765697087163192938512899314199518611673792654620563505951385713761697126872299805327673278710513763"

        LimitedPrecisionConstant.__init__(self, name, value,
                                          conversions=conversions,
                                          domain='positive')

khinchin = Khinchin().expression()

class TwinPrime(LimitedPrecisionConstant):
    r"""
    The Twin Primes constant is defined as
    `\prod 1 - 1/(p-1)^2` for primes `p > 2`.

    EXAMPLES::

        sage: float(twinprime)
        0.66016181584686962
        sage: R=RealField(200);R
        Real Field with 200 bits of precision
        sage: R(twinprime)
        0.66016181584686957392781211001455577843262336028473341331945
    """
    def __init__(self, name='twinprime'):
        """
        EXAMPLES::

            sage: loads(dumps(twinprime))
            twinprime
        """
        conversions = dict(maxima='twinprime')
        #digits come from http://www.gn-50uma.de/alula/essays/Moree/Moree-details.en.shtml

        value = "0.660161815846869573927812110014555778432623360284733413319448423335405642304495277143760031413839867911779005226693304002965847755123366227747165713213986968741097620630214153735434853131596097803669932135255299767199302474590593101082978291553834469297505205916657133653611991532464281301172462306379341060056466676584434063501649322723528968010934966475600478812357962789459842433655749375581854814173628678098705969498703841243363386589311969079150040573717814371081810615401233104810577794415613125444598860988997585328984038108718035525261719887112136382808782349722374224097142697441764455225265548994829771790977784043757891956590649994567062907828608828395990394287082529070521554595671723599449769037800675978761690802426600295711092099633708272559284672129858001148697941855401824639887493941711828528382365997050328725708087980662201068630474305201992394282014311102297265141514194258422242375342296879836738796224286600285358098482833679152235700192585875285961205994728621007171131607980572"

        LimitedPrecisionConstant.__init__(self, name, value,
                                          conversions=conversions,
                                          domain='positive')

twinprime = TwinPrime().expression()


class Merten(LimitedPrecisionConstant):
    """
    The Merten constant is related to the Twin Primes constant and
    appears in Merten's second theorem.

    EXAMPLES::

        sage: float(merten)
        0.26149721284764277
        sage: R=RealField(200);R
        Real Field with 200 bits of precision
        sage: R(merten)
        0.26149721284764278375542683860869585905156664826119920619206
    """
    def __init__(self, name='merten'):
        """
        EXAMPLES::

            sage: loads(dumps(merten))
            merten
        """
        conversions = dict(maxima='merten')

        # digits come from Sloane's tables at http://www.research.att.com/~njas/sequences/table?a=77761&fmt=0
        value = "0.261497212847642783755426838608695859051566648261199206192064213924924510897368209714142631434246651051617"

        LimitedPrecisionConstant.__init__(self, name, value,
                                          conversions=conversions,
                                          domain='positive')

merten = Merten().expression()

class Brun(LimitedPrecisionConstant):
    """
    Brun's constant is the sum of reciprocals of odd twin primes.

    It is not known to very high precision; calculating the number
    using twin primes up to `10^{16}` (Sebah 2002) gives the
    number `1.9021605831040`.

    EXAMPLES::

        sage: float(brun)
        Traceback (most recent call last):
        ...
        NotImplementedError: brun is only available up to 41 bits
        sage: R = RealField(41); R
        Real Field with 41 bits of precision
        sage: R(brun)
        1.90216058310
    """
    def __init__(self, name='brun'):
        """
        EXAMPLES::

            sage: loads(dumps(brun))
            brun
        """
        conversions = dict(maxima='brun')

        # digits come from Sloane's tables at http://www.research.att.com/~njas/sequences/table?a=65421&fmt=0
        value = "1.902160583104"

        LimitedPrecisionConstant.__init__(self, name, value,
                                          conversions=conversions,
                                          domain='positive')

brun = Brun().expression()

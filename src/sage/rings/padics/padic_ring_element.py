

from sage.libs.all import pari, pari_gen, PariError
from sage.structure.element import Element
from sage.misc.latex import latex
from sage.rings.infinity import infinity
import sage.rings.arith
import sage.rings.integer
import sage.rings.integer_mod
import sage.rings.ring_element
import sage.rings.padics.padic_ring
import sage.rings.arith
import sage.rings.rational
from sage.rings.rational_field import frac, QQ
import sage.misc.functional

class pAdicRingGenericElement(sage.rings.commutative_ring_element.CommutativeRingElement):
    def __init__(self, parent):
	commutative_ring_element.CommutativeRingElement.__init__(self,parent)
	self._p = parent.prime()

    def __cmp__(self, other):
	raise NotImplementedError

    def __getitem__(self, n):
        r"""
        Returns the coefficient of $p^n$ in the series expansion of self, as an integer in the range $0$ to $p-1$.

	EXAMPLE:
	    sage: R = pAdicRingCappedRelative(7,4); R.set_print_mode('series'); a = R(1/3); a
		5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
	    sage: a[0]
		5
	    sage: a[1]
		4
	"""

	if n >= self.precision_absolute():
	    raise ValueError, "element not known to enough precision."
	elif n < self.valuation():
	    return 0
	else:
	    return self.list()[n]

    def __getslice__(self, i, j):
	r"""
	Returns a list of the coefficients of the powers of $p$ from $p^i$ to $p^{j-1}$, as integers in the range $0$ to $p-1$.

	EXAMPLE:
	    sage: R = pAdicRingCappedRelative(7,4); R.set_print_mode('series'); a = R(1/3); a
		5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
	    sage: a[0:2]
		[5,4]
	    sage: a[1,4]
		[4,4,4]
	"""
	if j > self.precision_absolute():
	    raise ValueError, "element not known to enough precision."
	head = [0 for w in range(i,self.valuation())]
	return head + self.list()[max(i - self.valuation(), 0) : j - self.valuation()]

    def __invert__(self, prec=infinity):
	r"""
	Returns the multiplicative inverse of self.

	EXAMPLE:
	    sage: R = pAdicRingCappedRelative(7,4); R.set_print_mode('series'); a = R(3); a
		3 + O(7^4)
	    sage: ~a
		5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)

	NOTES:
	The element returned is an element of the fraction field.
	"""
	return self.parent().fraction_field()(self).__invert__()

    def __mod__(self, right):
	raise NotImplementedError

    def __neg__(self):
	raise NotImplementedError

    def __pos__(self):
	return self

    def __pow__(self, right):
	raise NotImplementedError

    def _add_(self, right):
	raise NotImplementedError

    def _div_(self, right):
	r"""
	Returns the quotient of self by right.

	EXAMPLE:
	    sage: R = pAdicRingCappedRelative(7,4); R.set_print_mode('series'); a = R(3); b = R(5); a / b
		2 + 4*7 + 5*7^2 + 2*7^3 + O(7^4)
	"""
	return self * right.__invert__()

    def _floordiv_(self, right):
	raise NotImplementedError

    def _integer_(self):
	raise NotImplementedError

    def _latex_(self):
	return self._repr_(do_latex=True)

    def _mul_(self):
	raise NotImplementedError

    def _pari_init_(self):
        raise NotImplementedError

    def _repr_(self, do_latex=False):
	r"""
	    Prints a string representation of the element.  See set_print_mode for more details.

	    EXAMPLE:
		sage: R = pAdicRingCappedRelative(7,4); R.set_print_mode('val-unit'); a = R(364); a
		    7^1 * 52 + O(7^4)
		sage: R.set_print_mode('integer'); a
		    364 + O(7^4)
		sage: R.set_print_mode('integer-p'); a
		    364 + O(p^4)
		sage: R.set_print_mode('val-unit-p'); a
		    p^1 * 52 + O(p^4)
		sage: R.set_print_mode('series'); a
		    3*7 + 7^3 + O(7^4)
		sage: R.set_print_mode('series-p'); a
		    3*p + p^3 + O(p^4)
	"""
	if self.parent().get_print_mode() == 'val-unit':
	    if do_latex:
		return "%s^{%s} \\cdot %s + O(%s^{%s})"%(self.prime(), self.valuation(), self.unit_part(), self.prime(), self.precision_absolute())
	    else:
		return "%s^%s * %s + O(%s^%s)"%(self.prime(), self.valuation(), self.unit_part(), self.prime(), self.precision_absolute())
	elif self.parent().get_print_mode() == 'integer':
	    if do_latex:
		return "%s + O(%s^{%s})"%(self.lift(), self.prime(), self.precision_absolute())
	    else:
		return "%s + O(%s^%s)"%(self.lift(), self.prime(), self.precision_absolute())
	elif self.parent().get_print_mode() == 'integer-p':
	    if do_latex:
		return "%s + O(p^{%s})"%(self.lift(), self.precision_absolute())
	    else:
		return "%s + O(p^%s)"%(self.lift(), self.precision_absolute())
	elif self.parent().get_print_mode() == 'val-unit-p':
	    if do_latex:
		return "p^{%s} * %s + O(p^{%s})"%(self.valuation(), self.unit_part(), self.precision_absolute())
	    else:
		return "p^%s * %s + O(p^%s)"%(self.valuation(), self.unit_part(), self.precision_absolute())
	else:
	    v = self.lift()
	    p = self.prime()
	    s = ""
	    exp = 0
	    while v != 0:
		coeff = v % p
		if coeff != 0:
		    if exp == 0:
			s += "%s + "%coeff
		    else:
			if self.parent().get_print_mode() == 'series':
			    var = "%s"%p
			else:
			    var = "p"
			if exp != 1:
			    if do_latex:
				var += "^{%s}"%exp
			    else:
				var += "^%s"%exp
			if coeff != 1:
			    if do_latex:
				s += "%s \\cdot %s + "%(coeff, var)
			    else:
				s += "%s*%s + "%(coeff, var)
			else:
			    s += "%s + "%var
		exp += 1
		v = (v - coeff) / p
	    if self.parent().get_print_mode() == 'series':
		s += "O(%s"%(p)
	    else:
		s += "O(p"
	    if self._modulus == 1:
		s += ")"
	    else:
		if do_latex:
		    s += "^{%s})"%self.precision_absolute()
		else:
		    s += "^%s)"%self.precision_absolute()
	    return s

    def _sub_(self, right):
	r"""
	    Returns the difference between self and right.

	    EXAMPLE:
		sage: R = pAdicRingCappedRelative(7,4); R.set_print_mode('series'); a = R(12); b = R(5); a - b
		    7 + O(7^4)
	"""
	return self + (-right)

    def add_bigoh(self, prec):
	raise NotImplementedError

    def additive_order(self, prec):
	r"""
	    Returns the additive order of self, where self is considered to be zero if it is zero modulo $p^{\mbox{prec}}$.
	"""
	if self.is_zero(prec):
	    return integer.Integer(1)
	else:
	    return infinity

    def algdep(self, n):
	"""
        Returns a polynomial of degree at most $n$ which is approximately
        satisfied by this number.  Note that the returned polynomial
        need not be irreducible, and indeed usually won't be if this number
        is a good approximation to an algebraic number of degree less than $n$.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLE:
        sage: R = pAdicRingCappedRelative(3,20); R.set_print_mode('series')
        sage: a = R(7/19); a
        1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
        sage: a.algdep(1)
        19*x - 7
        """
	return sage.rings.arith.algdep(self, n)

    algebraic_dependency = algdep

    def cache_set_precision(self, prec):
	raise NotImplementedError

    def cache_get_precision(self, prec):
	raise NotImplementedError

    def copy(self):
	raise NotImplementedError

    def exp(self):
	raise NotImplementedError

    def exp_artin_hasse(self):
	raise NotImplementedError

    def gamma(self):
	raise NotImplementedError

    def is_integral(self):
	"""
	Returns whether self is a p-adic integer, ie. true.
	"""
	return True

    def is_square(self): #should be overridden for lazy elements
	if self.valuation() == infinity:
	    return True
	elif self.precision_relative() < 1:
	    return True
	elif self.prime() != 2:
	    return (self.valuation() % 2 == 0) and (self.unit_part().reduce(self.prime()).is_square())
	else:
	    return (self.valuation() % 2 == 0) and (self.unit_part().residue(8) == 1)

    def is_unit(self):
	return self.valuation() == 0

    def is_zero(self, prec):
	raise NotImplementedError

    def is_equal_to(self, right, prec):
	raise NotImplementedError

    def lift(self):
	raise NotImplementedError

    def list(self):
	raise NotImplementedError

    def log(self):
	raise NotImplementedError

    def log_artin_hasse(self):
	raise NotImplementedError

    def minimal_polynomial(self, name):
	import sage.rings.polynomial_ring
	R = sage.rings.polynomial_ring.PolynomialRing(self.parent(), name)
	return R.gen() - R(self)

    def multiplicative_order(self, prec): #needs to be rewritten for lazy elements
	if self.valuation() > 0:
	    return infinity
	else:
	    return Mod(self.unit_part(), self.parent().prime_pow(min(prec, self.precision_relative()))).multiplicative_order()

    def norm(self, ground=None):
	if (ground != None) and (ground != self.parent()):
	    raise ValueError, "Ground Field not a subfield"
	else:
	    return self

    def ordp(self):
	return self.valuation()

    def padded_list(self, n):
	raise NotImplementedError

    def precision_absolute(self):
	raise NotImplementedError

    def precision_relative(self):
	raise NotImplementedError

    def prime(self):
	return self._p

    def prime_pow(self, n):
	return self._p ** n

    def rational_reconstruction(self):
	if self.is_zero():
	    return frac(0,1)
	p = self.prime()
	alpha = self._unit().lift()
	m = integer.Integer(p**self.precision_relative())
	r = arith.rational_reconstruction(alpha, m)
	return (frac(p, 1)**self.valuation())*r

    def residue(self, prec):
	raise NotImplementedError

    def sqrt(self):
	return square_root(self)

    def square_root(self):
	raise NotImplementedError

    def trace(self, ground=None):
	if (ground != None) and (ground != self.parent()):
	    raise ValueError, "Ground Field not a subfield"
	else:
	    return self

    def unit_part(self):
	raise NotImplementedError

    def valuation(self):
	raise NotImplementedError

class pAdicRingCappedRelativeElement(pAdicRingGenericElement):
    def __init__(self, parent, x, prec=None, construct=False):
	pAdicRingGenericElement.__init__(self, parent)
	if construct:
	    (self._ordp, self._unit, self._relprec) = x
	    return
	if prec == None or prec > parent.precision_cap():
	    prec = parent.precision_cap()
	if isinstance(x, pAdicRingGenericElement): #this construction doesn't quite work with lazy p-adics
	    if parent.prime() != x.prime():
		raise ValueError, "Cannot coerce between p-adic rings with different primes."
	    self._ordp = x.valuation()
	    if x.precision_relative() > parent.precision_cap():
		self._unit = x.unit_part().residue(x.prime_pow(parent.precision_cap()))
		self._relprec = parent.precision_cap()
	    else:
		self._unit = x.unit_part().residue(x.prime_pow(x.precision_relative()))
		self._relprec = x.precision_relative()
	    return
	if isinstance(x, pAdicFieldCappedRelativeElement):
	    if x.valuation() < 0:
		raise ValueError, "element not a p-adic integer."
	    else:
		if parent.prime() != x.prime():
		    raise ValueError, "Cannot coerce between p-adic rings with different primes."
		self._ordp = x.valuation()
		if x.precision_relative() > parent.precision_cap():
		    self._unit = x.unit_part().residue(x.prime_pow(parent.precision_cap()))
		    self._relprec = parent.precision_cap()
		else:
		    self._unit = x.unit_part().residue(x.prime_pow(x.precision_relative()))
		    self._relprec = x.precision_relative()
		return

	if isinstance(x, pari_gen) and x.type() == "t_PADIC":
	    t = x.lift()
	    big_oh = x.padicprec(self.prime())
	    if t.type() == 't_INT':
		x = int(t)
	    else:
		x = QQ(t)

	    # We now use the code, below, so don't make the next line elif
	if isinstance(x, (int, long)):
	    self._ordp = sage.rings.arith.valuation(x, self.prime())
	elif isinstance(x, (integer.Integer, rational.Rational)):
	    self._ordp = x.valuation(self.prime())
	else:
	    raise TypeError, "unable to compute ordp"
	if self._ordp < 0:
	    raise ValueError, "element not a p-adic integer."
	elif self._ordp == infinity:
	    self._unit = 1
	    self._relprec = parent.precision_cap()
	    return
	x = x // self.prime_pow(self._ordp)
	self._unit = Mod(x, self.prime_pow(parent.precision_cap()))
	self._relprec = parent.precision_cap()
	return

    def __cmp__(self, other):
	if not isinstance(other, pAdicRingGenericElement) or other.parent() != self.parent():
	    return coerce.cmp(self, other)
	m = min(self.precision_absolute(), other.precision_absolute())
	other_ordp = other.valuation()
        raise NotImplementedError

    def __mod__(self, right):
	val = self.valuation()
	rval = right.valuation()
	if rval > val:
	    raise ValueError, "not enough precision to reduce"
	if right == self.prime_pow(rval):
	    return integer_mod.Mod(self.lift(), right)
	else:
	    raise ValueError, "modulus must be a power of p"

    def __neg__(self):
	return pAdicRingCappedRelativeElement(self.parent(), (self.prime(), self._ordp, -self._unit, self._relprec), construct=True)

    def __pow__(self, right):
        raise NotImplementedError

    def _add_(self, right): #assumes both have the same parent (ie same p and same relative cap)
	if self.valuation() == right.valuation():
	    if self._relprec == right._relprec:
		u = self._unit + right._unit
		rprec = self._relprec
	    elif self._relprec < right._relprec:
		u = self._unit + Mod(right._unit, self.prime_pow(self._relprec))
		rprec = self._relprec
	    else:
		u = Mod(self._unit, self.prime_pow(right._relprec)) + right._unit
		rprec = right._relprec
	    ul = u.lift()
	    if ul == 0: #In this case we set the valuation of the sum to be the minimum possible valuation and say that we have zero precision.
		v = min(self._relprec, right._relprec)
		rprec = 0
		u = 1
	    else:
	    	v = ul.valuation(self.prime())
	    	if v > 0:
		    u = Mod(ul // self.prime_pow(v), self.prime_pow(rprec - v))
	            rprec = rprec - v
	else:
	    rprec = min(min(self.valuation() + self.precision_relative(), right.valuation + right.precision_relative()) - min(self.valuation(), right.valuation()), self.parent().precision_cap())
	    uself = self._unit.lift() * self.prime_pow(self.valuation() - min(self.valuation(), right.valuation()))
	    uright = right._unit.lift() * self.prime_pow(right.valuation() - min(self.valuation(), right.valuation()))
	    u = Mod(uself + uright, self.prime_pow(rprec))
	    v = 0
	return pAdicRingCappedRelativeElement(self.parent(), (min(self.valuation(), right.valuation()) + v, u, rprec), construct = True)

    def __floordiv__(self, right):
        raise NotImplementedError

    def _integer_(self):
	return self._unit.lift() * self.prime_pow(self.valuation())

    def _mul_(self):
	rprec = min(self._relprec, right._relprec)
	return pAdicRingCappedRelativeElement(self.parent(), (self.valuation() + right.valuation(), Mod(self._unit, self.prime_pow(rprec)) * Mod(right._unit, self.prime_pow(rprec)), rprec), construct = True)

    def add_bigoh(self, prec):
	rprec = min(self._relprec, prec - self.valuation())
	if rprec <= 0:
	    return pAdicRingCappedRelativeElement(self.parent(), (prec, 1, 0), construct = True) #need to think about whether setting relprec to 0 causes problems
	return pAdicRingCappedRelativeElement(self.parent(), (self.valuation(), Mod(self._unit, self.prime_pow(rprec)), rprec), construct = True)

    def cache_set_precision(self, prec):
	if prec <= self.precision_relative():
	    return true
	else:
	    return false

    def cache_get_precision(self):
	return self.precision_relative()

    def copy(self):
	return pAdicRingCappedRelativeElement(self.parent(), (self.valuation(), self._unit, self._relprec), construct = True)

    def exp(self):
        raise NotImplementedError

    def exp_artin_hasse(self):
        raise NotImplementedError
	#E_p(x) = exp(x + x^p/p + x^(p^2)/p^2 + ...)

    def gamma(self):
        raise NotImplementedError

    def is_zero(self, prec):
	return (self._relprec <= 0) or (self.valuation() >= prec)

    def is_equal_to(self, right, prec):
	return (self - right).is_zero(prec)

    def lift(self):
	if self.valuation() == infinity:
	    return 0
	else:
	    return self.prime_pow(self.valuation()) * self.unit_part()

    def list(self):
	if (self.valuation() == infinity) or (self.precision_relative() <= 0):
	    return []
	else:
	    def plist(n, p, prec):
		if prec == 0:
		    return []
		else:
		    return [n % p] + plist(n // p, p, prec - 1)
	    return [0 for w in range(self.valuation())] + plist(self._unit.lift(), self.prime(), self.precision_relative())

    def log(self):
        raise NotImplementedError

    def log_artin_hasse(self):
        raise NotImplementedError

    def padded_list(self, n):
	return self.list()[:n] + [0 for w in range(self.precision_absolute(), n)]

    def precision_absolute(self):
	return self._ordp + self._relprec

    def precision_relative(self):
	return self._relprec

    def residue(self, prec):
	return Mod(self.prime_pow(self.valuation()) * self._unit.lift(), self.prime_pow(self.precision_absolute()))

    def square_root(self): ############# Fix bug in integer_mod square root, make sure that precision is correct
	if self.is_square():
	    return pAdicRingCappedRelativeElement(self.parent(), (self.valuation() / 2, self._unit().sqrt(), self._relprec), construct = True)
	else:
	    raise ValueError, "element is not a square" # should eventually be changed to produce an element of a field extension

    def unit_part(self):
	return pAdicRingCappedRelativeElement(self.parent(), (0, self._unit, self._relprec), construct = True)

    def _unit(self):
	return self._unit

    def valuation(self):
	return self._ordp

class pAdicRingLazyElement(pAdicRingGenericElement):
    def __init__(self, parent, x, prec = None, construct = False):
	pAdicRingGenericElement.__init__(self, parent)
	if construct:
	    self._parent = parent
	    (self._p, self._shift, self._known, self._function, self._cache_size) = x
	    return
	if prec == None:
	    prec = parent.precision_cap()
	if isinstance(x, pAdicRingLazyElement):
	    self._parent = parent
	    if parent.prime() != x._p:
		raise ValueError, "Cannot coerce between p-adic rings with different primes."
	    self._p = x._p
	    self._shift = x._shift
	    self._known = x._known
	    self._function = x._function
	    self._cache_size = x._cache_size
	    self.cache_set_precision(prec)
	    return
	if isinstance(x, pAdicRingCappedRelativeElement):
	    self._parent = parent
	    if parent.prime() != x._p:
		raise ValueError, "Cannot coerce between p-adic rings with different primes."
	    self._p = x._p
	    if x.valuation() == infinity:
		self._shift = infinity
		self._known = []
		self._function = lambda x,y: None
		self._cache_size = infinity
	    else:
		self._shift = x.valuation()
		self._known = x.list()
		self._function = lambda x,y: None
		self._cache_size = x.precision_relative()
	    return
	if isinstance(x, pAdicRingFixedModElement):
	    raise NotImplementedError
	if isinstance(x, pAdicRingCappedAbsoluteElement):
	    raise NotImplementedError
	if isinstance(x, pAdicFieldCappedRelativeElement):
	    raise NotImplementedError
	if isinstance(x, pAdicFieldLazyElement):
	    raise NotImplementedError

	self._p = int(parent.prime())
	self._parent = parent

	if isinstance(x, pari_gen) and x.type() == "t_PADIC":
	    t = x.lift()
	    big_oh = x.padicprec(self.prime())
	    if t.type() == 't_INT':
		x = int(t)
	    else:
		x = QQ(t)

	    # We now use the code below, so don't make the next line elif
	if isinstance(x, (int, long)):
	    self._shift = sage.rings.arith.valuation(x, self.prime())
	elif isinstance(x, (integer.Integer, rational.Rational)):
	    self._shift = x.valuation(self.prime())
	else:
	    raise TypeError, "unable to compute valuation"
	if self._shift < 0:
	    raise ValueError, "element not a p-adic integer."
	elif self._shift == infinity:
	    self._known = []
	    self._function = lambda x,y: None
	    self._cache_size = infinity
	    return
	x = x // self.prime_pow(self._shift)


class pAdicRingFixedModElement(pAdicRingGenericElement):

    def __init__(self, parent, x, construct=False):
        r"""
	INPUT:
	    parent -- a pAdicRingFixedMod object.
        Types currently supported:
	    Integers
	    Rationals -- denominator must be relatively prime to p
	    FixedMod p-adics
	Types that should be supported:
	    Finite precision p-adics
	    Lazy p-adics
	    Elements of local extensions of THIS p-adic ring that actually lie in Zp
	    Elements of IntegerModRing(p^k) for k less than or equal to the modulus
    	"""
	pAdicRingGenericElement.__init__(self, parent)
	if construct:
	    (self._modulus, self._value) = x
	    return

	if isinstance(x, pAdicRingGenericElement): #this construction doesn't quite work with lazy p-adics
	    self._modulus = parent.precision_cap()
	    self._value = Mod(x.lift(), self.prime_pow(self._modulus))
	    return

	###########

	self._modulus = parent.precision_cap()

	if isinstance(x, pari_gen) and x.type() == "t_PADIC":
	    t = x.lift()
	    if t.type() == 't_INT':
		x = int(t)
	    else:
		x = QQ(t)

	#Now use the code below to convert from integer or rational, so don't make the next line elif

	if isinstance(x, integer.Integer):
	    self._value = integer_mod.Mod(x, self.prime_pow(self._modulus))
	elif isinstance(x, rational.Rational):
	    val = x.valuation(self.prime())
	    if val < 0:
		raise ValueError, "p divides the denominator"
	    else:
		self._value = integer_mod.Mod(x, self.prime_pow(self._modulus))
	elif isinstance(x, (int, long)):
	    self._value = integer_mod.Mod(x, self.prime_pow(self._modulus))
	else:
	    raise TypeError, "unable to create p-adic element"

    def __cmp__(self, other):
        raise NotImplementedError

    def __invert__(self):
	if self.valuation() > 0:
	    raise ValueError, "cannot invert non-unit"
	else:
	    inverse = 1 / self._value
	    return pAdicRingFixedModElement(self.parent(), (self._modulus, inverse), construct = True)

    def __mod__(self, right):
	val = self.valuation()
	rval = right.valuation()
	if rval > val:
	    raise ValueError, "not enough precision to reduce"
	if right == self.prime_pow(rval):
	    return integer_mod.Mod(self._value, right)
	else:
	    raise ValueError, "modulus must be a power of p" #or should we ignore units

    def __neg__(self):
	return pAdicRingFixedModElement(self.parent(), (self._modulus, -self._value), construct = True)

    def __pow__(self, right):
	right = integer.Integer(right) #Need to make sure that this works for p-adic exponents
	return pAdicRingFixedModElement(self.parent(), (self._modulus, self._value**right), construct = True)

    def _add_(self, right):
	return pAdicRingFixedModElement(self.parent(), (self._modulus, self._value + right._value), construct = True)

    def _div_(self, right):
	if right.valuation() > 0:
	    raise ValueError, "cannot invert non-unit"
	else:
	    return pAdicRingFixedModElement(self.parent(), (self._modulus, self._value / right.value), construct = True)

    def _floordiv_(self, right):
	ppow = self.prime_pow(right.valuation())
	quotient = Mod(self._value.lift() // ppow, self._modulus) / Mod(right._value.lift() / ppow, self._modulus)
	return pAdicRingFixedModElement(self.parent(), (self._modulus, quotient), construct = True)

    def _integer_(self):
	return self._value.lift()

    def _mul_(self, right):
	return pAdicRingFixedModElement(self.parent(), (self._modulus, self._value * right._value), construct = True)

    def _sub_(self, right):
	return pAdicRingFixedModElement(self.parent(), (self._modulus, self._value - right._value), construct = True)

    def add_bigoh(self, prec):
	return pAdicRingFixedModElement(self.parent(), (self._modulus, Mod(Mod(self._value, self.prime_pow(prec)), self.prime_pow(self._modulus))), construct = True)

    def cache_set_precision(self, prec):
	return false

    def cache_get_precision(self):
	return 0

    def copy(self):
	return pAdicRingFixedModElement(self.parent(), (self._modulus, self._value), construct = True)

    def exp(self):
        raise NotImplementedError

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def is_zero(self, prec):
	return Mod(self._value, self.prime_pow(prec)) == 0

    def is_equal_to(self, right, prec): #assumes they have the same parent
	return Mod(self._value, self.prime_pow(prec)) == Mod(right._value, self.prime_pow(prec))

    def lift(self):
	return self._value.lift()

    def list(self):
	def plist(n, p):
	    if n == 0:
		return []
	    else:
		return [n % p] + plist(n // p, p)
	rep = plist(self._value.lift(), self.prime())
	return rep + [0 for w in range(len(rep), self._modulus)]

    def log(self):
        raise NotImplementedError

    def log_artin_hasse(self):
        raise NotImplementedError

    def multiplicative_order(self):
	if self._value % self.prime() == 0:
	    return infinity
	if self._value == 1:
	    return integer.Integer(1)
	if self._value == -1:
	    return integer.Integer(2)
	else:
	    return self._value.multiplicative_order()

    def padded_list(self):
	return self.list()[:n] + [0 for w in range(self.precision_absolute(), n)]

    def precision_absolute(self):
	return self._modulus

    def precision_relative(self):
	return self._modulus - self.valuation()

    def residue(self, prec):
	return Mod(self._value, self.prime_pow(prec))

    def square_root(self):
	if self.is_square():
	    return pAdicRingFixedModElement(self.parent(), (self._modulus, self._value.square_root()), construct = True)
	else:
	    raise ValueError, "element is not a square" # should eventually change to return an element of an extension field

    def _unit(self):
	return Mod(self._value.lift() // self.prime_pow(self.valuation()), self.prime_pow(self._modulus))

    def unit_part(self):
	r"""
	Returns the unit part of self.

	EXAMPLE:
	    sage: R = pAdicRingFixedMod(17,4)
	    sage: a = R(18*17)
	    sage: a.unit_part()
		18
	    sage: type(a)
		<type 'sage.rings.padic_ring_element.pAdicRingFixedModElement'>
	AUTHOR:
	    --David Roe (2006-10)
	"""
	return pAdicRingFixedModElement(self.parent(), (self._modulus, self._unit()), construct = True)

    def valuation(self):
	"""
	Returns the valuation of self.

	EXAMPLE:
	    sage: R = pAdicRingFixedMod(17, 4)
	    sage: a = R(2*17^2)
	    sage: a.valuation()
		2
	AUTHOR:
	 -- David Roe (2006-10)
	"""
	p = self.prime()
	value = self._value
	if value == 0:
	    return self._modulus
	if value % p != 0:
	    return 0
	ppowhigh = p * p
	if value % ppowhigh != 0:
	    return 1
	low = 2
	high = 4
	ppowlow = ppowhigh
	while True:
	    if high >= self._modulus:
		high = self._modulus
		break
	    ppowhigh = ppowlow * ppowlow
	    if value % ppowhigh != 0:
		break
	    low = high
	    high = high + high
	    ppowlow = ppowhigh
	while high - low > 1:
	    mid = floor((high+low)/2)
	    ppowhigh = ppowlow * p**(mid - low)
	    if value % ppowhigh != 0:
		high = mid
	    else:
		low = mid
		ppowlow = ppowhigh
	return low

class pAdicRingCappedAbsoluteElement(pAdicRingFixedModElement):

    def __init__(self, parent, x, prec=None, construct=False):
	pAdicRingGenericElement.__init__(self, parent)
	if construct:
	    (self._modulus, self._value, self._absprec) = x
	    return

	if isinstance(x, pAdicRingGenericElement): #this construction doesn't quite work with lazy p-adics
	    self._modulus = parent.precision_cap()
	    if prec == None:
		self._absprec = x.precision_absolute()
	    else:
		self._absprec = min(x.precision_absolute(), prec)
	    self._value = Mod(Mod(x.lift(), self.prime_pow(self._absprec)), self.prime_pow(self._modulus))
	    return

	###########

	self._modulus = parent.precision_cap()
	self._absprec = parent.precision_cap()

	if isinstance(x, pari_gen) and x.type() == "t_PADIC":
	    t = x.lift()
	    if t.type() == 't_INT':
		x = int(t)
	    else:
		x = QQ(t)

	#Now use the code below to convert from integer or rational, so don't make the next line elif

	if isinstance(x, integer.Integer):
	    self._value = integer_mod.Mod(x, self.prime_pow(self._modulus))
	elif isinstance(x, rational.Rational):
	    val = x.valuation(self.prime())
	    if val < 0:
		raise ValueError, "p divides the denominator"
	    else:
		self._value = integer_mod.Mod(x, self.prime_pow(self._modulus))
	elif isinstance(x, (int, long)):
	    self._value = integer_mod.Mod(x, self.prime_pow(self._modulus))
	else:
	    raise TypeError, "unable to create p-adic element"

    def __invert__(self):
	if self.valuation() > 0:
	    raise ValueError, "cannot invert non-unit"
	else:
	    inverse = 1 / self._value
	    return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, inverse, self._absprec), construct = True)

    def __neg__(self):
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, -self._value, self._absprec), construct = True)

    def __pow__(self, right):
	new = integer.Integer(right) #Need to make sure that this works for p-adic exponents
	val = self.valuation()
	if (val > 0) and (isinstance(right, pAdicRingGenericElement) or isinstance(right, pAdicFieldGenericElement)):
	    raise ValueError, "Can only have p-adic exponent if base is a unit"
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, self._value**right, self._absprec - val + right * val), construct = True)

    def _add_(self, right):
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, self._value + right._value,
                                  min(self.precision_absolute(), right.precision_absolute())), construct = True)

    def _div_(self, right):
	if right.valuation() > 0:
	    raise ValueError, "cannot invert non-unit"
	else:
	    return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus,
                              self._value / right.value, self.valuation() +
                               min(self.precision_relative(), right.precision_relative())), construct = True)

    def _floordiv_(self, right):
	selfval = self.valuation()
	rightval = right.valuation()
	ppow = self.prime_pow(rightval)
	toppart = self._value.lift() // ppow
	bottompart = right._value.lift() // ppow
	if rightval > selfval:
	    val = arith.valuation(toppart, self.prime())
	    toprel = self.precision_relative() - (rightval - selfval) - val
	    bottomrel = right.precision_relative()
	    return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus,
                               Mod(toppart, self._modulus) / Mod(bottompart, self._modulus),
                                                  val + min(toprel, bottomrel)), construct = True)
	else:
	    return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus,
                               Mod(toppart, self._modulus) / Mod(bottompart, self._modulus),
                               selfval - rightval + min(self.precision_relative()), right.precision_relative()))

    def _mul_(self, right):
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus,
                            self._value * right._value, self.valuation() + right.valuation() + min(self.precision_relative(),
                                                  right.precision_relative())), construct = True)

    def _sub_(self, right):
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, self._value - right._value,
                                              min(self.precision_absolute(), right.precision_absolute())), construct = True)

    def add_bigoh(self, prec):
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus,
                   Mod(Mod(self._value, self.prime_pow(prec)), self.prime_pow(self._modulus)), min(prec, self._modulus)), construct = True)

    def cache_set_precision(self, prec):
	return prec >= self.precision_relative()

    def cache_get_precision(self):
	return self.precision_relative()

    def copy(self):
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, self._value, self._absprec), construct = True)

    def exp(self):
        raise NotImplementedError

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def list(self):
	def plist(n, p):
	    if n == 0:
		return []
	    else:
		return [n % p] + plist(n // p, p)
	rep = plist(self._value.lift(), self.prime())
	return rep + [0 for w in range(len(rep), self.precision_absolute())]

    def log(self):
        raise NotImplementedError

    def log_artin_hasse(self):
        raise NotImplementedError

    def multiplicative_order(self):
	if self._value % self.prime() == 0:
	    return infinity
	if self._value == 1:
	    return integer.Integer(1)
	if self._value == -1:
	    return integer.Integer(2)
	else:
	    return Mod(self._value, self.prime_pow(self.precision_absolute())).multiplicative_order()

    def padded_list(self):
	return self.list()[:n] + [0 for w in range(self.precision_absolute(), n)]

    def precision_absolute(self):
	return self._absprec

    def precision_relative(self):
	return self._absprec - self.valuation()

    def square_root(self):
	if self.is_square():
	    return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, self._value.square_root(), self.precision_absolute() - self.valuation() // 2), construct = True)
	else:
	    raise ValueError, "element is not a square" # should eventually change to return an element of an extension field

    def unit_part(self):
	r"""
	Returns the unit part of self.

	EXAMPLE:
	    sage: R = pAdicRingFixedMod(17,4)
	    sage: a = R(18*17)
	    sage: a.unit_part()
		18
	    sage: type(a)
		<type 'sage.rings.padic_ring_element.pAdicRingFixedModElement'>
	AUTHOR:
	    --David Roe (2006-10)
	"""
	return pAdicRingCappedAbsoluteElement(self.parent(), (self._modulus, self._unit(), self.precision_relative()), construct = True)





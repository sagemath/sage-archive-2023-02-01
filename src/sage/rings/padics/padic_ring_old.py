import weakref

import sage.rings.ring
from sage.rings.infinity import infinity
from sage.rings.integer_mod_ring import IntegerModRing
import sage.rings.padics06.padic_ring_element
import sage.rings.integer
import sage.rings.rational

class pAdicRingGeneric(ring.CommutativeRing):
    def characteristic(self):
	pass

class pAdicRing(pAdicRingGeneric): #, ring.DiscreteValuationRing):
    def __call__(self, x, prec=infinity):
	pass
    def __cmp__(self, other):
	pass
    def __init__(self, p, series_print = True, print_prec = infinity):
	pass
    def __call__(self, x, prec=infinity):
	pass
    def __cmp__(self, other):
	pass
    def _coerce_(self, x):
	pass
    def _repr_(self, do_latex=False):
	pass
    def get_default_prec(self):
	pass
    def get_print_mode(self):
	pass
    def set_print_mode(self, str):
	pass
    def prime(self):
	pass
    def residue_characteristic(self):
	pass
    def residue_class_field(self):
	pass
    def residue_system(self):
	pass
    def teichmuller(self, x):
	pass
    def teichmuller_system(self):
	pass
    def absolute_discriminant(self):
	pass
    def discriminant(self, K=None):
	pass
    def different(self):
	pass
    def defining_polynomial(self):
	pass
    def base_ring(self):
	pass
    def ground_ring(self):
	pass
    def degree(self):
	pass
    def ramification_index(self):
	pass
    def e(self):
	pass
    def inertia_degree(self):
	pass
    def f(self):
	pass
    def inertia_subring(self):
	pass
    def get_extension(self):
	pass
    def automorphisms(self):
	pass
    def galois_group(self):
	pass
    def is_abelian(self):
	pass
    def is_normal(self):
	pass
    def uniformizer(self):
	pass
    def fraction_field(self):
	pass
    def has_pth_root(self):
	pass
    def has_root_of_unity(self, n):
	pass
    def hasGNB(self):
	pass
    def hom(self, ring):
	pass
    def is_isomorphic(self, ring):
	pass
    def random_element(self):
	pass
    def unit_group(self):
	pass
    def unit_group_gens(self):
	pass
    def principal_unit_group(self):
	pass

class pAdicRingFixedMod(pAdicRingGeneric, ring.NoetherianRing):
    r"""
    An implementation of the p-adic integers using fixed modulus.
    """

    def __init__(self, p, mod = 20, print_mode = 'val-unit'):
        r"""
        INPUT:
            p -- a prime integer
            mod -- the fixed modulus of this ring, ie the number such that this is isomorphic to $\Z / p^mod \Z$
            print_mode -- the printing mode for this ring.  'val-unit' will print as p^k*u, 'integer' as just an integer and 'series' as a series in $p$.
        """
        self.__p = p
        self.__mod = mod
	self.set_print_mode(print_mode)

    def __call__(self, x):
        r"""
            Supports casting into this ring from other rings.  Note that any precision information from other types of p-adics will be lost.
            Types currently supported:
		Integers
		Rationals -- denominator must be relatively prime to p
            Types that will be supported:
		Finite precision p-adics
		Lazy p-adics
		Elements of local extensions of THIS p-adic ring that actually lie in here.
	"""
	return pAdicRingFixedModElement(self, x)

    def __cmp__(self, other):
        if not isinstance(other, pAdicRingGeneric):
	    return -1
	if self.__p < other.__p:
	    return -1
	elif self.__p > other.__p:
	    return 1
	return 0

    def _coerce_(self, x):
	if x.parent() is self:
	    return x
	if isinstance(x, (integer.Integer, rational.Rational)):
	    return self(x)
#need to add other canonical coercions
	raise TypeError

    def _repr_(self, do_latex=False):
	return "%s-adic Ring of fixed modulus %s"%(self.__p, self.__mod)

    def get_default_prec(self):
	"""
	Returns the modulus of this fixed modulus ring.

	Note that this is different from the precisions associated with other types of p-adic ring.
	EXAMPLES:
	    sage: R = pAdicRingFixedMod(7, 5)
	    sage: R.get_default_prec(5)
		5
     	"""
	return self.__mod

    def get_print_mode(self):
	"""
	Returns the current print mode as a string.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(7,5)
	    sage: R.get_print_mode()
		'val-unit'
	"""
	return self.__print_mode

    def set_print_mode(self, print_mode):
        """
	Sets the print mode.

	The options are:
	    'val-unit' -- elements are displayed as p^k*u
	    'integer' -- elements are displayed as an integer
	    'series' -- elements are displayed as series in p
	    'val-unit-p' -- same as val-unit, except that p is written as "p"
	    'integer-p' -- same as integer, except that p is written as "p"
	    'series-p' -- same as series, except that p is written as "p"
	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3,5)
	    sage: a = R(117); a
		3^2*13 + O(3^5)
	    sage: R.set_print_mode('integer'); a
		117 + O(3^5)
	    sage: R.set_print_mode('series'); a
		1*3^2 + 1*3^3 + 1*3^4 + O(3^5)
	    sage: R.set_print_mode('val-unit-p'); a
		p^2*13 + O(p^5)
	    sage: R.set_print_mode('integer-p'); a
		117 + O(p^5)
	    sage: R.set_print_mode('series-p'); a
		1*p^2 + 1*p^3 + 1*p^4 + O(p^5)
	"""
	if (print_mode in ['val-unit', 'integer', 'series', 'val-unit-p', 'integer-p', 'series-p']):
            self.__print_mode = print_mode
        else:
            raise ValueError, "print_mode must be either val-unit, integer, series, val-unit-p, integer-p, or series-p"

    def is_atomic_repr(self):
	"""
	Return False, since we want p-adics to be printed with parentheses around them
	when they are coefficients, e.g., in a polynomial.
	"""
	return False

    def is_field(self):
	"""
	Return False.
	"""
	return False

    def prime(self):
	"""
	Returns the prime, ie the characteristic of the residue field.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3,5)
	    sage: R.prime()
		3
	"""
	return self.__p

    def modulus(self):
	"""
	Returns the modulus of this fixed modulus p-adic ring.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3,5)
	    sage: R.modulus()
		5
	"""
	return self.__mod

    def residue_characteristic(self):
	"""
	Returns the prime, ie the characteristic of the residue field.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3,5)
	    sage: R.residue_characteristic()
		3
	"""
	return self.__p

    def residue_class_field(self):
	"""
	Returns the residue class field.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3,5)
	    sage: k = R.residue_class_field()
	    sage: k
		Finite Field with 3 elements
	"""
	return IntegerModRing(self.__p)

    def residue_system(self):
	"""
	Returns a list of elements representing all the residue classes.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3, 5)
	    sage: R.residue_system(self)
		[0 + O(3^5), 1 + O(3^5), 2 + O(3^5)]
	"""
	return [self(i) for i in range(self.__p)]

    def teichmuller(self, x):
	r"""
	Returns the teichmuller representative of x.

	INPUT:
	    x -- an integer or element of $\Z / p\Z$ that is not divisible by $p$
	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3, 5)
	    sage: R.teichmuller(2)
		????????????
	"""
	raise NotImplementedError

    def teichmuller_system(self):
	r"""
	Returns a set of teichmuller representatives for the invertible elements of $\Z / p\Z$.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3, 5)
	    sage: R.teichmuller_system()
		??????????
	"""
	return [self.teichmuller(i) for i in range(1, self.__p)]

    def absolute_discriminant(self):
	pass
    def discriminant(self, K=None):
	pass
    def different(self):
	pass
    def defining_polynomial(self):
	from polynomial_ring import PolynomialRing
	x = PolynomialRing(self).gen()
	self.__defining_polynomial = x - 1
	return self.__defining_polynomial

    def base_ring(self):
	r"""
	Returns self.

	For consistency with extensions.
	"""
	return self

    def ground_ring(self):
	r"""
	Returns self.

	For consistency with extensions.
	"""
	return self

    def degree(self):
	r"""
	Returns 1.
	"""

	return 1

    def ramification_index(self):
	"""
	Returns the ramification index, i.e. 1.
	"""
	return 1

    def e(self):
	"""
	Returns the ramification index, i.e. 1.
	"""
	return 1

    def inertia_degree(self):
	"""
	Returns the inertia degree, i.e. 1.
	"""
	return 1

    def f(self):
	"""
	Returns the inertia degree, i.e. 1.
	"""
	return 1

    def inertia_subring(self):
	"""
	Returns the inertia subring, i.e. self.
	"""
	return self

    def get_extension(self):
	"""
	Returns the trivial extension of self.
	"""
	raise NotImplementedError

    def automorphisms(self):
	r"""
	Returns the group of automorphisms of $\Z_p$, i.e. the trivial group.
	"""
	raise NotImplementedError

    def galois_group(self):
	r"""
	Returns the Galois group of $\Z_p$, i.e. the trivial group.
	"""
	raise NotImplementedError

    def is_abelian(self):
	"""
	Returns whether the Galois group is abelian, i.e. True.  #should this be automorphism group?
	"""
	return True

    def is_normal(self):
	"""
	Returns whether or not this is a normal extension, i.e. True.
	"""
	return True

    def uniformizer(self):
	"""
	Returns a uniformizer for this ring.

	EXAMPLES:
	    sage: R = pAdicRingFixedMod(3,5)
	    sage: R.uniformizer()
		3
	"""
	return self.__p

    def fraction_field(self):
	r"""
	Would normally return $\Q_p$, but there is no implementation of $Q_p$ matching this ring so this raises an error
	"""
	raise TypeError, "This implementation of the p-adic ring does not support fields of fractions."

    def has_pth_root(self):
	r"""
	Returns whether or not $\Z_p$ has a $p^{\mbox{th}}$ root of unity, ie false.
	"""
	return False

    def has_root_of_unity(self, n):
	r"""
	Returns whether or not $\Z_p$ has an $n^{\mbox{th}}$ root of unity.
	"""
	raise NotImplementedError

    def hasGNB(self):
	r"""
	Returns whether or not $\Z_p$ has a Gauss Normal Basis.
	"""
	raise NotImplementedError

    def hom(self, ring):
	r"""
	Returns the set of homomorphisms from $Z_p$ into another ring.
	"""
	raise NotImplementedError

    def is_isomorphic(self, ring):
	pass
    def random_element(self):
	pass
    def unit_group(self):
	pass
    def unit_group_gens(self):
	pass
    def principal_unit_group(self):
	pass

class pAdicRingLazy(pAdicRingGeneric): #, ring.DiscreteValuationRing):
    def __init__(self, p, series_print = True, print_prec = infinity):
	pass
    def __call__(self, x, prec=infinity):
	pass
    def __cmp__(self, other):
	pass
    def _coerce_(self, x):
	pass
    def _repr_(self, do_latex=False):
	pass
    def get_default_prec(self):
	pass
    def get_print_mode(self):
	pass
    def set_print_mode(self, str):
	pass
    def prime(self):
	pass
    def residue_characteristic(self):
	pass
    def residue_class_field(self):
	pass
    def residue_system(self):
	pass
    def teichmuller(self, x):
	pass
    def teichmuller_system(self):
	pass
    def absolute_discriminant(self):
	pass
    def discriminant(self, K=None):
	pass
    def different(self):
	pass
    def defining_polynomial(self):
	pass
    def base_ring(self):
	pass
    def ground_ring(self):
	pass
    def degree(self):
	pass
    def ramification_index(self):
	pass
    def e(self):
	pass
    def inertia_degree(self):
	pass
    def f(self):
	pass
    def inertia_subring(self):
	pass
    def get_extension(self):
	pass
    def automorphisms(self):
	pass
    def galois_group(self):
	pass
    def is_abelian(self):
	pass
    def is_normal(self):
	pass
    def uniformizer(self):
	pass
    def fraction_field(self):
	pass
    def has_pth_root(self):
	pass
    def has_root_of_unity(self, n):
	pass
    def hasGNB(self):
	pass
    def hom(self, ring):
	pass
    def is_isomorphic(self, ring):
	pass
    def random_element(self):
	pass
    def unit_group(self):
	pass
    def unit_group_gens(self):
	pass
    def principal_unit_group(self):
	pass

def Zp(p, default_prec = 20, cap_prec = infinity, mod = infinity, lazy = False):
    if (lazy == False):
        if (mod == infinity):
            if (cap_prec == infinity):
                return pAdicRing(p, default_prec)
            else:
                return pAdicRingCappedPrec(p, cap_prec)
        else:
	    return pAdicRingFixedMod(p, mod)
    else:
        return pAdicRingLazy(p)

#def is_pAdicRing(x):



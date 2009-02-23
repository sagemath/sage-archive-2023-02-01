

#############################################################
##                                                         ##
##  Wrappers for the Genus/Genus Symbol Code in ../genera/ ##
##                                                         ##
#############################################################

from sage.quadratic_forms.genera.genus import Genus, LocalGenusSymbol, \
        is_GlobalGenus, is_2_adic_genus, canonical_2_adic_compartments, \
        canonical_2_adic_trains, canonical_2_adic_reduction, \
        basis_complement, signature_of_matrix, p_adic_symbol, is_even, \
        split_odd,  trace_diag, two_adic_symbol, is_trivial_symbol
        #GenusSymbol_p_adic_ring, GenusSymbol_global_ring
## NOTE: Removed the signature routine here... and rewrote it for now.


from sage.rings.integer_ring import IntegerRing
from sage.rings.arith import is_prime, prime_divisors



def global_genus_symbol(self):
    """
    Returns the genus of a two times a quadratic form over ZZ.  These
    are defined by a collection of local genus symbols (a la Chapter
    15 of Conway-Sloane), and a signature.

    EXAMPLES:
    """
    ## Check that the form is defined over ZZ
    if not self.base_ring() == IntegerRing():
        raise TypeError, "Oops!  The quadratic form is not defined over the integers."

    ## Return the result
    try:
        return Genus(self.Hessian_matrix())
    except:
        raise TypeError, "Oops!  There is a problem computing the genus symbols for this form."



def local_genus_symbol(self, p):
    """
    Returns the Conway-Sloane genus symbol of 2 times a quadratic form defined
    over ZZ at a prime number p.

    INPUT:
        p -- a prime number > 0

    OUTPUT:
        Returns a Conway-Sloane genus symbol.

    EXAMPLES:
    """
    ## Check that p is prime and that the form is defined over ZZ.
    if not is_prime(p):
        raise TypeError, "Oops!  The number " + str(p) + " isn't prime."
    if not self.base_ring() == IntegerRing():
        raise TypeError, "Oops!  The quadratic form is not defined over the integers."

    ## Return the result
    try:
        M = self.Hessian_matrix()
        return LocalGenusSymbol(M, p)
    except:
        raise TypeError, "Oops!  There is a problem computing the local genus symbol at the prime " + str(p) + " for this form."






def CS_genus_symbol_list(self, force_recomputation=False):
    """
    Returns the list of Conway-Sloane genus symbols in increasing order of primes dividing 2*det.
    """
    ## Try to use the cached list
    if force_recomputation == False:
        try:
            return self.__CS_genus_symbol_list
        except:
            pass

    ## Otherwise recompute and cache the list
    list_of_CS_genus_symbols = [ ]

    for p in prime_divisors(2 * self.det()):
        list_of_CS_genus_symbols.append(self.local_genus_symbol(p))

    self.__CS_genus_symbol_list = list_of_CS_genus_symbols
    return list_of_CS_genus_symbols

include "sage/libs/linkages/padics/Polynomial_ram.pxi"
include "FM_template.pxi"

cdef class RelativeRamifiedFixedModElement(FMElement):
    def __cinit__(self, parent=None, x=None, absprec=infinity, relprec=infinity):
        # It's not possible to set self.value in cconstruct (because of the calling syntax)
        # so we do it here.
        cdef type t
        if parent is not None:
            t = type((<PowComputer_?>parent.prime_pow).modulus)
            self.value = t.__new__(t)

    cdef FMElement _new_c(self):
        """
        Creates a new element in this ring.

        This is meant to be the fastest way to create such an element; much
        faster than going through the usual mechanisms which involve
        ``__init__``.

        TESTS::

            sage: R = ZpFM(5); R(6) * R(7) #indirect doctest
            2 + 3*5 + 5^2 + O(5^20)

        """
        cdef type t = type(self)
        cdef type polyt = type(self.prime_pow.modulus)
        cdef FMElement ans = t.__new__(t)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        ans.value = polyt.__new__(polyt)
        cconstruct(ans.value, ans.prime_pow)
        return ans

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        EXAMPLES::

            sage: R.<a> = ZqFM(9)
            sage: S.<x> = ZZ[]
            sage: W.<w> = R.extension(x^2 - 3)
            sage: loads(dumps(w)) == w
            True
        """
        return unpickle_fme_rel_v2, (self.__class__, self.parent(), cpickle(self.value, self.prime_pow))

    def _poly_rep(self):
        return self.value

def unpickle_fme_rel_v2(cls, parent, value):
    """
    Unpickles a fixed mod element.

    EXAMPLES::

        sage: from sage.rings.padics.relative_ramified_FM import RelativeRamifiedFixedModElement, unpickle_fme_rel_v2
        sage: R.<a> = ZqFM(9)
        sage: S.<x> = ZZ[]
        sage: W.<w> = R.extension(x^2 - 3)
        sage: u = unpickle_fme_rel_v2(RelativeRamifiedFixedModElement, W, [a, 1]); u
        sage: u.parent() is W
        True
    """
    cdef RelativeRamifiedFixedModElement ans = cls.__new__(cls)
    ans._parent = parent
    ans.prime_pow = <PowComputer_?>parent.prime_pow
    cdef type polyt = type(ans.prime_pow.modulus)
    ans.value = polyt.__new__(polyt)
    cconstruct(ans.value, ans.prime_pow)
    cunpickle(ans.value, value, ans.prime_pow)
    return ans

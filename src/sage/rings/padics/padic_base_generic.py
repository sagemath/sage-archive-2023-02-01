"""
p-Adic Base Generic.

A superclass for implementations of Zp and Qp.

AUTHORS::

    - David Roe
"""

from padic_generic import pAdicGeneric
from sage.rings.infinity import infinity
from sage.rings.padics.pow_computer import PowComputer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

class pAdicBaseGeneric(pAdicGeneric):
    def __init__(self, p, prec, print_mode, names, element_class):
        """
        Initialization

        TESTS::

            sage: R = Zp(5) #indirect doctest
        """
        self.prime_pow = PowComputer(p, max(min(prec - 1, 30), 1), prec, self.is_field())
        pAdicGeneric.__init__(self, self, p, prec, print_mode, names, element_class)
        if self.is_field():
            self._populate_coercion_lists_(coerce_list=[QQ], element_constructor=element_class)
        else:
            self._populate_coercion_lists_(coerce_list=[ZZ], element_constructor=element_class)

#     def __reduce__(self):
#         """
#         For pickling.

#         This function is provided because prime_pow needs to be set
#         before _printer, so the standard unpickling fails.

#         TESTS::

#             sage: R = Zp(5, 15, print_mode='bars', print_sep='&'); S = loads(dumps(R))
#             sage: R == S
#             True
#             sage: S(2777)
#         """
#         from sage.rings.padics.factory import Zp, Qp
#         if self.is_field():
#             if self.is_capped_relative():
#                 return Qp, (self.prime(), self.precision_cap(), 'capped-rel', self.print_mode(), 40, self.variable_name())
#             else:
#                 return Qp, (self.prime(), self.precision_cap(), 'lazy', self.print_mode(), self.halting_parameter(), self.variable_name())
#         else:
#             if self.is_capped_relative():
#                 return Zp, (self.prime(), self.precision_cap(), 'capped-rel', self.print_mode(), 40, self.variable_name())
#             elif self.is_capped_absolute():
#                 return Zp, (self.prime(), self.precision_cap(), 'capped-abs', self.print_mode(), 40, self.variable_name())
#             elif self.is_fixed_mod():
#                 return Zp, (self.prime(), self.precision_cap(), 'fixed-mod', self.print_mode(), 40, self.variable_name())
#             else:
#                 return Zp, (self.prime(), self.precision_cap(), 'lazy', self.print_mode(), self.halting_parameter(), self.variable_name())

    def fraction_field(self, print_mode=None):
        r"""
        Returns the fraction field of self.

        INPUT::

            - print_mode - a dictionary containing print options.
              Defaults to the same options as this ring.

        OUTPUT::

            the fraction field of self.

        EXAMPLES::

            sage: R = Zp(5, print_mode='digits')
            sage: K = R.fraction_field(); repr(K(1/3))[3:]
            '31313131313131313132'
            sage: L = R.fraction_field({'max_ram_terms':4}); repr(L(1/3))[3:]
            '3132'
        """
        if self.is_field() and print_mode is None:
            return self
        from sage.rings.padics.factory import Qp
        return Qp(self.prime(), self.precision_cap(), 'capped-rel', print_mode=self._modified_print_mode(print_mode), names=self._uniformizer_print())

    def integer_ring(self, print_mode=None):
        r"""
        Returns the integer ring of self, possibly with print_mode changed.

        INPUT::

            - print_mode - a dictionary containing print options.
              Defaults to the same options as this ring.

        OUTPUT::

            - The ring of integral elements in self.

        EXAMPLES::

            sage: K = Qp(5, print_mode='digits')
            sage: R = K.integer_ring(); repr(R(1/3))[3:]
            '31313131313131313132'
            sage: S = K.integer_ring({'max_ram_terms':4}); repr(S(1/3))[3:]
            '3132'
        """
        if not self.is_field() and print_mode is None:
            return self
        from sage.rings.padics.factory import Zp
        return Zp(self.prime(), self.precision_cap(), self._prec_type(), print_mode=self._modified_print_mode(print_mode), names=self._uniformizer_print())

    def is_isomorphic(self, ring):
        r"""
        Returns whether self and ring are isomorphic, i.e. whether ring is an implementation of $\Z_p$ for the same prime as self.

        INPUT::

            self -- a p-adic ring
            ring -- a ring

        OUTPUT::

            boolean -- whether ring is an implementation of $\Z_p$ for the same prime as self.

        EXAMPLES::

            sage: R = Zp(5, 15, print_mode='digits'); S = Zp(5, 44, print_max_terms=4); R.is_isomorphic(S)
            True
        """
        return isinstance(ring, pAdicBaseGeneric) and self.prime() == ring.prime() and self.is_field() == ring.is_field()

    def gen(self, n=0):
        """
        Returns the nth generator of this extension.  For base
        rings/fields, we consider the generator to be the prime.

        EXAMPLES::

            sage: R = Zp(5); R.gen()
            5 + O(5^21)
        """
        if n != 0:
            raise IndexError, "only one generator"
        return self(self.prime())

    def absolute_discriminant(self):
        """
        Returns the absolute discriminant of this p-adic ring

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            integer -- the absolute discriminant of this ring, i.e., 1

        EXAMPLES::

            sage: Zp(5).absolute_discriminant()
            1
        """
        return 1

    def discriminant(self, K=None):
        """
        Returns the discriminant of this p-adic ring over K

        INPUT::

            self -- a p-adic ring
            K -- a sub-ring of self or None (default: None)

        OUTPUT::

            - integer -- the discriminant of this ring over K (or the
              absolute discriminant if K is None)

        EXAMPLES::

            sage: Zp(5).discriminant()
            1
        """
        if (K is None or K is self):
            return 1
        else:
            raise ValueError, "Ground Ring must be a subring of self"

    def is_abelian(self):
        """
        Returns whether the Galois group is abelian, i.e. True.  #should this be automorphism group?

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            boolean -- whether self is abelian, i.e., True

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.is_abelian()
            True
        """
        return True

    def is_normal(self):
        """
        Returns whether or not this is a normal extension, i.e. True.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            boolean -- whether self is normal, i.e., True

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.is_normal()
            True
        """
        return True

    def uniformizer(self):
        """
        Returns a uniformizer for this ring.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            element -- self's uniformizer

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod', 'series')
            sage: R.uniformizer()
            3 + O(3^5)
        """
        return self(self.prime_pow._prime())

    def uniformizer_pow(self, n):
        """
        Returns the nth power of the uniformizer of self (as an element of self).

        EXAMPLES::

            sage: R = Zp(5)
            sage: R.uniformizer_pow(5)
            5^5 + O(5^25)
            sage: R.uniformizer_pow(infinity)
            0
        """
        return self(self.prime_pow(n))

    def _uniformizer_print(self):
        """
        Returns how the uniformizer is supposed to print.

        EXAMPLES::

            sage: R = Zp(5, names='pi'); R._uniformizer_print()
            'pi'
        """
        return self.variable_name()

    def has_pth_root(self):
        r"""
        Returns whether or not $\Z_p$ has a primitive $p^{\mbox{th}}$ root of unity.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            boolean -- whether self has primitive $p^{\mbox{th}}$ root of unity.

        EXAMPLES::

            sage: Zp(2).has_pth_root()
            True
            sage: Zp(17).has_pth_root()
            False
        """
        return (self.prime() == 2)

    def has_root_of_unity(self, n):
        r"""
        Returns whether or not $\Z_p$ has a primitive $n^{\mbox{th}}$ root of unity.

        INPUT::

            - self -- a p-adic ring
            - n -- an integer

        OUTPUT::

            boolean -- whether self has primitive $n^{\mbox{th}}$ root of unity

        EXAMPLES::

            sage: R=Zp(37)
            sage: R.has_root_of_unity(12)
            True
            sage: R.has_root_of_unity(11)
            False
        """
	if (self.prime() == 2):
	    return n.divides(2)
        else:
	    return n.divides(self.prime() - 1)

    def zeta(self, n=None):
        r"""
        Returns a generator of the group of roots of unity.

        INPUT::

            - self -- a p-adic ring
            - n -- an integer ot None (default: None)

        OUTPUT::

            - element -- a generator of the $n^\mbox{th}$ roots of
              unity, or a generator of the full group of roots of
              unity if n is None

        EXAMPLES::

            sage: R = Zp(37,5)
            sage: R.zeta(12)
            8 + 24*37 + 37^2 + 29*37^3 + 23*37^4 + O(37^5)
        """
        if (self.prime() == 2):
            if (n is None) or (n == 2):
                return self(-1)
            if n == 1:
                return self(1)
            else:
                raise ValueError, "No, %sth root of unity in self"%n
        else:
            from sage.rings.finite_field import GF
            return self.teichmuller(GF(self.prime()).zeta(n).lift())

    def zeta_order(self):
        """
        Returns the order of the group of roots of unity.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            integer -- the order of the group of roots of unity

        EXAMPLES::

            sage: R = Zp(37); R.zeta_order()
            36
            sage: Zp(2).zeta_order()
            2
        """
        if (self.prime() == 2):
            return 2
        else:
            return self.prime() - 1

    def plot(self, max_points=2500, **args):
        r"""
        Creates a visualization of this p-adic ring as a fractal similar as
        a generalization of the the Sierpi\'nski triangle. The resulting image
        attempts to capture the algebraic and topological characteristics of $\Z_p$.

        INPUT::

            - point_count -- the maximum number or points to plot,
              which controls the depth of recursion (default 2500)

            - **args -- color, size, etc. that are passed to the
              underlying point graphics objects

        REFERENCES::

            Cuoco, A. ``Visualizing the $p$-adic Integers'', The American
                Mathematical Monthly, Vol. 98, No. 4 (Apr., 1991), pp. 355-364

        EXAMPLES::

            sage: Zp(3).plot()
            sage: Zp(5).plot(max_points=625)
            sage: Zp(23).plot(rgbcolor=(1,0,0))
        """
        if not args.has_key('pointsize'):
            args['pointsize'] = 1
        from sage.misc.mrange import cartesian_product_iterator
        from sage.rings.real_double import RDF
        from sage.plot.plot import points, circle, Graphics
        p = self.prime()
        phi = 2*RDF.pi()/p
        V = RDF**2
        vs = [V([(phi*t).sin(), (phi*t).cos()]) for t in range(p)]
        all = []
        depth = max(RDF(max_points).log(p).floor(), 1)
        scale = min(RDF(1.5/p), 1/RDF(3))
        pts = [vs]*depth
        if depth == 1 and 23 < p < max_points:
            extras = int(max_points/p)
            if p/extras > 5:
                pts = [vs]*depth + [vs[::extras]]
        for digits in cartesian_product_iterator(pts):
            p = sum([v * scale**n for n, v in enumerate(digits)])
            all.append(tuple(p))
        g = points(all, **args)
        # Set default plotting options
        g._Graphics__show_axes = False
        g._Graphics__aspect_ratio = 1
        return g

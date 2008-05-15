import padic_generic
import sage.rings.infinity
from sage.rings.padics.pow_computer import PowComputer

infinity = sage.rings.infinity.infinity

class pAdicBaseGeneric(padic_generic.pAdicGeneric):
    def __init__(self, p, prec, print_mode, names, element_class):
        self.prime_pow = PowComputer(p, max(min(prec - 1, 30), 1), prec, self.is_field())
        padic_generic.pAdicGeneric.__init__(self, self, p, prec, print_mode, names, element_class)

    def __reduce__(self):
        """
        For pickling.

        This function is provided because prime_pow needs to be set before _printer, so the standard unpickling fails.
        """
        from sage.rings.padics.factory import Zp, Qp
        if self.is_field():
            if self.is_capped_relative():
                return Qp, (self.prime(), self.precision_cap(), 'capped-rel', self.print_mode(), 40, self.variable_name())
            else:
                return Qp, (self.prime(), self.precision_cap(), 'lazy', self.print_mode(), self.halting_parameter(), self.variable_name())
        else:
            if self.is_capped_relative():
                return Zp, (self.prime(), self.precision_cap(), 'capped-rel', self.print_mode(), 40, self.variable_name())
            elif self.is_capped_absolute():
                return Zp, (self.prime(), self.precision_cap(), 'capped-abs', self.print_mode(), 40, self.variable_name())
            elif self.is_fixed_mod():
                return Zp, (self.prime(), self.precision_cap(), 'fixed-mod', self.print_mode(), 40, self.variable_name())
            else:
                return Zp, (self.prime(), self.precision_cap(), 'lazy', self.print_mode(), self.halting_parameter(), self.variable_name())

    def is_isomorphic(self, ring):
        r"""
        Returns whether self and ring are isomorphic, i.e. whether ring is an implementation of $\Z_p$ for the same prime as self.

        INPUT:
            self -- a p-adic ring
            ring -- a ring

        OUTPUT:
            boolean -- whether ring is an implementation of $\Z_p$ for the same prime as self.
        """
        return is_instance(ring, pAdicBaseGeneric) and self.prime() == ring.prime() and self.is_field() == ring.is_field()

    def uniformizer_pow(self, n):
        return self(self.prime_pow(n))

    uniformiser_pow = uniformizer_pow

    def _uniformizer_print(self):
        return self.variable_name()

    def plot(self, max_points=2500, **args):
        """
        Creates a visualization of this p-adic ring as a fractal similar as
        a generalization of the the Sierpi\'nski triangle. The resulting image
        attempts to capture the algebraic and topological characteristics of $\Z_p$.

        INPUT:
            point_count -- the maximum number or points to plot, which
                           controls the depth of recursion (default 2500)
            **args      -- color, size, etc. that are passed to the
                           underlying point graphics objects

        REFERENCES:
            Cuoco, A. ``Visualizing the $p$-adic Integers'', The American
                Mathematical Monthly, Vol. 98, No. 4 (Apr., 1991), pp. 355-364

        EXAMPLES:
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

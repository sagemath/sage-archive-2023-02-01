r"""
Rings
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.rngs import Rngs
from category import HomCategory

class Rings(CategoryWithAxiom):
    """
    The category of rings

    Associative rings with unit, not necessarily commutative

    EXAMPLES::

      sage: Rings()
      Category of rings
      sage: sorted(Rings().super_categories(), key=str)
      [Category of rngs, Category of semirings]
      sage: sorted(Rings().axioms())
      ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital', 'Associative', 'Unital']
      sage: Rings() is Rngs().Unital()
      True
      sage: Rings() is Semirings().AdditiveInverse()
      True

    TESTS::

        sage: TestSuite(Rings()).run()

    TODO (see: http://trac.sagemath.org/sage_trac/wiki/CategoriesRoadMap)

     - Make Rings() into a subcategory or alias of Algebras(ZZ);

     - A parent P in the category ``Rings()`` should automatically be
       in the category ``Algebras(P)``.
    """

    _base_category_class_and_axiom = [Rngs, "Unital"]

    class SubcategoryMethods:

        def NoZeroDivisors(self):
            """
            Returns the full subcategory of the objects of ``self`` having no zero divisors

            EXAMPLES::

                sage: Rings().NoZeroDivisors()
                Category of domains

            .. NOTE:: this could be generalized to DistributiveMagmasAndAdditiveMagmas.AdditiveUnital.Unital

            TESTS::

                sage: TestSuite(Rings().NoZeroDivisors()).run()
                sage: Algebras(QQ).NoZeroDivisors.__module__
                'sage.categories.rings'
            """
            return self._with_axiom('NoZeroDivisors')

        def Division(self):
            """
            Returns the full subcategory of the *division* objects of ``self``

            A ring satisfies the division axiom if all non-zero
            elements have multiplicative inverse.

            .. NOTE:: this could be generalized to DistributiveMagmasAndAdditiveMagmas.AdditiveUnital.Unital

            EXAMPLES::

                sage: Rings().Division()
                Category of division rings
                sage: Rings().Commutative().Division()
                Category of fields

            TESTS::

                sage: TestSuite(Rings().Division()).run()
                sage: Algebras(QQ).Division.__module__
                'sage.categories.rings'
            """
            return self._with_axiom('Division')


    NoZeroDivisors = LazyImport('sage.categories.domains', 'Domains', at_startup=True)
    Division       = LazyImport('sage.categories.division_rings', 'DivisionRings', at_startup=True)
    Commutative    = LazyImport('sage.categories.commutative_rings', 'CommutativeRings', at_startup=True)

    class ParentMethods:
        def is_ring(self):
            """
            Return True, since this in an object of the category of rings.

            EXAMPLES::

                sage: Parent(QQ,category=Rings()).is_ring()
                True

            """
            return True

        def bracket(self, x, y):
            """
            Returns the Lie bracket `[x, y] = x y - y x` of `x` and `y`.

            INPUT:

             - ``x``, ``y`` -- elements of ``self``

            EXAMPLES::

                sage: F = AlgebrasWithBasis(QQ).example()
                sage: F
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: a,b,c = F.algebra_generators()
                sage: F.bracket(a,b)
                B[word: ab] - B[word: ba]

            This measures the default of commutation between `x` and `y`.
            `F` endowed with the bracket operation is a Lie algebra;
            in particular, it satisfies Jacobi's identity::

                sage: F.bracket( F.bracket(a,b), c) + F.bracket(F.bracket(b,c),a) + F.bracket(F.bracket(c,a),b)
                0
            """
            return x*y - y*x

        def _Hom_(self, Y, category):
            r"""
            Returns the homset from ``self`` to ``Y`` in the category ``category``

            INPUT::

            - ``Y`` -- a ring
            - ``category`` -- a subcategory of :class:`Rings`() or None

            The sole purpose of this method is to construct the homset
            as a :class:`~sage.rings.homset.RingHomset`. If
            ``category`` is specified and is not a subcategory of
            :class:`Rings`, a ``TypeError`` is raised instead

            This method is not meant to be called directly. Please use
            :func:`sage.categories.homset.Hom` instead.

            EXAMPLES::

                sage: H = QQ._Hom_(QQ, category = Rings()); H
                Set of Homomorphisms from Rational Field to Rational Field
                sage: H.__class__
                <class 'sage.rings.homset.RingHomset_generic_with_category'>

            TESTS::

                sage: Hom(QQ, QQ, category = Rings()).__class__
                <class 'sage.rings.homset.RingHomset_generic_with_category'>

                sage: Hom(CyclotomicField(3), QQ, category = Rings()).__class__
                <class 'sage.rings.number_field.morphism.CyclotomicFieldHomset_with_category'>

                sage: TestSuite(Hom(QQ, QQ, category = Rings())).run() # indirect doctest

            """
            if category is not None and not category.is_subcategory(Rings()):
                raise TypeError, "%s is not a subcategory of Rings()"%category
            if Y not in Rings():
                raise TypeError, "%s is not a ring"%Y
            from sage.rings.homset import RingHomset
            return RingHomset(self, Y, category = category)

        # this is already in sage.rings.ring.Ring,
        # but not all rings descend from that class,
        # e.g., matrix spaces.
        def _mul_(self, x, switch_sides=False):
            """
            Multiplication of rings with, e.g., lists.

            NOTE:

            This method is used to create ideals. It is
            the same as the multiplication method for
            :class:`~sage.rings.ring.Ring`. However, not
            all parents that belong to the category of
            rings also inherits from the base class of
            rings. Therefore, we implemented a ``__mul__``
            method for parents, that calls a ``_mul_``
            method implemented here. See :trac:`7797`.

            INPUT:

            - `x`, an object to multiply with.
            - `switch_sides` (optional bool): If ``False``,
              the product is ``self*x``; if ``True``, the
              product is ``x*self``.

            EXAMPLE:

            As we mentioned above, this method is called
            when a ring is involved that does not inherit
            from the base class of rings. This is the case,
            e.g., for matrix algebras::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: isinstance(MS,Ring)
                False
                sage: MS in Rings()
                True
                sage: MS*2     # indirect doctest
                Left Ideal
                (
                  [2 0]
                  [0 2]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            In the next example, the ring and the other factor switch sides
            in the product::

                sage: [MS.2]*MS
                Right Ideal
                (
                  [0 0]
                  [1 0]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            AUTHOR:

            - Simon King (2011-03-22)

            """
            try:
                if self.is_commutative():
                    return self.ideal(x)
            except (AttributeError,NotImplementedError):
                pass
            try:
                side = x.side()
            except AttributeError:
                return self.ideal(x, side='right' if switch_sides else 'left')
            # presumably x is an ideal...
            try:
                x = x.gens()
            except (AttributeError, NotImplementedError):
                pass # ... not an ideal
            if switch_sides:
                if side in ['right','twosided']:
                    return self.ideal(x,side=side)
                elif side=='left':
                    return self.ideal(x,side='twosided')
            else:
                if side in ['left','twosided']:
                    return self.ideal(x,side=side)
                elif side=='right':
                    return self.ideal(x,side='twosided')
            # duck typing failed
            raise TypeError, "Don't know how to transform %s into an ideal of %s"%(x,self)

        @cached_method
        def ideal_monoid(self):
            """
            The monoid of the ideals of this ring.

            NOTE:

            The code is copied from the base class of rings.
            This is since there are rings that do not inherit
            from that class, such as matrix algebras.  See
            :trac:`7797`.

            EXAMPLE::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: isinstance(MS,Ring)
                False
                sage: MS in Rings()
                True
                sage: MS.ideal_monoid()
                Monoid of ideals of Full MatrixSpace of 2 by 2 dense matrices
                over Rational Field

            Note that the monoid is cached::

                sage: MS.ideal_monoid() is MS.ideal_monoid()
                True

            """
            try:
                from sage.rings.ideal_monoid import IdealMonoid
                return IdealMonoid(self)
            except TypeError:
                from sage.rings.noncommutative_ideals import IdealMonoid_nc
                return IdealMonoid_nc(self)

        def characteristic(self):
            """
            Return the characteristic of this ring.

            EXAMPLES::

                sage: QQ.characteristic()
                0
                sage: GF(19).characteristic()
                19
                sage: Integers(8).characteristic()
                8
                sage: Zp(5).characteristic()
                0
            """
            from sage.rings.infinity import infinity
            from sage.rings.integer_ring import ZZ
            order_1 = self.one().additive_order()
            return ZZ.zero() if order_1 is infinity else order_1

        def _test_characteristic(self, **options):
            """
            Run generic tests on the method :meth:`characteristic`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: ZZ._test_characteristic()
            """
            tester = self._tester(**options)
            try:
                characteristic = self.characteristic()
            except AttributeError:
                return # raised when self.one() does not have a additive_order()
            except NotImplementedError:
                return

            # test that #12988 is fixed
            from sage.rings.integer import Integer
            tester.assertIsInstance(characteristic, Integer)

        def ideal(self, *args, **kwds):
            """
            Create an ideal of this ring.

            NOTE:

            The code is copied from the base class
            :class:`~sage.rings.ring.Ring`. This is
            because there are rings that do not inherit
            from that class, such as matrix algebras.
            See :trac:`7797`.

            INPUT:

            - An element or a list/tuple/sequence of elements.
            - ``coerce`` (optional bool, default ``True``):
              First coerce the elements into this ring.
            - ``side``, optional string, one of ``"twosided"``
              (default), ``"left"``, ``"right"``: determines
              whether the resulting ideal is twosided, a left
              ideal or a right ideal.

            EXAMPLE::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: isinstance(MS,Ring)
                False
                sage: MS in Rings()
                True
                sage: MS.ideal(2)
                Twosided Ideal
                (
                  [2 0]
                  [0 2]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field
                sage: MS.ideal([MS.0,MS.1],side='right')
                Right Ideal
                (
                  [1 0]
                  [0 0],
                <BLANKLINE>
                  [0 1]
                  [0 0]
                )
                 of Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            """
            if kwds.has_key('coerce'):
                coerce = kwds['coerce']
                del kwds['coerce']
            else:
                coerce = True

            from sage.rings.ideal import Ideal_generic
            from types import GeneratorType
            if len(args) == 0:
                gens = [self(0)]
            else:
                gens = args
                while isinstance(gens, (list, tuple, GeneratorType)) and len(gens) == 1:
                    first = gens[0]
                    if isinstance(first, Ideal_generic):
                        R = first.ring()
                        m = self.convert_map_from(R)
                        if m is not None:
                            gens = [m(g) for g in first.gens()]
                            coerce = False
                        else:
                            m = R.convert_map_from(self)
                            if m is not None:
                                raise NotImplementedError
                            else:
                                raise TypeError
                        break
                    elif isinstance(first, (list, tuple, GeneratorType)):
                        gens = first
                    else:
                        try:
                            if self.has_coerce_map_from(first):
                                gens = first.gens() # we have a ring as argument
                            elif hasattr(first,'parent'):
                                gens = [first]
                            else:
                                raise ArithmeticError, "There is no coercion from %s to %s"%(first,self)
                        except TypeError: # first may be a ring element
                            pass
                        break
            if coerce:
                gens = [self(g) for g in gens]
            from sage.categories.principal_ideal_domains import PrincipalIdealDomains
            if self in PrincipalIdealDomains():
                # Use GCD algorithm to obtain a principal ideal
                g = gens[0]
                if len(gens) == 1:
                    try:
                        g = g.gcd(g) # note: we set g = gcd(g, g) to "canonicalize" the generator: make polynomials monic, etc.
                    except (AttributeError, NotImplementedError):
                        pass
                else:
                    for h in gens[1:]:
                        g = g.gcd(h)
                gens = [g]
            if kwds.has_key('ideal_class'):
                C = kwds['ideal_class']
                del kwds['ideal_class']
            else:
                C = self._ideal_class_(len(gens))
            if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
                gens = gens[0]
            return C(self, gens, **kwds)

        def _ideal_class_(self,n=0):
            """
            Return the class that is used to implement ideals of this ring.

            NOTE:

            We copy the code from :class:`~sage.rings.ring.Ring`. This is
            necessary because not all rings inherit from that class, such
            as matrix algebras.

            INPUT:

            - ``n`` (optional integer, default 0): The number of generators
              of the ideal to be created.

            OUTPUT:

            The class that is used to implement ideals of this ring with
            ``n`` generators.

            NOTE:

            Often principal ideals (``n==1``) are implemented via a different
            class.

            EXAMPLES::

                sage: MS = MatrixSpace(QQ,2,2)
                sage: MS._ideal_class_()
                <class 'sage.rings.noncommutative_ideals.Ideal_nc'>

            We don't know of a commutative ring in Sage that does not inherit
            from the base class of rings. So, we need to cheat in the next
            example::

                sage: super(Ring,QQ)._ideal_class_.__module__
                'sage.categories.rings'
                sage: super(Ring,QQ)._ideal_class_()
                <class 'sage.rings.ideal.Ideal_generic'>
                sage: super(Ring,QQ)._ideal_class_(1)
                <class 'sage.rings.ideal.Ideal_principal'>
                sage: super(Ring,QQ)._ideal_class_(2)
                <class 'sage.rings.ideal.Ideal_generic'>

            """
            from sage.rings.noncommutative_ideals import Ideal_nc
            try:
                if not self.is_commutative():
                    return Ideal_nc
            except (NotImplementedError,AttributeError):
                return Ideal_nc
            from sage.rings.ideal import Ideal_generic, Ideal_principal
            if n == 1:
                return Ideal_principal
            else:
                return Ideal_generic

        ##
        # Quotient rings
        # Again, this is defined in sage.rings.ring.pyx
        def quotient(self, I, names=None):
            """
            Quotient of a ring by a two-sided ideal.

            INPUT:

            - ``I``: A twosided ideal of this ring.
            - ``names``: a list of strings to be used as names
              for the variables in the quotient ring.

            EXAMPLES:

            Usually, a ring inherits a method :meth:`sage.rings.ring.Ring.quotient`.
            So, we need a bit of effort to make the following example work with the
            category framework::

                sage: F.<x,y,z> = FreeAlgebra(QQ)
                sage: from sage.rings.noncommutative_ideals import Ideal_nc
                sage: class PowerIdeal(Ideal_nc):
                ...    def __init__(self, R, n):
                ...        self._power = n
                ...        self._power = n
                ...        Ideal_nc.__init__(self,R,[R.prod(m) for m in CartesianProduct(*[R.gens()]*n)])
                ...    def reduce(self,x):
                ...        R = self.ring()
                ...        return add([c*R(m) for c,m in x if len(m)<self._power],R(0))
                ...
                sage: I = PowerIdeal(F,3)
                sage: Q = Rings().parent_class.quotient(F,I); Q
                Quotient of Free Algebra on 3 generators (x, y, z) over Rational Field by the ideal (x^3, x^2*y, x^2*z, x*y*x, x*y^2, x*y*z, x*z*x, x*z*y, x*z^2, y*x^2, y*x*y, y*x*z, y^2*x, y^3, y^2*z, y*z*x, y*z*y, y*z^2, z*x^2, z*x*y, z*x*z, z*y*x, z*y^2, z*y*z, z^2*x, z^2*y, z^3)
                sage: Q.0
                xbar
                sage: Q.1
                ybar
                sage: Q.2
                zbar
                sage: Q.0*Q.1
                xbar*ybar
                sage: Q.0*Q.1*Q.0
                0

            """
            from sage.rings.quotient_ring import QuotientRing
            return QuotientRing(self, I, names=names)
        def quo(self, I, names=None):
            """
            Quotient of a ring by a two-sided ideal.

            NOTE:

            This is a synonyme for :meth:`quotient`.

            EXAMPLE::

                sage: MS = MatrixSpace(QQ,2)
                sage: MS.full_category_initialisation()
                sage: I = MS*MS.gens()*MS

            ``MS`` is not an instance of :class:`~sage.rings.ring.Ring`.
            But since its category was fully initalised (which is not
            by default, by :trac:`11900`), it is an instance of
            the parent class of the category of rings. The quotient
            method is inherited from there::

                sage: isinstance(MS,sage.rings.ring.Ring)
                False
                sage: isinstance(MS,Rings().parent_class)
                True
                sage: MS.quo(I,names = ['a','b','c','d'])
                Quotient of Full MatrixSpace of 2 by 2 dense matrices over Rational Field by the ideal
                (
                  [1 0]
                  [0 0],
                <BLANKLINE>
                  [0 1]
                  [0 0],
                <BLANKLINE>
                  [0 0]
                  [1 0],
                <BLANKLINE>
                  [0 0]
                  [0 1]
                )

            """
            return self.quotient(I,names=names)

        def quotient_ring(self, I, names=None):
            """
            Quotient of a ring by a two-sided ideal.

            NOTE:

            This is a synonyme for :meth:`quotient`.

            EXAMPLE::

                sage: MS = MatrixSpace(QQ,2)
                sage: I = MS*MS.gens()*MS

            ``MS`` is not an instance of :class:`~sage.rings.ring.Ring`,
            but it is an instance of the parent class of the category of
            rings. The quotient method is inherited from there::

                sage: isinstance(MS,sage.rings.ring.Ring)
                False
                sage: isinstance(MS,Rings().parent_class)
                True
                sage: MS.quotient_ring(I,names = ['a','b','c','d'])
                Quotient of Full MatrixSpace of 2 by 2 dense matrices over Rational Field by the ideal
                (
                  [1 0]
                  [0 0],
                <BLANKLINE>
                  [0 1]
                  [0 0],
                <BLANKLINE>
                  [0 0]
                  [1 0],
                <BLANKLINE>
                  [0 0]
                  [0 1]
                )

            """
            return self.quotient(I,names=names)

        def __div__(self, I):
            """
            Since assigning generator names would not work properly,
            the construction of a quotient ring using division syntax
            is not supported.

            EXAMPLE::

                sage: MS = MatrixSpace(QQ,2)
                sage: I = MS*MS.gens()*MS
                sage: MS/I
                Traceback (most recent call last):
                ...
                TypeError: Use self.quo(I) or self.quotient(I) to construct the quotient ring.
            """
            raise TypeError, "Use self.quo(I) or self.quotient(I) to construct the quotient ring."

    class ElementMethods:
        def is_unit(self):
            r"""
            Return whether this element is a unit in the ring.

            .. NOTE::

                This is a generic implementation for (non-commutative) rings
                which only works for the one element, its additive inverse, and
                the zero element.  Most rings should provide a more specialized
                implementation.

            EXAMPLES::

                sage: MS = MatrixSpace(ZZ, 2)
                sage: MS.one().is_unit()
                True
                sage: MS.zero().is_unit()
                False
                sage: MS([1,2,3,4]).is_unit()
                Traceback (most recent call last):
                ...
                NotImplementedError

            """
            if self == 1 or self == -1:
                return True
            if self == 0: # now 0 != 1
                return False
            raise NotImplementedError

    class HomCategory(HomCategory):
        pass

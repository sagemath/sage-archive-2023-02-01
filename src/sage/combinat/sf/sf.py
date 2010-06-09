from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.misc.misc import attrcall
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.all import Rings, GradedHopfAlgebrasWithBasis, ModulesWithBasis
from sage.categories.homset import Hom
from sage.combinat.partition import Partition, Partition_class, Partitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ

#from sage.combinat.sf.categories import *
import sage.combinat.sf.sfa as sfa

class SymmetricFunctions(UniqueRepresentation, Parent):
    """
    The abstract algebra of commutative symmetric functions

    We construct the abstract algebra of commutative symmetric
    functions over the rational numbers::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Sym
        Symmetric Functions over Rational Field

    Todo: add one example of non trivial computation, and then proceed
    with the detailed explanations

    Todo: expand this tutorial, merging with that of ``MuPAD-Combinat``

        sage: h = Sym.h(); e = Sym.e(); s = Sym.s(); m = Sym.m(); p = Sym.p()
        sage: ( ( h[2,1] * ( 1 + 3 * h[2,1]) ) + s[2]. antipode()) . coproduct() # todo: not implemented

    Sym is the unique free commutative graded connected
    algebra, with one generator in each degree::

        sage: Sym.category()
        Category of abstract graded hopf algebras with basis over Rational Field

    We use the Sage standard renaming idiom to get shorter outputs::

        sage: Sym.rename("Sym")
        sage: Sym
        Sym

    Sym has many representations as a concrete algebra. Each of them
    has a distinguished basis, and its elements are expanded in this
    basis. Here is the p representation::

        sage: p = Sym.powersum()
        sage: p			# todo: not implemented (do we want this to be printed as Sym.p()?)
        Sym on the p basis

    Elements of p are linear combinations of such compositions::

      sage: p.an_element()
      1/2*p[] + 3*p[1, 1, 1] + 2*p[2, 1, 1]

    The basis itself is accessible through::

        sage: p.basis()
        Lazy family (Term map from Partitions to Symmetric Function Algebra over Rational Field, Power symmetric functions as basis(i))_{i in Partitions}
        sage: p.basis().keys()
        Partitions

    To construct an element one can therefore do::

        sage: p.basis()[Partition([2,1,1])]
        p[2, 1, 1]

    As this is rather cumbersome, the following abuses of notation are
    allowed::

        sage: p[Partition([2, 1, 1])]
        p[2, 1, 1]
        sage: p[[2, 1, 1]]
        p[2, 1, 1]
        sage: p[2, 1, 1]
        p[2, 1, 1]

    or even::

        sage: p[(i for i in [2, 1, 1])]
        p[2, 1, 1]

    Badly enough, due to a limitation in Python syntax, one cannot use::
        sage: p[]       # todo: not implemented

    Please use instead::
        sage: p[[]]
        p[]

    Now, we can construct linear combinations of basis elements::

        sage: p[2,1,1] + 2 * (p[4] + p[2,1])
        2*p[2, 1] + p[2, 1, 1] + 2*p[4]

    ..topic: Algebra structure

    Let us explore the other operations of p. First, we can ask for
    the mathematical properties of p:

        sage: p.categories() # todo: not implemented
        [The category of multiplicative bases on primitive elements of Sym,
         ...
         Category of graded hopf algebras with basis over Rational Field,
         ...
         Category of graded algebras with basis over Rational Field,
         ...
         Category of graded coalgebras with basis over Rational Field, ...]

    To start with, p is a graded algebra, the grading being induced
    by the size of compositions. Due to this, the one is the basis
    element indexed by the empty composition::

        sage: p.one()          # TODO: add the other bases
        p[]

    As we have seen above, the p basis is multiplicative; that is
    multiplication is induced by linearity from the concatenation of
    compositions::

        sage: p[3,1] * p[2,1]
        p[3, 2, 1, 1]

        sage: (p.one() + 2 * p[3,1]) * p[4, 2]
        p[4, 2] + 2*p[4, 3, 2, 1]

    ..topic: Hopf algebra structure

    p is further endowed with a coalgebra algebra structure (in
    fact, it is, up to isomorphism, the unique free algebra on
    primitive elements). The coproduct is an algebra morphism, and
    therefore determined by its values on the generators; those are
    primitive::

        sage: p[1].coproduct()		# todo: not implemented
        p[] # p[1] + p[1] # p[]
        sage: p[2].coproduct()		# todo: not implemented
        p[] # p[2] + p[2] # p[]

    The coproduct being cocommutative on the generators is cocommutative everywhere::
        sage: p[2, 1].coproduct()	# todo: not implemented
        p[] # p[2, 1] + p[1] # p[2] + p[2, 1] # p[] + p[2] # p[1]

    The antipode is an anti-algebra morphism (Todo: explain what it
    is); on the p basis, it sends the generators to their opposite:

        sage: p[3].antipode()		# todo: not implemented
        -p[3]
        sage: p[1,3,2].antipode()	# todo: not implemented
        -p[2, 3, 1]

    ..topic: Other concrete representations

    Todo: demonstrate how to customize the basis names

    Sym admits many other concrete representations::

        sage: s = Sym.schur()
        sage: h = Sym.complete()
        sage: e = Sym.elementary()
        sage: m = Sym.monomial()
        sage: f = Sym.forgotten() 	# todo: not implemented

    To change from one basis to another, one simply does:

        sage: m(p[3])
        m[3]
        sage: m(p[3,2])
        m[3, 2] + m[5]

    In general, one can mix up different basis in computations:

        sage: p( m[1] * ( e[3]*s[2] + 1 ))
        p[1] + 1/12*p[1, 1, 1, 1, 1, 1] - 1/6*p[2, 1, 1, 1, 1] - 1/4*p[2, 2, 1, 1] + 1/6*p[3, 1, 1, 1] + 1/6*p[3, 2, 1]

    Jack polynomials can be obtained as::

        sage: Sym = SymmetricFunctions(QQ['t'])
        sage: Jack = Sym.jack_polynomials()		# todo: not implemented
        sage: P = Jack.P(); J = Jack.J(); Q = Jack.Q()  # todo: not implemented
        sage: J(P[2,1])                                 # todo: not implemented

    t can be specialized as follow::

        sage: Sym = SymmetricFunctions(QQ)
        sage: Jack = Sym.jack_polynomials(t = 1)	# todo: not implemented
        sage: P = Jack.P(); J = Jack.J(); Q = Jack.Q()  # todo: not implemented
        sage: J(P[2,1])                                 # todo: not implemented

    Todo: introduce a field with degree 1 elements as in
    MuPAD-Combinat, to get proper plethysm.

    Similarly one can get Hall-Littlewood, Macdonald polynomials, etc::

        sage: Sym = SymmetricFunctions(QQ['q','t'])
        sage: Mcd = Sym.macdonald_polynomials()		# todo: not implemented
        sage: P = Mcd.P(); J = Mcd.J(); Q = Mcd.Q()     # todo: not implemented
        sage: J(P[2,1])                                 # todo: not implemented

    Further things to do:
     - Use UniqueRepresentation to get rid of all the manual cache handling for the bases
     - Devise a mechanism so that pickling bases of symmetric functions pickles
       the coercions which have a cache.

    """

    def __init__(self, R):
        """
        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)

        TESTS:

        There are a lot of missing features for this abstract parent. But some tests do pass::

            sage: TestSuite(Sym).run()
            Failure ...
            The following tests failed: _test_additive_associativity, _test_an_element, _test_associativity, _test_distributivity, _test_elements, _test_elements_eq, _test_not_implemented_methods, _test_one, _test_prod, _test_some_elements, _test_zero
        """
        assert(R in Rings())
        self._base = R # Won't be needed when CategoryObject won't override anymore base_ring
        Parent.__init__(self, category = GradedHopfAlgebrasWithBasis(R).abstract_category())

    def _repr_(self): # could be taken care of by the category
        """
        TESTS::

            sage: SymmetricFunctions(RR) # indirect doctest
            Symmetric Functions over Real Field with 53 bits of precision
        """
        return "Symmetric Functions over %s"%self.base_ring()

    # For Jason: all the functions below should call directly the corresponding classes:
    # 		return Schur(self))
    # Then, SymmetricFunctionAlgebra will just be a backward compatibility alias
    def schur(self):
        r"""
        EXAMPLES::

            sage: SymmetricFunctions(QQ).schur()
            Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis
        """
        return sfa.cache_s(self.base_ring())
    s = schur

    def powersum(self):
        r"""
        EXAMPLES::

            sage: SymmetricFunctions(QQ).powersum()
            Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
        """
        return sfa.cache_p(self.base_ring())
    p = powersum
    power = powersum # Todo: get rid of this one when it won't be needed anymore

    def complete(self):
        r"""
        EXAMPLES::

            sage: SymmetricFunctions(QQ).complete()
            Symmetric Function Algebra over Rational Field, Homogeneous symmetric functions as basis
        """
        return sfa.cache_h(self.base_ring())
    h = complete
    homogeneous = complete

    def elementary(self):
        r"""
        EXAMPLES::

            sage: SymmetricFunctions(QQ).elementary()
            Symmetric Function Algebra over Rational Field, Elementary symmetric functions as basis
        """
        return sfa.cache_e(self.base_ring())
    e = elementary

    def monomial(self):
        r"""
        EXAMPLES::

            sage: SymmetricFunctions(QQ).monomial()
            Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis
        """
        return sfa.cache_m(self.base_ring())
    m = monomial

    def from_polynomial(self, f):
        """
        This function converts a symmetric polynomial `f` in a polynomial ring in finitely
        many variables to a symmetric function in the monomial
        basis of the ring of symmetric functions over the same base ring.

        EXAMPLES::

            sage: P = PolynomialRing(QQ, 'x', 3)
            sage: x= P.gens()
            sage: f = x[0] + x[1] + x[2]
            sage: S = SymmetricFunctions(QQ)
            sage: S.from_polynomial(f)
            m[1]
        """
        return self.m().from_polynomial(f)

    def register_isomorphism(self, morphism):
        """
        Registers an isomorphism between two bases of self, as a canonical coercion

        EXAMPLES:

        We override the canonical coercion from the Schur basis to the
        powersum basis by a (stupid!) map `s_\lambda\mapsto 2p_\lambda`.

            sage: Sym = SymmetricFunctions(QQ['zorglub']) # make sure we are not going to screw up later tests
            sage: s = Sym.s(); p = Sym.p().dual_basis()
            sage: phi = s.module_morphism(diagonal = lambda t: 2, codomain = p)
            sage: phi(s[2, 1])
            2*d_p[2, 1]
            sage: Sym.register_isomorphism(phi)
            sage: p(s[2,1])
            2*d_p[2, 1]

        The map is supposed to implement the canonical isomorphism
        between the two basis. Otherwise, the results will be
        mathematically wrong, as above. Use with care!
        """
        morphism.codomain().register_coercion(morphism)

    _shorthands = set(['e', 'h', 'm', 'p', 's'])

    def inject_shorthands(self, shorthands = _shorthands):
        """
        Imports standard shorthands into the global namespace

        INPUT:

         - shorthands - a list (or iterable) of strings (default: ['e', 'h', 'm', 'p', 's'])

        EXAMPLES::

            sage: S = SymmetricFunctions(ZZ)
            sage: S.inject_shorthands()
            doctest:...: RuntimeWarning: redefining global value `e`
            doctest:...: RuntimeWarning: redefining global value `m`
            sage: s[1] + e[2] * p[1,1] + 2*h[3] + m[2,1]
            s[1] - 2*s[1, 1, 1] + s[1, 1, 1, 1] + s[2, 1] + 2*s[2, 1, 1] + s[2, 2] + 2*s[3] + s[3, 1]
            sage: e
            Symmetric Function Algebra over Integer Ring, Elementary symmetric functions as basis
            sage: p
            Symmetric Function Algebra over Integer Ring, Power symmetric functions as basis
            sage: s
            Symmetric Function Algebra over Integer Ring, Schur symmetric functions as basis

            sage: e == S.e(), h == S.h(), m == S.m(), p == S.p(), s == S.s()
            (True, True, True, True, True)

        One can also just import a subset of the shorthands::

            sage: S = SymmetricFunctions(QQ)
            sage: S.inject_shorthands(['p', 's'])
            doctest:...: RuntimeWarning: redefining global value `p`
            doctest:...: RuntimeWarning: redefining global value `s`
            sage: p
            Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
            sage: s
            Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis

        Note that ``e`` is left unchanged::

            sage: e
            Symmetric Function Algebra over Integer Ring, Elementary symmetric functions as basis
        """
        from sage.misc.misc import inject_variable
        for shorthand in shorthands:
            assert shorthand in self._shorthands
            inject_variable(shorthand, getattr(self, shorthand)())

    def __init_extra__(self):
        """
        Sets up the coercions between the different bases

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ) # indirect doctest
            sage: s = Sym.s(); p = Sym.p()
            sage: s.coerce_map_from(p)
            Generic morphism:
              From: Symmetric Function Algebra over Rational Field, Power symmetric functions as basis
              To:   Symmetric Function Algebra over Rational Field, Schur symmetric functions as basis
        """
        #category   = Bases(self)
        #category   = GradedHopfAlgebrasWithBasis(self.base_ring())
        category   = ModulesWithBasis(self.base_ring())
        powersum   = self.powersum  ()
        complete   = self.complete  ()
        elementary = self.elementary()
        schur      = self.schur     ()
        monomial   = self.monomial  ()

        iso = self.register_isomorphism

        from sage.combinat.sf.classical import conversion_functions
        from sage.categories.morphism import SetMorphism
        from sage.categories.homset import Hom

        for (basis1_name, basis2_name) in conversion_functions.keys():
            basis1 = getattr(self, basis1_name)()
            basis2 = getattr(self, basis2_name)()
            on_basis = SymmetricaConversionOnBasis(t = conversion_functions[basis1_name,basis2_name], domain = basis1, codomain = basis2)
            iso(basis1._module_morphism(on_basis, codomain = basis2, category = category))

        # Todo: fill in with other conversion functions on the classical bases

class SymmetricaConversionOnBasis:
    def __init__(self, t, domain, codomain):
        """
        INPUT:

         - t -- a function taking a monomial in CombinatorialFreeModule(QQ, Partitions()), and returning a
           (partition, coefficient) list.

         - ``domain``, ``codomain`` -- parents

        Construct a function mapping a partition to an element of ``codomain``.

        This is a temporary quick hack to wrap around the existing
        symmetrica conversions, without changing their specs.

        EXAMPLES:

            sage: Sym = SymmetricFunctions(QQ[x])
            sage: p = Sym.p(); s = Sym.s()
            sage: def t(x) : [(p,c)] = x; return [ (p,2*c), (p.conjugate(), c) ]
            sage: f = sage.combinat.sf.sf.SymmetricaConversionOnBasis(t, p, s)
            sage: f(Partition([3,1]))
            s[2, 1, 1] + 2*s[3, 1]
        """
        self._domain = domain
        self.fake_sym = CombinatorialFreeModule(QQ, Partitions())
        self._codomain = codomain
        self._t = t

    def __call__(self, partition):
        """
            sage: Sym = SymmetricFunctions(QQ[x])
            sage: p = Sym.p(); s = Sym.s()
            sage: p[1] + s[1]				# indirect doctest
            2*p[1]
        """
        R = self._codomain.base_ring()
        # TODO: use self._codomain.sum_of_monomials, when the later
        # will have an optional optimization for the case when there
        # is no repetition in the support
        return self._codomain._from_dict(dict(self._t(self.fake_sym.monomial(partition))), coerce = True)


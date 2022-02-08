r"""
Finitely presented modules over the Steenrod algebra

This package allows the user to define finitely presented modules
over the mod `p` Steenrod Algebra, elements of them, and
morphisms between them.  Methods are provided for doing basic homological
algebra, e.g. computing kernels and images of homomorphisms, and finding
free resolutions of modules.

======================
Theoretical background
======================

The category of finitely presented graded modules over an arbitrary non-Noetherian
graded ring `R` is not Abelian in general, since kernels of homomorphisms are
not necessarily finitely presented.

The mod `p` Steenrod algebra is non-Noetherian, but it is the
union of a countable set of finite sub-Hopf algebras
([Mar1983]_ Ch. 15, Sect. 1, Prop 7), and is therefore an example of a
`P`-algebra ([Mar1983]_ Ch. 13).

Any finitely presented module over the Steenrod algebra is
defined over one of these finite sub-Hopf algebras.  Similarly, any homomorphism
between finitely presented modules over the Steenrod algebra is
defined over a finite sub-Hopf algebra of the Steenrod algebra.
Further, tensoring up from the category of modules over a sub-Hopf
algebra is an exact functor, since the Steenrod algebra is free over
any sub-Hopf algebra.

It follows that kernels, cokernels, images, and, more generally, any finite
limits or colimits can be computed in the category of modules over the
Steenrod algebra, by working in the category of modules over an appropriate
finite sub-Hopf algebra.

It is also the case that presentations of modules and the images of the
generators of the domain of a homomorphism are the same over the sub-Hopf
algebra and over the whole Steenrod algebra, so that the tensoring up is
entirely implicit and requires no computation.

The definitions and computations carried out by this package are thus done
as follows.   First, find a small finite sub-Hopf algebra over which the
computation can be done.   Then, carry out the calculation there, where it
is a finite problem, and can be reduced to linear algebra over a finite
prime field.

==========
User guide
==========

Let `p` be a prime number.  The mod `p` Steenrod algebra `A`
is a connected algebra over `\GF{p}`, the finite field of `p` elements.
All modules presented here will be defined over `A`, or one of its sub-Hopf
algebras.  E.g.::

    sage: A = SteenrodAlgebra(p=2)

The constructor of the module class takes as arguments an ordered tuple of
degrees and the algebra over which the module is defined, together with an
optional set of relations::

    sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
    sage: F = SteenrodFPModule(A, [0, 1, 7]); F  # A free module
    Finitely presented left module on 3 generators and 0 relations over mod 2 Steenrod algebra, milnor basis

Denote the module generators of an `A`-module `M` by `g_1,\ldots, g_N`.
A homogeneous relation of degree `n` has the form

.. MATH::

    \sum_{i=1}^N a_i\cdot g_i = 0

where the homogeneous coefficients `a_1,\ldots a_N` in `A`, such that
`\deg(a_i) + \deg(g_i) = n` for `i=1\ldots N`.  To create a module with
relations, the coefficients for each relation is given to the constructor::

    sage: r1 = [Sq(8), Sq(7), 0]   # First relation
    sage: r2 = [Sq(7), 0, 1]       # Second relation
    sage: M = SteenrodFPModule(A, [0, 1, 7], relations=[r1, r2]); M
    Finitely presented left module on 3 generators and 2 relations over mod 2 Steenrod algebra, milnor basis

The resulting module will have three generators in the degrees we gave them::

    sage: M.generator_degrees()
    (0, 1, 7)

The connectivity of a module over a connected graded algebra is the minimum
degree of all its module generators.  Thus, if the module is non-trivial, the
connectivity is an integer::

    sage: M.connectivity()
    0

Each module is defined over a Steenrod algebra or some sub-Hopf algebra of it,
given by its base ring::

    sage: M.base_ring()
    mod 2 Steenrod algebra, milnor basis
    sage: SteenrodFPModule(SteenrodAlgebra(p=2,profile=(3,2,1)), [0]).base_ring()
    sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function
    [3, 2, 1]

.. NOTE::

    Calling :meth:`algebra` will not return the desired algebra. Users should
    use the :meth:`base_ring` method.

---------------
Module elements
---------------

Module elements are displayed in terms of generators, which by default
are called ``g[degree]]``::

    sage: M.an_element(n=5)
    Sq(2,1)*g[0] + Sq(4)*g[1]

    sage: e = M.an_element(n=15); e
    Sq(0,0,0,1)*g[0] + Sq(1,2,1)*g[1] + Sq(8)*g[7]
    sage: e.dense_coefficient_list()
    [Sq(0,0,0,1), Sq(1,2,1), Sq(8)]

The generators are themselves elements of the module::

    sage: gens = M.generators(); gens
    (g[0], g[1], g[7])
    sage: gens[0] in M
    True

Producing elements from a given set of algebra coefficients is possible using
the module class ()-method::

    sage: coeffs=[Sq(15), Sq(10)*Sq(1,1), Sq(8)]
    sage: x = M(coeffs); x
    Sq(15)*g[0] + (Sq(4,1,1)+Sq(7,0,1)+Sq(11,1))*g[1] + Sq(8)*g[7]

The module action produces new elements::

    sage: Sq(2)*x
    Sq(14,1)*g[0] + (Sq(7,1)+Sq(10))*g[7]

Each non-zero element has a well-defined degree::

    sage: x.degree()
    15

But the zero element does not::

    sage: zero = M.zero(); zero
    0
    sage: zero.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero element does not have a well-defined degree

Any two elements can be added as long as they have the same degree::

    sage: y = M.an_element(15); y
    Sq(0,0,0,1)*g[0] + Sq(1,2,1)*g[1] + Sq(8)*g[7]
    sage: x + y
    (Sq(0,0,0,1)+Sq(15))*g[0] + (Sq(1,2,1)+Sq(4,1,1)+Sq(7,0,1)+Sq(11,1))*g[1]

or when at least one of them is zero::

    sage: x + zero == x
    True

Finally, additive inverses exist::

    sage: x - x == 0
    True

At this point it may be useful to point out that elements are not reduced to
a minimal representation when they are created.  A normalization can be
enforced, however::

    sage: g3 = M([0, 0, 1]); g3
    g[7]
    sage: g3.normalize()
    Sq(7)*g[0]
    sage: g3 == g3.normalize()
    True

    sage: m = M([Sq(7), 0, 0])
    sage: sum = m + g3; sum     # m and g3 are related by `m = \operatorname{Sq}^7(g_1) = g_3`,
    Sq(7)*g[0] + g[7]
    sage: sum == 0              # so their sum should zero.
    True
    sage: sum.normalize()       # Its normalized form is more revealing.
    0

For every integer `n`, the set of module elements of degree `n` form a
vectorspace over the ground field `\GF{p}`.  A basis for this vectorspace can be
computed::

    sage: M.basis_elements(7)
    (Sq(0,0,1)*g[0],
     Sq(1,2)*g[0],
     Sq(4,1)*g[0],
     Sq(7)*g[0],
     Sq(0,2)*g[1],
     Sq(3,1)*g[1],
     Sq(6)*g[1])

Note that the third generator `g_3` of degree 7 having coordinates <0,0,1>
is apparently missing from the basis above.  This is because of the relation
`\operatorname{Sq}^7(g_1) = g_3`.

A vectorspace presentation can be produced::

    sage: M.vector_presentation(5)
    Vector space quotient V/W of dimension 4 over Finite Field of size 2 where
    V: Vector space of dimension 4 over Finite Field of size 2
    W: Vector space of degree 4 and dimension 0 over Finite Field of size 2
    Basis matrix:
    []

Given any element, its coordinates with respect to this basis can be computed::

    sage: x = M.an_element(7); x
    Sq(0,0,1)*g[0] + Sq(3,1)*g[1] + g[7]
    sage: v = x.vector_presentation(); v
    (1, 0, 0, 1, 0, 1, 0)

Going the other way, any element can be constructed by specifying its
coordinates::

    sage: x_ = M.element_from_coordinates((1, 0, 0, 1, 0, 1, 0), 7)
    sage: x_
    (Sq(0,0,1)+Sq(7))*g[0] + Sq(3,1)*g[1]
    sage: x_ == x
    True

-------------------------------------
Module homomorphisms
-------------------------------------

Homomorphisms of `A`-modules `M\to N` are linear maps of their
underlying `\GF{p}`-vectorspaces which commute with the `A`-module structure.

To create a homomorphism, first create the object modeling the set of all
such homomorphisms using the free function ``Hom``::

    sage: Hko = SteenrodFPModule(A, [0], [[Sq(2)], [Sq(1)]])
    sage: homspace = Hom(Hko, Hko); homspace
    Set of Morphisms from Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis to Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis in Category of finitely presented graded modules over mod 2 Steenrod algebra, milnor basis

Just as module elements, homomorphisms are created using the ()-method
of the homspace object.  The only argument is a list of module elements in the
codomain, corresponding to the module generators of the domain::

    sage: gen = Hko.generator(0)  # the generator of the codomain module.
    sage: values = [Sq(0, 0, 1)*gen]; values
    [Sq(0,0,1)*g[0]]
    sage: f = homspace(values)

The resulting homomorphism is the one sending the `i`-th generator of the
domain to the `i`-th codomain value given::

    sage: f
    Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> Sq(0,0,1)*g[0]

Homomorphisms can be evaluated on elements of the domain module::

    sage: v1 = f(Sq(4)*gen); v1
    Sq(4,0,1)*g[0]

    sage: v2 = f(Sq(2)*Sq(4)*gen); v2
    (Sq(3,1,1)+Sq(6,0,1))*g[0]

and they respect the module action::

    sage: v1 == Sq(4)*f(gen)
    True

    sage: v2 == Sq(2)*Sq(4)*f(gen)
    True

Convenience methods exist for creating the trivial morphism::

    sage: x = Sq(4)*Sq(7)*gen
    sage: x == 0
    False
    sage: zero_map = homspace.zero(); zero_map
    Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> 0
    sage: zero_map(x)
    0
    sage: zero_map(x).is_zero()
    True

as well as the identity endomorphism::

    sage: one = Hom(Hko, Hko).identity(); one
    Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> g[0]
    sage: one.is_endomorphism()
    True
    sage: one(x) == x
    True

Any non-trivial homomorphism has a well defined degree::

    sage: f.degree()
    7

but just as module elements, the trivial homomorphism does not::

    sage: zero_map = homspace.zero()
    sage: zero_map.degree()
    Traceback (most recent call last):
    ...
    ValueError: the zero morphism does not have a well-defined degree

Any two homomorphisms can be added as long as they are of the same degree::

    sage: f1 = homspace([Hko([Sq(0,0,3) + Sq(0,2,0,1)])])
    sage: f2 = homspace([Hko([Sq(8,2,1)])])
    sage: (f1 + f2).is_zero()
    False
    sage: f1 + f2
    Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> (Sq(0,0,3)+Sq(0,2,0,1)+Sq(8,2,1))*g[0]

or when at least one of them is zero::

    sage: f + zero_map == f
    True

Finally, additive inverses exist::

    sage: (f - f) == 0
    True

The restriction of a homomorphism to the vectorspace of `n`-dimensional module
elements is a linear transformation::

    sage: f_21 = f.vector_presentation(21); f_21
    Vector space morphism represented by the matrix:
    [1 0 0 0 0 0]
    [0 0 0 0 0 0]
    [1 0 0 0 0 0]
    Domain: Vector space quotient V/W of dimension 3 over Finite Field of size 2 where
    V: Vector space of dimension 20 over Finite Field of size 2
    W: Vector space of degree 20 and dimension 17 over Finite Field of size 2
    Basis matrix:
    [1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
    [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
    [0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1]
    [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1]
    [0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 1 0 0 0 1]
    [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1]
    [0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 1]
    [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1]
    [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1]
    Codomain: Vector space quotient V/W of dimension 6 over Finite Field of size 2 where
    V: Vector space of dimension 35 over Finite Field of size 2
    W: Vector space of degree 35 and dimension 29 over Finite Field of size 2
    Basis matrix:
    29 x 35 dense matrix over Finite Field of size 2

This is compatible with the vector presentations of its domain and codomain
modules::

    sage: f.domain() is Hko
    True
    sage: f.codomain() is Hko
    True
    sage: f_21.domain() is Hko.vector_presentation(21)
    True
    sage: f_21.codomain() is Hko.vector_presentation(21 + f.degree())
    True

Elements in the pre-image of homomorphism can be found::

    sage: f.solve(Sq(2)*Sq(4)*Sq(7)*gen)
    Sq(0,2)*g[0]

    sage: f.solve(Sq(8)*gen) is None
    True

Homomorphisms can be composed as expected::

    sage: g = homspace([Sq(0, 0, 0, 1)*gen]); g
    Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> Sq(0,0,0,1)*g[0]

    sage: g*f
    Module endomorphism of Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> Sq(0,0,1,1)*g[0]

    sage: one = homspace.identity()
    sage: f*one == f
    True

-------------------
Homological algebra
-------------------

The category of modules over `A` is Abelian, so kernels, images and cokernels
all exist and can be computed through the API belonging to the homomorphism class
:class:`sage.modules.fp_graded.steenrod.morphism.SteenrodFPModuleMorphism`.

.. NOTE:: Functions like :meth:`kernel`, :meth:`cokernel`, :meth:`image`,
    :meth:`homology` compute sub- and quotient modules related to homomorphisms,
    but they do not return instances of the module class.  Rather, they return
    the natural homomorphisms which connect these modules to the modules that
    gave rise to them.

    E.g. the function :meth:`kernel` returns an injective homomorphism
    which is onto the kernel submodule we asked it to compute.  On the other
    hand, the function :meth:`cokernel` provides a surjective homomorphism onto
    the cokernel module.

    In each case, getting a reference to the module instance requires calling
    :meth:`domain` or :meth:`codomain` on the returned homomorphism, depending
    on the case.

    Refer to each function's documentation for specific details.

.....................................
Cokernels
.....................................

In the following example, we define a cyclic module `H\ZZ` with one
relation in two ways: First explicitly, and then as the cokernel of a
homomorphism of free modules.  We then construct a candidate for an isomorphism
and check that it is both injective and surjective::

    sage: HZ = SteenrodFPModule(A, [0], [[Sq(1)]]); HZ
    Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis

    sage: F = SteenrodFPModule(A, [0])
    sage: j = Hom(F, F)([Sq(1)*F.generator(0)])
    sage: coker = j.cokernel_projection() # the natural quotient homomorphism onto the cokernel.
    sage: hz = coker.codomain(); hz
    Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis

    sage: a = Hom(HZ, hz)([hz.generator(0)])
    sage: a.is_injective()
    True
    sage: a.is_surjective()
    True

.......
Kernels
.......

When computing the kernel of a homomorphism `f`, the result is an injective homomorphism
into the domain of `f`::

    sage: k = f.kernel_inclusion(); k
    Module morphism:
      From: Finitely presented left module on 1 generator and 3 relations over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[7] |--> Sq(0,0,1)*g[0]
    sage: k.codomain() == f.domain()
    True
    sage: k.is_injective()
    True

    sage: ker = k.domain()
    sage: ker
    Finitely presented left module on 1 generator and 3 relations over mod 2 Steenrod algebra, milnor basis

We can check that the injective image of `k` is the kernel of `f` by
showing that `f` factors as `h\circ c`, where `c` is the quotient map
to the cokernel of `k`, and `h` is injective::

    sage: K = k.codomain()        # We want to check that this really is the kernel of f.
    sage: coker = k.cokernel_projection()    # coker is the natural map: Hko -> coker(f) with kernel K.
    sage: h = Hom(coker.codomain(), Hko)(f.values())
    sage: h*coker == f            # Is K contained in ker(f) ?
    True
    sage: h.is_injective()        # Is ker(f) contained in K ?
    True

......
Images
......

The method :meth:`image` behaves similarly, returning an injective
homomorphism with image equal to the submodule `\operatorname{im}(f)` ::

    sage: i = f.image(); i
    Module morphism:
      From: Finitely presented left module on 1 generator and 3 relations over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[7] |--> Sq(0,0,1)*g[0]
    sage: i.codomain() == f.codomain()
    True
    sage: i.is_injective()
    True

We can check that the injective image of `i` is the image of `f` by
lifting `f` over `i`, and showing that the lift is surjective::

    sage: f_ = f.lift(i); f_
    Module morphism:
      From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 3 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> g[7]
    sage: i*f_ == f             # Is im(i) contained in im(f) ?
    True
    sage: f_.is_surjective()    # Is im(f) contained in im(i) ?
    True

When a pair of composable homomorphisms `g\circ f: M\to N\to L` satisfy the condition
`g\circ f = 0`, the sub-quotient `\ker(g) / \operatorname{im}(f)` can be computed and
is given by the natural quotient homomorphism with domain `\ker(g)`::

    sage: f*f == 0          # Does the kernel of f contain the image of f ?
    True
    sage: K = f.kernel_inclusion()    # k: ker(f) -> Hko
    sage: h = f.homology(f) # h: ker(f) -> ker(f) / im(f)
    sage: h.codomain()      # This is the homology module.
    Finitely presented left module on 1 generator and 4 relations over mod 2 Steenrod algebra, milnor basis

................
Free resolutions
................

Finally, free resolutions can be computed.  These calculations usually take
some time to complete, so it is usually a good idea to raise the verbose flag
to output progress information.

The following examples are taken from
`Michael Catanzaro's thesis <https://digitalcommons.wayne.edu/oa_theses/602/>`_
where the first version of this software appeared::

    sage: res = Hko.resolution(6, verbose=True)
    Computing f_1 (1/6)
    Computing f_2 (2/6)
    Resolving the kernel in the range of dimensions [1, 8]: 1 2 3 4 5 6 7 8.
    Computing f_3 (3/6)
    Resolving the kernel in the range of dimensions [2, 10]: 2 3 4 5 6 7 8 9 10.
    Computing f_4 (4/6)
    Resolving the kernel in the range of dimensions [3, 13]: 3 4 5 6 7 8 9 10 11 12 13.
    Computing f_5 (5/6)
    Resolving the kernel in the range of dimensions [4, 18]: 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18.
    Computing f_6 (6/6)
    Resolving the kernel in the range of dimensions [5, 20]: 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20.

The result of the calculation is a list of all the maps in the resolution::

    sage: [f.domain() for f in res]
    [Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 3 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 4 generators over mod 2 Steenrod algebra, milnor basis,
     Free graded left module on 4 generators over mod 2 Steenrod algebra, milnor basis]

    sage: def is_complex(res):
    ....:     for i in range(len(res)-1):
    ....:         f = (res[i]*res[i+1])
    ....:         if not f.is_zero():
    ....:             return False
    ....:     return True
    ....:
    sage: is_complex(res)
    True

    sage: def is_exact(res):
    ....:     for i in range(len(res)-1):
    ....:         h = res[i].homology(res[i+1])
    ....:         if not h.codomain().is_trivial():
    ....:             return False
    ....:     return True

    sage: is_exact(res)
    True

    sage: [r.codomain().generator_degrees() for r in res]
    [(0,), (0,), (2, 1), (2, 4), (3, 7), (4, 8, 12), (5, 9, 13, 14)]

    sage: [r.values() for r in res]
    [(g[0],),
     (Sq(2)*g[0], Sq(1)*g[0]),
     (Sq(1)*g[1], Sq(3)*g[1] + Sq(2)*g[2]),
     (Sq(1)*g[2], Sq(2,1)*g[2] + Sq(3)*g[4]),
     (Sq(1)*g[3], Sq(2,1)*g[3] + Sq(1)*g[7], Sq(2,1)*g[7]),
     (Sq(1)*g[4],
      Sq(2,1)*g[4] + Sq(1)*g[8],
      Sq(2,1)*g[8] + Sq(1)*g[12],
      Sq(2)*g[12]),
     (Sq(1)*g[5],
      Sq(2,1)*g[5] + Sq(1)*g[9],
      Sq(2,1)*g[9] + Sq(1)*g[13],
      Sq(0,1)*g[13] + Sq(2)*g[14])]

-----------------------
Example: Counting lifts
-----------------------

In this more elaborate example we show how to find all possible lifts of a
particular homomorphism.  We will do this in two ways, and as a check of
validity, we will compare the results in the end.

We will work with the following modules over the mod 2
Steenrod algebra `A`:

.. MATH::

    \begin{align}
        H\ZZ &= A/A\cdot \operatorname{Sq}(1)\\
        Hko &= A/A\cdot \{\operatorname{Sq}(2), \operatorname{Sq}(1)\} \,.
    \end{align}

There is a natural projection `q: H\ZZ\to Hko`, and a non-trivial
endomorphism of degree 28, represented as a degree zero map `f: \Sigma^{28}Hko\to Hko`
which we define below.

The problem we will solve is to find all possible homomorphisms
`f': \Sigma^{28}Hko\to H\ZZ`, making the following diagram into a
commuting triangle:

.. MATH::

    H\ZZ\xrightarrow{q} Hko \xleftarrow{f} \Sigma^{28} Hko

We begin by defining the modules and the homomorphisms `f` and `q`.  In the following,
we let `L = \Sigma^{28}Hko`::

    sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
    sage: A = SteenrodAlgebra(2)
    sage: Hko = SteenrodFPModule(A, [0], [[Sq(2)],[Sq(1)]])
    sage: HZ = SteenrodFPModule(A, [0], [[Sq(1)]])
    sage: L = Hko.suspension(28)

    sage: # The projection:
    sage: q = Hom(HZ, Hko)([Hko.generator(0)])
    sage: q
    Module morphism:
      From: Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[0] |--> g[0]

    sage: # The map to lift over `q`:
    sage: f = Hom(L, Hko)([Sq(0,2,1,1)*Hko.generator(0)])
    sage: f
    Module morphism:
      From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      Defn: g[28] |--> Sq(0,2,1,1)*g[0]

    sage: f.is_zero()   # f is non-trivial.
    False

We will count the number of different lifts in two ways.  First, we will simply
compute the vectorspace of all possible maps `L \to H\ZZ`, and then check which
of those become `f` when composed with `q`::

    sage: basis = Hom(L, HZ).basis_elements(0)   # The basis for the vectorspace of degree 0 maps L -> HZ

    sage: from itertools import product

    sage: def from_coords(c):
    ....:     '''
    ....:     Create a linear combination of the three basis homomorphisms.
    ....:     '''
    ....:     return c[0]*basis[0] + c[1]*basis[1] + c[2]*basis[2]

    sage: for coords in product([0,1], repeat=3):
    ....:     print('%s: %s' % (coords, q*from_coords(coords) == f))
    (0, 0, 0): False
    (0, 0, 1): False
    (0, 1, 0): True
    (0, 1, 1): True
    (1, 0, 0): True
    (1, 0, 1): True
    (1, 1, 0): False
    (1, 1, 1): False

From this we conclude that four out of eight different homomorphisms
`L \to H\ZZ` are lifts of f::

    sage: lifts = [\
    ....:     from_coords((0,1,0)),\
    ....:     from_coords((0,1,1)),\
    ....:     from_coords((1,0,0)),\
    ....:     from_coords((1,0,1))]
    sage: lifts
    [Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> Sq(6,5,1)*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(6,5,1)+Sq(18,1,1))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> Sq(10,1,0,1)*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(10,1,0,1)+Sq(18,1,1))*g[0]]

Alternatively we can use left-exactness of the functor `\operatorname{Hom}_A(L, -)`
to enumerate all possible lifts of `f`.  Start by finding a single lift of `f`
over the projection `q`::

    sage: f_ = f.lift(q); f_
    Module morphism:
      From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
      To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
      Defn: g[28] |--> (Sq(4,3,0,1)+Sq(6,0,1,1)+Sq(7,2,0,1)+Sq(10,1,0,1))*g[0]
    sage: q*f_ == f  # Check that f_ is indeed a lift.
    True

There is an exact sequence

.. MATH::

    0\to \operatorname{Hom}_A(L, \ker(q)) \xrightarrow{iK_*} \operatorname{Hom}_A(L, H\ZZ) \xrightarrow{q_*} \operatorname{Hom}_A(L, Hko)\,,

which means that the indeterminacy of choosing a lift for
`f\in \operatorname{Hom}_A(L, Hko)` is represented by an element in
`\operatorname{Hom}_A(L,\ker(f))`.  Therefore, we can proceed to count the
number of lifts by computing this vectorspace of homomorphisms::

    sage: iK = q.kernel_inclusion()
    sage: K = iK.domain()
    sage: K.generator_degrees()
    (2,)
    sage: K.relations()
    (Sq(2)*g[2],)
    sage: ind = Hom(L, K).basis_elements(0); len(ind)
    2

So now we know that the vectorspace of indeterminacies is 2-dimensional over the
field of two elements.  This means that there are four distinct lifts of `f` over
`q`, and we can construct these by taking the one lift we already found, and add
to it all the different elements in the image of `iK_*`::

    sage: lifts_ = [f_,\
    ....:     f_ + iK*ind[0],\
    ....:     f_ + iK*ind[1],\
    ....:     f_ + iK*(ind[0] + ind[1])]
    sage: lifts_
    [Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(4,3,0,1)+Sq(6,0,1,1)+Sq(7,2,0,1)+Sq(10,1,0,1))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(0,7,1)+Sq(3,6,1)+Sq(4,1,3)+Sq(6,0,1,1)+Sq(6,5,1)+Sq(7,0,3))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(4,3,0,1)+Sq(6,0,1,1)+Sq(7,2,0,1)+Sq(10,1,0,1)+Sq(12,3,1)+Sq(15,2,1)+Sq(18,1,1))*g[0],
     Module morphism:
       From: Finitely presented left module on 1 generator and 2 relations over mod 2 Steenrod algebra, milnor basis
       To:   Finitely presented left module on 1 generator and 1 relation over mod 2 Steenrod algebra, milnor basis
       Defn: g[28] |--> (Sq(0,7,1)+Sq(3,6,1)+Sq(4,1,3)+Sq(6,0,1,1)+Sq(6,5,1)+Sq(7,0,3)+Sq(12,3,1)+Sq(15,2,1)+Sq(18,1,1))*g[0]]

As a test of correctness, we now compare the two sets of lifts.  As they stand,
it is not obvious that the lists ``lifts`` and ``lifts_`` are the same (up to a
re-ordering of list elements), so the following comparison is reassuring::

    sage: lifts_[0] == lifts[2]
    True
    sage: lifts_[1] == lifts[0]
    True
    sage: lifts_[2] == lifts[3]
    True
    sage: lifts_[3] == lifts[1]
    True

TESTS:

    sage: A = SteenrodAlgebra(2)
    sage: K = SteenrodFPModule(A, [1,3]); K
    Finitely presented left module on 2 generators ...
    sage: K.category()
    Category of finitely presented graded modules over mod 2 Steenrod algebra, milnor basis
    sage: L = SteenrodFPModule(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]]);L
    Finitely presented left module on 2 generators and 2 relations ...
    sage: M = SteenrodFPModule(A, [2,3], [[Sq(2),Sq(1)]]);M
    Finitely presented left module on 2 generators and 1 relation ...
    sage: K.element_class
    <class 'sage.modules.fp_graded.steenrod.module.SteenrodFPModule_with_category.element_class'>
    sage: m = M((0,1)); m
    g[3]
    sage: K.is_parent_of(m)
    False
    sage: L.is_parent_of(m)
    False
    sage: M.is_parent_of(m)
    True

    sage: SteenrodFPModule(ZZ, [0])
    Traceback (most recent call last):
    ...
    ValueError: the ground ring of the algebra must be a field

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

#*****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu>
#             and          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import Hom
from sage.rings.infinity import infinity
from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra_generic
from sage.modules.fp_graded.module import FPModule
from .profile import enveloping_profile_elements

class SteenrodFPModule(FPModule):
    r"""
    Create a finitely presented module over the Steenrod algebra.

    INPUT:

    - ``generator_degrees`` -- A tuple integer degrees.

    - ``algebra`` -- The Steenrod algebra over which the module is defined.

    - ``relations`` -- A tuple of relations.  A relation is a tuple of
      coefficients `(c_1, \ldots, c_n)` corresponding to the module
      generators.

    OUTPUT: The finitely presented module over ``algebra`` with
    presentation given by ``generator_degrees`` and ``relations``.

    TESTS:

        sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
        sage: SteenrodFPModule(SteenrodAlgebra(2), (0,))
        Finitely presented left module on 1 generator and 0 relations over mod 2 Steenrod algebra, milnor basis

    """
    def _Hom_(self, other, category=None):
        """
        The homset from ``self`` to ``other``

        TESTS::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: F = SteenrodFPModule(A, [1,3])
            sage: L = SteenrodFPModule(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]])

            sage: Hom(F, L)
            Set of Morphisms from Finitely presented left module on 2 generators ...

            sage: M = A.free_graded_module((2,1))
            sage: Hom(F, M)
            Set of Morphisms from Finitely presented left module on 2 generators ...
        """
        from .homspace import SteenrodFPModuleHomspace
        return SteenrodFPModuleHomspace(self, other, category=category)


    def change_ring(self, algebra):
        r"""
        Change the base ring of ``self``.

        INPUT:

        - ``algebra`` -- a connected graded algebra

        OUTPUT:

        The finitely presented module over ``algebra`` defined with the
        exact same number of generators of the same degrees and relations
        as ``self``.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
            sage: A = SteenrodAlgebra(2)
            sage: A1 = SteenrodAlgebra(2, profile=(2,1))

            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2), Sq(1)]])
            sage: N = M.change_ring(A1); N
            Finitely presented left module on 2 generators and 1 relation over
             sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [2, 1]

        Changing back yields the original module::

            sage: N.change_ring(A) is M
            True
        """
        return SteenrodFPModule(self._j.change_ring(algebra))


    def profile(self):
        r"""
        A finite profile over which this module can be defined.

        .. NOTE:: The profile produced by this function is reasonably small,
           but is not guaranteed to be minimal.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2),Sq(1)],[0,Sq(2)],[Sq(3),0]])
            sage: M.profile()
            (2, 1)

        TESTS:

            sage: from sage.modules.fp_graded.steenrod.module import *
            sage: A = SteenrodAlgebra(2)
            sage: X = SteenrodFPModule(A, [0])
            sage: X.profile()
            (1,)

        """
        elements = [coeffifient for value in self._j.values()\
                for coeffifient in value.dense_coefficient_list()]

        elements = [a for a in elements if a not in (0, 1)]

        profile = enveloping_profile_elements(elements)

        # Avoid returning the zero profile because it triggers a corner case
        # in FPModule.resolution().
        #
        # XXX todo: Fix FPModule_class.resolution().
        #
        return (1,) if profile == (0,) else profile


    def min_pres(self, verbose=False):
        r"""
        A minimal presentation of this module.

        OUTPUT: An isomorphism `M \to self`, where `M` has minimal
        presentation.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = SteenrodFPModule(A, [0,1], [[Sq(2),Sq(1)],[0,Sq(2)],[Sq(3),0]])

            sage: i = M.min_pres()
            sage: i.codomain() is M
            True

            sage: i.is_injective()
            True
            sage: i.is_surjective()
            True

            sage: i.domain().relations()
            (Sq(2)*g[0] + Sq(1)*g[1], Sq(2)*g[1])

            sage: i.codomain().relations()
            (Sq(2)*g[0] + Sq(1)*g[1], Sq(2)*g[1], Sq(3)*g[0])

        """
        return self._Hom_(self).identity().image(verbose=verbose)


    def resolution(self, k, top_dim=None, verbose=False):
        r"""
        A resolution of this module of the given length.

        INPUT:

        - ``k`` -- An non-negative integer.
        - ``verbose`` -- A boolean to control if log messages should be emitted.
          (optional, default: ``False``)

        OUTPUT: A list of homomorphisms `[\epsilon, f_1, \ldots, f_k]` such that

            `f_i: F_i \to F_{i-1}` for `1\leq i\leq k`,

            `\epsilon: F_0\to M`,

          where each `F_i` is a finitely generated free module, and the
          sequence

            `F_k \xrightarrow{f_k} F_{k-1} \xrightarrow{f_{k-1}} \ldots \rightarrow F_0 \xrightarrow{\epsilon} M \rightarrow 0`

          is exact.

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import *
            sage: A = SteenrodAlgebra(2)
            sage: Hko = SteenrodFPModule(A, [0], [[Sq(1)], [Sq(2)]])

            sage: res = Hko.resolution(5, verbose=True)
            Computing f_1 (1/5)
            Computing f_2 (2/5)
            Resolving the kernel in the range of dimensions [1, 8]: 1 2 3 4 5 6 7 8.
            Computing f_3 (3/5)
            Resolving the kernel in the range of dimensions [2, 10]: 2 3 4 5 6 7 8 9 10.
            Computing f_4 (4/5)
            Resolving the kernel in the range of dimensions [3, 13]: 3 4 5 6 7 8 9 10 11 12 13.
            Computing f_5 (5/5)
            Resolving the kernel in the range of dimensions [4, 18]: 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18.
            sage: [x.domain() for x in res]
            [Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 3 generators over mod 2 Steenrod algebra, milnor basis,
             Free graded left module on 4 generators over mod 2 Steenrod algebra, milnor basis]

        When there are no relations, the resolution is trivial::

            sage: M = SteenrodFPModule(A, [0])
            sage: M.resolution(4)
            [Module morphism:
               From: Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis
               To:   Finitely presented left module on 1 generator and 0 relations over mod 2 Steenrod algebra, milnor basis
               Defn: g[0] |--> g[0],
             Module morphism:
               From: Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis
               To:   Free graded left module on 1 generator over mod 2 Steenrod algebra, milnor basis,
             Module endomorphism of Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis,
             Module endomorphism of Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis,
             Module endomorphism of Free graded left module on 0 generators over mod 2 Steenrod algebra, milnor basis]
        """
        algebra = self.base_ring()
        finite_algebra = SteenrodAlgebra_generic(algebra.prime(), profile=self.profile())

        # Change rings to the finite algebra, and call the base class
        # implementation of this function.
        res = FPModule.resolution(
            self.change_ring(finite_algebra),
            k,
            top_dim=top_dim,
            verbose=verbose)

        # Change rings back to the original Steenrod algebra.
        return [j.change_ring(self.base_ring()) for j in res]


    def export_module_definition(self, powers_of_two_only=True):
        r"""
        Export the module to the input
        `format used by R. Bruner's Ext software <http://www.math.wayne.edu/~rrb/cohom/modfmt.html>`_.

        INPUT:

        - ``powers_of_two_only`` -- A boolean to control if the output should
          contain the action of all Steenrod squaring operations (restricted
          by the profile), or only the action of the operations of degree equal
          to a power of two. (optional, default: ``True``)

        EXAMPLES::

            sage: from sage.modules.fp_graded.steenrod.module import *
            sage: A1 = algebra=SteenrodAlgebra(p=2, profile=[2,1])
            sage: M = SteenrodFPModule(A1, [0])
            sage: M.export_module_definition()
            8 0 1 2 3 3 4 5 6
            0 1 1 1
            2 1 1 4
            3 1 1 5
            6 1 1 7
            0 2 1 2
            1 2 2 3 4
            2 2 1 5
            3 2 1 6
            4 2 1 6
            5 2 1 7
            sage: N = SteenrodFPModule(A1, [0], [[Sq(1)]])
            sage: N.export_module_definition()
            4 0 2 3 5
            1 1 1 2
            0 2 1 1
            2 2 1 3
            sage: N.export_module_definition(powers_of_two_only=False)
            4 0 2 3 5
            1 1 1 2
            0 2 1 1
            2 2 1 3
            0 3 1 2
            sage: A2 = SteenrodAlgebra(p=2, profile=[3,2,1])
            sage: Hko = SteenrodFPModule(A2, [0], [[Sq(1)], [Sq(2)]])
            sage: Hko.export_module_definition()
            8 0 4 6 7 10 11 13 17
            2 1 1 3
            4 1 1 5
            1 2 1 2
            5 2 1 6
            0 4 1 1
            2 4 1 4
            3 4 1 5
            6 4 1 7

        """
        if not self.base_ring().is_finite():
            raise (ValueError, 'this module is not defined over a finite algebra')
            return

        if self.base_ring().characteristic() != 2:
            raise (ValueError, 'this function is not implemented for odd primes')
            return

        n = self.connectivity()
        if n == infinity:
            print('The module connectivity is infinite, so there is ' +
                  'nothing to export.')
            return

        limit = self.base_ring().top_class().degree() + max(self.generator_degrees())

        # Create a list of bases, one for every module degree we consider.
        vector_space_basis = [self.basis_elements(i) for i in range(n, limit+1)]
        # print (vector_space_basis)

        additive_generator_degrees = []
        additive_generator_global_indices = [0]
        for dim, basis_vectors in enumerate(vector_space_basis):
            additive_generator_global_indices.append(
                len(basis_vectors) + additive_generator_global_indices[-1])
            additive_generator_degrees += len(basis_vectors)*[dim + n]

        # Print the degrees of the additive generators.
        print('%d %s' % (
            len(additive_generator_degrees),
            ' '.join(['%d' % x for x in additive_generator_degrees])))

        # A private function which transforms a vector in a given dimension
        # to a vector of global indices for the basis elements corresponding
        # to the non-zero entries in the vector.  E.g.
        # _GetIndices(dim=2, vec=(1,0,1)) will return a vector of length two,
        # (a, b), where a is the index of the first vector in the basis for
        # the 2-dimensional part of the module, and b is the index of the
        # last vector in the same part.
        def _GetIndices(dim, vec):
            if len(vector_space_basis[dim]) != len(vec):
                raise ValueError('the given vector\n%s\nhas the wrong size, it should be %d' % (str(vec), len(vector_space_basis[dim])))
            base_index = additive_generator_global_indices[dim]
            return [base_index + a for a,c in enumerate(vec) if c != 0]

        profile = self.base_ring()._profile
        powers = [2**i for i in range(profile[0])] if powers_of_two_only else\
            range(1, 2**profile[0])

        for k in powers:
            images = [[(self.base_ring().Sq(k)*x).vector_presentation() for x in D]\
                      for D in vector_space_basis]

            element_index = 0

            # Note that the variable dim is relative to the bottom dimension, n.
            for dim, image in enumerate(images):
                for im in image:
                    if im != 0 and im != None:
                        values = _GetIndices(dim + k, im)

                        print ("%d %d %d %s" % (
                            element_index,
                            k,
                            len(values),
                            " ".join(["%d" % x for x in values])))
                    element_index += 1

r"""
Combinatorics quickref
----------------------

Integer Sequences::

    sage: oeis([1,3,19,211])
    0: A000275: Coefficients of a Bessel function (reciprocal of J_0(z)); also pairs of permutations with rise/rise forbidden.
    sage: s.programs()
    0: (PARI) a(n)=if(n<0,0,n!^2*4^n*polcoeff(1/besselj(0,x+x*O(x^(2*n))),2*n)) /* Michael Somos May 17 2004 */

Combinatorial objects::

    sage: P = Partitions(10); P.cardinality(); P.<tab>
    sage: C = Combinations([1,3,7]); C.list()
    sage: Compositions(5, max_part = 3).unrank(3)

Constructions and Species::

    sage: for (p, c) in CartesianProduct(P, C): print p, c
    sage: DisjointUnionEnumeratedSets(Family(lambda n: IntegerVectors(n, 3), NonNegativeIntegers))

Words::

    sage: W=Words('abc')
    sage: Word('aabca').some_flashy_feature()

Polytopes::

    sage: L = LatticePolytope(random_matrix(ZZ, 3,6, x=7))
    sage: L.npoints(); L.plot3d()

Root systems, Coxeter and Weyl groups (:mod:`sage.combinat.root_system`)::

Crystals (:mod:`sage.combinat.crystals`)::

    sage: CrystalOfTableaux(["A",3], shape = [3,2])

Symmetric functions and combinatorial Hopf algebras (:mod:`sage.combinat.algebraic_combinatorics`)::

    sage: Sym = SymmetricFunctions(QQ); Sym.inject_shorthands()
    sage: s[3] * h[2] ...

Discrete groups, Permutation groups (:mod:`sage.groups`)::

    sage: S = SymmetricGroup(4)
    sage: M = MultivariatePolynomials('x0,x1,x2,x3')
    sage: M(...).action??? S.

Graph theory, posets, lattices (:class:`Graph`, :class:`Digraph`, :mod:`sage.combinat.posets`)::

    sage: Poset({1: [2,3], 2: [4], 3: [4]}).some_snappy_feature()
"""

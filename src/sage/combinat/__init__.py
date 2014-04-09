__doc__ = r"""
Combinatorics quickref
----------------------

 - :mod:`sage.combinat.demo`
 - sage.combinat.demo_short
 - sage.combinat.combinat?  (pretty outdated)
 - sage.combinat.root_systems?
 - sage.combinat.species?

See also:
 - :class:`EnumeratedSets`, :class:`FiniteEnumeratedSets`

* Integer Sequences

sloane_find(list), sloane.<tab>
s = sloane.find([1,3,19,211])[0]
s = sloane.i_am_lucky([1,3,19,211])
s(5), s.references()

* Combinatorial objects:

P = Partitions(10); P.count(); P.<tab>
C = Combinations([1,3,7]); C.list()
Compositions(5, max_part = 3, ...).unrank(3)
Tableau

* Constructions and Species

for (p, c) in CartesianProduct(P, C): ...

DisjointUnion(Family(lambda n: IntegerVectors(n, 3), NonNegativeIntegers))

* Words

W=Words('abc') W('aabca')

Franco: what about putting W('aabca').bar(), where bar would be some flashy feature?

* Posets

Posets: Poset([[1,2],[4],[3],[4],[]])

* Polytopes

L =LatticePolytope(random_matrix(ZZ, 3,6, x=7))
L.npoints() L.plot3d()

* Root systems, Coxeter and Weyl groups

See: sage.combinat.root_system?

* Crystals

CrystalOfTableaux(["A",3], shape = [3,2])

See sage.combinat.crystals?

* Symmetric functions and combinatorial Hopf algebras

 Sym = SymmetricFunctions(QQ)
 %from Sym.shortcuts() import *   /   %from Sym.shortcuts() import s, h, m
 Sym.import_shortcuts()           /   Sym.import_shortcuts("s,h,m")
 s[3] * h[2] ...

 NCSF
 QSym
 MultivariatePolynomials

 SymmetricGroupAlgebra

 HeckeAlgebra

* Discrete groups, Permutation groups

See sage.groups?

  S = SymmetricGroup(4)
  M = MultivariatePolynomials('x0,x1,x2,x3')
  M(...).action??? S.

* Lattices

* Graph theory and posets

See Graph?, Digraph?, graphs?

Poset({1: [2,3], 2: [4], 3: [4]}).some_snappy_feature()

"""

import demo
import demo_short
import demo_algebraic_combinatorics
import tutorial_enumerated_sets
import tutorial

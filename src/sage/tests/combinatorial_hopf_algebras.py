r"""
Tests For Combinatorial Hopf Algebras

There are various morphisms between the Hopf algebras below.
We test that the diagram of morphisms is commutative::

    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
    sage: QSym = QuasiSymmetricFunctions(QQ)
    sage: FQSym = algebras.FQSym(QQ)
    sage: FSym = algebras.FSym(QQ)
    sage: FSymDual = FSym.dual()
    sage: Sym = SymmetricFunctions(QQ)

    sage: def go(composition):
    ....:     x = NSym.a_realization()[composition]
    ....:     if QSym(FQSym(x)) != QSym(Sym(x)): return False
    ....:     if Sym(FSym(x)) != Sym(x): return False
    ....:     if FQSym(FSym(x)) != FQSym(x): return False
    ....:     if FSymDual(Sym(x)) != FSymDual(FQSym(x)): return False
    ....:     return True

    sage: go([2,1,2])  # not tested (needs more morphisms)
    True
    sage: all(all(go(comp) for comp in Compositions(n)) for n in range(5))  # not tested (needs more morphisms)
    True

    sage: def go2(n):
    ....:     for sigma in Permutations(n):
    ....:          x = FQSym.F()[sigma]
    ....:          if QSym(FSymDual(x)) != QSym(x): return False
    ....:     s = Sym.s()
    ....:     for mu in Partitions(n):
    ....:          x = s[mu]
    ....:          if QSym(FSymDual(x)) != QSym(x): return False
    ....:     return True

    sage: all(go2(n) for n in range(6))  # not tested (needs more morphisms)
    True

.. TODO::

    Add tests checking compatibility of the Hopf structure.
"""


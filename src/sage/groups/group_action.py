class FiniteGroupAction(SageObject):

    def __init__(self, A, X=None, G=None):
        r"""
        A finite group action.


        INPUT:

        - ``A`` -- can be one of the following:

          * a permutation group, if ``X`` and ``G`` are ``None``

          * a bijection on a finite set ``X``, if ``G`` is ``None``

          * a map `G\times X\to X`, otherwise.

        - ``X`` (optional, default ``None``) -- a finite set

        - ``G`` (optional, default ``None``) -- a finite group

        EXAMPLES::

            sage: G = PermutationGroup([['b','c','a']], domain=['a','b','c','d'])
            sage: a = FiniteGroupAction(G)
            sage: a.orbits()
            [['a', 'b', 'c'], ['d']]

            sage: A = lambda x: (2*x) % 6
            sage: X = [0,1,2,3,4,5]
            sage: a = FiniteGroupAction(A, X)
            sage: a.orbits()
            [[0], [1, 2, 4], [3], [5]]

            sage: A = lambda g, x: g*x
            sage: G = SymmetricGroup(3)
            sage: a = FiniteGroupAction(A, G, G)
            sage: a.orbits()
            [[1, 2, 3]]

            sage: A = lambda g, x: vector(g*x, immutable=True)
            sage: G = SL(2,3)
            sage: X = [vector(x, immutable=True) for x in GF(3)^2]
            sage: a = FiniteGroupAction(A, X, G)
            sage: a.orbits()
            [[(0, 0)], [(1, 0), (0, 2), (2, 2), (2, 0), (1, 2), (2, 1), (0, 1), (1, 1)]]

            sage: G = graphs.ButterflyGraph().automorphism_group()
            sage: a = FiniteGroupAction(G)
            sage: a.orbits() 
            [[0, 1, 2, 3], [4]]
        """
        if X is None and G is None:
            assert A in PermutationGroups
            self._G = A
            self._X = A.domain()

        elif G is None:
            # make this lazy!
            self._X = X
            from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
            gens = [tuple(o) for o in orbit_decomposition(X, A)]
            self._G = PermutationGroup(gens, domain = X)

        else:
            # make this lazy!            
            self._X = X
            from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition            
            gens = [[tuple(o) for o in orbit_decomposition(X, lambda x: A(g, x))]
                    for g in G.gens()]
            self._G = PermutationGroup(gens, domain = X)

    def orbits(self):
        return self._G.orbits()

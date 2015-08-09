r"""
Sets of morphisms between differentiable manifolds

The class :class:`DiffManifoldHomset` implements sets of morphisms between
two differentiable manifolds over the same topological field `K` (in most
applications, `K = \RR` or `K = \CC`), a morphism being a *differentiable map*
for the category of differentiable manifolds.

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

.. [1] J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
   (New York) (2013)
.. [2] S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
   Interscience Publishers (New York) (1963)

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.manifolds.manifold_homset import TopManifoldHomset
from sage.manifolds.differentiable.diff_map import DiffMap

class DiffManifoldHomset(TopManifoldHomset):
    r"""
    Set of differentiable maps between two differentiable manifolds.

    Given two differentiable manifolds `M` and `N` over a topological field `K`,
    the class :class:`DiffManifoldHomset` implements the set
    `\mathrm{Hom}(U,V)` of morphisms (i.e. differentiable maps)
    `U\rightarrow V`, where `U` is an open subset of `M` and `V` an open
    subset of `N`. Note that, as open subsets of differentiable manifolds, `U`
    and `V` are differentiable manifolds by themselves.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.manifolds.differentiable.diff_map.DiffMap`.

    INPUT:

    - ``domain`` -- open subset `U\subset M` (domain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.differentiable.manifold.DiffManifold`
    - ``codomain`` -- open subset `V\subset N` (codomain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.differentiable.manifold.DiffManifold`
    - ``name`` -- (default: ``None``) string; name given to the homset; if
      none is provided, Hom(U,V) will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      homset; if none is provided, `\mathrm{Hom}(U,V)` will be used

    EXAMPLES:

    Set of differentiable maps between a 2-dimensional differentiable manifold
    and a 3-dimensional one::

        sage: M = DiffManifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: N = DiffManifold(3, 'N')
        sage: Y.<u,v,w> = N.chart()
        sage: H = Hom(M, N) ; H
        Set of Morphisms from 2-dimensional differentiable manifold M to
         3-dimensional differentiable manifold N in Category of sets
        sage: type(H)
        <class 'sage.manifolds.differentiable.manifold_homset.DiffManifoldHomset_with_category'>
        sage: H.category()
        Category of homsets of sets
        sage: latex(H)
        \mathrm{Hom}\left(M,N\right)
        sage: H.domain()
        2-dimensional differentiable manifold M
        sage: H.codomain()
        3-dimensional differentiable manifold N

    An element of ``H`` is a differentiable map from ``M`` to ``N``::

        sage: H.Element
        <class 'sage.manifolds.differentiable.diff_map.DiffMap'>
        sage: f = H.an_element() ; f
        Differentiable map from the 2-dimensional differentiable manifold M to the
         3-dimensional differentiable manifold N
        sage: f.display()
        M --> N
           (x, y) |--> (u, v, w) = (0, 0, 0)

    The test suite is passed::

        sage: TestSuite(H).run()

    When the codomain coincides with the domain, the homset is a set of
    *endomorphisms* in the category of differentiable manifolds::

        sage: E = Hom(M, M) ; E
        Set of Morphisms from 2-dimensional differentiable manifold M to
         2-dimensional differentiable manifold M in Category of sets
        sage: E.category()
        Category of endsets of sets
        sage: E.is_endomorphism_set()
        True
        sage: E is End(M)
        True

    In this case, the homset is a monoid for the law of morphism composition::

        sage: E in Monoids()
        True

    This was of course not the case for ``H = Hom(M, N)``::

        sage: H in Monoids()
        False

    The identity element of the monoid is of course the identity map of ``M``::

        sage: E.one()
        Identity map Id_M of the 2-dimensional differentiable manifold M
        sage: E.one() is M.identity_map()
        True
        sage: E.one().display()
        Id_M: M --> M
           (x, y) |--> (x, y)

    The test suite is passed by ``E``::

        sage: TestSuite(E).run()

    This test suite includes more tests than in the case of ``H``, since ``E``
    has some extra structure (monoid)::

        sage: TestSuite(H).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass

    ::

        sage: TestSuite(E).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass

    """

    Element = DiffMap

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        TESTS::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = DiffManifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M, N) ; H
            Set of Morphisms from 2-dimensional differentiable manifold M to
             3-dimensional differentiable manifold N in Category of sets
            sage: TestSuite(H).run()

        Test for an endomorphism set::

            sage: E = Hom(M, M) ; E
            Set of Morphisms from 2-dimensional differentiable manifold M to
             2-dimensional differentiable manifold M in Category of sets
            sage: TestSuite(E).run()

        """
        from sage.manifolds.differentiable.manifold import DiffManifold
        if not isinstance(domain, DiffManifold):
            raise TypeError("domain = {} is not an ".format(domain) +
                            "instance of DiffManifold")
        if not isinstance(codomain, DiffManifold):
            raise TypeError("codomain = {} is not an ".format(codomain) +
                            "instance of DiffManifold")
        TopManifoldHomset.__init__(self, domain, codomain, name=name,
                                   latex_name=latex_name)

    #### Parent methods ####

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        EXAMPLE::

            sage: M = DiffManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = DiffManifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M,N)
            sage: H._coerce_map_from_(ZZ)
            False
            sage: H._coerce_map_from_(M)
            False
            sage: H._coerce_map_from_(N)
            False

        """
        #!# for the time being:
        return False

    #### End of parent methods ####

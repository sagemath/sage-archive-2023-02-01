r"""
Sets of morphisms between topologial manifolds

The class :class:`TopManifoldHomset` implements sets of morphisms between
two topological manifolds over the same topological field `K`, a morphism
being a *continuous map* for the category of topological manifolds.

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed., Springer (New
  York) (2011)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)
- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
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

from sage.categories.homset import Homset
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.manifolds.continuous_map import ContinuousMap

class TopManifoldHomset(Homset, UniqueRepresentation):
    r"""
    Set of continuous maps between two topological manifolds.

    Given two topological manifolds `M` and `N` over a topological field `K`,
    the class :class:`TopManifoldHomset` implements the set
    `\mathrm{Hom}(U,V)` of morphisms (i.e. continuous maps) `U\rightarrow V`,
    where `U` is an open subset of `M` and `V` an open subset of `N`.
    Note that, as open subsets of topological manifolds, `U` and `V` are
    topological manifolds by themselves.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.manifolds.continuous_map.ContinuousMap`.

    INPUT:

    - ``domain`` -- open subset `U\subset M` (domain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.manifold.TopManifold`
    - ``codomain`` -- open subset `V\subset N` (codomain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.manifold.TopManifold`
    - ``name`` -- (default: ``None``) string; name given to the homset; if
      none is provided, Hom(U,V) will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      homset; if none is provided, `\mathrm{Hom}(U,V)` will be used

    EXAMPLES:

    Set of continuous maps between a 2-dimensional manifold and a
    3-dimensional one::

        sage: M = TopManifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: N = TopManifold(3, 'N')
        sage: Y.<u,v,w> = N.chart()
        sage: H = Hom(M, N) ; H
        Set of Morphisms from 2-dimensional topological manifold M to
         3-dimensional topological manifold N in Category of sets
        sage: type(H)
        <class 'sage.manifolds.manifold_homset.TopManifoldHomset_with_category'>
        sage: H.category()
        Category of homsets of sets
        sage: latex(H)
        \mathrm{Hom}\left(M,N\right)
        sage: H.domain()
        2-dimensional topological manifold M
        sage: H.codomain()
        3-dimensional topological manifold N

    An element of ``H`` is a continuous map from ``M`` to ``N``::

        sage: H.Element
        <class 'sage.manifolds.continuous_map.ContinuousMap'>
        sage: f = H.an_element() ; f
        Continuous map from the 2-dimensional topological manifold M to the
         3-dimensional topological manifold N
        sage: f.display()
        M --> N
           (x, y) |--> (u, v, w) = (0, 0, 0)

    The test suite is passed::

        sage: TestSuite(H).run()

    When the codomain coincides with the domain, the homset is a set of
    *endomorphisms* in the category of topological manifolds::

        sage: E = Hom(M, M) ; E
        Set of Morphisms from 2-dimensional topological manifold M to
         2-dimensional topological manifold M in Category of sets
        sage: E.category()
        Category of endsets of sets
        sage: E.is_endomorphism_set()
        True
        sage: E is End(M)
        True

    In this case, the homset is a monoid for the law of morphism composition::

        sage: E in Monoids()
        True

    This was of course not the case of ``H = Hom(M, N)``::

        sage: H in Monoids()
        False

    The identity element of the monoid is of course the identity map of ``M``::

        sage: E.one()
        Identity map Id_M of the 2-dimensional topological manifold M
        sage: E.one() is M.identity_map()
        True
        sage: E.one().display()
        Id_M: M --> M
           (x, y) |--> (x, y)

    The test suite is passed by ``E``::

        sage: TestSuite(E).run()

    This test suite includes more tests than in the case of ``H``, since ``E``
    has some extra structure (monoid).

    """

    Element = ContinuousMap

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = TopManifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M, N) ; H
            Set of Morphisms from 2-dimensional topological manifold M to
             3-dimensional topological manifold N in Category of sets
            sage: TestSuite(H).run()

        Test for an endomorphism set::

            sage: E = Hom(M, M) ; E
            Set of Morphisms from 2-dimensional topological manifold M to
             2-dimensional topological manifold M in Category of sets
            sage: TestSuite(E).run()

        """
        from sage.manifolds.manifold import TopManifold
        if not isinstance(domain, TopManifold):
            raise TypeError("domain = {} is not an ".format(domain) +
                            "instance of TopManifold")
        if not isinstance(codomain, TopManifold):
            raise TypeError("codomain = {} is not an ".format(codomain) +
                            "instance of TopManifold")
        Homset.__init__(self, domain, codomain)
        if name is None:
            self._name = "Hom(" + domain._name + "," + codomain._name + ")"
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = \
                    r"\mathrm{Hom}\left(" + domain._latex_name + "," + \
                    codomain._latex_name + r"\right)"
        else:
            self._latex_name = latex_name
        self._one = None # to be set by self.one() if self is a set of
                         # endomorphisms (codomain = domain)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        EXAMPLE::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = TopManifold(3, 'N')
            sage: H = Hom(M, N)
            sage: H._latex_()
            '\\mathrm{Hom}\\left(M,N\\right)'

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name


    #### Parent methods ####

    def _element_constructor_(self, coord_functions, name=None, latex_name=None,
                              is_isomorphism=False, is_identity=False):
        r"""
        Construct an element of ``self``, i.e. a continuous map
        U --> V, where U is the domain of ``self`` and V its codomain.

        INPUT:

        - ``coord_functions`` -- a dictionary of the coordinate expressions
          (as lists (or tuples) of the coordinates of the image expressed in
          terms of the coordinates of the considered point) with the pairs of
          charts (chart1, chart2) as keys (chart1 being a chart on `U` and
          chart2 a chart on `V`). If the dimension of the arrival manifold
          is 1, a single coordinate expression can be passed instead of a tuple
          with a single element
        - ``name`` -- (default: ``None``) name given to the continuous map
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          continuous map; if none is provided, the LaTeX symbol is set to
          ``name``
        - ``is_isomorphism`` -- (default: ``False``) determines whether the
          constructed object is a isomorphism (i.e. a homeomorphism); if set to
          ``True``, then the manifolds `M` and `N` must have the same dimension.
        - ``is_identity`` -- (default: ``False``) determines whether the
          constructed object is the identity map; if set to ``True``,
          then `V` must be `U` and the entry ``coord_functions`` is not used.

        .. NOTE::

            If the information passed by means of the argument
            ``coord_functions`` is not sufficient to fully specify the
            continuous map, further coordinate expressions, in other charts,
            can be subsequently added by means of the method :meth:`add_expr`

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.continuous_map.ContinuousMap`

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = TopManifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M, N)
            sage: f = H._element_constructor_({(X, Y): [x+y, x-y, x*y]}, name='f') ; f
            Continuous map f from the 2-dimensional topological manifold M to
             the 3-dimensional topological manifold N
            sage: f.display()
            f: M --> N
               (x, y) |--> (u, v, w) = (x + y, x - y, x*y)
            sage: id = Hom(M, M)._element_constructor_({}, is_identity=True) ; id
            Identity map Id_M of the 2-dimensional topological manifold M
            sage: id.display()
            Id_M: M --> M
               (x, y) |--> (x, y)

        """
        # Standard construction
        return self.element_class(self, coord_functions=coord_functions,
                                  name=name, latex_name=latex_name,
                                  is_isomorphism=is_isomorphism,
                                  is_identity=is_identity)

    def _an_element_(self):
        r"""
        Construct some element.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.continuous_map.ContinuousMap`

        EXAMPLE::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = TopManifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M,N)
            sage: f = H._an_element_() ; f
            Continuous map from the 2-dimensional topological manifold M to the
             3-dimensional topological manifold N
            sage: f.display()
            M --> N
               (x, y) |--> (u, v, w) = (0, 0, 0)
            sage: p = M((-2,3)) ; p
            Point on the 2-dimensional topological manifold M
            sage: f(p)
            Point on the 3-dimensional topological manifold N
            sage: f(p).coord(Y)
            (0, 0, 0)
            sage: TestSuite(f).run()

        """
        dom = self.domain()
        codom = self.codomain()
        # A constant map is constructed:
        chart2 = codom.default_chart()
        target_point = chart2.domain().an_element()
        target_coord = target_point.coord(chart2)
        coord_functions = {}
        for chart in dom.atlas():
            coord_functions[(chart, chart2)] = target_coord
        return self.element_class(self, coord_functions=coord_functions)

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent.

        EXAMPLE::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = TopManifold(3, 'N')
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

    #!# check
    def __call__(self, *args, **kwds):
        r"""
        To bypass Homset.__call__, enforcing Parent.__call__ instead.

        EXAMPLE::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = TopManifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M,N)
            sage: f = H.__call__({(X, Y): [x+y, x-y, x*y]}, name='f') ; f
            Continuous map f from the 2-dimensional topological manifold M to
             the 3-dimensional topological manifold N
            sage: f.display()
            f: M --> N
               (x, y) |--> (u, v, w) = (x + y, x - y, x*y)

        """
        return Parent.__call__(self, *args, **kwds)

    #### End of parent methods ####

    #### Monoid methods (case of an endomorphism set) ####

    def one(self):
        r"""
        Return the identity element of ``self`` considered as a monoid (case of
        an endomorphism set).

        This applies only when the codomain of ``self`` is equal to its domain,
        i.e. when ``self`` is of the type `\mathrm{Hom}(U,U)` where `U` is
        an open subset of some manifold. `\mathrm{Hom}(U,U)` equiped with the
        law of morphisms composition is then a monoid, whose identity element
        is nothing but the identity map of `U`.

        OUTPUT:

        - the identity map of `U`, as an instance of
          :class:`~sage.manifolds.continuous_map.ContinuousMap`

        EXAMPLE:

        The identity map of a 2-dimensional manifold::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: H = Hom(M, M) ; H
            Set of Morphisms from 2-dimensional topological manifold M to
             2-dimensional topological manifold M in Category of sets
            sage: H in Monoids()
            True
            sage: H.one()
            Identity map Id_M of the 2-dimensional topological manifold M
            sage: H.one().parent() is H
            True
            sage: H.one().display()
            Id_M: M --> M
               (x, y) |--> (x, y)

        The identity map is cached::

            sage: H.one() is H.one()
            True

        If ``self`` is not a set of endomorphisms, the identity element is
        meaningless::

            sage: N = TopManifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: Hom(M, N).one()
            Traceback (most recent call last):
            ...
            TypeError: the Set of Morphisms from 2-dimensional topological
             manifold M to 3-dimensional topological manifold N in Category of
             sets is not a monoid

        """
        if self._one is None:
            if self.codomain() != self.domain():
                raise TypeError("the {} is not a monoid".format(self))
            self._one = self.element_class(self, is_identity=True)
        return self._one

    #### End of monoid methods ####

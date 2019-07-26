r"""
Topological Vector Bundle

Let `K` be a topological field. A *vector bundle* of rank `n` over the field
`K` and over a topological manifold `B` (base space) is a topological manifold
`E` (total space) together with a continuous and surjective projection
`\pi: E \to B` such that for every point `x \in B`
- the set `E_x=\pi^{-1}(x)` has the vector space structure of `K^n`,
- there is a neighborhood `U \subset B` of `x` and a homeomorphism
  `\varphi: \pi^{-1}(x) \to U \times K^n` such that
  `v \mapsto \varphi^{-1}(y,v)` is a linear isomorphism for any `y \in U`.

AUTHORS:

- Michael Jung (2019) : initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.category_object import CategoryObject
from sage.categories.vector_bundles import VectorBundles
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import CC
from sage.rings.real_mpfr import RR, RealField_class
from sage.rings.complex_field import ComplexField_class
from sage.rings.integer import Integer
from sage.manifolds.vector_bundle_fiber import VectorBundleFiber
from sage.manifolds.section_module import SectionModule

class TopologicalVectorBundle(CategoryObject, UniqueRepresentation):
    r"""
    An instance of this class is a topological vector bundle `E \to B` over a
    topological field `K`.

    INPUT:

    - ``rank`` -- positive integer; rank of the vector bundle
    - ``name```-- string representation given to the total space
    - ``base_space`` -- the base space (topological manifold) over which the
      vector bundle shall be defined
    - ``field`` -- field `K` which gives the fibers the structure of a
     vector space over `K`; allowed values are

      - ``'real'`` or an object of type ``RealField`` (e.g., ``RR``) for
        a vector bundle over `\RR`
      - ``'complex'`` or an object of type ``ComplexField`` (e.g., ``CC``)
        for a vector bundle over `\CC`
      - an object in the category of topological fields (see
        :class:`~sage.categories.fields.Fields` and
        :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
        for other types of topological fields

    - ``latex_name`` -- (default: ``None``) latex representation given to the
      total space
    - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``VectorBundles(base_space, c_field)`` is assumed (see the
      category :class:`~sage.categories.vector_bundles.VectorBundles`)

    EXAMPLES:

    A rational line bundle over some 4-dimensional topological manifold::

        sage: M = Manifold(4, 'M', structure='top')
        sage: E = M.vector_bundle(1, 'E', field=QQ); E
        Topological neither_real_nor_complex vector bundle E -> M of rank 1 over
         the base space 4-dimensional topological manifold M
        sage: E.base_space()
        4-dimensional topological manifold M
        sage: E.base_ring()
        Rational Field
        sage: E.rank()
        1

    """
    def __init__(self, rank, name, base_space, field='real',
                 latex_name=None, category=None, unique_tag=None):
        r"""
        Construct a topological vector bundle.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: from sage.manifolds.vector_bundle import TopologicalVectorBundle
            sage: TopologicalVectorBundle(2, 'E', M)
            Topological real vector bundle E -> M of rank 2 over the base space
             2-dimensional topological manifold M

        """
        if base_space is None:
            raise ValueError("a base space must be provided")
        # Handle the field:
        if field == 'real':
            self._field = RR
            self._field_type = field
        elif field == 'complex':
            self._field = CC
            self._field_type = field
        else:
            self._field = field
            if isinstance(field, RealField_class):
                self._field_type = 'real'
            elif isinstance(field, ComplexField_class):
                self._field_type = 'complex'
            else:
                self._field_type = 'neither_real_nor_complex'
        # Get the category:
        if category is None:
            category = VectorBundles(base_space, self._field)
        CategoryObject.__init__(self, base=self._field,
                                category=category)
        # Check rank:
        if not isinstance(rank, (int, Integer)):
            raise TypeError("the rank must be an integer")
        if rank < 1:
            raise ValueError("the rank must be strictly positive")
        # Define left attributes:
        self._rank = rank
        self._base_space = base_space
        # Set names:
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        ###
        # Now, the trivializations come into play:
        self._atlas = []  # list of trivializations defined on self
        self._transitions = {} # dictionary of transition maps (key: pair of
                    # of trivializations)
        self._frames = []  # list of local frames for self
        self._def_frame = None  # default frame
        self._coframes = [] # list of local coframes for self
        self._trivial_parts = set() # subsets of base space on which self is
                    # trivial

    def base_space(self):
        r"""
        Return the base space of the vector bundle

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.base_space()
            2-dimensional topological manifold M

        """
        return self._base_space

    def base_field_type(self):
        r"""
        Return the type of topological field on which the fibers are defined.

        OUTPUT:

        - a string describing the field, with three possible values:

          - ``'real'`` for the real field `\RR`
          - ``'complex'`` for the complex field `\CC`
          - ``'neither_real_nor_complex'`` for a field different
            from `\RR`, `\CC` and `\QQ`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E', field=QQ)
            sage: E.base_field_type()
            'neither_real_nor_complex'

        """
        return self._field_type

    def base_field(self):
        r"""
        Return the field on which the fibers are defined.

        OUTPUT:

        - a topological field

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: E = M.vector_bundle(2, 'E', field=QQ)
            sage: E.base_field()
            Rational Field

        """
        return self._field

    def rank(self):
        r"""
        Return the rank of the vector bundle.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(3, 'E')
            sage: E.rank()
            3

        """
        return self._rank

    def _repr_object_name(self):
        r"""
        String name of the object without structure.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: E._repr_object_name()
            'real vector bundle E -> M of rank 1 over the base space
             2-dimensional topological manifold M'

        """
        desc = self.base_field_type() + " "
        desc += "vector bundle "
        desc += self._name + " -> " + self.base_space()._name + " "
        desc += "of rank {} ".format(self._rank)
        desc += "over the base space {}".format(self.base_space())
        return desc

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: E._repr_()
            'Topological real vector bundle E -> M of rank 1 over the base space
             2-dimensional topological manifold M'

        """
        desc = "Topological "
        desc += self._repr_object_name()
        return desc

    def _latex_(self):
        r"""
        Return the LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: E._latex_()
            'E\\to M'

        """
        latex = self._latex_name
        latex += r'\to '
        latex += self.base_space()._latex_name
        return latex

    def trivialization(self, domain=None, name=None, latex_name=None):
        r"""
        Return a trivialization of ``self`` over the domain ``domain``.

        OUTPUT:

        - a :class:`~sage.manifolds.trivialization.Trivialization` representing
          a trivialization of `E`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization(U, name='phi'); phi
            Trivialization (phi:E|_U -> U)

        """
        from .trivialization import Trivialization
        return Trivialization(self, domain=domain, name=name,
                              latex_name=latex_name)

    def transitions(self):
        r"""
        Return the transition maps defined over subsets of the manifold.

        OUTPUT:

        - dictionary of transition maps, with pairs of trivializations as keys

        EXAMPLES::



        """
        return self._transitions

    def transition(self, triv1, triv2):
        r"""
        Return the transition map between two trivializations defined over the
        manifold.

        The transition map must have been defined previously, for instance by
        the method
        :meth:`~sage.manifolds.trivialization.Trivialization.transition_map`.

        INPUT:

        - ``triv1`` -- trivialization 1
        - ``triv2`` -- trivialization 2

        OUTPUT:

        - instance of :class:`~sage.manifolds.trivialization.TransitionMap`
          representing the transition map from trivialization 1 to
          trivialization 2

        EXAMPLES::



        """
        if (triv1, triv2) not in self._transitions:
            raise TypeError("the transition map from " +
                            "{} to {}".format(triv1, triv2) + " has not " +
                            "been defined on the {}".format(self))
        return self._transitions[(triv1, triv2)]

    def atlas(self):
        r"""
        Return the list of trivializations that have been defined for ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: U = M.open_subset('U')
            sage: triv_U = E.trivialization(U); triv_U
            Trivialization (E|_U -> U)
            sage: V = M.open_subset('V')
            sage: triv_V = E.trivialization(V); triv_V
            Trivialization (E|_V -> V)
            sage: triv_M E.trivialization(M); triv_M
            sage: E.atlas()
            [Trivialization (E|_U -> U),
             Trivialization (E|_V -> V),
             Trivialization (E|_M -> M)]

        """
        return list(self._atlas) # Make a (shallow) copy

    def is_manifestly_trivial(self):
        r"""
        Return ``True`` if ``self`` is manifestly a trivial bundle, i.e. there
        exists a trivialization defined on the whole base space.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: U = M.open_subset('U')
            sage: triv_U = E.trivialization(U); triv_U
            Trivialization (E|_U -> U)
            sage: V = M.open_subset('V')
            sage: triv_V = E.trivialization(V); triv_V
            Trivialization (E|_V -> V)
            sage: E.is_manifestly_trivial()
            False
            sage: E.trivialization(M)
            Trivialization (E|_M -> M)
            sage: E.is_manifestly_trivial()
            True

        """
        return self.base_space() in self._trivial_parts

    def section_module(self, domain=None, force_free=False):
        r"""
        Return the section module of continuous section on ``self``.

        See :class:`~sage.manifolds.section_module.SectionModule` for a complete
        documentation.

        INPUT:

        - ``force_free`` -- (default: ``False``) if set to ``True``, force
          the construction of a *free* module (this implies that `E` is trivial)

        OUTPUT:

        - a
          :class:`~sage.manifolds.section_module.SectionModule`
          (or if `E` is trivial, a
          :class:`~sage.manifolds.section_module.SectionFreeModule`)
          representing the module `\Gamma(E)` of continuous sections on
          `M` taking values on `E`

        """
        from sage.manifolds.section_module import \
                                            SectionModule, SectionFreeModule
        if domain is None:
            domain = self._base_space
        for triv in self._atlas:
            if triv.domain() is domain:
                return SectionFreeModule(domain)
        for frame in self._frames:
            if frame.domain() is domain:
                return True
        return False

    def fiber(self, point):
        r"""
        Return the vector bundle fiber at ``point``.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: p = M((0,2,1), name='p'); p
            Point p on the 3-dimensional topological manifold M
            sage: E = M.vector_bundle(2, 'E'); E
            Topological real vector bundle E -> M of rank 2 over the base space
             3-dimensional topological manifold NoneM
            sage: E.fiber(p)
            Fiber of E at Point p on the 3-dimensional topological manifold M

        """
        return VectorBundleFiber(self, point)

    def section(self, name=None, latex_name=None):
        r"""
        Return a continuous section of ``self``.

        INPUT:

        - ``domain`` -- (default: ````) domain on which the section shall be
          defined; if ``None``, the base space is assumed

        OUTPUT:

        - an instance of :class:`~sage.manifolds.section.Section` representing
          a continuous section of `M` with values on `E`

        """
        # TODO: Implement
        pass

    def whitney_sum(self, vector_bundle):
        r"""
        Return the Whitney sum ``self`` with ``vector_bundle``.

        INPUT:

        - ``vector_bundle`` --

        OUTPUT:

        -

        """
        # TODO: Implement
        pass

    def tensor_product(self, vector_bundle):
        r"""
        Return the tensor product of ``self```with ``vector_bundle``.

        INPUT:

        - ``vector_bundle`` --

        OUTPUT:

        -

        """
        # TODO: Implement
        pass

    def total_space(self):
        r"""
        Return the total space of ``self``.

        OUTPUT:

        - the total space of ``self`` as an instance of :class:`~sage.manifolds.manifold.TopologicalManifold`

        """
        from sage.manifolds.manifold import Manifold
        base_space = self._base_space
        dim = base_space._dim * self._rank
        sindex = base_space.start_index()
        total_space = Manifold(dim, self._name, latex_name=self._latex_name,
                               field=self._field, structure='topological',
                               start_index=sindex)
        # TODO: if self._atlas not empty, introduce charts
        return total_space

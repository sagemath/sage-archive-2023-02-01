r"""
Topological Vector Bundle

Let `K` be a topological field. A *vector bundle* of rank `n` over the field
`K` and over a topological manifold `B` (base space) is a topological manifold
`E` (total space) together with a continuous and surjective projection
`\pi: E \to B` such that for every point `p \in B`
- the set `E_x=\pi^{-1}(x)` has the vector space structure of `K^n`,
- there is a neighborhood `U \subset B` of `p` and a homeomorphism
  `\varphi: \pi^{-1}(p) \to U \times K^n` such that
  `v \mapsto \varphi^{-1}(q,v)` is a linear isomorphism for any `q \in U`.

AUTHORS:

- Michael Jung (2019) : initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
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
      vector bundle is defined
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

    A real line bundle over some 4-dimensional topological manifold::

        sage: M = Manifold(4, 'M', structure='top')
        sage: E = M.vector_bundle(1, 'E'); E
        Topological real vector bundle E -> M of rank 1 over the base space
         4-dimensional topological manifold M
        sage: E.base_space()
        4-dimensional topological manifold M
        sage: E.base_ring()
        Real Field with 53 bits of precision
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
        bs_field = base_space.base_field()
        if not bs_field.is_subring(self._field):
            raise ValueError("for concrete implementation, manifold's base "
                             "field must be a subfield of the vector bundle's "
                             "base field")
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
        # Define remaining attributes:
        self._rank = rank
        self._diff_degree = 0
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
        self._frame_changes = {}  # dictionary of changes of frames
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
            sage: E = M.vector_bundle(2, 'E', field=CC)
            sage: E.base_field_type()
            'complex'

        """
        return self._field_type

    def base_field(self):
        r"""
        Return the field on which the fibers are defined.

        OUTPUT:

        - a topological field

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: E = M.vector_bundle(2, 'E', field=CC)
            sage: E.base_field()
            Complex Field with 53 bits of precision

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

    def _add_local_frame(self, frame):
        r"""
        Helper method to add local frames to the vector bundle.

        INPUT:

        - ``frame`` -- the local frame that shall be added

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: E._add_local_frame(e)
            sage: E._frames
            [Local frame (E|_M, (e_0,e_1)), Local frame (E|_M, (e_0,e_1))]

        """
        self._trivial_parts.add(frame.domain())
        self._frames.append(frame)

    def trivialization(self, name, domain=None, latex_name=None):
        r"""
        Return a trivialization of ``self`` over the domain ``domain``.

        INPUT:

        - ``domain`` -- (default: ``None``) domain on which the trivialization
          is defined; if ``None`` the base space is assumed
        - ``name`` -- (default: ``None``) name given to the trivialization
        - ``latex_name`` -- (default: ``None``) latex name given to the
          trivialization

        OUTPUT:

        - a :class:`~sage.manifolds.trivialization.Trivialization` representing
          a trivialization of `E`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi', domain=U); phi
            Trivialization (phi, E|_U)

        """
        if domain is None:
            domain = self._base_space
        from .trivialization import Trivialization
        return Trivialization(self, name, domain=domain, latex_name=latex_name)

    def transitions(self):
        r"""
        Return the transition maps defined over subsets of the base space.

        OUTPUT:

        - dictionary of transition maps, with pairs of trivializations as keys

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: X_UV = X.restrict(U.intersection(V))
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_U', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, 1)
            sage: E.transitions()
            {(Trivialization (phi_U, E|_U),
             Trivialization (phi_U, E|_V)): Transition map from Trivialization
             (phi_U, E|_U) to Trivialization (phi_U, E|_V),
             (Trivialization (phi_U, E|_V),
             Trivialization (phi_U, E|_U)): Transition map from Trivialization
             (phi_U, E|_V) to Trivialization (phi_U, E|_U)}

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

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: X_UV = X.restrict(U.intersection(V))
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_U_to_phi_V = phi_U.transition_map(phi_V, 1)
            sage: E.transition(phi_V, phi_U)
            Transition map from Trivialization (phi_V, E|_V) to Trivialization
             (phi_U, E|_U)

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

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: phi_V = E.trivialization('phi_V', domain=V)
            sage: phi_M = E.trivialization('phi_M')
            sage: E.atlas()
            [Trivialization (phi_U, E|_U),
             Trivialization (phi_V, E|_V),
             Trivialization (phi_M, E|_M)]

        """
        return list(self._atlas) # Make a (shallow) copy

    def is_manifestly_trivial(self):
        r"""
        Return ``True`` if ``self`` is manifestly a trivial bundle, i.e. there
        exists a frame or a trivialization defined on the whole base space.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: phi_U = E.trivialization('phi_U', domain=U); phi_U
            Trivialization (phi_U, E|_U)
            sage: phi_V = E.trivialization('phi_V', domain=V); phi_V
            Trivialization (phi_V, E|_V)
            sage: E.is_manifestly_trivial()
            False
            sage: E.trivialization('phi_M', M)
            Trivialization (phi_M, E|_M)
            sage: E.is_manifestly_trivial()
            True

        """
        return self.base_space() in self._trivial_parts

    def section_module(self, domain=None, force_free=False):
        r"""
        Return the section module of continuous sections on ``self``.

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

        EXAMPLES:

        Module of sections on the Moebius bundle::

            sage: M = Manifold(1, 'S^1', structure='top', start_index=1)
            sage: U = M.open_subset('U')  # the complement of one point
            sage: c_t.<t> =  U.chart('t:(0,2*pi)') # the standard angle coordinate
            sage: V = M.open_subset('V') # the complement of the point t=pi
            sage: M.declare_union(U,V)   # S^1 is the union of U and V
            sage: c_u.<u> = V.chart('u:(0,2*pi)') # the angle t-pi
            sage: t_to_u = c_t.transition_map(c_u, (t-pi,),
            ....:                             intersection_name='W',
            ....:                             restrictions1 = t!=pi,
            ....:                             restrictions2 = u!=pi)
            sage: u_to_t = t_to_u.inverse()
            sage: W = U.intersection(V)
            sage: E = M.vector_bundle(1, 'E')
            sage: phi_U = E.trivialization('phi_U', latex_name=r'\varphi_U',
            ....:                          domain=U)
            sage: phi_V = E.trivialization('phi_V', latex_name=r'\varphi_V',
            ....:                          domain=V)
            sage: transf = phi_U.transition_map(phi_V, [[-1]])
            sage: C0 = E.section_module(); C0
            Module C^0(S^1;E) of sections on the 1-dimensional topological
             manifold S^1 with values in the real vector bundle E of rank 1

        `C^0(S^1;E)` is a module over the algebra `C^0(M)`::

            sage: C0.category()
            Category of modules over Algebra of scalar fields on the 1-dimensional
             topological manifold S^1
            sage: C0.base_ring() is M.scalar_field_algebra()
            True

        However, `C^0(S^1;E)` is not a free module::

            sage: isinstance(C0, FiniteRankFreeModule)
            False

        since the Moebius bundle is not trivial::

            sage: E.is_manifestly_trivial()
            False

        The section module over `U`, on the other hand, is a free module since
        `E|_U` admits a trivialization and therefore has a local frame::

            sage: C0_U = E.section_module(domain=U)
            sage: isinstance(C0_U, FiniteRankFreeModule)
            True

        The elements of `C^0(U)` are sections on `U`::

            sage: C0_U.an_element()
            Section on the Open subset U of the 1-dimensional topological
             manifold S^1 with values in the real vector bundle E of rank 1
            sage: C0_U.an_element().display(phi_U.frame())
            2 (phi_U^*e_1)


        """
        if domain is None:
            domain = self._base_space
        from sage.manifolds.section_module import \
                                            SectionModule, SectionFreeModule
        if force_free or domain in self._trivial_parts:
            return SectionFreeModule(self, domain)
        return SectionModule(self, domain)

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
             3-dimensional topological manifold M
            sage: E.fiber(p)
            Fiber of E at Point p on the 3-dimensional topological manifold M

        """
        return VectorBundleFiber(self, point)

    def local_frame(self, symbol, latex_symbol=None, indices=None,
                    latex_indices=None, symbol_dual=None,
                    latex_symbol_dual=None, domain=None):
        r"""
        Define a local frame on ``self``.

        A *local frame* is a section on a subset `U \subset M` in `E` that
        provides, at each point `p` of the base space, a vector basis of the
        fiber `E_p` at `p`.

        .. SEEALSO::

            :class:`~sage.manifolds.local_frame.LocalFrame` for complete
            documentation.

        INPUT:

        - ``symbol`` -- (default: ``None``) either a string, to be used as a
          common base for the symbols of the sections constituting the
          local frame, or a list/tuple of strings, representing the individual
          symbols of the sections; can be ``None`` only if ``from_frame``
          is not ``None`` (see below)
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the sections
          constituting the local frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the sections;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the sections;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe
        - ``domain`` -- (default: ``None``) domain on which the local frame
          is defined; if ``None``, the whole base space is assumed

        OUTPUT:

        - a :class:`~sage.manifolds.local_frame.LocalFrame` representing the
          defined local frame

        EXAMPLES:

        Setting a local frame on a real rank-2 vector bundle::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U); e
            Local frame (E|_U, (e_0,e_1))

        For a global frame, the argument ``domain`` is omitted::

            sage: f = E.local_frame('f'); f
            Local frame (E|_M, (f_0,f_1))

        .. SEEALSO::

            For more options, in particular for the choice of symbols and
            indices, see :class:`~sage.manifolds.local_frame.LocalFrame`.

        """
        from sage.manifolds.local_frame import LocalFrame
        sec_module = self.section_module(domain=domain, force_free=True)
        return LocalFrame(sec_module, symbol=symbol, latex_symbol=latex_symbol,
                          indices=indices, latex_indices=latex_indices,
                          symbol_dual=symbol_dual,
                          latex_symbol_dual=latex_symbol_dual)

    def section(self, *comp, **kwargs):
        r"""
        Return a continuous section of ``self``.

        INPUT:

        - ``domain`` -- (default: ``None``) domain on which the section shall be
          defined; if ``None``, the base space is assumed
        - ``name`` -- (default: ``None``) name of the local section
        - ``latex_name`` -- (default``None``) latex representation of the local
          section

        OUTPUT:

        - an instance of :class:`~sage.manifolds.section.Section` representing
          a continuous section of `M` with values on `E`

        EXAMPLES::

        TODO

        """
        domain = kwargs.pop('domain', self._base_space)
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        smodule = self.section_module(domain=domain)  # the parent
        resu = smodule.element_class(smodule, name=name, latex_name=latex_name)
        if comp:
            # Some components are to be initialized
            resu._init_components(*comp, **kwargs)
        return resu

    # def whitney_sum(self, vector_bundle):
    #     r"""
    #     Return the Whitney sum ``self`` with ``vector_bundle``.
    #
    #     INPUT:
    #
    #     - ``vector_bundle`` --
    #
    #     OUTPUT:
    #
    #     -
    #
    #     """
    #     # TODO: Implement
    #     pass
    #
    # def tensor_product(self, vector_bundle):
    #     r"""
    #     Return the tensor product of ``self```with ``vector_bundle``.
    #
    #     INPUT:
    #
    #     - ``vector_bundle`` --
    #
    #     OUTPUT:
    #
    #     -
    #
    #     """
    #     # TODO: Implement
    #     pass

    def total_space(self):
        r"""
        Return the total space of ``self``.

        OUTPUT:

        - the total space of ``self`` as an instance of
          :class:`~sage.manifolds.manifold.TopologicalManifold`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: E.total_space()
            6-dimensional topological manifold E

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

    def set_change_of_frame(self, frame1, frame2, change_of_frame,
                            compute_inverse=True):
        r"""
        Relate two vector frames by an automorphism.

        This updates the internal dictionary ``self._frame_changes``.

        INPUT:

        - ``frame1`` -- frame 1, denoted `(e_i)` below
        - ``frame2`` -- frame 2, denoted `(f_i)` below
        - ``change_of_frame`` -- instance of class
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          describing the automorphism `P` that relates the basis `(e_i)` to
          the basis `(f_i)` according to `f_i = P(e_i)`
        - ``compute_inverse`` (default: True) -- if set to True, the inverse
          automorphism is computed and the change from basis `(f_i)` to `(e_i)`
          is set to it in the internal dictionary ``self._frame_changes``

        EXAMPLES:

            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: f = E.local_frame('f')
            sage: a = E.section_module().automorphism()
            sage: a[e,:] = [[1,2],[0,3]]
            sage: E.set_change_of_frame(e, f, a)
            sage: f[0].display(e)
            f_0 = e_0
            sage: f[1].display(e)
            f_1 = 2 e_0 + 3 e_1
            sage: e[0].display(f)
            e_0 = f_0
            sage: e[1].display(f)
            e_1 = -2/3 f_0 + 1/3 f_1
            sage: M.change_of_frame(e,f)[e,:]
            [1 2]
            [0 3]

        """
        from sage.tensor.modules.free_module_automorphism import \
            FreeModuleAutomorphism
        sec_module = frame1._fmodule
        if frame2._fmodule != sec_module:
            raise ValueError("the two frames are not defined on the same " +
                             "section module")
        if isinstance(change_of_frame, FreeModuleAutomorphism):
            auto = change_of_frame
        else: # Otherwise try to coerce the input
            auto_group = sec_module.general_linear_group()
            auto = auto_group(change_of_frame, basis=frame1)
        sec_module.set_change_of_basis(frame1, frame2, auto,
                                       compute_inverse=compute_inverse)
        self._frame_changes[(frame1, frame2)] = auto
        if compute_inverse:
            self._frame_changes[(frame2, frame1)] = ~auto

    def change_of_frame(self, frame1, frame2):
        r"""
        Return a change of local frames defined on ``self``.

        INPUT:

        - ``frame1`` -- local frame 1
        - ``frame2`` -- local frame 2

        OUTPUT:

        - a
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`
          representing, at each point, the vector space automorphism `P`
          that relates frame 1, `(e_i)` say, to frame 2, `(f_i)` say,
          according to `f_i = P(e_i)`

        EXAMPLES:

            sage: M = Manifold(3, 'M', structure='top')
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: a = E.section_module().automorphism() # Now, the section module is free
            sage: a[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
            sage: f = e.new_frame(a, 'f')
            sage: E.change_of_frame(e, f)
            Automorphism of the Free module C^0(M;E) of sections on the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2
            sage: a == E.change_of_frame(e, f)
            True
            sage: a.inverse() == E.change_of_frame(f, e)
            True

        """
        if (frame1, frame2) not in self._frame_changes:
            raise ValueError("the change of frame from {} to {}".format(frame1, frame2) +
                             " has not been defined on the {}".format(self))
        return self._frame_changes[(frame1, frame2)]

    def changes_of_frame(self):
        r"""
        Return all the changes of local frames defined on ``self``.

        OUTPUT:

        - dictionary of vector bundle automorphisms representing
          the changes of frames, the keys being the pair of frames

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: auto_group = E.section_module().general_linear_group()
            sage: e_to_f = auto_group([[0,1],[1,0]]); e_to_f
            Automorphism of the Free module C^0(M;E) of sections on the
             3-dimensional topological manifold M with values in the real vector
             bundle E of rank 2
            sage: f_in_e = auto_group([[0,1],[1,0]])
            sage: f = e.new_frame(f_in_e, 'f'); f
            Local frame (E|_M, (f_0,f_1))
            sage: E.changes_of_frame()
            {(Local frame (E|_M, (f_0,f_1)),
             Local frame (E|_M, (e_0,e_1))): Automorphism of the Free module
             C^0(M;E) of sections on the 3-dimensional topological manifold M
             with values in the real vector bundle E of rank 2,
             (Local frame (E|_M, (e_0,e_1)),
             Local frame (E|_M, (f_0,f_1))): Automorphism of the Free module
             C^0(M;E) of sections on the 3-dimensional topological manifold M
             with values in the real vector bundle E of rank 2}

        """
        return self._frame_changes

    def frames(self):
        r"""
        Return the list of local frames defined on ``self``.

        OUTPUT:

        - list of local frames defined on ``self``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: e = E.local_frame('e', domain=V)
            sage: E.frames()
            [Trivialization frame (E|_U, ((phi_U^*e_1),(phi_U^*e_2))),
             Local frame (E|_V, (e_0,e_1))]

        """
        return self._frames

    def coframes(self):
        r"""
        Return the list of coframes defined on ``self``.

        OUTPUT:

        - list of coframes defined on ``self``

        EXAMPLES:

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: e = E.local_frame('e', domain=V)
            sage: E.coframes()
            [Trivialization coframe (E|_U, ((phi_U^*e^1),(phi_U^*e^2))),
             Local coframe (E|_V, (e^0,e^1))]

        """
        return self._coframes

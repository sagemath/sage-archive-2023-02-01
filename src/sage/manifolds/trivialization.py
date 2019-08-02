r"""
Trivializations

AUTHORS:

- Michael Jung (2019) : initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.manifolds.local_frame import TrivializationFrame

class Trivialization(UniqueRepresentation, SageObject):
    r"""
    A local trivialization of a given vector bundle.

    Let `\pi:E \to M` be a topological vector bundle of rank `n` over the field
    `K` (see: :class:~sage.manifolds.vector_bundle.TopologicalVectorBundle`). A
    *local trivialization* over an open subset `U \subset M` is a homeomorphism
    `\varphi: \pi^{-1}(U) \to U \times K^n` such that `(p, v) \mapsto \pi^{-1}(p)`
    is a linear isomorphism.

    EXAMPLES:

    """
    def __init__(self, vector_bundle, domain, name=None, latex_name=None):
        r"""
        Construct a local trivialization of the vector bundle ``vector_bundle``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: E = M.vector_bundle(2, 'E')
            sage: triv = E.trivialization(M)
            sage: TestSuite(triv).run()

        """
        if domain is None:
            domain = vector_bundle.base_space()
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._base_space = vector_bundle.base_space()
        self._vbundle = vector_bundle
        self._bdl_rank = vector_bundle.rank()
        self._base_field = vector_bundle.base_field()
        self._sindex = self._base_space.start_index()
        self._domain = domain
        # Add this trivialization to the atlas of the vector bundle:
        vector_bundle._atlas.append(self)
        self._frame = TrivializationFrame(self)
        self._coframe = self._frame._coframe

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: triv = E.trivialization(M)
            sage: triv._repr_()
            'Trivialization (E|_M -> M)'

        """
        desc = "Trivialization ("
        if self._name is not None:
            desc += self._name + ":"
        desc += "{}|_{} -> {})".format(self._vbundle._name,
                                               self._base_space._name,
                                               self._base_space._name)
        return desc

    def _latex_(self):
        r"""
        Return the LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(1, 'E')
            sage: triv = E.trivialization(M)
            sage: triv._latex_()
            'E |_{M} \\to M \\times \\Bold{R}^1'

        """
        latex = str()
        if self._latex_name is not None:
            latex += self._latex_name + r':'
        latex += r'{} |_{{{}}} \to {} \times {}^{}'.format(self._vbundle._latex_name,
                            self._domain._latex_(), self._domain._latex_(),
                            self._base_field._latex_(), self._bdl_rank)
        return latex

    def base_space(self):
        r"""
        Return the manifold on which the trivialization is defined.

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization(U)
            sage: phi_U.base_space()
            2-dimensional topological manifold M

        """
        return self._base_space

    def transition_map(self, other, transf, compute_inverse=True):
        r"""
        Return the transition map between ``self`` and ``other``.

        INPUT:

        - ``other`` -- the trivialization where the transition map from ``self``
          goes to
        - ``transf`` -- transformation of the transition map

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization(U, name='phi_U')
            sage: phi_V = E.trivialization(V, name='phi_V')
            sage: phi_U.transition_map(phi_V, 1)
            Transition map from Trivialization (phi_U:E|_U -> U) to
             Trivialization (phi_V:E|_V -> V)


        """
        return TransitionMap(self, other, transf,
                             compute_inverse=compute_inverse)

    def vector_bundle(self):
        r"""
        Return the vector bundle on which the trivialization is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization(U)
            sage: phi_U.vector_bundle()
            Topological real vector bundle E -> M of rank 2 over the base space
             2-dimensional topological manifold M

        """
        return self._vbundle

    def domain(self):
        r"""
        Return the domain on which the trivialization is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization(U)
            sage: phi_U.domain()
            Open subset U of the 2-dimensional topological manifold M

        """
        return self._domain

    def frame(self):
        r"""

        """
        return self._frame

    def coframe(self):
        r"""

        """
        return self._frame._coframe

# *****************************************************************************

class TransitionMap(SageObject):
    r"""

    """
    def __init__(self, triv1, triv2, transf, compute_inverse=True):
        r"""
        Construct a transition map between two trivializations.

        TESTS::



        """
        bs1 = triv1.base_space()
        bs2 = triv2.base_space()
        if bs1 is not bs2:
            raise ValueError("base spaces must coincide")
        self._base_space = bs1

        vb1 = triv1.vector_bundle()
        vb2 = triv2.vector_bundle()
        if vb1 is not vb2:
            raise ValueError("vector bundles must coincide")
        self._vbundle = vb1
        self._bdl_rank = self._vbundle.rank()

        mfd_field = self._base_space.base_field()
        vb_field = self._vbundle.base_field()
        if not mfd_field.is_subring(vb_field):
            raise ValueError("for concrete implementation, manifold's base "
                             "field must be a subfield of the vector bundle's "
                             "base field")
        dom1 = triv1.domain()
        dom2 = triv2.domain()
        dom = dom1.intersection(dom2)
        self._domain = dom
        self._frame1 = triv1._frame.restrict(dom)
        self._frame2 = triv2._frame.restrict(dom)
        self._triv1 = triv1
        self._triv2 = triv2
        self._inverse = None
        self._vbundle._transitions[(triv1, triv2)] = self
        ###
        # Define the automorphism
        sec_module = self._vbundle.section_module(dom, force_free=True)
        auto_group = sec_module.general_linear_group()
        auto = auto_group(transf, basis=self._frame1)
        self._automorphism = auto
        # Add this change of basis to the basis changes
        self._vbundle.set_change_of_frame(self._frame1, self._frame2, auto,
                                          compute_inverse=compute_inverse)
        if compute_inverse:
            self._inverse = type(self)(self._triv2, self._triv1, ~auto)
            self._inverse._inverse = self
        else:
            self._inverse = None

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::



        """
        desc = "Transition map from {} to {}".format(self._triv1, self._triv2)
        return desc

    def _latex_(self):
        r"""
        Return the LaTeX representation of the object.

        TESTS::



        """
        from sage.misc.latex import latex
        return r'(' + latex(self._triv1) + r') \curvearrowright (' + \
               latex(self._triv2) + r')'

    def automorphism(self):
        r"""

        """
        return self._automorphism

    def __mul__(self, other):
        # TODO: Coming soon...
        pass

    def inverse(self):
        r"""

        """
        if self._inverse is None:
            self._vbundle.set_change_of_frame(self._frame2, self._frame12,
                                              ~auto)
            if compute_inverse:
                self._inverse = type(self)(self._triv2, self._triv1, ~auto)
                self._inverse._inverse = self
        return self._inverse

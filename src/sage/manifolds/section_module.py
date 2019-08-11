r"""
Section Modules

The set of continuous local sections over a topological vector bundle `E \to M`
on a domain `U \in M` is a module over the algebra `C^0(U)` of continuous scalar
fields on `U`.

Depending on the domain, there are two classes of section modules:

- :class:`SectionModule` for local sections over a non-trivial part of a
  topological vector bundle
- :class:`SectionFreeModule` for local sections over a trivial part of a
  topological vector bundle

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version from
  vectorfield_module.py
- Michael Jung (2019): Generalization to vector bundles

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.categories.modules import Modules
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.manifolds.section import Section, TrivialSection

class SectionModule(UniqueRepresentation, Parent):
    r"""
    Module of continuous local sections over a topological vector bundle
    `E \to M` on a domain `U \in M`.

    The *section module* `C^0(U;E)` is the set of all local sections of the type

    .. MATH::

        s: U \longrightarrow E

    such that

    .. MATH::

        \forall p \in U,\ s(p) \in E_p,

    where `E_p` is the vector bundle fiber of `E` at the point `p`.

    EXAMPLES:

    TODO

    """
    Element = Section

    def __init__(self, vbundle, domain):
        r"""
        Construct the module of continuous sections over a topological vector
        bundle.

        TESTS::

            TODO

        """
        base_space = vbundle.base_space()
        if not domain.is_subset(base_space):
            raise ValueError("domain must be a subset of base space")
        name = "C^0({};{})".format(domain._name, vbundle._name)
        latex_name = r'C^0({};{})'.format(domain._latex_name,
                                          vbundle._latex_name)
        self._name = name
        self._latex_name = latex_name
        self._vbundle = vbundle
        self._domain = domain
        self._base_space = vbundle.base_space()
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring,
                        category=Modules(self._ring))
        # Dictionary of the tensor modules built on self
        #   (keys = (k,l) --the tensor type)
        # This dictionary is to be extended on need by the method tensor_module
        self._tensor_modules = {(1,0): self} # self is considered as the set
                                             # of tensors of type (1,0)
        # Dictionaries of exterior powers of self and of its dual
        #   (keys = p --the power degree)
        # These dictionaries are to be extended on need by the methods
        # exterior_power and dual_exterior_power
        #self._exterior_powers = {1: self} <- coming soon
        #self._dual_exterior_powers = {} <- coming soon

    #### Begin of parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of the module.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: c_xyz = M.chart()
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U,V)
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module()
            sage: e = E.local_frame('e', domain=U)
            sage: f = E.local_frame('f', domain=V)
            sage: s = C0([-x,y], frame=e, name='s'); s

            TODO

        """
        if isinstance(comp, (int, Integer)) and comp == 0:
            return self.zero()
        if isinstance(comp, Section):
            if self._domain.is_subset(comp._domain):
                return comp.restrict(self._domain)
            else:
                raise ValueError("cannot convert the {} ".format(comp) +
                                 "to a local section in {}".format(self))
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""

        """
        # TODO: Implement
        pass

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent.

        TESTS::

            TODO

        """
        if isinstance(other, (SectionModule, SectionFreeModule)):
            return self._domain.is_subset(other._domain)
        else:
            return False

    #### End of parent methods

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            TODO

        """
        desc = "Module "
        desc += self._name + " "
        desc += "of continuous sections "
        desc += "on " + self._domain._name + " "
        desc += "with values in " + self._vbundle._name
        return desc

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            TODO

        """
        return self._latex_name

    def base_space(self):
        r"""
        Return the base space of the sections in this module.

        OUTPUT:

        - TODO

        EXAMPLES::

            TODO

        """
        return self._base_space

    def domain(self):
        r"""

        """
        return self._domain

    def vector_bundle(self):
        r"""
        TODO
        """
        return self._vbundle

    @cached_method
    def zero(self):
        """
        Return the zero of ``self``.

        EXAMPLES::

            TODO

        """
        elt = self.element_class(self, name='zero', latex_name='0')
        for frame in self._vbundle._frames:
            if frame._domain.is_subset(self._domain):
                elt.add_comp(frame)
                # (since new components are initialized to zero)
        return elt

#******************************************************************************

class SectionFreeModule(FiniteRankFreeModule):
    r"""
    TODO
    """
    Element = TrivialSection

    def __init__(self, vbundle, domain):
        r"""
        Construct the free module of sections over a trivialized vector bundle.

        TESTS::

            TODO

        """
        from .scalarfield import ScalarField
        self._domain = domain
        name = "C^0({};{})".format(domain._name, vbundle._name)
        latex_name = r'C^0({};{})'.format(domain._latex_name,
                                          vbundle._latex_name)
        base_space = vbundle.base_space()
        self._base_space = base_space
        self._vbundle = vbundle
        cat = Modules(domain.scalar_field_algebra()).FiniteDimensional()
        FiniteRankFreeModule.__init__(self, domain.scalar_field_algebra(),
                               vbundle.rank(), name=name,
                               latex_name=latex_name,
                               start_index=base_space._sindex,
                               output_formatter=ScalarField.coord_function,
                               category=cat)

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of ``self``.

        TESTS::

            TODO

        """
        if isinstance(comp, (int, Integer)) and comp == 0:
            return self.zero()
        if isinstance(comp, Section):
            if self._domain.is_subset(comp._domain):
                return comp.restrict(self._domain)
            else:
                raise ValueError("cannot convert the {}".format(comp) +
                                 "to a local section in {}".format(self))
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp != []:
            resu.set_comp(basis)[:] = comp
        return resu

    # Rem: _an_element_ is declared in the superclass FiniteRankFreeModule

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from parent ``other``.

        TESTS::

            TODO

        """
        if isinstance(other, (SectionModule, SectionFreeModule)):
            return self._domain.is_subset(other._domain)
        else:
            return False

    #### End of parent methods

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            TODO

        """
        desc = "Free module "
        desc += self._name + " "
        desc += "of continuous sections "
        desc += "on " + self._domain._name + " "
        desc += "with values in " + self._vbundle._name
        return desc

    def domain(self):
        r"""
        Return the domain on which ``self`` is defined.

        TODO

        """
        return self._domain

    def base_space(self):
        r"""
        TODO
        """
        return self._base_space

    def vector_bundle(self):
        r"""
        TODO
        """
        return self._vbundle

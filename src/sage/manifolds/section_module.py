r"""
Section Modules



AUTHORS:

- Michael Jung (2019): initial version

"""

#******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.modules import Modules
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.manifolds.section import Section, TrivialSection

class SectionModule(UniqueRepresentation, Parent):
    r"""

    """
    Element = Section
    
    def __init__(self, vbundle, domain=None):
        r"""

        """
        base_space = vbundle.base_space()
        if domain is None:
            domain = base_space
        elif not domain.is_subset(base_space):
            raise ValueError("domain must be a subset of base space")
        name = "C^0({};{})".format(domain._name, vbundle._name)
        latex_name = r'C^0({};{})'.format(domain._latex_name, vbundle._latex_name)
        self._name = name
        self._latex_name = latex_name
        self._vbundle = vbundle
        self._domain = domain
        self._base_space = vbundle.base_space()
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring,
                        category=Modules(self._ring))

    #### Begin of parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""

        """
        # TODO: Implement
        pass

    def _an_element_(self):
        r"""

        """
        # TODO: Implement
        pass

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent.

        TESTS::



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



        """
        vbundle = self._vbundle
        desc = "Module "
        desc += self._name + " "
        desc += "of {} sections ".format(self._domain._structure.name)
        desc += "over the {} ".format(self._domain)
        desc += "with values on the " + self._domain._structure.name + " "
        desc += vbundle.base_field_type() + " "
        desc += "vector bundle "
        desc += vbundle._name + " -> " + self.base_space()._name + " "
        desc += "of rank {} ".format(vbundle.rank())
        return desc

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

        """
        return self._latex_name

    def base_space(self):
        r"""
        Return the base space of the sections in this module.

        OUTPUT:

        -

        EXAMPLES::


        """
        return self._base_space

    def domain(self):
        r"""

        """
        return self._domain

    def vector_bundle(self):
        r"""

        """
        return self._vbundle

#******************************************************************************

class SectionFreeModule(FiniteRankFreeModule):
    r"""

    """

    Element = TrivialSection

    def __init__(self, vbundle, domain=None):
        r"""
        Construct the free module of sections over a trivialized vector bundle.

        TESTS::



        """
        from .scalarfield import ScalarField
        self._domain = domain
        name = "C^0({};{})".format(domain._name, vbundle._name)
        latex_name = r'C^0({};{})'.format(domain._latex_name, vbundle._latex_name)
        base_space = vbundle.base_space()
        self._base_space = base_space
        self._vbundle = vbundle
        cat = Modules(domain.scalar_field_algebra()).FiniteDimensional()
        FiniteRankFreeModule.__init__(self, domain.scalar_field_algebra(),
                               base_space._dim, name=name,
                               latex_name=latex_name,
                               start_index=base_space._sindex,
                               output_formatter=None,
                               category=cat)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::



        """
        desc = "Free module "
        desc += self._name + " "
        desc += "of {} sections ".format(self._domain._structure.name)
        desc += "over the {} ".format(self._domain)
        desc += "with values on the {}".format(self._vbundle)
        return desc

    def domain(self):
        r"""
        Return the domain on which ``self`` is defined.


        """
        return self._domain

    def base_space(self):
        r"""

        """
        return self._base_space

    def vector_bundle(self):
        r"""

        """
        return self._vbundle
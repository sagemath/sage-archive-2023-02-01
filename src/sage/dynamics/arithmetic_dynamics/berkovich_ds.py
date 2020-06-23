r"""
Dynamical systmes on Berkovich space over \CC_p.

A dynamical system on Berkovich space over \CC_p is 
determined by a dynamical system on A^1(\CC_p) or P^1(\CC_p),
which naturally induces a dynamical system on affine or
projective Berkovich space.
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Element
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sage.misc.classcall_metaclass import typecall
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.rings.padics.generic_nodes import is_pAdicField
from sage.schemes.berkovich.berkovich_space import Berkovich_Cp_Affine, Berkovich_Cp_Projective
from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective

class DynamicalSystem_Berkovich(Element):

    @staticmethod
    def __classcall_private__(cls, system_morphism_polys, domain=None):
        r"""
        Base class for dynamical systems on Berkovich space over \CC_p.

        INPUT::

        - ``system_morphism_polys`` -- a list of polynomials or rational functions,
          or a dynamical system. In any case, this input should be
          defined over affine or projective space of a finite extension of \QQ_p.

        - ``domain`` -- (optional) affine or projective Berkovich space over \CC_p
        """
        if not(domain is None) or 
        if isinstance(system_morphism_polys,(list,tuple)):
            polys = list(system_morphism_polys)
            if len(system_morphism_polys) not in [1,2]:
                raise ValueError('system_morphism_polys must be length 1 or 2')
            from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing
            from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
            test = lambda x: is_PolynomialRing(x) or is_MPolynomialRing(x)
            if not all(test(poly.parent()) for poly in polys):
                try:
                    polys = [poly.lift() for poly in polys]
                except AttributeError:
                    raise ValueError('{} must be elements of a polynomial ring'.format(system_morphism_polys))
            if domain is None:
                PR = get_coercion_model().common_parent(*polys)
                if (not is_pAdicField(PR)) or :
                    raise ValueError('common parent for polynomials not a padic field')
                system_morphism_polys = DynamicalSystem_projective(polys)
            else:
                
                if isinstance(domain, Berkovich_Cp_Affine):


        if isinstance(system_morphism_polys, DynamicalSystem):
            system_morphism_polys = system_morphism_polys.as_scheme_morphism()

        if isinstance(system_morphism_polys, SchemeMorphism_polynomial):
            R = system_morphism_polys.base_ring()
            morphism_domain = system_morphism_polys.domain()
            polys = list(system_morphism_polys)
            if morphism_domain != system_morphism_polys.codomain():
                raise ValueError('domain and codomain do not agree')
            if not is_pAdicField(R):
                if domain is None:
                    raise ValueError('system_morphism_polys not defined over padic field and domain is None')
                try:
                    system_morphism_polys = system_morphism_polys.change_ring()
            if morphism_domain != morphism_domain.ambient_space():
                raise ValueError('morphism must be defined on the ambient space')
            if domain == None:
                domain = Berkovich_Cp_Projective(morphism_domain)
            else:
                from sage.schemes.berkovich.berkovich_space import Berkovich_Cp
                if not isinstance(domain, Berkovich_Cp):
                    raise ValueError('domain must be Berkovich space over Cp')
            return typecall(cls,polys,domain)

        


"""
Structure maps for finite fields

This module provides classes for isomorphisms between finite fields and vector spaces.

AUTHORS:

- Kwankyu Lee (2017-11-07): initial version

"""

#*****************************************************************************
#       Copyright (C) 2017 Kwankyu <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism

class FiniteFieldIsomorphism(SetMorphism):
    """
    Base class of the vector space isomorphism between a finite field
    and a vector space over a subfield of the finite field.
    """
    def _repr_(self):
        """
        Return the string representation of this isomorphism
        between a finite field and a vector space.

        EXAMPLES::

            sage: E = GF(16)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi
            Isomorphism:
              From: Vector space of dimension 1 over Finite Field in z4 of size 2^4
              To:   Finite Field in z4 of size 2^4

        """
        s = "Isomorphism:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s

    def is_injective(self):
        """
        EXAMPLES::

            sage: E = GF(9)
            sage: F = GF(3)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        EXAMPLES::

            sage: E = GF(9)
            sage: F = GF(3)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi.is_surjective()
            True
        """
        return True

class MorphismVectorSpaceToFiniteField(FiniteFieldIsomorphism):
    """
    Isomorphisms from vector spaces to finite fields.
    """
    def __init__(self, V, K, function):
        """
        Initialize.

        INPUT:

        - ``V`` -- vector space

        - ``K`` -- finite field

        - ``function`` -- Python function that inputs a vector and outputs
          an element of the finite field

        EXAMPLES::

            sage: E = GF(16)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: phi
            Isomorphism:
              From: Vector space of dimension 1 over Finite Field in z4 of size 2^4
              To:   Finite Field in z4 of size 2^4
        """
        FiniteFieldIsomorphism.__init__(self, Hom(V, K), function)

class MorphismFiniteFieldToVectorSpace(FiniteFieldIsomorphism):
    """
    Isomorphisms from finite fields to vector spaces
    """
    def __init__(self, K, V, function):
        """
        Initialize.

        INPUT:

        - ``K`` -- finite field

        - ``V`` -- vector space

        - ``function`` -- Python function that inputs an element of the finite field
          and outputs a vector of the vector space

        EXAMPLES::

            sage: E = GF(16)
            sage: F = GF(4)
            sage: V, phi, psi = E.vector_space(E, map=True)
            sage: psi
            Isomorphism:
              From: Finite Field in z4 of size 2^4
              To:   Vector space of dimension 1 over Finite Field in z4 of size 2^4
        """
        FiniteFieldIsomorphism.__init__(self, Hom(K, V), function)

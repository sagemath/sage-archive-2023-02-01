r"""
Field embedding

Considering a *big field* `F_{q^m}` and a *small_field* `F_q`, with
`q = p^s`, `p` being a prime and `s, m` being integers, this file
contains a method to take care of the representation of `F_{q^m}`-elements
as `F_q`-elements.
"""

#*****************************************************************************
#       Copyright (C) 2016 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.functions.all import log
from sage.structure.sage_object import SageObject
from sage.categories.homset import Hom
from sage.matrix.constructor import column_matrix
from sage.modules.free_module_element import vector

class FieldEmbedding(SageObject):
    r"""
    Represents the embedding of a big non prime field into a smaller
    non prime field.

    INPUT:

    - ``big_field``, ``small_field`` -- two finite fields, ``small_field``
      being a subfield of ``big_field``

    - ``embedding`` -- (default: ``None``) an homomorphism from ``small_field`` to
      ``big_field``. If ``None`` is provided, it will default to the first
      homomorphism of the list of homomorphisms Sage can build.

    EXAMPLES::

        sage: from sage.coding.field_embedding import *
        sage: Fqm.<aa> = GF(16)
        sage: Fq.<a> = GF(4)
        sage: FieldEmbedding(Fqm, Fq)
        Embedding between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2

    It is possible to specify the embedding to use
    from ``big_field`` to ``small_field``::

        sage: Fqm.<aa> = GF(16)
        sage: Fq.<a> = GF(4)
        sage: FieldEmbedding(Fqm, Fq, embedding=Hom(Fq, Fqm)[1])
        Embedding between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
    """

    def __init__(self, big_field, small_field, embedding=None):
        r"""
        TESTS:

        If ``big_field`` is not a finite field, an error is raised::

            sage: from sage.coding.field_embedding import *
            sage: Fqm = RR
            sage: Fq.<a> = GF(4)
            sage: FieldEmbedding(Fqm, Fq)
            Traceback (most recent call last):
            ...
            ValueError: big_field has to be a finite field

        Same for ``small_field``::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq = RR
            sage: FieldEmbedding(Fqm, Fq)
            Traceback (most recent call last):
            ...
            ValueError: small_field has to be a finite field

        If ``small_field`` is not a subfield of ``big_field``, an exception
        is raised::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(8)
            sage: FieldEmbedding(Fqm, Fq)
            Traceback (most recent call last):
            ...
            ValueError: small_field has to be a subfield of big_field
        """
        if not big_field.is_finite():
            raise ValueError("big_field has to be a finite field")
        if not small_field.is_finite():
            raise ValueError("small_field has to be a finite field")
        p = small_field.characteristic()
        s = log(small_field.order(), p)
        sm = log(big_field.order(), p)
        if not s.divides(sm):
            raise ValueError("small_field has to be a subfield of big_field")
        H = Hom(small_field, big_field)
        if embedding is not None and not embedding in H:
            raise ValueError("embedding has to be an embedding from small_field to big_field")
        elif embedding is not None:
            self._phi = embedding
        else:
            self._phi = H[0]
        self._prime_field = small_field.base_ring()
        self._small_field = small_field
        self._big_field = big_field
        alpha = small_field.gen()
        beta = big_field.gen()
        self._alphas = [alpha ** i for i in range(s)]
        self._betas = [beta ** i for i in range(sm)]
        self._small_field_power = s
        self._big_field_power = sm

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FieldEmbedding(Fqm, Fq)
            Embedding between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
        """
        return "Embedding between %s and %s" % (self.big_field(), self.small_field())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: latex(FieldEmbedding(Fqm, Fq))
            \textnormal{Embedding between \Bold{F}_{2^{4}} and \Bold{F}_{2^{2}}}
        """
        return "\\textnormal{Embedding between %s and %s}" % (self.big_field()._latex_(),
                self.small_field()._latex_())

    @cached_method
    def representation_matrix(self):
        r"""
        Returns the matrix used to represent `b` as an element of the
        small field.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.representation_matrix()
            [1 0 0 0]
            [0 0 1 1]
            [0 1 1 1]
            [0 0 0 1]
        """
        s = self.small_field_power()
        m = self.big_field_power() / s
        betas = self.big_field_basis()
        phi_alphas = [ self._phi(self._alphas[i]) for i in range(s) ]
        A = column_matrix([vector(betas[i] * phi_alphas[j])
            for i in range(m) for j in range(s)])
        return A.inverse()

    def small_field_vector_representation(self, b):
        r"""
        Returns a vector representation of ``b`` in the basis of
        the small field over the base field.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: b = aa^3 + aa^2 + aa + 1
            sage: FE.small_field_vector_representation(b)
            (1, 0, 1, 1)
        """
        return self.representation_matrix() * vector(b)

    def small_field_polynomial_representation(self, b):
        r"""
        Returns a polynomial representation of ``b`` in the basis of
        the small field over the base field.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: b = aa^3 + aa^2 + aa + 1
            sage: FE.small_field_polynomial_representation(b)
            a
        """
        Fq = self.small_field()
        vect = self.representation_matrix() * vector(b)
        pol = Fq.zero()
        s = self.small_field_power()
        sm = self.big_field_power()
        for i in range(0, sm, s):
            pol += Fq(vect[i:i+s])
        return pol

    def big_field_representation(self, a):
        r"""
        Returns a polynomial representation of ``a`` over the big field.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: v = vector(GF(2), [1, 0, 1, 1])
            sage: FE.big_field_representation(v)
            aa^3 + aa^2 + aa + 1
        """
        alphas = self.small_field_basis()
        betas = self.big_field_basis()
        phi = self.embedding()
        s = self.small_field_power()
        m = self.big_field_power() / s
        b = self.big_field().zero()
        for i in range(m):
            b += betas[i] * phi(sum([a[j] * alphas[j%s] for j in range(i*s, i*s + s)]))
        return b

    def embedding(self):
        r"""
        Returns the embedding which is used to go from the
        small field to the big field.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.embedding()
            Ring morphism:
             From: Finite Field in a of size 2^2
             To:   Finite Field in aa of size 2^4
             Defn: a |--> aa^2 + aa
        """
        return self._phi

    def small_field_basis(self):
        r"""
        Returns a basis of the small field over the prime field.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.small_field_basis()
            [1, a]
        """
        return self._alphas

    def big_field_basis(self):
        r"""
        Returns a basis of the big field over the prime field.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.big_field_basis()
            [1, aa, aa^2, aa^3]
        """
        return self._betas

    def small_field_power(self):
        r"""
        Let `F_p` be the base field of our small field `F_q`.
        Returns `s` where `p^s = q`

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.small_field_power()
            2
        """
        return self._small_field_power

    def big_field_power(self):
        r"""
        Let `F_p` be the base field of our big field `F_{q^m}`.
        Returns `sm` where `p^{sm} = q^{m}`

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.big_field_power()
            4
        """
        return self._big_field_power

    def prime_field(self):
        r"""
        Returns the base field of our big and small fields.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.prime_field()
            Finite Field of size 2
        """
        return self._prime_field

    def small_field(self):
        r"""
        Returns the small field of ``self``.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.small_field()
            Finite Field in a of size 2^2
        """
        return self._small_field

    def big_field(self):
        r"""
        Returns the big field of ``self``.

        EXAMPLES::

            sage: from sage.coding.field_embedding import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = FieldEmbedding(Fqm, Fq)
            sage: FE.big_field()
            Finite Field in aa of size 2^4
        """
        return self._big_field

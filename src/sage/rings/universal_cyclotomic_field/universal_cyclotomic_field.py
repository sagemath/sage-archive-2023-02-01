r"""
The universal cyclotomic field (UCF)

Implementation of the universal cyclotomic field using the :meth:`Zumbroich basis<UniversalCyclotomicField.zumbroich_basis_indices>`.
The universal cyclotomic field is the smallest subfield of the complex field
containing all roots of unity.

REFERENCES:

.. [Bre97] T. Breuer "Integral bases for subfields of cyclotomic fields" AAECC 8, 279--289 (1997).

AUTHORS:

- Christian Stump

.. NOTE::

    - This function behaves exactly like the *Cyclotomics* in *GAP*.
    - The universal cyclotomic field is used to work with non-crystallographic
      reflection groups. E.g., to work with elements as matrices, computing
      *reflecting hyperplanes*, and *characters*.
    - To multiply matrices over the universal cyclotomic field, it is still
      *much* faster to coerce it to a cyclotomic field and to the
      computation there.

.. TODO::

    - implementation of matrices over the universal cyclotomic field.
    - speed improvements of the cythonized methods.
    - speed improvements for scalar multiples.
    - Remove the inheritance from Field and FieldElement as soon as
      the methods ``is_field(proof=True)`` is implemented in the Fields category.


EXAMPLES:

The universal cyclotomic field is constructed using::

    sage: UCF = UniversalCyclotomicField(); UCF
    Universal Cyclotomic Field

One can as well construct it through :func:`~sage.ring.number_field.number_field.CyclotomicField`::

    sage: UCF = CyclotomicField(); UCF
    Universal Cyclotomic Field

The cyclotomics themselves are accessable through::

    sage: UCF.gen(5)
    E(5)
    sage: UCF.gen(5,2)
    E(5)^2

or the alias::

    sage: UCF.gen(5)
    E(5)
    sage: UCF.gen(5,2)
    E(5)^2

One can as well access the universal cyclotomic field using::

    sage: UCF.<E> = UniversalCyclotomicField();
    sage: E(5)
    E(5)

Other names are supported as well::

    sage: UCF.<zeta> = UniversalCyclotomicField();
    sage: zeta(5)
    zeta(5)

As are other bracketings::

    sage: UCF.<E> = UniversalCyclotomicField(bracket='');
    sage: E(5)
    E5

    sage: UCF.<E> = UniversalCyclotomicField(bracket="[]");
    sage: E(5)
    E[5]

    sage: UCF.<E> = UniversalCyclotomicField(bracket="(ABCXYZ)");
    sage: E(5)
    E(ABC5XYZ)

We use the generator "E" and the standard bracketing throughout this file::

    sage: UCF.<E> = UniversalCyclotomicField();

Some very first examples::

    sage: E(2)
    -1
    sage: E(3)
    E(3)
    sage: E(6)
    -E(3)^2

Equality and inequality checks::

    sage: E(6,2) == E(6)^2 == E(3)
    True

    sage: E(6)^2 != E(3)
    False

Addition and multiplication::

    sage: E(2) * E(3)
    -E(3)
    sage: f = E(2) + E(3); f
    2*E(3) + E(3)^2

Inverses::

    sage: f^-1
    1/3*E(3) + 2/3*E(3)^2
    sage: f.inverse()
    1/3*E(3) + 2/3*E(3)^2
    sage: f * f.inverse()
    1

Complex conjugation::

    sage: f.conjugate()
    E(3) + 2*E(3)^2

Galois conjugation::

    sage: f.galois_conjugates()
    [2*E(3) + E(3)^2, E(3) + 2*E(3)^2]
    sage: f.norm_of_galois_extension()
    3

Coercion to the algebraic field :class:`QQbar<sage.rings.qqbar.AlgebraicField>`::

    sage: QQbar(E(3))
    -0.500000000000000? + 0.866025403784439?*I
    sage: QQbar(f)
    -1.500000000000000? + 0.866025403784439?*I

Partial conversion to the real algebraic field :class:`AA<sage.rings.qqbar.AlgebraicRealField>`::

    sage: AA(E(5)+E(5).conjugate())
    0.618033988749895?

    sage: AA(E(5))
    Traceback (most recent call last):
    ...
    TypeError: No conversion of E(5) to the real algebraic field AA.

One can as well define the universal cyclotomic field without any embedding::

    sage: UCF.<E> = UniversalCyclotomicField(embedding=None); UCF
    Universal Cyclotomic Field

    sage: UCF.<E> = UniversalCyclotomicField(embedding=False); UCF
    Universal Cyclotomic Field

    sage: QQbar(E(5))
    Traceback (most recent call last):
    ...
    TypeError: Illegal initializer for algebraic number

Conversion to :class:`CyclotomicField<sage.rings.number_field.number_field.CyclotomicField>`:

.. WARNING::

    This is only possible if ``self`` has the standard embedding

::

    sage: UCF.<E> = UniversalCyclotomicField()

    sage: E(5).to_cyclotomic_field()
    zeta5

    sage: f = E(2) + E(3)
    sage: f.to_cyclotomic_field()
    zeta3 - 1

    sage: CF = CyclotomicField(5)
    sage: CF(E(5))
    zeta5

    sage: CF = CyclotomicField(7)
    sage: CF(E(5))
    Traceback (most recent call last):
    ...
    TypeError: The element E(5) cannot be converted to Cyclotomic Field of order 7 and degree 6

    sage: CF = CyclotomicField(10)
    sage: CF(E(5))
    zeta10^2

Conversions to and from GAP::

    sage: a = gap('E(6)'); a
    -E(3)^2
    sage: a.parent()
    Gap

    sage: b = UCF.from_gap(a); b
    -E(3)^2
    sage: b.parent()
    Universal Cyclotomic Field

    sage: gap(b)
    -E(3)^2

Conversions to and from the *cyclotomic field*::

    sage: a = E(6).to_cyclotomic_field(); a
    zeta3 + 1

    sage: UCF.from_cyclotomic_field(a)
    -E(3)^2

One can also do basic arithmetics with matrices over the universal cyclotomic field::

    sage: m = matrix(2,[E(3),1,1,E(4)]); m
    [E(3)    1]
    [   1 E(4)]
    sage: m.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Universal Cyclotomic Field

    sage: m^2
    [                       -E(3) E(12)^4 - E(12)^7 - E(12)^11]
    [E(12)^4 - E(12)^7 - E(12)^11                            0]

    sage: -m
    [-E(3)    -1]
    [   -1 -E(4)]

And compute its *characteristic polynomial*, *echelon form*, *pivots*, and thus its *rank*::

    sage: m.charpoly()
    x^2 + (-E(12)^4 + E(12)^7 + E(12)^11)*x + E(12)^4 + E(12)^7 + E(12)^8

    sage: m.echelon_form()
    [1 0]
    [0 1]

    sage: m.pivots()
    (0, 1)

    sage: m.rank()
    2

The eigenvalues do not (yet) work::

    sage: m.eigenvalues() # not implemented
    ...
    NotImplementedError:

A long real life test. Computing ``N3`` is much faster than computing
``N2`` which is again 3 times faster than computing ``N1``::

    sage: W = gap3.ComplexReflectionGroup(14)       #optional - gap3 # long time
    sage: UC = W.UnipotentCharacters()              #optional - gap3 # long time
    sage: UCF.<E> = UniversalCyclotomicField();     #optional - gap3 # long time
    sage: M = matrix(UCF,UC.families[2].fourierMat) #optional - gap3 # long time
    sage: N1 = M*M                                  #optional - gap3 # long time

    sage: N2 = UCF._matrix_mult(M,M)                #optional - gap3 # long time
    sage: CF = CyclotomicField(24)                  #optional - gap3 # long time
    sage: M = matrix(CF,M)                          #optional - gap3 # long time
    sage: N3 = matrix(UCF,M*M)                      #optional - gap3 # long time
    sage: N1 == N2 == N3                            #optional - gap3 # long time
    True

TESTS:

As an indication that everything works, we start with a test that we
obtain the same answers as in GAP::

    sage: all(str(E(n,k)).translate(None,' ') == gap.execute('E('+str(n)+')^'+str(k)).translate(None,'\n ') for n in range(1,15) for k in range(n))
    True

The following didn't work first::

    sage: str(-E(9)^4-E(9)^7).translate(None,' ') == gap.execute('-E(9)^4-E(9)^7').translate(None,'\n ')
    True
    sage: str(-E(9)^5-E(9)^8).translate(None,' ') == gap.execute('-E(9)^5-E(9)^8').translate(None,'\n ')
    True
"""
#*****************************************************************************
#       Copyright (C) 2012 Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from random import randint, randrange, sample, choice

import sage.structure.parent_base
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import FieldElement, Element
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.sage_object import have_same_parent

from sage.categories.morphism import SetMorphism
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
from sage.categories.sets_cat import Sets
from sage.categories.homset import Hom

from sage.rings.all import ZZ, QQ
from sage.rings.ring import Field
from sage.rings.qqbar import QQbar, AA
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.integer import GCD_list, LCM_list

from sage.rings.real_mpfr import RealField, mpfr_prec_min
from sage.rings.complex_field import ComplexField
from sage.rings.complex_double import CDF

from sage.rings.real_lazy import RLF, CLF
from sage.combinat.dict_addition import dict_linear_combination, dict_addition
from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field_c import \
    ZumbroichBasisCython, push_down_cython, ZumbroichDecomposition, galois_conjugates_cython, \
    push_to_higher_field, dict_multiplication, dict_vector_multiplication

class UniversalCyclotomicField(UniqueRepresentation, Field):
    r"""
    The *universal cyclotomic field*, which is the smallest field containing
    the rational numbers together with all roots of unity.

    Its elements are represented as linear combinations of the so-called
    *Zumbroich basis*.

    EXAMPLES::

        sage: UCF.<E> = UniversalCyclotomicField(); UCF
        Universal Cyclotomic Field
        sage: E(12)
        -E(12)^7
        sage: E(12) in UCF
        True

    One also has access to the universal cyclotomic field using the function :func:`CyclotomicField`::

        sage: UCF = CyclotomicField(); UCF
        Universal Cyclotomic Field

    One can also construct a vector space over the universal cyclotomic field::

        sage: UCF^3
        Vector space of dimension 3 over Universal Cyclotomic Field
    """
    def __init__(self,names="E",bracket="()",embedding=True):
        r"""

        :param names: The string name for the root of unity
        :type names: optional, default:"E"
        :param bracket: The bracket used in string representations. Can be any even length string
        :type bracket: optional, default:"()"
        :param embedding: The given embedding of ``self`` into the complex numbers. Can be
            - ``True``: The standard embedding
            - ``None`` or ``False``: No embedding
            - An indexable `f(n) \in \mathbb{N}` such that `0 \leq f(n) < n` sending the generator ``E(n)`` to `e^{ 2 \pi i k / n }`.
        :param embedding: optional, default:``True``

        .. WARNING::

            The optional argument ``embedding`` is not checked for consistency.

        TESTS::

            sage: F = UniversalCyclotomicField()
            sage: TestSuite(F).run()
        """
        # the data is stored as linear combinations of elements in the Zumbroich Basis
        from sage.combinat.free_module import CombinatorialFreeModule
        from sage.categories.fields import Fields
        from sage.categories.algebras import Algebras

        # getting the optional argument "names" right
        if names is None:
            names = "E"
        if not isinstance(names,str):
            if not isinstance(names,(list,tuple)):
                raise ValueError("The given name %s is not valid."%names)
            if len(names) != 1:
                raise ValueError("The given name %s is not valid."%names)
            names = names[0]
        if not isinstance(names,str):
            raise ValueError("The given name %s is not valid."%names)

        # getting the optional argument "bracket" right
        if isinstance(bracket,str) and len(bracket) % 2 == 0:
            bracket_len = len(bracket)
            bracket = (bracket[:bracket_len//2],bracket[bracket_len//2:])
        else:
            raise ValueError("The given bracket %s is not a string of even length."%bracket)

        self._data = CombinatorialFreeModule(QQ, ZumbroichBasisIndices(),prefix=names,bracket=bracket)
        Parent.__init__(self, base = QQ, category = (Fields(), Algebras(QQ)))

        self._has_standard_embedding = False
        if embedding not in [False,None]:
            if embedding is True:
                embedding = lambda n: QQbar.zeta(n)
                self._has_standard_embedding = True

            P = get_parent_of_embedding(embedding)

           # embedding from self into P
            def on_basis(x):
                return embedding(x[0])**(int(x[1]))

            mor = self._data.module_morphism(on_basis, codomain=P)
            H = SetMorphism(Hom(self, QQbar), lambda z: mor(z.value))
            self.register_embedding(H)

            # define partial conversions to AA, if possible
            SetMorphism(
                    Hom(self, AA, SetsWithPartialMaps()),
                    lambda elem: elem._real_()
               ).register_as_conversion()

        # string representations of elements
        def repr_term(m,type="repr"):
            if m[1] == 0:
                return '1'
            elif 2*m[1] == m[0]:
                return '-1'
            elif m[1] == 1:
                if type == "repr":
                    return '%s%s%s%s'%(names,bracket[0],m[0],bracket[1])
                elif type == "latex":
                    return '\\zeta_{%s}'%m[0]
            else:
                if type == "repr":
                    return '%s%s%s%s^%s'%(names,bracket[0],m[0],bracket[1],m[1])
                elif type == "latex":
                    return '\\zeta_{%s}^{%s}'%(m[0],m[1])

        self._data._repr_term = lambda m: repr_term(m,type="repr")
        self._data._latex_term = lambda m: repr_term(m,type="latex")

        # setting zero and one
        self._zero = self._from_dict({}, remove_zeros=False)
        self._one = self._from_dict({(1,0):QQ(1)}, remove_zeros=False)

    @cached_method
    def gen(self,n, k=1):
        r"""
        Returns `\zeta^k` living in :class:`UniversalCyclotomicField`, where `\zeta` denotes the primitive `n`-th root of unity `\zeta = exp(2 \pi i / n)`.

        :param n: positive integer.
        :param k: positive integer.
        :type n: integer
        :type k: integer; optional, default ``1``

        .. NOTE::

            - For the mathematical description of the Zumbroich basis and the
              algorithmic behind, see [Bre97]_.

            - This function behaves exactly like the *Cyclotomics* in *GAP*.

        EXAMPLES::

            sage: UCF.<E> = UniversalCyclotomicField()

            sage: E(3) # indirect doctest
            E(3)
            sage: E(6) # indirect doctest
            -E(3)^2
            sage: E(12) # indirect doctest
            -E(12)^7
            sage: E(6,2) # indirect doctest
            E(3)
            sage: E(6)^2 # indirect doctest
            E(3)
        """
        if n != ZZ(n) or not n > 0 or k != ZZ(k):
            raise TypeError('The argument for a root of unity is not correct.')
        else:
            g = GCD_list([n,k])
            n = ZZ(n/g)
            k = ZZ(k/g) % n

        return self._from_dict(push_down_cython(n, ZumbroichDecomposition(n,k)), coerce=True, remove_zeros=False)

    def _first_ngens(self,n):
        r"""
        Returns the method :meth:`gen` if ``n=1``, and raises an error otherwise.

        This method is needed to make the following work::

            sage: UCF.<E> = UniversalCyclotomicField() # indirect doctest
        """
        if n == 1:
            return (self.gen,)
        else:
            raise ValueError("This ring has only a single generator method.")

    def _element_constructor_(self, arg):
        r"""
        The only way to get here is if there was no coercion found.

        In this case, we give other parents the option to define a conversion
        using the method ``_universal_cyclotomic_field_``, or raise a TypeError otherwise.

        TESTS::

            sage: UCF = UniversalCyclotomicField()

            sage: UCF(CC(5)) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError:  No coercion to the universal cyclotomic field found for the input 5.00000000000000 with parent Complex Field with 53 bits of precision.

            sage: p = Permutation([2,1])
            sage: UCF(p)
            Traceback (most recent call last):
            ...
            TypeError: No coercion to the universal cyclotomic field found for the input [2, 1] with parent Standard permutations.
        """
        if hasattr(arg,"_universal_cyclotomic_field_"):
            return arg._universal_cyclotomic_field_()
        elif isinstance(arg,(sage.interfaces.gap3.GAP3Element,sage.interfaces.gap.GapElement)):
            return self.from_gap(arg)

        error_str = "No coercion to the universal cyclotomic field found for the input %s"%str(arg)
        if hasattr(arg,"parent"):
            error_str ="%s with parent %s."%(error_str,str(arg.parent()))
        else:
            error_str ="%s."%(error_str)
        raise TypeError(error_str)

    def _coerce_map_from_(self, other):
        r"""
        If ``self`` has the standard embedding and
        - ``other`` is a cyclotomic number field: returns
          the coercion thereof, taking non-standard embeddings
          of ``other.parent()`` into account.
        - ``other`` is also a universal cyclotomic field with the
          standard embedding: returns the obvious morphism

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: zeta = CyclotomicField(5).gen()
            sage: UCF(zeta) # indirect doctest
            E(5)

            sage: zeta = CyclotomicField(5,embedding=CC(exp(2*pi*I/5))).gen()
            sage: UCF(zeta) # indirect doctest
            E(5)

            sage: zeta = CyclotomicField(5,embedding=CC(exp(4*pi*I/5))).gen()
            sage: UCF(zeta) # indirect doctest
            E(5)^2

            sage: UCF.<E> = UniversalCyclotomicField();
            sage: UCF2.<E2> = UniversalCyclotomicField();
            sage: UCF2(E(5))
            E2(5)

            sage: UCF = UniversalCyclotomicField(embedding=None)
            sage: UCF(zeta) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: No coercion to the universal cyclotomic field found for the input zeta5 with parent Cyclotomic Field of order 5 and degree 4.
        """
        from sage.rings.number_field.number_field import NumberField_cyclotomic
        from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
        if self._has_standard_embedding:
            if isinstance(other,NumberField_cyclotomic):
                return  NumberFieldEmbedding(other, self, self.from_cyclotomic_field(other.gen()))
            elif isinstance(other,UniversalCyclotomicField) and other._has_standard_embedding:
                return SetMorphism(Hom(other,self), lambda z: self._from_dict(z._dict_()))

    def __pow__(self,n):
        r"""
        Returns the ``n``-th power of self as a vector space.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF^3
            Vector space of dimension 3 over Universal Cyclotomic Field
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self,n)

    def is_finite(self):
        r"""
        Returns False as ``self`` is not finite.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.is_finite()
            False
        """
        return False

    def is_subring(self,other):
        r"""
        Returns True if ``self`` is a subring of ``other``.


        .. WARNING::

            Currently, it is only checked if ``self is other``!

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: UCF.is_subring(UCF)
            True

            sage: UCF.is_subring(CC)
            False
        """
        return other is self

    def _repr_(self):
        r"""
        Returns the string representation of ``self``.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF._repr_()
            'Universal Cyclotomic Field'
        """
        return "Universal Cyclotomic Field"

    def _gap_init_(self):
        r"""
        Returns gap string representation of ``self``.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF._gap_init_()
            'Cyclotomics'
        """
        return 'Cyclotomics'

    def degree(self):
        r"""
        Returns the *degree* of ``self`` as a field extension over the Rationals.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.degree()
            +Infinity
        """
        from sage.rings.infinity import infinity
        return infinity

    def characteristic(self):
        r"""
        Returns ``0`` which is the *characteristic* of ``self``.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.characteristic()
            0
        """
        return ZZ(0)

    def prime_subfield(self):
        r"""
        Returns `\QQ` which is the *prime subfield* of ``self``.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.prime_subfield()
            Rational Field
        """
        return QQ

    def is_prime_field(self):
        r"""
        Returns False since ``self`` is not a prime field.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.is_prime_field()
            False
        """
        return False

    def an_element(self, order = 3):
        r"""
        Returns an element of ``self`` of order ``order``.

        :param order: a positive integer.
        :type order: integer; optional, default:``3``

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: UniversalCyclotomicField().an_element()
            E(3)

            sage: UniversalCyclotomicField().an_element(order=6)
            -E(3)^2

            sage: UniversalCyclotomicField().an_element(order=10)
            -E(5)^3
        """
        return self.gen(order)

    def random_element(self, order=None):
        r"""
        Returns a (often non-trivial) pseudo-random element of ``self``.

        :param order:
        :type order: integer or None; optional, default:``None``

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: UCF.random_element() # random
            3*E(7)^2 + E(7)^3 + 2*E(7)^4 - 5*E(7)^5

            sage: UCF.random_element(order=4) # random
            -3*E(4)

            sage: UCF.random_element(order=12) # random
            E(12)^7 - 4*E(12)^8 + E(12)^11
        """
        F = self
        seq = [-5,-4,-3,-2,-1,1,2,3,4,5]
        if order is not None:
            n = order
        else:
            n = randint(1, 17)
        B = ZumbroichBasisIndices().indices(n)
        # TODO: could this be written in a more conceptual way
        # by having appropriate constructors?
        k = randrange(len(B))
        B = sample(B, k)
        dict_in_basis = {}
        for key in B:
            dict_in_basis[ key.value ] = QQ(choice(seq))
        return F._from_dict(push_down_cython(n,dict_in_basis), remove_zeros=False)

    def _from_dict(self, D, coerce=True, remove_zeros=True):
        r"""
        Returns the element in ``self`` from the given dictionary ``D``.

        :param D: a dictionary with keys being elements in the Zumbroich basis and values being Rationals.
        :param coerce: if True, the values are coerced to the Rationals.
        :type coerce: Boolean; optional, default:``False``
        :param remove_zeros: if True, zeros are removed from the dict first. Should be ``True`` unless it is clear that ``D`` doesn't contain zeros.
        :type remove_zeros: Boolean; optional, default:``True``

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: D = dict([((1,0),2)])
            sage: UCF._from_dict(D)
            2

            sage: D = dict([((1,0),2)],remove_zeros=False)
            sage: UCF._from_dict(D)
            2

            sage: D = dict([((1,0),2),((3,1),1),((3,2),0)])
            sage: UCF._from_dict(D)
            2 + E(3)
        """
        if coerce:
            for X,a in D.iteritems():
                D[X] = QQ(a)
        elem = self.element_class(self, self._data._from_dict(D, remove_zeros=remove_zeros))
        return elem

    def from_base_ring(self,coeff):
        r"""
        Returns the base ring element ``coeff`` as an element in ``self``.

        :param coeff: A rational number.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: x = UCF.from_base_ring(2); x
            2
            sage: x.parent()
            Universal Cyclotomic Field
        """
        return self._from_dict({ (1,0):coeff })

    def zero(self):
        r"""
        Returns the zero in ``self``.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.zero()
            0
        """
        return self._zero

    def one(self):
        r"""
        Returns the one in ``self``.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.one()
            1
        """
        return self._one

    def monomial(self, mon, check=True):
        r"""
        Returns the monomial in ``self`` associated to ``mon`` in the Zumbroich basis.

        :param mon: an element in the Zumbroich basis
        :param check: if True, the monomial is checked to be in the Zumbroich basis
        :type check: Boolean; optional, default:``True``

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: UCF.monomial((1,0))
            1

            sage: UCF.monomial((4,2))
            Traceback (most recent call last):
            ...
            ValueError: The given data is not a monomial of the universal cyclotomic field.
        """
        if check:
            if not mon in ZumbroichBasisIndices():
                raise ValueError("The given data is not a monomial of the universal cyclotomic field.")
        return self._from_dict({ mon : QQ(1) }, remove_zeros=False)

    def sum(self, L):
        r"""
        Returns the sum of all elements (which must be coerceable into ``self``) in ``L``.

        :param L: list or tuple of elements in ``self``

        .. NOTE::

            Faster than the usual sum as operated directly on dictionaries, as all steps are done together.

        EXAMPLES::

            sage: UCF.<E> = UniversalCyclotomicField()

            sage: UCF.sum([ E(i) for i in range(1,5) ])
            E(12)^4 - E(12)^7 - E(12)^11
        """
        l = LCM_list([ other.field_order() for other in L ])
        large_dict_list = [ push_to_higher_field(other.value._monomial_coefficients, other.field_order(), l) for other in L ]
        return self._from_dict(push_down_cython(l,dict_addition(large_dict_list)), remove_zeros=False)

    def _matrix_mult(self,M1,M2,order=None):
        r"""
        Returns the product ``M1`` `\times` ``M2`` of the two matrices ``M1`` and ``M2`` over ``self``.

        .. WARNING::

            This method is not for public use, but only to provide a quick test how fast we can multiply matrices.

        EXAMPLES::

            sage: UCF.<E> = UniversalCyclotomicField()

            sage: M = matrix(UCF,[[E(3),E(4)],[E(5),E(6)]]); M
            [   E(3)    E(4)]
            [   E(5) -E(3)^2]

            sage: M2 = UCF._matrix_mult(M,M); M2
            [-E(60)^4 - E(60)^7 - E(60)^16 - E(60)^28 - E(60)^47 - E(60)^52                                             E(12)^7 - E(12)^11]
            [                                            E(15)^8 - E(15)^13 -E(60)^7 - E(60)^8 - E(60)^32 - E(60)^44 - E(60)^47 - E(60)^56]

            sage: M2 == M*M
            True
        """
        from sage.matrix.all import zero_matrix
        if not M1.nrows() == M2.ncols():
            raise ValueError("The given matrices cannot be multiplied.")
        dim1, dim, dim2 = M1.ncols(), M1.nrows(), M2.nrows()
        m_rows = M1.rows()
        m_cols = M2.columns()
        rows,cols = [],[]
        for i in xrange(dim):
            rows.append(tuple(x.value._monomial_coefficients for x in m_rows[i]))
            cols.append(tuple(x.value._monomial_coefficients for x in m_cols[i]))

        M_new = zero_matrix(self,dim1,dim2)
        if order:
            n = order
        else:
            LCM = [ x.field_order() for x in set(M1.list()).union(M2.list()) ]
            n = LCM_list(LCM)
        for i in xrange(dim1):
            for j in xrange(dim2):
                M_new[i,j] = self._from_dict(push_down_cython(n,dict_vector_multiplication(n,rows[i],cols[j])))
        return M_new

    def zumbroich_basis_indices(self, n):
        r"""
        Returns the indices of the *Zumbroich basis* of order ``n``.

        The Zumbroich basis is a linear basis of the universal cyclotomic field
        that behaves very well with considering an primitive `d`-th root of unity
        as a (non primitive) `kd`-th root. See [Bre97]_ for further details.

        :param n: positive integer

        OUTPUT:

        - a set of tuples `(n,k)` of all elements in the Zumbroich basis of order `n`.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.zumbroich_basis_indices(8)
            {(8, 1), (8, 3), (8, 0), (8, 2)}
        """
        return ZumbroichBasisIndices().indices(n)

    def zumbroich_basis(self,n):
        r"""
        Returns the *Zumbroich basis* of order ``n``.

        The Zumbroich basis is a linear basis of the universal cyclotomic field
        that behaves very well with considering an primitive `d`-th root of unity
        as a (non primitive) `kd`-th root. See [Bre97]_ for further details.

        :param n: positive integer

        OUTPUT:

        - the set of elements in the universal cyclotomic field forming the Zumbroich basis of order `n`.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.zumbroich_basis(8)
            {E(8)^3, E(4), E(8), 1}

            sage: UCF.zumbroich_basis(9)
            {E(9)^5, E(9)^4, E(3)^2, E(3), E(9)^7, E(9)^2}
        """
        return set(self.gen(n,k) for n,k in self.zumbroich_basis_indices(n))

    def from_gap(self, elem):
        r"""
        Returns the element in ``self`` obtained from the gap by executing ``string``.

        :param string: string representing an element in the universal cyclotomic field

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.from_gap(gap("-E(3)^2"))
            -E(3)^2

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.from_gap(gap("E(3)^2"))
            E(3)^2

            sage: UCF.from_gap(gap("1/6*E(3)")) # testing a former bug
            1/6*E(3)
        """
        if not isinstance(elem,(sage.interfaces.gap3.GAP3Element,sage.interfaces.gap.GapElement)):
            raise ValueError("The input %s is not a GAP object."%elem)
        if hasattr(elem,"sage"):
            try:
                return self(elem.sage())
            except NotImplementedError:
                pass

        string = str(elem)
        terms = string.replace('\n','').replace('-','+-').split('+')
        if terms[0] == '':
            del terms[0]
        for i in range(len(terms)):
            if '^' in terms[i]:
                terms[i] = terms[i].replace(')^',',') + ')'
            terms[i] = terms[i].replace("E","self.gen")
            if '*' in terms[i]:
                pos = terms[i].index('*')
                fac = QQ(terms[i][:pos])
                terms[i] = terms[i][pos+1:]
            else:
                fac = None
            exec('terms[%s]='%i + terms[i])
            if fac is not None:
                terms[i] = fac*terms[i]
        return self.sum(terms)

    def from_cyclotomic_field(self, elem):
        r"""
        Returns the element in ``self`` coming from the element in NumberField_cyclotomic.

        :param elem: an element of NumberField_cyclotomic

        .. WARNING::

            This method raises an error if self does not have the standard embedding.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: a = CyclotomicField(6).gen(); a
            zeta6
            sage: UCF.from_cyclotomic_field(a)
            -E(3)^2

        An example with another embedding::

            sage: a = CyclotomicField(5,embedding=CC(exp(4*pi*I/5))).gen(); a
            zeta5
            sage: UCF.from_cyclotomic_field(a)
            E(5)^2

        TESTS::

            sage: UCF.from_cyclotomic_field(4)
            Traceback (most recent call last):
            ...
            TypeError: The given data (4) is not a cyclotomic field element.

            sage: UCF = UniversalCyclotomicField(embedding=None);
            sage: a = CyclotomicField(5).gen()
            sage: UCF.from_cyclotomic_field(a)
            Traceback (most recent call last):
            ...
            TypeError: This method can only be used if Universal Cyclotomic Field uses the standard embedding.
        """
        from sage.rings.number_field.number_field import NumberField_cyclotomic
        if not self._has_standard_embedding:
            raise TypeError("This method can only be used if %s uses the standard embedding."%self)
        if not hasattr(elem,'parent') or not isinstance(elem.parent(), NumberField_cyclotomic):
            raise TypeError("The given data (%s) is not a cyclotomic field element."%elem)
        n = elem.parent()._n()
        CF = CyclotomicField(n)
        elem = CF(elem)
        coeff_list = elem.list()
        return self._from_dict(push_down_cython(n,dict_linear_combination((ZumbroichDecomposition(n, k), coeff_list[ k ]) for k in range(len(coeff_list)) if coeff_list[k] != 0)), remove_zeros=False)

    # This cannot inherit from FieldElement and ElementWrapper, this causes the _add_()
    #   _sub_(), _mul_(), etc. to not call the corresponding functions in this class,
    #   leading to a segfault.
    class Element(FieldElement):
        r"""
        An element of the universal cyclotomic field.

        .. SEEALSO::

            - :class:`UniversalCyclotomicField`
        """
        def __init__(self, parent, value):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: elt = E(6)
                sage: TestSuite(elt).run()
            """
            self.value = value
            FieldElement.__init__(self, parent)

        def _repr_(self):
            """
            Return a string representation of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(6)
                -E(3)^2
            """
            return repr(self.value)

        def _latex_(self):
            r"""
            Return a `\LaTeX` representation of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: latex(E(6))
                -\zeta_{3}^{2}
            """
            from sage.misc.latex import latex
            return latex(self.value)

        def __lt__(self,other):
            r"""
            Pushes the method forward to ``QQbar``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3).__lt__(E(4))
                True
            """
            return QQbar(self).__lt__(other)

        def __gt__(self,other):
            r"""
            Pushes the method forward to ``QQbar``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3).__gt__(E(4))
                False
            """
            return QQbar(self).__gt__(other)

        def __le__(self,other):
            r"""
            Pushes the method forward to ``QQbar``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3).__le__(E(4))
                True
            """
            return QQbar(self).__le__(other)

        def __ge__(self,other):
            r"""
            Pushes the method forward to ``QQbar``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3).__ge__(E(4))
                False
            """
            return QQbar(self).__ge__(other)

        def __eq__(self, other):
            r"""
            Returns ``True`` if self and other are equal.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3) == E(3) # indirect doctest
                True

                sage: E(3) == E(4) # indirect doctest
                False
            """
            if have_same_parent(self, other):
                return self.value._monomial_coefficients == other.value._monomial_coefficients
            from sage.structure.element import get_coercion_model
            import operator
            try:
                return get_coercion_model().bin_op(self, other, operator.eq)
            except TypeError:
                return False

        def __ne__(left, right):
            r"""
            Returns ``True`` if self and other are not equal.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3) != E(3) # indirect doctest
                False

                sage: E(3) != E(4) # indirect doctest
                True
            """
            return not left.__eq__(right)

        # only implemented to survive the equality check for matrices
        def __cmp__(self, other):
            r"""
            The ordering is the one on the underlying sorted list of (monomial,coefficients) pairs.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: cmp(E(3),E(4)) # indirect doctest
                -1

                sage: cmp(UCF(3), 1) # indirect doctest
                1
            """
            if self.__eq__(other):
                return 0
            if self.field_order() > other.field_order():
                return 1
            if self.field_order() < other.field_order():
                return -1
            return cmp(sorted(self), sorted(other))

        def __hash__(self):
            r"""
            Returns the hash of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: L = set( E(n) for n in [1..1000] )
                sage: all( E(n) in L for n in [1..1000] )
                True
            """
            return hash(frozenset(self.value._monomial_coefficients.items()))

        def _gap_init_(self):
            r"""
            Returns gap string representation of ``self``.

            EXAMPLES::

            sage: UCF.<X> = UniversalCyclotomicField();
            sage: X(4)
            X(4)
            sage: X(4)._gap_init_()
            'E(4)'
            """
            UCF = UniversalCyclotomicField()
            return str(UCF._from_dict(self.value._monomial_coefficients))

        def _real_(self):
            r"""
            Returns ``self`` in the real algebraic field, if possible. Raises an error otherwise.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: x = E(5)
                sage: AA(E(5)+E(5).conjugate()) # indirect doctest
                0.618033988749895?

                sage: AA(E(5)) # indirect doctest
                Traceback (most recent call last):
                ...
                TypeError: No conversion of E(5) to the real algebraic field AA.
            """
            if self.is_real():
                P = self.parent().coerce_embedding().codomain()
                return AA(P(self))
            raise TypeError("No conversion of %s to the real algebraic field AA."%str(self))

        def __float__(self):
            r"""
            Returns ``self`` as a float if ``self`` is real. Raises an error otherwise.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: float(E(5)+E(5)^(-1))
                0.6180339887498949

                sage: float(E(5))
                Traceback (most recent call last):
                ...
                ValueError: E(5) is not real
            """
            if self.is_real():
                return float(CDF(self).real_part())
            raise ValueError("{} is not real".format(self))

        def _rational_(self):
            r"""
            Returns ``self`` in the Rationals, if possible. Raises an error otherwise.

            EXAMPLES::

                sage: UCF = UniversalCyclotomicField()

                sage: x = UCF(4/5)._rational_(); x
                4/5
                sage: x.parent()
                Rational Field

                sage: UCF(0)._rational_()
                0

                sage: UCF.gen(5)._rational_()
                Traceback (most recent call last):
                ...
                TypeError:  No conversion of E(5) to the rational field QQ.
            """
            if self.is_zero():
                return QQ.zero()
            elif self.is_rational():
                return self.value._monomial_coefficients[(1,0)]
            else:
                raise TypeError("No conversion of %s to the rational field QQ."%str(self))

        def _integer_(self):
            r"""
            Returns ``self`` in the Integers, if possible. Raises an error otherwise.

            EXAMPLES::

                sage: UCF = UniversalCyclotomicField()

                sage: x = UCF(5)._integer_(); x
                5
                sage: x.parent()
                Integer Ring
            """
            if self.is_rational():
                coeff = self.value._monomial_coefficients[(1,0)]
                return ZZ(coeff)
            raise TypeError("No conversion of %s to the integer ring ZZ."%str(self))

        def __nonzero__(self):
            r"""
            Returns True if ``self`` is not zero.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).is_zero()
                False

                sage: (E(3)-E(3)).is_zero()
                True
            """
            return bool(self.value._monomial_coefficients)

        def is_one(self):
            r"""
            Returns True if ``self`` is one.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).is_one()
                False

                sage: UCF(1).is_one()
                True
            """
            try:
                x = self.value._monomial_coefficients[(1,0)]
                return x == 1
            except KeyError:
                return False

        def is_rational(self):
            r"""
            Returns True if ``self`` is rational.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).is_rational()
                False

                sage: UCF(1/3).is_rational()
                True

                sage: UCF(0).is_rational()
                True
            """
            return (1,0) in self.value._monomial_coefficients or self.is_zero()

        def is_real(self):
            r"""
            Returns True if ``self`` is real.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(5).is_real()
                False

                sage: (E(5)^2 + E(5)^3).is_real()
                True

                sage: (E(5)^4 + E(5)^3).is_real()
                False
            """
            return self == self.conjugate()

        def is_real_positive(self):
            r"""
            Returns True if ``self`` is real and positive.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: (E(5)^2 + E(5)^3).is_real_positive()
                False

                sage: (-E(5)^4 - E(5)^3).is_real_positive()
                False
            """
            return self.is_real() and AA(self) > 0

        def __neg__(self):
            r"""
            Returns the negative of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3).__neg__()
                -E(3)
            """
            F = self.parent()
            elem = F.element_class(F, self.value.__neg__())
            return elem

        def __invert__(self):
            r"""
            Returns the inverse of ``self``.

            .. NOTE::

                We make use of the fact that the norm of ``self``, i.e., the product of all Galois conjugates of ``self``, is rational.
                Thus the inverse of ``self`` is the product of all Galois conjugates except ``self`` multiplied by the inverse of the norm of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).__invert__()
                E(3)^2

                sage: E(6).__invert__()
                -E(3)

                sage: UCF(0).__invert__()
                Traceback (most recent call last):
                ...
                ZeroDivisionError: The given element is zero.
            """
            if self.is_zero():
                raise ZeroDivisionError("The given element is zero.")
            elif self.is_rational():
                return self.parent().from_base_ring(QQ(1)/self.value._monomial_coefficients[(1,0)])
            else:
                inv = self.parent().prod(self.galois_conjugates()[1:])
                self_norm_inv = QQ(1)/(self*inv).value._monomial_coefficients[(1,0)]
                return inv*self_norm_inv

        def __iter__(self):
            r"""
            Returns an iterator of the monomial coefficients of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: a = E(3)+E(4); a
                E(12)^4 - E(12)^7 - E(12)^11

                sage: for x in a: print x
                ((12, 7), -1)
                ((12, 4), 1)
                ((12, 11), -1)
            """
            return self.value.__iter__()

        def _dict_(self):
            r"""
            Returns the *monomial coefficients* of ``self`` as a dictionary.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(6)._dict_()
                {(3, 2): -1}
            """
            return self.value._monomial_coefficients

        def inverse(self):
            r"""
            Returns the inverse of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: f = 2 * E(3) + E(4); f
                2*E(12)^4 - E(12)^7 - E(12)^11

                sage: f.inverse()
                2/13*E(12)^4 - 3/13*E(12)^7 + 8/13*E(12)^8 + 1/13*E(12)^11

                sage: f * f.inverse()
                1
            """
            return self.__invert__()

        def __div__(self, other):
            r"""
            Returns ``self `\times` ``other.inverse()``.

            EXAMPLES::

            sage: UCF.<E> = UniversalCyclotomicField()

            sage: E(6).__div__(1)
            -E(3)^2

            sage: E(6).__div__(-1)
            E(3)^2

            sage: E(6).__div__(E(3))
            -E(3)
            """
            return self * other**-1

        def __pow__(self, k):
            r"""
            Returns ``self`` `{}^k`.

            :param k: an integer

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).__pow__(3)
                1
                sage: E(6).__pow__(4)
                E(3)^2

                sage: E(6).__pow__(-2)
                E(3)^2
                sage: E(6).__pow__(0)
                1
            """
            if k == 0:
                return self.parent().one()
            elif k == 1:
                return self
            elif self.is_zero():
                return self.parent().zero()
            elif self.is_rational():
                return self.parent()._from_dict({ (1,0) : self.value._monomial_coefficients[(1,0)]**k }, remove_zeros=False)
            elif len(self.value._monomial_coefficients) == 1:
                mon,coeff = self.value._monomial_coefficients.iteritems().next()
                n = self.field_order()
                return self.parent()._from_dict(push_down_cython(n,dict_linear_combination([ (ZumbroichDecomposition(n, k*mon[1] % n), coeff**k,) ])), remove_zeros=False)
            elif k < 0:
                return self.__invert__().__pow__(-k)
            else:
                if k % 2 == 0:
                    return (self*self).__pow__(k/2)
                else:
                    return self * self.__pow__(k-1)

        def _acted_upon_(self, scalar, self_on_left = False):
            r"""
            Returns the action of a scalar on ``self``.

            :param scalar: the rational factor by which ``self`` gets multiplied
            :param self_on_left: if ``True``, the multiplication is done on the left
            :type self_on_left: Boolean; optional, default:``True``

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3)._acted_upon_(2)
                2*E(3)
            """
            # With the current design, the coercion model does not have
            # enough information to detect apriori that this method only
            # accepts scalars; so it tries on some elements(), and we need
            # to make sure to report an error.
            if hasattr(scalar, 'parent') and scalar.parent() != self.base_ring():
                # Temporary needed by coercion (see Polynomial/FractionField tests).
                if self.base_ring().has_coerce_map_from(scalar.parent()):
                    scalar = self.base_ring()(scalar)
                else:
                    return None

            F = self.parent()
            elem = F.element_class(F, self.value._acted_upon_(scalar, self_on_left = self_on_left))
            return elem

        # For backward compatibility
        _lmul_ = _acted_upon_
        _rmul_ = _acted_upon_

        def _mul_(self, other):
            r"""
            Returns ``self`` `\times` ``other``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3)._mul_(E(4))
                E(12)^7
            """
            F = self.parent()
            if other.parent() is not F:
                other = F(other)

            if self.is_zero():
                return self
            if other.is_zero():
                return other
            if self.is_one():
                return other
            if other.is_one():
                return self

            D_self = self.value._monomial_coefficients
            D_other = other.value._monomial_coefficients

            if (1,0) in D_self:
                coeff = D_self[(1,0)]
                return F._from_dict(dict_linear_combination([ (D_other, coeff) ]), remove_zeros=False)
            elif (1,0) in D_other:
                coeff = D_other[(1,0)]
                return F._from_dict(dict_linear_combination([ (D_self, coeff) ]), remove_zeros=False)

            n1,n2 = self.field_order(),other.field_order()
            n = LCM_list([n1,n2])
            return F._from_dict(push_down_cython(n,dict_multiplication(D_self, D_other, n1, n2, n)), remove_zeros=False)

        def _sub_(self, other):
            r"""
            Returns ``self`` `+` ``other.__neg__()``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3).__sub__(E(4))
                E(12)^4 + E(12)^7 + E(12)^11
            """
            return self + other.__neg__()

        def _add_(self, other):
            r"""
            Returns ``self`` `+` ``other``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(3) + E(4) # indirect doctest
                E(12)^4 - E(12)^7 - E(12)^11
            """
            n,m = self.field_order(),other.field_order()
            l = LCM_list([ n, m ])

            D_self = push_to_higher_field(self.value._monomial_coefficients, n, l)
            D_other = push_to_higher_field(other.value._monomial_coefficients, m, l)

            return self.parent()._from_dict(push_down_cython(l,dict_addition([D_self,D_other])), remove_zeros=False)

        def minpoly(self, var='x'):
            r"""
            The minimal polynomial of ``self`` element over `\QQ`.

            :param var: the minimal polynomial is defined over a polynomial ring
               in a variable with this name

            :type var: optional, default:``'x'``

            .. SEEALSO::

                - :meth:`~sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic.minpoly`

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: UCF(4).minpoly()
                x - 4

                sage: UCF(4).minpoly(var='y')
                y - 4

                sage: E(3).minpoly()
                x^2 + x + 1

                sage: E(3).minpoly(var='y')
                y^2 + y + 1

            TESTS::

                sage: x = UCF(4)
                sage: x.minpoly() == x.to_cyclotomic_field().minpoly()
                True

                sage: x = E(3)
                sage: x.minpoly() == x.to_cyclotomic_field().minpoly()
                True

                sage: x = E(3)
                sage: x.minpoly(var='y') == x.to_cyclotomic_field().minpoly(var='y')
                True
            """
            if self.is_rational():
                R = QQ[var]
                return R([-self._rational_(), 1])
            else:
                return self.to_cyclotomic_field().minpoly(var=var)

        def abs(self):
            r"""
            Return the absolute value of this element.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).abs()
                1

                sage: x = 2*E(3); x.abs()
                2

                sage: x = E(3)+E(4); x.abs()
                1.931851652578137?

            If no embedding is given, an error is raised::

                sage: UCF.<E> = UniversalCyclotomicField(embedding=None)
                sage: E(3).abs()
                Traceback (most recent call last):
                ...
                ValueError: Universal Cyclotomic Field has no embedding defined.
            """
            P = self.parent().coerce_embedding()
            if P is not None:
                P = P.codomain()
                return P(self).abs()
            raise ValueError("%s has no embedding defined."%str(self.parent()))

        def conjugate(self):
            r"""
            Returns the complex conjugate of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).conjugate()
                E(3)^2

                sage: E(4).conjugate()
                -E(4)

            .. NOTE::

                the conjugate of a monomial is always a monomial or the negation thereof.
            """
            n = self.field_order()
            return self.parent()._from_dict(dict_linear_combination((ZumbroichDecomposition(n, n - key[1]), c) for (key,c) in self), remove_zeros=False)

        def galois_conjugates(self, m=None):
            r"""
            Returns all Galois conjugates of ``self``.

            Those are the elements in the universal cyclotomic field obtained
            from self by substituting `\zeta_n` by `\zeta_n^k` for all `{\rm
            gcd}(n,k)=1`. Remark that for odd primes, the Galois conjugates
            permutes the Zumbroich basis. The first Galois conjugate in the
            list is ``self``.

            :param m: if given, it must be a multiple of :meth:`field_order`;
                the Galois conjugates are then computed with respect to the cyclotomics of order ``m``

            OUTPUT:

            - a list `[p_{i_1},...,p_{i_{max}}]`, where `p_{i_j}` is obtained
              from ``self`` by substituting `E(n)` by `E(n)^{i_j}` and where
              `i_j` is the `j`-th integer coprime to n

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(6).galois_conjugates()
                [-E(3)^2, -E(3)]

                sage: E(6).galois_conjugates(6)
                [-E(3)^2, -E(3)]

                sage: E(6).galois_conjugates(12)
                [-E(3)^2, -E(3), -E(3)^2, -E(3)]

                sage: E(8).galois_conjugates()
                [E(8), E(8)^3, -E(8), -E(8)^3]

                sage: E(8).galois_conjugates(16)
                [E(8), E(8)^3, -E(8), -E(8)^3, E(8), E(8)^3, -E(8), -E(8)^3]

                sage: E(9).galois_conjugates()
                [-E(9)^4 - E(9)^7, E(9)^2, E(9)^4, E(9)^5, E(9)^7, -E(9)^2 - E(9)^5]

                sage: E(11).galois_conjugates()
                [E(11), E(11)^2, E(11)^3, E(11)^4, E(11)^5, E(11)^6, E(11)^7, E(11)^8, E(11)^9, E(11)^10]

                sage: E(6).galois_conjugates(5)
                Traceback (most recent call last):
                ...
                ValueError: The given integer (5) is not a multiple of the field order of -E(3)^2.
            """
            n = self.field_order()
            if m is None:
                m = n
            else:
                if not m%n == 0:
                    raise ValueError("The given integer (%s) is not a multiple of the field order of %s."%(m,self))
            coprimes = [1] + [ i for i in range(2,m) if GCD_list([m,i])==1 ]
            conjugates = galois_conjugates_cython(self.value._monomial_coefficients, n, m, coprimes)
            return [ self.parent()._from_dict(conjugate, remove_zeros=False) for conjugate in conjugates ]

        def support(self):
            r"""
            Returns the support of ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()
                sage: E(6).support()
                [(3, 2)]
            """
            return self.value.support()

        def coefficient(self, mon):
            r"""
            Returns the coefficient of ``mon`` in ``self``.

            :param mon: an element in the Zumbroich basis

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(6)
                -E(3)^2

                sage: E(6).coefficient((3,2))
                -1

                sage: E(6).coefficient((3,1))
                0

            Alternatively, one can use indexed access::

                sage: E(6)[(3,2)]
                -1
            """
            try:
                return self.value._monomial_coefficients[mon]
            except KeyError:
                return 0

        __getitem__ = coefficient

        def field_order(self):
            r"""
            Returns the order of the smallest field containing ``self``.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(4).field_order()
                4

                sage: E(6).field_order()
                3
            """
            if bool(self.value._monomial_coefficients):
                return self.value._monomial_coefficients.iterkeys().next()[0]
            else:
                return 1

        def norm_of_galois_extension(self):
            r"""
            Returns the norm as a Galois extension of `\QQ`, which is
            given by the product of all galois_conjugates.

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(3).norm_of_galois_extension()
                1

                sage: E(6).norm_of_galois_extension()
                1

                sage: (E(2) + E(3)).norm_of_galois_extension()
                3
            """
            return self.parent().prod(self.galois_conjugates()).value._monomial_coefficients[(1,0)]

        def to_cyclotomic_field(self):
            r"""
            Returns ``self`` in :class:`CyclotomicField`.

            .. WARNING::

                This method raises an error if ``self.parent()`` does not
                have the standard embedding

            EXAMPLES::

                sage: UCF.<E> = UniversalCyclotomicField()

                sage: E(5).to_cyclotomic_field()
                zeta5

            This method is as well used to convert to a cyclotomic field::

                sage: CF = CyclotomicField(5)
                sage: CF(E(5))
                zeta5

            .. SEEALSO::

                :class:`CyclotomicField`
            """
            if not self.parent()._has_standard_embedding:
                raise TypeError("This method can only be used if %s uses the standard embedding."%self)
            CF = CyclotomicField(self.field_order())
            zeta = CF.gen()
            def on_basis(x):
                return zeta**(x[1])
            return self.parent()._data._apply_module_morphism(self.value, on_basis, codomain=CF)

class ZumbroichBasisIndices(UniqueRepresentation, Parent):
    def __init__(self):
        r"""
        This class is a thin wrapper to work with indices in the Zumbroich basis.

        EXAMPLES::

            sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import ZumbroichBasisIndices

            sage: ZumbroichBasisIndices()
            The indices of the Zumbroich basis

        One can ask for an element to play with::

            sage: a = ZumbroichBasisIndices().an_element(); a
            (12, 4)

        The element ``a`` is indeed an element of this class::

            sage: a.parent()
            The indices of the Zumbroich basis

        And one can check if an element is indeed contained in the Zumbroich basis::

            sage: a in ZumbroichBasisIndices()
            True

            sage: (12,4) in ZumbroichBasisIndices()
            True
        """
        Parent.__init__(self, category=Sets())

    def _repr_(self):
        r"""
        Returns the string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import ZumbroichBasisIndices
            sage: ZumbroichBasisIndices() # indirect doctest
            The indices of the Zumbroich basis
        """
        return "The indices of the Zumbroich basis"

    def an_element(self):
        r"""
        Returns an element of the Zumbroich basis.

        EXAMPLES::

            sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import ZumbroichBasisIndices
            sage: a = ZumbroichBasisIndices().an_element(); a
            (12, 4)
        """
        return self.element_class(self, (12, 4))

    def __contains__(self, x):
        r"""
        Returns ``True`` if ``x`` is contained in ``self``.

        EXAMPLES::

            sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import ZumbroichBasisIndices
            sage: a = ZumbroichBasisIndices().an_element(); a
            (12, 4)
            sage: a in ZumbroichBasisIndices() # indirect doctest
            True
            sage: (12,4) in ZumbroichBasisIndices()
            True
            sage: (12,5) in ZumbroichBasisIndices()
            False
        """
        if isinstance(x,tuple) and len(x) == 2:
            n,i = int(x[0]),int(x[1])
            x = self.element_class( self, (n, int(i)) )
        return x in self.indices(x[0])

    def indices(self, n, m=1):
        r"""
        Returns the list of tuples `(n,k)` such that the set `\zeta_n^k` form a Zumbroich basis for `QQ(\zeta_n)` over `QQ(\zeta_m)`.

        :param n: positive integer
        :param m: positive integer dividing ``n``
        :type m: optional, default:``1``

        EXAMPLES::

            sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import ZumbroichBasisIndices

            sage: ZumbroichBasisIndices().indices(6)
            {(6, 4), (6, 2)}
            sage: ZumbroichBasisIndices().indices(12)
            {(12, 7), (12, 4), (12, 11), (12, 8)}
            sage: ZumbroichBasisIndices().indices(24)
            {(24, 19), (24, 8), (24, 17), (24, 16), (24, 14), (24, 1), (24, 22), (24, 11)}
        """
        if not n%m == 0:
            raise ValueError('%s does not divide %s.'%(m,n))
        B = ZumbroichBasisCython(n, m)
        return set([self.element_class( self, (n, int(i)) ) for i in B])

    class Element(ElementWrapper):
        def __getitem__(self,i):
            r"""
            Passes the indexing to its value.

            EXAMPLES::

            sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import ZumbroichBasisIndices
            sage: a = ZumbroichBasisIndices().an_element()
            sage: a[0],a[1] # indirect doctest
            (12, 4)
            """
            return self.value[i]

def get_parent_of_embedding(embedding):
    r"""
    Returns the parent of an element in the image of ``embedding``.

    :param embedding: A function from the positive integers `\{1,2,3,\ldots\}` into a common parent

    If the images are in a real or complex field, then
    it creates an image into a lazy field.

    EXAMPLES::

        sage: from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import get_parent_of_embedding

        sage: get_parent_of_embedding(lambda n: QQbar.zeta()^n)
        Algebraic Field

        sage: get_parent_of_embedding(lambda n: CC(exp(2*pi*I/n)))
        Complex Lazy Field
    """
    if callable(embedding) and isinstance(embedding(3), Element):
        P = embedding(3).parent()
        if not P.is_exact():
            RR = RealField(mpfr_prec_min())
            CC = ComplexField(mpfr_prec_min())
            if RR.has_coerce_map_from(P):
                P = RLF
            elif CC.has_coerce_map_from(P):
                P = CLF
        return P
    else:
        raise TypeError("Embedding (type %s) must be an element." % type(embedding))

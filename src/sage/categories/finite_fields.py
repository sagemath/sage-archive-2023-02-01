r"""
Finite Fields
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom

class FiniteFields(CategoryWithAxiom):
    """
    The category of finite fields.

    EXAMPLES::

        sage: K = FiniteFields()
        sage: K
        Category of finite fields

    A finite field is a finite monoid with the structure of a field::

        sage: K.super_categories()
        [Category of fields, Category of finite commutative rings]

    Some examples of membership testing and coercion::

        sage: FiniteField(17) in K
        True
        sage: RationalField() in K
        False
        sage: K(RationalField())
        Traceback (most recent call last):
        ...
        TypeError: unable to canonically associate a finite field to Rational Field

    TESTS::

        sage: TestSuite(FiniteFields()).run()
        sage: FiniteFields().is_subcategory(FiniteEnumeratedSets())
        True
    """

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in FiniteFields()
            True
            sage: QQ in FiniteFields()
            False
            sage: IntegerModRing(4) in FiniteFields()
            False
        """
        from sage.categories.fields import Fields
        return x in Fields() and x.is_finite()

    # As is, this does no more than the usual __call__ of Category, but for the error message
    def _call_(self, x):
        """
        EXAMPLES::

            sage: FiniteFields()(GF(4, "a"))
            Finite Field in a of size 2^2
            sage: FiniteFields()(RationalField())   # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unable to canonically associate a finite field to Rational Field
        """
        raise TypeError("unable to canonically associate a finite field to %s"%x)
        # TODO: local dvr ring?

    #@lazy_attribute
    #def element_class(self):
    #    """
    #    A common super class for all elements of finite fields
    #
    #    EXAMPLES::
    #
    #        sage: C = FiniteFields().element_class; C
    #        <type 'sage.rings.finite_rings.element_base.FiniteFieldElement'>
    #        sage: type(C)
    #        <type 'type'>
    #    """
    #    from sage.rings.finite_rings.element_base import FiniteFieldElement
    #    return FiniteFieldElement

    class ParentMethods:
        pass

    class ElementMethods:
        pass

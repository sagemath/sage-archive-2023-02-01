r"""
Monoid algebras
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def MonoidAlgebras(base_ring):
    """
    The category of monoid algebras over ``base_ring``

    EXAMPLES::

        sage: C = MonoidAlgebras(QQ); C
        Category of monoid algebras over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of algebras with basis over Rational Field,
         Category of semigroup algebras over Rational Field,
         Category of unital magma algebras over Rational Field]

    This is just an alias for::

        sage: C is Monoids().Algebras(QQ)
        True

    TESTS::

        sage: TestSuite(MonoidAlgebras(ZZ)).run()
    """
    from sage.categories.all import Monoids
    return Monoids().Algebras(base_ring)

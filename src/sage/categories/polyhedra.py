r"""
Polyhedral subsets of free ZZ, QQ or RR-modules.
"""
#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring

class PolyhedralSets(Category_over_base_ring):
    r"""
    The category of polyhedra over a ring.

    EXAMPLES:

    We create the category of polyhedra over `\QQ`::

        sage: PolyhedralSets(QQ)
        Category of polyhedral sets over Rational Field

    TESTS::

        sage: TestSuite(PolyhedralSets(RDF)).run()

        sage: P = Polyhedron()
        sage: P.parent().category().element_class
        <class 'sage.categories.polyhedra.PolyhedralSets.element_class'>
        sage: P.parent().category().element_class.mro()
        [<class 'sage.categories.polyhedra.PolyhedralSets.element_class'>,
         <class 'sage.categories.magmas.Magmas.element_class'>,
         <class 'sage.categories.additive_magmas.AdditiveMagmas.element_class'>,
         <class 'sage.categories.sets_cat.Sets.element_class'>,
         <class 'sage.categories.sets_with_partial_maps.SetsWithPartialMaps.element_class'>,
         <class 'sage.categories.objects.Objects.element_class'>,
         <type 'object'>]
        sage: isinstance(P, P.parent().category().element_class)
        True
    """

    def __init__(self, R):
        """
        TESTS::

            sage: PolyhedralSets((1,2,3))
            Traceback (most recent call last):
            ...
            TypeError: base ring R (=(1, 2, 3)) must be ZZ, QQ, or RDF.
        """
        from sage.rings.all import ZZ, QQ, RDF
        if R not in [ZZ, QQ, RDF]:
            raise TypeError, 'base ring R (='+str(R)+') must be ZZ, QQ, or RDF.'
        Category_over_base_ring.__init__(self, R)

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: PolyhedralSets(QQ).super_categories()
            [Category of magmas, Category of additive magmas]
        """
        from sage.categories.all import Magmas, AdditiveMagmas
        return [Magmas(), AdditiveMagmas()]





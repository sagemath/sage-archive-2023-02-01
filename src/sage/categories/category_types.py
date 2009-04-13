"""
Specific category classes.

This is placed in a separate file from categories.py to avoid circular imports
(as morphisms must be very low in the hierarchy with the new coercion model).
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from category import *
from sage.rings.ring import is_Ring
from sage.rings.field import is_Field
from sage.algebras.algebra import is_Algebra




####################################################################
#   Different types of categories
####################################################################

#############################################################
# Category of elements of some object
#############################################################
class Elements(uniq1, Category):
    """
    The category of all elements of a given object.

    EXAMPLES:
        sage: a = IntegerRing()(5)
        sage: C = a.category(); C
        Category of elements of Integer Ring
        sage: a in C
        True
        sage: 2/3 in C
        False
        sage: loads(C.dumps()) == C
        True
    """
    def __init__(self, object):
        Category.__init__(self)
        self.__object = object

    def _call_(self, x):
        """
        EXAMPLES:
            sage: V = VectorSpace(QQ,3)
            sage: x = V.0
            sage: C = x.category()
            sage: C
            Category of elements of Vector space of dimension 3 over Rational Field
            sage: w = C([1,2,3]); w
            (1, 2, 3)
            sage: w.category()
            Category of elements of Vector space of dimension 3 over Rational Field
        """
        return self.__object(x)

    def object(self):
        return self.__object

    def __reduce__(self):
        return Elements, (self.__object, )

    def _repr_(self):
        return "Category of elements of %s"%self.object()

    def _latex_(self):
        """
        EXAMPLES:
            sage: V = VectorSpace(QQ,3)
            sage: x = V.0
            sage: latex(x.category())
            \mathbf{Elt}_{\mathbf{Q}^{3}}
        """
        return "\\mathbf{Elt}_{%s}"%latex(self.__object)


#############################################################
# Category of sequences of elements of some objects
#############################################################
class Sequences(uniq1, Category):
    """
    The category of all elements of a given object.

    EXAMPLES:
        sage: v = Sequence([1,2,3]); v
        [1, 2, 3]
        sage: C = v.category(); C
        Category of sequences in Integer Ring
        sage: loads(C.dumps()) == C
        True
        sage: Sequences(ZZ) is C
        True
    """
    def __init__(self, object):
        Category.__init__(self)
        self.__object = object

    def _call_(self, x):
        """
        EXAMPLES:
            sage: v = Sequence([1,2,3]); v
            [1, 2, 3]
            sage: C = v.category(); C
            Category of sequences in Integer Ring
            sage: w = C([2/1, 3/1, GF(2)(5)]); w
            [2, 3, 1]
            sage: w.category()
            Category of sequences in Integer Ring
        """
        from sage.structure.sequence import Sequence
        return Sequence(x, self.__object)

    def object(self):
        return self.__object

    def __reduce__(self):
        return Sequences, (self.__object, )

    def _repr_(self):
        return "Category of sequences in %s"%self.object()

    def _latex_(self):
        """
        EXAMPLES:
            sage: v = Sequence([1,2,3])
            sage: latex(v.category())
            \mathbf{Seq}_{\mathbf{Z}}
        """

        return "\\mathbf{Seq}_{%s}"%latex(self.__object)

#############################################################
# Category of objects over some base object
#############################################################
class Category_over_base(uniq1, Category):
    def __init__(self, base, name=None):
        Category.__init__(self, name)
        self.__base = base

    def base(self):
        """
        Return the base over which elements of this category are
        defined.
        """
        return self.__base

    def _repr_(self):
        return Category._repr_(self) + " over %s"%self.__base

    def _latex_(self):
        return "\\mathbf{%s}_{%s}"%(self.__label, latex(self.__base))


#############################################################
# Category of objects over some base ring
#############################################################
class Category_over_base_ring(Category_over_base):
    def __init__(self, base, name=None):
        Category_over_base.__init__(self, base, name)
        if not is_Ring(base):
            raise TypeError, "base must be a ring"
        self.__base = base

    def base_ring(self):
        """
        Return the base ring over which elements of this category are
        defined.
        """
        return self.base()

    def _parameters(self):
        return set([self.__base])


class AbelianCategory:
    def is_abelian(self):
        return True


#############################################################
# Category of objects in some ambient object
#############################################################
class Category_in_ambient(uniq1, Category):
    def __init__(self, ambient, name=None):
        Category.__init__(self, name)
        self.__ambient = ambient

    def ambient(self):
        """
        Return the ambient object in which objects of this category are
        embedded.
        """
        return self.__ambient

    def _repr_(self):
        return Category._repr_(self) + " in %s"%self.__ambient

    def _parameters(self):
        return set([self.__ambient])


class Category_module(Category_over_base_ring, AbelianCategory):
    def __init__(self,  base, name=None):
        Category_over_base_ring.__init__(self, base, name)


class Category_ideal(Category_in_ambient):
    def __init__(self, ring, name=None):
        Category_in_ambient.__init__(self, ring, name)

    def ring(self):
        return self.ambient()

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: C = Ideals(IntegerRing())
            sage: IntegerRing().zero_ideal() in C
            True
        """
        import sage.rings.all
        if sage.rings.all.is_Ideal(x) and x.ring() == self.ring():
            return True
        return False

    def __call__(self, v):
        if v in self:
            return v
        return self.ring().ideal(v)




####################################################################
# Actual categories (below)
####################################################################


#############################################################
# Generic category (default when requesting category of
# an object using misc.functional.category
#############################################################
class Objects(Category_uniq):
    """
    The category of all Sage objects.

    EXAMPLES:
        sage: Objects()
        Category of objects
    """
    def __reduce__(self):
        return Objects, tuple([])

    def __contains__(self, x):
        return True

#############################################################
# Pointed Sets
#############################################################

class PointedSets(Category_uniq):
    """
    The category of pointed sets.

    EXAMPLES:
        sage: PointedSets()
        Category of pointed sets
    """
    def __call__(self, X, pt):
        import sage.sets.all
        return sage.sets.all.Set(X, pt)

    def __reduce__(self):
        return PointedSets, tuple([])

#############################################################
# Sets with Partial Maps
#############################################################

class SetsWithPartialMaps(Category_uniq):
    """
    The category whose objects are sets and whose morphisms are
    maps that are allowed to raise a ValueError on some inputs.

    This category is equivalent to the category of pointed sets,
    via the equivalence sending an object X to X union {error},
    a morphism f to the morphism of pointed sets that sends x
    to f(x) if f does not raise an error on x, or to error if it
    does.

    EXAMPLES:
        sage: SetsWithPartialMaps()
        Category with objects Sets and morphisms partially defined maps
    """
    def __call__(self, X, pt):
        import sage.sets.all
        return sage.sets.all.Set(X, pt)

    def __reduce__(self):
        return SetsWithPartialMaps, tuple([])

    def __repr__(self):
        return "Category with objects Sets and morphisms partially defined maps"

#############################################################
# GSets
#     $G$-Sets play an important role in permutation groups.
#############################################################
class GSets(uniq1, Category):
    """
    The category of $G$-sets, for a group $G$.

    EXAMPLES:
        sage: S = SymmetricGroup(3)
        sage: GSets(S)
        Category of G-sets for SymmetricGroup(3)
    """
    def __init__(self, G):
        Category.__init__(self, "G-sets")
        self.__G = G

    def __repr__(self):
        return "Category of G-sets for %s"%self.__G

    def __reduce__(self):
        """
        EXAMPLES:
            sage: S8 = SymmetricGroup(8)
            sage: C = GSets(S8)
            sage: loads(C.dumps()) == C
            True
        """
        return GSets, (self.__G, )


#############################################################
# SimplicialComplex
#############################################################
class SimplicialComplexes(Category_uniq):
    """
    The category of simplicial complexes.

    EXAMPLES::

        sage: SimplicialComplexes()
        Category of simplicial complexes
    """
    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = SimplicialComplexes()
            sage: loads(C.dumps()) == C
            True
        """
        return SimplicialComplexes, tuple([])

#############################################################
# Semigroup
#############################################################

class Semigroups(Category_uniq):
    """
    The category of all semigroups.

    EXAMPLES:
        sage: Semigroups()
        Category of semigroups
    """
    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = Semigroups()
            sage: loads(C.dumps()) == C
            True
        """
        return Semigroups, tuple([])


#############################################################
# Monoid
#############################################################
class Monoids(Category_uniq):
    """
    The category of monoids.

    EXAMPLES:
        sage: Monoids()
        Category of monoids
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = Monoids()
            sage: loads(C.dumps()) == C
            True
        """
        return Monoids, tuple([])


#############################################################
# Groupoid
#############################################################
class Groupoid(uniq1, Category):
    """
    The category of groupoids, for a set (usually a group) $G$.

    EXAMPLES:
        sage: Groupoid(DihedralGroup(3))
        Groupoid with underlying set Dihedral group of order 6 as a permutation group
    """
    def __init__(self, G):
        Category.__init__(self, "Groupoid")
        self.__G = G

    def __repr__(self):
        return "Groupoid with underlying set %s"%self.__G

    def __reduce__(self):
        """
        EXAMPLES:
            sage: S8 = SymmetricGroup(8)
            sage: C = Groupoid(S8)
            sage: loads(C.dumps()) == C
            True
        """
        return Groupoid, (self.__G, )


#############################################################
# Group
#############################################################
class Groups(Category_uniq):
    """
    The category of groups.

    EXAMPLES:
        sage: Groups()
        Category of groups
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = Groups()
            sage: loads(C.dumps()) == C
            True
        """
        return Groups, tuple([])



#############################################################
# AbelianSemigroup
#############################################################
class AbelianSemigroups(Category_uniq):
    """
    The category of all abelian semigroups.


    EXAMPLES:
        sage: AbelianSemigroups()
        Category of abelian semigroups
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = AbelianSemigroups()
            sage: loads(C.dumps()) == C
            True
        """
        return AbelianSemigroups, tuple([])


#############################################################
# AbelianMonoid
#############################################################
class AbelianMonoids(Category_uniq):
    """
    The category of all monoids.

    EXAMPLES:
        sage: AbelianMonoids()
        Category of abelian monoids
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = AbelianMonoids()
            sage: loads(C.dumps()) == C
            True
        """
        return AbelianMonoids, tuple([])


#############################################################
# AbelianGroup
#############################################################
class AbelianGroups(Category_uniq, AbelianCategory):
    """
    The category of all abelian groups.

    EXAMPLES:
        sage: AbelianGroups()
        Category of abelian groups
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = AbelianGroups()
            sage: loads(C.dumps()) == C
            True
        """
        return AbelianGroups, tuple([])


#############################################################
# Ring
#############################################################
class Rings(Category_uniq):
    """
    The category of all rings.

    EXAMPLES:
        sage: Rings()
        Category of rings
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = Rings()
            sage: loads(C.dumps()) == C
            True
        """
        return Rings, tuple([])


#############################################################
# CommutativeRing
#############################################################
class CommutativeRings(Category_uniq):
    """
    The category of commutative rings.

    EXAMPLES:
        sage: CommutativeRings()
        Category of commutative rings
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = CommutativeRings()
            sage: loads(C.dumps()) == C
            True
        """
        return CommutativeRings, tuple([])


#############################################################
# Field
#############################################################
class Fields(Category_uniq):
    """
    The category of fields.

    EXAMPLES:
        sage: K = Fields()
        sage: K
        Category of fields
        sage: K(IntegerRing())
        Rational Field
        sage: K(PolynomialRing(GF(3), 'x'))
        Fraction Field of Univariate Polynomial Ring in x over
        Finite Field of size 3
        sage: K(RealField())
        Real Field with 53 bits of precision
    """
    def __contains__(self, x):
        return is_Field(x)

    def __call__(self, x):
        if x in self:
            return x
        try:
            return x.fraction_field()
        except AttributeError:
            raise TypeError, "unable to associate a field to %s"%x

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = Fields()
            sage: loads(C.dumps()) == C
            True
        """
        return Fields, tuple([])


#############################################################
# FiniteField
#############################################################
class FiniteFields(Category_uniq):
    """
    The category of all finite fields.

    EXAMPLES:
        sage: K = FiniteFields()
        sage: K
        Category of finite fields
        sage: FiniteField(17) in K
        True
        sage: RationalField() in K
        False
        sage: K(RationalField())
        Traceback (most recent call last):
        ...
        TypeError: unable to canonically associate a finite field to Rational Field
    """
    def __contains__(self, x):
        return is_Field(x) and x.is_finite()

    def __call__(self, x):
        if x in self:
            return x
        raise TypeError, "unable to canonically associate a finite field to %s"%x
        # TODO: local dvr ring?

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = FiniteFields()
            sage: loads(C.dumps()) == C
            True
        """
        return FiniteFields, tuple([])


#############################################################
# NumberField
#############################################################
class NumberFields(Category_uniq):
    r"""
    The category of number fields.

    EXAMPLES:
    We create the category of number fields.
        sage: C = NumberFields()
        sage: C
        Category of number fields

    Notice that the rational numbers $\Q$ *are* considered as
    an object in this category.
        sage: RationalField() in C
        True

    However, we can define a degree 1 extension of $\Q$, which is of
    course also in this category.
        sage: x = PolynomialRing(RationalField(), 'x').gen()
        sage: K = NumberField(x - 1, 'a'); K
        Number Field in a with defining polynomial x - 1
        sage: K in C
        True

    Number fields all lie in this category, regardless of the name
    of the variable.
        sage: K = NumberField(x^2 + 1, 'a')
        sage: K in C
        True
    """
    def __contains__(self, x):
        import sage.rings.all
        return sage.rings.all.is_NumberField(x)

    def __call__(self, x):
        if x in self:
            return x
        try:
            return x.number_field()
        except AttributeError:
            raise TypeError, "unable to canonically associate a number field to %s"%x

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = NumberFields()
            sage: loads(C.dumps()) == C
            True
        """
        return NumberFields, tuple([])


#############################################################
# Algebra
#############################################################
class Algebras(Category_over_base_ring):
    """
    The category of algebras over a fixed base ring.


    """
    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = Algebras(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return Algebras, (self.base(), )

    def __contains__(self, x):
        return is_Algebra(x) and x.base_field() == self.base_field()

    def base_field(self):
        """
        Return the base field over which the algebras of this category
        are all defined.
        """
        return self.base_ring()


#############################################################
# CommutativeAlgebra
#############################################################
class CommutativeAlgebras(Category_over_base_ring):
    """
    The category of commutative algebras over a given base ring.

    EXAMPLES:
        sage: M = CommutativeAlgebras(GF(19))
        sage: M
        Category of commutative algebras over Finite Field of size 19
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = CommutativeAlgebras(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return CommutativeAlgebras, (self.base(), )

#############################################################
# MonoidAlgebra
#############################################################
class MonoidAlgebras(Category_over_base_ring):
    """
    The category of all monoid algebras over a given base ring.
    EXAMPLES:
        sage: MonoidAlgebras(GF(2))
        Category of monoid algebras over Finite Field of size 2
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = MonoidAlgebras(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return MonoidAlgebras, (self.base(), )


#############################################################
# GroupAlgebra
#############################################################
class GroupAlgebras(Category_over_base_ring):
    """
    EXAMPLES:
        sage: GroupAlgebras(IntegerRing())
        Category of group algebras over Integer Ring
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = GroupAlgebras(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return GroupAlgebras, (self.base(), )



#############################################################
# MatrixAlgebra
#############################################################
class MatrixAlgebras(Category_over_base_ring):
    """
    The category of matrix algebras over a field.

    EXAMPLES:
        sage: MatrixAlgebras(RationalField())
        Category of matrix algebras over Rational Field
    """

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = MatrixAlgebras(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return MatrixAlgebras, (self.base(), )



#############################################################
# RingModule
#############################################################
class RingModules(Category_module):
    """
    The category of all modules over a base ring.

    EXAMPLES:
        sage: Modules(RationalField())
        Category of ring modules over Rational Field

        sage: Modules(Integers(9))
        Category of ring modules over Ring of integers modulo 9
    """


    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = RingModules(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return RingModules, (self.base(), )


# Synonym

Modules = RingModules



#############################################################
# FreeModule
#############################################################
class FreeModules(Category_module):
    """
    The category of free modules over a base ring.

    EXAMPLES:
        sage: FreeModules(IntegerRing())
        Category of free modules over Integer Ring
    """
    def __init__(self, R):
        if not is_Ring(R):
            raise TypeError, "R (=%s) must be a ring"%R
        Category_module.__init__(self, R)


    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = FreeModules(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return FreeModules, (self.base(), )

    def __call__(self, x):
        try:
            M = x.free_module()
            if M.base_ring() != self.base_ring():
                M = M.change_ring(self.base_ring())
        except (TypeError, AttributeError), msg:
            raise TypeError, "%s\nunable to coerce x (=%s) into %s"%(msg,x,self)
        return M

    def is_abelian(self):
        return self.base_ring().is_field()


#############################################################
# VectorSpace
#############################################################
class VectorSpaces(Category_module):
    """
    The category of vector spaces over a specified field,
    with an embedding in an ambient vector space.

    EXAMPLES:
        sage: VectorSpaces(RationalField())
        Category of vector spaces over Rational Field
    """
    def __init__(self, K):
        if not is_Field(K):
            raise TypeError, "K (=%s) must be a field"%K
        Category_module.__init__(self, K)

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = QQ^10      # vector space
            sage: loads(C.dumps()) == C
            True
        """
        return VectorSpaces, (self.base(), )

    def __call__(self, x):
        try:
            V = x.vector_space(self.base_field())
            if V.base_field() != self.base_field():
                V = V.change_ring(self.base_field())
        except (TypeError, AttributeError), msg:
            raise TypeError, "%s\nunable to coerce x (=%s) into %s"%(msg,x,self)
        return V

    def base_field(self):
        return self.base_ring()


#############################################################
# HeckeModules
#############################################################
class HeckeModules(Category_module):
    r"""
    The category of Hecke modules.

    A Hecke module is a module $M$ over the \emph{anemic} Hecke
    algebra, i.e., the Hecke algebra generated by Hecke operators
    $T_n$ with $n$ coprime to the level of $M$.  (Every Hecke module
    defines a level function, which is a positive integer.)  The
    reason we require that $M$ only be a module over the anemic Hecke
    algebra is that many natural maps, e.g., degeneracy maps,
    Atkin-Lehner operators, etc., are $\T$-module homomorphisms; but
    they are homomorphisms over the anemic Hecke algebra.

    EXAMPLES:
    We create the category of Hecke modules over $\Q$.
        sage: C = HeckeModules(RationalField()); C
        Category of Hecke modules over Rational Field

    Note that the base ring can be an arbitrary commutative ring.
        sage: HeckeModules(IntegerRing())
        Category of Hecke modules over Integer Ring
        sage: HeckeModules(FiniteField(5))
        Category of Hecke modules over Finite Field of size 5

    The base ring doesn't have to be a principal ideal domain.
        sage: HeckeModules(PolynomialRing(IntegerRing(), 'x'))
        Category of Hecke modules over Univariate Polynomial Ring in x over Integer Ring
    """
    def __init__(self, R):
        if not (is_Ring(R) and R.is_commutative()):
            raise TypeError, "R (=%s) must be a commutative ring"%R
        Category_module.__init__(self, R, "Hecke modules")

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = HeckeModules(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return HeckeModules, (self.base(), )


#############################################################
# RingIdeal
#############################################################
class RingIdeals(Category_ideal):
    """
    The category of all ideals in a fixed ring.

    EXAMPLES:
        sage: Ideals(Integers(200))
        Category of ring ideals in Ring of integers modulo 200
        sage: C = Ideals(IntegerRing()); C
        Category of ring ideals in Integer Ring
        sage: I = C([8,12,18])
        sage: I
        Principal ideal (2) of Integer Ring
    """
    def __init__(self, R):
        if not is_Ring(R):
            raise TypeError, "R (=%s) must be a ring"%R
        Category_ideal.__init__(self, R)

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = RingIdeals(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return RingIdeals, (self.ring(), )


Ideals = RingIdeals


#############################################################
# CommutativeRingIdeal
#############################################################
class CommutativeRingIdeals(Category_ideal):
    """
    The category of ideals in a fixed commutative ring.

    EXAMPLES:

        sage: C = CommutativeRingIdeals(IntegerRing())
        sage: C
        Category of commutative ring ideals in Integer Ring
    """
    def __init__(self, R):
        import sage.rings.all
        if not sage.rings.all.is_CommutativeRing(R):
            raise TypeError, "R (=%s) must be a commutative ring"%R
        Category_ideal.__init__(self, R)


    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = CommutativeRingIdeals(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return CommutativeRingIdeals, (self.ring(), )



#############################################################
# AlgebraModule
#############################################################
class AlgebraModules(Category_module):
    """
    The category of modules over a fixed algebra $A$.
    """
    def __init__(self, A):
        if not is_Algebra(A):
            raise TypeError, "A (=%s) must be an algebra"%A
        Category_module.__init__(self, A)

    def algebra(self):
        return self.base_ring()

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = AlgebraModules(FreeAlgebra(QQ,2,'a,b'))
            sage: loads(C.dumps()) == C
            True
        """
        return AlgebraModules, (self.algebra(),)



#############################################################
# AlgebraIdeal
#############################################################
class AlgebraIdeals(Category_ideal):
    """
    The category of ideals in a fixed algebra $A$.

    EXAMPLES:
        sage: C = AlgebraIdeals(FreeAlgebra(QQ,2,'a,b'))
        sage: loads(C.dumps()) == C
        True
    """
    def __init__(self, A):
        if not is_Algebra(A):
            raise TypeError, "A (=%s) must be an algebra"%A
        Category_ideal.__init__(self, A)

    def algebra(self):
        return self.ambient()

    def __reduce__(self):
        return AlgebraIdeals, (self.algebra(),)



#############################################################
# CommutativeAlgebraIdeal
#############################################################
class CommutativeAlgebraIdeals(Category_ideal):
    """
    The category of ideals in a fixed commutative algebra $A$.
    """
    def __init__(self, A):
        if not is_CommutativeAlgebra(A):
            raise TypeError, "A (=%s) must be a commutative algebra"%A
        Category_in_ambient.__init__(self, A)

    def algebra(self):
        return self.ambient()


    def __reduce__(self):
        """
        """
        return CommutativeAlgebraIdeals, (self.algebra(),)

#############################################################
# ChainComplex
#############################################################
class ChainComplexes(Category_module):
    """
    The category of all chain complexes over a base ring.

    EXAMPLES::

        sage: ChainComplexes(RationalField())
        Category of chain complexes over Rational Field

        sage: ChainComplexes(Integers(9))
        Category of chain complexes over Ring of integers modulo 9
    """

    def __reduce__(self):
        """
        EXAMPLES::

            sage: C = ChainComplexes(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return ChainComplexes, (self.base(), )

#############################################################
# Schemes over a given base scheme.
#############################################################
class Schemes_over_base(Category_over_base):
    """
    The category of schemes over a given base scheme.

    EXAMPLES:
        sage: Schemes(Spec(ZZ))
        Category of schemes over Spectrum of Integer Ring
    """
    def __init__(self, Y):
        Category_over_base.__init__(self, Y)

    def base_scheme(self):
        return self.base()

    def _repr_(self):
        return "Category of schemes over %s"%self.base_scheme()

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = Schemes(ZZ)
            sage: loads(C.dumps()) == C
            True
        """
        return Schemes_over_base, (self.base_scheme(), )

class Schemes_abstract(Category_uniq):
    """
    The category of all abstract schemes.

    EXAMPLES:
        sage: Schemes()
        Category of Schemes
    """
    def __init__(self):
        Category_uniq.__init__(self, "Schemes")

    def __reduce__(self):
        return Schemes_abstract, ()

    def __call__(self, x):
        """
        EXAMPLES:
        Fetch the category of schemes.
            sage: S = Schemes(); S
            Category of Schemes

        We create a scheme from a ring.
            sage: X = S(ZZ); X
            Spectrum of Integer Ring

        We create a scheme from a scheme (do nothing).
            sage: S(X)
            Spectrum of Integer Ring

        We create a scheme morphism from a ring homomorphism.x
            sage: phi = ZZ.hom(QQ); phi
            Ring Coercion morphism:
              From: Integer Ring
              To:   Rational Field
            sage: f = S(phi); f
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

            sage: f.domain()
            Spectrum of Rational Field
            sage: f.codomain()
            Spectrum of Integer Ring
            sage: S(f)
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

        """
        from sage.rings.all import is_CommutativeRing, is_RingHomomorphism
        from sage.schemes.all import is_Scheme, Spec, is_SchemeMorphism
        if is_Scheme(x) or is_SchemeMorphism(x):
            return x
        elif is_CommutativeRing(x):
            return Spec(x)
        elif is_RingHomomorphism(x):
            A = Spec(x.codomain())
            return A.hom(x)
        else:
            raise TypeError, "No way to create an object or morphism in %s from %s"%(self, x)


def Schemes(X=None):
    if X is None:
        return Schemes_abstract()
    from sage.schemes.all import is_Scheme
    if not is_Scheme(X):
        X = Schemes()(X)
    return Schemes_over_base(X)


class ModularAbelianVarieties(Category_over_base):
    """
    The category of modular abelian varieties over a given field.

    EXAMPLES:
        sage: ModularAbelianVarieties(QQ)
        Category of modular abelian varieties over Rational Field
    """
    def __init__(self, Y):
        assert Y.is_field()
        Category_over_base.__init__(self, Y)

    def base_field(self):
        return self.base()

    def _repr_(self):
        return "Category of modular abelian varieties over %s"%self.base_field()

    def __reduce__(self):
        """
        EXAMPLES:
            sage: C = ModularAbelianVarieties(QQ)
            sage: loads(C.dumps()) == C
            True
        """
        return ModularAbelianVarieties, (self.base_field(), )

###################################################################
#
# Natural inclusions between categories.
#
###################################################################

category_hierarchy = {\
    Objects                : [Sets],\
    Sets                   : [],\
    GSets                  : [Sets],\
    Semigroups             : [Sets],\
    Monoids                : [Semigroups, Sets],\
    Groups                 : [Monoids, Semigroups, Sets],\
    AbelianSemigroups      : [Semigroups, Sets],\
    AbelianMonoids         : [Monoids, Sets],\
    AbelianGroups          : [Groups, AbelianMonoids, AbelianSemigroups, Sets],\
    Rings                  : [AbelianGroups, Sets],\
    CommutativeRings       : [Rings, Sets],\
    Fields                 : [CommutativeRings, Rings, Sets],\
    FiniteFields           : [Fields, CommutativeRings, Rings, Sets],\
    NumberFields           : [CommutativeRings, Rings, Sets, Fields],\
    Algebras               : [Rings, RingModules, Sets],\
    CommutativeAlgebras    : [Algebras, CommutativeRings, Rings, Sets],\
    MonoidAlgebras         : [Algebras, Sets],\
    GroupAlgebras          : [MonoidAlgebras, Algebras, Sets],\
    MatrixAlgebras         : [Algebras, Sets],\
    RingModules            : [AbelianGroups, Sets],\
    FreeModules            : [RingModules, AbelianGroups, Sets],\
    VectorSpaces           : [FreeModules, RingModules, AbelianGroups, Sets],\
    HeckeModules           : [FreeModules, RingModules, AbelianGroups, Sets],\
    AlgebraModules         : [RingModules, AbelianGroups, Sets],\
    RingIdeals             : [RingModules, AbelianGroups, Sets],\
    CommutativeRingIdeals  : [RingIdeals, AbelianGroups, Sets],\
    AlgebraIdeals          : [AlgebraModules, AbelianGroups, Sets],\
    CommutativeAlgebraIdeals:[CommutativeRingIdeals, RingIdeals, AbelianGroups, Sets], \
    Elements               : [Sets], \
    Sequences              : [Sets], \
    Schemes_abstract       : [Sets], \
    Schemes_over_base      : [Schemes_abstract, Sets], \
    ModularAbelianVarieties: [Sets], \
    }

for k in category_hierarchy:
    category_hierarchy[k].append(Objects)


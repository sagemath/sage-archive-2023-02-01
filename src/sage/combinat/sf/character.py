"""
character bases of the symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2015 Mike Zabrocki <zabrocki@mathstat.yorku.ca
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.sf.sfa import SymmetricFunctionAlgebra_generic as SFA_generic
from sage.misc.cachefunc import cached_method
from sage.categories.homset import Hom
from sage.combinat.partition import Partition
from sage.rings.arith import divisors, moebius
from sage.functions.other import binomial

class generic_character(SFA_generic):
    def _my_key(self, la):
        r"""
        A rank function for partitions

        The leading term of a homogeneous expression will
        be the partition with the largest key value.

        This key value is `|\lambda|^2+\lambda_0` and
        using the ``max`` function on a list of Partitions.

        Of course it is possible that this rank function
        is equal for some partitions, but the leading
        term should be any one of these partitions.

        INPUT:

        - ``la`` -- a partition

        OUTPUT:

        - an integer

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: ht._my_key(Partition([2,1,1]))
            18
            sage: ht._my_key(Partition([2,2]))
            18
            sage: ht._my_key(Partition([3,1]))
            19
            sage: ht._my_key(Partition([1,1,1,1]))
            17
        """
        if la:
            return la.size()**2+la[0]
        else:
            return 0

    def _other_to_self(self, sexpr):
        r"""
        Convert an expression the target basis to the character-basis

        We use triangularity to determine the expansion
        by subtracting off the leading term.  The target basis
        is specified by the method ``self.target_basis()``.

        INPUT:

        - ``sexpr`` -- an element of ``self.target_basis()`` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: h = Sym.h()
            sage: ht._other_to_self(h[2] + h([]))
            ht[] + ht[1] + ht[2]
            sage: st = SymmetricFunctions(QQ).st()
            sage: s = Sym.s()
            sage: st._other_to_self(s[1] + s([]))
            2*st[] + st[1]
        """
        if sexpr==0:
            return self(0)
        if sexpr.support()==[[]]:
            return self([])
        out = self.zero()
        while sexpr:
            mup = max(sexpr.support(),key=self._my_key)
            out += sexpr.coefficient(mup)*self(mup)
            sexpr -= sexpr.coefficient(mup)*self._self_to_other_on_basis(mup)
        return out

    class Element(SFA_generic.Element):
        pass

class character_basis(generic_character):
    r"""
    General code for a character basis (irreducible and induced trivial)

    This is a basis of the symmetric functions that has the
    property that ``self(la).character_to_frobenius_image(n)``
    is equal to ``other([n-sum(la)]+la)``
    It should also have the property that the (outer) structure
    constants are the analogue of the stable kronecker
    coefficients on the ``other`` basis (where ``other`` is either the
    Schur or homogeneous bases).

    REFERENCES:

    .. [OZ2015] R. Orellana, M. Zabrocki, *Symmetric group characters
       as symmetric functions*, :arxiv:`1000.0000v1`.

    EXAMPLES::

        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.s()
        sage: h = Sym.h()
        sage: ht = SymmetricFunctions(QQ).ht()
        sage: st = SymmetricFunctions(QQ).st()
        sage: ht(s[2,1])
        ht[1, 1] + ht[2, 1] - ht[3]
        sage: s(ht[2,1])
        s[1] - 2*s[1, 1] - 2*s[2] + s[2, 1] + s[3]
        sage: ht(h[2,1])
        ht[1] + 2*ht[1, 1] + ht[2, 1]
        sage: h(ht[2,1])
        h[1] - 2*h[1, 1] + h[2, 1]
        sage: st(ht[2,1])
        st[] + 2*st[1] + st[1, 1] + 2*st[2] + st[2, 1] + st[3]
        sage: ht(st[2,1])
        ht[1] - ht[1, 1] + ht[2, 1] - ht[3]
        sage: ht[2]*ht[1,1]
        ht[1, 1] + 2*ht[1, 1, 1] + ht[2, 1, 1]
        sage: h[4,2].kronecker_product(h[4,1,1])
        h[2, 2, 1, 1] + 2*h[3, 1, 1, 1] + h[4, 1, 1]
        sage: s(st[2,1])
        3*s[1] - 2*s[1, 1] - 2*s[2] + s[2, 1]
        sage: st(s[2,1])
        st[] + 3*st[1] + 2*st[1, 1] + 2*st[2] + st[2, 1]
        sage: st[2]*st[1]
        st[1] + st[1, 1] + st[2] + st[2, 1] + st[3]
        sage: s[4,2].kronecker_product(s[5,1])
        s[3, 2, 1] + s[3, 3] + s[4, 1, 1] + s[4, 2] + s[5, 1]
    """
    def __init__(self, Sym, other_basis, bname, pfix):
        r"""
        Initialize the basis and register coercions

        The coercions are set up between the ``other_basis``

        INPUT:

        - ``Sym`` -- an instance of the symmetric function algebra
        - ``other_basis`` -- a basis of Sym
        - ``bname`` -- the name for this basis (convention: ends in "character")
        - ``pfix`` -- a prefix to use for the basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht(); ht
            Symmetric Functions over Rational Field in the induced trivial character basis
            sage: st = SymmetricFunctions(QQ).st(); st
            Symmetric Functions over Rational Field in the irreducible symmetric group character basis
        """
        SFA_generic.__init__(self, Sym, basis_name=bname, prefix=pfix)
        self._other = other_basis
        from sage.categories.graded_hopf_algebras_with_basis \
          import GradedHopfAlgebrasWithBasis
        categ = GradedHopfAlgebrasWithBasis(self.base_ring())
        self.module_morphism(self._self_to_other_on_basis,
          codomain=self._other, category=categ).register_as_coercion()
        from sage.categories.morphism import SetMorphism
        self.register_coercion(
          SetMorphism(Hom(self._other, self, categ), self._other_to_self))

    def target_basis(self):
        r"""
        Return the basis that was passed as the ``other_basis``
        in the intialization.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: st = SymmetricFunctions(QQ).st()
            sage: ht.target_basis()
            Symmetric Functions over Rational Field in the homogeneous basis
            sage: st.target_basis()
            Symmetric Functions over Rational Field in the Schur basis
        """
        return self._other

    @cached_method
    def _self_to_other_on_basis(self, lam):
        r"""
        Convert a character-basis element to the ``self.target_basis()`` basis

        This is a recursive procedure that is calculated
        by the assumption that the leading term of ``self(lam)``
        is ``other(lam)`` and ``evalsf(self(lam),n) == other([n-sum(lam)]+lam)``

        INPUT:

        - ``lam`` -- a partition

        OUTPUT:

        - an expression in the ``self.target_basis()`` basis

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: ht = SymmetricFunctions(QQ).ht()
            sage: ht._self_to_other_on_basis(Partition([2,1]))
            h[1] - 2*h[1, 1] + h[2, 1]
            sage: st = SymmetricFunctions(QQ).st()
            sage: st._self_to_other_on_basis(Partition([2,1]))
            3*s[1] - 2*s[1, 1] - 2*s[2] + s[2, 1]
        """
        if lam==[]:
            return self._other([])
        n = sum(lam)+lam[0]
        sim = self._other(self._other(lam).character_to_frobenius_image(n))
        return self._other(lam)-sum(c*self._self_to_other_on_basis( \
          Partition(mu[1:])) for (mu,c) in sim if mu[1:]!=lam)

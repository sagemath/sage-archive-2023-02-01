"""
Generic dual bases symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
import sfa, multiplicative, classical
import sage.combinat.partition
from sage.matrix.all import matrix

class SymmetricFunctionAlgebra_dual(classical.SymmetricFunctionAlgebra_classical):
    def __init__(self, dual_basis, scalar, scalar_name="", prefix=None):
        """
        TESTS::

            sage: e = SFAElementary(QQ)
            sage: f = e.dual_basis()
            sage: f == loads(dumps(f))
            True
        """
        self._dual_basis = dual_basis
        self._scalar = scalar
        self._scalar_name = scalar_name


        #Set up the cache
        self._to_self_cache = {}
        self._from_self_cache = {}
        self._transition_matrices = {}
        self._inverse_transition_matrices = {}


        #
        scalar_target = scalar(sage.combinat.partition.Partition_class([1])).parent()
        scalar_target = (scalar_target(1)*dual_basis.base_ring()(1)).parent()


        self._p = sfa.SFAPower(scalar_target)

        if prefix is None:
            prefix = 'd_'+dual_basis.prefix()

        classical.SymmetricFunctionAlgebra_classical.__init__(self, scalar_target, "dual_"+dual_basis.basis_name(), SymmetricFunctionAlgebraElement_dual, prefix)


    def dual_basis(self, scalar=None, scalar_name="", prefix=None):
        """
        Return the dual basis to self. If a the scalar option is not
        passed, then it returns the dual basis with respect to the scalar
        product used to define self.

        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: h.dual_basis()
            Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis
            sage: m2 = h.dual_basis(zee, prefix='m2')
            sage: m([2])^2
            2*m[2, 2] + m[4]
            sage: m2([2])^2
            2*m2[2, 2] + m2[4]
        """
        if scalar is None:
            return self._dual_basis
        else:
            return SymmetricFunctionAlgebra_dual(self, scalar, scalar_name, prefix)

    def __repr__(self):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee); h #indirect doctests
            Dual basis to Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis
            sage: h = m.dual_basis(scalar=zee, scalar_name='Hall scalar product'); h #indirect doctest
            Dual basis to Symmetric Function Algebra over Rational Field, Monomial symmetric functions as basis with respect to the Hall scalar product
        """
        if self._scalar_name:
            return "Dual basis to %s"%self._dual_basis + " with respect to the " + self._scalar_name
        else:
            return "Dual basis to %s"%self._dual_basis


    def __call__(self, x):
        """
        Coerce things into self. We do this by coercing things into self's
        dual basis, and then converting that to self.

        EXAMPLES::

            sage: s = SFASchur(QQ['t'].fraction_field())
            sage: t = QQ['t'].fraction_field().gen()
            sage: zee_hl = lambda x: x.centralizer_size(t=t)
            sage: S = s.dual_basis(zee_hl)
            sage: S(s([2,1]))
            (-t/(t^5-2*t^4+t^3-t^2+2*t-1))*d_s[1, 1, 1] + ((-t^2-1)/(t^5-2*t^4+t^3-t^2+2*t-1))*d_s[2, 1] + (-t/(t^5-2*t^4+t^3-t^2+2*t-1))*d_s[3]
        """

        if sfa.is_SymmetricFunction(x):
            if x.parent() is self:
                return x
            else:
                #Coerce x into the dual basis
                d_x = self._dual_basis(x)
                return self._element_class(self, dual=d_x)


        return classical.SymmetricFunctionAlgebra_classical.__call__(self, x)


    def _precompute(self, n):
        """
        Computes the transition matrix between self and its dual basis for
        the homogenous component of size n.

        EXAMPLES::

            sage: e = SFAElementary(QQ)
            sage: f = e.dual_basis()
            sage: f._precompute(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(f._to_self_cache)
            [([1, 1], [([1, 1], 2), ([2], 1)]), ([2], [([1, 1], 1), ([2], 1)])]
            sage: l(f._from_self_cache)
            [([1, 1], [([1, 1], 1), ([2], -1)]), ([2], [([1, 1], -1), ([2], 2)])]
            sage: f._transition_matrices[2]
            [1 1]
            [1 2]
            sage: f._inverse_transition_matrices[2]
            [ 2 -1]
            [-1  1]
        """
        base_ring = self.base_ring()
        zero = base_ring(0)

        #Handle the n == 0 and n == 1 cases separately
        if n == 0 or n == 1:
            part = sage.combinat.partition.Partition_class([1]*n)
            self._to_self_cache[ part ] = { part: base_ring(1) }
            self._from_self_cache[ part ] = { part: base_ring(1) }
            self._transition_matrices[n] = matrix(base_ring, [[1]])
            self._inverse_transition_matrices[n] = matrix(base_ring, [[1]])
            return

        #Get all the basis elements of the n^th homogeneous component
        #of the dual basis and express them in the power-sum basis
        partitions_n = sage.combinat.partition.Partitions_n(n).list()
        d = {}
        for part in partitions_n:
            d[part] = self._p(self._dual_basis(part))._monomial_coefficients

        #This contains the data for the transition matrix from the
        #dual basis to self.
        transition_matrix_n = matrix(base_ring, len(partitions_n), len(partitions_n))

        #This first section calculates how the basis elements of the
        #dual basis are expressed in terms of self's basis.

        #For every partition p of size n, compute self(p) in
        #terms of the dual basis using the scalar product.
        i = 0
        for s_part in partitions_n:
            #s_part corresponds to self(dual_basis(part))
            #s_mcs  corresponds to self(dual_basis(part))._monomial_coefficients
            s_mcs = {}

            #We need to compute the scalar product of d[s_part] and
            #all of the d[p_part]'s
            j = 0
            for p_part in partitions_n:
                #Compute the scalar product of d[s_part] and d[p_part]
                sp = zero
                for ds_part in d[s_part]:
                    if ds_part in d[p_part]:
                        sp += d[s_part][ds_part]*d[p_part][ds_part]*self._scalar(ds_part)
                if sp != zero:
                    s_mcs[p_part] = sp
                    transition_matrix_n[i,j] = sp

                j += 1

            self._to_self_cache[ s_part ] = s_mcs
            i += 1

        #Save the transition matrix
        self._transition_matrices[n] = transition_matrix_n

        #This second section calculates how the basis elements of
        #self in terms of the dual basis.  We do this by computing
        #the inverse of the matrix obtained above.
        inverse_transition = ~transition_matrix_n

        for i in range(len(partitions_n)):
            d_mcs = {}
            for j in range(len(partitions_n)):
                if inverse_transition[i,j] != zero:
                    d_mcs[ partitions_n[j] ] = inverse_transition[i,j]

            self._from_self_cache[ partitions_n[i] ] = d_mcs

        self._inverse_transition_matrices[n] = inverse_transition

    def transition_matrix(self, basis, n ):
        r"""
        Returns the transition matrix between the `n^{th}`
        homogeneous component of self and basis.

        EXAMPLES::

            sage: s = SFASchur(QQ)
            sage: e = SFAElementary(QQ)
            sage: f = e.dual_basis()
            sage: f.transition_matrix(s, 5)
            [ 1 -1  0  1  0 -1  1]
            [-2  1  1 -1 -1  1  0]
            [-2  2 -1 -1  1  0  0]
            [ 3 -1 -1  1  0  0  0]
            [ 3 -2  1  0  0  0  0]
            [-4  1  0  0  0  0  0]
            [ 1  0  0  0  0  0  0]
            sage: e.transition_matrix(s, 5).inverse().transpose()
            [ 1 -1  0  1  0 -1  1]
            [-2  1  1 -1 -1  1  0]
            [-2  2 -1 -1  1  0  0]
            [ 3 -1 -1  1  0  0  0]
            [ 3 -2  1  0  0  0  0]
            [-4  1  0  0  0  0  0]
            [ 1  0  0  0  0  0  0]
        """
        if n not in self._transition_matrices:
            self._precompute(n)

        if basis is self._dual_basis:
            return self._inverse_transition_matrices[n]
        else:
            return self._inverse_transition_matrices[n]*self._dual_basis.transition_matrix(basis, n)


    def _multiply(self, left, right):
        """
        Multiplication is done by performing the multiplication in self's
        dual basis and then converting back to self.

        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: a = h([2])
            sage: b = a*a; b # indirect doctest
            d_m[2, 2]
            sage: b.dual()
            6*m[1, 1, 1, 1] + 4*m[2, 1, 1] + 3*m[2, 2] + 2*m[3, 1] + m[4]
        """

        #Do the multiplication in the dual basis
        #and then convert back to self.
        eclass = left.__class__
        d_product = left.dual()*right.dual()

        return eclass(self, dual=d_product)

class SymmetricFunctionAlgebraElement_dual(classical.SymmetricFunctionAlgebraElement_classical):
    def __init__(self, A, dictionary=None, dual=None):
        """
        Create an element of a dual basis.

        INPUT: At least one of the following must be specified. The one (if
        any) which is not provided will be computed.


        -  ``dictionary`` - an internal dictionary for the
           monomials and coefficents of self

        -  ``dual`` - self as an element of the dual basis.


        TESTS::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee, prefix='h')
            sage: a = h([2])
            sage: ec = h._element_class
            sage: ec(h, dual=m([2]))
            -h[1, 1] + 2*h[2]
            sage: h(m([2]))
            -h[1, 1] + 2*h[2]
            sage: h([2])
            h[2]
            sage: h([2])._dual
            m[1, 1] + m[2]
            sage: m(h([2]))
            m[1, 1] + m[2]
        """
        if dictionary is None and dual is None:
            raise ValueError, "you must specify either x or dual"

        parent = A
        base_ring = parent.base_ring()
        zero = base_ring(0)

        if dual is None:
            #We need to compute the dual
            dual_dict = {}
            from_self_cache = parent._from_self_cache

            #Get the underlying dictionary for self
            s_mcs = dictionary

            #Make sure all the conversions from self to
            #to the dual basis have been precomputed
            for part in s_mcs:
                if part not in from_self_cache:
                    parent._precompute(sum(part))

            #Create the monomial coefficient dictionary from the
            #the monomial coefficient dictionary of dual
            for s_part in s_mcs:
                from_dictionary = from_self_cache[s_part]
                for part in from_dictionary:
                    dual_dict[ part ] = dual_dict.get(part, zero) + base_ring(s_mcs[s_part]*from_dictionary[part])

            dual = parent._dual_basis._from_dict(dual_dict)


        if dictionary is None:
            #We need to compute the monomial coefficients dictionary
            dictionary = {}
            to_self_cache = parent._to_self_cache

            #Get the underlying dictionary for the
            #dual
            d_mcs = dual._monomial_coefficients

            #Make sure all the conversions from the dual basis
            #to self have been precomputed
            for part in d_mcs:
                if part not in to_self_cache:
                    parent._precompute(sum(part))

            #Create the monomial coefficient dictionary from the
            #the monomial coefficient dictionary of dual
            for d_part in d_mcs:
                to_dictionary = to_self_cache[d_part]
                for part in to_dictionary:
                    dictionary[ part ] = dictionary.get(part, zero) + base_ring(d_mcs[d_part]*to_dictionary[part])


        #Initialize self
        self._dual = dual
        classical.SymmetricFunctionAlgebraElement_classical.__init__(self, A, dictionary)


    def dual(self):
        """
        Returns self in the dual basis.

        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: a = h([2,1])
            sage: a.dual()
            3*m[1, 1, 1] + 2*m[2, 1] + m[3]
        """
        return self._dual

    def omega(self):
        """
        Returns the image of self under the Frobenius / omega
        automorphism.

        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(zee)
            sage: hh = SFAHomogeneous(QQ)
            sage: hh([2,1]).omega()
            h[1, 1, 1] - h[2, 1]
            sage: h([2,1]).omega()
            d_m[1, 1, 1] - d_m[2, 1]
        """
        eclass = self.__class__
        return eclass(self.parent(), dual=self._dual.omega() )

    def scalar(self, x):
        """
        Returns the standard scalar product of self and x.

        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: a = h([2,1])
            sage: a.scalar(a)
            2
        """
        return self._dual.scalar(x)

    def scalar_hl(self, x):
        """
        Returns the Hall-Littlewood scalar product of self and x.

        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: a = h([2,1])
            sage: a.scalar_hl(a)
            (t + 2)/(-t^4 + 2*t^3 - 2*t + 1)
        """
        return self._dual.scalar_hl(x)



    def _add_(self, y):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(zee)
            sage: a = h([2,1])+h([3]); a # indirect doctest
            d_m[2, 1] + d_m[3]
            sage: a.dual()
            4*m[1, 1, 1] + 3*m[2, 1] + 2*m[3]
        """
        eclass = self.__class__
        return eclass(self.parent(), dual=(self.dual()+y.dual()))

    def _neg_(self):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(zee)
            sage: a = h([2,1])+h([3])
            sage: -h([2,1]) # indirect doctest
            -d_m[2, 1]
        """
        eclass = self.__class__
        return eclass(self.parent(), dual=self.dual()._neg_())

    def _sub_(self, y):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(zee)
            sage: a = h([2,1])+h([3])
            sage: h([2,1])-h([3]) # indirect doctest
            d_m[2, 1] - d_m[3]
        """
        eclass = self.__class__
        return eclass(self.parent(), dual=(self.dual()-y.dual()))

    def _div_(self, y):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(zee)
            sage: a = h([2,1])+h([3])
            sage: a/2 # indirect doctest
            1/2*d_m[2, 1] + 1/2*d_m[3]
        """
        return self*(~y)

    def __invert__(self):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(zee)
            sage: a = h(2); a
            2*d_m[]
            sage: ~a
            1/2*d_m[]
        """
        eclass = self.__class__
        return eclass(self.parent(), dual=~self.dual())

    def expand(self, n, alphabet='x'):
        """
        EXAMPLES::

            sage: m = SFAMonomial(QQ)
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(zee)
            sage: a = h([2,1])+h([3])
            sage: a.expand(2)
            2*x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + 2*x1^3
            sage: a.dual().expand(2)
            2*x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + 2*x1^3
            sage: a.expand(2,alphabet='y')
            2*y0^3 + 3*y0^2*y1 + 3*y0*y1^2 + 2*y1^3
            sage: a.expand(2,alphabet='x,y')
            2*x^3 + 3*x^2*y + 3*x*y^2 + 2*y^3
        """
        return self._dual.expand(n, alphabet)

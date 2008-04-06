from functor import Functor
from category_types import *

# TODO, think through the rankings, and override pushout where necessary.

class ConstructionFunctor(Functor):
    def __mul__(self, other):
        if not isinstance(self, ConstructionFunctor) and not isinstance(other, ConstructionFunctor):
            raise TypeError, "Non-constructive product"
        return CompositConstructionFunctor(other, self)

    def pushout(self, other):
        if self.rank > other.rank:
            return self * other
        else:
            return other * self

    def __cmp__(self, other):
        """
        Equality here means that they are mathematically equivalent, though they may have specific implementation data.
        See the \code{merge} function.
        """
        return cmp(type(self), type(other))

    def __str__(self):
        s = str(type(self))
        import re
        return re.sub("<.*'.*\.([^.]*)'>", "\\1", s)

    def __repr__(self):
        return str(self)

    def merge(self, other):
        if self == other:
            return self
        else:
            return None

    def commutes(self, other):
        return False

class CompositeConstructionFunctor(ConstructionFunctor):
    def __init__(self, first, second):
        Functor.__init__(self, first.domain(), second.codomain())
        self._first = first
        self._second = second

    def __call__(self, R):
        return self._second(self._first(R))

    def __cmp__(self, other):
        c = cmp(self._first, other._first)
        if c == 0:
            c = cmp(self._second, other._second)
        return c

    def __str__(self):
        return "%s(%s)" % (self._second, self._first)

class IdentityConstructionFunctor(ConstructionFunctor):
    def __init__(self):
        Functor.__init__(self, Sets(), Sets())
        self.rank = -100
    def __call__(self, R):
        return R
    def __mul__(self, other):
        if isinstance(self, IdentityConstructionFunctor):
            return other
        else:
            return self

class PolynomialFunctor(ConstructionFunctor):
    def __init__(self, var, multi_variate=False):
        Functor.__init__(self, Rings(), Rings())
        self.var = var
        self.multi_variate = multi_variate
        self.rank = 9
    def __call__(self, R):
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.polynomial.multi_polynomial_ring_generic import is_MPolynomialRing
        if self.multi_variate and (is_MPolynomialRing(R) or is_PolynomialRing(R)):
            return PolynomialRing(R.base_ring(), (list(R.variable_names()) + [self.var]))
        else:
            return PolynomialRing(R, self.var)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.var, other.var)
        return c
    def merge(self, other):
        if self == other:
            return PolynomialFunctor(self.var, (self.multi_variate or other.multi_variate))
        else:
            return None
#    def __str__(self):
#        return "Poly(%s)" % self.var

class MatrixFunctor(ConstructionFunctor):
    def __init__(self, nrows, ncols, is_sparse=False):
#        if nrows == ncols:
#            Functor.__init__(self, Rings(), RingModules()) # takes a basering
#        else:
#            Functor.__init__(self, Rings(), MatrixAlgebras()) # takes a basering
        Functor.__init__(self, Rings(), Rings())
        self.nrows = nrows
        self.ncols = ncols
        self.is_sparse = is_sparse
        self.rank = 10
    def __call__(self, R):
        from sage.matrix.matrix_space import MatrixSpace
        return MatrixSpace(R, self.nrows, self.ncols, sparse=self.is_sparse)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp((self.nrows, self.ncols), (other.nrows, other.ncols))
        return c
    def merge(self, other):
        if self != other:
            return None
        else:
            return MatrixFunctor(self.nrows, self.ncols, self.is_sparse and other.is_sparse)

class VectorFunctor(ConstructionFunctor):
    def __init__(self, n, is_sparse=False, inner_product_matrix=None):
#        if nrows == ncols:
#            Functor.__init__(self, Rings(), RingModules()) # takes a basering
#        else:
#            Functor.__init__(self, Rings(), MatrixAlgebras()) # takes a basering
        Functor.__init__(self, Rings(), Rings())
        self.n = n
        self.is_sparse = is_sparse
        self.inner_product_matrix = inner_product_matrix
        self.rank = 10 # ranking of functor, not rank of module
    def __call__(self, R):
        from sage.modules.free_module import FreeModule
        return FreeModule(R, self.n, sparse=self.is_sparse, inner_product_matrix=self.inner_product_matrix)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.n, other.n)
        return c
    def merge(self, other):
        if self != other:
            return None
        else:
            return VectorFunctor(self.n, self.is_sparse and other.is_sparse)

class SubspaceFunctor(ConstructionFunctor):
    def __init__(self, basis):
        self.basis = basis
        self.rank = 11 # ranking of functor, not rank of module
    def __call__(self, ambient):
        return ambient.span_of_basis(self.basis)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.basis, other.basis)
        return c
    def merge(self, other):
        if isinstance(other, SubspaceFunctor):
            return SubspaceFunctor(self.basis + other.basis) # TODO: remove linear dependancies
        else:
            return None


class FractionField(ConstructionFunctor):
    def __init__(self):
        Functor.__init__(self, Rings(), Fields())
        self.rank = 5
    def __call__(self, R):
        return R.fraction_field()

class LocalizationFunctor(ConstructionFunctor):
    def __init__(self, t):
        Functor.__init__(self, Rings(), Rings())
        self.t = t
        self.rank = 6
    def __call__(self, R):
        return R.localize(t)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.t, other.t)
        return c

class CompletionFunctor(ConstructionFunctor):
    def __init__(self, p, prec, extras=None):
        Functor.__init__(self, Rings(), Rings())
        self.p = p
        self.prec = prec
        self.extras = extras
        self.rank = 4
    def __call__(self, R):
        return R.completion(self.p, self.prec, self.extras)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.p, other.p)
        return c
    def merge(self, other):
        if self.p == other.p:
            if self.prec == other.prec:
                extras = self.extras.copy()
                extras.update(other.extras)
                return CompletionFunctor(self.p, self.prec, extras)
            elif self.prec < other.prec:
                return self
            else: # self.prec > other.prec
                return other
        else:
            return None

class QuotientFunctor(ConstructionFunctor):
    def __init__(self, I):
        Functor.__init__(self, Rings(), Rings()) # much more general...
        self.I = I
        self.rank = 7
    def __call__(self, R):
        I = self.I
        if I.ring() != R:
            I.base_extend(R)
        return R.quo(I)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.I, other.I)
        return c
    def merge(self, other):
        if self == other:
            return self
        try:
            gcd = self.I + other.I
        except (TypeError, NotImplementedError):
            return None
        if gcd.is_trivial() and not gcd.is_zero():
            # quotient by gcd would result in the trivial ring/group/...
            # Rather than create the zero ring, we claim they can't be merged
            # TODO: Perhaps this should be detected at a higher level...
            raise TypeError, "Trivial quotient intersection."
        return QuotientFunctor(gcd)

class AlgebraicExtensionFunctor(ConstructionFunctor):
    def __init__(self, poly, name, elt=None):
        Functor.__init__(self, Rings(), Rings())
        self.poly = poly
        self.name = name
        self.elt = elt
        self.rank = 3
    def __call__(self, R):
        return R.extension(self.poly, self.name)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.poly, other.poly)
        return c

class AlgebraicClosureFunctor(ConstructionFunctor):
    def __init__(self):
        Functor.__init__(self, Rings(), Rings())
        self.rank = 3
    def __call__(self, R):
        return R.algebraic_closure()
    def merge(self, other):
        # Algebraic Closure subsumes Algebraic Extension
        return self

def BlackBoxConstructionFunctor(ConstructionFunctor):
    def __init__(self, box):
        self.box = box
        self.rank = 100
    def __call__(self, R):
        return box(R)
    def __cmp__(self, other):
        return self.box == other.box


def pushout(R, S):
    """
    Given a pair of Objects R and S, try and construct a
    reasonable object $Y$ and return maps such that
    cannonically $R \leftarrow Y \rightarrow S$.

    ALGORITHM:
       This incorperates the idea of functors discussed SAGE Days 4.
       Every object $R$ can be viewed as an initial object and
       a series of functors (e.g. polynomial, quotient, extension,
       completion, vector/matrix, etc.) Call the series of
       increasingly-simple rings (with the associated functors)
       the "tower" of $R$. The \code{construction} method is used to
       create the tower.

       Given two objects $R$ and $S$, try and find a common initial
       object $Z$. If the towers of $R$ and $S$ meet, let $Z$ be their
       join. Otherwise, see if the top of one coerces naturally into
       the other.

       Now we have an initial object and two \emph{ordered} lists of
       functors to apply. We wish to merge these in an unambiguous order,
       popping elements off the top of one or the other tower as we
       apply them to $Z$.

       - If the functors are distinct types, there is an absolute ordering
           given by the rank attribute. Use this.
       - Otherwise:
          - If the tops are equal, we (try to) merge them.
          - If \emph{exactly} one occurs lower in the other tower
              we may unambiguously apply the other (hoping for a later merge).
          - If the tops commute, we can apply either first.
          - Otherwise fail due to ambiguity.

    EXAMPLES:
        Here our "towers" are $R = Complete_7(Frac(\Z)$ and $Frac(Poly_x(\Z))$, which give us $Frac(Poly_x(Complete_7(Frac(\Z)))$
            sage: from sage.categories.pushout import pushout
            sage: pushout(Qp(7), Frac(ZZ['x']))
            Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20

        Note we get the same thing with
            sage: pushout(Zp(7), Frac(QQ['x']))
            Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20
            sage: pushout(Zp(7)['x'], Frac(QQ['x']))
            Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20

        Note that polynomial variable ordering must be unambiguously determined.
            sage: pushout(ZZ['x,y,z'], QQ['w,z,t'])
            Traceback (most recent call last):
            ...
            TypeError: Ambiguous Base Extension
            sage: pushout(ZZ['x,y,z'], QQ['w,x,z,t'])
            Multivariate Polynomial Ring in w, x, y, z, t over Rational Field

        Some other examples
            sage: pushout(Zp(7)['y'], Frac(QQ['t'])['x,y,z'])
            Multivariate Polynomial Ring in x, y, z over Fraction Field of Univariate Polynomial Ring in t over 7-adic Field with capped relative precision 20
            sage: pushout(ZZ['x,y,z'], Frac(ZZ['x'])['y'])
            Multivariate Polynomial Ring in y, z over Fraction Field of Univariate Polynomial Ring in x over Integer Ring
            sage: pushout(MatrixSpace(RDF, 2, 2), Frac(ZZ['x']))
            Full MatrixSpace of 2 by 2 dense matrices over Fraction Field of Univariate Polynomial Ring in x over Real Double Field
            sage: pushout(ZZ, MatrixSpace(ZZ[['x']], 3, 3))
            Full MatrixSpace of 3 by 3 dense matrices over Power Series Ring in x over Integer Ring
            sage: pushout(QQ['x,y'], ZZ[['x']])
            Univariate Polynomial Ring in y over Power Series Ring in x over Rational Field
            sage: pushout(Frac(ZZ['x']), QQ[['x']])
            Laurent Series Ring in x over Rational Field

    AUTHORS:
       -- Robert Bradshaw
    """
    if R is S or R == S:
        return R

    if isinstance(R, type):
        R = type_to_parent(R)

    if isinstance(S, type):
        S = type_to_parent(S)

    R_tower = construction_tower(R)
    S_tower = construction_tower(S)
    Rs = [c[1] for c in R_tower]
    Ss = [c[1] for c in S_tower]

    if R in Ss:
        return S
    elif S in Rs:
        return R

#    print Rs
#    print Ss

    if R_tower[-1][1] in Ss:
      Rs, Ss = Ss, Rs
      R_tower, S_tower = S_tower, R_tower

    # look for join
    if Ss[-1] in Rs:
        if Rs[-1] == Ss[-1]:
            while Rs[-1] == Ss[-1]:
                Rs.pop()
                Z = Ss.pop()
        else:
            Rs = Rs[:Rs.index(Ss[-1])]
            Z = Ss.pop()

    # look for topmost coercion
    elif S.has_coerce_map_from(Rs[-1]):
        while not Ss[-1].has_coerce_map_from(Rs[-1]):
            Ss.pop()
        while len(Rs) > 0 and Ss[-1].has_coerce_map_from(Rs[-1]):
            Rs.pop()
        Z = Ss.pop()

    elif R.has_coerce_map_from(Ss[-1]):
        while not Rs[-1].has_coerce_map_from(Ss[-1]):
            Rs.pop()
        while len(Ss) > 0 and Rs[-1].has_coerce_map_from(Ss[-1]):
            Ss.pop()
        Z = Rs.pop()

    else:
        raise TypeError, "No common base"

    # Rc is a list of functors from Z to R and Sc is a list of functors from Z to S
    Rc = [c[0] for c in R_tower[1:len(Rs)+1]]
    Sc = [c[0] for c in S_tower[1:len(Ss)+1]]

    while len(Rc) > 0 or len(Sc) > 0:
        # print Z
        # if we are out of functors in either tower, there is no ambiguity
        if len(Sc) == 0:
            c = Rc.pop()
            Z = c(Z)
        elif len(Rc) == 0:
            c = Sc.pop()
            Z = c(Z)
        # if one of the functors has lower rank, do it first
        elif Rc[-1].rank < Sc[-1].rank:
            c = Rc.pop()
            Z = c(Z)
        elif Sc[-1].rank < Rc[-1].rank:
            c = Sc.pop()
            Z = c(Z)
        else:
            # the ranks are the same, so things are a bit subtler
            if Rc[-1] == Sc[-1]:
                # If they are indeed the same operation, we only do it once.
                # The \code{merge} function here takes into account non-mathematical
                # distinctions (e.g. single vs. multivariate polynomials)
                cR = Rc.pop()
                cS = Sc.pop()
                c = cR.merge(cS) or cS.merge(cR)
                if c:
                    Z = c(Z)
                else:
                    raise TypeError, "Incompatable Base Extension %r, %r (on %r, %r)" % (R, S, cR, cS)
            else:
                # Now we look ahead to see if either top functor is
                # applied later on in the other tower.
                # If this is the case for exactly one of them, we unambiguously
                # postpone that operation, but if both then we abort.
                if Rc[-1] in Sc:
                    if Sc[-1] in Rc:
                        raise TypeError, "Ambiguous Base Extension"
                    else:
                        c = Sc.pop()
                        Z = c(Z)
                elif Sc[-1] in Rc:
                    c = Rc.pop();
                    Z = c(Z)
                # If, perchance, the two functors commute, then we may do them in any order.
                elif Rc[-1].commutes(Sc[-1]):
                    c = Rc.pop()
                    Z = c(Z)
                    c = Sc.pop()
                    Z = c(Z)
                else:
                    # try and merge (default merge is failure for unequal functors)
                    cR = Rc.pop()
                    cS = Sc.pop()
                    c = cR.merge(cS) or cS.merge(cR)
                    if c is not None:
                        Z = c(Z)
                    else:
                        # Otherwise, we cannot proceed.
                        raise TypeError, "Ambiguous Base Extension"
    return Z



def pushout_lattice(R, S):
    """
    Given a pair of Objects R and S, try and construct a
    reasonable object $Y$ and return maps such that
    cannonically $R \leftarrow Y \rightarrow S$.

    ALGORITHM:
       This is based on the model that arose from much discussion at SAGE Days 4.
       Going up the tower of constructions of $R$ and $S$ (e.g. the reals
       come from the rationals come from the integers) try and find a
       common parent, and then try and fill in a lattice with these
       two towers as sides with the top as the common ancestor and
       the bottom will be the desired ring.

       See the code for a specific worked-out example.

    EXAMPLES:
        sage: from sage.categories.pushout import pushout_lattice
        sage: A, B = pushout_lattice(Qp(7), Frac(ZZ['x']))
        sage: A.codomain()
        Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20
        sage: A.codomain() is B.codomain()
        True
        sage: A, B = pushout_lattice(ZZ, MatrixSpace(ZZ[['x']], 3, 3))
        sage: B
        Identity endomorphism of Full MatrixSpace of 3 by 3 dense matrices over Power Series Ring in x over Integer Ring

    AUTHOR:
       -- Robert Bradshaw
    """
    R_tower = construction_tower(R)
    S_tower = construction_tower(S)
    Rs = [c[1] for c in R_tower]
    Ss = [c[1] for c in S_tower]

    # look for common ancestor
    start = None
    for Z in Rs:
        if Z in Ss:
            start = Z
    if start is None:
        # Should I test for a map between the tops of the towers?
        # Or, if they're both not ZZ, is it hopeless?
        return None

    # truncate at common ancestor
    R_tower = list(reversed(R_tower[:Rs.index(start)+1]))
    S_tower = list(reversed(S_tower[:Ss.index(start)+1]))
    Rs = [c[1] for c in R_tower] # the list of objects
    Ss = [c[1] for c in S_tower]
    Rc = [c[0] for c in R_tower] # the list of functors
    Sc = [c[0] for c in S_tower]

    # Here we try and construct a 2-dimensional lattice as follows.
    # Suppose our towers are Z -> Q -> Qp = R and Z -> Z[t] -> Frac(Z[t]) = S
    lattice = {}
    # First we fill in the sides
    #
    #         Z
    #       /   \
    #      Q    Z[t]
    #    /         \
    #   Qp       Frac(Z[t])
    #
    for i in range(len(Rs)):
        lattice[i,0] = Rs[i]
    for j in range(len(Ss)):
        lattice[0,j] = Ss[j]

    # Now we attempt to fill in the center, one (diagonal) row at a time,
    # one commuting square at a time.
    #
    #          Z
    #       /    \
    #      Q     Z[t]
    #    /   \  /    \
    #   Qp   Q[t]   Frac(Z[t])
    #    \   /
    #    Qp[t]
    #
    # There is always exactly one "correct" path/order in which to apply operations
    # from the top to the bottom. In our example, this is down the far left side.
    # We keep track of which that is by clearing out Rc and Sc as we go along.
    #
    # Note that when applying the functors in the correct order, base extension
    # is not needed (though it may occur in the resulting morphisms).
    #
    for i in range(len(Rc)-1):
        for j in range(len(Sc)-1):
            try:
                if lattice[i,j+1] == lattice[i+1,j]:
                    # In this case we have R <- S -> R
                    # We don't want to perform the operation twice
                    # and all subsequent squares will come from objects
                    # where the operation was already performed (either
                    # to the left or right)
                    Rc[i] = Sc[j] = None # IdentityConstructionFunctor()
                    lattice[i+1,j+1] = lattice[i,j+1]
                elif Rc[i] is None and Sc[j] is None:
                    lattice[i+1,j+1] = lattice[i,j+1]
                elif Rc[i] is None:
                    lattice[i+1,j+1] = Sc[j](lattice[i+1,j])
                elif Sc[j] is None:
                    lattice[i+1,j+1] = Rc[i](lattice[i,j+1])
                else:
                    # For now, we just look at the rank.
                    # TODO: be more sophisticated and query the functors themselves
                    if Rc[i].rank < Sc[j].rank:
                        lattice[i+1,j+1] = Sc[j](lattice[i+1,j])
                        Rc[i] = None # force us to use pre-applied Rc[i]
                    else:
                        lattice[i+1,j+1] = Rc[i](lattice[i,j+1])
                        Sc[j] = None # force us to use pre-applied Sc[i]
            except (AttributeError, NameError):
                print i, j
                pp(lattice)
                raise TypeError, "%s does not support %s" % (lattice[i,j], 'F')

    # If we are successful, we should have something that looks like this.
    #
    #          Z
    #       /    \
    #      Q     Z[t]
    #    /   \  /    \
    #   Qp   Q[t]   Frac(Z[t])
    #    \   /  \    /
    #    Qp[t]  Frac(Q[t])
    #      \      /
    #     Frac(Qp[t])
    #
    R_loc = len(Rs)-1
    S_loc = len(Ss)-1

    # Find the composition coercion morphisms along the bottom left...
    if S_loc > 0:
        R_map = lattice[R_loc,1].coerce_map_from(R)
        for i in range(1, S_loc):
            map = lattice[R_loc, i+1].coerce_map_from(lattice[R_loc, i]) # The functor used is implicit here, should it be?
            R_map = map * R_map
    else:
        R_map = R.coerce_map_from(R) # id

    # ... and bottom right
    if R_loc > 0:
        S_map = lattice[1, S_loc].coerce_map_from(S)
        for i in range(1, R_loc):
            map = lattice[i+1, S_loc].coerce_map_from(lattice[i, S_loc])
            S_map = map * S_map
    else:
        S_map = S.coerce_map_from(S) # id

    return R_map, S_map


def pp(lattice):
    """
    Used in debugging to print the current lattice.
    """
    for i in range(100):
        for j in range(100):
            try:
                R = lattice[i,j]
                print i, j, R
            except KeyError:
                break

def construction_tower(R):
    tower = [(None, R)]
    c = R.construction()
    while c is not None:
        f, R = c
        if not isinstance(f, ConstructionFunctor):
            f = BlackBoxConstructionFunctor(f)
        tower.append((f,R))
        c = R.construction()
    return tower



def type_to_parent(P):
    import sage.rings.all
    if P in [int, long]:
        return sage.rings.all.ZZ
    elif P is float:
        return sage.rings.all.RDF
    elif P is complex:
        return sage.rings.all.CDF
    else:
        raise TypeError, "Not a scalar type."

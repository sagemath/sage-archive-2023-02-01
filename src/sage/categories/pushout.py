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

class CompositConstructionFunctor(ConstructionFunctor):
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
        from sage.rings.polynomial.polynomial_ring import PolynomialRing, is_PolynomialRing
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
    def __init__(self, nrows, ncols):
#        if nrows == ncols:
#            Functor.__init__(self, Rings(), RingModules()) # takes a basering
#        else:
#            Functor.__init__(self, Rings(), MatrixAlgebras()) # takes a basering
        Functor.__init__(self, Rings(), Rings())
        self.nrows = nrows
        self.ncols = ncols
        self.rank = 10
    def __call__(self, R):
        from sage.matrix.matrix_space import MatrixSpace
        return MatrixSpace(R, self.nrows, self.ncols)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp((self.nrows, self.ncols), (other.nrows, other.ncols))
        return c

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
        self.rank = 7
    def __call__(self, R):
        return R.completion(self.p, self.prec, self.extras)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.p, other.p)
        return c

class QuotientFunctor(ConstructionFunctor):
    def __init__(self, I):
        Functor.__init__(self, Rings(), Rings()) # much more general...
        self.I = I
        self.rank = 2
    def __call__(self, R):
        I = self.I
        if I.base_ring != R:
            I.base_extend(R)
        return R.quo(I)
    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.I, other.I)
        return c

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

def pushout(R, S):

    if R == S:
        return R

    R_tower = construction_tower(R)
    S_tower = construction_tower(S)
    Rs = [c[1] for c in R_tower]
    Ss = [c[1] for c in S_tower]

    if R in Ss:
        return S
    elif S in Rs:
        return R

    if R_tower[-1][1] in Ss:
      Rs, Ss = Ss, Rs
      R_tower, S_tower = S_tower, R_tower

#    print "Rs", Rs
#    print "Ss", Ss

    if Ss[-1] in Rs:
        if Rs[-1] == Ss[-1]:
            while Rs[-1] == Ss[-1]:
                Rs.pop()
                Z = Ss.pop()
        else:
            Rs = Rs[:Rs.index(Ss[-1])]
            Z = Ss.pop()
    elif Rs[-1].has_coerce_map_from(Ss[-1]):
        Z = Rs[-1]

    elif Ss[-1].has_coerce_map_from(Rs[-1]):
        Z = Ss[-1]

    else:
        raise TypeError, "No common base"

    # Rc is a list of functors from Z to R and Sc is a list of functors from Z to S
    Rc = [c[0] for c in R_tower[1:len(Rs)+1]]
    Sc = [c[0] for c in S_tower[1:len(Ss)+1]]

    while len(Rc) > 0 or len(Sc) > 0:
        print Z
        if len(Sc) == 0:
            c = Rc.pop()
            Z = c(Z)
        elif len(Rc) == 0:
            c = Sc.pop()
            Z = c(Z)
        elif Rc[-1].rank < Sc[-1].rank:
            c = Rc.pop()
            Z = c(Z)
        elif Sc[-1].rank < Rc[-1].rank:
            c = Sc.pop()
            Z = c(Z)
        else:
            if Rc[-1] == Sc[-1]:
                cR = Rc.pop()
                cS = Sc.pop()
                c = cR.merge(cS)
                if c:
                    Z = c(Z)
                else:
                    raise TypeError, "Incompatable Base Extension %r, %r (on %r, %r)" % (R, S, cR, cS)
            else:
                if Rc[-1] in Sc:
                    if Sc[-1] in Rc:
                        raise TypeError, "Ambiguous Base Extension"
                    else:
                        c = Sc.pop()
                        Z = c(Z)
                elif Sc[-1] in Rc:
                    c = Rc.pop();
                    Z = c(Z)
                elif Rc[-1].commutes(Sc[-1]):
                    c = Rc.pop()
                    Z = c(Z)
                    c = Sc.pop()
                    Z = c(Z)
                else:
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
        sage: from sage.categories.pushout import pushout
        sage: A, B = pushout(Qp(7), Frac(ZZ['x']))
        sage: A.codomain()
        Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20
        sage: A.codomain() is B.codomain()
        True
        sage: A, B = pushout(ZZ, MatrixSpace(ZZ[['x']], 3, 3))
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
        tower.append(c)
        R = c[1]
        c = R.construction()
    return tower

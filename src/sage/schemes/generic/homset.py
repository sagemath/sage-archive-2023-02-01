"""
Set of homomorphisms between two schemes
"""

import sage.rings.integer_ring
from sage.rings.arith import gcd
Z = sage.rings.integer_ring.ZZ

# Some naive point enumeration routines for default.
# AUTHOR: David R. Kohel <kohel@maths.usyd.edu.au>

def enum_projective_rational_field(X,B):
    n = X.codomain().ambient_space().ngens()
    Q = [ k+1 for k in range(B) ]
    R = [ 0 ] + [ s*k for k in Q for s in [1,-1] ]
    pts = []
    i = int(n-1)
    while not i < 0:
        P = [ 0 for _ in range(n) ]; P[i] = 1
        try:
            pts.append(X(P))
        except:
            pass
        iters = [ iter(R) for _ in range(i) ]
        [ iters[j].next() for j in range(i) ]
        j = 0
        while j < i:
            try:
                aj = ZZ(iters[j].next())
                P[j] = aj
                for ai in Q:
                    P[i] = ai
                    if gcd(P) == 1:
                        try:
                            pts.append(X(P))
                        except:
                            pass
                j = 0
            except StopIteration:
                iters[j] = iter(R)  # reset
                P[j] = 0
                P[j] = iters[j].next() # reset P[j] to 0 and increment
                j += 1
        i -= 1
    return pts

def enum_affine_rational_field(X,B):
    n = X.codomain().ambient_space().ngens()
    if X.value_ring() is ZZ:
        Q = [ 1 ]
    else: # rational field
        Q = [ k+1 for k in range(B) ]
    R = [ 0 ] + [ s*k for k in range(1,B+1) for s in [1,-1] ]
    pts = []
    P = [ 0 for _ in range(n) ]
    m = ZZ(0)
    try:
        pts.append(X(P))
    except:
        pass
    iters = [ iter(R) for _ in range(n) ]
    [ iters[j].next() for j in range(n) ]
    i = 0
    while i < n:
        try:
            a = ZZ(iters[i].next())
            m = m.gcd(a)
            P[i] = a
            for b in Q:
                if m.gcd(b) == 1:
                    try:
                        pts.append(X([ num/b for num in P ]))
                    except:
                        pass
            i = 0
            m = ZZ(0)
        except StopIteration:
            iters[i] = iter(R) # reset
	    P[i] = iters[i].next() # reset P[i] to 0 and increment
            i += 1
    return pts

def enum_projective_finite_field(X):
    n = X.codomain().ambient_space().ngens()
    R = X.value_ring()
    pts = []
    i = int(n-1)
    while not i < 0:
        P = [ 0 for _ in range(n) ]; P[i] = 1
        try:
            pts.append(X(P))
        except:
            pass
	# define some iterators and increment them:
        iters = [ iter(R) for _ in range(i) ]
        [ iters[j].next() for j in range(i) ]
        j = 0
        while j < i:
            try:
                P[j] = iters[j].next()
                try:
                    pts.append(X(P))
                except:
                    pass
                j = 0
            except StopIteration:
                iters[j] = iter(R) # reset iterator at j
                P[j] = iters[j].next() # reset P[j] to 0 and increment
                j += 1
        i -= 1
    return pts

def enum_affine_finite_field(X):
    n = X.codomain().ambient_space().ngens()
    R = X.value_ring()
    pts = []
    zero = R(0)
    P = [ zero for _ in range(n) ]
    pts.append(X(P))
    iters = [ iter(R) for _ in range(n) ]
    for x in iters: x.next() # put at zero
    i = 0
    while i < n:
        try:
            P[i] = iters[i].next()
            try:
                pts.append(X(P))
            except:
                pass
            i = 0
        except StopIteration:
            iters[i] = iter(R)  # reset
            iters[i].next() # put at zero
            P[i] = zero
            i += 1
    return pts

#*****************************************************************************
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.structure.parent_old as parent_old

from sage.categories.all import HomsetWithBase, Schemes
from sage.rings.all      import (
    is_FiniteField, is_RationalField, is_RingHomomorphism, ZZ)
import spec

import morphism

SCH = Schemes()

def is_SchemeHomset(H):
    return isinstance(H, SchemeHomset_generic)

def SchemeHomset(R, S, cat=None, check=True):
    if spec.is_Spec(R) and spec.is_Spec(S):
        return SchemeHomset_spec(R, S, cat=cat, check=check)
    else:
        return SchemeHomset_generic(R, S, cat=cat, check=check)

class SchemeHomset_generic(parent_old.Parent, HomsetWithBase):
    def __init__(self, X, Y, cat=None, check=True, base=ZZ):
        HomsetWithBase.__init__(self, X, Y, cat=cat, check=check, base=base)

    def has_coerce_map_from_impl(self, S):
        if self == S:   # an obvious case
            return True
        # Todo -- implement more cases.
        return False

    def _repr_(self):
        try:
            return "Set of points of %s defined over %s"%(self.codomain(), self.domain().coordinate_ring())
        except ValueError:
            return "Set of morphisms from %s to %s"%(self.domain(), self.codomain())

    def natural_map(self):
        X = self.domain()
        Y = self.codomain()
        if spec.is_Spec(Y) and Y.coordinate_ring() == X.base_ring():
            return morphism.SchemeMorphism_structure_map(self)
        raise NotImplementedError

    def __call__(self, x, check=True):
        """
        EXAMPLES::

            sage: f = ZZ.hom(QQ); f
            Ring Coercion morphism:
              From: Integer Ring
              To:   Rational Field

        ::

            sage: H = Hom(Spec(QQ,ZZ), Spec(ZZ)); H
            Set of points of Spectrum of Integer Ring defined over Rational Field

        ::

            sage: phi = H(f); phi
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field
        """
        if isinstance(x, (list, tuple)):
            return self.codomain()._point_morphism_class(self, x, check=check)

        if is_RingHomomorphism(x):
            return morphism.SchemeMorphism_spec(self, x, check=check)

class SchemeHomset_spec(SchemeHomset_generic):
    pass

class SchemeHomset_coordinates(SchemeHomset_generic):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __init__(self, X, S):
        R = X.base_ring()
        if R != S:
            X = X.base_extend(S)
        SchemeHomset_generic.__init__(self, spec.Spec(S, R), X)

    def _repr_(self):
        S = self.domain()
        if S == self.codomain().base_scheme():
            return "Set of Rational Points of %s"%self.codomain()
        if spec.is_Spec(S):
            S = S.coordinate_ring()
        return "Set of Rational Points over %s of %s"%(S, self.codomain())

    def value_ring(self):
        """
        Returns S for a homset X(T) where T = Spec(S).
        """
        T = self.domain()
        if spec.is_Spec(T):
            return T.coordinate_ring()
        else:
            raise TypeError, "Domain of argument must be of the form Spec(S)."

class SchemeHomset_affine_coordinates(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __call__(self, *v):
        if len(v) == 1:
            v = v[0]
        return morphism.SchemeMorphism_affine_coordinates(self, v)

    def points(self, B=0):
        try:
            R = self.value_ring()
        except TypeError:
            raise TypeError, "Domain of argument must be of the form Spec(S)."
        if is_RationalField(R) or R == Z:
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            return enum_affine_rational_field(self,B)
        elif is_FiniteField(R):
            return enum_affine_finite_field(self)
        else:
            raise TypeError, "Unable to enumerate points over %s."%R

class SchemeHomset_projective_coordinates_field(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __call__(self, *v):
        if len(v) == 1:
            v = v[0]
        X = self.codomain()
        try:
            return X._point_class(X, v)
        except AttributeError:  # should be very rare
            return morphism.SchemeMorphism_projective_coordinates_field(self, v)

    def points(self, B=0):
        try:
            R = self.value_ring()
        except TypeError:
            raise TypeError, "Domain of argument must be of the form Spec(S)."
        if is_RationalField(R):
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            return enum_projective_rational_field(self,B)
        elif is_FiniteField(R):
            return enum_projective_finite_field(self)
        else:
            raise TypeError, "Unable to enumerate points over %s."%R

class SchemeHomset_projective_coordinates_ring(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __call__(self, *v):
        raise NotImplementedError

    def points(self, B=0):
        raise NotImplementedError # fixed when __call__ is defined.
        try:
            R = self.value_ring()
        except TypeError:
            raise TypeError, "Domain of argument must be of the form Spec(S)."
        if R == Z:
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            return enum_projective_rational_field(self,B)
        else:
            raise TypeError, "Unable to enumerate points over %s."%R

class SchemeHomsetModule_abelian_variety_coordinates_field(SchemeHomset_projective_coordinates_field):
    def __init__(self, X, S, cat=None, check=True):
        r"""
        EXAMPLES: The bug reported at trac #1785 is fixed::

            sage: K.<a> = NumberField(x^2 + x - (3^3-3))
            sage: E = EllipticCurve('37a')
            sage: X = E(K)
            sage: X
            Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + (-1)*x over Number Field in a with defining polynomial x^2 + x - 24
            sage: P = X([3,a])
            sage: P
            (3 : a : 1)
            sage: P in E
            False
            sage: P in E.base_extend(K)
            True
        """
        R = X.base_ring()
        if R != S:
            X = X.base_extend(S)
        Y = spec.Spec(S, R)
        HomsetWithBase.__init__(self, Y, X, cat=cat,
                                check = check,
                                base = sage.rings.integer_ring.ZZ)

    def _repr_(self):
        return "Abelian group of points on %s"%self.codomain()

    def base_extend(self, R):
        if R != sage.rings.integer_ring.ZZ:
            raise NotImplementedError, "Abelian variety point sets not implemented as modules over rings other than ZZ."
        return self


r"""
Projective plane conics over `\QQ`

AUTHORS:

- Marco Streng (2010-07-20)

- Nick Alexander (2008-01-08)

"""
#*****************************************************************************
#       Copyright (C) 2008 Nick Alexander <ncalexander@gmail.com>
#       Copyright (C) 2009/2010 Marco Streng <marco.streng@gmail.com>
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

from sage.rings.all import (PolynomialRing, ZZ, QQ)

from sage.rings.morphism import is_RingHomomorphism
from sage.rings.real_mpfr import is_RealField

from sage.structure.sequence import Sequence
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.matrix.constructor import Matrix

from sage.quadratic_forms.qfsolve import qfsolve, qfparam

from con_number_field import ProjectiveConic_number_field

from sage.structure.element import is_InfinityElement

from sage.rings.arith import (lcm, hilbert_symbol)

class ProjectiveConic_rational_field(ProjectiveConic_number_field):
    r"""
    Create a projective plane conic curve over `\QQ`.
    See ``Conic`` for full documentation.

    EXAMPLES::

        sage: P.<X, Y, Z> = QQ[]
        sage: Conic(X^2 + Y^2 - 3*Z^2)
        Projective Conic Curve over Rational Field defined by X^2 + Y^2 - 3*Z^2

    TESTS::

        sage: Conic([2, 1, -1])._test_pickling()
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        EXAMPLES::

            sage: Conic([1, 1, 1])
            Projective Conic Curve over Rational Field defined by x^2 + y^2 + z^2
        """
        ProjectiveConic_number_field.__init__(self, A, f)


    def has_rational_point(self, point = False, obstruction = False,
                           algorithm = 'default', read_cache = True):
        r"""
        Returns True if and only if ``self`` has a point defined over `\QQ`.

        If ``point`` and ``obstruction`` are both False (default), then
        the output is a boolean ``out`` saying whether ``self`` has a
        rational point.

        If ``point`` or ``obstruction`` is True, then the output is
        a pair ``(out, S)``, where ``out`` is as above and the following
        holds:

         - if ``point`` is True and ``self`` has a rational point,
           then ``S`` is a rational point,

         - if ``obstruction`` is True and ``self`` has no rational point,
           then ``S`` is a prime such that no rational point exists
           over the completion at ``S`` or `-1` if no point exists over `\RR`.

        Points and obstructions are cached, whenever they are found.
        Cached information is used if and only if ``read_cache`` is True.

        ALGORITHM:

        The parameter ``algorithm``
        specifies the algorithm to be used:

         - ``'qfsolve'`` -- Use Denis Simon's GP script ``qfsolve``
           (see ``sage.quadratic_forms.qfsolve.qfsolve``)

         - ``'rnfisnorm'`` -- Use PARI's function rnfisnorm
           (cannot be combined with ``obstruction = True``)

         - ``'local'`` -- Check if a local solution exists for all primes
           and infinite places of `\QQ` and apply the Hasse principle
           (cannot be combined with ``point = True``)

         - ``'default'`` -- Use ``'qfsolve'``

         - ``'magma'`` (requires Magma to be installed) --
           delegates the task to the Magma computer algebra
           system.

        EXAMPLES::

            sage: C = Conic(QQ, [1, 2, -3])
            sage: C.has_rational_point(point = True)
            (True, (1 : 1 : 1))
            sage: D = Conic(QQ, [1, 3, -5])
            sage: D.has_rational_point(point = True)
            (False, 3)
            sage: P.<X,Y,Z> = QQ[]
            sage: E = Curve(X^2 + Y^2 + Z^2); E
            Projective Conic Curve over Rational Field defined by X^2 + Y^2 + Z^2
            sage: E.has_rational_point(obstruction = True)
            (False, -1)

        The following would not terminate quickly with
        ``algorithm = 'rnfisnorm'`` ::

            sage: C = Conic(QQ, [1, 113922743, -310146482690273725409])
            sage: C.has_rational_point(point = True)
            (True, (-76842858034579/5424 : -5316144401/5424 : 1))
            sage: C.has_rational_point(algorithm = 'local', read_cache = False)
            True
            sage: C.has_rational_point(point=True, algorithm='magma', read_cache=False) # optional - magma
            (True, (30106379962113/7913 : 12747947692/7913 : 1))

        TESTS:

        Create a bunch of conics over `\QQ`, check if ``has_rational_point`` runs without errors
        and returns consistent answers for all algorithms. Check if all points returned are valid. ::

            sage: l = Sequence(cartesian_product_iterator([[-1, 0, 1] for i in range(6)]))
            sage: c = [Conic(QQ, a) for a in l if a != [0,0,0] and a != (0,0,0,0,0,0)]
            sage: d = []
            sage: d = [[C]+[C.has_rational_point(algorithm = algorithm, read_cache = False, obstruction = (algorithm != 'rnfisnorm'), point = (algorithm != 'local')) for algorithm in ['local', 'qfsolve', 'rnfisnorm']] for C in c[::10]] # long time: 7 seconds
            sage: assert all([e[1][0] == e[2][0] and e[1][0] == e[3][0] for e in d])
            sage: assert all([e[0].defining_polynomial()(Sequence(e[i][1])) == 0 for e in d for i in [2,3] if e[1][0]])
        """
        if read_cache:
            if self._rational_point is not None:
                if point or obstruction:
                    return True, self._rational_point
                else:
                    return True
            if self._local_obstruction is not None:
                if point or obstruction:
                    return False, self._local_obstruction
                else:
                    return False
            if (not point) and self._finite_obstructions == [] and \
               self._infinite_obstructions == []:
                if obstruction:
                    return True, None
                return True
        if self.has_singular_point():
            if point:
                return self.has_singular_point(point = True)
            if obstruction:
                return True, None
            return True
        if algorithm == 'default' or algorithm == 'qfsolve':
            M = self.symmetric_matrix()
            M *= lcm([ t.denominator() for t in M.list() ])
            pt = qfsolve(M)
            if pt in ZZ:
                if self._local_obstruction is None:
                    self._local_obstruction = pt
                if point or obstruction:
                    return False, pt
                return False
            pt = self.point([pt[0], pt[1], pt[2]])
            if point or obstruction:
                return True, pt
            return True
        ret = ProjectiveConic_number_field.has_rational_point( \
                                           self, point = point, \
                                           obstruction = obstruction, \
                                           algorithm = algorithm, \
                                           read_cache = read_cache)
        if point or obstruction:
            if is_RingHomomorphism(ret[1]):
                ret[1] = -1
        return ret


    def is_locally_solvable(self, p):
        r"""
        Returns True if and only if ``self`` has a solution over the
        `p`-adic numbers. Here `p` is a prime number or equals
        `-1`, infinity, or `\RR` to denote the infinite place.

        EXAMPLES::

            sage: C = Conic(QQ, [1,2,3])
            sage: C.is_locally_solvable(-1)
            False
            sage: C.is_locally_solvable(2)
            False
            sage: C.is_locally_solvable(3)
            True
            sage: C.is_locally_solvable(QQ.hom(RR))
            False
            sage: D = Conic(QQ, [1, 2, -3])
            sage: D.is_locally_solvable(infinity)
            True
            sage: D.is_locally_solvable(RR)
            True

        """
        D, T = self.diagonal_matrix()
        abc = [D[j, j] for j in range(3)]
        if abc[2] == 0:
            return True
        a = -abc[0]/abc[2]
        b = -abc[1]/abc[2]
        if is_RealField(p) or is_InfinityElement(p):
            p = -1
        elif is_RingHomomorphism(p):
            if p.domain() is QQ and is_RealField(p.codomain()):
                p = -1
            else:
                raise TypeError("p (=%s) needs to be a prime of base field " \
                                 "B ( =`QQ`) in is_locally_solvable" % p)
        if hilbert_symbol(a, b, p) == -1:
            if self._local_obstruction is None:
                self._local_obstruction = p
            return False
        return True


    def local_obstructions(self, finite = True, infinite = True, read_cache = True):
        r"""
        Returns the sequence of finite primes and/or infinite places
        such that self is locally solvable at those primes and places.

        The infinite place is denoted `-1`.

        The parameters ``finite`` and ``infinite`` (both True by default) are
        used to specify whether to look at finite and/or infinite places.
        Note that ``finite = True`` involves factorization of the determinant
        of ``self``, hence may be slow.

        Local obstructions are cached. The parameter ``read_cache`` specifies
        whether to look at the cache before computing anything.

        EXAMPLES ::

            sage: Conic(QQ, [1, 1, 1]).local_obstructions()
            [2, -1]
            sage: Conic(QQ, [1, 2, -3]).local_obstructions()
            []
            sage: Conic(QQ, [1, 2, 3, 4, 5, 6]).local_obstructions()
            [41, -1]

        """
        obs0 = []
        obs1 = []
        if infinite:
            if read_cache and self._infinite_obstructions is not None:
                obs0 = self._infinite_obstructions
            else:
                if not self.is_locally_solvable(-1):
                    obs0 = [-1]
                self._infinite_obstructions = obs0
        if finite:
            if read_cache and self._finite_obstructions is not None:
                obs1 = self._finite_obstructions
            else:
                candidates = []
                if self.determinant() != 0:
                    for a in self.symmetric_matrix().list():
                        if a != 0:
                            for f in a.factor():
                                if f[1] < 0 and not f[0] in candidates:
                                    candidates.append(f[0])
                    for f in (2*self.determinant()).factor():
                        if f[1] > 0 and not f[0] in candidates:
                            candidates.append(f[0])
                for b in candidates:
                    if not self.is_locally_solvable(b):
                       obs1.append(b)
                self._infinite_obstructions = obs1
        obs = obs1 + obs0
        if finite and infinite:
            assert len(obs) % 2 == 0
        return obs


    def parametrization(self, point=None, morphism=True):
        r"""
        Return a parametrization `f` of ``self`` together with the
        inverse of `f`.

        If ``point`` is specified, then that point is used
        for the parametrization. Otherwise, use ``self.rational_point()``
        to find a point.

        If ``morphism`` is True, then `f` is returned in the form
        of a Scheme morphism. Otherwise, it is a tuple of polynomials
        that gives the parametrization.

        ALGORITHM:

        Uses Denis Simon's GP script ``qfparam``.
        See ``sage.quadratic_forms.qfsolve.qfparam``.

        EXAMPLES ::

            sage: c = Conic([1,1,-1])
            sage: c.parametrization()
            (Scheme morphism:
              From: Projective Space of dimension 1 over Rational Field
              To:   Projective Conic Curve over Rational Field defined by x^2 + y^2 - z^2
              Defn: Defined on coordinates by sending (x : y) to
                    (2*x*y : -x^2 + y^2 : x^2 + y^2),
             Scheme morphism:
              From: Projective Conic Curve over Rational Field defined by x^2 + y^2 - z^2
              To:   Projective Space of dimension 1 over Rational Field
              Defn: Defined on coordinates by sending (x : y : z) to
                    (1/2*x : 1/2*y + 1/2*z))

        An example with ``morphism = False`` ::

            sage: R.<x,y,z> = QQ[]
            sage: C = Curve(7*x^2 + 2*y*z + z^2)
            sage: (p, i) = C.parametrization(morphism = False); (p, i)
            ([-2*x*y, 7*x^2 + y^2, -2*y^2], [-1/2*x, -1/2*z])
            sage: C.defining_polynomial()(p)
            0
            sage: i[0](p) / i[1](p)
            x/y

        A ``ValueError`` is raised if ``self`` has no rational point ::

            sage: C = Conic(x^2 + 2*y^2 + z^2)
            sage: C.parametrization()
            Traceback (most recent call last):
            ...
            ValueError: Conic Projective Conic Curve over Rational Field defined by x^2 + 2*y^2 + z^2 has no rational points over Rational Field!

        A ``ValueError`` is raised if ``self`` is not smooth ::

            sage: C = Conic(x^2 + y^2)
            sage: C.parametrization()
            Traceback (most recent call last):
            ...
            ValueError: The conic self (=Projective Conic Curve over Rational Field defined by x^2 + y^2) is not smooth, hence does not have a parametrization.
        """
        if (not self._parametrization is None) and not point:
            par = self._parametrization
        else:
            if not self.is_smooth():
                raise ValueError("The conic self (=%s) is not smooth, hence does not have a parametrization." % self)
            if point is None:
                point = self.rational_point()
            point = Sequence(point)
            Q = PolynomialRing(QQ, 'x,y')
            [x, y] = Q.gens()
            gens = self.ambient_space().gens()
            M = self.symmetric_matrix()
            M *= lcm([ t.denominator() for t in M.list() ])
            par1 = qfparam(M, point)
            B = Matrix([[par1[i][j] for j in range(3)] for i in range(3)])
            # self is in the image of B and does not lie on a line,
            # hence B is invertible
            A = B.inverse()
            par2 = [sum([A[i,j]*gens[j] for j in range(3)]) for i in [1,0]]
            par = ([Q(pol(x/y)*y**2) for pol in par1], par2)
            if self._parametrization is None:
                self._parametrization = par
        if not morphism:
            return par
        P1 = ProjectiveSpace(self.base_ring(), 1, 'x,y')
        return P1.hom(par[0],self), self.Hom(P1)(par[1], check = False)


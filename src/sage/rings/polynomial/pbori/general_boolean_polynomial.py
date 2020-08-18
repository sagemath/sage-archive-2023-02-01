from __future__ import print_function

import sys
import resource
from .gbcore import *
from time import time, clock
from .PyPolyBoRi import *
from .blocks import declare_ring, Block
from math import sqrt
from .parallel import groebner_basis_first_finished
from optparse import OptionParser
from .interred import interred


# Just for debugging
def print_matrix(A):
    res = ""
    for i in range(len(A)):
        for j in range(len(A[i])):
            res += str(A[i][j]) + " "
        res += "\n"
    return res


# TODO: Implement constructor to convert Polynomials to
# GeneralBooleanPolynomial (with e_1 + ... + e_k as coefficient)
class GeneralBooleanPolynomial:
    r"""
    Class to represent Boolean polynomials over F_2^k
    """

    def __init__(self, k, coeff, polynomial):
        r"""
        Construct a GeneralBooleanPolynomial given by coeff * polynomial
        
        INPUT:
        
        - ``k`` -- Number of factors of F_2
        - ``coeff`` -- Array containing natural numbers in {0, ..., k-1} representing an element of F_2^k as a set
        - ``polynomial`` -- Polynomial
        """
        #print "type of polynomials", type(polynomial)
        #print "polynomials is int?", isinstance(polynomial, int)
        assert isinstance(polynomial, Polynomial) or isinstance(polynomial,
            Monomial) or isinstance(polynomial, Variable) or isinstance(
            polynomial, int)
        self.polys = []
        self.k = k
        self.ring = polynomial.ring()
        for i in range(k):
            if i in coeff:
                self.polys.append(polynomial)
            else:
                self.polys.append(Polynomial(0, self.ring))

    def __len__(self):
        r"""
        Returns the number of factors k in the product of the underlying ring F_2^k
        """
        return self.k

    def __getitem__(self, k):
        r"""
        Return the k-th component (i.e. the projection to the k-th factor)
        """
        return self.polys[k]

    def __setitem__(self, k, value):
        r"""
        Sets the k-th component (i.e. the projection to the k-th factor)
        """
        self.polys[k] = value

    def __eq__(self, other):
        r"""
        Tests equality by testing that
        - both objects are defined over the same ring (i.e. the number of factors is the same)
        - the objects are equal in each component
        """
        if not len(self) == len(other):
            return False
        for i in range(len(self)):
            if self[i] != other[i]:
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        r"""
        Returns a representation of the polynomial as string
        """
        res = ""
        # Build a list of all terms occurring
        terms = set([])
        for i in range(self.k):
            if not isinstance(self.polys[i], Polynomial):
                assert isinstance(self.polys[i], Monomial) or isinstance(self.
                    polys[i], Variable)
                self[i] = Polynomial(self[i])
            terms = terms.union(set(self[i].terms()))
        # Sort the terms
        terms = list(terms)
        terms.sort()
        # Determine the coefficient for each term and build up the string
        for t in terms:
            comps = [j for j in range(self.k) if t in set(Polynomial(self[j])
                .terms())]
            if len(comps) == self.k:
                # We use this simplified notation of the coefficient it 1
                res += str(t) + " + "
            else:
                res += str(comps) + " * " + str(t) + " + "
        res = res[0:len(res) - 3]
        if res == "":
            return "0"
        else:
            return res

    def __add__(self, other):
        r"""
        Addition of two GeneralBooleanPolynomial
        """
        if not len(self) == len(other):
            print("Cannot add polynomials defined over different rings")
            print("Len(self)=", len(self))
            print("Len(other)=", len(other))
            assert len(self) == len(other)
        sum = GeneralBooleanPolynomial(self.k, [], Polynomial(0, self.ring))
        for i in range(self.k):
            sum[i] = self[i] + other[i]
        return sum

    def __mul__(self, other):
        r"""
        Multiplication of two GeneralBooleanPolynomial
        """
        if not len(self) == len(other):
            print("Cannot multiply polynomials defined over different rings")
            print("Len(self)=", len(self))
            print("Len(other)=", len(other))
            assert len(self) == len(other)
        prod = GeneralBooleanPolynomial(self.k, [], Polynomial(0, self.ring))
        for i in range(self.k):
            prod[i] = Polynomial(self[i]) * Polynomial(other[i])
        return prod

    def __sub__(self, other):
        r"""
        Subtraction of two GeneralBooleanPolynomial
        """
        if not len(self) == len(other):
            print("Cannot subtract polynomials defined over different rings")
            print("Len(self)=", len(self))
            print("Len(other)=", len(other))
            assert len(self) == len(other)
        sub = GeneralBooleanPolynomial(self.k, [], Polynomial(1, self.ring))
        for i in range(self.k):
            sub[i] = self[i] - other[i]
        return sub

    def lc(self):
        r"""
        Returns leading coefficient as constant GeneralBooleanPolynomial
        """
        return GeneralBooleanPolynomial(self.k, self.lc_as_set(),
                        Polynomial(1, self.ring))

    def lc_as_set_array(self):
        r"""
        Returns leading coefficient as array containing the indices of the
        non-zero components of the leading coefficient.
        """
        max_term = 1
        for i in range(len(self)):
            if not self[i].is_zero():
                if self[i].lead() > max_term:
                    max_term = self[i].lead()
        comps = [j for j in range(len(self)) if max_term
                 in set(self[j].terms())]
        return comps

    def lc_as_set(self):
        r"""
        Returns leading coefficient as set containing the indices of the
        non-zero components of the leading coefficient.
        """
        return set(self.lc_as_set_array())

    def lc_binary(self):
        r"""
        Returns leading coefficient as array containing the integers 0 or 1
        representing the leading coefficient in a binary form.
        """
        lc_set = self.lc_as_set()
        lc_binary = [0] * len(self)
        for i in range(len(self)):
            if i in lc_set:
                lc_binary[i] = 1
        return lc_binary

    def lt(self):
        r"""
        Leading term in form of a GeneralBooleanPolynomial
        """
        max_term = 1
        for i in range(self.k):
            if not self[i].is_zero():
                if self[i].lead() > max_term:
                    max_term = self[i].lead()
        comps = [j for j in range(self.k) if max_term in set(self[j].terms())]
        return GeneralBooleanPolynomial(self.k, comps, max_term)

    def lm(self):
        r"""
        Leading monomial in form of a GeneralBooleanPolynomial
        """
        max_term = 1
        contains_non_zero_term = False
        for i in range(len(self)):
            if not self[i].is_zero():
                contains_non_zero_term = True
                if self[i].lead() > max_term:
                    max_term = self[i].lead()
        if not contains_non_zero_term:
            raise ValueError("lm of zero polynomial is not defined")
        return GeneralBooleanPolynomial(self.k, [i for i in range(self.k)],
            max_term)

    def constant_part_binary(self):
        r"""
        Constant part as binary tuple indicading which coefficients are non-zero
        """
        comps = []
        for i in range(len(self)):
            if self[i].has_constant_part():
                comps.append(1)
            else:
                comps.append(0)
        return comps

    def constant_part_as_set_array(self):
        r"""
        Constant part as array containing the indices of the non-zero coefficients of the constant part (sorted increasingly)
        """
        res = []
        for i in range(len(self)):
            if self[i].has_constant_part():
                res.append(i)
        return res

    def constant_part_as_set(self):
        r"""
        Constant part as set containing the indices of the non-zero coefficients of the constant part
        """
        return set(self.constant_part_as_set_array())

    def constant_part(self):
        r"""
        Constant part as GeneralBoolenPolynomial
        """
        res = GeneralBooleanPolynomial(len(self), [], Polynomial(0, self.ring))
        for i in range(len(self)):
            if self[i].has_constant_part():
                res[i] = Polynomial(1, self.ring)
            else:
                res[i] = Polynomial(0, self.ring)
        return res

    def to_expanded_polynomial_ring(self, new_variables):
        r"""
        Returns a representation in form of a Polynomial over a ring with
        additional variables, one for each factor in the product of fields
        F_2^k
        """
        assert self.k == len(new_variables)
        return sum(new_variables[i] * self.polys[i] for i in range(self.k))

    def is_monomial(self):
        r"""
        Test if self is a Monomial
        """
        # Build a list of all terms occurring
        terms = set([])
        for i in range(self.k):
            if not isinstance(self.polys[i], Polynomial):
                assert isinstance(self.polys[i], Monomial) or isinstance(self.
                    polys[i], Variable)
                self[i] = Polynomial(self[i])
            terms = terms.union(set(self[i].terms()))
        if len(terms) > 1:
            return False
        else:
            return True

    def is_zero(self):
        r"""
        Tests if self is zero
        """
        for i in range(len(self)):
            if self[i] != 0:
                return False
        return True

    def monomial(self):
        r"""
        Returns a PolyBoRi Monomial representing the leading monomial of self, where self should be a monomial
        """
        assert self.is_monomial()
        for i in range(self.k):
            if self.polys[i] != Polynomial(0, self.ring):
                return self.polys[i].lead()
        return Polynomial(0, self.ring)

    def divides(self, other):
        r"""
        Tests if self divides other
        """
        assert len(self) == len(other)
        assert (self.is_monomial() and other.is_monomial())
        self_divides_other = True
        for i in range(len(self)):
            if self[i] == 0:
                if other[i] != 0:
                    return False
                else:
                    continue
            if self[i] == 1:
                continue
            else:
                if other[i] == 1:
                    return False

            if other[i] == 0:
                continue

            assert other[i] == other[i].lead()
            other[i] = other[i].lead()
            if not other[i] in (Polynomial(self.polys[i]).lead()).divisors():
                return False
        return True


# Simple Gaus algorithm
# Should really be replaced by something faster (maybe from PolyBoRi)
def triangulate_over_F2(A, b):
    assert len(A) == len(b)
    n = len(b)
    m = len(A[0])
    print("n: ", n, "m: ", m)
    print("A, b", A, b)
    for i in range(0, min(n, m)):    # Row
        print("i", i)
        if A[i][i] == 0:
            # permutate rows
            changed = False
            for l in range(i, n):
                if A[l][i] != 0:
                    for k in range(n):
                        A[l][k], A[i][k] = A[i][k], A[l][k]
                    b[l], b[i] = b[i], b[l]
                changed = True
            if not changed:
                return -1
        for j in range(0, i):
            if A[i][j] != 0:
                for k in range(j, n):
                    A[i][k] += A[i - 1][k]
                b[i] += b[i - 1]
    res_A = [[A[i][j] % 2 for j in range(m)] for i in range(n)]
    res_b = [b[i] % 2 for i in range(n)]
    return (res_A, res_b)


def projection_of_expanded_polynomial(f, e, e_vars):
    r"""
    Compute the projection of the expanded polynomial f to the component
    corresponding to the variable e (which is part of e_vars)
    """
    for v in e_vars:
        if v == e:
            # Substitute v by 1
            f = Polynomial(f.set().subset0(v.index())) + Polynomial(f.set().
                subset1(v.index()))
        else:
            # Substitute v by 0
            f = Polynomial(f.set().subset0(v.index()))
    return f


def expanded_polynomial2general_polynomial(polynomial, new_variables, ring):
    r"""
    Returns a GeneralBooleanPolynomial associated to a Polynomial (polynomial) in
    additional variables (new_variables), where the
    GeneralBooleanPolynomial in obtained from the projections of polynomial
    to the factors of F_2^k (by specialization of the additional variables)
    """
    comps = [projection_of_expanded_polynomial(polynomial, e, new_variables)
        for e in new_variables]
    sum2 = GeneralBooleanPolynomial(len(new_variables), [0] * len(
        new_variables), Polynomial(0, ring))
    for i in range(len(new_variables)):
        sum2[i] = comps[i]
    return sum2


def reduce_general_boolean_polynomial(F, polynomial):
    r"""
    Computes the reduction of polynomial via the ideal given by F
    """
    r = polynomial
    s = len(F)
    k = len(polynomial)
    h = [GeneralBooleanPolynomial(len(polynomial), [0] * len(polynomial),
        Polynomial(0, self.ring))] * s
    # Indices i where leading monomial of F[i] divided leading monomial of r

    def calc_Is():
        ret = []
        for i in range(s):
            if (F[i].is_zero() and not r.is_zero()):
                continue
            if (F[i].lm()).divides(r.lm()):
                ret.append(i)
        return ret

    def calc_X():
        ret = []
        for i in range(len(Is)):
            ret.append((r.lm()).monomial() / (F[Is[i]].lm()).monomial())
        return ret

    def included(i, set):
        if i in set:
            return 1
        else:
            return 0

    Is = calc_Is()
    X = calc_X()

    for f in F:
        assert len(f) == len(polynomial)

    lm_polynomial = (polynomial.lm()).monomial()
    lc_polynomial_binary = polynomial.lc_binary()

    while len(Is) != 0:
        exp = [F[i].lm() for i in range(s)]
        matrix = [[included(i, F[j].lc_as_set()) for j in Is] for i in range(
            k)]

        # Compute solution
        coeff = [[0] * len(Is) for j in range(k)]
        for j in range(k):
            if lc_polynomial_binary[j]:
                coeff[j][0] = 1

        sum = GeneralBooleanPolynomial(k, [0] * k, Polynomial(0, self.ring))
        for i in range(len(Is)):
            c = [coeff[l][i] for l in range(k)]
            c_set = [l for l in range(k) if coeff[l][i] == 1]
            poly1 = GeneralBooleanPolynomial(k, c_set, X[i])

            sum += GeneralBooleanPolynomial(k, c_set, X[i]) * F[Is[i]].lt()

        if polynomial.lt() == sum:
            r -= sum
        else:
            break

        if r.is_zero():
            return r

        Is = calc_Is()
        X = calc_X()
    return r


def stratified(I):
    r"""
    Tests if I does no contain two polynomials with the same leading monomials
    """
    leading_monomials = []
    for p in I:
        if p.is_zero():
            continue
        lm = p.lm()
        if lm in leading_monomials:
            return False
        leading_monomials.append(lm)
    return True


def build_dict_from_array_of_extended_polynomials(A, e_vars):
    new_dict = {}
    for p in A:
        new_dict[p] = expanded_polynomial2general_polynomial(p, e_vars)
    return new_dict


def stratify_dict_I_gb_I(dict, e_vars, debug=0):
    r"""
    Wrapper (calls either stratify_dict_I_gb_I_our_alg or stratify_dict_I_gb_I_Inoue
    """
    #return stratify_dict_I_gb_I_Inoue(dict, e_vars, debug)
    return stratify_dict_I_gb_I_our_alg(dict, e_vars, debug)


def stratify_dict_I_gb_I_our_alg(dict, e_vars, debug=0):
    r"""
    Build a stratified Groebner bases for dict
    """
    # Ideal for the polynomials of the new basis
    A = []
    # and their leading monomials
    LMs = []
    new_dict = {}


    while len(dict) > 0:
        # We copy the keys of dict into I in order to sort them.
        # This way the elements are traversed in a unique order.
        I = sorted(dict.keys())
        p = I[0]

        p_gb = dict[p]

        del dict[p]

        if debug > 1:
            print("Processing p", p)
            print("A before proceeding", A)

        if p_gb.is_zero():
            LMs.append(Polynomial(0, self.ring))
            A.append(p)
            if debug > 1:
                print("Adding p that becomes zero")
            continue

        p_gb_lm = p_gb.lm()

        # If the leading monomial of p does not coincide with any of
        # polynomials already in A, then add p to A
        if not p_gb_lm in LMs:
            A.append(p)
            LMs.append(p_gb_lm)
            if debug > 1:
                print("Appending", p, "since its lm is not contained in A yet")
            continue

        # Index of p_gb_lm in LMs
        i = LMs.index(p_gb_lm)

        # Polynomial in A with the same lm as p
        b = A[i]

        # The Polynomial we want to add to A while keeping A stratified
        r = p

        b_gb = expanded_polynomial2general_polynomial(b, e_vars)
        r_gb = expanded_polynomial2general_polynomial(r, e_vars)

        # Leading coefficients as GeneralBooleanPolynomial
        lc_b_gb = b_gb.lc()
        lc_r_gb = r_gb.lc()
        unit = GeneralBooleanPolynomial(len(e_vars), [o for o in range(len(
            e_vars))], Polynomial(1, p.ring()))

        b1_gb = b_gb * (unit + lc_r_gb) + r_gb
        r2_gb = r_gb * (unit + lc_r_gb + lc_r_gb * lc_b_gb) + lc_r_gb * b_gb

        b1 = b1_gb.to_expanded_polynomial_ring(e_vars)
        r2 = r2_gb.to_expanded_polynomial_ring(e_vars)

        A[i] = b1
        if debug > 1:
            print("New polynomial in A (replaced)", A[i])

        if r2 != 0 and r2 not in A:
            dict[r2] = r2_gb

    return build_dict_from_array_of_extended_polynomials(A, e_vars)


def stratify_dict_I_gb_I_Inoue(dict, e_vars, debug=0):
    r"""
    Reimplementation of a simple algorithm of Inoue from BGSet
    """
    # Ideal for the polynomials of the new basis
    A = []
    new_dict = {}
    while len(dict) > 0:
        p = dict.keys()[0]
        p_gb = dict[p]
        del dict[p]

        if p_gb.is_zero():
            A.append(p)
            continue

        p_gb_lm = p_gb.lm()

        for f in dict.keys():
            f_gb = dict[f]

            if f_gb.is_zero():
                continue

            if p_gb_lm == f_gb.lm():
                p = p + f
                del dict[f]
        A.append(p)
    return build_dict_from_array_of_extended_polynomials(A, e_vars)

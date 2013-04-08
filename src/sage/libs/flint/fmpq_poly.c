///////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2010 Sebastian Pancratz                                     //
//                                                                           //
// Distributed under the terms of the GNU General Public License as          //
// published by the Free Software Foundation; either version 2 of the        //
// License, or (at your option) any later version.                           //
//                                                                           //
// http://www.gnu.org/licenses/                                              //
///////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
    extern "C" {
#endif

#include "fmpq_poly.h"

/**
 * \file     fmpq_poly.c
 * \brief    Fast implementation of the rational function field
 * \author   Sebastian Pancratz
 * \date     Mar 2010 -- July 2010
 * \version  0.1.8
 *
 * \mainpage
 *
 * <center>
 * A FLINT-based implementation of the rational polynomial ring.
 * </center>
 *
 * \section Overview
 *
 * The module <tt>fmpq_poly</tt> provides functions for performing
 * arithmetic on rational polynomials in \f$\mathbf{Q}[t]\f$, represented as
 * quotients of integer polynomials of type <tt>fmpz_poly_t</tt> and
 * denominators of type <tt>fmpz_t</tt>.  These functions start with the
 * prefix <tt>fmpq_poly_</tt>.
 *
 * Rational polynomials are stored in objects of type #fmpq_poly_t,
 * which is an array of #fmpq_poly_struct's of length one.  This
 * permits passing parameters of type #fmpq_poly_t by reference.
 * We also define the type #fmpq_poly_ptr to be a pointer to
 * #fmpq_poly_struct's.
 *
 * The representation of a rational polynomial as the quotient of an integer
 * polynomial and an integer denominator can be made canonical by demanding
 * the numerator and denominator to be coprime and the denominator to be
 * positive.  As the only special case, we represent the zero function as
 * \f$0/1\f$.  All arithmetic functions assume that the operands are in this
 * canonical form, and canonicalize their result.  If the numerator or
 * denominator is modified individually, for example using the methods in the
 * group \ref NumDen, it is the user's responsibility to canonicalize the
 * rational function using the function #fmpq_poly_canonicalize() if necessary.
 *
 * All methods support aliasing of their inputs and outputs \e unless
 * explicitly stated otherwise, subject to the following caveat.  If
 * different rational polynomials (as objects in memory, not necessarily in
 * the mathematical sense) share some of the underlying integer polynomials
 * or integers, the behaviour is undefined.
 *
 * \section Changelog  Version history
 * - 0.1.8
 *   - Extra checks for <tt>NULL</tt> returns of <tt>malloc</tt> calls
 * - 0.1.7
 *   - Explicit casts for return values of <tt>malloc</tt>
 *   - Changed a few calls to GMP and FLINT methods from <tt>_si</tt>
 *     to <tt>_ui</tt>
 * - 0.1.6
 *   - Added tons of testing code in the subdirectory <tt>test</tt>
 *   - Made a few changes throughout the code base indicated by the tests
 * - 0.1.5
 *   - Re-wrote the entire code to always initialise the denominator
 * - 0.1.4
 *   - Fixed a bug in #fmpq_poly_divrem()
 *   - Swapped calls to #fmpq_poly_degree to calls to #fmpq_poly_length in
 *     many places
 * - 0.1.3
 *   - Changed #fmpq_poly_inv()
 *   - Another sign check to <tt>fmpz_abs</tt> in #fmpq_poly_content()
 * - 0.1.2
 *   - Introduce a function #fmpq_poly_monic() and use this to simplify the
 *     code for the gcd and xgcd functions
 *   - Make further use of #fmpq_poly_is_zero()
 * - 0.1.1
 *   - Replaced a few sign checks and negations by <tt>fmpz_abs</tt>
 *   - Small changes to comments and the documentation
 *   - Moved some function bodies from #fmpq_poly.h to #fmpq_poly.c
 * - 0.1.0
 *   - First draft, based on the author's Sage code
 */

/**
 * \defgroup Definitions        Type definitions
 * \defgroup MemoryManagement   Memory management
 * \defgroup NumDen             Accessing numerator and denominator
 * \defgroup Assignment         Assignment and basic manipulation
 * \defgroup Coefficients       Setting and retrieving individual coefficients
 * \defgroup PolyParameters     Polynomial parameters
 * \defgroup Comparison         Comparison
 * \defgroup Addition           Addition and subtraction
 * \defgroup ScalarMul          Scalar multiplication and division
 * \defgroup Multiplication     Multiplication
 * \defgroup Division           Euclidean division
 * \defgroup Powering           Powering
 * \defgroup GCD                Greatest common divisor
 * \defgroup Derivative         Derivative
 * \defgroup Evaluation         Evaluation
 * \defgroup Content            Gaussian content
 * \defgroup Resultant          Resultant
 * \defgroup Composition        Composition
 * \defgroup Squarefree         Square-free test
 * \defgroup Subpolynomials     Subpolynomials
 * \defgroup StringConversions  String conversions and I/O
 */

///////////////////////////////////////////////////////////////////////////////
// Auxiliary functions                                                       //
///////////////////////////////////////////////////////////////////////////////

/**
 * Returns the number of digits in the decimal representation of \c n.
 */
unsigned int fmpq_poly_places(ulong n)
{
    unsigned int count;
    if (n == 0)
        return 1u;
    count = 0;
    while (n > 0)
    {
        n /= 10ul;
        count++;
    }
    return count;
}

///////////////////////////////////////////////////////////////////////////////
// Implementation                                                            //
///////////////////////////////////////////////////////////////////////////////

/**
 * \brief    Puts <tt>f</tt> into canonical form.
 * \ingroup  Definitions
 *
 * This method ensures that the denominator is coprime to the content
 * of the numerator polynomial, and that the denominator is positive.
 *
 * The optional parameter <tt>temp</tt> is the temporary variable that this
 * method would otherwise need to allocate.  If <tt>temp</tt> is provided as
 * an initialized <tt>fmpz_t</tt>, it is assumed \e without further checks
 * that it is large enough.  For example,
 * <tt>max(f-\>num-\>limbs, fmpz_size(f-\>den))</tt> limbs will certainly
 * suffice.
 */
void fmpq_poly_canonicalize(fmpq_poly_ptr f, fmpz_t temp)
{
    if (fmpz_is_one(f->den))
        return;

    if (fmpq_poly_is_zero(f))
        fmpz_set_ui(f->den, 1ul);

    else if (fmpz_is_m1(f->den))
    {
        fmpz_poly_neg(f->num, f->num);
        fmpz_set_ui(f->den, 1ul);
    }

    else
    {
        int tcheck;  /* Whether temp == NULL */
        //if (temp == NULL)
        //{
            tcheck = 1;
            temp = fmpz_init(FLINT_MAX(f->num->limbs, fmpz_size(f->den)));
        //}
        //else
        //    tcheck = 0;

        fmpz_poly_content(temp, f->num);
        fmpz_abs(temp, temp);

        if (fmpz_is_one(temp))
        {
            if (fmpz_sgn(f->den) < 0)
            {
                fmpz_poly_neg(f->num, f->num);
                fmpz_neg(f->den, f->den);
            }
        }
        else
        {
            fmpz_t tempgcd;
            tempgcd = fmpz_init(FLINT_MAX(f->num->limbs, fmpz_size(f->den)));

            fmpz_gcd(tempgcd, temp, f->den);
            fmpz_abs(tempgcd, tempgcd);

            if (fmpz_is_one(tempgcd))
            {
                if (fmpz_sgn(f->den) < 0)
                {
                    fmpz_poly_neg(f->num, f->num);
                    fmpz_neg(f->den, f->den);
                }
            }
            else
            {
                if (fmpz_sgn(f->den) < 0)
                    fmpz_neg(tempgcd, tempgcd);
                fmpz_poly_scalar_div_fmpz(f->num, f->num, tempgcd);
                fmpz_tdiv(temp, f->den, tempgcd);
                fmpz_set(f->den, temp);
            }

            fmpz_clear(tempgcd);
        }
        if (tcheck)
            fmpz_clear(temp);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Assignment and basic manipulation

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the additive inverse of \c op.
 */
void fmpq_poly_neg(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    if (rop == op)
    {
        fmpz_poly_neg(rop->num, op->num);
    }
    else
    {
        fmpz_poly_neg(rop->num, op->num);
        rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
        fmpz_set(rop->den, op->den);
    }
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the multiplicative inverse of \c op.
 *
 * Assumes that the element \c op is a unit.  Otherwise, an exception
 * is raised in the form of an <tt>abort</tt> statement.
 */
void fmpq_poly_inv(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    /* Assertion! */
    if (fmpq_poly_length(op) != 1ul)
    {
        printf("ERROR (fmpq_poly_inv).  Element is not a unit.\n");
        abort();
    }

    if (rop == op)
    {
        fmpz_t t;
        t = fmpz_init(rop->num->limbs);
        fmpz_set(t, rop->num->coeffs);
        fmpz_poly_zero(rop->num);
        if (fmpz_sgn(t) > 0)
        {
            fmpz_poly_set_coeff_fmpz(rop->num, 0, rop->den);
            rop->den = fmpz_realloc(rop->den, fmpz_size(t));
            fmpz_set(rop->den, t);
        }
        else
        {
            fmpz_neg(rop->den, rop->den);
            fmpz_poly_set_coeff_fmpz(rop->num, 0, rop->den);
            rop->den = fmpz_realloc(rop->den, fmpz_size(t));
            fmpz_neg(rop->den, t);
        }
        fmpz_clear(t);
    }
    else
    {
        fmpz_poly_zero(rop->num);
        fmpz_poly_set_coeff_fmpz(rop->num, 0, op->den);
        rop->den = fmpz_realloc(rop->den, fmpz_size(op->num->coeffs));
        fmpz_set(rop->den, op->num->coeffs);
        if (fmpz_sgn(rop->den) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_neg(rop->den, rop->den);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Setting/ retrieving individual coefficients

/**
 * \ingroup  Coefficients
 *
 * Obtains the <tt>i</tt>th coefficient from the polynomial <tt>op</tt>.
 */
void fmpq_poly_get_coeff_mpq(mpq_t rop, const fmpq_poly_ptr op, ulong i)
{
    fmpz_poly_get_coeff_mpz(mpq_numref(rop), op->num, i);
    fmpz_to_mpz(mpq_denref(rop), op->den);
    mpq_canonicalize(rop);
}

/**
 * \ingroup  Coefficients
 *
 * Sets the <tt>i</tt>th coefficient of the polynomial <tt>op</tt> to \c x.
 */
void fmpq_poly_set_coeff_fmpz(fmpq_poly_ptr rop, ulong i, const fmpz_t x)
{
    fmpz_t t;
    int canonicalize;

    /* Is the denominator 1?  This includes the case when rop == 0. */
    if (fmpz_is_one(rop->den))
    {
        fmpz_poly_set_coeff_fmpz(rop->num, i, x);
        return;
    }

    t = fmpz_poly_get_coeff_ptr(rop->num, i);
    canonicalize = !(t == NULL || fmpz_is_zero(t));

    t = fmpz_init(fmpz_size(x) + fmpz_size(rop->den));
    fmpz_mul(t, x, rop->den);
    fmpz_poly_set_coeff_fmpz(rop->num, i, t);
    fmpz_clear(t);
    if (canonicalize)
        fmpq_poly_canonicalize(rop, NULL);
}

/**
 * \ingroup  Coefficients
 *
 * Sets the <tt>i</tt>th coefficient of the polynomial <tt>op</tt> to \c x.
 */
void fmpq_poly_set_coeff_mpz(fmpq_poly_ptr rop, ulong i, const mpz_t x)
{
    mpz_t z;
    fmpz_t t;
    int canonicalize;

    /* Is the denominator 1?  This includes the case when rop == 0. */
    if (fmpz_is_one(rop->den))
    {
        fmpz_poly_set_coeff_mpz(rop->num, i, x);
        return;
    }

    t = fmpz_poly_get_coeff_ptr(rop->num, i);
    canonicalize = !(t == NULL || fmpz_is_zero(t));

    mpz_init(z);
    fmpz_to_mpz(z, rop->den);
    mpz_mul(z, z, x);
    fmpz_poly_set_coeff_mpz(rop->num, i, z);
    if (canonicalize)
        fmpq_poly_canonicalize(rop, NULL);
    mpz_clear(z);
}

/**
 * \ingroup  Coefficients
 *
 * Sets the <tt>i</tt>th coefficient of the polynomial <tt>op</tt> to \c x.
 */
void fmpq_poly_set_coeff_mpq(fmpq_poly_ptr rop, ulong i, const mpq_t x)
{
    fmpz_t oldc;
    mpz_t den, gcd;
    int canonicalize;

    mpz_init(den);
    mpz_init(gcd);

    fmpz_to_mpz(den, rop->den);
    mpz_gcd(gcd, den, mpq_denref(x));

    oldc = fmpz_poly_get_coeff_ptr(rop->num, i);
    canonicalize = !(oldc == NULL || fmpz_is_zero(oldc));

    if (mpz_cmp(mpq_denref(x), gcd) == 0)
    {
        mpz_divexact(den, den, gcd);
        mpz_mul(gcd, den, mpq_numref(x));
        fmpz_poly_set_coeff_mpz(rop->num, i, gcd);
    }
    else
    {
        mpz_t t;

        mpz_init(t);
        mpz_divexact(t, mpq_denref(x), gcd);
        fmpz_poly_scalar_mul_mpz(rop->num, rop->num, t);

        mpz_divexact(gcd, den, gcd);
        mpz_mul(gcd, gcd, mpq_numref(x));

        fmpz_poly_set_coeff_mpz(rop->num, i, gcd);

        mpz_mul(den, den, t);
        rop->den = fmpz_realloc(rop->den, mpz_size(den));
        mpz_to_fmpz(rop->den, den);

        mpz_clear(t);
    }

    if (canonicalize)
        fmpq_poly_canonicalize(rop, NULL);

    mpz_clear(den);
    mpz_clear(gcd);
}

/**
 * \ingroup  Coefficients
 *
 * Sets the <tt>i</tt>th coefficient of the polynomial <tt>op</tt> to \c x.
 */
void fmpq_poly_set_coeff_si(fmpq_poly_ptr rop, ulong i, long x)
{
    mpz_t z;
    fmpz_t t;
    int canonicalize;

    /* Is the denominator 1?  This includes the case when rop == 0. */
    if (fmpz_is_one(rop->den))
    {
        fmpz_poly_set_coeff_si(rop->num, i, x);
        return;
    }

    t = fmpz_poly_get_coeff_ptr(rop->num, i);
    canonicalize = !(t == NULL || fmpz_is_zero(t));

    mpz_init(z);
    fmpz_to_mpz(z, rop->den);
    mpz_mul_si(z, z, x);
    fmpz_poly_set_coeff_mpz(rop->num, i, z);
    if (canonicalize)
        fmpq_poly_canonicalize(rop, NULL);
    mpz_clear(z);
}

///////////////////////////////////////////////////////////////////////////////
// Comparison

/**
 * \brief    Returns whether \c op1 and \c op2 are equal.
 * \ingroup  Comparison
 *
 * Returns whether the two elements \c op1 and \c op2 are equal.
 */
int fmpq_poly_equal(const fmpq_poly_ptr op1, const fmpq_poly_ptr op2)
{
    return (op1->num->length == op2->num->length)
        && (fmpz_equal(op1->den, op2->den))
        && (fmpz_poly_equal(op1->num, op2->num));
}

/**
 * \brief   Compares the two polynomials <tt>left</tt> and <tt>right</tt>.
 * \ingroup Comparison
 *
 * Compares the two polnomials <tt>left</tt> and <tt>right</tt>, returning
 * <tt>-1</tt>, <tt>0</tt>, or <tt>1</tt> as <tt>left</tt> is less than,
 * equal to, or greater than <tt>right</tt>.
 *
 * The comparison is first done by degree and then, in case of a tie, by
 * the individual coefficients, beginning with the highest.
 */
int fmpq_poly_cmp(const fmpq_poly_ptr left, const fmpq_poly_ptr right)
{
    long degdiff, i;

    /* Quick check whether left and right are the same object. */
    if (left == right)
        return 0;

    i = fmpz_poly_degree(left->num);
    degdiff = i - fmpz_poly_degree(right->num);

    if (degdiff < 0)
        return -1;
    else if (degdiff > 0)
        return 1;
    else
    {
        int ans;
        ulong limbs;
        fmpz_t lcoeff, rcoeff;

        if (fmpz_equal(left->den, right->den))
        {
            if (i == -1)  /* left and right are both zero */
                return 0;
            limbs = FLINT_MAX(left->num->limbs, right->num->limbs) + 1;
            lcoeff = fmpz_init(limbs);
            while (fmpz_equal(fmpz_poly_get_coeff_ptr(left->num, i),
                fmpz_poly_get_coeff_ptr(right->num, i)) && i > 0)
                i--;
            fmpz_sub(lcoeff, fmpz_poly_get_coeff_ptr(left->num, i),
                fmpz_poly_get_coeff_ptr(right->num, i));
            ans = fmpz_sgn(lcoeff);
            fmpz_clear(lcoeff);
            return ans;
        }
        else if (fmpz_is_one(left->den))
        {
            limbs = FLINT_MAX(left->num->limbs + fmpz_size(right->den),
                right->num->limbs) + 1;
            lcoeff = fmpz_init(limbs);
            fmpz_mul(lcoeff, fmpz_poly_get_coeff_ptr(left->num, i), right->den);
            while (fmpz_equal(lcoeff, fmpz_poly_get_coeff_ptr(right->num, i))
                && i > 0)
            {
                i--;
                fmpz_mul(lcoeff, fmpz_poly_get_coeff_ptr(left->num, i),
                    right->den);
            }
            fmpz_sub(lcoeff, lcoeff, fmpz_poly_get_coeff_ptr(right->num, i));
            ans = fmpz_sgn(lcoeff);
            fmpz_clear(lcoeff);
            return ans;
        }
        else if (fmpz_is_one(right->den))
        {
            limbs = FLINT_MAX(left->num->limbs,
                right->num->limbs + fmpz_size(left->den)) + 1;
            rcoeff = fmpz_init(limbs);
            fmpz_mul(rcoeff, fmpz_poly_get_coeff_ptr(right->num, i), left->den);
            while (fmpz_equal(fmpz_poly_get_coeff_ptr(left->num, i), rcoeff)
                && i > 0)
            {
                i--;
                fmpz_mul(rcoeff, fmpz_poly_get_coeff_ptr(right->num, i),
                    left->den);
            }
            fmpz_sub(rcoeff, fmpz_poly_get_coeff_ptr(left->num, i), rcoeff);
            ans = fmpz_sgn(rcoeff);
            fmpz_clear(rcoeff);
            return ans;
        }
        else
        {
            limbs = FLINT_MAX(left->num->limbs + fmpz_size(right->den),
                right->num->limbs + fmpz_size(left->den)) + 1;
            lcoeff = fmpz_init(limbs);
            rcoeff = fmpz_init(right->num->limbs + fmpz_size(left->den));
            fmpz_mul(lcoeff, fmpz_poly_get_coeff_ptr(left->num, i), right->den);
            fmpz_mul(rcoeff, fmpz_poly_get_coeff_ptr(right->num, i), left->den);
            while (fmpz_equal(lcoeff, rcoeff) && i > 0)
            {
                i--;
                fmpz_mul(lcoeff, fmpz_poly_get_coeff_ptr(left->num, i),
                    right->den);
                fmpz_mul(rcoeff, fmpz_poly_get_coeff_ptr(right->num, i),
                    left->den);
            }
            fmpz_sub(lcoeff, lcoeff, rcoeff);
            ans = fmpz_sgn(lcoeff);
            fmpz_clear(lcoeff);
            fmpz_clear(rcoeff);
            return ans;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Addition/ subtraction

/**
 * \ingroup  Addition
 *
 * Sets \c rop to the sum of \c rop and \c op.
 *
 * \todo  This is currently implemented by creating a copy!
 */
void _fmpq_poly_add_in_place(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    if (rop == op)
    {
        fmpq_poly_scalar_mul_si(rop, rop, 2l);
        return;
    }

    if (fmpq_poly_is_zero(rop))
    {
        fmpq_poly_set(rop, op);
        return;
    }
    if (fmpq_poly_is_zero(op))
    {
        return;
    }

    /* Now we may assume that rop and op refer to distinct objects in        */
    /* memory and that both polynomials are non-zero.                        */
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_poly_add(t, rop, op);
    fmpq_poly_swap(rop, t);
    fmpq_poly_clear(t);
}

/**
 * \ingroup  Addition
 *
 * Sets \c rop to the sum of \c op1 and \c op2.
 */
void fmpq_poly_add(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2)
{
    if (op1 == op2)
    {
        fmpq_poly_scalar_mul_si(rop, op1, 2l);
        return;
    }

    if (rop == op1)
    {
        _fmpq_poly_add_in_place(rop, op2);
        return;
    }
    if (rop == op2)
    {
        _fmpq_poly_add_in_place(rop, op1);
        return;
    }

    /* From here on, we may assume that rop, op1 and op2 all refer to        */
    /* distinct objects in memory, although they may still be equal.         */

    if (fmpz_is_one(op1->den))
    {
        if (fmpz_is_one(op2->den))
        {
            fmpz_poly_add(rop->num, op1->num, op2->num);
            fmpz_set_ui(rop->den, 1ul);
        }
        else
        {
            fmpz_poly_scalar_mul_fmpz(rop->num, op1->num, op2->den);
            fmpz_poly_add(rop->num, rop->num, op2->num);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op2->den));
            fmpz_set(rop->den, op2->den);
        }
    }
    else
    {
        if (fmpz_is_one(op2->den))
        {
            fmpz_poly_scalar_mul_fmpz(rop->num, op2->num, op1->den);
            fmpz_poly_add(rop->num, op1->num, rop->num);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op1->den));
            fmpz_set(rop->den, op1->den);
        }
        else
        {
            fmpz_poly_t tpoly;
            fmpz_t tfmpz;
            ulong limbs;

            fmpz_poly_init(tpoly);

            limbs = fmpz_size(op1->den) + fmpz_size(op2->den);
            rop->den = fmpz_realloc(rop->den, limbs);

            limbs = FLINT_MAX(limbs, fmpz_size(op2->den) + op1->num->limbs);
            limbs = FLINT_MAX(limbs, fmpz_size(op1->den) + op2->num->limbs);
            tfmpz = fmpz_init(limbs);

            fmpz_gcd(rop->den, op1->den, op2->den);
            fmpz_tdiv(tfmpz, op2->den, rop->den);
            fmpz_poly_scalar_mul_fmpz(rop->num, op1->num, tfmpz);
            fmpz_tdiv(tfmpz, op1->den, rop->den);
            fmpz_poly_scalar_mul_fmpz(tpoly, op2->num, tfmpz);
            fmpz_poly_add(rop->num, rop->num, tpoly);
            fmpz_mul(rop->den, tfmpz, op2->den);

            fmpq_poly_canonicalize(rop, tfmpz);

            fmpz_poly_clear(tpoly);
            fmpz_clear(tfmpz);
        }
    }
}

/**
 * \ingroup  Addition
 *
 * Sets \c rop to the difference of \c rop and \c op.
 *
 * \note  This is implemented using the methods #fmpq_poly_neg() and
 *        #_fmpq_poly_add_in_place().
 */
void _fmpq_poly_sub_in_place(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    if (rop == op)
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpq_poly_neg(rop, rop);
    _fmpq_poly_add_in_place(rop, op);
    fmpq_poly_neg(rop, rop);
}

/**
 * \ingroup  Addition
 *
 * Sets \c rop to the difference of \c op1 and \c op2.
 */
void fmpq_poly_sub(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2)
{
    if (op1 == op2)
    {
        fmpq_poly_zero(rop);
        return;
    }
    if (rop == op1)
    {
        _fmpq_poly_sub_in_place(rop, op2);
        return;
    }
    if (rop == op2)
    {
        _fmpq_poly_sub_in_place(rop, op1);
        fmpq_poly_neg(rop, rop);
        return;
    }

    /* From here on, we know that rop, op1 and op2 refer to distinct objects */
    /* in memory, although as rational functions they may still be equal     */

    if (fmpz_is_one(op1->den))
    {
        if (fmpz_is_one(op2->den))
        {
            fmpz_poly_sub(rop->num, op1->num, op2->num);
            fmpz_set_ui(rop->den, 1ul);
        }
        else
        {
            fmpz_poly_scalar_mul_fmpz(rop->num, op1->num, op2->den);
            fmpz_poly_sub(rop->num, rop->num, op2->num);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op2->den));
            fmpz_set(rop->den, op2->den);
        }
    }
    else
    {
        if (fmpz_is_one(op2->den))
        {
            fmpz_poly_scalar_mul_fmpz(rop->num, op2->num, op1->den);
            fmpz_poly_sub(rop->num, op1->num, rop->num);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op1->den));
            fmpz_set(rop->den, op1->den);
        }
        else
        {
            fmpz_poly_t tpoly;
            fmpz_t tfmpz;
            ulong limbs;

            fmpz_poly_init(tpoly);

            limbs = fmpz_size(op1->den) + fmpz_size(op2->den);
            rop->den = fmpz_realloc(rop->den, limbs);

            limbs = FLINT_MAX(limbs, fmpz_size(op2->den) + op1->num->limbs);
            limbs = FLINT_MAX(limbs, fmpz_size(op1->den) + op2->num->limbs);
            tfmpz = fmpz_init(limbs);

            fmpz_gcd(rop->den, op1->den, op2->den);
            fmpz_tdiv(tfmpz, op2->den, rop->den);
            fmpz_poly_scalar_mul_fmpz(rop->num, op1->num, tfmpz);
            fmpz_tdiv(tfmpz, op1->den, rop->den);
            fmpz_poly_scalar_mul_fmpz(tpoly, op2->num, tfmpz);
            fmpz_poly_sub(rop->num, rop->num, tpoly);
            fmpz_mul(rop->den, tfmpz, op2->den);

            fmpq_poly_canonicalize(rop, tfmpz);

            fmpz_poly_clear(tpoly);
            fmpz_clear(tfmpz);
        }
    }
}

/**
 * \ingroup  Addition
 *
 * Sets \c rop to <tt>rop + op1 * op2</tt>.
 *
 * Currently, this method refers to the methods #fmpq_poly_mul() and
 * #fmpq_poly_add() to form the result in the naive way.
 *
 * \todo  Implement this method more efficiently.
 */
void fmpq_poly_addmul(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2)
{
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_poly_mul(t, op1, op2);
    fmpq_poly_add(rop, rop, t);
    fmpq_poly_clear(t);
}

/**
 * \ingroup  Addition
 *
 * Sets \c rop to <tt>rop - op1 * op2</tt>.
 *
 * Currently, this method refers to the methods #fmpq_poly_mul() and
 * #fmpq_poly_sub() to form the result in the naive way.
 *
 * \todo  Implement this method more efficiently.
 */
void fmpq_poly_submul(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2)
{
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_poly_mul(t, op1, op2);
    fmpq_poly_sub(rop, rop, t);
    fmpq_poly_clear(t);
}

///////////////////////////////////////////////////////////////////////////////
// Scalar multiplication and division

/**
 * \ingroup  ScalarMul
 *
 * Sets \c rop to the scalar product of \c op and the integer \c x.
 */
void fmpq_poly_scalar_mul_si(fmpq_poly_ptr rop, const fmpq_poly_ptr op, long x)
{
    if (fmpz_is_one(op->den))
    {
        fmpz_poly_scalar_mul_si(rop->num, op->num, x);
        fmpz_set_ui(rop->den, 1ul);
    }
    else
    {
        fmpz_t fx, g;

        g  = fmpz_init(fmpz_size(op->den));
        fx = fmpz_init(1);
        fmpz_set_si(fx, x);
        fmpz_gcd(g, op->den, fx);
        fmpz_abs(g, g);
        if (fmpz_is_one(g))
        {
            fmpz_poly_scalar_mul_si(rop->num, op->num, x);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
            fmpz_set(rop->den, op->den);
        }
        else
        {
            if (rop == op)
            {
                fmpz_t t;
                t = fmpz_init(fmpz_size(op->den));
                fmpz_tdiv(t, fx, g);
                fmpz_poly_scalar_mul_fmpz(rop->num, op->num, t);
                fmpz_tdiv(t, op->den, g);
                fmpz_clear(rop->den);  // FIXME  Fixed
                rop->den = t;
            }
            else
            {
                rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
                fmpz_tdiv(rop->den, fx, g);
                fmpz_poly_scalar_mul_fmpz(rop->num, op->num, rop->den);
                fmpz_tdiv(rop->den, op->den, g);
            }
        }
        fmpz_clear(g);
        fmpz_clear(fx);
    }
}

/**
 * \ingroup  ScalarMul
 *
 * Sets \c rop to the scalar multiple of \c op with the \c mpz_t \c x.
 */
void fmpq_poly_scalar_mul_mpz(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpz_t x)
{
    fmpz_t g, y;

    if (fmpz_is_one(op->den))
    {
        fmpz_poly_scalar_mul_mpz(rop->num, op->num, x);
        fmpz_set_ui(rop->den, 1UL);
        return;
    }

    if (mpz_cmpabs_ui(x, 1) == 0)
    {
        if (mpz_sgn(x) > 0)
            fmpq_poly_set(rop, op);
        else
            fmpq_poly_neg(rop, op);
        return;
    }
    if (mpz_sgn(x) == 0)
    {
        fmpq_poly_zero(rop);
        return;
    }

    g = fmpz_init(FLINT_MAX(fmpz_size(op->den), mpz_size(x)));
    y = fmpz_init(mpz_size(x));
    mpz_to_fmpz(y, x);

    fmpz_gcd(g, op->den, y);

    if (fmpz_is_one(g))
    {
        fmpz_poly_scalar_mul_fmpz(rop->num, op->num, y);
        if (rop != op)
        {
            rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
            fmpz_set(rop->den, op->den);
        }
    }
    else
    {
        fmpz_t t;
        t = fmpz_init(FLINT_MAX(fmpz_size(y), fmpz_size(op->den)));
        fmpz_divides(t, y, g);
        fmpz_poly_scalar_mul_fmpz(rop->num, op->num, t);
        fmpz_divides(t, op->den, g);
        fmpz_clear(rop->den);
        rop->den = t;
    }

    fmpz_clear(g);
    fmpz_clear(y);
}

/**
 * \ingroup  ScalarMul
 *
 * Sets \c rop to the scalar multiple of \c op with the \c mpq_t \c x.
 */
void fmpq_poly_scalar_mul_mpq(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpq_t x)
{
    fmpz_poly_scalar_mul_mpz(rop->num, op->num, mpq_numref(x));
    if (fmpz_is_one(op->den))
    {
        rop->den = fmpz_realloc(rop->den, mpz_size(mpq_denref(x)));
        mpz_to_fmpz(rop->den, mpq_denref(x));
    }
    else
    {
        fmpz_t s, t;
        s = fmpz_init(mpz_size(mpq_denref(x)));
        t = fmpz_init(fmpz_size(op->den));
        mpz_to_fmpz(s, mpq_denref(x));
        fmpz_set(t, op->den);
        rop->den = fmpz_realloc(rop->den, fmpz_size(s) + fmpz_size(t));
        fmpz_mul(rop->den, s, t);
        fmpz_clear(s);
        fmpz_clear(t);
    }
    fmpq_poly_canonicalize(rop, NULL);
}

/**
 * /ingroup  ScalarMul
 *
 * Divides \c rop by the integer \c x.
 *
 * Assumes that \c x is non-zero.  Otherwise, an exception is raised in the
 * form of an <tt>abort</tt> statement.
 */
void _fmpq_poly_scalar_div_si_in_place(fmpq_poly_ptr rop, long x)
{
    fmpz_t cont, fx, gcd;

    /* Assertion! */
    if (x == 0l)
    {
        printf("ERROR (_fmpq_poly_scalar_div_si_in_place).  Division by zero.\n");
        abort();
    }

    if (x == 1l)
    {
        return;
    }

    cont = fmpz_init(rop->num->limbs);
    fmpz_poly_content(cont, rop->num);
    fmpz_abs(cont, cont);

    if (fmpz_is_one(cont))
    {
        if (x > 0l)
        {
            if (fmpz_is_one(rop->den))
            {
                fmpz_set_si(rop->den, x);
            }
            else
            {
                fmpz_t t;
                t = fmpz_init(fmpz_size(rop->den) + 1);
                fmpz_mul_ui(t, rop->den, (ulong) x);
                fmpz_clear(rop->den);  // FIXME  Fixed
                rop->den = t;
            }
        }
        else
        {
            fmpz_poly_neg(rop->num, rop->num);
            if (fmpz_is_one(rop->den))
            {
                fmpz_set_si(rop->den, -x);
            }
            else
            {
                fmpz_t t;
                t = fmpz_init(fmpz_size(rop->den) + 1);
                fmpz_mul_ui(t, rop->den, (ulong) -x);
                fmpz_clear(rop->den);  // FIXME  Fixed
                rop->den = t;
            }
        }
        fmpz_clear(cont);
        return;
    }

    fx = fmpz_init(1);
    fmpz_set_si(fx, x);

    gcd = fmpz_init(FLINT_MAX(rop->num->limbs, fmpz_size(fx)));
    fmpz_gcd(gcd, cont, fx);
    fmpz_abs(gcd, gcd);

    if (fmpz_is_one(gcd))
    {
        if (x > 0l)
        {
            if (fmpz_is_one(rop->den))
            {
                fmpz_set_si(rop->den, x);
            }
            else
            {
                fmpz_t t;
                t = fmpz_init(fmpz_size(rop->den) + 1);
                fmpz_mul_ui(t, rop->den, (ulong) x);
                fmpz_clear(rop->den);  // FIXME  Fixed
                rop->den = t;
            }
        }
        else
        {
            fmpz_poly_neg(rop->num, rop->num);
            if (fmpz_is_one(rop->den))
            {
                fmpz_set_si(rop->den, -x);
            }
            else
            {
                fmpz_t t;
                t = fmpz_init(fmpz_size(rop->den) + 1);
                fmpz_mul_ui(t, rop->den, (ulong) -x);
                fmpz_clear(rop->den);  // FIXME  Fixed
                rop->den = t;
            }
        }
    }
    else
    {
        fmpz_poly_scalar_div_fmpz(rop->num, rop->num, gcd);
        if (fmpz_is_one(rop->den))
        {
            fmpz_tdiv(rop->den, fx, gcd);
        }
        else
        {
            fmpz_tdiv(cont, fx, gcd);
            fx = fmpz_realloc(fx, fmpz_size(rop->den));
            fmpz_set(fx, rop->den);
            rop->den = fmpz_realloc(rop->den, fmpz_size(rop->den) + 1);
            fmpz_mul(rop->den, fx, cont);
        }
        if (x < 0l)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_neg(rop->den, rop->den);
        }
    }

    fmpz_clear(cont);
    fmpz_clear(fx);
    fmpz_clear(gcd);
}


/**
 * \ingroup  ScalarMul
 *
 * Sets \c rop to the scalar multiple of \c op with the multiplicative inverse
 * of the integer \c x.
 *
 * Assumes that \c x is non-zero.  Otherwise, an exception is raised in the
 * form of an <tt>abort</tt> statement.
 */
void fmpq_poly_scalar_div_si(fmpq_poly_ptr rop, const fmpq_poly_ptr op, long x)
{
    fmpz_t cont, fx, gcd;

    /* Assertion! */
    if (x == 0l)
    {
        printf("ERROR (fmpq_poly_scalar_div_si).  Division by zero.\n");
        abort();
    }

    if (rop == op)
    {
        _fmpq_poly_scalar_div_si_in_place(rop, x);
        return;
    }

    /* From here on, we may assume that rop and op denote two different      */
    /* rational polynomials (as objects in memory).                          */

    if (x == 1l)
    {
        fmpq_poly_set(rop, op);
        return;
    }

    cont = fmpz_init(op->num->limbs);
    fmpz_poly_content(cont, op->num);
    fmpz_abs(cont, cont);

    if (fmpz_is_one(cont))
    {
        if (x > 0l)
        {
            fmpz_poly_set(rop->num, op->num);
            if (fmpz_is_one(op->den))
            {
                fmpz_set_si(rop->den, x);
            }
            else
            {
                rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + 1);
                fmpz_mul_ui(rop->den, op->den, (ulong) x);
            }
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            if (fmpz_is_one(op->den))
            {
                fmpz_set_si(rop->den, -x);
            }
            else
            {
                rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + 1);
                fmpz_mul_ui(rop->den, op->den, (ulong) -x);
            }
        }
        fmpz_clear(cont);
        return;
    }

    fx = fmpz_init(1);
    fmpz_set_si(fx, x);

    gcd = fmpz_init(FLINT_MAX(op->num->limbs, fmpz_size(fx)));
    fmpz_gcd(gcd, cont, fx);
    fmpz_abs(gcd, gcd);

    if (fmpz_is_one(gcd))
    {
        if (x > 0l)
        {
            fmpz_poly_set(rop->num, op->num);
            if (fmpz_is_one(op->den))
            {
                fmpz_set_si(rop->den, x);
            }
            else
            {
                rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + 1);
                fmpz_mul_ui(rop->den, op->den, (ulong) x);
            }
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            if (fmpz_is_one(op->den))
            {
                fmpz_set_si(rop->den, -x);
            }
            else
            {
                rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + 1);
                fmpz_mul_ui(rop->den, op->den, (ulong) -x);
            }
        }
    }
    else
    {
        fmpz_poly_scalar_div_fmpz(rop->num, op->num, gcd);
        if (fmpz_is_one(op->den))
        {
            fmpz_tdiv(rop->den, fx, gcd);
        }
        else
        {
            rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + 1);
            fmpz_tdiv(cont, fx, gcd);  /* fx and gcd are word-sized */
            fmpz_mul(rop->den, op->den, cont);
        }
        if (x < 0l)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_neg(rop->den, rop->den);
        }
    }

    fmpz_clear(cont);
    fmpz_clear(fx);
    fmpz_clear(gcd);
}

/**
 * \ingroup  ScalarMul
 *
 * Sets \c rop to the scalar multiple of \c op with the multiplicative inverse
 * of the integer \c x.
 *
 * Assumes that \c x is non-zero.  Otherwise, an exception is raised in the
 * form of an <tt>abort</tt> statement.
 */
void _fmpq_poly_scalar_div_mpz_in_place(fmpq_poly_ptr rop, const mpz_t x)
{
    fmpz_t cont, fx, gcd;

    /* Assertion! */
    if (mpz_sgn(x) == 0)
    {
        printf("ERROR (_fmpq_poly_scalar_div_mpz_in_place).  Division by zero.\n");
        abort();
    }

    if (mpz_cmp_ui(x, 1) == 0)
        return;

    cont = fmpz_init(rop->num->limbs);
    fmpz_poly_content(cont, rop->num);
    fmpz_abs(cont, cont);

    if (fmpz_is_one(cont))
    {
        if (fmpz_is_one(rop->den))
        {
            rop->den = fmpz_realloc(rop->den, mpz_size(x));
            mpz_to_fmpz(rop->den, x);
        }
        else
        {
            fmpz_t t;
            fx = fmpz_init(mpz_size(x));
            mpz_to_fmpz(fx, x);
            t = fmpz_init(fmpz_size(rop->den) + mpz_size(x));
            fmpz_mul(t, rop->den, fx);
            fmpz_clear(rop->den);  // FIXME  Fixed
            rop->den = t;
            fmpz_clear(fx);
        }
        if (mpz_sgn(x) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_neg(rop->den, rop->den);
        }
        fmpz_clear(cont);
        return;
    }

    fx = fmpz_init(mpz_size(x));
    mpz_to_fmpz(fx, x);

    gcd = fmpz_init(FLINT_MAX(rop->num->limbs, fmpz_size(fx)));
    fmpz_gcd(gcd, cont, fx);
    fmpz_abs(gcd, gcd);

    if (fmpz_is_one(gcd))
    {
        if (fmpz_is_one(rop->den))
        {
            rop->den = fmpz_realloc(rop->den, mpz_size(x));
            mpz_to_fmpz(rop->den, x);
        }
        else
        {
            fmpz_t t;
            t = fmpz_init(fmpz_size(rop->den) + mpz_size(x));
            fmpz_mul(t, rop->den, fx);
            fmpz_clear(rop->den);  // FIXME  Fixed
            rop->den = t;
        }
        if (mpz_sgn(x) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_neg(rop->den, rop->den);
        }
    }
    else
    {
        fmpz_poly_scalar_div_fmpz(rop->num, rop->num, gcd);
        if (fmpz_is_one(rop->den))
        {
            rop->den = fmpz_realloc(rop->den, fmpz_size(fx));
            fmpz_tdiv(rop->den, fx, gcd);
        }
        else
        {
            cont = fmpz_realloc(cont, fmpz_size(fx));
            fmpz_tdiv(cont, fx, gcd);
            fx = fmpz_realloc(fx, fmpz_size(rop->den));
            fmpz_set(fx, rop->den);
            rop->den = fmpz_realloc(rop->den, fmpz_size(fx) + fmpz_size(cont));
            fmpz_mul(rop->den, fx, cont);
        }
        if (mpz_sgn(x) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_neg(rop->den, rop->den);
        }
    }

    fmpz_clear(cont);
    fmpz_clear(fx);
    fmpz_clear(gcd);
}

/**
 * \ingroup  ScalarMul
 *
 * Sets \c rop to the scalar multiple of \c op with the multiplicative inverse
 * of the integer \c x.
 *
 * Assumes that \c x is non-zero.  Otherwise, an exception is raised in the
 * form of an <tt>abort</tt> statement.
 */
void fmpq_poly_scalar_div_mpz(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpz_t x)
{
    fmpz_t cont, fx, gcd;

    /* Assertion! */
    if (mpz_sgn(x) == 0)
    {
        printf("ERROR (fmpq_poly_scalar_div_mpz).  Division by zero.\n");
        abort();
    }

    if (rop == op)
    {
        _fmpq_poly_scalar_div_mpz_in_place(rop, x);
        return;
    }

    /* From here on, we may assume that rop and op denote two different      */
    /* rational polynomials (as objects in memory).                          */

    if (mpz_cmp_ui(x, 1) == 0)
    {
        fmpq_poly_set(rop, op);
        return;
    }

    cont = fmpz_init(op->num->limbs);
    fmpz_poly_content(cont, op->num);
    fmpz_abs(cont, cont);

    if (fmpz_is_one(cont))
    {
        if (fmpz_is_one(op->den))
        {
            rop->den = fmpz_realloc(rop->den, mpz_size(x));
            mpz_to_fmpz(rop->den, x);
        }
        else
        {
            fmpz_t t;
            t = fmpz_init(mpz_size(x));
            mpz_to_fmpz(t, x);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + fmpz_size(t));
            fmpz_mul(rop->den, op->den, t);
            fmpz_clear(t);
        }
        if (mpz_sgn(x) > 0)
        {
            fmpz_poly_set(rop->num, op->num);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            fmpz_neg(rop->den, rop->den);
        }
        fmpz_clear(cont);
        return;
    }

    fx = fmpz_init(mpz_size(x));
    mpz_to_fmpz(fx, x);

    gcd = fmpz_init(FLINT_MAX(op->num->limbs, fmpz_size(fx)));
    fmpz_gcd(gcd, cont, fx);
    fmpz_abs(gcd, gcd);

    if (fmpz_is_one(gcd))
    {
        if (fmpz_is_one(op->den))
        {
            rop->den = fmpz_realloc(rop->den, mpz_size(x));
            mpz_to_fmpz(rop->den, x);
        }
        else
        {
            fmpz_t t;
            t = fmpz_init(mpz_size(x));
            mpz_to_fmpz(t, x);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + fmpz_size(t));
            fmpz_mul(rop->den, op->den, t);
            fmpz_clear(t);
        }
        if (mpz_sgn(x) > 0)
        {
            fmpz_poly_set(rop->num, op->num);
        }
        else
        {
            fmpz_poly_neg(rop->num, op->num);
            fmpz_neg(rop->den, rop->den);
        }
    }
    else
    {
        fmpz_poly_scalar_div_fmpz(rop->num, op->num, gcd);
        if (fmpz_is_one(op->den))
        {
            rop->den = fmpz_realloc(rop->den, fmpz_size(fx));
            fmpz_tdiv(rop->den, fx, gcd);
        }
        else
        {
            cont = fmpz_realloc(cont, fmpz_size(fx));
            fmpz_tdiv(cont, fx, gcd);
            rop->den = fmpz_realloc(rop->den, fmpz_size(op->den) + fmpz_size(cont));
            fmpz_mul(rop->den, op->den, cont);
        }
        if (mpz_sgn(x) < 0)
        {
            fmpz_poly_neg(rop->num, rop->num);
            fmpz_neg(rop->den, rop->den);
        }
    }

    fmpz_clear(cont);
    fmpz_clear(fx);
    fmpz_clear(gcd);
}

/**
 * \ingroup  ScalarMul
 *
 * Sets \c rop to the scalar multiple of \c op with the multiplicative inverse
 * of the rational \c x.
 *
 * Assumes that the rational \c x is in lowest terms and non-zero.  If the
 * rational is not in lowest terms, the resulting value of <tt>rop</tt> is
 * undefined.  If <tt>x</tt> is zero, an exception is raised in the form
 * of an <tt>abort</tt> statement.
 */
void fmpq_poly_scalar_div_mpq(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpq_t x)
{
    /* Assertion! */
    if (mpz_sgn(mpq_numref(x)) == 0)
    {
        printf("ERROR (fmpq_poly_scalar_div_mpq).  Division by zero.\n");
        abort();
    }

    fmpz_poly_scalar_mul_mpz(rop->num, op->num, mpq_denref(x));
    if (fmpz_is_one(op->den))
    {
        rop->den = fmpz_realloc(rop->den, mpz_size(mpq_numref(x)));
        mpz_to_fmpz(rop->den, mpq_numref(x));
    }
    else
    {
        fmpz_t s, t;
        s = fmpz_init(mpz_size(mpq_numref(x)));
        t = fmpz_init(fmpz_size(op->den));
        mpz_to_fmpz(s, mpq_numref(x));
        fmpz_set(t, op->den);
        rop->den = fmpz_realloc(rop->den, fmpz_size(s) + fmpz_size(t));
        fmpz_mul(rop->den, s, t);
        fmpz_clear(s);
        fmpz_clear(t);
    }
    fmpq_poly_canonicalize(rop, NULL);
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication

/**
 * \ingroup  Multiplication
 *
 * Multiplies <tt>rop</tt> by <tt>op</tt>.
 */
void _fmpq_poly_mul_in_place(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    fmpq_poly_t t;

    if (rop == op)
    {
        fmpz_poly_power(rop->num, op->num, 2ul);
        if (!fmpz_is_one(rop->den))
        {
            rop->den = fmpz_realloc(rop->den, 2 * fmpz_size(rop->den));
            fmpz_pow_ui(rop->den, op->den, 2ul);
        }
        return;
    }

    if (fmpq_poly_is_zero(rop) || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    /* From here on, rop and op point to two different objects in memory,    */
    /* and these are both non-zero rational polynomials.                     */

    fmpq_poly_init(t);
    fmpq_poly_mul(t, rop, op);
    fmpq_poly_swap(rop, t);
    fmpq_poly_clear(t);
}

/**
 * \ingroup  Multiplication
 *
 * Sets \c rop to the product of \c op1 and \c op2.
 */
void fmpq_poly_mul(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2)
{
    fmpz_t gcd1, gcd2;

    if (op1 == op2)
    {
        fmpz_poly_power(rop->num, op1->num, 2ul);
        if (fmpz_is_one(op1->den))
        {
            fmpz_set_ui(rop->den, 1);
        }
        else
        {
            if (rop == op1)
            {
                fmpz_t t;
                t = fmpz_init(2 * fmpz_size(op1->den));
                fmpz_pow_ui(t, op1->den, 2ul);
                fmpz_clear(rop->den);  // FIXME  Fixed
                rop->den = t;
            }
            else
            {
                rop->den = fmpz_realloc(rop->den, 2 * fmpz_size(op1->den));
                fmpz_pow_ui(rop->den, op1->den, 2ul);
            }
        }
        return;
    }

    if (rop == op1)
    {
        _fmpq_poly_mul_in_place(rop, op2);
        return;
    }
    if (rop == op2)
    {
        _fmpq_poly_mul_in_place(rop, op1);
        return;
    }

    /* From here on, we may assume that rop, op1 and op2 all refer to        */
    /* distinct objects in memory, although they may still be equal.         */

    gcd1 = NULL;
    gcd2 = NULL;

    if (!fmpz_is_one(op2->den))
    {
        gcd1 = fmpz_init(FLINT_MAX(op1->num->limbs, fmpz_size(op2->den)));
        fmpz_poly_content(gcd1, op1->num);
        if (!fmpz_is_one(gcd1))
        {
            fmpz_t t;
            t = fmpz_init(FLINT_MAX(fmpz_size(gcd1), fmpz_size(op2->den)));
            fmpz_gcd(t, gcd1, op2->den);
            fmpz_clear(gcd1);
            gcd1 = t;
        }
    }
    if (!fmpz_is_one(op1->den))
    {
        gcd2 = fmpz_init(FLINT_MAX(op2->num->limbs, fmpz_size(op1->den)));
        fmpz_poly_content(gcd2, op2->num);
        if (!fmpz_is_one(gcd2))
        {
            fmpz_t t;
            t = fmpz_init(FLINT_MAX(fmpz_size(gcd2), fmpz_size(op1->den)));
            fmpz_gcd(t, gcd2, op1->den);
            fmpz_clear(gcd2);
            gcd2 = t;
        }
    }

    /* TODO:  If gcd1 and gcd2 are very large compared to the degrees of  */
    /* poly1 and poly2, we might want to create copies of the polynomials */
    /* and divide out the common factors *before* the multiplication.     */

    fmpz_poly_mul(rop->num, op1->num, op2->num);
    rop->den = fmpz_realloc(rop->den, fmpz_size(op1->den)
             + fmpz_size(op2->den));
    fmpz_mul(rop->den, op1->den, op2->den);

    if (gcd1 == NULL || fmpz_is_one(gcd1))
    {
        if (gcd2 != NULL && !fmpz_is_one(gcd2))
        {
            fmpz_t h;
            h = fmpz_init(fmpz_size(rop->den));
            fmpz_poly_scalar_div_fmpz(rop->num, rop->num, gcd2);
            fmpz_divides(h, rop->den, gcd2);
            fmpz_clear(rop->den);
            rop->den = h;
        }
    }
    else
    {
        if (gcd2 == NULL || fmpz_is_one(gcd2))
        {
            fmpz_t h;
            h = fmpz_init(fmpz_size(rop->den));
            fmpz_poly_scalar_div_fmpz(rop->num, rop->num, gcd1);
            fmpz_divides(h, rop->den, gcd1);
            fmpz_clear(rop->den);
            rop->den = h;
        }
        else
        {
            fmpz_t g, h;
            g = fmpz_init(fmpz_size(gcd1) + fmpz_size(gcd2));
            h = fmpz_init(fmpz_size(rop->den));
            fmpz_mul(g, gcd1, gcd2);
            fmpz_poly_scalar_div_fmpz(rop->num, rop->num, g);
            fmpz_divides(h, rop->den, g);
            fmpz_clear(rop->den);
            rop->den = h;
            fmpz_clear(g);
        }
    }

    if (gcd1 != NULL) fmpz_clear(gcd1);
    if (gcd2 != NULL) fmpz_clear(gcd2);
}

///////////////////////////////////////////////////////////////////////////////
// Division

/**
 * \ingroup  Division
 *
 * Returns the quotient of the Euclidean division of \c a by \c b.
 *
 * Assumes that \c b is non-zero.  Otherwise, an exception is raised in the
 * form of an <tt>abort</tt> statement.
 */
void fmpq_poly_floordiv(fmpq_poly_ptr q, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    ulong m;
    fmpz_t lead;

    /* Assertion! */
    if (fmpq_poly_is_zero(b))
    {
        printf("ERROR (fmpq_poly_floordiv).  Division by zero.\n");
        abort();
    }

    /* Catch the case when a and b are the same objects in memory. */
    /* As of FLINT version 1.5.1, there is a bug in this case.     */
    if (a == b)
    {
        fmpq_poly_set_si(q, 1);
        return;
    }

    /* Deal with the various other cases of aliasing. */
    if (q == a | q == b)
    {
        fmpq_poly_t tpoly;
        fmpq_poly_init(tpoly);
        fmpq_poly_floordiv(tpoly, a, b);
        fmpq_poly_swap(q, tpoly);
        fmpq_poly_clear(tpoly);
        return;
    }

    /* Deal separately with the case deg(b) = 0. */
    if (fmpq_poly_length(b) == 1ul)
    {
        lead = fmpz_poly_get_coeff_ptr(b->num, 0);
        if (fmpz_is_one(a->den))
        {
            if (fmpz_is_one(b->den))  /* a->den == b->den == 1 */
                fmpz_poly_set(q->num, a->num);
            else  /* a->den == 1, b->den > 1 */
                fmpz_poly_scalar_mul_fmpz(q->num, a->num, b->den);

            q->den = fmpz_realloc(q->den, fmpz_size(lead));
            fmpz_set(q->den, lead);
            fmpq_poly_canonicalize(q, NULL);
        }
        else
        {
            if (fmpz_is_one(b->den))  /* a->den > 1, b->den == 1 */
                fmpz_poly_set(q->num, a->num);
            else  /* a->den, b->den > 1 */
                fmpz_poly_scalar_mul_fmpz(q->num, a->num, b->den);

            q->den = fmpz_realloc(q->den, fmpz_size(a->den) + fmpz_size(lead));
            fmpz_mul(q->den, a->den, lead);
            fmpq_poly_canonicalize(q, NULL);
        }
        return;
    }

    /* General case..                            */
    /* Set q to b->den q->num / (a->den lead^m). */

    /* Since the fmpz_poly_pseudo_div method is broken as of FLINT 1.5.0, we */
    /* use the fmpz_poly_pseudo_divrem method instead.                       */
    /* fmpz_poly_pseudo_div(q->num, &m, a->num, b->num);                     */
    {
        fmpz_poly_t r;
        fmpz_poly_init(r);
        fmpz_poly_pseudo_divrem(q->num, r, &m, a->num, b->num);
        fmpz_poly_clear(r);
    }

    lead = fmpz_poly_lead(b->num);

    if (!fmpz_is_one(b->den))
        fmpz_poly_scalar_mul_fmpz(q->num, q->num, b->den);

    /* Case 1:  lead^m is 1 */
    if (fmpz_is_one(lead) || m == 0ul || (fmpz_is_m1(lead) & m % 2 == 0))
    {
        q->den = fmpz_realloc(q->den, fmpz_size(a->den));
        fmpz_set(q->den, a->den);
        fmpq_poly_canonicalize(q, NULL);
    }
    /* Case 2:  lead^m is -1 */
    else if (fmpz_is_m1(lead) & m % 2)
    {
        fmpz_poly_neg(q->num, q->num);
        q->den = fmpz_realloc(q->den, fmpz_size(a->den));
        fmpz_set(q->den, a->den);
        fmpq_poly_canonicalize(q, NULL);
    }
    /* Case 3:  lead^m is not +-1 */
    else
    {
        ulong limbs;
        limbs = m * fmpz_size(lead);

        if (fmpz_is_one(a->den))
        {
            q->den = fmpz_realloc(q->den, limbs);
            fmpz_pow_ui(q->den, lead, m);
        }
        else
        {
            fmpz_t t;
            t = fmpz_init(limbs);
            q->den = fmpz_realloc(q->den, limbs + fmpz_size(a->den));
            fmpz_pow_ui(t, lead, m);
            fmpz_mul(q->den, t, a->den);
            fmpz_clear(t);
        }
        fmpq_poly_canonicalize(q, NULL);
    }
}

/**
 * \ingroup  Division
 *
 * Sets \c r to the remainder of the Euclidean division of \c a by \c b.
 *
 * Assumes that \c b is non-zero.  Otherwise, an exception is raised in the
 * form of an <tt>abort</tt> statement.
 */
void fmpq_poly_mod(fmpq_poly_ptr r, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    ulong m;
    fmpz_t lead;

    /* Assertion! */
    if (fmpq_poly_is_zero(b))
    {
        printf("ERROR (fmpq_poly_mod).  Division by zero.\n");
        abort();
    }

    /* Catch the case when a and b are the same objects in memory. */
    /* As of FLINT version 1.5.1, there is a bug in this case.     */
    if (a == b)
    {
        fmpq_poly_set_si(r, 0);
        return;
    }

    /* Deal with the various other cases of aliasing. */
    if (r == a | r == b)
    {
        fmpq_poly_t tpoly;
        fmpq_poly_init(tpoly);
        fmpq_poly_mod(tpoly, a, b);
        fmpq_poly_swap(r, tpoly);
        fmpq_poly_clear(tpoly);
        return;
    }

    /* Deal separately with the case deg(b) = 0. */
    if (fmpq_poly_length(b) == 1ul)
    {
        fmpz_poly_zero(r->num);
        fmpz_set_ui(r->den, 1);
        return;
    }

    /* General case..                     */
    /* Set r to r->num / (a->den lead^m). */


    /* As of FLINT 1.5.0, the fmpz_poly_pseudo_mod method is not             */
    /* asymptotically fast and hence we swap to the fmpz_poly_pseudo_divrem  */
    /* method if one of the operands' lengths is at least 32.                */
    if (fmpq_poly_length(a) < 32 && fmpq_poly_length(b) < 32)
        fmpz_poly_pseudo_rem(r->num, &m, a->num, b->num);
    else
    {
        fmpz_poly_t q;
        fmpz_poly_init(q);
        fmpz_poly_pseudo_divrem(q, r->num, &m, a->num, b->num);
        fmpz_poly_clear(q);
    }

    lead = fmpz_poly_lead(b->num);

    /* Case 1:  lead^m is 1 */
    if (fmpz_is_one(lead) || m == 0ul || (fmpz_is_m1(lead) & m % 2 == 0))
    {
        r->den = fmpz_realloc(r->den, fmpz_size(a->den));
        fmpz_set(r->den, a->den);
    }
    /* Case 2:  lead^m is -1 */
    else if (fmpz_is_m1(lead) & m % 2)
    {
        r->den = fmpz_realloc(r->den, fmpz_size(a->den));
        fmpz_neg(r->den, a->den);
    }
    /* Case 3:  lead^m is not +-1 */
    else
    {
        ulong limbs;
        limbs = m * fmpz_size(lead);

        if (fmpz_is_one(a->den))
        {
            r->den = fmpz_realloc(r->den, limbs);
            fmpz_pow_ui(r->den, lead, m);
        }
        else
        {
            fmpz_t t;
            t = fmpz_init(limbs);
            r->den = fmpz_realloc(r->den, limbs + fmpz_size(a->den));
            fmpz_pow_ui(t, lead, m);
            fmpz_mul(r->den, t, a->den);
            fmpz_clear(t);
        }
    }
    fmpq_poly_canonicalize(r, NULL);
}

/**
 * \ingroup  Division
 *
 * Sets \c q and \c r to the quotient and remainder of the Euclidean
 * division of \c a by \c b.
 *
 * Assumes that \c b is non-zero, and that \c q and \c r refer to distinct
 * objects in memory.
 */
void fmpq_poly_divrem(fmpq_poly_ptr q, fmpq_poly_ptr r, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    ulong m;
    fmpz_t lead;

    /* Assertion! */
    if (fmpq_poly_is_zero(b))
    {
        printf("ERROR (fmpq_poly_divrem).  Division by zero.\n");
        abort();
    }
    if (q == r)
    {
        printf("ERROR (fmpq_poly_divrem).  Output arguments aliased.\n");
        abort();
    }

    /* Catch the case when a and b are the same objects in memory.           */
    /* As of FLINT version 1.5.1, there is a bug in this case.               */
    if (a == b)
    {
        fmpq_poly_set_si(q, 1);
        fmpq_poly_set_si(r, 0);
        return;
    }

    /* Deal with the various other cases of aliasing.                        */
    if (r == a | r == b)
    {
        if (q == a | q == b)
        {
            fmpq_poly_t tempq, tempr;
            fmpq_poly_init(tempq);
            fmpq_poly_init(tempr);
            fmpq_poly_divrem(tempq, tempr, a, b);
            fmpq_poly_swap(q, tempq);
            fmpq_poly_swap(r, tempr);
            fmpq_poly_clear(tempq);
            fmpq_poly_clear(tempr);
            return;
        }
        else
        {
            fmpq_poly_t tempr;
            fmpq_poly_init(tempr);
            fmpq_poly_divrem(q, tempr, a, b);
            fmpq_poly_swap(r, tempr);
            fmpq_poly_clear(tempr);
            return;
        }
    }
    else
    {
        if (q == a | q == b)
        {
            fmpq_poly_t tempq;
            fmpq_poly_init(tempq);
            fmpq_poly_divrem(tempq, r, a, b);
            fmpq_poly_swap(q, tempq);
            fmpq_poly_clear(tempq);
            return;
        }
    }

    // TODO: Implement the case `\deg(b) = 0` more efficiently!

    /* General case..                                                        */
    /* Set q to b->den q->num / (a->den lead^m)                              */
    /* and r to r->num / (a->den lead^m).                                    */

    fmpz_poly_pseudo_divrem(q->num, r->num, &m, a->num, b->num);

    lead = fmpz_poly_lead(b->num);

    if (!fmpz_is_one(b->den))
        fmpz_poly_scalar_mul_fmpz(q->num, q->num, b->den);

    /* Case 1.  lead^m is 1 */
    if (fmpz_is_one(lead) || m == 0ul || (fmpz_is_m1(lead) & m % 2 == 0))
    {
        q->den = fmpz_realloc(q->den, fmpz_size(a->den));
        fmpz_set(q->den, a->den);
    }
    /* Case 2.  lead^m is -1 */
    else if (fmpz_is_m1(lead) & m % 2)
    {
        fmpz_poly_neg(q->num, q->num);
        q->den = fmpz_realloc(q->den, fmpz_size(a->den));
        fmpz_set(q->den, a->den);

        fmpz_poly_neg(r->num, r->num);
    }
    /* Case 3.  lead^m is not +-1 */
    else
    {
        ulong limbs;
        limbs = m * fmpz_size(lead);

        if (fmpz_is_one(a->den))
        {
            q->den = fmpz_realloc(q->den, limbs);
            fmpz_pow_ui(q->den, lead, m);
        }
        else
        {
            fmpz_t t;
            t = fmpz_init(limbs);
            q->den = fmpz_realloc(q->den, limbs + fmpz_size(a->den));
            fmpz_pow_ui(t, lead, m);
            fmpz_mul(q->den, t, a->den);
            fmpz_clear(t);
        }
    }
    r->den = fmpz_realloc(r->den, fmpz_size(q->den));
    fmpz_set(r->den, q->den);

    fmpq_poly_canonicalize(q, NULL);
    fmpq_poly_canonicalize(r, NULL);
}

///////////////////////////////////////////////////////////////////////////////
// Powering

/**
 * \ingroup  Powering
 *
 * Sets \c rop to the <tt>exp</tt>th power of \c op.
 *
 * The corner case of <tt>exp == 0</tt> is handled by setting \c rop to the
 * constant function \f$1\f$.  Note that this includes the case \f$0^0 = 1\f$.
 */
void fmpq_poly_power(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong exp)
{
    if (exp == 0ul)
    {
        fmpq_poly_one(rop);
    }
    else
    {
        fmpz_poly_power(rop->num, op->num, exp);
        if (fmpz_is_one(op->den))
        {
            fmpz_set_ui(rop->den, 1);
        }
        else
        {
            if (rop == op)
            {
                fmpz_t t;
                t = fmpz_init(exp * fmpz_size(op->den));
                fmpz_pow_ui(t, op->den, exp);
                fmpz_clear(rop->den);  // FIXME  Fixed
                rop->den = t;
            }
            else
            {
                rop->den = fmpz_realloc(rop->den, exp * fmpz_size(op->den));
                fmpz_pow_ui(rop->den, op->den, exp);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Greatest common divisor

/**
 * \ingroup  GCD
 *
 * Returns the (monic) greatest common divisor \c res of \c a and \c b.
 *
 * Corner cases:  If \c a and \c b are both zero, returns zero.  If only
 * one of them is zero, returns the other polynomial, up to normalisation.
 */
void fmpq_poly_gcd(fmpq_poly_ptr rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    fmpz_t lead;
    fmpz_poly_t num;

    /* Deal with aliasing */
    if (rop == a | rop == b)
    {
        fmpq_poly_t tpoly;
        fmpq_poly_init(tpoly);
        fmpq_poly_gcd(tpoly, a, b);
        fmpq_poly_swap(rop, tpoly);
        fmpq_poly_clear(tpoly);
        return;
    }

    /* Deal with corner cases */
    if (fmpq_poly_is_zero(a))
    {
        if (fmpq_poly_is_zero(b))  /* a == b == 0 */
            fmpq_poly_zero(rop);
        else                       /* a == 0, b != 0 */
            fmpq_poly_monic(rop, b);
        return;
    }
    else
    {
        if (fmpq_poly_is_zero(b))  /* a != 0, b == 0 */
        {
            fmpq_poly_monic(rop, a);
            return;
        }
    }

    /* General case.. */
    fmpz_poly_init(num);
    fmpz_poly_primitive_part(rop->num, a->num);
    fmpz_poly_primitive_part(num, b->num);

    /* Since rop.num is the greatest common divisor of the primitive parts   */
    /* of a.num and b.num, it is also primitive.  But as of FLINT 1.4.0, the */
    /* leading term *might* be negative.                                     */
    fmpz_poly_gcd(rop->num, rop->num, num);

    lead = fmpz_poly_lead(rop->num);
    rop->den = fmpz_realloc(rop->den, fmpz_size(lead));
    if (fmpz_sgn(lead) < 0)
    {
        fmpz_poly_neg(rop->num, rop->num);
        fmpz_neg(rop->den, lead);
    }
    else
        fmpz_set(rop->den, lead);

    fmpz_poly_clear(num);
}

/**
 * \ingroup  GCD
 *
 * Returns polynomials \c s, \c t and \c rop such that \c rop is
 * (monic) greatest common divisor of \c a and \c b, and such that
 * <tt>rop = s a + t b</tt>.
 *
 * Corner cases:  If \c a and \c b are zero, returns zero polynomials.
 * Otherwise, if only \c a is zero, returns <tt>(res, s, t) = (b, 0, 1)</tt>
 * up to normalisation, and similarly if only \c b is zero.
 *
 * Assumes that the output parameters \c rop, \c s, and \c t refer to
 * distinct objects in memory.  Otherwise, an exception is raised in the
 * form of an <tt>abort</tt> statement.
 */
void fmpq_poly_xgcd(fmpq_poly_ptr rop, fmpq_poly_ptr s, fmpq_poly_ptr t, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    fmpz_t lead, temp;
    ulong bound, limbs;

    /* Assertion! */
    if (rop == s | rop == t | s == t)
    {
        printf("ERROR (fmpq_poly_xgcd).  Output arguments aliased.\n");
        abort();
    }

    /* Deal with aliasing */
    if (rop == a | rop == b)
    {
        if (s == a | s == b)
        {
            /* We know that t does not coincide with a or b, since otherwise */
            /* two of rop, s, and t coincide, too.                           */
            fmpq_poly_t tempg, temps;
            fmpq_poly_init(tempg);
            fmpq_poly_init(temps);
            fmpq_poly_xgcd(tempg, temps, t, a, b);
            fmpq_poly_swap(rop, tempg);
            fmpq_poly_swap(s, temps);
            fmpq_poly_clear(tempg);
            fmpq_poly_clear(temps);
            return;
        }
        else
        {
            if (t == a | t == b)
            {
                fmpq_poly_t tempg, tempt;
                fmpq_poly_init(tempg);
                fmpq_poly_init(tempt);
                fmpq_poly_xgcd(tempg, s, tempt, a, b);
                fmpq_poly_swap(rop, tempg);
                fmpq_poly_swap(t, tempt);
                fmpq_poly_clear(tempg);
                fmpq_poly_clear(tempt);
                return;
            }
            else
            {
                fmpq_poly_t tempg;
                fmpq_poly_init(tempg);
                fmpq_poly_xgcd(tempg, s, t, a, b);
                fmpq_poly_swap(rop, tempg);
                fmpq_poly_clear(tempg);
                return;
            }
        }
    }
    else
    {
        if (s == a | s == b)
        {
            if (t == a | t == b)
            {
                fmpq_poly_t temps, tempt;
                fmpq_poly_init(temps);
                fmpq_poly_init(tempt);
                fmpq_poly_xgcd(rop, temps, tempt, a, b);
                fmpq_poly_swap(s, temps);
                fmpq_poly_swap(t, tempt);
                fmpq_poly_clear(temps);
                fmpq_poly_clear(tempt);
                return;
            }
            else
            {
                fmpq_poly_t temps;
                fmpq_poly_init(temps);
                fmpq_poly_xgcd(rop, temps, t, a, b);
                fmpq_poly_swap(s, temps);
                fmpq_poly_clear(temps);
                return;
            }
        }
        else
        {
            if (t == a | t == b)
            {
                fmpq_poly_t tempt;
                fmpq_poly_init(tempt);
                fmpq_poly_xgcd(rop, s, tempt, a, b);
                fmpq_poly_swap(t, tempt);
                fmpq_poly_clear(tempt);
                return;
            }
        }
    }

    /* From here on, we may assume that none of the output variables are */
    /* aliases for the input variables.                                  */

    /* Deal with the following three corner cases: */
    /*   a == 0, b == 0                            */
    /*   a == 0, b =! 0                            */
    /*   a != 0, b == 0                            */
    if (fmpq_poly_is_zero(a))
    {
        if (fmpq_poly_is_zero(b))  /* Case 1.  a == b == 0 */
        {
            fmpq_poly_zero(rop);
            fmpq_poly_zero(s);
            fmpq_poly_zero(t);
            return;
        }
        else  /* Case 2.  a == 0, b != 0 */
        {
            fmpz_t blead = fmpz_poly_lead(b->num);

            fmpq_poly_monic(rop, b);
            fmpq_poly_zero(s);

            fmpz_poly_zero(t->num);
            fmpz_poly_set_coeff_fmpz(t->num, 0, b->den);
            if (fmpz_is_one(blead))
                fmpz_set_ui(t->den, 1);
            else
            {
                fmpz_t g;
                g = fmpz_init(FLINT_MAX(b->num->limbs, fmpz_size(b->den)));
                fmpz_gcd(g, blead, b->den);
                if (fmpz_sgn(g) != fmpz_sgn(blead))
                    fmpz_neg(g, g);
                fmpz_poly_scalar_div_fmpz(t->num, t->num, g);
                t->den = fmpz_realloc(t->den, fmpz_size(blead));
                fmpz_divides(t->den, blead, g);
                fmpz_clear(g);
            }
            return;
        }
    }
    else
    {
        if (fmpq_poly_is_zero(b))  /* Case 3.  a != 0, b == 0 */
        {
            fmpq_poly_xgcd(rop, t, s, b, a);
            return;
        }
    }

    /* We are now in the general case where a and b are non-zero. */

    s->den = fmpz_realloc(s->den, a->num->limbs);
    t->den = fmpz_realloc(t->den, b->num->limbs);
    fmpz_poly_content(s->den, a->num);
    fmpz_poly_content(t->den, b->num);
    fmpz_poly_scalar_div_fmpz(s->num, a->num, s->den);
    fmpz_poly_scalar_div_fmpz(t->num, b->num, t->den);

    /* Note that, since s->num and t->num are primitive, rop->num is         */
    /* primitive, too.  In fact, it is the rational greatest common divisor  */
    /* of a and b.  As of FLINT 1.4.0, the leading coefficient of res.num    */
    /* *might* be negative.                                                  */

    fmpz_poly_gcd(rop->num, s->num, t->num);
    if (fmpz_sgn(fmpz_poly_lead(rop->num)) < 0)
        fmpz_poly_neg(rop->num, rop->num);
    lead = fmpz_poly_lead(rop->num);

    /* Now rop->num is a (primitive) rational greatest common divisor of */
    /* a and b.                                                          */

    if (fmpz_poly_degree(rop->num) > 0)
    {
        fmpz_poly_div(s->num, s->num, rop->num);
        fmpz_poly_div(t->num, t->num, rop->num);
    }

    bound = fmpz_poly_resultant_bound(s->num, t->num);
    if (fmpz_is_one(lead))
        bound = bound / FLINT_BITS + 2;
    else
        bound = FLINT_MAX(bound / FLINT_BITS + 2, fmpz_size(lead));
    rop->den = fmpz_realloc(rop->den, bound);

    fmpz_poly_xgcd(rop->den, s->num, t->num, s->num, t->num);

    /* Now the following equation holds:                                     */
    /*   rop->den rop->num ==                                                */
    /*             (s->num a->den / s->den) a +  (t->num b->den / t->den) b. */

    limbs = FLINT_MAX(s->num->limbs, t->num->limbs);
    limbs = FLINT_MAX(limbs, fmpz_size(s->den));
    limbs = FLINT_MAX(limbs, fmpz_size(t->den) + fmpz_size(rop->den) + fmpz_size(lead));
    temp = fmpz_init(limbs);

    s->den = fmpz_realloc(s->den, fmpz_size(s->den) + fmpz_size(rop->den)
                                                    + fmpz_size(lead));
    if (!fmpz_is_one(a->den))
        fmpz_poly_scalar_mul_fmpz(s->num, s->num, a->den);
    fmpz_mul(temp, s->den, rop->den);
    fmpz_mul(s->den, temp, lead);

    t->den = fmpz_realloc(t->den, fmpz_size(t->den) + fmpz_size(rop->den)
                                                    + fmpz_size(lead));
    if (!fmpz_is_one(b->den))
        fmpz_poly_scalar_mul_fmpz(t->num, t->num, b->den);
    fmpz_mul(temp, t->den, rop->den);
    fmpz_mul(t->den, temp, lead);

    fmpq_poly_canonicalize(s, temp);
    fmpq_poly_canonicalize(t, temp);

    fmpz_set(rop->den, lead);

    fmpz_clear(temp);
}

/**
 * \ingroup  GCD
 *
 * Computes the monic (or zero) least common multiple of \c a and \c b.
 *
 * If either of \c a and \c b is zero, returns zero.  This behaviour ensures
 * that the relation
 * \f[
 * \text{lcm}(a,b) \gcd(a,b) \sim a b
 * \f]
 * holds, where \f$\sim\f$ denotes equality up to units.
 */
void fmpq_poly_lcm(fmpq_poly_ptr rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    /* Handle aliasing */
    if (rop == a | rop == b)
    {
        fmpq_poly_t tpoly;
        fmpq_poly_init(tpoly);
        fmpq_poly_lcm(tpoly, a, b);
        fmpq_poly_swap(rop, tpoly);
        fmpq_poly_clear(tpoly);
        return;
    }

    if (fmpq_poly_is_zero(a))
        fmpq_poly_zero(rop);
    else if (fmpq_poly_is_zero(b))
        fmpq_poly_zero(rop);
    else
    {
        fmpz_poly_t prod;
        fmpz_t lead;

        fmpz_poly_init(prod);
        fmpq_poly_gcd(rop, a, b);
        fmpz_poly_mul(prod, a->num, b->num);
        fmpz_poly_primitive_part(prod, prod);
        fmpz_poly_div(rop->num, prod, rop->num);

        /* Note that GCD returns a monic rop and so a primitive rop.num.   */
        /* Dividing the primitive prod by this yields a primitive quotient */
        /* rop->num.                                                       */

        lead = fmpz_poly_lead(rop->num);
        rop->den = fmpz_realloc(rop->den, fmpz_size(lead));
        if (fmpz_sgn(lead) < 0)
            fmpz_poly_neg(rop->num, rop->num);
        fmpz_set(rop->den, fmpz_poly_lead(rop->num));

        fmpz_poly_clear(prod);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Derivative

/**
 * \ingroup  Derivative
 *
 * Sets \c rop to the derivative of \c op.
 *
 * \todo  The second argument should be declared \c const, but as of
 *        FLINT 1.5.0 this generates a compile-time warning.
 */
void fmpq_poly_derivative(fmpq_poly_ptr rop, fmpq_poly_ptr op)
{
    if (fmpq_poly_length(op) < 2ul)
        fmpq_poly_zero(rop);
    else
    {
        fmpz_poly_derivative(rop->num, op->num);
        rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
        fmpz_set(rop->den, op->den);
        fmpq_poly_canonicalize(rop, NULL);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Evaluation

/**
 * \ingroup  Evaluation
 *
 * Evaluates the integer polynomial \c f at the rational \c a using Horner's
 * method.
 */
void _fmpz_poly_evaluate_mpq_horner(mpq_t rop, const fmpz_poly_t f, const mpq_t a)
{
    mpq_t temp;
    ulong n;

    /* Handle aliasing */
    if (rop == a)
    {
        mpq_t tempr;
        mpq_init(tempr);
        _fmpz_poly_evaluate_mpq_horner(tempr, f, a);
        mpq_swap(rop, tempr);
        mpq_clear(tempr);
        return;
    }

    n = fmpz_poly_length(f);

    if (n == 0ul)
    {
        mpq_set_ui(rop, 0, 1);
    }
    else if (n == 1ul)
    {
        fmpz_poly_get_coeff_mpz(mpq_numref(rop), f, 0);
        mpz_set_ui(mpq_denref(rop), 1);
    }
    else
    {
        n--;
        fmpz_poly_get_coeff_mpz(mpq_numref(rop), f, n);
        mpz_set_ui(mpq_denref(rop), 1);
        mpq_init(temp);
        do {
            n--;
            mpq_mul(temp, rop, a);
            fmpz_poly_get_coeff_mpz(mpq_numref(rop), f, n);
            mpz_set_ui(mpq_denref(rop), 1);
            mpq_add(rop, rop, temp);
        } while (n);
        mpq_clear(temp);
    }
}

/**
 * \ingroup  Evaluation
 *
 * Evaluates the rational polynomial \c f at the integer \c a.
 *
 * Assumes that the numerator and denominator of the <tt>mpq_t</tt>
 * \c rop are distinct (as objects in memory) from the <tt>mpz_t</tt>
 * \c a.
 */
void fmpq_poly_evaluate_mpz(mpq_t rop, fmpq_poly_ptr f, const mpz_t a)
{
    fmpz_t num, t;
    ulong limbs, max;

    if (fmpq_poly_is_zero(f))
    {
        mpq_set_ui(rop, 0, 1);
        return;
    }

    /* Establish a bound on the size of f->num evaluated at a */
    max = (f->num->length) * (f->num->limbs + f->num->length * mpz_size(a));

    /* Compute the result */
    num = fmpz_init(max);
    t   = fmpz_init(mpz_size(a));
    mpz_to_fmpz(t, a);
    fmpz_poly_evaluate(num, f->num, t);
    fmpz_to_mpz(mpq_numref(rop), num);
    if (fmpz_is_one(f->den))
    {
        mpz_set_ui(mpq_denref(rop), 1);
    }
    else
    {
        fmpz_to_mpz(mpq_denref(rop), f->den);
        mpq_canonicalize(rop);
    }

    fmpz_clear(num);
    fmpz_clear(t);
}

/**
 * \ingroup  Evaluation
 *
 * Evaluates the rational polynomial \c f at the rational \c a.
 */
void fmpq_poly_evaluate_mpq(mpq_t rop, fmpq_poly_ptr f, const mpq_t a)
{
    if (rop == a)
    {
        mpq_t tempr;
        mpq_init(tempr);
        fmpq_poly_evaluate_mpq(tempr, f, a);
        mpq_swap(rop, tempr);
        mpq_clear(tempr);
        return;
    }

    _fmpz_poly_evaluate_mpq_horner(rop, f->num, a);
    if (!fmpz_is_one(f->den))
    {
        mpz_t den;
        mpz_init(den);
        fmpz_to_mpz(den, f->den);
        mpz_mul(mpq_denref(rop), mpq_denref(rop), den);
        mpq_canonicalize(rop);
        mpz_clear(den);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Gaussian content

/**
 * \ingroup  Content
 *
 * Returns the non-negative content of \c op.
 *
 * The content of \f$0\f$ is defined to be \f$0\f$.
 */
void fmpq_poly_content(mpq_t rop, const fmpq_poly_ptr op)
{
    fmpz_t numc;

    numc = fmpz_init(op->num->limbs);
    fmpz_poly_content(numc, op->num);
    fmpz_abs(numc, numc);
    fmpz_to_mpz(mpq_numref(rop), numc);
    if (fmpz_is_one(op->den))
        mpz_set_ui(mpq_denref(rop), 1);
    else
        fmpz_to_mpz(mpq_denref(rop), op->den);
    fmpz_clear(numc);
}

/**
 * \ingroup  Content
 *
 * Returns the primitive part (with non-negative leading coefficient) of
 * \c op as an element of type #fmpq_poly_t.
 */
void fmpq_poly_primitive_part(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    if (fmpq_poly_is_zero(op))
        fmpq_poly_zero(rop);
    else
    {
        fmpz_poly_primitive_part(rop->num, op->num);
        if (fmpz_sgn(fmpz_poly_lead(rop->num)) < 0)
            fmpz_poly_neg(rop->num, rop->num);
        fmpz_set_ui(rop->den, 1);
    }
}

/**
 * \brief    Returns whether \c op is monic.
 * \ingroup  Content
 *
 * Returns whether \c op is monic.
 *
 * By definition, the zero polynomial is \e not monic.
 */
int fmpq_poly_is_monic(const fmpq_poly_ptr op)
{
    fmpz_t lead;

    if (fmpq_poly_is_zero(op))
        return 0;
    else
    {
        lead = fmpz_poly_lead(op->num);
        if (fmpz_is_one(op->den))
            return fmpz_is_one(lead);
        else
            return (fmpz_sgn(lead) > 0) && (fmpz_cmpabs(lead, op->den) == 0);
    }
}

/**
 * Sets \c rop to the unique monic scalar multiple of \c op.
 *
 * As the only special case, if \c op is the zero polynomial, \c rop is set
 * to zero, too.
 */
void fmpq_poly_monic(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    fmpz_t lead;

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpz_poly_primitive_part(rop->num, op->num);
    lead = fmpz_poly_lead(rop->num);
    rop->den = fmpz_realloc(rop->den, fmpz_size(lead));
    if (fmpz_sgn(lead) < 0)
        fmpz_poly_neg(rop->num, rop->num);
    fmpz_set(rop->den, fmpz_poly_lead(rop->num));
}

///////////////////////////////////////////////////////////////////////////////
// Resultant

/**
 * \brief    Returns the resultant of \c a and \c b.
 * \ingroup  Resultant
 *
 * Returns the resultant of \c a and \c b.
 *
 * Enumerating the roots of \c a and \c b over \f$\bar{\mathbf{Q}}\f$ as
 * \f$r_1, \dotsc, r_m\f$ and \f$s_1, \dotsc, s_n\f$, respectively, and
 * letting \f$x\f$ and \f$y\f$ denote the leading coefficients, the resultant
 * is defined as
 * \f[
 * x^{\deg(b)} y^{\deg(a)} \prod_{1 \leq i, j \leq n} (r_i - s_j).
 * \f]
 *
 * We handle special cases as follows:  if one of the polynomials is zero,
 * the resultant is zero.  Note that otherwise if one of the polynomials is
 * constant, the last term in the above expression is the empty product.
 */
void fmpq_poly_resultant(mpq_t rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    fmpz_t rest, t1, t2;
    fmpz_poly_t anum, bnum;
    fmpz_t acont, bcont;
    long d1, d2;
    ulong bound, denbound, numbound;
    fmpz_poly_t g;

    d1 = fmpq_poly_degree(a);
    d2 = fmpq_poly_degree(b);

    /* We first handle special cases. */
    if (d1 < 0l | d2 < 0l)
    {
        mpq_set_ui(rop, 0, 1);
        return;
    }
    if (d1 == 0l)
    {
        if (d2 == 0l)
            mpq_set_ui(rop, 1, 1);
        else if (d2 == 1l)
        {
            fmpz_to_mpz(mpq_numref(rop), fmpz_poly_lead(a->num));
            fmpz_to_mpz(mpq_denref(rop), a->den);
        }
        else
        {
            if (fmpz_is_one(a->den))
                bound = a->num->limbs;
            else
                bound = FLINT_MAX(a->num->limbs, fmpz_size(a->den));
            t1 = fmpz_init(d2 * bound);
            fmpz_pow_ui(t1, fmpz_poly_lead(a->num), d2);
            fmpz_to_mpz(mpq_numref(rop), t1);
            if (fmpz_is_one(a->den))
                mpz_set_ui(mpq_denref(rop), 1);
            else
            {
                fmpz_pow_ui(t1, a->den, d2);
                fmpz_to_mpz(mpq_denref(rop), t1);
            }
            fmpz_clear(t1);
        }
        return;
    }
    if (d2 == 0l)
    {
        fmpq_poly_resultant(rop, b, a);
        return;
    }

    /* We are now in the general case, with both polynomials of degree at */
    /* least 1.                                                           */

    /* We set a->num =: acont anum with acont > 0 and anum primitive. */
    acont = fmpz_init(a->num->limbs);
    fmpz_poly_content(acont, a->num);
    fmpz_abs(acont, acont);
    fmpz_poly_init(anum);
    fmpz_poly_scalar_div_fmpz(anum, a->num, acont);

    /* We set b->num =: bcont bnum with bcont > 0 and bnum primitive. */
    bcont = fmpz_init(b->num->limbs);
    fmpz_poly_content(bcont, b->num);
    fmpz_abs(bcont, bcont);
    fmpz_poly_init(bnum);
    fmpz_poly_scalar_div_fmpz(bnum, b->num, bcont);

    fmpz_poly_init(g);
    fmpz_poly_gcd(g, anum, bnum);

    /* If the gcd has positive degree, the resultant is zero. */
    if (fmpz_poly_degree(g) > 0)
    {
        mpq_set_ui(rop, 0, 1);

        /* Clean up */
        fmpz_clear(acont);
        fmpz_clear(bcont);
        fmpz_poly_clear(anum);
        fmpz_poly_clear(bnum);
        fmpz_poly_clear(g);
        return;
    }

    /* Set some bounds */
    if (fmpz_is_one(a->den))
    {
        if (fmpz_is_one(b->den))
        {
            numbound = FLINT_MAX(d2 * fmpz_size(acont), d1 * fmpz_size(bcont));
            denbound = 1;
        }
        else
        {
            numbound = FLINT_MAX(d2 * fmpz_size(acont), d1 * fmpz_size(bcont));
            denbound = d1 * fmpz_size(b->den);
        }
    }
    else
    {
        if (fmpz_is_one(b->den))
        {
            numbound = FLINT_MAX(d2 * fmpz_size(acont), d1 * fmpz_size(bcont));
            denbound = d2 * fmpz_size(a->den);
        }
        else
        {
            numbound = FLINT_MAX(d2 * fmpz_size(acont), d1 * fmpz_size(bcont));
            denbound = FLINT_MAX(d2 * fmpz_size(a->den), d1 * fmpz_size(b->den));
        }
    }

    /* Now anum and bnum are coprime and we compute their resultant using */
    /* the method from the fmpz_poly module.                              */
    bound = fmpz_poly_resultant_bound(anum, bnum);
    bound = bound/FLINT_BITS + 2 + d2 * fmpz_size(acont)
        + d1 * fmpz_size(bcont);
    rest = fmpz_init(FLINT_MAX(bound, denbound));
    fmpz_poly_resultant(rest, anum, bnum);

    /* Finally, w take into account the factors acont/a.den and bcont/b.den. */
    t1 = fmpz_init(FLINT_MAX(bound, denbound));
    t2 = fmpz_init(FLINT_MAX(numbound, denbound));

    if (!fmpz_is_one(acont))
    {
        fmpz_pow_ui(t2, acont, d2);
        fmpz_set(t1, rest);
        fmpz_mul(rest, t1, t2);
    }
    if (!fmpz_is_one(bcont))
    {
        fmpz_pow_ui(t2, bcont, d1);
        fmpz_set(t1, rest);
        fmpz_mul(rest, t1, t2);
    }

    fmpz_to_mpz(mpq_numref(rop), rest);

    if (fmpz_is_one(a->den))
    {
        if (fmpz_is_one(b->den))
            fmpz_set_ui(rest, 1);
        else
            fmpz_pow_ui(rest, b->den, d1);
    }
    else
    {
        if (fmpz_is_one(b->den))
            fmpz_pow_ui(rest, a->den, d2);
        else
        {
            fmpz_pow_ui(t1, a->den, d2);
            fmpz_pow_ui(t2, b->den, d1);
            fmpz_mul(rest, t1, t2);
        }
    }

    fmpz_to_mpz(mpq_denref(rop), rest);
    mpq_canonicalize(rop);

    /* Clean up */
    fmpz_clear(rest);
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_poly_clear(anum);
    fmpz_poly_clear(bnum);
    fmpz_clear(acont);
    fmpz_clear(bcont);
    fmpz_poly_clear(g);
}

/**
 * \brief    Computes the discriminant of \c a.
 * \ingroup  Discriminant
 *
 * Computes the discriminant of the polynomial \f$a\f$.
 *
 * The discriminant \f$R_n\f$ is defined as
 * \f[
 * R_n = a_n^{2 n-2} \prod_{1 \le i < j \le n} (r_i - r_j)^2,
 * \f]
 * where \f$n\f$ is the degree of this polynomial, \f$a_n\f$ is the leading
 * coefficient and \f$r_1, \ldots, r_n\f$ are the roots over
 * \f$\bar{\mathbf{Q}}\f$ are.
 *
 * The discriminant of constant polynomials is defined to be \f$0\f$.
 *
 * This implementation uses the identity
 * \f[
 * R_n(f) := (-1)^(n (n-1)/2) R(f,f') / a_n,
 * \f]
 * where \f$n\f$ is the degree of this polynomial, \f$a_n\f$ is the leading
 * coefficient and \f$f'\f$ is the derivative of \f$f\f$.
 *
 * \see  #fmpq_poly_resultant()
 */
void fmpq_poly_discriminant(mpq_t disc, fmpq_poly_t a)
{
    fmpq_poly_t der;
    mpq_t t;
    long n;

    n = fmpq_poly_degree(a);
    if (n < 1L)
    {
        mpq_set_ui(disc, 0, 1);
        return;
    }

    fmpq_poly_init(der);
    fmpq_poly_derivative(der, a);
    fmpq_poly_resultant(disc, a, der);
    mpq_init(t);

    fmpz_to_mpz(mpq_numref(t), a->den);
    fmpz_to_mpz(mpq_denref(t), fmpz_poly_lead(a->num));
    mpq_canonicalize(t);
    mpq_mul(disc, disc, t);

    if (n % 4 > 1)
        mpz_neg(mpq_numref(disc), mpq_numref(disc));

    fmpq_poly_clear(der);
    mpq_clear(t);
}

///////////////////////////////////////////////////////////////////////////////
// Composition

/**
 * \brief    Returns the composition of \c a and \c b.
 * \ingroup  Composition
 *
 * Returns the composition of \c a and \c b.
 *
 * To be clear about the order of the composition, denoting the polynomials
 * \f$a(T)\f$ and \f$b(T)\f$, returns the polynomial \f$a(b(T))\f$.
 */
void fmpq_poly_compose(fmpq_poly_ptr rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b)
{
    mpq_t x;            /* Rational to hold the inverse of b->den */

    if (fmpz_is_one(b->den))
    {
        fmpz_poly_compose(rop->num, a->num, b->num);
        rop->den = fmpz_realloc(rop->den, fmpz_size(a->den));
        fmpz_set(rop->den, a->den);
        fmpq_poly_canonicalize(rop, NULL);
        return;
    }

    /* Aliasing.                                                            */
    /*                                                                      */
    /* Note that rop and a, as well as a and b may be aliased, but rop and  */
    /* b may not be aliased.                                                */
    if (rop == b)
    {
        fmpq_poly_t tempr;
        fmpq_poly_init(tempr);
        fmpq_poly_compose(tempr, a, b);
        fmpq_poly_swap(rop, tempr);
        fmpq_poly_clear(tempr);
        return;
    }

    /* Set x = 1/b.den, and note this is already in canonical form. */
    mpq_init(x);
    mpz_set_ui(mpq_numref(x), 1);
    fmpz_to_mpz(mpq_denref(x), b->den);

    /* First set rop = a(T / b.den) and then use FLINT's composition to */
    /* set rop->num = rop->num(b->num).                                 */
    fmpq_poly_rescale(rop, a, x);
    fmpz_poly_compose(rop->num, rop->num, b->num);
    fmpq_poly_canonicalize(rop, NULL);

    mpq_clear(x);
}

/**
 * \brief    Rescales the co-ordinate.
 * \ingroup  Composition
 *
 * Denoting this polynomial \f$a(T)\f$, computes the polynomial \f$a(x T)\f$.
 */
void fmpq_poly_rescale(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpq_t x)
{
    ulong numsize, densize;
    ulong limbs;
    fmpz_t num, den;
    fmpz_t coeff, power, t;
    long i, n;

    fmpq_poly_set(rop, op);

    if (fmpq_poly_length(rop) < 2ul)
        return;

    num = fmpz_init(mpz_size(mpq_numref(x)));
    den = fmpz_init(mpz_size(mpq_denref(x)));
    mpz_to_fmpz(num, mpq_numref(x));
    mpz_to_fmpz(den, mpq_denref(x));
    numsize = fmpz_size(num);
    densize = fmpz_size(den);

    n = fmpz_poly_degree(rop->num);

    if (fmpz_is_one(den))
    {
        coeff = fmpz_init(rop->num->limbs + n * numsize);
        power = fmpz_init(n * numsize);
        t = fmpz_init(rop->num->limbs + n * numsize);

        fmpz_set(power, num);

        fmpz_poly_get_coeff_fmpz(t, rop->num, 1);
        fmpz_mul(coeff, t, power);
        fmpz_poly_set_coeff_fmpz(rop->num, 1, coeff);

        for (i = 2; i <= n; i++)
        {
            fmpz_set(t, power);
            fmpz_mul(power, t, num);
            fmpz_poly_get_coeff_fmpz(t, rop->num, i);
            fmpz_mul(coeff, t, power);
            fmpz_poly_set_coeff_fmpz(rop->num, i, coeff);
        }
    }
    else
    {
        coeff = fmpz_init(rop->num->limbs + n * (numsize + densize));
        power = fmpz_init(n * (numsize + densize));
        limbs = rop->num->limbs + n * (numsize + densize);
        limbs = FLINT_MAX(limbs, fmpz_size(rop->den));
        t = fmpz_init(limbs);

        fmpz_pow_ui(power, den, n);

        if (fmpz_is_one(rop->den))
        {
            rop->den = fmpz_realloc(rop->den, n * densize);
            fmpz_set(rop->den, power);
        }
        else
        {
            fmpz_set(t, rop->den);
            limbs = n * densize + fmpz_size(rop->den);
            rop->den = fmpz_realloc(rop->den, limbs);
            fmpz_mul(rop->den, power, t);
        }

        fmpz_set_ui(power, 1);
        for (i = n - 1; i >= 0; i--)
        {
            fmpz_set(t, power);
            fmpz_mul(power, t, den);
            fmpz_poly_get_coeff_fmpz(t, rop->num, i);
            fmpz_mul(coeff, t, power);
            fmpz_poly_set_coeff_fmpz(rop->num, i, coeff);
        }

        fmpz_set_ui(power, 1);
        for (i = 1; i <= n; i++)
        {
            fmpz_set(t, power);
            fmpz_mul(power, t, num);
            fmpz_poly_get_coeff_fmpz(t, rop->num, i);
            fmpz_mul(coeff, t, power);
            fmpz_poly_set_coeff_fmpz(rop->num, i, coeff);
        }
    }
    fmpq_poly_canonicalize(rop, NULL);
    fmpz_clear(num);
    fmpz_clear(den);
    fmpz_clear(coeff);
    fmpz_clear(power);
    fmpz_clear(t);
}

///////////////////////////////////////////////////////////////////////////////
// Square-free

/**
 * \brief    Returns whether \c op is squarefree.
 * \ingroup  Squarefree
 *
 * Returns whether \c op is squarefree.
 *
 * By definition, a polynomial is square-free if it is not a multiple of a
 * non-unit factor.  In particular, polynomials up to degree 1 (including)
 * are square-free.
 */
int fmpq_poly_is_squarefree(const fmpq_poly_ptr op)
{
    if (fmpq_poly_length(op) < 3ul)
        return 1;
    else
    {
        int ans;
        fmpz_poly_t prim;

        fmpz_poly_init(prim);
        fmpz_poly_primitive_part(prim, op->num);
        ans = fmpz_poly_is_squarefree(prim);
        fmpz_poly_clear(prim);
        return ans;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Subpolynomials

/**
 * \brief    Returns a contiguous subpolynomial.
 * \ingroup  Subpolynomials
 *
 * Returns the slice with coefficients from \f$x^i\f$ (including) to
 * \f$x^j\f$ (excluding).
 */
void fmpq_poly_getslice(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong i, ulong j)
{
    ulong k;

    j = FLINT_MIN(fmpq_poly_length(op), j);

    /* Aliasing */
    if (rop == op)
    {
        for (k = 0; k < i; k++)
            fmpz_poly_set_coeff_ui(rop->num, k, 0);
        for (k = j; k < fmpz_poly_length(rop->num); k++)
            fmpz_poly_set_coeff_ui(rop->num, k, 0);
        fmpq_poly_canonicalize(rop, NULL);
        return;
    }

    fmpz_poly_zero(rop->num);
    if (i < j)
    {
        for (k = j; k != i; )
        {
            k--;
            fmpz_poly_set_coeff_fmpz(rop->num, k,
                                     fmpz_poly_get_coeff_ptr(op->num, k));
        }
        rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
        fmpz_set(rop->den, op->den);
        fmpq_poly_canonicalize(rop, NULL);
    }
    else
    {
        fmpz_set_ui(rop->den, 1);
    }
}

/**
 * \brief    Shifts this polynomial to the left by \f$n\f$ coefficients.
 * \ingroup  Subpolynomials
 *
 * Notionally multiplies \c op by \f$t^n\f$ and stores the result in \c rop.
 */
void fmpq_poly_left_shift(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n)
{
    /* XXX:  As a workaround for a bug in FLINT 1.5.1, we need to handle */
    /* the zero polynomial separately.                                   */
    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    if (n == 0ul)
    {
        fmpq_poly_set(rop, op);
        return;
    }

    fmpz_poly_left_shift(rop->num, op->num, n);
    if (rop != op)
    {
        rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
        fmpz_set(rop->den, op->den);
    }
}

/**
 * \brief    Shifts this polynomial to the right by \f$n\f$ coefficients.
 * \ingroup  Subpolynomials
 *
 * Notationally returns the quotient of floor division of \c rop by \c op.
 */
void fmpq_poly_right_shift(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n)
{
    /* XXX:  As a workaround for a bug in FLINT 1.5.1, we need to handle */
    /* the zero polynomial separately.                                   */
    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }

    if (n == 0ul)
    {
        fmpq_poly_set(rop, op);
        return;
    }

    fmpz_poly_right_shift(rop->num, op->num, n);
    if (rop != op)
    {
        rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
        fmpz_set(rop->den, op->den);
    }
    fmpq_poly_canonicalize(rop, NULL);
}

/**
 * \brief    Truncates this polynomials.
 * \ingroup  Subpolynomials
 *
 * Truncates <tt>op</tt>  modulo \f$x^n\f$ whenever \f$n\f$ is positive.
 * Returns zero otherwise.
 */
void fmpq_poly_truncate(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n)
{
    fmpq_poly_set(rop, op);
    if (fmpq_poly_length(rop) > n)
    {
        fmpz_poly_truncate(rop->num, n);
        fmpq_poly_canonicalize(rop, NULL);
    }
}

/**
 * \brief    Reverses this polynomial.
 * \ingroup  Subpolynomials
 *
 * Reverses the coefficients of \c op - thought of as a polynomial of
 * length \c n - and places the result in <tt>rop</tt>.
 */
void fmpq_poly_reverse(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n)
{
    ulong len;
    len = fmpq_poly_length(op);

    if (n == 0ul || len == 0ul)
    {
        fmpq_poly_zero(rop);
        return;
    }

    fmpz_poly_reverse(rop->num, op->num, n);
    if (rop != op)
    {
        rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
        fmpz_set(rop->den, op->den);
    }

    if (n < len)
        fmpq_poly_canonicalize(rop, NULL);
}

///////////////////////////////////////////////////////////////////////////////
// String conversion

/**
 * \addtogroup StringConversions
 *
 * The following three methods enable users to construct elements of type
 * \c fmpq_poly_ptr from strings or to obtain string representations of such
 * elements.
 *
 * The format used is based on the FLINT format for integer polynomials of
 * type <tt>fmpz_poly_t</tt>, which we recall first:
 *
 * A non-zero polynomial \f$a_0 + a_1 X + \dotsb + a_n X^n\f$ of length
 * \f$n + 1\f$ is represented by the string <tt>n+1  a_0 a_1 ... a_n</tt>,
 * where there are two space characters following the length and single space
 * characters separating the individual coefficients.  There is no leading or
 * trailing white-space.  In contrast, the zero polynomial is represented by
 * <tt>0</tt>.
 *
 * We adapt this notation for rational polynomials by using the <tt>mpq_t</tt>
 * notation for the coefficients without any additional white-space.
 *
 * There is also a <tt>_pretty</tt> variant available.
 *
 * Note that currently these functions are not optimized for performance and
 * are intended to be used only for debugging purposes or one-off input and
 * output, rather than as a low-level parser.
 */

/**
 * \ingroup StringConversions
 * \brief Constructs a polynomial from a list of rationals.
 *
 * Given a list of length <tt>n</tt> containing rationals of type
 * <tt>mpq_t</tt>, this method constructs a polynomial with these
 * coefficients, beginning with the constant term.
 */
void _fmpq_poly_from_list(fmpq_poly_ptr rop, mpq_t * a, ulong n)
{
    mpz_t den, t;
    ulong i;

    mpz_init(t);
    mpz_init_set_ui(den, 1);
    for (i = 0; i < n; i++)
        mpz_lcm(den, den, mpq_denref(a[i]));

    for (i = 0; i < n; i++)
    {
        mpz_divexact(t, den, mpq_denref(a[i]));
        mpz_mul(mpq_numref(a[i]), mpq_numref(a[i]), t);
    }

    fmpz_poly_realloc(rop->num, n);
    for (i = n - 1ul; i != -1ul; i--)
        fmpz_poly_set_coeff_mpz(rop->num, i, mpq_numref(a[i]));

    if (mpz_cmp_ui(den, 1) == 0)
        fmpz_set_ui(rop->den, 1);
    else
    {
        rop->den = fmpz_realloc(rop->den, mpz_size(den));
        mpz_to_fmpz(rop->den, den);
    }
    mpz_clear(den);
    mpz_clear(t);
}

/**
 * \ingroup  StringConversions
 *
 * Sets the rational polynomial \c rop to the value specified by the
 * null-terminated string \c str.
 *
 * The behaviour is undefined if the format of the string \c str does not
 * conform to the specification.
 *
 * \todo  In a future version, it would be nice to have this function return
 *        <tt>0</tt> or <tt>1</tt> depending on whether the input format
 *        was correct.  Currently, the method always returns <tt>1</tt>.
 */
int fmpq_poly_from_string(fmpq_poly_ptr rop, const char * str)
{
    mpq_t * a;
    char * strcopy;
    ulong strcopy_len;
    ulong i, j, k, n;

    n = atoi(str);

    /* Handle the case of the zero polynomial */
    if (n == 0ul)
    {
        fmpq_poly_zero(rop);
        return 1;
    }

    /* Compute the maximum length that the copy buffer has to be */
    strcopy_len = 0;
    for (j = 0; str[j] != ' '; j++);
    j += 2;
    for (i = 0; i < n; i++)
    {
        for (k = j; !(str[k] == ' ' || str[k] == '\0'); k++);
        strcopy_len = FLINT_MAX(strcopy_len, k - j + 1);
        j = k + 1;
    }

    strcopy = (char *) malloc(strcopy_len * sizeof(char));
    if (!strcopy)
    {
        printf("ERROR (fmpq_poly_from_string).  Memory allocation failed.\n");
        abort();
    }

    /* Read the data into the array a of mpq_t's */
    a = (mpq_t *) malloc(n * sizeof(mpq_t));
    for (j = 0; str[j] != ' '; j++);
    j += 2;
    for (i = 0; i < n; i++)
    {
        for (k = j; !(str[k] == ' ' || str[k] == '\0'); k++);
        memcpy(strcopy, str + j, k - j + 1);
        strcopy[k - j] = '\0';

        mpq_init(a[i]);
        mpq_set_str(a[i], strcopy, 10);
        mpq_canonicalize(a[i]);
        j = k + 1;
    }

    _fmpq_poly_from_list(rop, a, n);

    /* Clean-up */
    free(strcopy);
    for (i = 0; i < n; i++)
        mpq_clear(a[i]);
    free(a);

    return 1;
}

/**
 * \ingroup  StringConversions
 *
 * Returns the string representation of the rational polynomial \c op.
 */
char * fmpq_poly_to_string(const fmpq_poly_ptr op)
{
    ulong i, j;
    ulong len;     /* Upper bound on the length          */
    ulong denlen;  /* Size of the denominator in base 10 */
    mpz_t z;
    mpq_t q;

    char * str;

    if (fmpq_poly_is_zero(op))
    {
        str = (char *) malloc(2 * sizeof(char));
        if (!str)
        {
            printf("ERROR (fmpq_poly_to_string).  Memory allocation failed.\n");
            abort();
        }
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    mpz_init(z);
    if (fmpz_is_one(op->den))
    {
        denlen = 0;
    }
    else
    {
        fmpz_to_mpz(z, op->den);
        denlen = mpz_sizeinbase(z, 10);
    }
    len = fmpq_poly_places(fmpq_poly_length(op)) + 2;
    for (i = 0; i < fmpq_poly_length(op); i++)
    {
        fmpz_poly_get_coeff_mpz(z, op->num, i);
        len += mpz_sizeinbase(z, 10) + 1;
        if (mpz_sgn(z) != 0)
            len += 2 + denlen;
    }

    mpq_init(q);
    str = (char *) malloc(len * sizeof(char));
    if (!str)
    {
        printf("ERROR (fmpq_poly_to_string).  Memory allocation failed.\n");
        abort();
    }
    sprintf(str, "%lu", fmpq_poly_length(op));
    for (j = 0; str[j] != '\0'; j++);
    str[j++] = ' ';
    for (i = 0; i < fmpq_poly_length(op); i++)
    {
        str[j++] = ' ';
        fmpz_poly_get_coeff_mpz(mpq_numref(q), op->num, i);
        if (op->den == NULL)
            mpz_set_ui(mpq_denref(q), 1);
        else
            fmpz_to_mpz(mpq_denref(q), op->den);
        mpq_canonicalize(q);
        mpq_get_str(str + j, 10, q);
        for ( ; str[j] != '\0'; j++);
    }

    mpq_clear(q);
    mpz_clear(z);

    return str;
}

/**
 * \ingroup  StringConversions
 *
 * Returns the pretty string representation of \c op.
 *
 * Returns the pretty string representation of the rational polynomial \c op,
 * using the string \c var as the variable name.
 */
char * fmpq_poly_to_string_pretty(const fmpq_poly_ptr op, const char * var)
{
    long i;
    ulong j;
    ulong len;     /* Upper bound on the length          */
    ulong denlen;  /* Size of the denominator in base 10 */
    ulong varlen;  /* Length of the variable name        */
    mpz_t z;       /* op->den (if this is not 1)         */
    mpq_t q;
    char * str;

    if (fmpq_poly_is_zero(op))
    {
        str = (char *) malloc(2 * sizeof(char));
        if (!str)
        {
            printf("ERROR (fmpq_poly_to_string_pretty).  Memory allocation failed.\n");
            abort();
        }
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    if (fmpq_poly_length(op) == 1ul)
    {
        mpq_init(q);
        fmpz_to_mpz(mpq_numref(q), fmpz_poly_lead(op->num));
        if (fmpz_is_one(op->den))
            mpz_set_ui(mpq_denref(q), 1);
        else
        {
            fmpz_to_mpz(mpq_denref(q), op->den);
            mpq_canonicalize(q);
        }
        str = mpq_get_str(NULL, 10, q);
        mpq_clear(q);
        return str;
    }

    varlen = strlen(var);

    /* Copy the denominator into an mpz_t */
    mpz_init(z);
    if (fmpz_is_one(op->den))
    {
        denlen = 0;
    }
    else
    {
        fmpz_to_mpz(z, op->den);
        denlen = mpz_sizeinbase(z, 10);
    }

    /* Estimate the length */
    len = 0;
    for (i = 0; i < fmpq_poly_length(op); i++)
    {
        fmpz_poly_get_coeff_mpz(z, op->num, i);
        len += mpz_sizeinbase(z, 10);            /* Numerator                */
        if (mpz_sgn(z) != 0)
            len += 1 + denlen;                   /* Denominator and /        */
        len += 3;                                /* Operator and whitespace  */
        len += 1 + varlen + 1;                   /* *, x and ^               */
        len += fmpq_poly_places(i);              /* Exponent                 */
    }

    mpq_init(q);
    str = (char *) malloc(len * sizeof(char));
    if (!str)
    {
        printf("ERROR (fmpq_poly_to_string_pretty).  Memory allocation failed.\n");
        abort();
    }
    j = 0;

    /* Print the leading term */
    fmpz_to_mpz(mpq_numref(q), fmpz_poly_lead(op->num));
    fmpz_to_mpz(mpq_denref(q), op->den);
    mpq_canonicalize(q);

    if (mpq_cmp_ui(q, 1, 1) != 0)
    {
        if (mpq_cmp_si(q, -1, 1) == 0)
            str[j++] = '-';
        else
        {
            mpq_get_str(str, 10, q);
            for ( ; str[j] != '\0'; j++);
            str[j++] = '*';
        }
    }
    sprintf(str + j, "%s", var);
    j += varlen;
    if (fmpz_poly_degree(op->num) > 1)
    {
        str[j++] = '^';
        sprintf(str + j, "%li", fmpz_poly_degree(op->num));
        for ( ; str[j] != '\0'; j++);
    }

    for (i = fmpq_poly_degree(op) - 1; i >= 0; i--)
    {
        if (fmpz_is_zero(fmpz_poly_get_coeff_ptr(op->num, i)))
            continue;

        fmpz_poly_get_coeff_mpz(mpq_numref(q), op->num, i);
        fmpz_to_mpz(mpq_denref(q), op->den);
        mpq_canonicalize(q);

        str[j++] = ' ';
        if (mpq_sgn(q) < 0)
        {
            mpq_abs(q, q);
            str[j++] = '-';
        }
        else
            str[j++] = '+';
        str[j++] = ' ';

        mpq_get_str(str + j, 10, q);
        for ( ; str[j] != '\0'; j++);

        if (i > 0)
        {
            str[j++] = '*';
            sprintf(str + j, "%s", var);
            j += varlen;
            if (i > 1)
            {
                str[j++] = '^';
                sprintf(str + j, "%li", i);
                for ( ; str[j] != '\0'; j++);
            }
        }
    }

    mpq_clear(q);
    mpz_clear(z);
    return str;
}

#ifdef __cplusplus
    }
#endif


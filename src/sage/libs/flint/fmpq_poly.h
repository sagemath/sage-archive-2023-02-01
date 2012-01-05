///////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2010 Sebastian Pancratz                                     //
//                                                                           //
// Distributed under the terms of the GNU General Public License as          //
// published by the Free Software Foundation; either version 2 of the        //
// License, or (at your option) any later version.                           //
//                                                                           //
// http://www.gnu.org/licenses/                                              //
///////////////////////////////////////////////////////////////////////////////

#ifndef _FMPQ_POLY_H_
#define _FMPQ_POLY_H_

#ifdef __cplusplus
    extern "C" {
#endif

#include <stdio.h>
#include <gmp.h>

#include "fmpz.h"
#include "fmpz_poly.h"

/**
 * \file     fmpq_poly.h
 * \brief    Fast implementation of the rational polynomial ring
 * \author   Sebastian Pancratz
 * \date     Mar 2010 -- July 2010
 * \version  0.1.8
 */

/**
 * \ingroup  Definitions
 * \brief    Data type for a rational polynomial.
 * \details  We represent a rational polynomial as the quotient of an integer
 *           polynomial of type <tt>fmpz_poly_t</tt> and a denominator of type
 *           <tt>fmpz_t</tt>, enforcing coprimality and that the denominator
 *           is positive.  The zero polynomial is represented as \f$0/1\f$.
 */
typedef struct
{
    fmpz_poly_t num;  //!< Numerator
    fmpz_t den;       //!< Denominator
} fmpq_poly_struct;

/**
 * \ingroup  Definitions
 * \brief    Array of #fmpq_poly_struct of length one.
 * \details  This allows passing <em>by reference</em> without having to
 *           explicitly allocate memory for the structure, as one would have
 *           to when using pointers.
 */
typedef fmpq_poly_struct fmpq_poly_t[1];

/**
 * \ingroup  Definitions
 * \brief    Pointer to #fmpq_poly_struct.
 */
typedef fmpq_poly_struct * fmpq_poly_ptr;

void fmpq_poly_canonicalize(fmpq_poly_ptr f, fmpz_t temp);

///////////////////////////////////////////////////////////////////////////////
// fmpq_poly_*                                                               //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Accessing numerator and denominator

/**
 * \def      fmpq_poly_numref(op)
 * \brief    Returns a reference to the numerator of \c op.
 * \ingroup  NumDen
 * \details  The <tt>fmpz_poly</tt> functions can be used on the result of
 *           this macro.
 */
#define fmpq_poly_numref(op)  ((op)->num)

/**
 * \def      fmpq_poly_denref(op)
 * \brief    Returns a reference to the denominator of \c op.
 * \ingroup  NumDen
 * \details  The <tt>fmpz</tt> functions can be used on the result of
 *           this macro.
 */
#define fmpq_poly_denref(op)  ((op)->den)

///////////////////////////////////////////////////////////////////////////////
// Memory management

/**
 * \ingroup  MemoryManagement
 *
 * Initializes the element \c rop to zero.
 *
 * This function should be called once for the #fmpq_poly_ptr \c rop, before
 * using it with other <tt>fmpq_poly</tt> functions, or following a
 * preceeding call to #fmpq_poly_clear().
 */
static inline
void fmpq_poly_init(fmpq_poly_ptr rop)
{
    fmpz_poly_init(rop->num);
    rop->den = fmpz_init(1ul);
    fmpz_set_ui(rop->den, 1ul);
}

/**
 * \ingroup  MemoryManagement
 *
 * Clears the element \c rop.
 *
 * This function should only be called on an element which has been
 * initialised.
 */
static inline
void fmpq_poly_clear(fmpq_poly_ptr rop)
{
    fmpz_poly_clear(rop->num);
    fmpz_clear(rop->den);
}

///////////////////////////////////////////////////////////////////////////////
// Assignment and basic manipulation

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the same value as the element \c op.
 */
static inline
void fmpq_poly_set(fmpq_poly_ptr rop, const fmpq_poly_ptr op)
{
    if (rop != op)
    {
        fmpz_poly_set(rop->num, op->num);

        rop->den = fmpz_realloc(rop->den, fmpz_size(op->den));
        fmpz_set(rop->den, op->den);
    }
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the same value as the element \c op.
 */
static inline
void fmpq_poly_set_fmpz_poly(fmpq_poly_ptr rop, const fmpz_poly_t op)
{
    fmpz_poly_set(rop->num, op);
    fmpz_set_ui(rop->den, 1);
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the value given by the \c int \c op.
 */
static inline
void fmpq_poly_set_si(fmpq_poly_ptr rop, long op)
{
    fmpz_poly_zero(rop->num);
    fmpz_poly_set_coeff_si(rop->num, 0, op);
    fmpz_set_ui(rop->den, 1ul);
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the integer \c x.
 */
static inline
void fmpq_poly_set_fmpz(fmpq_poly_ptr rop, const fmpz_t x)
{
    fmpz_poly_zero(rop->num);
    fmpz_poly_set_coeff_fmpz(rop->num, 0, x);
    fmpz_set_ui(rop->den, 1ul);
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the integer \c x.
 */
static inline
void fmpq_poly_set_mpz(fmpq_poly_ptr rop, const mpz_t x)
{
    fmpz_poly_zero(rop->num);
    fmpz_poly_set_coeff_mpz(rop->num, 0, x);
    fmpz_set_ui(rop->den, 1ul);
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to the rational \c x, assumed to be given
 * in lowest terms.
 */
static inline
void fmpq_poly_set_mpq(fmpq_poly_ptr rop, const mpq_t x)
{
    fmpz_poly_zero(rop->num);
    fmpz_poly_set_coeff_mpz(rop->num, 0, mpq_numref(x));
    rop->den = fmpz_realloc(rop->den, mpz_size(mpq_denref(x)));
    mpz_to_fmpz(rop->den, mpq_denref(x));
}

/**
 * \ingroup  Assignment
 *
 * Swaps the elements \c op1 and \c op2.
 *
 * This is done efficiently by swapping pointers.
 */
static inline
void fmpq_poly_swap(fmpq_poly_ptr op1, fmpq_poly_ptr op2)
{
    if (op1 != op2)
    {
        fmpz_t t;
        fmpz_poly_swap(op1->num, op2->num);
        t        = op1->den;
        op1->den = op2->den;
        op2->den = t;
    }
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to zero.
 */
static inline
void fmpq_poly_zero(fmpq_poly_ptr rop)
{
    fmpz_poly_zero(rop->num);
    fmpz_set_ui(rop->den, 1ul);
}

/**
 * \ingroup  Assignment
 *
 * Sets the element \c rop to one.
 */
static inline
void fmpq_poly_one(fmpq_poly_ptr rop)
{
    fmpz_poly_zero(rop->num);
    fmpz_poly_set_coeff_si(rop->num, 0, 1);
    fmpz_set_ui(rop->den, 1ul);
}

void fmpq_poly_neg(fmpq_poly_ptr rop, const fmpq_poly_ptr op);
void fmpq_poly_inv(fmpq_poly_ptr rop, const fmpq_poly_ptr op);

///////////////////////////////////////////////////////////////////////////////
// Setting/ retrieving individual coefficients

void fmpq_poly_get_coeff_mpq(mpq_t rop, const fmpq_poly_ptr op, ulong i);

void fmpq_poly_set_coeff_fmpz(fmpq_poly_ptr rop, ulong i, const fmpz_t x);
void fmpq_poly_set_coeff_mpz(fmpq_poly_ptr rop, ulong i, const mpz_t x);
void fmpq_poly_set_coeff_mpq(fmpq_poly_ptr rop, ulong i, const mpq_t x);
void fmpq_poly_set_coeff_si(fmpq_poly_ptr rop, ulong i, long x);

///////////////////////////////////////////////////////////////////////////////
// Comparison

/**
 * \brief    Returns whether \c op is zero.
 * \ingroup  Comparison
 *
 * Returns whether the element \c op is zero.
 */
static inline
int fmpq_poly_is_zero(const fmpq_poly_ptr op)
{
    return op->num->length == 0ul;
}

/**
 * \brief    Returns whether \c op is equal to \f$1\f$.
 * \ingroup  Comparison
 *
 * Returns whether the element \c op is equal to the constant polynomial
 * \f$1\f$.
 */
static inline
int fmpq_poly_is_one(const fmpq_poly_ptr op)
{
    return (op->num->length == 1ul) && fmpz_is_one(op->num->coeffs)
                                    && fmpz_is_one(op->den);
}

int fmpq_poly_equal(const fmpq_poly_ptr op1, const fmpq_poly_ptr op2);
int fmpq_poly_cmp(const fmpq_poly_ptr left, const fmpq_poly_ptr right);

///////////////////////////////////////////////////////////////////////////////
// Polynomial parameters

/**
 * \brief    Returns the length of \c op.
 * \ingroup  PolyParameters
 * \details  Returns the length of the polynomial \c op, which is one greater
 *           than its degree.
 */
static inline
ulong fmpq_poly_length(const fmpq_poly_ptr op)
{
    return op->num->length;
}

/**
 * \brief    Returns the degree of \c op.
 * \ingroup  Parameters
 * \details  Returns the degree of the polynomial \c op.
 */
static inline
long fmpq_poly_degree(const fmpq_poly_ptr op)
{
    return (long) (op->num->length - 1ul);
}

///////////////////////////////////////////////////////////////////////////////
// Addition/ subtraction

void fmpq_poly_add(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2);
void fmpq_poly_sub(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2);

void fmpq_poly_addmul(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2);
void fmpq_poly_submul(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2);

///////////////////////////////////////////////////////////////////////////////
// Scalar multiplication and division

void fmpq_poly_scalar_mul_si(fmpq_poly_ptr rop, const fmpq_poly_ptr op, long x);
void fmpq_poly_scalar_mul_mpz(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpz_t x);
void fmpq_poly_scalar_mul_mpq(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpq_t x);

void fmpq_poly_scalar_div_si(fmpq_poly_ptr rop, const fmpq_poly_ptr op, long x);
void fmpq_poly_scalar_div_mpz(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpz_t x);
void fmpq_poly_scalar_div_mpq(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpq_t x);

///////////////////////////////////////////////////////////////////////////////
// Multiplication

void fmpq_poly_mul(fmpq_poly_ptr rop, const fmpq_poly_ptr op1, const fmpq_poly_ptr op2);

///////////////////////////////////////////////////////////////////////////////
// Division

void fmpq_poly_floordiv(fmpq_poly_ptr q, const fmpq_poly_ptr a, const fmpq_poly_ptr b);
void fmpq_poly_mod(fmpq_poly_ptr r, const fmpq_poly_ptr a, const fmpq_poly_ptr b);
void fmpq_poly_divrem(fmpq_poly_ptr q, fmpq_poly_ptr r, const fmpq_poly_ptr a, const fmpq_poly_ptr b);

///////////////////////////////////////////////////////////////////////////////
// Powering

void fmpq_poly_power(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong exp);

///////////////////////////////////////////////////////////////////////////////
// Greatest common divisor

void fmpq_poly_gcd(fmpq_poly_ptr rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b);
void fmpq_poly_xgcd(fmpq_poly_ptr rop, fmpq_poly_ptr s, fmpq_poly_ptr t, const fmpq_poly_ptr a, const fmpq_poly_ptr b);
void fmpq_poly_lcm(fmpq_poly_ptr rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b);

///////////////////////////////////////////////////////////////////////////////
// Derivative

void fmpq_poly_derivative(fmpq_poly_ptr rop, fmpq_poly_ptr op);

///////////////////////////////////////////////////////////////////////////////
// Evaluation

void fmpq_poly_evaluate_mpz(mpq_t rop, fmpq_poly_ptr f, const mpz_t a);
void fmpq_poly_evaluate_mpq(mpq_t rop, fmpq_poly_ptr f, const mpq_t a);

///////////////////////////////////////////////////////////////////////////////
// Gaussian content

void fmpq_poly_content(mpq_t rop, const fmpq_poly_ptr op);
void fmpq_poly_primitive_part(fmpq_poly_ptr rop, const fmpq_poly_ptr op);
int fmpq_poly_is_monic(const fmpq_poly_ptr op);
void fmpq_poly_monic(fmpq_poly_ptr rop, const fmpq_poly_ptr op);

///////////////////////////////////////////////////////////////////////////////
// Resultant

void fmpq_poly_resultant(mpq_t rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b);
void fmpq_poly_discriminant(mpq_t disc, fmpq_poly_t a);

///////////////////////////////////////////////////////////////////////////////
// Composition

void fmpq_poly_compose(fmpq_poly_ptr rop, const fmpq_poly_ptr a, const fmpq_poly_ptr b);
void fmpq_poly_rescale(fmpq_poly_ptr rop, const fmpq_poly_ptr op, const mpq_t x);

///////////////////////////////////////////////////////////////////////////////
// Square-free

int fmpq_poly_is_squarefree(const fmpq_poly_ptr op);

///////////////////////////////////////////////////////////////////////////////
// Subpolynomials

void fmpq_poly_getslice(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong i, ulong j);
void fmpq_poly_left_shift(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n);
void fmpq_poly_right_shift(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n);
void fmpq_poly_truncate(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n);
void fmpq_poly_reverse(fmpq_poly_ptr rop, const fmpq_poly_ptr op, ulong n);

///////////////////////////////////////////////////////////////////////////////
// String conversion

void _fmpq_poly_from_list(fmpq_poly_ptr rop, mpq_t * a, ulong n);
int fmpq_poly_from_string(fmpq_poly_ptr rop, const char * s);
char * fmpq_poly_to_string(const fmpq_poly_ptr op);
char * fmpq_poly_to_string_pretty(const fmpq_poly_ptr op, const char * x);

#ifdef __cplusplus
    }
#endif

#endif


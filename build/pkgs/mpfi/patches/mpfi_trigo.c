/* mpfi_trigo.c -- Trigonometric functions in mpfi

Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005,
                     Spaces project, Inria Lorraine
                     and Salsa project, INRIA Rocquencourt,
                     and Arenaire project, Inria Rhone-Alpes, France
                     and Lab. ANO, USTL (Univ. of Lille),  France


This file is part of the MPFI Library.

The MPFI Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFI Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFI Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */


#include "mpfi.h"
#include "mpfi-impl.h"

/* Conversion of a mpfr which has an integer value into a mpz */
static void mpfi_mpz_exact_set_fr (mpz_ptr z, mpfr_srcptr x)
{
  mp_exp_t expo;
  int sign;

  sign = MPFR_SIGN(x);
  expo = mpfr_get_z_exp(z, x);
  if (expo >= 0)
    mpz_mul_2exp(z, z, expo);
  else
    mpz_fdiv_q_2exp(z, z, -expo);

}

/* returns in quad the integer part of the division of x by Pi/2        */
/* the result is exact                                                  */
/* the returned value is the precision required to perform the division */
static mp_prec_t mpfi_quadrant (mpz_ptr quad, mpfr_srcptr x)
{
/* Assumption: x is neither a NaN nor an Infinite */
  int ok=0;
  mp_prec_t prec;
  mpfi_t two_over_pi, tmp;

  prec = mpfr_get_prec(x);

  if (MPFR_IS_ZERO(x)) {
    mpz_set_ui(quad, 0);
  }
  else {
    mpfi_init2(two_over_pi, prec);
    mpfi_init2(tmp, prec);

    do {
      mpfi_const_pi(two_over_pi);
      mpfi_ui_div(two_over_pi, 2, two_over_pi);

      mpfi_mul_fr(tmp, two_over_pi, x);
      mpfr_floor(&(tmp->left), &(tmp->left));
      mpfr_floor(&(tmp->right), &(tmp->right));

      ok =  mpfr_cmp(&(tmp->left), &(tmp->right));
      if (ok != 0) {
        prec += BITS_PER_MP_LIMB;
        mpfi_set_prec(two_over_pi, prec);
        mpfi_set_prec(tmp, prec);
      }
    } while (ok != 0);

    mpfi_mpz_exact_set_fr(quad, &(tmp->left));

    mpfi_clear(two_over_pi);
    mpfi_clear(tmp);
  }
 return prec;
}

/* compares z * Pi/2 - x and y where z is an integer (mpz)              */
/* the result is exact                                                  */
static int mpfi_cmp_sym_pi (mpz_srcptr z, mpfr_srcptr x, mpfr_srcptr y, mp_prec_t prec_init)
{
/* Assumption: x and y are neither NaN nor Infinite */
  mp_prec_t prec;
  mpfi_t pi_over_two, tmp;
  int not_ok;


  if (mpz_cmp_ui(z, 0) == 0) {
    /* We could do this faster by negating x in place, but mutating
       an input variable isn't thread-safe. */
    mpfr_t neg_x;
    int cmp_result;

    mpfr_init2(neg_x, mpfr_get_prec(x));
    mpfr_neg(neg_x, x, GMP_RNDN);
    cmp_result = mpfr_cmp(neg_x, y);
    mpfr_clear(neg_x);
    return cmp_result;
  }

  prec = prec_init;
  mpfi_init2(pi_over_two, prec);
  mpfi_init2(tmp, prec);

  do {
    mpfi_const_pi(pi_over_two);
    mpfi_div_2exp(pi_over_two, pi_over_two, 1);

    mpfi_mul_z(tmp, pi_over_two, z);
    mpfi_sub_fr(tmp, tmp, x);

    not_ok = mpfi_is_inside_fr(y, tmp);
    if (not_ok) {
      prec += BITS_PER_MP_LIMB;
      mpfi_set_prec(pi_over_two, prec);
      mpfi_set_prec(tmp, prec);
    }
  } while (not_ok);

  not_ok = mpfi_cmp_fr_default(tmp, y);
  mpfi_clear(pi_over_two);
  mpfi_clear(tmp);

  return not_ok;
}



int mpfi_sin(mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  mp_prec_t prec, prec_left, prec_right;
  mpfr_t tmp;
  mpz_t z, zmod4;
  mpz_t quad_left, quad_right;
  int ql_mod4, qr_mod4;

  if (MPFI_NAN_P(b)) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if (MPFI_INF_P(b)) {
    /* the two endpoints are the same infinite */
    if ( mpfr_cmp(&(b->left), &(b->right)) == 0) {
      mpfr_set_nan(&(a->left));
      mpfr_set_nan(&(a->right));
      MPFR_RET_NAN;
    }
    mpfr_set_si(&(a->left), -1, GMP_RNDD);
    mpfr_set_si(&(a->right), 1, GMP_RNDU);
    return 0;
  }

  mpz_init(quad_left);
  mpz_init(quad_right);
  mpz_init(z);
  /* quad_left gives the quadrant where the left endpoint of b is */
  /* quad_left = floor(2 b->left / pi)                            */
  /* idem for quad_right and b->right                             */
  prec_left = mpfi_quadrant(quad_left, &(b->left));
  prec_right = mpfi_quadrant(quad_right, &(b->right));


  /* if there is at least one period in b, then a = [-1, 1] */
  mpz_sub(z, quad_right, quad_left);
  if (mpz_cmp_ui(z, 4) >= 0) {
    mpfr_set_si(&(a->left), -1, GMP_RNDD);
    mpfr_set_si(&(a->right), 1, GMP_RNDU);
    inexact = 0;
  }
  else {
  /* there is less than one period in b */
  /* let us discuss according to the position (quadrant) of the endpoints of b    */

    /* computing precision = maximal precision required to determine the          */
    /* relative position of the endpoints of b and of integer multiples of Pi / 2 */
    prec = mpfi_get_prec(a);
    mpfr_init2(tmp, prec);

    if (prec_left > prec) prec = prec_left;
    if (prec_right > prec) prec = prec_right;

    mpz_add(z, quad_left, quad_right);
    /* z = quad_right + quad_left */

    mpz_init(zmod4);

    /* qr_mod4 gives the quadrant where the right endpoint of b is */
    /* qr_mod4 = floor(2 b->right / pi) mod 4 */
    mpz_mod_ui(zmod4, quad_right, 4);
    qr_mod4 = mpz_get_ui(zmod4);

    /* quad_left gives the quadrant where the left endpoint of b is */
    /* quad_left = floor(2 b->left / pi) mod 4 */
    mpz_mod_ui(zmod4, quad_left, 4);
    ql_mod4 = mpz_get_ui(zmod4);


    switch (qr_mod4) {
      case 0:
        switch (ql_mod4) {
          case 0:
          case 3:
            inexact_left = mpfr_sin(&(a->left), &(b->left), GMP_RNDD);
            inexact_right = mpfr_sin(&(a->right), &(b->right), GMP_RNDU);
            break;
          case 1:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->left), &(b->right), prec) >= 0)
              inexact_right = mpfr_sin(&(a->right), &(b->left), GMP_RNDU);
            else
              inexact_right = mpfr_sin(&(a->right), &(b->right), GMP_RNDU);
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            break;
          case 2:
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            inexact_right = mpfr_sin(&(a->right), &(b->right), GMP_RNDU);
            break;
        }
        break;
      case 1:
        switch (ql_mod4) {
          case 0:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->right), &(b->left), prec) >= 0)
              inexact_left = mpfr_sin(&(a->left), &(b->left), GMP_RNDD);
            else
              inexact_left = mpfr_sin(&(a->left), &(b->right), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
          case 1:
            inexact_left = mpfr_sin(tmp, &(b->right), GMP_RNDD);
            inexact_right = mpfr_sin(&(a->right), &(b->left), GMP_RNDU);
            mpfr_set(&(a->left), tmp, GMP_RNDD);
            break;
          case 2:
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
          case 3:
            inexact_left = mpfr_sin(&(a->left), &(b->left), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
        }
        break;
      case 2:
        switch (ql_mod4) {
          case 0:
            inexact_left = mpfr_sin(&(a->left), &(b->right), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
          case 1:
          case 2:
            inexact_left = mpfr_sin(tmp, &(b->right), GMP_RNDD);
            inexact_right = mpfr_sin(&(a->right), &(b->left), GMP_RNDU);
            mpfr_set(&(a->left), tmp, GMP_RNDD);
            break;
          case 3:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->left), &(b->right), prec) >= 0)
              inexact_left = mpfr_sin(&(a->left), &(b->left), GMP_RNDD);
            else
              inexact_left = mpfr_sin(&(a->left), &(b->right), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
        }
        break;
      case 3:
        switch (ql_mod4) {
          case 0:
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
          case 1:
            inexact_right = mpfr_sin(&(a->right), &(b->left), GMP_RNDU);
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            break;
          case 2:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->right), &(b->left), prec) >= 0)
              inexact_right = mpfr_sin(&(a->right), &(b->left), GMP_RNDU);
            else
              inexact_right = mpfr_sin(&(a->right), &(b->right), GMP_RNDU);
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            break;
          case 3:
            inexact_left = mpfr_sin(&(a->left), &(b->left), GMP_RNDD);
            inexact_right = mpfr_sin(&(a->right), &(b->right), GMP_RNDU);
            break;
        }
        break;
      }

    if (inexact_left) inexact = 1;
    if (inexact_right) inexact += 2;

    mpz_clear(zmod4);
    mpfr_clear (tmp);
  }

  mpz_clear(quad_left);
  mpz_clear(quad_right);
  mpz_clear(z);

  return inexact;
}



int mpfi_cos(mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  mp_prec_t prec, prec_left, prec_right;
  mpfr_t tmp;
  mpz_t z, zmod4;
  mpz_t quad_left, quad_right;
  int ql_mod4, qr_mod4;

  if (MPFI_NAN_P(b)) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if (MPFI_INF_P(b)) {
    /* the two endpoints are the same infinite */
    if ( mpfr_cmp(&(b->left), &(b->right)) == 0) {
      mpfr_set_nan(&(a->left));
      mpfr_set_nan(&(a->right));
      MPFR_RET_NAN;
    }
    mpfr_set_si(&(a->left), -1, GMP_RNDD);
    mpfr_set_si(&(a->right), 1, GMP_RNDU);
    return 0;
  }

  mpz_init(quad_left);
  mpz_init(quad_right);
  mpz_init(z);
  /* quad_left gives the quadrant where the left endpoint of b is */
  /* quad_left = floor(2 b->left / pi)                            */
  /* idem for quad_right and b->right                             */
  prec_left = mpfi_quadrant(quad_left, &(b->left));
  prec_right = mpfi_quadrant(quad_right, &(b->right));

  /* if there is at least one period in b, then a = [-1, 1] */
  mpz_sub(z, quad_right, quad_left);
  if (mpz_cmp_ui(z, 4) >= 0) {
    mpfr_set_si(&(a->left), -1, GMP_RNDD);
    mpfr_set_si(&(a->right), 1, GMP_RNDU);
    inexact = 0;
  }
  else {
  /* there is less than one period in b */
  /* let us discuss according to the position (quadrant) of the endpoints of b    */

    /* computing precision = maximal precision required to determine the          */
    /* relative position of the endpoints of b and of integer multiples of Pi / 2 */
    prec = mpfi_get_prec(a);
    mpfr_init2(tmp, prec);

    if (prec_left > prec) prec = prec_left;
    if (prec_right > prec) prec = prec_right;

    mpz_add(z, quad_left, quad_right);
    /* z = quad_right + quad_left */

    mpz_init(zmod4);

    /* qr_mod4 gives the quadrant where the right endpoint of b is */
    /* qr_mod4 = floor(2 b->right / pi) mod 4 */
    mpz_mod_ui(zmod4, quad_right, 4);
    qr_mod4 = mpz_get_ui(zmod4);

    /* quad_left gives the quadrant where the left endpoint of b is */
    /* quad_left = floor(2 b->left / pi) mod 4 */
    mpz_mod_ui(zmod4, quad_left, 4);
    ql_mod4 = mpz_get_ui(zmod4);


    switch (qr_mod4) {
      case 0:
        switch (ql_mod4) {
          case 0:
            inexact_left = mpfr_cos(tmp, &(b->right), GMP_RNDD);
            inexact_right = mpfr_cos(&(a->right), &(b->left), GMP_RNDU);
            mpfr_set(&(a->left), tmp, GMP_RNDD);
            break;
          case 1:
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
          case 2:
            inexact_left = mpfr_cos(&(a->left), &(b->left), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
          case 3:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->right), &(b->left), prec) >= 0)
              inexact_left = mpfr_cos(&(a->left), &(b->left), GMP_RNDD);
            else
              inexact_left = mpfr_cos(&(a->left), &(b->right), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
        }
        break;
      case 1:
        switch (ql_mod4) {
          case 0:
          case 1:
            inexact_left = mpfr_cos(tmp, &(b->right), GMP_RNDD);
            inexact_right = mpfr_cos(&(a->right), &(b->left), GMP_RNDU);
            mpfr_set(&(a->left), tmp, GMP_RNDD);
            break;
          case 2:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->right), &(b->left), prec) >= 0)
              inexact_left = mpfr_cos(&(a->left), &(b->left), GMP_RNDD);
            else
              inexact_left = mpfr_cos(&(a->left), &(b->right), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
          case 3:
            inexact_left = mpfr_cos(&(a->left), &(b->right), GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
        }
        break;
      case 2:
        switch (ql_mod4) {
          case 0:
            inexact_right = mpfr_cos(&(a->right), &(b->left), GMP_RNDU);
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            break;
          case 1:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->left), &(b->right), prec) >= 0)
              inexact_right = mpfr_cos(&(a->right), &(b->left), GMP_RNDU);
            else
              inexact_right = mpfr_cos(&(a->right), &(b->right), GMP_RNDU);
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            break;
          case 2:
            inexact_left = mpfr_cos(&(a->left), &(b->left), GMP_RNDD);
            inexact_right = mpfr_cos(&(a->right), &(b->right), GMP_RNDU);
            break;
          case 3:
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            inexact_right = mpfr_set_si(&(a->right), 1, GMP_RNDU);
            break;
        }
        break;
      case 3:
        switch (ql_mod4) {
          case 0:
            mpz_add_ui(z, z, 1);
            if (mpfi_cmp_sym_pi(z, &(b->right), &(b->left), prec) >= 0)
              inexact_right = mpfr_cos(&(a->right), &(b->left), GMP_RNDU);
            else
              inexact_right = mpfr_cos(&(a->right), &(b->right), GMP_RNDU);
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            break;
          case 1:
            inexact_right = mpfr_cos(&(a->right), &(b->right), GMP_RNDU);
            inexact_left = mpfr_set_si(&(a->left), -1, GMP_RNDD);
            break;
          case 2:
          case 3:
            inexact_left = mpfr_cos(&(a->left), &(b->left), GMP_RNDD);
            inexact_right = mpfr_cos(&(a->right), &(b->right), GMP_RNDU);
            break;
        }
        break;
      }

    if (inexact_left) inexact = 1;
    if (inexact_right) inexact += 2;

    mpz_clear(zmod4);
    mpfr_clear (tmp);
  }

  mpz_clear(quad_left);
  mpz_clear(quad_right);
  mpz_clear(z);

  return inexact;
}


int mpfi_tan(mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  mpz_t z_left, z_right, tmp;

  if (MPFI_NAN_P(b)) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if (MPFI_INF_P(b)) {
    /* the two endpoints are the same infinite */
    if ( mpfr_cmp(&(b->left), &(b->right)) == 0) {
      mpfr_set_nan(&(a->left));
      mpfr_set_nan(&(a->right));
      MPFR_RET_NAN;
    }
    mpfr_set_inf(&(a->left), -1);
    mpfr_set_inf(&(a->right), 1);
    return 0;
  }

  mpz_init(z_left);
  mpz_init(z_right);
  mpz_init(tmp);

  mpfi_quadrant(z_left, &(b->left));
  mpfi_quadrant(z_right, &(b->right));

  /* if there is at least one period in b or if b contains a Pi/2 + k*Pi, */
  /* then a = ]-oo, +oo[ */
  mpz_sub(tmp, z_right, z_left);
  if ( (mpz_cmp_ui(tmp, 2) >= 0) ||
       (mpz_even_p(z_left) && mpz_odd_p(z_right)) ) {
    mpfr_set_inf(&(a->left), -1);
    mpfr_set_inf(&(a->right), 1);
    inexact = 0;
  }

  else { /* within one period, tan is increasing */
    inexact_left = mpfr_tan(&(a->left), &(b->left), GMP_RNDD);
    inexact_right = mpfr_tan(&(a->right), &(b->right), GMP_RNDU);
    if (inexact_left) inexact += 1;
    if (inexact_right) inexact += 2;
  }

  mpz_clear(z_left);
  mpz_clear(z_right);
  mpz_clear(tmp);

  return inexact;
}



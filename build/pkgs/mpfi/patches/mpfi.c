/* mpfi.c -- Implementation mpfi functions.

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

/***************************************************************/
/*                   MPFI                                      */
/***************************************************************/


/************************************************/
/* Error                                        */
/************************************************/

static int mpfi_error=0;

#define MPFI_ERROR(s) \
{\
if(!mpfi_error) mpfi_error=1;fprintf(stderr,"\n%s\n",s);\
}

void mpfi_set_error(const int i)
{
  mpfi_error=i;
}

void mpfi_reset_error()
{
  mpfi_error=0;
}

int mpfi_is_error()
{
  return(mpfi_error==1);
}

/* Miscellaneous utilities */

#define MPFI_IS_POS(x) ((MPFR_SIGN((&(x->left)))>=0) && (MPFR_SIGN((&(x->right)))>0))
#define MPFI_IS_STRICTLY_POS(x) ((MPFR_SIGN((&(x->left)))>0) && (MPFR_SIGN((&(x->right)))>0))
#define MPFI_IS_NONNEG(x) ((MPFR_SIGN((&(x->left)))>=0) && (MPFR_SIGN((&(x->right)))>=0))
#define MPFI_IS_NEG(x) ((MPFR_SIGN((&(x->left)))<0) && (MPFR_SIGN((&(x->right)))<=0))
#define MPFI_IS_STRICTLY_NEG(x) ((MPFR_SIGN((&(x->left)))<0) && (MPFR_SIGN((&(x->right)))<0))
#define MPFI_IS_NONPOS(x) ((MPFR_SIGN((&(x->left)))<=0) && (MPFR_SIGN((&(x->right)))<=0))
#define MPFI_IS_NULL(x) ((MPFR_SIGN((&(x->left)))==0) && (MPFR_SIGN((&(x->right)))==0))
#define MPFI_HAS_ZERO(x) ((MPFR_SIGN((&(x->left)))<0) && (MPFR_SIGN((&(x->right)))>0))
#define MPFI_HAS_ZERO_NONSTRICT(x) ((MPFR_SIGN((&(x->left)))<=0) && (MPFR_SIGN((&(x->right)))>=0))

/* Default sign testing functions                  */

int mpfi_is_pos_default(mpfi_srcptr a)
{
  if ( mpfi_nan_p(a) )
    return 1;
  return(MPFI_IS_POS(a));
}

int mpfi_is_strictly_pos_default(mpfi_srcptr a)
{
  if ( mpfi_nan_p(a) )
    return 1;
  return(MPFI_IS_STRICTLY_POS(a));
}

int mpfi_is_nonneg_default(mpfi_srcptr a)
{
  if ( mpfi_nan_p(a) )
    return 1;
  return(MPFI_IS_NONNEG(a));
}

int mpfi_is_neg_default(mpfi_srcptr a)
{
  if ( mpfi_nan_p(a) )
    return 1;
  return(MPFI_IS_NEG(a));
}

int mpfi_is_strictly_neg_default(mpfi_srcptr a)
{
  if ( mpfi_nan_p(a) )
    return 1;
  return(MPFI_IS_STRICTLY_NEG(a));
}

int mpfi_is_nonpos_default(mpfi_srcptr a)
{
  if ( mpfi_nan_p(a) )
    return 1;
  return(MPFI_IS_NONPOS(a));
}

int mpfi_is_zero_default(mpfi_srcptr a)
{
  if ( mpfi_nan_p(a) )
    return 0;
  return((MPFR_SIGN(&(a->right))==0) &&
         (MPFR_SIGN(&(a->left))==0));
}


/* Customizable sign testing functions */

int (*mpfi_is_pos)  (mpfi_srcptr)=mpfi_is_pos_default;
int (*mpfi_is_strictly_pos)  (mpfi_srcptr)=mpfi_is_strictly_pos_default;
int (*mpfi_is_nonneg)  (mpfi_srcptr)=mpfi_is_nonneg_default;
int (*mpfi_is_neg)  (mpfi_srcptr)=mpfi_is_neg_default;
int (*mpfi_is_strictly_neg)  (mpfi_srcptr)=mpfi_is_strictly_neg_default;
int (*mpfi_is_nonpos)  (mpfi_srcptr)=mpfi_is_nonpos_default;
int (*mpfi_is_zero)  (mpfi_srcptr)=mpfi_is_zero_default;

int mpfi_has_zero (mpfi_srcptr a)
{
  return(   !MPFI_NAN_P(a)
         && (mpfr_cmp_ui(&(a->left),0)<=0)
         && (mpfr_cmp_ui(&(a->right),0)>=0) );
}



/************************************************/
/* Rounding                                     */
/************************************************/

/* The precision of the interval x is set to prec,     */
/* the previous value of x is kept.                    */
int mpfi_round_prec(mpfi_ptr x, mp_prec_t prec)
  {
  int inexact_left, inexact_right, inexact = 0;
  inexact_left = mpfr_round_prec (&(x->left), MPFI_RNDD, prec);
  inexact_right = mpfr_round_prec (&(x->right), MPFI_RNDU, prec);

  if ( MPFI_NAN_P(x) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
  }

/************************************************/
/* Initialization, destruction and assignment   */
/************************************************/

/* Interval initialization and destruction */
void mpfi_init (mpfi_t x)
{
  mpfr_init(&(x->left));
  mpfr_init(&(x->right));
}

void mpfi_init2 (mpfi_t x, mp_prec_t p)
{
  mpfr_init2(&(x->left),p);
  mpfr_init2(&(x->right),p);
}

/* ensures that the result is [a,b] with a<=b */
/* should be useless but bugs are not excluded... */
/* Returns 1 if endpoints have been exchanged     */

int mpfi_revert_if_needed (mpfi_ptr a)
{
  if ( MPFI_NAN_P(a) )
    return 0;

  if (mpfr_cmp(&(a->right),&(a->left))<0) {
    mpfr_swap(&(a->left), &(a->right));
    return 1;
  }
  else
    return 0;
}

/*
#define mpfi_revert_if_needed(a) 0
*/

void mpfi_clear(mpfi_ptr a)
{
  /* There is no test to check that the two endpoints are different
     and thus are not cleared twice. They should be different if
      only mpfi functions have been used...                        */
  mpfr_clear(&(a->right));
  mpfr_clear(&(a->left));
}

/* Returns the largest precision of the endpoints of x */
/* Reminder: the endpoints' precisions are supposed to be the same */
mp_prec_t mpfi_get_prec(mpfi_srcptr x)
{
  mp_prec_t prec_left, prec_right;
  prec_left=mpfr_get_prec(&(x->left));
  prec_right=mpfr_get_prec(&(x->right));
  return (prec_left>prec_right ? prec_left : prec_right);
}

/* The precision of the interval x is set to prec,     */
/* the previous value of x is lost.                    */
void mpfi_set_prec(mpfi_ptr x,mp_prec_t prec)
{
  mpfr_set_prec(&(x->right),prec);
  mpfr_set_prec(&(x->left),prec);
}

/************************************************/
/* Assignment functions                         */
/************************************************/

/* Interval assignments */

int mpfi_set(mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_set(&(a->left),&(b->left),MPFI_RNDD);
  inexact_right = mpfr_set(&(a->right),&(b->right),MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_set: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_set_si(mpfi_ptr a,const long b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_set_si(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_set_si(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_set_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_set_ui(mpfi_ptr a,const unsigned long b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_set_ui(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_set_ui(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_set_ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_set_d(mpfi_ptr a, const double b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_set_d(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_set_d(&(a->right),b,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_set_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_set_z(mpfi_ptr a, mpz_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_set_z(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_set_z(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_set_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_set_q(mpfi_ptr a, mpq_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_set_q(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_set_q(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_set_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_set_fr(mpfi_ptr a,mpfr_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_set(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_set(&(a->right),b,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_set_fr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

/* Combined initialization and assignment      */
int mpfi_init_set(mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_init_set(&(a->left),&(b->left),MPFI_RNDD);
  inexact_right = mpfr_init_set(&(a->right),&(b->right),MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_init_set: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_init_set_ui(mpfi_ptr a, const unsigned long b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_init_set_ui (&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_init_set_ui(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_init_set_ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_init_set_si(mpfi_ptr a, const long b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_init_set_si(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_init_set_si(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_init_set_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_init_set_d(mpfi_ptr a, const double b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_init_set_d(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_init_set_d(&(a->right),b,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_init_set_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_init_set_z(mpfi_ptr a, mpz_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_init_set_z(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_init_set_z(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_init_set_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_init_set_q(mpfi_ptr a, mpq_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_init_set_q(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_init_set_q(&(a->right),b,MPFI_RNDU);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_init_set_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int mpfi_init_set_fr(mpfi_ptr a, mpfr_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;
  inexact_left = mpfr_init_set(&(a->left),b,MPFI_RNDD);
  inexact_right = mpfr_init_set(&(a->right),b,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if ( mpfi_revert_if_needed(a) ) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_init_set_fr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

/* Swapping the two arguments */
void mpfi_swap (mpfi_ptr a, mpfi_ptr b)
{
  mpfr_swap(&(a->left), &(b->left));
  mpfr_swap(&(a->right), &(b->right));
}


/************************************************/
/* Various useful interval functions            */
/* with scalar results                          */
/************************************************/

/* Absolute diameter                            */
int mpfi_diam_abs(mpfr_ptr diam, mpfi_srcptr interv)
{
  return (mpfr_sub(diam,&(interv->right),&(interv->left), GMP_RNDU));
}

/* Relative diameter                                     */
/* Always computes an overestimation                     */
/* return value: 0 if the result is exact, > 0 otherwise */
int mpfi_diam_rel(mpfr_ptr diam, mpfi_srcptr interv)
{
  mpfr_t centre;
  int inexact_sub, inexact_mid, inexact_neg, inexact_sub2=0, inexact=0;

  inexact_sub = mpfr_sub(diam,&(interv->right),&(interv->left), GMP_RNDU);

  if (mpfi_bounded_p(interv)) {
    mpfr_init2(centre, mpfr_get_prec(diam));
    inexact_mid = mpfi_mid(centre, interv);

    if (mpfr_cmp_ui(centre, 0) <0)
      {
      inexact_neg = mpfr_neg(centre, centre, GMP_RNDD);
      if ( (!inexact_neg) || (inexact_mid<0) )
        mpfr_sub_one_ulp(centre, GMP_RNDD);
      }

    if (mpfr_cmp_ui(centre,0))
      inexact = mpfr_div(diam, diam, centre, GMP_RNDU);

    mpfr_clear(centre);
  } /* if interv is bounded, then a relative diameter can be computed */

  if ( mpfr_nan_p(diam) )
    MPFR_RET_NAN;

  if ( inexact || inexact_sub || inexact_sub2 )
    return 1;
  else
    return 0;
}

/* Diameter: relative if the interval does not contain 0 */
/* absolute otherwise                                    */
int mpfi_diam(mpfr_ptr diam, mpfi_srcptr interv)
{

  if (!mpfi_has_zero(interv)) {
    return (mpfi_diam_rel(diam, interv));
  }
  else {
    return (mpfi_diam_abs(diam, interv));
  }
}

/* Magnitude: the largest absolute value of any element */
int mpfi_mag(mpfr_ptr m, mpfi_srcptr x)
{
  int inexact;

  if (mpfi_is_nonneg_default(x))
    inexact = mpfr_set(m, &(x->right), GMP_RNDU);
  else if (mpfi_is_nonpos_default(x))
    inexact = mpfr_neg(m,&(x->left), GMP_RNDU);
  else { /* x contains 0 */
    inexact = mpfr_neg(m, &(x->left), GMP_RNDU);
    if (mpfr_nan_p(m)) {
      MPFR_RET_NAN;
    }
    else {
      if (mpfr_cmp(m, &(x->right)) < 0)
        inexact = mpfr_set(m, &(x->right), GMP_RNDU);
      if (mpfr_nan_p(m))
        MPFR_RET_NAN;
    }
  }
  return inexact;
}

/* Mignitude: the smallest absolute value of any element */
int mpfi_mig(mpfr_ptr m, mpfi_srcptr x)
{
  int inexact;

  if (mpfi_is_nonneg_default(x))
    inexact = mpfr_set(m, &(x->left), GMP_RNDD);
  else if (mpfi_is_nonpos_default(x))
    inexact = mpfr_neg(m,&(x->right), GMP_RNDD);
  else { /* x contains 0 */
    inexact = mpfr_set_ui(m, 0, GMP_RNDD);
  }
  return inexact;
}


/* Middle of y                                              */
/* With an IEEE-compliant arithmetic, this      */
/* formula provides more accurate results than  */
/* x->left + 0.5*(x->right - x->left) or        */
/* x->right - 0.5*(x->right - x->left)          */

int mpfi_mid (mpfr_ptr m, mpfi_srcptr y)
  {
  int inexact_add, inexact_div2;

  inexact_add = mpfr_add(m, &(y->left), &(y->right), GMP_RNDN);
  inexact_div2 = mpfr_div_2ui (m, m, 1, GMP_RNDN);

  /* Hope it copes correctly with an underflow in the division by 2... */
  if (inexact_div2)
    return inexact_div2;
  else
    return inexact_add;
  }

/* Picks randomly a point m in y */
void mpfi_alea (mpfr_ptr m, mpfi_srcptr y)
  {
  mp_prec_t prec;
  mpfr_t diam, fact;
  int dummy;
  prec=mpfr_get_prec(m);

  if ( MPFI_NAN_P(y) ) {
    mpfr_set_nan(m);
  }

  mpfr_init2(diam, prec);
  mpfr_init2(fact, prec);

  dummy = mpfi_diam_abs(diam, y);
  mpfr_random(fact);

  /* fact lies between 0 and 1, the picked point lies at a relative
     distance "fact" of the left endpoint:  m = inf + (sup - inf)*fact  */
  dummy = mpfr_mul(diam, diam, fact, GMP_RNDD);

  dummy = mpfr_add(m, &(y->left), diam, GMP_RNDD);
  if (!mpfr_nan_p(m)) {  /* Ensure that m belongs to y */
    if (mpfr_cmp(m, &(y->left)) < 0)
      dummy = mpfr_set(m, &(y->left), GMP_RNDU);

    if (mpfr_cmp(&(y->right), m) < 0)
      dummy = mpfr_set(m, &(y->right), GMP_RNDD);
    }

  mpfr_clear(diam);
  mpfr_clear(fact);
  }





/************************************************/
/* Conversions                                  */
/************************************************/

double mpfi_get_d (mpfi_srcptr a)
{
  double res;
  mpfr_t tmp;
  int dummy;

  mpfr_init2(tmp, 53);
  dummy = mpfi_mid(tmp, a);
  res = mpfr_get_d(tmp, GMP_RNDN);
  mpfr_clear(tmp);
  return res;
}

void mpfi_get_fr (mpfr_ptr m, mpfi_srcptr a)
{
  int dummy;
  dummy = mpfi_mid(m, a);
}

/************************************************/
/* Basic arithmetic operations                  */
/************************************************/

/* Arithmetic operations between two interval operands */

int  mpfi_add(mpfi_ptr a,mpfi_srcptr b,mpfi_srcptr c)
    /* relies on the fact that each operand is OK, */
    /* ie its left endpoint is <= its right one.   */
    /* This should always hold.                    */
{
  int inexact_left, inexact_right, inexact=0;

  if (MPFI_IS_ZERO(c)) {
    return (mpfi_set(a,b) );
  }
  else if (MPFI_IS_ZERO(b)) {
    return ( mpfi_set(a,c) );
  }
  else {
    inexact_left = mpfr_add(&(a->left),&(b->left),&(c->left),MPFI_RNDD);
    inexact_right = mpfr_add(&(a->right),&(b->right),&(c->right),MPFI_RNDU);
    if (MPFI_NAN_P(a))
      MPFR_RET_NAN;
    if (inexact_left)
      inexact += 1;
    if (inexact_right)
      inexact += 2;

    if (mpfi_revert_if_needed(a)) {
      /*
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_add\n");
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_add, 1st operand: ");
      mpfi_out_str(stdout, 10, 0, b);
      fprintf(stdout, "\n");
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_add, 2nd operand: ");
      mpfi_out_str(stdout, 10, 0, c);
      fprintf(stdout, "\n");
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_add, result: ");
      mpfi_out_str(stdout, 10, 0, a);
      fprintf(stdout, "\n");
      */
      inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
    }
    return inexact;
  }
}

int     mpfi_sub(mpfi_ptr a,mpfi_srcptr b,mpfi_srcptr c)
{
  mpfr_t tmp;
  int inexact_left, inexact_right, inexact=0;

  if (MPFI_IS_ZERO(c)) {
    return(mpfi_set(a,b));
   }
  else if (MPFI_IS_ZERO(b)) {
    return (mpfi_neg(a,c));
  }
  else {
    mpfr_init2(tmp, mpfi_get_prec(a));
    inexact_left = mpfr_sub(tmp, &(b->left), &(c->right), MPFI_RNDD);
    inexact_right = mpfr_sub(&(a->right), &(b->right), &(c->left), MPFI_RNDU);
    inexact_left |= mpfr_set(&(a->left), tmp, MPFI_RNDD);
    mpfr_clear(tmp);
    if (MPFI_NAN_P(a))
      MPFR_RET_NAN;
    if (inexact_left)
      inexact += 1;
    if (inexact_right)
      inexact += 2;

    if (mpfi_revert_if_needed(a)) {
      /*
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_sub\n");
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_sub, 1st operand: ");
      mpfi_out_str(stdout, 10, 0, b);
      fprintf(stdout, "\n");
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_sub, 2nd operand: ");
      mpfi_out_str(stdout, 10, 0, c);
      fprintf(stdout, "\n");
      fprintf(stdout, "Pb endpoints in reverse order in mpfi_sub, result: ");
      mpfi_out_str(stdout, 10, 0, a);
      fprintf(stdout, "\n");
      */
      inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
    }
    return inexact;
  }
}

int     mpfi_mul(mpfi_ptr a,mpfi_srcptr u,mpfi_srcptr c)
{
  mpfr_t t1;
  mpfr_t t2;
  int inexact_left, inexact_right, inexact_set_left=0, inexact_set_right=0, inexact=0;

  /* Handling the NaN cases */
  if ( MPFI_NAN_P(u) || MPFI_NAN_P(c) )
    {
      mpfr_set_nan(&(a->left));
      mpfr_set_nan(&(a->right));
      MPFR_RET_NAN;
    }

  /* Handling the case where one operand is 0, in order */
  /* to avoid problems with 0 * an infinite interval    */
  if (MPFI_IS_ZERO(u)) {
    return (mpfi_set(a,u));
  }
  if (MPFI_IS_ZERO(c)) {
    return (mpfi_set(a,c));
  }

  /* In the following, double rounding can occur: in order to cope with a result equal
     to one argument, a multiplication is performed and stored in a temporary variable
     and then assigned to the corresponding endpoint.                                  */
  if (MPFR_SIGN(&(u->left))>=0) {
    if (MPFR_SIGN(&(c->left))>=0) {
      /* u nonnegative and c nonnegative */
      inexact_left = mpfr_mul(&(a->left),&(u->left),&(c->left),MPFI_RNDD);
      inexact_right = mpfr_mul(&(a->right),&(u->right),&(c->right),MPFI_RNDU);
    }
    else {
      if (MPFR_SIGN(&(c->right))<=0) {
	/* u nonnegative and c non-positive */
        mpfr_init2(t1, mpfi_get_prec(a));
        inexact_left = mpfr_mul(t1, &(u->right), &(c->left), MPFI_RNDD);
        inexact_right = mpfr_mul(&(a->right), &(u->left), &(c->right), MPFI_RNDU);
        inexact_set_left = mpfr_set(&(a->left), t1, MPFI_RNDD);
        mpfr_clear(t1);
      }
      else {
	/* u nonnegative and c overlapping 0 */
        mpfr_init2(t1, mpfi_get_prec(a));
	inexact_left = mpfr_mul(t1,&(u->right),&(c->left),MPFI_RNDD);
	inexact_right = mpfr_mul(&(a->right),&(u->right),&(c->right),MPFI_RNDU);
        inexact_set_left = mpfr_set(&(a->left), t1, MPFI_RNDD);
        mpfr_clear(t1);
      }
    }
  }
  else {
    if (MPFR_SIGN(&(u->right))<=0) {
      /* u non-positive */
      if (MPFR_SIGN(&(c->left))>=0) {
        /* u non-positive and c nonnegative */
	mpfr_init2(t1,mpfi_get_prec(a));
        inexact_left = mpfr_mul(t1,&(u->left),&(c->right),MPFI_RNDD);
        inexact_right = mpfr_mul(&(a->right),&(u->right),&(c->left),MPFI_RNDU);
        inexact_set_left = mpfr_set(&(a->left), t1, MPFI_RNDD);
	mpfr_clear(t1);
      }
      else {
        if (MPFR_SIGN(&(c->right))<=0) {
	  /* u non-positive and c non-positive */
          mpfr_init2(t1, mpfi_get_prec(a));
          inexact_left = mpfr_mul(t1, &(u->right), &(c->right), MPFI_RNDD);
          inexact_right = mpfr_mul(&(a->right), &(u->left), &(c->left), MPFI_RNDU);
          inexact_set_left = mpfr_set(&(a->left), t1, MPFI_RNDD);
          mpfr_clear(t1);
        }
        else {
	  /* u non-positive and c overlapping 0 */
          mpfr_init2(t1, mpfi_get_prec(a));
	  inexact_left = mpfr_mul(t1,&(u->left),&(c->right),MPFI_RNDD);
	  inexact_right = mpfr_mul(&(a->right),&(u->left),&(c->left),MPFI_RNDU);
          inexact_set_left = mpfr_set(&(a->left), t1, MPFI_RNDD);
          mpfr_clear(t1);
        }
      }
    }
    else {
      /* u contains 0 */
      if (MPFR_SIGN(&(c->left))>=0) {
	/* u overlapping 0 and c nonnegative  */
        mpfr_init2(t1, mpfi_get_prec(a));
	inexact_left = mpfr_mul(t1,&(u->left),&(c->right),MPFI_RNDD);
	inexact_right = mpfr_mul(&(a->right),&(u->right),&(c->right),MPFI_RNDU);
        inexact_set_left = mpfr_set(&(a->left), t1, MPFI_RNDD);
        mpfr_clear(t1);
      }
      else {
	if (MPFR_SIGN(&(c->right))<=0) {
	  /* u overlapping 0 and c non-positive */
          mpfr_init2(t1, mpfi_get_prec(a));
	  inexact_left = mpfr_mul(t1,&(u->right),&(c->left),MPFI_RNDD);
	  inexact_right = mpfr_mul(&(a->right),&(u->left),&(c->left),MPFI_RNDU);
          inexact_set_left = mpfr_set(&(a->left), t1, MPFI_RNDD);
          mpfr_clear(t1);
	}
	else {
	  /* u overlapping 0 and c overlapping 0
	     Beware the case where the result is one of the operands! */
	  mpfr_init2(t1,mpfi_get_prec(a));
	  mpfr_init2(t2,mpfi_get_prec(a));
	  inexact_right = mpfr_mul(t1,&(u->left),&(c->right),MPFI_RNDD);
	  inexact_left = mpfr_mul(t2,&(u->right),&(c->left),MPFI_RNDD);
	  if (mpfr_cmp(t1,t2)<0) {
	    mpfr_swap(t2,t1);
            inexact_left=inexact_right;
	  }
	  inexact_right = mpfr_mul(t1,&(u->left),&(c->left),MPFI_RNDU);
	  inexact_set_right = mpfr_mul(&(a->right),&(u->right),&(c->right),MPFI_RNDU);
	  if (mpfr_cmp(t1,&(a->right))>0) {
	    inexact_set_right = mpfr_set(&(a->right),t1,MPFI_RNDU);
	  }
          else
            inexact_right = 0;
	  inexact_set_left = mpfr_set(&(a->left),t2,MPFI_RNDD);
	  mpfr_clear(t1);
	  mpfr_clear(t2);
	}
      }
    }
  }

  if (inexact_left || inexact_set_left)
    inexact += 1;
  if (inexact_right || inexact_set_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int     mpfi_div(mpfi_ptr a,mpfi_srcptr b,mpfi_srcptr c)
{
  mpfr_t tmp;
  int inexact_left=0, inexact_right=0, inexact=0;

  if ( MPFI_NAN_P(b) || MPFI_NAN_P(c) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if (MPFI_HAS_ZERO_NONSTRICT(c)) {
    if ( MPFI_HAS_ZERO(c) || MPFI_HAS_ZERO(b) ) { /* a = ]-oo, +oo [ */
      mpfr_set_inf(&(a->left), -1);
      mpfr_set_inf(&(a->right), 1);
    }
    else if (MPFI_IS_NONNEG(c)) {  /* c >= 0 and its left endpoint is 0 */
      if (MPFI_IS_NONNEG(b)) {                    /* a = [ bl/cr, +oo [ */
        inexact_left = mpfr_div(&(a->left), &(b->left), &(c->right), MPFI_RNDD);
        mpfr_set_inf(&(a->right), 1);
      }
      else { /* b <= 0 */                         /* a = ] -oo, br/cr ] */
        inexact_right = mpfr_div(&(a->right), &(b->right), &(c->right), MPFI_RNDU);
        mpfr_set_inf(&(a->left), -1);
      }
    }
    else { /* c <= 0 and its right endpoint is 0 */
      if (MPFI_IS_NONNEG(b)) {                    /* a = ] -oo, bl/cl ] */
        inexact_right = mpfr_div(&(a->right), &(b->left), &(c->left), MPFI_RNDU);
        mpfr_set_inf(&(a->left), -1);
      }
      else { /* b <= 0 */                         /* a = [ br/cl, +oo [ */
        inexact_left = mpfr_div(&(a->left), &(b->right), &(c->left), MPFI_RNDD);
        mpfr_set_inf(&(a->right), 1);
      }
    }
  }
  else if (MPFI_IS_POS(c)) {
    mpfr_init2(tmp, mpfi_get_prec(a));
    if (MPFI_IS_NONNEG(b)) {                       /* a = [ bl/cr, br/cl ] */
      inexact_left  = mpfr_div(tmp, &(b->left), &(c->right), MPFI_RNDD);
      inexact_right = mpfr_div(&(a->right), &(b->right), &(c->left), MPFI_RNDU);
      inexact_left |= mpfr_set(&(a->left), tmp, MPFI_RNDD);
    }
    else if (MPFI_IS_NONPOS(b)) {                 /* a = [ bl/cl, br/cr ] */
      inexact_left  = mpfr_div(&(a->left), &(b->left), &(c->left), MPFI_RNDD);
      inexact_right = mpfr_div(&(a->right), &(b->right), &(c->right), MPFI_RNDU);
    }
    else { /* b contains 0 in its interior */     /* a = [ bl/cl, br/cl ] */
      inexact_right = mpfr_div(tmp, &(b->right), &(c->left), MPFI_RNDU);
      inexact_left  = mpfr_div(&(a->left), &(b->left), &(c->left), MPFI_RNDD);
      inexact_right|= mpfr_set(&(a->right), tmp, MPFI_RNDU);
    }
    mpfr_clear(tmp);
  }
  else { /* c < 0 */
    mpfr_init2(tmp, mpfi_get_prec(a));
    if (MPFI_IS_NONNEG(b)) {                       /* a = [ br/cr, bl/cl ] */
      inexact_left  = mpfr_div(tmp, &(b->right), &(c->right), MPFI_RNDD);
      inexact_right = mpfr_div(&(a->right), &(b->left), &(c->left), MPFI_RNDU);
      inexact_left |= mpfr_set(&(a->left), tmp, MPFI_RNDD);
    }
    else if (MPFI_IS_NONPOS(b)) {                 /* a = [ br/cl, bl/cr ] */
      inexact_left  = mpfr_div(tmp, &(b->right), &(c->left), MPFI_RNDD);
      inexact_right = mpfr_div(&(a->right), &(b->left), &(c->right), MPFI_RNDU);
      inexact_left |= mpfr_set(&(a->left), tmp, MPFI_RNDD);
    }
    else { /* b contains 0 in its interior */     /* a = [ br/cr, bl/cr ] */
      inexact_left  = mpfr_div(tmp, &(b->right), &(c->right), MPFI_RNDD);
      inexact_right = mpfr_div(&(a->right), &(b->left), &(c->right), MPFI_RNDU);
      inexact_left |= mpfr_set(&(a->left), tmp, MPFI_RNDD);
    }
    mpfr_clear(tmp);
  }

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}


/* Arithmetic operations between an interval operand */
/* and a double prec. floating-point                 */

int   mpfi_add_d(mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  int inexact_set, inexact_add, inexact=0;

  mpfi_init2(tmp, 53);
  inexact_set = mpfi_set_d(tmp, c);
  inexact_add = mpfi_add(a, b, tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_add)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_add))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_add) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_add))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_add_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int   mpfi_sub_d(mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp, 53);
  inexact_set = mpfi_set_d(tmp,c);
  inexact_sub = mpfi_sub(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_sub_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int   mpfi_d_sub(mpfi_ptr a, const double b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp, 53);
  inexact_set = mpfi_set_d(tmp,b);
  inexact_sub = mpfi_sub(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_d_sub: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }
  return inexact;
}

int   mpfi_mul_d(mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  int inexact_set, inexact_mul, inexact=0;

  mpfi_init2(tmp, 53);
  inexact_set = mpfi_set_d(tmp,c);
  inexact_mul = mpfi_mul(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_mul))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_mul))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_mul;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_div_d(mpfi_ptr a, mpfi_srcptr b, const double c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp, 53);
  inexact_set = mpfi_set_d(tmp,c);
  inexact_div = mpfi_div(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_d_div(mpfi_ptr a, const double b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp, 53);
  inexact_set = mpfi_set_d(tmp,b);
  inexact_div = mpfi_div(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of b into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_d_div: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* Arithmetic operations between an interval operand */
/* and a long integer                                */

int   mpfi_add_ui(mpfi_ptr a, mpfi_srcptr b, const unsigned long c)
{
  int inexact_left, inexact_right, inexact=0;

  if (c==0) {
    return (mpfi_set(a,b) );
  }
  else if (MPFI_IS_ZERO(b)) {
    return ( mpfi_set_ui(a,c) );
  }
  else {
    inexact_left = mpfr_add_ui(&(a->left),&(b->left),c,MPFI_RNDD);
    inexact_right = mpfr_add_ui(&(a->right),&(b->right),c,MPFI_RNDU);
    if (MPFI_NAN_P(a))
      MPFR_RET_NAN;
    if (inexact_left)
      inexact += 1;
    if (inexact_right)
      inexact += 2;

    if (mpfi_revert_if_needed(a)) {
      inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
    }
    return inexact;
  }
}

int   mpfi_sub_ui(mpfi_ptr a, mpfi_srcptr b, const unsigned long c)
{

  int inexact_left, inexact_right, inexact=0;

  if (c==0) {
    return(mpfi_set(a,b));
   }
  else if (MPFI_IS_ZERO(b)) {
    mpfi_set_ui(a,c);
    return (mpfi_neg(a,a));
  }
  else {
    inexact_left  = mpfr_sub_ui(&(a->left), &(b->left),c, MPFI_RNDD);
    inexact_right = mpfr_sub_ui(&(a->right), &(b->right), c, MPFI_RNDU);
    if (MPFI_NAN_P(a))
      MPFR_RET_NAN;
    if (inexact_left)
      inexact += 1;
    if (inexact_right)
      inexact += 2;

    if (mpfi_revert_if_needed(a)) {
      inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
    }
    return inexact;
  }
}

int   mpfi_ui_sub(mpfi_ptr a, const unsigned long b, mpfi_srcptr c)
{

  int inexact_left, inexact_right, inexact=0;

  if (MPFI_IS_ZERO(c)) {
    return(mpfi_set_ui(a,b));
   }
  else if (b==0) {
    return (mpfi_sub(a,a,c));
  }
  else {
    inexact_left  = mpfr_ui_sub(&(a->left), b,&(c->right), MPFI_RNDD);
    inexact_right = mpfr_ui_sub(&(a->right),b,&(c->left), MPFI_RNDU);
    if (MPFI_NAN_P(a))
      MPFR_RET_NAN;
    if (inexact_left)
      inexact += 1;
    if (inexact_right)
      inexact += 2;

    if (mpfi_revert_if_needed(a)) {
      inexact=MPFI_REVERT_INEXACT_FLAGS(inexact);
    }
    return inexact;
  }
}

int   mpfi_mul_ui(mpfi_ptr a, mpfi_srcptr b, const unsigned long c)
{
  mpfi_t tmp;
  int inexact_set, inexact_mul, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_ui(tmp,c);
  inexact_mul = mpfi_mul(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_mul))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_mul))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_mul;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_div_ui(mpfi_ptr a, mpfi_srcptr b, const unsigned long c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_ui(tmp,c);
  inexact_div = mpfi_div(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_ui_div(mpfi_ptr a, const unsigned long b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_ui(tmp,b);
  inexact_div = mpfi_div(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of b into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_ui_div: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* Arithmetic operations between an interval operand */
/* and a long integer                                */

int   mpfi_add_si(mpfi_ptr a, mpfi_srcptr b, const long c)
{
  mpfi_t tmp;
  int inexact_set, inexact_add, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_si(tmp,c);
  inexact_add = mpfi_add(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_add)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_add))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_add) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_add))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_add_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_sub_si(mpfi_ptr a, mpfi_srcptr b, const long c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_si(tmp,c);
  inexact_sub = mpfi_sub(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_sub_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_si_sub(mpfi_ptr a, const long b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_si(tmp,b);
  inexact_sub = mpfi_sub(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_si_sub: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_mul_si(mpfi_ptr a, mpfi_srcptr b, const long c)
{
  mpfi_t tmp;
  int inexact_set, inexact_mul, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_si(tmp,c);
  inexact_mul = mpfi_mul(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_mul))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_mul))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_mul;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_div_si(mpfi_ptr a, mpfi_srcptr b, const long c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_si(tmp,c);
  inexact_div = mpfi_div(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_si_div(mpfi_ptr a, const long b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_si(tmp,b);
  inexact_div = mpfi_div(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of b into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_si_div: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* Arithmetic operations between an interval operand */
/* and a multiple precision integer                  */

int   mpfi_add_z(mpfi_ptr a, mpfi_srcptr b, mpz_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_add, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_z(tmp,c);
  inexact_add = mpfi_add(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_add)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_add))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_add) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_add))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_add_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_sub_z(mpfi_ptr a, mpfi_srcptr b, mpz_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_z(tmp,c);
  inexact_sub = mpfi_sub(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_sub_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_z_sub(mpfi_ptr a, mpz_srcptr b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_z(tmp,b);
  inexact_sub = mpfi_sub(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_z_sub: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_mul_z(mpfi_ptr a, mpfi_srcptr b, mpz_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_mul, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_z(tmp,c);
  inexact_mul = mpfi_mul(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_mul))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_mul))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_mul;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_div_z(mpfi_ptr a, mpfi_srcptr b, mpz_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_z(tmp,c);
  inexact_div = mpfi_div(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_z_div(mpfi_ptr a, mpz_srcptr b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_z(tmp,b);
  inexact_div = mpfi_div(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of b into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_z_div: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* Arithmetic operations between an interval operand */
/* and a multiple precision rational                 */

int   mpfi_add_q(mpfi_ptr a, mpfi_srcptr b, mpq_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_add, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_q(tmp,c);
  inexact_add = mpfi_add(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_add)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_add))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_add) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_add))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_add_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_sub_q(mpfi_ptr a, mpfi_srcptr b, mpq_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_q(tmp,c);
  inexact_sub = mpfi_sub(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_sub_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_q_sub(mpfi_ptr a, mpq_srcptr b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_q(tmp,b);
  inexact_sub = mpfi_sub(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_q_sub: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_mul_q(mpfi_ptr a, mpfi_srcptr b, mpq_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_mul, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_q(tmp,c);
  inexact_mul = mpfi_mul(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_mul))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_mul))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_mul;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_div_q(mpfi_ptr a, mpfi_srcptr b, mpq_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_q(tmp,c);
  inexact_div = mpfi_div(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_q_div(mpfi_ptr a, mpq_srcptr b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_q(tmp,b);
  inexact_div = mpfi_div(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of b into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_q_div: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* Arithmetic operations between an interval operand */
/* and a multiple precision floating-point nb        */
/* To do without converting to a mpfi first,         */
/* to avoid double roundings and save time           */

int   mpfi_add_fr(mpfi_ptr a, mpfi_srcptr b, mpfr_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_add, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_fr(tmp,c);
  inexact_add = mpfi_add(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_add)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_add))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_add) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_add))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_add: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_sub_fr(mpfi_ptr a, mpfi_srcptr b, mpfr_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_fr(tmp,c);
  inexact_sub = mpfi_sub(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_sub_fr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_fr_sub(mpfi_ptr a, mpfr_srcptr b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_sub, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_fr(tmp,b);
  inexact_sub = mpfi_sub(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if ( mpfr_inf_p(&(a->left)) ) {
    if  (MPFI_LEFT_IS_INEXACT(inexact_sub)) /* overflow */
    inexact += 1;
  }
  else if (MPFI_LEFT_IS_INEXACT(inexact_set) || MPFI_LEFT_IS_INEXACT(inexact_sub))
    inexact += 1;
  if ( mpfr_inf_p(&(a->right)) ) {
    if (MPFI_RIGHT_IS_INEXACT(inexact_sub) )  /* overflow */
    inexact += 2;
  }
  else if (MPFI_RIGHT_IS_INEXACT(inexact_set) || MPFI_RIGHT_IS_INEXACT(inexact_sub))
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_fr_sub: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_mul_fr(mpfi_ptr a, mpfi_srcptr b, mpfr_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_mul, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_fr(tmp,c);
  inexact_mul = mpfi_mul(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_mul))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_mul))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_mul;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_fr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_div_fr(mpfi_ptr a, mpfi_srcptr b, mpfr_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_fr(tmp,c);
  inexact_div = mpfi_div(a,b,tmp);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of c into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_fr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int   mpfi_fr_div(mpfi_ptr a, mpfr_srcptr b, mpfi_srcptr c)
{
  mpfi_t tmp;
  int inexact_set, inexact_div, inexact=0;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_fr(tmp,b);
  inexact_div = mpfi_div(a,tmp,c);
  MPFI_CLEAR(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (mpfr_inf_p(&(a->left))  && MPFI_LEFT_IS_INEXACT(inexact_div))   /* overflow */
    inexact += 1;
  if (mpfr_inf_p(&(a->right)) && MPFI_RIGHT_IS_INEXACT(inexact_div))  /* overflow */
    inexact += 2;
  if (mpfi_bounded_p(a)) {
    if (inexact_set) /* if the conversion of b into a mpfi is inexact,
                        then so are both endpoints of the result.      */
      inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
    else
      inexact = inexact_div;
  }

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_fr_div: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* Arithmetic operations taking a single interval operand */

int mpfi_neg(mpfi_ptr a, mpfi_srcptr b)
{
  mpfr_t tmp;
  int inexact_left, inexact_right, inexact=0;

  mpfr_init2(tmp,mpfi_get_prec(a));
  inexact_right = mpfr_set(tmp,&(b->left),MPFI_RNDD);
  inexact_left = mpfr_neg(&(a->left),&(b->right),MPFI_RNDD);
  inexact_right |= mpfr_neg(&(a->right),tmp,MPFI_RNDU);
  mpfr_clear(tmp);

  if (MPFI_NAN_P(a))
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stdout, "Pb endpoints in reverse order in mpfi_neg, operand: ");
    mpfi_out_str(stdout, 10, 0, b);
    fprintf(stdout, "Pb endpoints in reverse order in mpfi_neg, result: ");
    mpfi_out_str(stdout, 10, 0, a);
    fprintf(stdout, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int     mpfi_sqr(mpfi_ptr a,mpfi_srcptr u)
{
  mpfr_t t1;
  int inexact_left, inexact_right, inexact=0;

  if ( MPFI_NAN_P(u) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }
  if (MPFR_SIGN(&(u->left))>=0) {
      /* u nonnegative */
      inexact_left = mpfr_mul(&(a->left),&(u->left),&(u->left),MPFI_RNDD);
      inexact_right = mpfr_mul(&(a->right),&(u->right),&(u->right),MPFI_RNDU);
  }
  else {
    if (MPFR_SIGN(&(u->right))<=0) {
      /* u non-positive -> beware the case where a = u */
      mpfr_init2(t1,mpfi_get_prec(a));
      inexact_right = mpfr_mul(t1, &(u->left), &(u->left), MPFI_RNDU);
      inexact_left = mpfr_mul(&(a->left), &(u->right), &(u->right), MPFI_RNDD);
      inexact_right |= mpfr_set(&(a->right), t1, MPFI_RNDU);
      mpfr_clear(t1);
    }
    else {
      /* inf = 0, sup = max of the squares of the endpoints of u */
      mpfr_init2(t1,mpfi_get_prec(u));
      inexact_right = mpfr_abs(t1, &(u->left), MPFI_RNDU);
      if (mpfr_cmp(t1, &(u->right)) <=0) {
        inexact_right = mpfr_mul(&(a->right), &(u->right), &(u->right), MPFI_RNDU);
      }
      else {
        inexact_right = mpfr_mul(&(a->right), &(u->left), &(u->left), MPFI_RNDU);
      }
      mpfr_clear(t1);
      inexact_left = mpfr_set_si(&(a->left), (long)0, GMP_RNDZ);
    }
  }

  /* The NaN case has already been handled */
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_sqr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_inv(mpfi_ptr a, mpfi_srcptr b)
{
  mpfr_t tmp;
  int inexact_left, inexact_right, inexact=0;

  if ( MPFI_NAN_P(b) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }
  else if ( (MPFR_SIGN(&(b->left))<0) && (MPFR_SIGN(&(b->right))>0) ) {
    /* b has endpoints of opposite signs */
    /* b can even be [0-, 0+] */
    /* The result is [-oo, +oo] and both endpoints are exact. */
    /*MPFI_ERROR("mpfi_inv: inversion by an interval that contains 0");*/
    mpfr_set_inf( &(a->left), -1);
    mpfr_set_inf( &(a->right), 1);
  }
  else { /* signed zeroes and Infinity are properly handled by MPFR */
      mpfr_init2(tmp,mpfi_get_prec(a));
      inexact_right = mpfr_ui_div(tmp,1,&(b->left),MPFI_RNDU);
      inexact_left = mpfr_ui_div(&(a->left),1,&(b->right),MPFI_RNDD);
      inexact_right |= mpfr_set(&(a->right),tmp,MPFI_RNDU);
      if (inexact_left)
        inexact += 1;
      if (inexact_right)
        inexact += 2;
      if (mpfi_revert_if_needed(a) ) {
        /* Not an error, but due to lazy programming
        fprintf(stderr, "Pb endpoints in reverse order in mpfi_inv: ");
        mpfi_out_str(stderr, 10, 0, a);
        fprintf(stderr, "\n");
        */
        inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
      }
      mpfr_clear(tmp);
  }

  return inexact;
}

int mpfi_sqrt(mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  /* if b is (partially) negative, the left bound will be a NaN */
  /* it is handled by MPFR */
  inexact_left = mpfr_sqrt(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_sqrt(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_abs(mpfi_ptr a, mpfi_srcptr b)
{
  /* the result a contains the absolute values of every element of b */
  mpfr_t tmp;
  int inexact_left, inexact_right, inexact=0;

  if (MPFI_NAN_P(b)) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if (mpfi_is_nonneg_default(b))
    inexact = mpfi_set(a, b);
  else if (mpfi_is_nonpos_default(b))
    inexact = mpfi_neg(a,b);
  else { /* b contains 0 */
    mpfr_init2(tmp, mpfi_get_prec(a));
    inexact_right = mpfr_neg(tmp, &(b->left), MPFI_RNDU);
    if (mpfr_cmp(tmp, &(b->right)) < 0) {
      inexact_right = mpfr_set(&(a->right), &(b->right), MPFI_RNDU);
    }
    else {
      inexact_right |= mpfr_set(&(a->right), tmp, MPFI_RNDU);
    }
    inexact_left = mpfr_set_si(&(a->left), 0, MPFI_RNDD);
    mpfr_clear(tmp);

    if (inexact_left)
      inexact += 1;
    if (inexact_right)
      inexact += 2;
  }

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;
  else
    return inexact;
}

/* Miscellaneous arithmetic operations                    */

int mpfi_mul_2exp(mpfi_ptr a, mpfi_srcptr b,unsigned long c)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_mul_2exp(&(a->left),&(b->left),c,MPFI_RNDD);
  inexact_right = mpfr_mul_2exp(&(a->right),&(b->right),c,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_2exp: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


int mpfi_mul_2ui(mpfi_ptr a, mpfi_srcptr b,unsigned long c)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_mul_2ui(&(a->left),&(b->left),c,MPFI_RNDD);
  inexact_right = mpfr_mul_2ui(&(a->right),&(b->right),c,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_2ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


int mpfi_mul_2si(mpfi_ptr a, mpfi_srcptr b, long c)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_mul_2si(&(a->left),&(b->left),c,MPFI_RNDD);
  inexact_right = mpfr_mul_2si(&(a->right),&(b->right),c,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_mul_2si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_div_2exp(mpfi_ptr a, mpfi_srcptr b,unsigned long c)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_div_2exp(&(a->left),&(b->left),c,MPFI_RNDD);
  inexact_right = mpfr_div_2exp(&(a->right),&(b->right),c,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_2exp: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_div_2ui(mpfi_ptr a, mpfi_srcptr b,unsigned long c)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_div_2ui(&(a->left),&(b->left),c,MPFI_RNDD);
  inexact_right = mpfr_div_2ui(&(a->right),&(b->right),c,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_2ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_div_2si (mpfi_ptr a, mpfi_srcptr b, long c)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_div_2si (&(a->left),&(b->left),c,MPFI_RNDD);
  inexact_right = mpfr_div_2si (&(a->right),&(b->right),c,MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_div_2si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}





/************************************************/
/* Special functions                            */
/************************************************/

/* Computes the log of an interval              */
/* Special cases are handled by mpfr_log        */
int mpfi_log (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_log(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_log(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

/* Computes the exp of an interval              */
int mpfi_exp (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_exp(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_exp(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_exp2 (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_exp2(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_exp2(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_acos (mpfi_ptr a, mpfi_srcptr b)
{
  mpfr_t tmp;
  int inexact_left, inexact_right, inexact=0;

  mpfr_init2(tmp, mpfi_get_prec(a));
  inexact_left = mpfr_acos(tmp, &(b->right), MPFI_RNDD);
  inexact_right = mpfr_acos(&(a->right), &(b->left), MPFI_RNDU);
  inexact_left |= mpfr_set(&(a->left), tmp, MPFI_RNDD);
  mpfr_clear(tmp);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_asin (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_asin(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_asin(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_atan (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_atan(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_atan(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_cosh (mpfi_ptr a, mpfi_srcptr b)
{
  mpfr_t tmp;
  int inexact_left, inexact_right, inexact=0;

  if ( MPFI_NAN_P(b) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if ( MPFI_IS_NONNEG(b) ) {
    inexact_left = mpfr_cosh(&(a->left), &(b->left), MPFI_RNDD);
    inexact_right = mpfr_cosh(&(a->right), &(b->right), MPFI_RNDU);
  }
  else if ( MPFI_HAS_ZERO(b) ) {
    mpfr_init2(tmp, mpfi_get_prec(a));
    inexact_right = mpfr_neg(tmp, &(b->left), MPFI_RNDU);
    if (mpfr_cmp(tmp, &(b->right)) > 0)
      inexact_right |= mpfr_cosh( &(a->right), tmp, MPFI_RNDU);
    else
      inexact_right = mpfr_cosh( &(a->right), &(b->right), MPFI_RNDU);
    inexact_left = mpfr_set_ui(&(a->left), 1, MPFI_RNDD);
    mpfr_clear(tmp);
  }
  else { /* b <= 0 */
    mpfr_init2(tmp, mpfi_get_prec(a));
    inexact_left = mpfr_cosh(tmp, &(b->right), MPFI_RNDD);
    inexact_right = mpfr_cosh(&(a->right), &(b->left), MPFI_RNDU);
    inexact_left |= mpfr_set(&(a->left), tmp, MPFI_RNDD);
    mpfr_clear(tmp);
  }

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_sinh (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_sinh(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_sinh(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_tanh (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_tanh(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_tanh(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_acosh (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_acosh(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_acosh(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_asinh (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_asinh(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_asinh(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_atanh (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_atanh(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_atanh(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

/* Computes the log of 1 plus an interval       */
/* Special cases are handled by mpfr_log        */
int mpfi_log1p (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_log1p(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_log1p(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

/* Computes the exp of an interval, minus 1     */
int mpfi_expm1 (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_expm1(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_expm1(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

/* Computes the logarithm in base 2 of an interval     */
int mpfi_log2 (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_log2(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_log2(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

/* Computes the logarithm in base 10 of an interval     */
int mpfi_log10 (mpfi_ptr a, mpfi_srcptr b)
{
  int inexact_left, inexact_right, inexact=0;

  inexact_left = mpfr_log10(&(a->left), &(b->left), MPFI_RNDD);
  inexact_right = mpfr_log10(&(a->right), &(b->right), MPFI_RNDU);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  return inexact;
}

int mpfi_const_log2(mpfi_ptr a)
{
  mpfr_const_log2(&(a->left), MPFI_RNDD);
  mpfr_const_log2(&(a->right), MPFI_RNDU);

  return MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
}

int mpfi_const_pi(mpfi_ptr a)
{
  mpfr_const_pi(&(a->left), MPFI_RNDD);
  mpfr_const_pi(&(a->right), MPFI_RNDU);

  return MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
}

int mpfi_const_euler(mpfi_ptr a)
{
  mpfr_const_euler(&(a->left), MPFI_RNDD);
  mpfr_const_euler(&(a->right), MPFI_RNDU);

  return MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
}


/************************************************/
/* Comparison functions                         */
/* Warning: there is no canonical order =>      */
/* interval comparison is not clearly defined   */
/************************************************/

/* Default comparison functions */
/* They return 1 if one (at least) of their operands is invalid (contains NaN) */

int mpfi_cmp_default (mpfi_srcptr a,mpfi_srcptr b)
{
  if ( MPFI_NAN_P(a) || MPFI_NAN_P(b) )
    return 1;
  return((mpfr_cmp(&(a->right),&(b->left))<0) ? -1 :(mpfr_cmp(&(b->right),&(a->left))<0));
}

int (*mpfi_cmp)(mpfi_srcptr,mpfi_srcptr)=mpfi_cmp_default;


int mpfi_cmp_d_default (mpfi_srcptr a,const double b)
{
  int dummy, res=0;
  mpfi_t tmp;
  mpfi_init2(tmp,mpfi_get_prec(a));
  dummy = mpfi_set_d(tmp,b);
  res=mpfi_cmp(a,tmp);
  MPFI_CLEAR(tmp);
  return(res);
}

int mpfi_cmp_ui_default (mpfi_srcptr a,const unsigned long b)
{
  int dummy, res=0;
  mpfi_t tmp;
  mpfi_init2(tmp,mpfi_get_prec(a));
  dummy = mpfi_set_ui(tmp,b);
  res=mpfi_cmp(a,tmp);
  MPFI_CLEAR(tmp);
  return(res);
}

int mpfi_cmp_si_default (mpfi_srcptr a,const long b)
{
  int dummy, res=0;
  mpfi_t tmp;
  mpfi_init2(tmp,mpfi_get_prec(a));
  dummy = mpfi_set_si(tmp,b);
  res=mpfi_cmp(a,tmp);
  MPFI_CLEAR(tmp);
  return(res);
}

int mpfi_cmp_z_default (mpfi_srcptr a,mpz_srcptr b)
{
  int dummy, res=0;
  mpfi_t tmp;
  mpfi_init2(tmp,mpfi_get_prec(a));
  dummy = mpfi_set_z(tmp,b);
  res=mpfi_cmp(a,tmp);
  MPFI_CLEAR(tmp);
  return(res);
}

int mpfi_cmp_q_default (mpfi_srcptr a,mpq_srcptr b)
{
  int dummy, res=0;
  mpfi_t tmp;
  mpfi_init2(tmp,mpfi_get_prec(a));
  dummy = mpfi_set_q(tmp,b);
  res=mpfi_cmp(a,tmp);
  MPFI_CLEAR(tmp);
  return(res);
}

int mpfi_cmp_fr_default (mpfi_srcptr a,mpfr_srcptr b)
{
  int dummy, res=0;
  mpfi_t tmp;
  mpfi_init2(tmp,mpfi_get_prec(a));
  dummy = mpfi_set_fr(tmp,b);
  res=mpfi_cmp(a,tmp);
  MPFI_CLEAR(tmp);
  return(res);
}

/* Customizable comparison functions */
/* Since the mpfi_cmp_* are based on mpfi_cmp, only mpfi_cmp needs to be modified */

int    (*mpfi_cmp_d)   (mpfi_srcptr,const double)=mpfi_cmp_d_default;
int    (*mpfi_cmp_ui)  (mpfi_srcptr,const unsigned long)=mpfi_cmp_ui_default;
int    (*mpfi_cmp_si)  (mpfi_srcptr,const long)=mpfi_cmp_si_default;
int    (*mpfi_cmp_z)   (mpfi_srcptr,mpz_srcptr)=mpfi_cmp_z_default;
int    (*mpfi_cmp_q)   (mpfi_srcptr,mpq_srcptr)=mpfi_cmp_q_default;
int    (*mpfi_cmp_fr)   (mpfi_srcptr,mpfr_srcptr)=mpfi_cmp_fr_default;

/* Predicates */

int mpfi_nan_p (mpfi_srcptr a)
{
  return ( mpfr_nan_p(&(a->left)) || mpfr_nan_p(&(a->right)) );
}

int mpfi_inf_p (mpfi_srcptr a)
{
  return ( mpfr_inf_p(&(a->left)) || mpfr_inf_p(&(a->right)) );
}

int mpfi_bounded_p (mpfi_srcptr a)
{
  return ( mpfr_number_p(&(a->left)) && mpfr_number_p(&(a->right)) );
}


/************************************************/
/* Interval manipulation                        */
/************************************************/

/* Operations related to the internal representation */
/* by endpoints                                      */

int   mpfi_get_left (mpfr_ptr b,mpfi_srcptr a)
{
  return mpfr_set(b,&(a->left),MPFI_RNDD);
}

int   mpfi_get_right (mpfr_ptr b,mpfi_srcptr a)
{
  return mpfr_set(b,&(a->right),MPFI_RNDU);
}


/* The result is the first interval extended so */
/* that it contains the second (scalar) one.    */
/* Also called "convex hull".                   */

int mpfi_put(mpfi_ptr a,mpfi_srcptr b)
{
  int inexact_left=0, inexact_right=0, inexact=0;

  if ( MPFI_NAN_P(b) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }
  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if ( mpfr_cmp(&(a->left),&(b->left))>0 ) {
    inexact_left = mpfr_set(&(a->left),&(b->left),MPFI_RNDD);
  }
  if (mpfr_cmp(&(a->right),&(b->right))<0) {
    inexact_right = mpfr_set(&(a->right),&(b->right),MPFI_RNDU);
  }
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_put: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_put_d(mpfi_ptr a,const double b)
{
  mpfi_t tmp;
  int inexact_set=0, inexact_left=0, inexact_right=0, inexact=0;

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_d(tmp,b);
  if ( MPFI_NAN_P(tmp) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
  }
  else {
    if (mpfr_cmp(&(a->left), &(tmp->left)) > 0 ) {
      inexact_left = mpfr_set(&(a->left), &(tmp->left), MPFI_RNDD);
      if ( MPFI_LEFT_IS_INEXACT(inexact_set) )
        inexact_left = 1;
      }
    if (mpfr_cmp(&(a->right), &(tmp->right)) < 0 ) {
      inexact_right = mpfr_set(&(a->right), &(tmp->right), MPFI_RNDD);
      if ( MPFI_RIGHT_IS_INEXACT(inexact_set) )
        inexact_right = 1;
      }
  }
  MPFI_CLEAR(tmp);

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_put_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_put_si(mpfi_ptr a,const long b)
{
  int inexact_left=0, inexact_right=0, inexact=0;

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if ( mpfr_cmp_si(&(a->left),b)>0 ) {
    inexact_left = mpfr_set_si(&(a->left),b,MPFI_RNDD);
  }
  if (mpfr_cmp_si(&(a->right),b)<0) {
    inexact_right = mpfr_set_si(&(a->right),b,MPFI_RNDU);
  }
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_put_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_put_ui(mpfi_ptr a,const unsigned long b)
{
  int inexact_left=0, inexact_right=0, inexact=0;

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if ( mpfr_cmp_ui(&(a->left),b)>0 ) {
    inexact_left = mpfr_set_ui(&(a->left),b,MPFI_RNDD);
  }
  if (mpfr_cmp_ui(&(a->right),b)<0) {
    inexact_right = mpfr_set_ui(&(a->right),b,MPFI_RNDU);
  }
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_put_ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_put_z(mpfi_ptr a,mpz_srcptr b)
{
  mpfi_t tmp;
  int inexact_set=0, inexact_left=0, inexact_right=0, inexact=0;

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_z(tmp,b);
  /* No need to test if tmp is a NaN, it is not possible with mpfr_set_z */
  if (mpfr_cmp(&(a->left), &(tmp->left)) > 0 ) {
    inexact_left = mpfr_set(&(a->left), &(tmp->left), MPFI_RNDD);
    if ( MPFI_LEFT_IS_INEXACT(inexact_set) )
      inexact_left = 1;
    }
  if (mpfr_cmp(&(a->right), &(tmp->right)) < 0 ) {
    inexact_right = mpfr_set(&(a->right), &(tmp->right), MPFI_RNDD);
    if ( MPFI_RIGHT_IS_INEXACT(inexact_set) )
      inexact_right = 1;
    }
  MPFI_CLEAR(tmp);

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_put_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_put_q(mpfi_ptr a,mpq_srcptr b)
{
  mpfi_t tmp;
  int inexact_set=0, inexact_left=0, inexact_right=0, inexact=0;

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  mpfi_init2(tmp,mpfi_get_prec(a));
  inexact_set = mpfi_set_q(tmp,b);
  /* No need to test if tmp is a NaN, it is not possible with mpfr_set_z */
  if (mpfr_cmp(&(a->left), &(tmp->left)) > 0 ) {
    inexact_left = mpfr_set(&(a->left), &(tmp->left), MPFI_RNDD);
    if ( MPFI_LEFT_IS_INEXACT(inexact_set) )
      inexact_left = 1;
    }
  if (mpfr_cmp(&(a->right), &(tmp->right)) < 0 ) {
    inexact_right = mpfr_set(&(a->right), &(tmp->right), MPFI_RNDD);
    if ( MPFI_RIGHT_IS_INEXACT(inexact_set) )
      inexact_right = 1;
    }
  MPFI_CLEAR(tmp);

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_put_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_put_fr(mpfi_ptr a,mpfr_srcptr b)
{
  int inexact_left=0, inexact_right=0, inexact=0;

  if ( mpfr_nan_p(b) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }
  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if ( mpfr_cmp(&(a->left),b)>0 ) {
    inexact_left = mpfr_set(&(a->left),b,MPFI_RNDD);
  }
  if (mpfr_cmp(&(a->right),b)<0) {
    inexact_right = mpfr_set(&(a->right),b,MPFI_RNDU);
  }
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_put_fr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* The result is an interval containing the     */
/* two scalar arguments.                        */

int mpfi_interv_d(mpfi_ptr a,const double b,const double c)
{
  int inexact_left, inexact_right, inexact=0;

  if (b<=c) {
    inexact_left = mpfr_set_d(&(a->left), b, MPFI_RNDD);
    inexact_right = mpfr_set_d(&(a->right), c, MPFI_RNDU);
  }
  else {
    inexact_left = mpfr_set_d(&(a->left), c, MPFI_RNDD);
    inexact_right = mpfr_set_d(&(a->right), b, MPFI_RNDU);
  }

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_interv_d: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_interv_si(mpfi_ptr a,const long b,const long c)
{
  int inexact_left, inexact_right, inexact=0;
  mpfi_set_si(a,b);
  mpfi_put_si(a,c);

  if (b<=c) {
    inexact_left = mpfr_set_si(&(a->left), b, MPFI_RNDD);
    inexact_right = mpfr_set_si(&(a->right), c, MPFI_RNDU);
  }
  else {
    inexact_left = mpfr_set_si(&(a->left), c, MPFI_RNDD);
    inexact_right = mpfr_set_si(&(a->right), b, MPFI_RNDU);
  }
  /* a cannot be a NaN, mpfr_set_si never returns a NaN */
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_interv_si: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_interv_ui(mpfi_ptr a,const unsigned long b,const unsigned long c)
{
  int inexact_left, inexact_right, inexact=0;

  if (b<=c) {
    inexact_left = mpfr_set_ui(&(a->left), b, MPFI_RNDD);
    inexact_right = mpfr_set_ui(&(a->right), c, MPFI_RNDU);
  }
  else {
    inexact_left = mpfr_set_ui(&(a->left), c, MPFI_RNDD);
    inexact_right = mpfr_set_ui(&(a->right), b, MPFI_RNDU);
  }
  /* a cannot be a NaN, mpfr_set_ui never returns a NaN */
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_interv_ui: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_interv_z(mpfi_ptr a,mpz_srcptr b,mpz_srcptr c)
{
  int inexact_left, inexact_right, inexact=0;

  if ( mpz_cmp(b,c) <= 0 ) {
    inexact_left = mpfr_set_z(&(a->left), b, MPFI_RNDD);
    inexact_right = mpfr_set_z(&(a->right), c, MPFI_RNDU);
  }
  else {
    inexact_left = mpfr_set_z(&(a->left), c, MPFI_RNDD);
    inexact_right = mpfr_set_z(&(a->right), b, MPFI_RNDU);
  }
  /* a cannot be a NaN, mpfr_set_z never returns a NaN */
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_interv_z: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_interv_q(mpfi_ptr a,mpq_srcptr b,mpq_srcptr c)
{
  int inexact_left, inexact_right, inexact=0;

  if ( mpq_cmp(b,c) <= 0) {
    inexact_left = mpfr_set_q(&(a->left), b, MPFI_RNDD);
    inexact_right = mpfr_set_q(&(a->right), c, MPFI_RNDU);
  }
  else {
    inexact_left = mpfr_set_q(&(a->left), c, MPFI_RNDD);
    inexact_right = mpfr_set_q(&(a->right), b, MPFI_RNDU);
  }
  /* a cannot be a NaN, mpfr_set_q never returns a NaN */
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_interv_q: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}

int mpfi_interv_fr(mpfi_ptr a,mpfr_srcptr b,mpfr_srcptr c)
{
  int inexact_left, inexact_right, inexact=0;

  if ( mpfr_nan_p(b) || mpfr_nan_p(c) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if ( mpfr_cmp(b,c) <= 0 ) {
    inexact_left = mpfr_set(&(a->left), b, MPFI_RNDD);
    inexact_right = mpfr_set(&(a->right), c, MPFI_RNDU);
  }
  else {
    inexact_left = mpfr_set(&(a->left), c, MPFI_RNDD);
    inexact_right = mpfr_set(&(a->right), b, MPFI_RNDU);
  }
  /* a cannot be a NaN, it has been tested before */
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    /*
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_interv_fr: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    */
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}



/* Inclusion tests: tests whether the 1st op is in the 2nd one      */

int   mpfi_is_strictly_inside    (mpfi_srcptr a,mpfi_srcptr b)
{
  /* Returns 1 if one of the operands is a NaN */
  if ( MPFI_NAN_P(a) || MPFI_NAN_P(b) )
    return 1;
  return(   ((mpfr_cmp(&(b->left),&(a->left))<0) &&
	     (mpfr_cmp(&(a->right),&(b->right))<0)) );
}

int   mpfi_is_inside    (mpfi_srcptr a,mpfi_srcptr b)
{
  /* Returns 1 if one of the operands is a NaN */
  if ( MPFI_NAN_P(a) || MPFI_NAN_P(b) )
    return 1;
  return((mpfr_cmp(&(b->left),&(a->left))<=0) &&
	 (mpfr_cmp(&(a->right),&(b->right))<=0));
}

int   mpfi_is_inside_d     (const double a, mpfi_srcptr b)
{
  int dummy, res;
  mpfi_t tmp;
  mpfi_init2(tmp,mpfi_get_prec(b));
  dummy = mpfi_set_d(tmp,a);
  res=mpfi_is_inside(tmp,b);
  MPFI_CLEAR(tmp);
  return(res);
}

int   mpfi_is_inside_ui     (const unsigned long a, mpfi_srcptr b)
{
  /* Returns 1 if one of the operands is a NaN */
  if ( MPFI_NAN_P(b) )
    return 1;
  return ( (mpfr_cmp_ui(&(b->left),a) <= 0) &&
           (mpfr_cmp_ui(&(b->right),a) >= 0) );
}

int   mpfi_is_inside_si     (const long a, mpfi_srcptr b)
{
  /* Returns 1 if one of the operands is a NaN */
  if ( MPFI_NAN_P(b) )
    return 1;
  return ( (mpfr_cmp_si(&(b->left),a) <= 0) &&
           (mpfr_cmp_si(&(b->right),a) >= 0) );
}

int   mpfi_is_inside_z     (mpz_srcptr a, mpfi_srcptr b)
{
  int dummy, res;
  mpfi_t tmp;
  /* Returns 1 if one of the operands is a NaN */
  if ( MPFI_NAN_P(b) )
    return 1;
  mpfi_init2(tmp,mpfi_get_prec(b));
  dummy = mpfi_set_z(tmp,a);
  res=mpfi_is_inside(tmp,b);
  MPFI_CLEAR(tmp);
  return(res);
}

int   mpfi_is_inside_q     (mpq_srcptr a, mpfi_srcptr b)
{
  int dummy, res;
  mpfi_t tmp;
  /* Returns 1 if one of the operands is a NaN */
  if ( MPFI_NAN_P(b) )
    return 1;
  mpfi_init2(tmp,mpfi_get_prec(b));
  dummy = mpfi_set_q(tmp,a);
  res=mpfi_is_inside(tmp,b);
  MPFI_CLEAR(tmp);
  return(res);
}

int   mpfi_is_inside_fr     (mpfr_srcptr a, mpfi_srcptr b)
{
  /* Returns 1 if one of the operands is a NaN */
  if ( mpfr_nan_p(a) || MPFI_NAN_P(b) )
    return 1;
  return ( (mpfr_cmp(a,&(b->left)) >= 0) &&
           (mpfr_cmp(a,&(b->right)) <= 0) );
}

/* Set operations                                 */
/* Be careful with empty results: nothing is done */
/* to handle empty intervals in math functions    */
int   mpfi_is_empty         (mpfi_srcptr a)
{
  return ( MPFI_NAN_P(a) || (mpfr_cmp(&(a->left), &(a->right)) > 0) );
}

int  mpfi_intersect        (mpfi_ptr a, mpfi_srcptr b, mpfi_srcptr c)
{
/* If intersection = empty set, a will have its left bound > right one */
  int inexact_left, inexact_right, inexact=0;

  if ( MPFI_NAN_P(b) || MPFI_NAN_P(c) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if (mpfr_cmp(&(b->left), &(c->left)) <= 0)
    inexact_left = mpfr_set(&(a->left), &(c->left), MPFI_RNDD);
  else
    inexact_left = mpfr_set(&(a->left), &(b->left), MPFI_RNDD);

  if (mpfr_cmp(&(c->right), &(b->right)) <= 0)
    inexact_right = mpfr_set(&(a->right), &(c->right), MPFI_RNDU);
  else
    inexact_right = mpfr_set(&(a->right), &(b->right), MPFI_RNDU);

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;
  if (mpfr_cmp(&(a->left), &(a->right)) > 0) /* a is the empty set */
    inexact = MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;

  return inexact;
}

int  mpfi_union            (mpfi_ptr a, mpfi_srcptr b, mpfi_srcptr c)
{
/* Union of two intervals, without gap: convex hull of the union */
  int inexact_left, inexact_right, inexact=0;

  if ( MPFI_NAN_P(b) || MPFI_NAN_P(c) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  if (mpfr_cmp(&(b->left), &(c->left)) <= 0)
    inexact_left = mpfr_set(&(a->left), &(b->left), MPFI_RNDD);
  else
    inexact_left = mpfr_set(&(a->left), &(c->left), MPFI_RNDD);

  if (mpfr_cmp(&(c->right), &(b->right)) <= 0)
    inexact_right = mpfr_set(&(a->right), &(b->right), MPFI_RNDU);
  else
    inexact_right = mpfr_set(&(a->right), &(c->right), MPFI_RNDU);

  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_union: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\n");
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}



/************************************************/
/* Miscellaneous                                */
/************************************************/

int mpfi_increase(mpfi_ptr a,mpfr_srcptr e)
{
  int inexact_left, inexact_right, inexact=0;

  if ( MPFI_NAN_P(a) )
    MPFR_RET_NAN;
  else if ( mpfr_nan_p(e) ) {
    mpfr_set_nan(&(a->left));
    mpfr_set_nan(&(a->right));
    MPFR_RET_NAN;
  }

  inexact_right = mpfr_add(&(a->right),&(a->right),e,MPFI_RNDU);
  inexact_left = mpfr_sub(&(a->left),&(a->left),e,MPFI_RNDD);
  if (inexact_left)
    inexact += 1;
  if (inexact_right)
    inexact += 2;

  if (mpfi_revert_if_needed(a)) {
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_increase: ");
    mpfi_out_str(stderr, 10, 0, a);
    fprintf(stderr, "\nincreased by (a possibly negative value): ");
    mpfr_out_str(stderr, 10, 0, e, GMP_RNDN);
    fprintf(stderr, "\n");
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* keeps the same center and multiply the radius by (1+fact) */
int mpfi_blow(mpfi_ptr y, mpfi_srcptr x, double fact)
/* if c = mid(x) and r = rad(x), y = [c - (1+fact)*r , c + (1+fact)*r] */
{
  mp_prec_t prec;
  mpfr_t radius, factor;
  mpfr_t centre;
  int inexact_diam, inexact_div, inexact_factor, inexact_rad, inexact_centre;
  int inexact_left, inexact_right, inexact=0;

  if ( MPFI_NAN_P(x) ) {
    mpfr_set_nan(&(y->left));
    mpfr_set_nan(&(y->right));
    MPFR_RET_NAN;
  }

  prec=mpfi_get_prec(x);
  mpfr_init2(radius, prec);
  mpfr_init2(factor, prec);
  mpfr_init2(centre, prec);

  inexact_diam = mpfi_diam_abs(radius, x);
  inexact_div = mpfr_div_2exp(radius, radius, 1, GMP_RNDU);
  if (fact > 0.0)
    inexact_factor = mpfr_set_d(factor, (1.0+fact), GMP_RNDU);
  else
    inexact_factor = mpfr_set_d(factor, (1.0-fact), GMP_RNDU);
  inexact_rad = mpfr_mul(radius, radius, factor, GMP_RNDU);
  inexact_centre = mpfi_mid(centre, x);
  inexact_left = mpfr_sub(&(y->left), centre, radius, GMP_RNDD);
  inexact_right = mpfr_add(&(y->right), centre, radius, GMP_RNDU);

  mpfr_clear(radius);
  mpfr_clear(factor);
  mpfr_clear(centre);

  if ( MPFI_NAN_P(y) )
    MPFR_RET_NAN;

  if ( MPFI_LEFT_IS_INEXACT(inexact_diam)   || MPFI_LEFT_IS_INEXACT(inexact_div) ||
       MPFI_LEFT_IS_INEXACT(inexact_factor) || MPFI_LEFT_IS_INEXACT(inexact_rad) ||
       MPFI_LEFT_IS_INEXACT(inexact_centre) || inexact_left  )
    inexact += 1;
  if ( MPFI_RIGHT_IS_INEXACT(inexact_diam)   || MPFI_RIGHT_IS_INEXACT(inexact_div) ||
       MPFI_RIGHT_IS_INEXACT(inexact_factor) || MPFI_RIGHT_IS_INEXACT(inexact_rad) ||
       MPFI_RIGHT_IS_INEXACT(inexact_centre) || inexact_right  )
    inexact += 2;

  if (mpfi_revert_if_needed(y)) {
    fprintf(stderr, "Pb endpoints in reverse order in mpfi_blow: ");
    mpfi_out_str(stderr, 10, 0, y);
    fprintf(stderr, "\n");
    inexact = MPFI_REVERT_INEXACT_FLAGS(inexact);
  }

  return inexact;
}


/* Splits an interval into 2 halves */
int mpfi_bisect(mpfi_ptr y1, mpfi_ptr y2, mpfi_srcptr y)
{
  mp_prec_t prec, prec1, prec2;
  mpfr_t centre;
  int inexact_centre, dummy;

  if ( MPFI_NAN_P(y) ) {
    mpfr_set_nan(&(y1->left));
    mpfr_set_nan(&(y1->right));
    mpfr_set_nan(&(y2->left));
    mpfr_set_nan(&(y2->right));
    MPFR_RET_NAN;
  }
  else if ( !mpfi_bounded_p(y) ) {
    dummy = mpfi_set(y1, y);
    mpfr_set_nan(&(y2->left));
    mpfr_set_nan(&(y2->right));
    MPFR_RET_NAN;
  }

  prec=mpfi_get_prec(y);
  prec1=mpfi_get_prec(y1);
  prec2=mpfi_get_prec(y2);
  if ((prec1>=prec2) && (prec1>=prec))
    prec = prec1;
  else if ((prec2>=prec1) && (prec2>=prec))
    prec = prec2;
  mpfr_init2(centre, prec);

  inexact_centre = mpfi_mid(centre, y);

  dummy = mpfr_set(&(y1->left), &(y->left), MPFI_RNDD);
  dummy = mpfr_set(&(y2->right), &(y->right), MPFI_RNDU);
  dummy = mpfr_set(&(y1->right), centre, MPFI_RNDU);
  dummy = mpfr_set(&(y2->left), centre, MPFI_RNDD);

  mpfr_clear(centre);
  return inexact_centre;
}

char * mpfi_get_version()
{
  return(VERSION);
}

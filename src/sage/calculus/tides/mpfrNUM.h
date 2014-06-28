/****************************************************************************
 libTIDES. 
 This file is part of TIDES.
 
 Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
 TIDES is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TIDES is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with TIDES.  If not, see <http://www.gnu.org/licenses/>.
 
 *****************************************************************************/

/*  V-20  */

#ifndef _mpfrNUM_H
#define _mpfrNUM_H

#include <stdio.h>
#include <stdlib.h>
#include "mpfr.h"


/* DEFINITION OF ROUND TYPE*/
//#define TIDES_RND GMP_RNDN

int		binary_precision(int prec) ;
int		mpfrts_get_prec (void); 
void 	mpfrts_set_prec (int dig);

void 	mpfrts_set_rnd (mpfr_rnd_t trnd);

void	mpfrts_init (mpfr_t *rop); /* INITIALIZE ROP */
void	mpfrts_set_i(mpfr_t  *rop, long op); /* SET INTEGER NUMBER */
void	mpfrts_set_d (mpfr_t *rop, double op); /* SET DOUBLE NUMBER */
void	mpfrts_set_str (mpfr_t *rop, char *op); /* SET STRING AS NUMBER */
void	mpfrts_set (mpfr_t *rop, mpfr_t op); /* SET REALNUM AS NUMBER */
double	mpfrts_get_d (mpfr_t op); /* TRANSFORM TO DOUBLE */
long    mpfrts_get_i(mpfr_t op);


void	mpfrts_add (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 + OP2 */
void	mpfrts_add_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_sub (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 - OP2 */
void 	mpfrts_sub_i (mpfr_t *rop, mpfr_t op1, long int op2); 
void 	mpfrts_i_sub (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_mul (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 x OP2 */
void	mpfrts_mul_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_div (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 / OP2 */
void 	mpfrts_div_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_i_div (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_pow (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = OP1 ^ OP2 */
void	mpfrts_pow_i (mpfr_t *rop, mpfr_t op1, long int op2);
void 	mpfrts_i_pow (mpfr_t *rop, unsigned long int op1, mpfr_t op2);
void	mpfrts_abs(mpfr_t  *rop, mpfr_t op); /* ROP = |OP| */
void	mpfrts_neg(mpfr_t  *rop, mpfr_t op);

int		mpfrts_greater(mpfr_t op1, mpfr_t op2); /* 1 IF OP1 > OP2; 0 OTHERWHISE */
int		mpfrts_greaterequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_less(mpfr_t op1, mpfr_t op2); 
int		mpfrts_lessequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_equal(mpfr_t op1, mpfr_t op2); /* 1 IF OP1 == OP2; 0 OTHERWHISE */

void	mpfrts_sqrt(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log(mpfr_t  *rop, mpfr_t op); /* ROP = LOG (OP) */
void	mpfrts_log2(mpfr_t  *rop, mpfr_t op); /* ROP = LOG_2 (OP) */
void	mpfrts_log10(mpfr_t  *rop, mpfr_t op); /* ROP = LOG_1(OP) */
void	mpfrts_exp(mpfr_t  *rop, mpfr_t op); /* ROP = e ^ OP */
void	mpfrts_exp2(mpfr_t  *rop, mpfr_t op); /* ROP = 2 ^ OP */
void	mpfrts_exp10(mpfr_t  *rop, mpfr_t op); /* ROP = 10 ^ OP */
void	mpfrts_cos(mpfr_t  *rop, mpfr_t op); /* ROP = COS (OP) */
void	mpfrts_sin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sin_cos(mpfr_t *rsin, mpfr_t *rcos, mpfr_t op); /* RSIN = SIN (OP); RCOS = COS (OP) */
void	mpfrts_tan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sec(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csc(mpfr_t  *rop, mpfr_t op);
void	mpfrts_cot(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acos(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan2(mpfr_t  *rop, mpfr_t op1, mpfr_t op2); /* ROP = ATAN2 (OP1, OP2), THE SAME AS IN DOUBLE */
void	mpfrts_cosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_tanh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sech(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csch(mpfr_t  *rop, mpfr_t op);
void	mpfrts_coth(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atanh(mpfr_t  *rop, mpfr_t op);


void	mpfrts_write_var(mpfr_t op);
void	mpfrts_write (char *c, mpfr_t op);
void	mpfrts_fread (FILE *file, mpfr_t rop);
void	mpfrts_fwrite (FILE *file, mpfr_t op, int prec); 

typedef  mpfr_t*	Array1MP;
typedef  mpfr_t**	Array2MP;

void	Array1MP_init(Array1MP *vec, long dim);
void	Array1MP_clear(Array1MP *vec, long dim);
void	Array2MP_init(Array2MP *vec, long rows, long columns);
void	Array2MP_clear(Array2MP *vec, long rows, long columns);
void	Array1MP_set(Array1MP rop, Array1MP op, long dim);
void	Array2MP_set(Array2MP rop, Array2MP op, long rows, long columns);
void	Array2MP_column_set(Array2MP rop, Array1MP op, long c, long dim);
void	Array2MP_row_set(Array2MP rop, Array1MP op, long r, long dim);


#endif



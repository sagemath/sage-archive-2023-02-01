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

#ifndef _doubleNUM_H
#define _doubleNUM_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void	double_init(double *rop); /* INITIALIZE */
void	double_set_d(double *rop, double op); /* SET DOUBLE NUMBER */
void	double_set_str(double *rop, char *op); /* SET STRING */
void	double_set(double *rop, double op); /* SET REALNUM */
double	double_get_d(double op); /* TRANSFORM TO DOUBLE */

void	double_clear (double op);


void	double_add(double *rop, double op1, double op2); /* ROP = OP1 + OP2 */
void	double_sub(double *rop, double op1, double op2); /* ROP = OP1 - OP2 */
void	double_mul(double *rop, double op1, double op2); /* ROP = OP1 x OP2 */
void	double_div(double *rop, double op1, double op2); /* ROP = OP1 / OP2 */
void	double_pow(double *rop, double op1, double op2); /* ROP = OP1 ^ OP2 */
void	double_abs(double *rop, double op);  /* ROP = |OP| */

void	double_add_i (double *rop, double op1, long   op2); /* ROP = OP1 + OP2 */
void	double_sub_i (double *rop, double op1, long   op2); /* ROP = OP1 - OP2 */
void	double_i_sub (double *rop, long   op1, double op2);
void	double_mul_i (double *rop, double op1, long   op2); /* ROP = OP1 x OP2 */
void	double_div_i (double *rop, double op1, long  op2); /* ROP = OP1 / OP2 */
void 	double_i_div (double *rop, long   op1, double op2); 
void 	double_pow_i (double *rop, double op1, long   op2);
void	double_i_pow (double *rop, unsigned long  op1, double op2);

int		double_greater(double op1, double op2); /* 1 IF OP1 > OP2; 0 OTHERWHISE */
int		double_greaterequal(double op1, double op2);
int		double_less(double op1, double op2); 
int		double_lessequal(double op1, double op2);
int		double_equal(double op1, double op2); /* 1 IF OP1 == OP2; 0 OTHERWHISE */

void	double_log(double *rop, double op); /* ROP = LOG (OP) */
void	double_log10(double *rop, double op); /* ROP = LOG_1(OP) */
void	double_exp(double *rop, double op); /* ROP = e ^OP */
void	double_exp2(double *rop, double op); /* ROP = 2 ^ OP */
void	double_exp10(double *rop, double op); /* ROP = 10 ^ OP */
void	double_cos(double *rop, double op); /* ROP = COS (OP) */
void	double_sin(double *rop, double op);
void	double_sin_cos(double *rsin, double *rcos, double op); /* RSIN = SIN (OP); RCOS = COS (OP) */
void	double_tan(double *rop, double op);
void	double_sec(double *rop, double op);
void	double_csc(double *rop, double op);
void	double_cot(double *rop, double op);
void	double_acos(double *rop, double op);
void	double_asin(double *rop, double op);
void	double_atan(double *rop, double op);
void	double_atan2(double *rop, double op1, double op2); /* ROP = ATAN2 (OP1, OP2); SAME AS IN DOUBLE */
void	double_cosh(double *rop, double op);
void	double_sinh(double *rop, double op);
void	double_tanh(double *rop, double op);
void	double_sech(double *rop, double op);
void	double_csch(double *rop, double op);
void	double_coth(double *rop, double op);
void	double_acosh(double *rop, double op);
void	double_asinh(double *rop, double op);
void	double_atanh(double *rop, double op);

typedef  double*	Array1DB;
typedef  double**	Array2DB;

void	Array1DB_init(Array1DB *vec, long dim); /* Initialize 1-dimensional arrays */
void	Array1DB_clear(Array1DB *vec, long dim);
void	Array2DB_init(Array2DB *vec, long rows, long columns); /* Initialize 2-dimensional arrays */
void	Array2DB_clear(Array2DB *vec, long rows, long columns);
void	Array3DB_init(Array2DB *vec, long dim, long rows, long columns); /* Initialize one array of 2-dimensional arrays */

void	Array1DB_set(Array1DB rop, Array1DB op, long dim);
void	Array2DB_set(Array2DB rop, Array2DB op, long rows, long columns);
void	Array2DB_column_set(Array2DB rop, Array1DB op, long c, long dim);
void	Array2DB_row_set(Array2DB rop, Array1DB op, long r, long dim);


void	double_write (char *c, double op);


#endif


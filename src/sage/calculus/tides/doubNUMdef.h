/****************************************************************************
 libTIDES. Version 1.0.0.
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


#ifndef _doubNUMdef_H
#define _doubNUMdef_H
#ifdef real_Double
	typedef double		realNUM;
	typedef double* 	realVEC;
	typedef double**	realMAT;

	#define  tides_init		double_init
	#define  tides_set_d	double_set_d
	#define  tides_set_str  double_set_str
	#define  tides_set		double_set
	#define  tides_get_d	double_get_d
	#define	 tides_clear	double_clear

	#define  tides_add  double_add
	#define  tides_add_i double_add_i
	#define  tides_sub  double_sub
	#define  tides_sub_i double_sub_i
	#define  tides_i_sub double_i_sub
	#define  tides_mul  double_mul
	#define  tides_mul_i double_mul_i
	#define  tides_div  double_div
	#define  tides_div_i double_div_i
	#define  tides_i_div double_i_div
	#define  tides_pow  double_pow
	#define  tides_pow_i double_pow_i
	#define  tides_i_pow double_i_pow
	#define  tides_abs  double_abs

	#define  tides_greater			double_greater
	#define  tides_greaterequal		double_greaterequal
	#define  tides_less				double_less
	#define  tides_lessequal		double_lessequal
	#define  tides_equal			double_equal

	#define  tides_log		double_log
	#define  tides_log10	double_log10
	#define  tides_exp		double_exp
	#define  tides_exp2		double_exp2
	#define  tides_exp10	double_exp10
	#define  tides_cos		double_cos
	#define  tides_sin		double_sin
	#define  tides_tan		double_tan
	#define  tides_sin_cos  double_sin_cos
	#define  tides_sec		double_sec
	#define  tides_csc		double_csc
	#define  tides_cot		double_cot
	#define  tides_acos		double_acos
	#define  tides_asin		double_asin
	#define  tides_atan		double_atan
	#define  tides_atan2	double_atan2
	#define  tides_cosh		double_cosh
	#define  tides_sinh		double_sinh
	#define  tides_tanh		double_tanh
	#define  tides_sech		double_sech
	#define  tides_csch		double_csch
	#define  tides_coth		double_coth
	#define  tides_acosh	double_acosh
	#define  tides_asinh	double_asinh
	#define  tides_atanh	double_atanh


	#define  realMAT_init		Array2DB_init
	#define  realVEC_init		Array1DB_init

	#define  variables_init		varDB_init
	#define  parameters_init	parDB_init
	#define  links_init			linkDB_init
	#define  derivatives_init	derDB_init
	#define  write_solution		write_solution_DB

	#define  variables_free		varDB_free
	#define  parameters_free	parDB_free
	#define  links_free			linkDB_free

	#define  set_precision_digits	double_set_prec

	#define  var_t				double_var_t
	#define  var_t_c			double_var_t_c
	#define  add_t				double_add_t
	#define  add_t_c			double_add_t_c
	#define  sub_t				double_sub_t
	#define  sub_t_c			double_sub_t_c
	#define  mul_t				double_mul_t
	#define  mul_t_c			double_mul_t_c
	#define  divide_t			double_div_t

	#define  inv_t				double_inv_t
	#define  exp_t				double_exp_t
	#define  pow_t_c			double_pow_t_c
	#define  sin_t				double_sin_t
	#define  sincos_t			double_sin_cos_t
	#define  sincosh_t			double_sinh_cosh_t
	#define  cos_t				double_cos_t
	#define  sinh_t				double_sinh_t
	#define  cosh_t				double_cosh_t

	#define  asin_t				double_asin_t
	#define  acos_t				double_acos_t
	#define  atan_t				double_atan_t
	#define  asinh_t			double_asinh_t
	#define  acosh_t			double_acosh_t
	#define  atanh_t			double_atanh_t
	#define  log_t				double_log_t

	#define tides_write			double_write

#endif
#endif

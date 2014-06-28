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


#ifndef _mpfrNUMdef_H
#define _mpfrNUMdef_H
#ifdef real_MP
	typedef mpfr_t 		realNUM;
	typedef mpfr_t* 	realVEC;
	typedef mpfr_t**	realMAT;

	#define  tides_init		mpfrts_init
	#define  tides_set_d	mpfrts_set_d
	#define  tides_set_str	mpfrts_set_str
	#define  tides_set		mpfrts_set
	#define  tides_get_d	mpfrts_get_d
	#define	 tides_clear	mpfr_clear

	#define	 tides_add	mpfrts_add
	#define  tides_add_i	mpfrts_add_i
	#define  tides_sub	mpfrts_sub
	#define  tides_sub_i	mpfrts_sub_i
	#define  tides_i_sub	mpfrts_i_sub
	#define  tides_mul	mpfrts_mul
	#define  tides_mul_i	mpfrts_mul_i
	#define  tides_div	mpfrts_div
	#define	 tides_div_i	mpfrts_div_i
	#define  tides_i_div	mpfrts_i_div
	#define  tides_pow	mpfrts_pow
	#define  tides_pow_i	mpfrts_pow_i
	#define	 tides_i_pow	mpfrts_i_pow
	#define  tides_abs	mpfrts_abs

	#define  tides_greater		mpfrts_greater           
	#define  tides_greaterequal mpfrts_greaterequal
	#define  tides_less			mpfrts_less                 
	#define  tides_lessequal	mpfrts_lessequal 	     
	#define  tides_equal		mpfrts_equal            

	#define  tides_log		mpfrts_log
	#define  tides_log2		mpfrts_log2
	#define  tides_log10		mpfrts_log10
	#define  tides_exp		mpfrts_exp
	#define  tides_exp2		mpfrts_exp2
	#define  tides_exp10		mpfrts_exp10
	#define  tides_cos		mpfrts_cos
	#define  tides_sin		mpfrts_sin
	#define  tides_sin_cos		mpfrts_sin_cos
	#define  tides_tan		mpfrts_tan
	#define  tides_sec		mpfrts_sec
	#define  tides_csc		mpfrts_csc
	#define  tides_cot		mpfrts_cot
	#define  tides_acos		mpfrts_acos
	#define  tides_asin		mpfrts_asin
	#define  tides_atan		mpfrts_atan
	#define  tides_atan2	mpfrts_atan2
	#define  tides_cosh		mpfrts_cosh
	#define  tides_sinh		mpfrts_sinh
	#define  tides_tanh		mpfrts_tanh
	#define  tides_sech		mpfrts_sech
	#define  tides_csch		mpfrts_csch
	#define  tides_coth		mpfrts_coth
	#define  tides_acosh	mpfrts_acosh
	#define  tides_asinh	mpfrts_asinh
	#define  tides_atanh	mpfrts_atanh

	#define tides_write	mpfrts_write


	#define  realMAT_init		Array2MP_init
	#define	 realVEC_init		Array1MP_init

	#define  variables_init		varMP_init
	#define  parameters_init	parMP_init
	#define  links_init		linkMP_init
	#define  derivatives_init	derMP_init

	#define  variables_free		varMP_free
	#define  parameters_free	parMP_free
	#define  links_free		linkMP_free
	#define  write_solution		write_solution_MP

	#define  set_precision_digits	mpfrts_set_prec

	#define  var_t				mpfrts_var_t
	#define  var_t_c			mpfrts_var_t_c
	#define  add_t				mpfrts_add_t
	#define  add_t_c			mpfrts_add_t_c
	#define  sub_t				mpfrts_sub_t
	#define  sub_t_c			mpfrts_sub_t_c
	#define  mul_t				mpfrts_mul_t
	#define  mul_t_c			mpfrts_mul_t_c
	#define  divide_t			mpfrts_div_t

	#define  inv_t				mpfrts_inv_t
	#define  exp_t				mpfrts_exp_t
	#define  pow_t_c			mpfrts_pow_t_c
	#define  sincos_t			mpfrts_sin_cos_t
	#define  sincosh_t			mpfrts_sinh_cosh_t
	#define  sin_t				mpfrts_sin_t
	#define  cos_t				mpfrts_cos_t
	#define  sinh_t				mpfrts_sinh_t
	#define  cosh_t				mpfrts_cosh_t

	#define  asin_t				mpfrts_asin_t
	#define  acos_t				mpfrts_acos_t
	#define  atan_t				mpfrts_atan_t
	#define  asinh_t			mpfrts_asinh_t
	#define  acosh_t			mpfrts_acosh_t
	#define  atanh_t			mpfrts_atanh_t
	#define  log_t				mpfrts_log_t

#endif

#endif

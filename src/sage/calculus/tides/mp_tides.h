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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mpfr.h"

#ifndef Header_DP_TIDES_h
#define Header_DP_TIDES_h

//Extern

extern int		_info_steps_taylor_;
extern int		num_etapas;
extern int		order_series;
extern int		order_estimator;
extern double   fac1;
extern double   fac2;
extern double   fac3;
extern double   rmaxstep; 
extern double   rminstep; 
extern int      nitermax; 
extern int      nordinc; 
extern int      minord; 
extern int      defect_error_control;
extern int      kahan_summation;
extern int      compensated_horner;

extern int		test_relminstep ;

extern int		printError ;

//typedef
typedef struct it_dt{
	int		MAX_ORDER;
	long	NDER;
	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
	int		*PARTIAL_LIST, *FUNCTION_LIST;
	long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
	long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;
	long	LINK_ELEMENTS;
	int		clearPartials;
} iteration_data;


//CommonITER
int		same_der(int dim, int *der1, int *der2);
void	string_to_der(int dim, char *sder, int *der);
long	position_derivative(char* sder, int *pdd);
long	position_variable(int v, char* der, int NVARS, int NFUNS, int *pdd);
long	position_function(int f, char* der, int NVARS, int NFUNS, int *pdd);
int		is_variable(iteration_data *itd, int num);
void	set_iteration_parameters(iteration_data *itd, int v, int p, int f, int o, int *pdd);
void	delete_iteration_parameters(iteration_data *itd);
void	set_links(iteration_data *itd, int l, int *flst);
void	check_iteration_data_parameters(int t, int pa, int pb);
void	set_info_error(void);
void	unset_info_error(void);

void	check_relminstep(void);
void	uncheck_relminstep(void);
void	use_default_step_estimator(void);
void	set_maxnumsteps(unsigned long val);


#endif



#ifndef Header_MP_TIDES_h
#define Header_MP_TIDES_h

//Extern

extern mpfr_t   mp_relminstep;

extern mpfr_t	etapa_mp_minima;
extern mpfr_t	etapa_mp_maxima;
extern mpfr_t	etapa_mp_total;

extern mpfr_rnd_t	TIDES_RND;
extern int		TIDES_PREC;
extern int		BINARY_PRECISION;
extern int		DECIMAL_PRECISION;


//Defines
#define initialize_mp_case() \
long  NUM_COLUMNS;\
set_links(itd, LINKS, POS_FUNCTIONS);\
NUM_COLUMNS = (itd->NVARS+itd->NFUNS)*itd->NDER;\
if(ORDER < 0) return NUM_COLUMNS;\
check_iteration_data_parameters(0, itd->NVARS, VARIABLES);\
check_iteration_data_parameters(1, itd->NPARS, PARAMETERS);\
check_iteration_data_parameters(2, itd->NFUNS, FUNCTIONS);\
itd->LINK_ELEMENTS = itd->NDER*itd->MAX_ORDER;\
mpfr_t var[itd->NVARS+1][itd->LINK_ELEMENTS];\
mpfr_t par[itd->NPARS][itd->LINK_ELEMENTS];\
mpfr_t link[LINKS][itd->LINK_ELEMENTS];\
varMP_init(itd,(mpfr_t*)var,v,t);\
parMP_init(itd,(mpfr_t*)par,p);\
linkMP_init(itd,(mpfr_t*)link);\
derMP_init(itd,(mpfr_t*)var, (mpfr_t*)par, v);

#define write_mp_solution()	write_sol_MP(itd,cvfd,(mpfr_t*)var,(mpfr_t*)link);

#define clear_vpl()	\
clear_mpfr_vpl(itd, (mpfr_t*)var, (mpfr_t*)par, (mpfr_t*)link);

#define clear_cts()	\
for(i = 0; i < NCONST ; i++ ) mpfr_clear(ct[i]);

//typedef

typedef long (*MPLinkedFunction)(iteration_data *itd, mpfr_t t, mpfr_t *v, mpfr_t *p, int orden, mpfr_t *coef);

typedef struct mp_DM {
	int rows;
	int columns;
	mpfr_t **data;
} mp_data_matrix;


typedef struct mp_node_event {
	mpfr_t *data;
	struct mp_node_event *next;
} mp_event_node;

typedef struct mp_linked_events_nodes {
	int total;
	int dim;
	mp_event_node *first;
} mp_event_list;

typedef  mpfr_t*	Array1MP;
typedef  mpfr_t**	Array2MP;

//doubNUMdef
#define TIDES_RND GMP_RNDN

#define  set_precision_digits	mpfrts_set_prec
#define  set_rounding_mode  	mpfrts_set_rnd
	


//mpfrITER		
void	varMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *v, mpfr_t t);
void	parMP_init(iteration_data *itd, mpfr_t *par, mpfr_t *p);
void	linkMP_init(iteration_data *itd, mpfr_t *lk);
void	derMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *par, mpfr_t *v);
void 	clear_mpfr_vpl (iteration_data *itd, mpfr_t *var, mpfr_t *par, mpfr_t *link);
void	write_sol_MP(iteration_data *itd, mpfr_t *cvfd, mpfr_t *var, mpfr_t *link);
void	mpfrts_htilde(iteration_data *itd, mpfr_t *h, long j, long v, long i, mpfr_t *ht, int ORDER_INDEX);
void	mpfrts_var_t(iteration_data *itd, mpfr_t *f, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_var_t_c(iteration_data *itd, char* cs, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_var_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_div_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_div_t_vc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_div_t_cv(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_inv_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_exp_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_pow_t_c(iteration_data *itd, mpfr_t *u, char* cs, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_pow_t_cc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sct_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i);
void	mpfrts_sct_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX);
void	mpfrts_sin_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_cos_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_sin_cos_t (iteration_data *itd, mpfr_t *f, mpfr_t *s, mpfr_t *c, int ORDER_INDEX);
void	mpfrts_sinh_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_cosh_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_sinh_cosh_t (iteration_data *itd, mpfr_t *f,  mpfr_t *s, mpfr_t *c, int ORDER_INDEX);
void	mpfrts_fgt_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i);
void	mpfrts_fgt_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX);
void	mpfrts_log_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_asin_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_acos_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_atan_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_asinh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_acosh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_atanh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);



//mpfrNUM
int		binary_precision(int prec) ;
void	mpfrts_init (mpfr_t *rop); 
void	mpfrts_set_i(mpfr_t  *rop, long op); 
void	mpfrts_set_d (mpfr_t *rop, double op); 
void	mpfrts_set_str (mpfr_t *rop, char *op); 
void	mpfrts_set (mpfr_t *rop, mpfr_t op); 
double	mpfrts_get_d (mpfr_t op); 
long    mpfrts_get_i(mpfr_t op);
int		mpfrts_get_prec (void); 
void 	mpfrts_set_prec (int dig);
void	mpfrts_add (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_add_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_sub (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void 	mpfrts_sub_i (mpfr_t *rop, mpfr_t op1, long int op2); 
void 	mpfrts_i_sub (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_mul (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_mul_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_div (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void 	mpfrts_div_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_i_div (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_pow (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_pow_i (mpfr_t *rop, mpfr_t op1, long int op2);
void 	mpfrts_i_pow (mpfr_t *rop, unsigned long int op1, mpfr_t op2);
void	mpfrts_abs(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_neg(mpfr_t  *rop, mpfr_t op);
int		mpfrts_greater(mpfr_t op1, mpfr_t op2); 
int		mpfrts_greaterequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_less(mpfr_t op1, mpfr_t op2); 
int		mpfrts_lessequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_equal(mpfr_t op1, mpfr_t op2); 
void	mpfrts_sqrt(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log2(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log10(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp2(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp10(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_cos(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_sin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sin_cos(mpfr_t *rsin, mpfr_t *rcos, mpfr_t op); 
void	mpfrts_tan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sec(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csc(mpfr_t  *rop, mpfr_t op);
void	mpfrts_cot(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acos(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan2(mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_cosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_tanh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sech(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csch(mpfr_t  *rop, mpfr_t op);
void	mpfrts_coth(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atanh(mpfr_t  *rop, mpfr_t op);
void 	mpfrts_write_var(mpfr_t op);
void 	mpfrts_write (char *c, mpfr_t op);
void 	mpfrts_fread (FILE *file, mpfr_t rop);
void 	mpfrts_fwrite (FILE *file, mpfr_t op, int prec); 
void	Array1MP_init(Array1MP *vec, long dim);
void	Array1MP_clear(Array1MP *vec, long dim);
void	Array2MP_init(Array2MP *vec, long rows, long columns);
void	Array2MP_clear(Array2MP *vec, long rows, long columns);
void	Array1MP_set(Array1MP rop, Array1MP op, long dim);
void	Array2MP_set(Array2MP rop, Array2MP op, long rows, long columns);
void	Array2MP_column_set(Array2MP rop, Array1MP op, long c, long dim);
void	Array2MP_row_set(Array2MP rop, Array1MP op, long r, long dim);

//mpfrPOL
void mp_pol_derivative(mpfr_t *pol, int grado, mpfr_t *dpol);
void mp_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_horner_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_classic_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_twosum(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y);
void mp_split_real(mpfr_t a, mpfr_t *ah, mpfr_t *al);
void mp_twoproduct(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y);
void mp_compensated_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_ch_evaluate_poly_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *vpol, mpfr_t *vder);
int  mp_one_zero(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz);
int  mp_one_extremum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_one_maximum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_one_minimum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_brent(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz);
int  mp_rtsafe(mpfr_t *pol, int grado, mpfr_t x1, mpfr_t x2, mpfr_t tol, mpfr_t *raiz);

//mpfrTODE
void mp_set_relminstep(mpfr_t val);
void mp_set_info_taylor(void);
void mp_unset_info_taylor(void);
void mp_str_info_taylor(void);
void mp_add_info_step(mpfr_t tstep);
int  mp_taylor_order(mpfr_t eps);
void mp_norm_inf_vec(mpfr_t *rop, int nvar, mpfr_t *y) ;
void mp_norm_inf(mpfr_t *rop, int n, int k, mpfr_t *coef, int MAX_ORDER);
void mp_compute_step(mpfr_t *rop, double *hant, mpfr_t tol, int n, int ord, mpfr_t *coef, int MAX_ORDER);
void mp_compute_tol (mpfr_t *tol, mpfr_t tolrel, mpfr_t tolabs, int nvar, mpfr_t *y); 
void mp_taylor_horner(int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t *x, int MAX_ORDER);
void mp_taylor_horner_der(int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t *x, int MAX_ORDER);
void mp_write_taylor_solution( int n, int j, mpfr_t tini, mpfr_t *x,  mpfr_t*** mat, FILE* fileout);	
int  mp_valid_step (MPLinkedFunction fcn, iteration_data *itd, mpfr_t *step, mpfr_t tip, mpfr_t eps, int nvar, int ncol, int order, mpfr_t *cvfd, mpfr_t *p, int MAX_ORDER);
int  mp_tides(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t *lt, int ntes, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_point(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t tf, mpfr_t dt, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_list(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t *lt, int ntes, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout) ;
int  mp_tides_delta(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t dt, int ntot, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_delta_ft(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t tf, mpfr_t dt, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_kernel(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t dt,  mpfr_t *dtl, int ntot, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout, mp_data_matrix *der);
long mp_number_of_columns(MPLinkedFunction fcn);
void init_mp_data_matrix(mp_data_matrix *dm, int r, int c);
void delete_mp_data_matrix(mp_data_matrix *dm);
int  getOrder (void);
int  getNsteps (void);

//mpfrEVENTS
void mp_init_event_list(mp_event_list *lista, int ncol);
void mp_add_to_event_list(mp_event_list *lista, mpfr_t *val);
void mp_event_list_to_array(mp_event_list lista, mp_data_matrix *array);

int  mp_tides_find_zeros(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout) ;
int  mp_tides_find_extrema(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout);
int  mp_tides_find_minimum(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout);
int  mp_tides_find_maximum(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout);
int  mp_tides_events(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout,int evcase);

#endif








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
extern double   dp_relminstep;


extern double	etapa_dp_minima;
extern double	etapa_dp_maxima;
extern double	etapa_dp_total;

extern int		printError ;

//Defines
#define initialize_dp_case() \
long  NUM_COLUMNS;\
set_links(itd, LINKS, POS_FUNCTIONS);\
NUM_COLUMNS = (itd->NVARS+itd->NFUNS)*itd->NDER;\
if(ORDER < 0) return NUM_COLUMNS;\
check_iteration_data_parameters(0, itd->NVARS, VARIABLES);\
check_iteration_data_parameters(1, itd->NPARS, PARAMETERS);\
check_iteration_data_parameters(2, itd->NFUNS, FUNCTIONS);\
itd->LINK_ELEMENTS = itd->NDER*itd->MAX_ORDER;\
double var[itd->NVARS+1][itd->LINK_ELEMENTS];\
double par[itd->NPARS][itd->LINK_ELEMENTS];\
double link[LINKS][itd->LINK_ELEMENTS];\
varDB_init(itd,(double*)var,v,t);\
parDB_init(itd,(double*)par,p);\
linkDB_init(itd,(double*)link);\
derDB_init(itd,(double*)var, (double*)par, v);

#define write_dp_solution()	write_sol_DB(itd,cvfd,(double*)var,(double*)link);

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

typedef long (*DBLinkedFunction)(iteration_data *itd, double t, double *v, double *p, int orden, double *cvfd);

typedef struct dp_DM {
	int rows;
	int columns;
	double **data;
} dp_data_matrix;


typedef struct dp_node_event {
	double *data;
	struct dp_node_event *next;
} dp_event_node;

typedef struct dp_linked_events_nodes {
	int total;
	int dim;
	dp_event_node *first;
} dp_event_list;

typedef  double*	Array1DB;
typedef  double**	Array2DB;

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
void	set_maxnumsteps(unsigned long val);

//doubITER
void	varDB_init(iteration_data *itd, double *var,  double *v, double t);
void	parDB_init(iteration_data *itd, double *par, double *p);
void	linkDB_init(iteration_data *itd, double *lk);
void	derDB_init(iteration_data *itd, double *var, double *par, double *v);
void	write_sol_DB(iteration_data *itd, double *cvfd, double *var, double *link);
void	double_htilde(iteration_data *itd, double *h, long j, long v, long i, double *ht, int ORDER_INDEX);
void	double_var_t(iteration_data *itd, double *f, double *u, int ORDER_INDEX);
void	double_var_t_c(iteration_data *itd, char* cs, double *w, int ORDER_INDEX);
void	double_var_t_cc(iteration_data *itd, double c, double *w, int ORDER_INDEX);
void	double_add_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_sub_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_add_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_sub_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_add_t_cc(iteration_data *itd, double c, double *u,  double *w, int ORDER_INDEX);
void	double_sub_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);
void	double_mul_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_mul_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_mul_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);
void	double_div_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_div_t_vc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX);
void	double_div_t_cv(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);
void	double_inv_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_exp_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_pow_t_c(iteration_data *itd, double *u, char* cs, double *w, int ORDER_INDEX);
void	double_pow_t_cc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX);
void	double_sct_0(iteration_data *itd, double *f, double *g, double *h, long i);
void	double_sct_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX);
void	double_sin_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX);
void	double_cos_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX);
void    double_sin_cos_t (iteration_data *itd, double *f, double *s, double *c,  int ORDER_INDEX);
void	double_sinh_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX);
void	double_cosh_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX);
void    double_sinh_cosh_t (iteration_data *itd, double *f, double *s, double *c, int ORDER_INDEX);
void	double_fgt_0(iteration_data *itd, double *f, double *g, double *h, long i);
void	double_fgt_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX);
void	double_log_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_asin_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_acos_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_atan_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_asinh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_acosh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_atanh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);

//doubNUM
void	double_init(double *rop); 
void	double_set_d(double *rop, double op); 
void	double_set_str(double *rop, char *op); 
void	double_set(double *rop, double op); 
double	double_get_d(double op); 
void	double_clear (double op);
void	double_add(double *rop, double op1, double op2); 
void	double_sub(double *rop, double op1, double op2); 
void	double_mul(double *rop, double op1, double op2); 
void	double_div(double *rop, double op1, double op2); 
void	double_pow(double *rop, double op1, double op2); 
void	double_abs(double *rop, double op);  
void	double_add_i (double *rop, double op1, long   op2); 
void	double_sub_i (double *rop, double op1, long   op2); 
void	double_i_sub (double *rop, long   op1, double op2);
void	double_mul_i (double *rop, double op1, long   op2); 
void	double_div_i (double *rop, double op1, long  op2); 
void 	double_i_div (double *rop, long   op1, double op2); 
void 	double_pow_i (double *rop, double op1, long   op2);
void	double_i_pow (double *rop, unsigned long  op1, double op2);
int		double_greater(double op1, double op2); 
int		double_greaterequal(double op1, double op2);
int		double_less(double op1, double op2); 
int		double_lessequal(double op1, double op2);
int		double_equal(double op1, double op2); 
void	double_log(double *rop, double op); 
void	double_log10(double *rop, double op); 
void	double_exp(double *rop, double op); 
void	double_exp2(double *rop, double op); 
void	double_exp10(double *rop, double op); 
void	double_cos(double *rop, double op); 
void	double_sin(double *rop, double op);
void	double_sin_cos(double *rsin, double *rcos, double op); 
void	double_tan(double *rop, double op);
void	double_sec(double *rop, double op);
void	double_csc(double *rop, double op);
void	double_cot(double *rop, double op);
void	double_acos(double *rop, double op);
void	double_asin(double *rop, double op);
void	double_atan(double *rop, double op);
void	double_atan2(double *rop, double op1, double op2); 
void	double_cosh(double *rop, double op);
void	double_sinh(double *rop, double op);
void	double_tanh(double *rop, double op);
void	double_sech(double *rop, double op);
void	double_csch(double *rop, double op);
void	double_coth(double *rop, double op);
void	double_acosh(double *rop, double op);
void	double_asinh(double *rop, double op);
void	double_atanh(double *rop, double op);
void	Array1DB_init(Array1DB *vec, long dim); 
void	Array1DB_clear(Array1DB *vec, long dim);
void	Array2DB_init(Array2DB *vec, long rows, long columns); 
void	Array2DB_clear(Array2DB *vec, long rows, long columns);
void	Array3DB_init(Array2DB *vec, long dim, long rows, long columns); 
void	Array1DB_set(Array1DB rop, Array1DB op, long dim);
void	Array2DB_set(Array2DB rop, Array2DB op, long rows, long columns);
void	Array2DB_column_set(Array2DB rop, Array1DB op, long c, long dim);
void	Array2DB_row_set(Array2DB rop, Array1DB op, long r, long dim);
void	double_write (char *c, double op);

//doubPOL
void dp_pol_derivative(double *pol, int grado, double *dpol);
void dp_horner(double *pol, int grado, double t, double *eval);
void dp_horner_der(double *pol, int grado, double t, double *eval);
void dp_classic_horner(double *pol, int grado, double t, double *eval);
void dp_twosum(double a, double b, double *x, double *y);
void dp_split_real(double a, double *ah, double *al);
void dp_twoproduct(double a, double b, double *x, double *y);
void dp_compensated_horner(double *pol, int grado, double t, double *eval);
void dp_ch_evaluate_poly_der(double *pol, int grado, double t, double *vpol, double *vder);
int  dp_one_zero(double *pol, int grado, double ta, double tb, double tol, double *raiz);
int  dp_one_extremum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_one_maximum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_one_minimum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_brent(double *pol, int grado, double ta, double tb, double tol, double *raiz);
int  dp_rtsafe(double *pol, int grado, double x1, double x2, double tol, double *raiz);

//doubTODE
void	check_relminstep(void);
void	uncheck_relminstep(void);
void	dp_set_relminstep(double val);

void use_default_step_estimator(void);
void dp_set_info_taylor(void);
void dp_unset_info_taylor(void);
void dp_str_info_taylor(void);
void dp_add_info_step(double tstep);
int  dp_taylor_order(double eps);
void dp_norm_inf_var(double *rop, int nvar, double *y) ;
void dp_norm_inf(double *rop, int n, int k, double *coef, int MAX_ORDER);
void dp_compute_step(double *rop, double *hant, double tol, int n, int ord, double *cvfd, int MAX_ORDER);
void dp_compute_tol (double *tol, double tolrel, double tolabs, int nvar, double *y) ;
void dp_taylor_horner(int n, int ord, double *cvfd, double t, double *x, int MAX_ORDER);
void dp_taylor_horner_der (int n, int ord, double *cvfd, double t, double *x, int MAX_ORDER);
void dp_write_taylor_solution(int n,  int j, double tini, double *x, double ***mat, FILE* fileout);
int  dp_valid_step (DBLinkedFunction fcn,  iteration_data *itd, double *step, double tip,double tol, int nvar, int ncol, int order, double *cvfd,double *p, int MAX_ORDER);
int  dp_tides(DBLinkedFunction fcn,int *pdd, int nvar,  int npar, int nfun, double *x, double *p, double *lt, int ntes, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_point(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double tf, double dt, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_list(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double *lt, int ntes, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout) ;
int  dp_tides_delta(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double dt, int ntot, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_delta_ft(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double tf, double dt, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_kernel(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double dt, double *dtl, int ntot, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout, dp_data_matrix *der) ;
long dp_number_of_columns(DBLinkedFunction fcn);
void init_dp_data_matrix(dp_data_matrix *dm, int r, int c);
void delete_dp_data_matrix(dp_data_matrix *dm);

//doubEVENTS
void dp_init_event_list(dp_event_list *lista, int ncol);
void dp_add_to_event_list(dp_event_list *lista, double *val);
void dp_event_list_to_array(dp_event_list lista, dp_data_matrix *array);
int  dp_tides_find_zeros(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout) ;
int  dp_tides_find_extrema(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout) ;
int  dp_tides_find_minimum(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout) ;
int  dp_tides_find_maximum(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout);
int  dp_tides_events(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout, int evcase) ;


#endif





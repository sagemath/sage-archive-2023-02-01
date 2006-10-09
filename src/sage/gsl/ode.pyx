include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'

include 'gsl.pxi'

cdef struct param_container:
   int param_n
   int y_n
   double *params



ode_function = None
ode_jacobian= None

cdef int c_jac(double t,double *y,double *dfdy,double *dfdt,void *params):
   """This is a C function designed to wrap a user defined python function representing the jacabian of a function.The global variable ode_jacobian is set to reference the user defined functio   and that function is called from within thisc function. The user defined jacobian should accept a time value a list of y_values and a list of paramters. It should return a two component list or tuple whose first component is a dictionary which associates to the key (i,j) the (i,j) component of the jacobian. The second component of the list should be another list that contines the derivatives with respect to time"""



   cdef int i
   cdef int j
   cdef int y_n
   cdef int param_n
   cdef param_container* container
   container = <param_container *> params
   y_n=container.y_n
   param_n=container.param_n
   y_list=[]
   param_list=[]
   for i from 0<=i<y_n:
      y_list.append(y[i])
   if param_n >0:
      for i from 0 <=i<param_n:
         param_list.append(container.params[i])
   jac_list=ode_jacobian(t,y_list,param_list)
   for i from 0<=i<n:
      for j from 0<=j<n:
         dfdy[i*n+j]=jac_list[0][(i,j)]
   for i from 0 <=i<n:
      dfdt[i]=jac_list[1][i]
   return GSL_SUCCESS

cdef int c_f(double t,double* y, double* dydt,void *params):
   cdef int i
   cdef int y_n
   cdef int param_n
   cdef param_container* container
   container = <param_container *> params
   y_n= container.y_n
   param_n = container.param_n
   y_list=[]
   param_list=[]
   for i from 0<=i<y_n:
       y_list.append(y[i])
   if param_n >0:
      for i from 0 <= i<param_n:
         param_list.append(container.params[i])
   dydt_list=ode_function(t,y_list,param_list)
   for i from 0<=i<y_n:
       dydt[i]=dydt_list[i]
   return GSL_SUCCESS



def ode_solver_wrapper(g,jac,n,initial_y,params,double t_start,double t_end,double h,double error_abs,double error_rel):
   cdef int i
   global ode_function
   global ode_jacobian
   ode_function = g
   ode_jacobian = jac
   cdef double t
   cdef double *y
   cdef param_container parameters
   param_n=int(len(params))
   cdef int c_test_int
   parameters.param_n=param_n
   parameters.y_n=n
   y= <double*> PyMem_Malloc(sizeof(double)*(n))
   cdef double *c_params
   test_int=len([1,2,3,4])
   c_test_int = test_int

   if param_n > 0:
      c_params=<double *> PyMem_Malloc(sizeof(double)*(parameters.param_n))
      for i from 0<=i < parameters.param_n:
         c_params[i]=params[i]
   else:
      c_params=NULL
   for i from 0 <=i< n:
       y[i]=initial_y[i]


   parameters.params = <double *> c_params

   t = t_start
   result=[]
   v=[0]*n
   cdef gsl_odeiv_step_type * T

#Include for different types of solvers add a parameter

   T = gsl_odeiv_step_rk8pd
   cdef gsl_odeiv_step * s
   s  = gsl_odeiv_step_alloc (T, n)
   cdef gsl_odeiv_control * c
   c  = gsl_odeiv_control_y_new (error_abs, error_rel)
   cdef gsl_odeiv_evolve * e
   e  = gsl_odeiv_evolve_alloc(n)
   cdef gsl_odeiv_system sys
   sys.function = c_f

#change this
   sys.jacobian = c_jac
   sys.dimension = n
   sys.params = &parameters

   cdef int status
   import copy
   while (t < t_end):
    _sig_on
    status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t_end, &h, y)
    _sig_off
    if (status != GSL_SUCCESS):
      break
    else:
      for i from 0<=i<n:
          v[i]=<double> y[i]

   gsl_odeiv_evolve_free (e)
   gsl_odeiv_control_free (c)
   gsl_odeiv_step_free (s)
   PyMem_Free(y)
   PyMem_Free(c_params)
   return v



def f(t,y,params):
    return[y[1],-y[0]-params[0]*y[1]*(y[0]**2-1)]

def j(t,y,params):
   return[{(0,0):0,(0,1):1.0,(1,0):-2.0*params[0]*y[0]*y[1]-1.0,(1,1):-params[0]*(y[0]*y[0]-1.0)},[0,0]]

#def test():
#   return ode_solver_wrapper(f, j,1,[1],[1],0,10,1e-2,1e-10,1e-10)

def test():
    result=[(0,(1,0))]
    for i from 1 <=i<1000:
       result.append((i*100.0/1000,ode_solver_wrapper(f,j,2,[result[i-1][1][0],result[i-1][1][1]],[10],(i-1)*100.0/1000,i*(100.0/1000),1e-2,1e-10,1e-10)) )
    return result


#(g,jac,n,initial_y,params,double t_start,double t_end,double h,double error_abs,double error_rel)

#def f(t,y,params):
#    v=[]
#    v.append(y[1])
#    v.append(-y[0]-params[0]*y[1]*(y[0]*y[0]-1))
#    return v
#
#def test():
#    simple_ode_solve(f,f,2,(1,0),(10),0,1,1e-3,1e-3,1e-3)

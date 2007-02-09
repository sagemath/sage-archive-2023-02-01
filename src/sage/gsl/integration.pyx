##############################################################################
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################
#TODO: Make sure this works with compiled functions
#TODO: Expose fourier integral functionality

include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
include 'gsl.pxi'

import sage.rings.complex_double
import sage.plot.plot
import sage.gsl.interpolation

cdef class PyFunctionWrapper:
   cdef object the_function
   cdef object the_parameters

cdef class compiled_integrand:
   cdef int c_f(self,double t):  #void *params):
      return 0

cdef double c_f_compiled(double t, void *params):
   cdef double value
   cdef compiled_integrand wrapper
   wrapper = <compiled_integrand> params
   status =  wrapper.c_f(t)  #Could add paramters
   return value

cdef double c_f(double t,void *params):
   cdef double value
   cdef PyFunctionWrapper wrapper
   wrapper = <PyFunctionWrapper> params
   try:
      if len(wrapper.the_parameters)!=0:
         value=wrapper.the_function(t,wrapper.the_parameters)
      else:
         value=wrapper.the_function(t)
   except:
      return NaN

   return value



def NumericalIntegration(func,interval=[0,1],algorithm='qag',max_points=87,params=[],eps_abs=1e-6,eps_rel=1e-6,rule=6):
   """
   This uses GSL to numerically interate a function. The keyword options are

   interval: The interval of integration, specified as tuple/list with the first element the lower bound and the second element the upper bound.

   algorithm: valid choices are
       1.'qag' for an adaptive integration
       2. 'qng' for a non-adaptive gauss kronrod (samples at a maximum of 87pts)

   max_points: Sets the maximum number of sample points
   params: used to pass parameters to your function
   eps_abs, eps_rel:Absolute and relative error tolerances

   rule: This controls the Gauss-Kronrod rule used in the adaptive integration
   rule=1: 15 pt rule
   rule=2: 21 pt rule
   rule=3: 31 pt rule
   rule=4: 41 pt rule
   rule=5: 51 pt rule
   rule=6: 61 pt sule

   Higher key values are more accurate for smooth functions but lower key values deal better with discontinuities.

   NumericalIntegration returns a tuple whose first component is the answer and whose second component is an error estimate.

   For example to integrate the function x^2 from 0 to 1, we can do
   sage: f= lambda x:x**2
   sage: NumericalIntegration(f,interval=[0,1],max_points=100)

   If we want to change the error tolerances and gauss rule used
   sage:NumericalIntegration(f,interval=[0,1],max_points=200,err_abs=1e-7,arr_rel=1e-7,key=4)

   For a function with parameters
   sage: f=lambda x,a:1.0/(a[0]+x**2)
   sage: [NumericalIntegration(f,[1,2],max_points=100,params=[n]) for n in range(10)]
   Note the parameters are always a tuple everen if they have one component

   It is possible to perform on infinite intervals as well by using 'inf' or
   '-inf' in the interval argument. For example
   sage:f=lambda x: float(exp(RR(-x)))
   sage:NumericalIntegration(f,interval=[0,'inf'])
   Note the RR, this prevents underflow which can

   sage:f=lambda x:float(exp(RR(-x**2)))
   sage:NumericalIntegration(f,interval=['-inf','inf'])

   """
   import inspect
   cdef double abs_err # step size
   cdef double result
   cdef int i
   cdef int j
   cdef double a
   cdef double b
   cdef int type
   cdef PyFunctionWrapper wrapper #struct to pass information into GSL C function
   #self.params=params
   type = isinstance(func,compiled_integrand)
 #  if interval!=False:
 #     self.interval=interval

   if type == 0:
      wrapper = PyFunctionWrapper()
      if func!=None:
         wrapper.the_function = func
      else:
         raise ValueError, "No integrand defined"
      if params==[] and len(inspect.getargspec(wrapper.the_function)[0])==1:
         wrapper.the_parameters=[]
      elif params==[] and len(inspect.getargspec(wrapper.the_function)[0])>1:
         raise ValueError, "Integrand has parameters but no parameters specified"
      elif params!=[]:
         wrapper.the_parameters = params


   cdef gsl_function F
   cdef gsl_integration_workspace* W
   W=NULL
   F.function=c_f
   F.params=<void *> wrapper
   cdef size_t n
   n=max_points
   if algorithm=="qng":
      a=interval[0]
      b=interval[1]
      _sig_on
      gsl_integration_qng(&F,a,b,eps_abs,eps_rel,&result,&abs_err,&n)
      _sig_off

   elif algorithm=="qag":
      if interval[0]!="-inf" and interval[1]!="inf":
         a=interval[0]
         b=interval[1]
         W = <gsl_integration_workspace*> gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qag(&F,a,b,eps_abs,eps_rel,n,rule,W,&result,&abs_err)
         _sig_off

      if interval[0]=="-inf" and interval[1] =="inf":
#         a=interval[0]
#         b=interval[1]
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qagi(&F,eps_abs,eps_rel,n,W,&result,&abs_err)
         _sig_off
      elif interval[0]!="-inf" and interval[1]=="inf":
         a=interval[0]
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qagiu(&F,a,eps_abs,eps_rel,n,W,&result,&abs_err)
         _sig_off
      elif interval[0]=="-inf" and interval[1]!="inf":
         b=interval[1]
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qagil(&F,b,eps_abs,eps_rel,n,W,&result,&abs_err)
         _sig_off
   else:
      raise TypeError, "invalid integration algorithm"

   if W!=NULL:
      gsl_integration_workspace_free(W)
   return [result,abs_err]



## class FourierTransform:
##    def __init__(self,f):
##       self.f=f
##    def neg(self,x):
##       return self.f(-x)
##    def __call__(self,w):
##       a_1=sage.rings.complex_double.CDF(FourierIntegral(self.f,interval=[0,'inf'],omega=-w,type='cos')[0],FourierIntegral(self.f,interval=[0,'inf'],omega=-w,type='sin')[0])
##       a_2=sage.rings.complex_double.CDF(FourierIntegral(self.neg,interval=[0,'inf'],omega=w,type='cos')[0],FourierIntegral(self.neg,interval=[0,'inf'],omega=w,type='cos')[0])
##       return a_1+a_2



## def FourierIntegral(f,interval=[0,'inf'],omega=1,type='sin',eps_abs=1e-6,eps_rel=1e-6,max_points=100):
##    import inspect
##    cdef double abs_err
##    cdef double result
##    cdef double a
##    cdef PyFunctionWrapper wrapper
##    cdef w
##    cdef gsl_function F
##    cdef gsl_integration_workspace* W
##    cdef gsl_integration_workspace* CW
##    cdef gsl_integration_qawo_table* wf
##    W=NULL
##    CW=NULL
##    wf=NULL
##    w=omega
##    wrapper=PyFunctionWrapper()
##    wrapper.the_function=f
##    wrapper.the_parameters=[]
##    F.function=c_f
##    F.params=<void *> wrapper
##    cdef size_t n
##    n=max_points
##    cdef size_t sin_or_cos
##    if type=='sin':
##       sin_or_cos=GSL_INTEG_SINE
##    else:
##       sin_or_cos=GSL_INTEG_COSINE

##    if interval[0]!='-inf' and interval[1]=='inf':
##       a=interval[0]
##       _sig_on
##       W= <gsl_integration_workspace*> gsl_integration_workspace_alloc(n)
##       CW=<gsl_integration_workspace*> gsl_integration_workspace_alloc(n)
##       wf = <gsl_integration_qawo_table*> gsl_integration_qawo_table_alloc(w,1,sin_or_cos,n)
##       gsl_integration_qawf(&F,a,eps_abs,n,W,CW,wf,&result,&abs_err)
##       _sig_off
##    gsl_integration_workspace_free(W)
##    gsl_integration_workspace_free(CW)
##    gsl_integration_qawo_table_free(wf)
##    return [result,abs_err]





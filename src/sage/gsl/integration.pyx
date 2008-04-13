"""
Numerical Integration

AUTHORS:
    -- Josh Kantor (2007-02): first version
    -- William Stein (2007-02): rewrite of docs, conventions, etc.
"""

##############################################################################
#       Copyright (C) 2004,2005,2006,2007 Joshua Kantor <kantor.jm@gmail.com>
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
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
   status =  wrapper.c_f(t)  #Could add parameters
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
   except Exception, msg:
      print msg
      return 0

   return value



def numerical_integral(func, a, b=None,
                       algorithm='qag',
                       max_points=87, params=[], eps_abs=1e-6,
                       eps_rel=1e-6, rule=6):
   r"""
   Returns the numerical integral of the function on the interval
   from a to b and an error bound.

   EXAMPLES:
      To integrate the function $x^2$ from 0 to 1, we do
          sage: numerical_integral(lambda x: x^2, 0, 1, max_points=100)
          (0.33333333333333331, 3.7007434154171879e-15)

      To integrate the function $\sin(x)^3 + \sin(x)$ we do
         sage: numerical_integral(lambda x: sin(x)^3 + sin(x),  0, pi)
         (3.333333333333333, 3.7007434154171883e-14)

      We check this with a symbolic integration:
         sage: (sin(x)^3+sin(x)).integral(x,0,pi)
         10/3

   INPUT:
      -- a, b: The interval of integration, specified as two numbers
               or as a tuple/list with the first element the lower bound
               and the second element the upper bound.  Use
               +Infinity and -Infinity for plus or minus infinity.
      -- algorithm: valid choices are
                    'qag' -- for an adaptive integration
                    'qng' -- for a non-adaptive gauss kronrod (samples at a maximum of 87pts)
      -- max_points: sets the maximum number of sample points
      -- params: used to pass parameters to your function
      -- eps_abs, eps_rel: absolute and relative error tolerances
      -- rule: This controls the Gauss-Kronrod rule used in the adaptive integration
               rule=1: 15 pt rule
               rule=2: 21 pt rule
               rule=3: 31 pt rule
               rule=4: 41 pt rule
               rule=5: 51 pt rule
               rule=6: 61 pt rule
         Higher key values are more accurate for smooth functions but lower
         key values deal better with discontinuities.

   OUTPUT:
       numerical_integral returns a tuple whose first component is
       the answer and whose second component is an error estimate.

   REMARK:
       There is also a method \code{nintegral} on symbolic expressions
       that implements numerical integration using Maxima.  It is potentially
       very useful for symbolic expressions.


   MORE EXAMPLES:
   If we want to change the error tolerances and gauss rule used
       sage: f = lambda x: x^2
       sage: numerical_integral(f, 0, 1, max_points=200, eps_abs=1e-7, eps_rel=1e-7, rule=4)
       (0.33333333333333331, 3.7007434154171879e-15)

   For a Python function with parameters:
      sage: f = lambda x, a:1.0/(a[0]+x**2)
      sage: [numerical_integral(f, 1, 2, max_points=100, params=[n]) for n in range(10)]   # slightly random output (architecture and os dependent)
      [(0.49999999999998657, 5.5511151231256336e-15),
       (0.32175055439664557, 3.5721487367706477e-15),
       (0.24030098317249229, 2.6678768435816325e-15),
       (0.19253082576711697, 2.1375215571674764e-15),
       (0.16087527719832367, 1.7860743683853337e-15),
       (0.13827545676349412, 1.5351659583939151e-15),
       (0.12129975935702741, 1.3466978571966261e-15),
       (0.10806674191683065, 1.1997818507228991e-15),
       (0.09745444625548845, 1.0819617008493815e-15),
       (0.088750683050217577, 9.8533051773561173e-16)]

   Note the parameters are always a tuple even if they have one component.

   It is possible to perform on infinite intervals as well by using
   +Infinity or -Infinity in the interval argument. For example,
       sage: f = lambda x: float(exp(RR(-x)))
       sage: numerical_integral(f, 0, +Infinity)       # slightly random output
       (0.99999999999957279, 1.8429811298996553e-07)

   Note the coercion to the real field RR, which prevents underflow.

       sage: f = lambda x: float(exp(RR(-x**2)))
       sage: numerical_integral(f, -Infinity, +Infinity)           # slightly random output
       (1.7724538509060035, 3.4295192165889879e-08)

   One can integrate any real-valued callable function:
       sage: numerical_integral(lambda x: abs(zeta(x)), [1.1,1.5])           # slightly random output
       (1.8488570602160455, 2.052643677492633e-14)

   We can also numerically integrate symbolic expressions using either this
   function (which uses GSL) or the native integration (which uses Maxima):
       sage: exp(-1/x).nintegral(x, 1, 2)   # via maxima
       (0.50479221787318396, 5.6043194293440752e-15, 21, 0)
       sage: numerical_integral(exp(-1/x), 1, 2)
       (0.50479221787318407, 5.6043194293440744e-15)

   IMPLEMENTATION NOTES:
       Uses calls to the GSL -- the GNU Scientific Library -- C library.

   AUTHORS:
       -- Josh Kantor
       -- William Stein
   """

   import inspect
   cdef double abs_err # step size
   cdef double result
   cdef int i
   cdef int j
   cdef double _a, _b
   cdef PyFunctionWrapper wrapper  # struct to pass information into GSL C function

   if b is None or isinstance(a, (list, tuple)):
      b = a[1]
      a = a[0]

   if not isinstance(func, compiled_integrand):
      wrapper = PyFunctionWrapper()
      if not func is None:
         wrapper.the_function = func
      else:
         raise ValueError, "No integrand defined"
      try:
         if params==[] and len(inspect.getargspec(wrapper.the_function)[0])==1:
            wrapper.the_parameters=[]
         elif params==[] and len(inspect.getargspec(wrapper.the_function)[0])>1:
            raise ValueError, "Integrand has parameters but no parameters specified"
         elif params!=[]:
            wrapper.the_parameters = params
      except TypeError:
         wrapper.the_function = eval("lambda x: func(x)", {'func':func})
         wrapper.the_parameters = []


   cdef gsl_function F
   cdef gsl_integration_workspace* W
   W=NULL
   F.function=c_f
   F.params=<void *> wrapper
   cdef size_t n
   n=max_points

   gsl_set_error_handler_off()

   if algorithm=="qng":
      _a=a
      _b=b
      _sig_on
      gsl_integration_qng(&F,_a,_b,eps_abs,eps_rel,&result,&abs_err,&n)
      _sig_off

   elif algorithm=="qag":
      from sage.rings.infinity import Infinity
      if a is -Infinity and b is +Infinity:
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qagi(&F,eps_abs,eps_rel,n,W,&result,&abs_err)
         _sig_off

      elif a is -Infinity:
         _b=b
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qagil(&F,_b,eps_abs,eps_rel,n,W,&result,&abs_err)
         _sig_off

      elif b is +Infinity:
         _a=a
         W=<gsl_integration_workspace*>gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qagiu(&F,_a,eps_abs,eps_rel,n,W,&result,&abs_err)
         _sig_off

      else:
         _a=a
         _b=b
         W = <gsl_integration_workspace*> gsl_integration_workspace_alloc(n)
         _sig_on
         gsl_integration_qag(&F,_a,_b,eps_abs,eps_rel,n,rule,W,&result,&abs_err)
         _sig_off

   else:
      raise TypeError, "invalid integration algorithm"

   if W != NULL:
      gsl_integration_workspace_free(W)

   return result, abs_err


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





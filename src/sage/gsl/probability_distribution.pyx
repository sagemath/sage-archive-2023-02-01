"""nodoctest
Probability Distributions

AUTHORS:
    -- Josh Kantor (2007-02): first version
    -- William Stein (2007-02): rewrite of docs, conventions, etc.
"""


##############################################################################
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

import sage.plot.plot
#include '../ext/cdefs.pxi'
#include '../ext/interrupt.pxi'
#include 'gsl.pxi'
#cimport sage.rings.real_double
import sage.rings.real_double
import random
import integration
from sage.modules.free_module_element import vector


#TODO: Add more distributions available in gsl
#available but not currently wrapped are exponential,laplace,cauchy,landau,gamma,
#gamma,beta logistic.

cdef enum:
   uniform
   gaussian
   rayleigh
   lognormal
   pareto
   t
   chisquared
   exppow
   weibull
   beta

cdef class ProbabilityDistribution:
   def __init__(self):
      pass

   def get_random_element(self):
      raise NotImplementedError,"implement in derived class"

cdef class SphericalDistribution(ProbabilityDistribution):
   r"""
   This class is capable of producing random points uniformly distributed on the surface of an
   n-1 sphere in n dimensional euclidean space. The dimension, n is selected via the
   keyword dimension. The random number generator which drives it can be selected using the keyword
   rng. Valid choices are 'default' which uses the Mersenne-Twister, 'luxury' which uses RANDLXS,
   and 'taus' which uses the tausworth generator. The default dimension is 3.
   sage: T=SphericalDistribution()
   sage: T.get_random_element()
   sage: T=SphericalDistribution(dimension=4,rng='luxury')
   sage: T.get_random_element()
   """
   def __init__(self,dimension=3,rng='default',seed=None):
      gsl_rng_env_setup()
      self.set_random_number_generator(rng)
      self.r=gsl_rng_alloc(self.T)
      if seed==None:
         self.seed=random.randint(1,2^32)
      self.set_seed(self.seed)
      self.dimension=dimension
      self.vec=<double *>malloc(self.dimension*(sizeof(double)))

   def set_seed(self,seed):
      gsl_rng_set(self.r,seed)
      self.seed=seed

   def set_random_number_generator(self,rng='default'):
      if rng=='default':
         self.T=gsl_rng_default
      elif rng=='luxury':
         self.T=gsl_rng_ranlxd2
      elif rng=='taus':
         self.T=gsl_rng_taus2
      else:
         raise TypeError,"Not a valid random number generator"

   def __dealloc__(self):
      if self.r!=NULL:
         gsl_rng_free(self.r)
      free(self.vec)

   def get_random_element(self):
      cdef int i
      v=[0]*self.dimension
      gsl_ran_dir_nd(self.r,self.dimension,self.vec)
      for i from 0<=i<self.dimension:
         v[i]=self.vec[i]
      return vector(sage.rings.real_double.RDF,v) #This could be made more efficient by directly constructing the vector, TODO.

   def reset_distribution(self):
      """
      This method resets the distribution.

      EXAMPLES:
         sage: T = SphericalDistribution()
         sage: v = [T.get_random_element() for _ in range(10)]
         sage: T.reset_distribution()
         sage: w = [T.get_random_element() for _ in range(10)]
      """
      if self.r!=NULL:
         gsl_rng_free(self.r)
      self.r = gsl_rng_alloc(self.T)
      self.set_seed(self.seed)
#      gsl_rng_env_setup()



cdef class RealDistribution(ProbabilityDistribution):
   r"""
   The RealDistribution provides a number of routines for sampling
   from and analyzing and visualizing probability distributions.  For
   precise definitions of the distributions and their parameters see
   the gsl reference manuals chapter on random number generators and
   probability distributions for precise definitions.  The probability
   distributions available currently are uniform.

   EXAMPLES:
   To create a uniform distribution on the interval [a,b]:
      sage: a=0
      sage: b=2
      sage: T=RealDistribution('uniform',[a,b])
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   gaussian:
   The gaussian distribuiton takes 1 parameters sigma, sigma =1 gives the standard
   gaussian distribution

      sage: sigma=1
      sage: T=RealDistribution('gaussian',sigma)
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   rayleigh>
   The rayleigh distribution has 1 parameter sigma.

      sage: sigma=3
      sage: T=RealDistribution('rayleigh',sigma)
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   lognormal>
   The lognormal distribution has two parameters sigma and zeta:
      sage: zeta=0
      sage: sigma=1
      sage: T=RealDistribution('lognormal',[zeta,sigma])
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   pareto>
   The pareto distribution has two parameters a, and b:
      sage: a=1
      sage: b=1
      sage: T=RealDistribution('pareto',[a,b])
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   student's t>
   The student's t distribution has one parameter nu.
      sage: nu=1
      sage: T=RealDistribution('t',nu)
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   chi squared>
      The chi squared distribution has one parameter nu
      sage: nu=1
      sage: T=RealDistribution('chisquared',nu)
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   exponential power>
      The exponential power distribution has two parameters a and b
      sage: a=1
      sage: b=2.5
      sage: T=RealDistribution('exppow',[a,b])
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   beta>
      sage: a=2
      sage: b=2
      sage: T=RealDistribution('beta',[a,b])
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)

   weibull>
      sage: a=1
      sage: b=1
      sage: T=RealDistribution('weibull',[a,b])
      sage: T.get_random_element()
      sage: T.distribution_function(0)
      sage: T.cum_distribution_function(1)
      sage: T.cum_distribution_function_inv(.5)

   It is possible to select which random number generator drives the
   sampling as well as the seed.  The default is the Mersenne
   twister. Also available are the RANDLXS algorithm and the
   Tausworthe generator (see the GSL reference manual for more
   details). These are all supposed to be simulation quality
   generators.   For RANDLXS us rng='luxury' for tausworth rng='taus'

       sage: T = RealDistribution('gaussian',1,rng='luxury',seed=10)

   To change the seed at a later time use \code{set_seed}:
       sage: T.set_seed(100)
   """
   def __init__(self,type='uniform',parameters=[],rng='default',seed=None):
      gsl_rng_env_setup()
      self.parameters = NULL
      self.set_random_number_generator(rng)
      self.r=gsl_rng_alloc(self.T)
      if seed==None:
         self.seed=random.randint(1,2^32)
      self.set_seed(self.seed)
      self.name =" "
      self.set_distribution(type,parameters)

   def set_seed(self,seed):
      gsl_rng_set(self.r,seed)
      self.seed=seed

   def set_random_number_generator(self,rng='default'):
      if rng=='default':
         self.T=gsl_rng_default
      elif rng=='luxury':
         self.T=gsl_rng_ranlxd2
      elif rng=='taus':
         self.T=gsl_rng_taus2
      else:
         raise TypeError,"Not a valid random number generator"

   def __dealloc__(self):
      if self.r!=NULL:
         gsl_rng_free(self.r)
      if self.parameters!=NULL:
         free(self.parameters)

   def __str__(self):
      return self.name


   def get_random_element(self):
      """
      This method generates a random sample from the current probability distribution
      sage: T=RealDistribution('gaussian',1)
      sage: T.get_random_element()
      """
      cdef double result
      if self.distribution_type==uniform:
         result = gsl_ran_flat(self.r,self.parameters[0],self.parameters[1])
#         result = gsl_rng_uniform(self.r)
      elif self.distribution_type==gaussian:
         result=gsl_ran_gaussian(self.r,self.parameters[0])
      elif self.distribution_type==rayleigh:
         result=gsl_ran_rayleigh(self.r,self.parameters[0])
      elif self.distribution_type==lognormal:
         result=gsl_ran_lognormal(self.r,self.parameters[0],self.parameters[1])
      elif self.distribution_type==pareto:
         result=gsl_ran_pareto(self.r,self.parameters[0],self.parameters[1])
      elif self.distribution_type==t:
         result=gsl_ran_tdist(self.r,self.parameters[0])
      elif self.distribution_type==chisquared:
         result=gsl_ran_chisq(self.r,self.parameters[0])
      elif self.distribution_type==exppow:
         result=gsl_ran_exppow(self.r,self.parameters[0],self.parameters[1])
      elif self.distribution_type==weibull:
         result=gsl_ran_weibull(self.r,self.parameters[0],self.parameters[1])
      elif self.distribution_type==beta:
         result=gsl_ran_beta(self.r,self.parameters[0],self.parameters[1])


      return sage.rings.real_double.RDF(result)
   def set_distribution(self,name='uniform',parameters = []):
      """
      This method can be called to change the current probability distribution.

      EXAMPLES:
         sage: T=RealDistribution('gaussian', 1)
         sage: T.set_distribution('gaussian', 1)
         sage: T.set_distribution('pareto', [0,1])
      """
      if self.parameters!=NULL:
         free(self.parameters)

      if name =='uniform':
        self.distribution_type = uniform
        for x in parameters:
           try:
              float(x)
           except:
              raise TypeError, "Uniform distribution requires parameters coercible to float"
        self.parameters=<double*>malloc(sizeof(double)*2)
        self.parameters[0] = parameters[0]
        self.parameters[1] = parameters[1]
      elif name == 'gaussian':
         try:
            float(parameters)
         except:
            raise TypeError,"gaussian distribution requires parameter sigma coercible to float"
         self.parameters=<double*>malloc(sizeof(double))
         self.parameters[0]=float(parameters)
         self.distribution_type = gaussian
      elif name == 'pareto':
         if len(parameters)!=2:
            raise TypeError, "pareto distribution has two prameters"
         try:
            map(float,parameters)
         except:
            raise TypeError, "parameters must be coercible to float"
         self.parameters=<double*>malloc(sizeof(double)*2)
         self.parameters[0]=float(parameters[0])
         self.parameters[1]=float(parameters[1])
         self.distribution_type=pareto
      elif name =='rayleigh':
         self.distribution_type=rayleigh
         try:
            float(parameters)
         except:
            raise TypeError,"rayleigh distribution requires parameter sigma coercible to float"
         self.parameters=<double*>malloc(sizeof(double))
         self.parameters[0]=float(parameters)
         self.distribution_type = rayleigh
      elif name=='lognormal':
         if len(parameters)!=2:
            raise TypeError, "Lognormal distribution requires two parameters"
         for x in parameters:
            try:
               float(x)
            except:
               raise TypeError,"Lognormal distributions requires real parameters"

         self.parameters=<double*>malloc(sizeof(double)*2)
         self.parameters[0]=float(parameters[0])
         self.parameters[1]=float(parameters[1])
         self.distribution_type=lognormal
      elif name=='t':
         try:
            float(parameters)
         except:
            raise TypeError,"parameter to t distribution must be coercible to float"
         self.parameters=<double*>malloc(sizeof(double))
         self.parameters[0]=float(parameters)
         self.distribution_type=t

      elif name=='chisquared':
         try:
            float(parameters)
         except:
            raise TypeError,"parameters to t distribution must be coercible to float"
         self.parameters=<double *>malloc(sizeof(double))
         self.parameters[0]=float(parameters)
         self.distribution_type=chisquared
      elif name=='exppow':
         if len(parameters)!=2:
            raise TypeError, "exponential power distribution requires two parameters"
         for x in parameters:
            try:
               float(x)
            except:
               raise TypeError,"exponential power distributions requires real parameters"

         self.parameters=<double*>malloc(sizeof(double)*2)
         self.parameters[0]=float(parameters[0])
         self.parameters[1]=float(parameters[1])
         self.distribution_type=exppow
      elif name=='weibull':
         if len(parameters)!=2:
            raise TypeError, "weibull distribution requires two real parameters"
         try:
            map(float,parameters)
         except:
            raise TypeError,"weibull distribution requires real parameters"

         self.parameters=<double *>malloc(sizeof(double)*2)
         self.parameters[0]=float(parameters[0])
         self.parameters[1]=float(parameters[1])
         self.distribution_type=weibull

      elif name=='beta':
         if len(parameters)!=2:
            raise TypeError, "beta distribution requires two real parameters"
         try:
            map(float,parameters)
         except:
            raise TypeError,"beta distribution requires real parameters"

         self.parameters=<double *>malloc(sizeof(double)*2)
         self.parameters[0]=float(parameters[0])
         self.parameters[1]=float(parameters[1])
         self.distribution_type=beta


      else:
         raise TypeError,"Not a valid probability distribution"

      self.name=name

  # def _get_random_element_c():


   def reset_distribution(self):
      """
      This method resets the distribution.

      EXAMPLES:
         sage: T = RealDistribution('gaussian',1,seed=10)
         sage: v = [T.get_random_element() for _ in range(10)]
         sage: T.reset_distribution()
         sage: w = [T.get_random_element() for _ in range(10)]
      """
      if self.r!=NULL:
         gsl_rng_free(self.r)
      self.r = gsl_rng_alloc(self.T)
      self.set_seed(self.seed)
#      gsl_rng_env_setup()

   def distribution_function(self,x):
      """
      This methods evaluates the distribution function of the
      probability distribution at $x$.
      """
      if self.distribution_type ==uniform:
         return sage.rings.real_double.RDF(gsl_ran_flat_pdf(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==gaussian:
         return sage.rings.real_double.RDF(gsl_ran_gaussian_pdf(x,self.parameters[0]))
      elif self.distribution_type==rayleigh:
         return sage.rings.real_double.RDF(gsl_ran_rayleigh_pdf(x,self.parameters[0]))
      elif self.distribution_type==lognormal:
         return sage.rings.real_double.RDF(gsl_ran_lognormal_pdf(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==pareto:
         return sage.rings.real_double.RDF(gsl_ran_pareto_pdf(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==t:
         return sage.rings.real_double.RDF(gsl_ran_tdist_pdf(x,self.parameters[0]))
      elif self.distribution_type==chisquared:
         return sage.rings.real_double.RDF(gsl_ran_chisq_pdf(x,self.parameters[0]))
      elif self.distribution_type==exppow:
         return sage.rings.real_double.RDF(gsl_ran_exppow_pdf(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==weibull:
         return sage.rings.real_double.RDF(gsl_ran_weibull_pdf(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==beta:
         return sage.rings.real_double.RDF(gsl_ran_beta_pdf(x,self.parameters[0],self.parameters[1]))


   def cum_distribution_function(self,x):
      """
      This method evaluates the cumulative distribution function of
      the probability distribution at $x$.
      """
      if self.distribution_type ==uniform:
         return sage.rings.real_double.RDF(gsl_cdf_flat_P(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==gaussian:
         return sage.rings.real_double.RDF(gsl_cdf_gaussian_P(x,self.parameters[0]))
      elif self.distribution_type==rayleigh:
         return sage.rings.real_double.RDF(gsl_cdf_rayleigh_P(x,self.parameters[0]))
      elif self.distribution_type==lognormal:
         return sage.rings.real_double.RDF(gsl_cdf_lognormal_P(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==pareto:
         return sage.rings.real_double.RDF(gsl_cdf_pareto_P(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==t:
         return sage.rings.real_double.RDF(gsl_cdf_tdist_P(x,self.parameters[0]))
      elif self.distribution_type==chisquared:
         return sage.rings.real_double.RDF(gsl_cdf_chisq_P(x,self.parameters[0]))
      elif self.distribution_type==exppow:
         return sage.rings.real_double.RDF(gsl_cdf_exppow_P(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==weibull:
         return sage.rings.real_double.RDF(gsl_cdf_weibull_P(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==beta:
         return sage.rings.real_double.RDF(gsl_cdf_beta_P(x,self.parameters[0],self.parameters[1]))


   def cum_distribution_function_inv(self,x):
      """
      This method evaluates the inverse of the cumulative distribution
      distribution function of the probability distribution at $x$.
      """
      if self.distribution_type ==uniform:
         return sage.rings.real_double.RDF(gsl_cdf_flat_Pinv(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==gaussian:
         return sage.rings.real_double.RDF(gsl_cdf_gaussian_Pinv(x,self.parameters[0]))
      elif self.distribution_type==rayleigh:
         return sage.rings.real_double.RDF(gsl_cdf_rayleigh_Pinv(x,self.parameters[0]))
      elif self.distribution_type==lognormal:
         return sage.rings.real_double.RDF(gsl_cdf_lognormal_Pinv(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==pareto:
         return sage.rings.real_double.RDF(gsl_cdf_pareto_Pinv(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==t:
         return sage.rings.real_double.RDF(gsl_cdf_tdist_Pinv(x,self.parameters[0]))
      elif self.distribution_type==chisquared:
         return sage.rings.real_double.RDF(gsl_cdf_chisq_Pinv(x,self.parameters[0]))
      elif self.distribution_type==exppow:
         raise NotImplementedError,"gsl does not provide inverse for exponential power"
#         return sage.rings.real_double.RDF(gsl_cdf_exppow_Pinv(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==weibull:
         return sage.rings.real_double.RDF(gsl_cdf_weibull_Pinv(x,self.parameters[0],self.parameters[1]))
      elif self.distribution_type==beta:
         return sage.rings.real_double.RDF(gsl_cdf_beta_Pinv(x,self.parameters[0],self.parameters[1]))

   def plot(self, *args, **kwds):
      return sage.plot.plot.plot(self.distribution_function, *args, **kwds)

   def generate_histogram_data(self,num_samples=1000,bins=50):
      import pylab
      l=[float(self.get_random_element()) for _ in range(num_samples)]
      S=pylab.hist(l,bins,normed=True,hold=False)
      return [list(S[0]),list(S[1])]

   def generate_histogram_plot(self,name,num_samples=1000,bins=50):
      import pylab
      l = [float(self.get_random_element()) for _ in range(num_samples)]
      pylab.hist(l,bins,normed=True,hold=False)
      pylab.savefig(name)






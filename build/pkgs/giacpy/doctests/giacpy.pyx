#*****************************************************************************
#    AUTHOR:  Han Frederic <frederic.han@imj-prg.fr>
#  2012
#  Distributed under the terms of the GNU General Public License (GPL)
#  version 2 or higher.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
:Name: giacpy
:Summary: A Cython frontend to the c++ library giac. (Computer Algebra System)
:License:  GPL v2 or above
:Home-page: http://www.math.jussieu.fr/~han/xcas/giacpy/
:Author: Frederic Han
:Author-email: frederic.han@imj-prg.fr


This is the sage version of giacpy.

- Giacpy is an interface to the c++ giac library. This interface is built with cython, and the resulting speed is similar to giac's one.

- Giac is a general purpose Computer algebra system by Bernard Parisse released under GPLv3.

      * http://www-fourier.ujf-grenoble.fr/~parisse/giac.html
      * It is build on C and C++ libraries:
        NTL (arithmetic), GSL (numerics), GMP (big integers), MPFR (bigfloats)
      * It  provides fast  algorithms  for multivariate polynomial operations (product, GCD, factorisation) and
      * symbolic  computations: solver, simplifications, limits/series, integration, summation...
      * Linear Algebra with numerical or symbolic coefficients.


AUTHORS:

- Frederic Han (2013-09-23): initial version


EXAMPLES:

- The class Pygen is the main tool to interact from python/sage with the c++ library giac via cython.

  The initialisation of a Pygen just create an object in giac, but the mathematical computation  is not done. This class is mainly for cython users.


  Here A is a Pygen element, and it is ready for any giac function.

  ::

        sage: from giacpy import *    # random
        //...
        sage: A=Pygen('2+2');A
        2+2
        sage: A.eval()
        4

  In general, you may prefer to directly create a Pygen and execute the evaluation in giac. This is exactly  the meaning of the libgiac() function.

  ::

        sage: a=libgiac('2+2');a;isinstance(a,Pygen)
        4
        True

- Most common usage of this package in sage will be with the libgiac() function. This function is just the composition of the Pygen initialisation and the evaluation of this object in giac.

  ::

        sage: from giacpy import libgiac,giacsettings
        sage: x,y,z=libgiac('x,y,z');  # add some giac objects
        sage: f=(x+3*y)/(x+z+1)^2 -(x+z+1)^2/(x+3*y)
        sage: f.factor()
        (3*y-x^2-2*x*z-x-z^2-2*z-1)*(3*y+x^2+2*x*z+3*x+z^2+2*z+1)/((x+z+1)^2*(3*y+x))
        sage: f.normal()
        (-x^4-4*x^3*z-4*x^3-6*x^2*z^2-12*x^2*z-5*x^2+6*x*y-4*x*z^3-12*x*z^2-12*x*z-4*x+9*y^2-z^4-4*z^3-6*z^2-4*z-1)/(x^3+3*x^2*y+2*x^2*z+2*x^2+6*x*y*z+6*x*y+x*z^2+2*x*z+x+3*y*z^2+6*y*z+3*y)

- To obtain more hints on giacpy consider the help of the :func:`libgiac<_giac>` function.

  ::

        sage: libgiac?             # doctest: +SKIP

- Some settings of giac are available via the ``giacsettings`` element. (Ex: maximal number of threads in computations, allowing probabilistic algorithms or not...

  ::

        sage: R=PolynomialRing(QQ,8,'x')
        sage: I=sage.rings.ideal.Katsura(R,8)
        sage: giacsettings.proba_epsilon=1e-15;
        sage: Igiac=libgiac(I.gens());
        sage: time Bgiac=Igiac.gbasis([R.gens()],'revlex')  # doctest: +SKIP
        Running a probabilistic check for the reconstructed Groebner basis. If successfull, error probability is less than 1e-15 and is estimated to be less than 10^-109. Use proba_epsilon:=0 to certify (this takes more time).
        Time: CPU 0.46 s, Wall: 0.50 s
        sage: giacsettings.proba_epsilon=0;
        sage: Igiac=libgiac(I.gens());
        sage: time Bgiac=Igiac.gbasis([R.gens()],'revlex') # doctest: +SKIP
        Time: CPU 2.74 s, Wall: 2.75 s

  ::

        sage: from giacpy import *
        sage: x=libgiac('x');f=1/(2+sin(5*x))
        sage: f.int()
        2/5/sqrt(3)*(atan((2*tan(5*x/2)+1)/sqrt(3))+pi*floor(5*x/2/pi+1/2))
        sage: f.series(x,0,3)
        1/2-5/4*x+25/8*x^2-125/48*x^3+x^4*order_size(x)
        sage: libgiac(sqrt(5)+pi).approx(100)
        5.377660631089582934871817052010779119637787758986631545245841837718337331924013898042449233720899343


SEEALSO:

    ``libgiac``, ``giacsettings``, ``Pygen``,``loadgiacgen``



GETTING HELP:

- To obtain some help on a giac keyword use the help() method. In sage the htmlhelp() method for Pygen element is disabled. Just use the ? or .help() method.

  ::

        sage: libgiac.gcd?             # doctest: +SKIP
        "Returns the greatest common divisor of 2 polynomials of several variables or of 2 integers or of 2 rationals.
        (Intg or Poly),(Intg or Poly)
        gcd(45,75);gcd(15/7,50/9);gcd(x^2-2*x+1,x^3-1);gcd(t^2-2*t+1,t^2+t-2);gcd((x^2-1)*(y^2-1)*z^2,x^3*y^3*z+(-(y^3))*z+x^3*z-z)
        lcm,euler,modgcd,ezgcd,psrgcd,heugcd,Gcd"

- You can find full html documentation about the **giac** functions  at:

      * http://www-fourier.ujf-grenoble.fr/~parisse/giac/doc/en/cascmd_en/

      * http://www-fourier.ujf-grenoble.fr/~parisse/giac/doc/fr/cascmd_fr/

      * http://www-fourier.ujf-grenoble.fr/~parisse/giac/doc/el/cascmd_el/

      * or in :doc:`$SAGE_LOCAL/share/giac/doc/en/cascmd_en/index.html`



REMARK:

- Graphics 2D Output via qcas (the qt frontend to giac) is removed in the sage version of giacpy.

"""
########################################################


cimport giacpy

# sage includes
from sage.libs.gmp.mpz cimport mpz_t, mpz_init_set
from sage.rings.integer cimport Integer
from sage.rings.all import ZZ, QQ, IntegerModRing
from sage.symbolic.ring import SR
from sage.rings.rational cimport Rational
from sage.matrix.matrix import is_Matrix

from sage.plot.line import line
from sage.plot.scatter_plot import scatter_plot


from sys import maxsize as Pymaxint
import os
import math


#### Python3 compatibility  #############################
from sys import version_info as Pyversioninfo
if Pyversioninfo[0]>2 :
   PythonVersion3 = True
else:
   PythonVersion3 = False

def decstring23(s):
  if PythonVersion3 :
      return s.decode()
  else:
      return s


def encstring23(s):
  if PythonVersion3 :
      return bytes(s,'UTF-8')
  else:
      return s


if PythonVersion3:
    listrange=list,range
else:
    listrange=list

####End of Python3 compatibility########################

#
include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'

########################################################


########################################################
# A global context pointer. One by giac session.
########################################################
cdef context * context_ptr=new context()
#
# Some global variables for optimisation
GIACNULL=Pygen('NULL')
#



# Create a giac setting instance
giacsettings=GiacSetting()
Pygen('I:=sqrt(-1)').eval()


#######################################################
# The wrapper to eval with giac
#######################################################
# in sage we don't export the giac function. We replace it by libgiac
def _giac(s):
   """
   This function evaluate a python/sage object with the giac library. It creates in python/sage a Pygen element and evaluate it with giac:

   - **First Example**:

     ::

       sage: from giacpy import libgiac
       sage: x,y=libgiac('x,y')
       sage: (x+2*y).cos().texpand()
       cos(x)*(2*cos(y)^2-1)-sin(x)*2*cos(y)*sin(y)



   - **Coercion, Pygen and internal giac variables**:

     The most usefull objects will be the Python object of type Pygen.

     ::

       sage: from giacpy import *
       sage: x,y,z=libgiac('x,y,z')
       sage: f=sum([x[i] for i in range(5)])^15/(y+z);f.coeff(x[0],12)
       (455*(x[1])^3+1365*(x[1])^2*x[2]+1365*(x[1])^2*x[3]+1365*(x[1])^2*x[4]+1365*x[1]*(x[2])^2+2730*x[1]*x[2]*x[3]+2730*x[1]*x[2]*x[4]+1365*x[1]*(x[3])^2+2730*x[1]*x[3]*x[4]+1365*x[1]*(x[4])^2+455*(x[2])^3+1365*(x[2])^2*x[3]+1365*(x[2])^2*x[4]+1365*x[2]*(x[3])^2+2730*x[2]*x[3]*x[4]+1365*x[2]*(x[4])^2+455*(x[3])^3+1365*(x[3])^2*x[4]+1365*x[3]*(x[4])^2+455*(x[4])^3)/(y+z)

     Warning: The complex number sqrt(-1) is exported in python as I. (But it may appears as i)

     ::

       sage: libgiac((1+I*sqrt(3))^3).normal(); libgiac(1+I)
       -8
       1+i

     Python integers and reals can be directly converted to giac.

     ::

       sage: from giacpy import *
       sage: a=libgiac(2^1024);a.nextprime();(libgiac(1.234567)).erf().approx(10)
       179769313486231590772930519078902473361797697894230657273430081157732675805500963132708477322407536021120113879871393357658789768814416622492847430639474124377767893424865485276302219601246094119453082952085005768838150682342462881473913110540827237163350510684586298239947245938479716304835356329624224137859
       0.9191788641

     The Python object y defined above is of type Pygen. It is not an internal giac variable. (Most of the time you won't need to use internal giac variables).

     ::

       sage: from giacpy import *
       sage: libgiac('y:=1');y
       1
       y
       sage: libgiac.purge('y')
       1
       sage: libgiac('y')
       y

     There are some natural coercion to Pygen elements:

     ::

       sage: from giacpy import *
       sage: libgiac(pi)>3.14 ; libgiac(pi) >3.15 ; libgiac(3)==3
       True
       False
       True

   - **Lists of Pygen and Giac lists**:

     Here l1 is a giac list and l2 is a python list of Pygen type objects.

     ::

       sage: from giacpy import *
       sage: l1=libgiac(range(10)); l2=[1/(i^2+1) for i in l1]
       sage: sum(l2)
       33054527/16762850

     So l1+l1 is done in giac and means a vector addition. But l2+l2 is done in Python so it is the list concatenation.

     ::

       sage: from giacpy import *
       sage: l1+l1
       [0,2,4,6,8,10,12,14,16,18]
       sage: l2+l2
       [1, 1/2, 1/5, 1/10, 1/17, 1/26, 1/37, 1/50, 1/65, 1/82, 1, 1/2, 1/5, 1/10, 1/17, 1/26, 1/37, 1/50, 1/65, 1/82]

     Here V is not a Pygen element. We need to push it to giac to use a giac method like dim, or we need to use an imported function.

     ::

       sage: from giacpy import *
       sage: V=[ [x[i]^j for i in range(8)] for j in range(8)]
       sage: libgiac(V).dim()
       [8,8]
       sage: libgiac.det_minor(V).factor()
       (x[6]-(x[7]))*(x[5]-(x[7]))*(x[5]-(x[6]))*(x[4]-(x[7]))*(x[4]-(x[6]))*(x[4]-(x[5]))*(x[3]-(x[7]))*(x[3]-(x[6]))*(x[3]-(x[5]))*(x[3]-(x[4]))*(x[2]-(x[7]))*(x[2]-(x[6]))*(x[2]-(x[5]))*(x[2]-(x[4]))*(x[2]-(x[3]))*(x[1]-(x[7]))*(x[1]-(x[6]))*(x[1]-(x[5]))*(x[1]-(x[4]))*(x[1]-(x[3]))*(x[1]-(x[2]))*(x[0]-(x[7]))*(x[0]-(x[6]))*(x[0]-(x[5]))*(x[0]-(x[4]))*(x[0]-(x[3]))*(x[0]-(x[2]))*(x[0]-(x[1]))


   - **Modular objects with %**

     ::

       sage: from giacpy import *
       sage: V=libgiac.ranm(5,6) % 2;
       sage: V.ker().rowdim()+V.rank()
       6
       sage: a=libgiac(7)%3;a;a%0;7%3
       1 % 3
       1
       1

     Do not confuse with the  python integers:

     ::

       sage: type(7%3)==type(a);type(a)==type(7%3)
       False
       False

   - **Syntaxes with reserved or unknown Python/sage symbols**:

     In general equations needs symbols such as ``=`` ``<`` ``>`` that have another meaning in Python or Sage.
     So those objects must be quoted.

     ::

       sage: from giacpy import *
       sage: x=libgiac('x')
       sage: (1+2*sin(3*x)).solve(x).simplify()
       Warning, argument is not an equation, solving 1+2*sin(3*x)=0
       list[-pi/18,7*pi/18]

       sage: libgiac.solve('sin(3*x)>2*sin(x)',x)        # doctest: +SKIP
       ...
       RuntimeError: Unable to find numeric values solving equation. For trigonometric equations this may be solved using assumptions, e.g. assume(x>-pi && x<pi) Error: Bad Argument Value


     You can also add some hypothesis to a giac symbol:

     ::

       sage: libgiac.assume('x>-pi && x<pi')
       x
       sage: libgiac.solve('sin(3*x)>2*sin(x)',x)
       list[((x>(-5*pi/6)) and (x<(-pi/6))),((x>0) and (x<(pi/6))),((x>(5*pi/6)) and (x<pi))]

     To remove those hypothesis use the giac function: purge

     ::

       sage: libgiac.purge('x')
       assume[[],[line[-pi,pi]],[-pi,pi]]
       sage: libgiac.solve('x>0')
       list[x>0]

     Same problems with the ..

     ::

       sage: from giacpy import *
       sage: x=libgiac('x')
       sage: f=1/(5+cos(4*x));f.int(x)
       1/2/(2*sqrt(6))*(atan(2*tan(4*x/2)/sqrt(6))+pi*floor(4*x/2/pi+1/2))
       sage: libgiac.fMax(f,'x=-0..pi').simplify()
       pi/4,3*pi/4
       sage: libgiac.fMax.help()  # doctest: +SKIP
       "Returns the abscissa of the maximum of the expression.
       Expr,[Var]
       fMax(-x^2+2*x+1,x)
       fMin"
       sage: libgiac.sum(1/(1+x^2),'x=0..infinity').simplify()
       (pi*exp(pi)^2+pi+exp(pi)^2-1)/(2*exp(pi)^2-2)

   - **From giac to sage**:

     One can convert a Pygen element to sage with the sage method.
     Get more details with:

     ::

       sage: Pygen.sage?          # doctest: +SKIP

     ::

       sage: from giacpy import *
       sage: L=libgiac('[1,sqrt(5),[1.3,x]]')
       sage: L.sage()       # All entries are converted recursively
       [1, sqrt(5), [1.3, x]]

     To obtain matrices and vectors, use the :meth:`matrix<Pygen._matrix_>` and :meth:`vector<Pygen._vector_>` commands. Get more details with:

     ::

       sage: Pygen._matrix_?          # doctest: +SKIP
       sage: Pygen._vector_?          # doctest: +SKIP

     ::

       sage: n=var('n');A=matrix([[1,2],[-1,1]])
       sage: B=libgiac(A).matpow(n)    # We compute the symbolic power on A via libgiac
       sage: C=matrix(B,SR); C         # We convert B to sage
       [                     1/2*(I*sqrt(2) + 1)^n + 1/2*(-I*sqrt(2) + 1)^n -1/2*I*sqrt(2)*(I*sqrt(2) + 1)^n + 1/2*I*sqrt(2)*(-I*sqrt(2) + 1)^n]
       [ 1/4*I*sqrt(2)*(I*sqrt(2) + 1)^n - 1/4*I*sqrt(2)*(-I*sqrt(2) + 1)^n                      1/2*(I*sqrt(2) + 1)^n + 1/2*(-I*sqrt(2) + 1)^n]
       sage: (C.subs(n=3)-A^3).expand()
       [0 0]
       [0 0]


   **MEMENTO of usual GIAC functions**:

   - *Expand with simplification*

         * ``ratnormal``, ``normal``, ``simplify``   (from the fastest to the most sophisticated)

         *  NB: ``expand`` function doesn't regroup nor cancel terms, so it could be slow. (pedagogical purpose only?)

   - *Factor/Regroup*

         * ``factor``, ``factors``, ``regroup``, ``cfactor``, ``ifactor``

   - *Misc*

         * ``unapply``, ``op``, ``subst``

   - *Polynomials/Fractions*

         * ``coeff``,  ``gbasis``, ``greduce``, ``lcoeff``, ``pcoeff``, ``canonical_form``,

         * ``proot``,  ``poly2symb``,  ``symb2poly``, ``posubLMQ``, ``poslbdLMQ``, ``VAS``, ``tcoeff``,  ``valuation``

         * ``gcd``, ``egcd``, ``lcm``, ``quo``, ``rem``, ``quorem``, ``abcuv``, ``chinrem``,

         * ``peval``, ``horner``, ``lagrange``, ``ptayl``, ``spline``,  ``sturm``,  ``sturmab``

         * ``partfrac``, ``cpartfrac``

   - *Memory/Variables*

         * ``assume``, ``about``, ``purge``, ``ans``

   - *Calculus/Exact*

         * ``linsolve``,  ``solve``,  ``csolve``,  ``desolve``,  ``seqsolve``, ``reverse_rsolve``, ``matpow``

         * ``limit``, ``series``, ``sum``, ``diff``, ``fMax``, ``fMin``,

         * ``integrate``, ``subst``, ``ibpdv``, ``ibpu``, ``preval``

   - *Calculus/Exp, Log, powers*

         * ``exp2pow``, ``hyp2exp``, ``expexpand``, ``lin``, ``lncollect``, ``lnexpand``, ``powexpand``, ``pow2exp``

   - *Trigo*

         * ``trigexpand``, ``tsimplify``, ``tlin``, ``tcollect``,

         * ``halftan``, ``cos2sintan``, ``sin2costan``, ``tan2sincos``, ``tan2cossin2``, ``tan2sincos2``, ``trigcos``, ``trigsin``, ``trigtan``, ``shift_phase``

         * ``exp2trig``, ``trig2exp``

         * ``atrig2ln``, ``acos2asin``, ``acos2atan``, ``asin2acos``, ``asin2atan``, ``atan2acos``, ``atan2asin``

   - *Linear Algebra*

         * ``identity``, ``matrix``, ``makemat``, ``syst2mat``, ``matpow``

         * ``det``,  ``det_minor``, ``rank``, ``ker``, ``image``, ``rref``, ``simplex_reduce``,

         * ``egv``, ``egvl``,  ``eigenvalues``, ``pcar``, ``pcar_hessenberg``, ``pmin``,

         * ``jordan``, ``adjoint_matrix``, ``companion``, ``hessenberg``, ``transpose``,

         * ``cholesky``, ``lll``,  ``lu``, ``qr``, ``svd``, ``a2q``, ``gauss``, ``gramschmidt``,
           ``q2a``, ``isom``, ``mkisom``


   - *Finite Fieds*

         * ``%``, ``% 0``, ``mod``, ``GF``, ``powmod``


   - *Integers*

         * ``gcd``, ``iabcuv``, ``ichinrem``, ``idivis``, ``iegcd``,

         * ``ifactor``, ``ifactors``, ``iquo``, ``iquorem``, ``irem``,

         * ``is_prime, is_pseudoprime``, ``lcm``, ``mod``, ``nextprime``, ``pa2b2``, ``prevprime``,
           ``smod``, ``euler``, ``fracmod``

   - *List*

         * ``append``, ``accumulate_head_tail``, ``concat``, ``head``, ``makelist``, ``member``, ``mid``, ``revlist``, ``rotate``, ``shift``, ``size``, ``sizes``, ``sort``, ``suppress``, ``tail``

   - *Set*

         * ``intersect``, ``minus``, ``union``, ``is_element``, ``is_included``

   """
   return Pygen(s).eval()


#######################################
# A class to adjust giac configuration
#######################################
cdef class GiacSetting(Pygen):
     """
     A class to customise the Computer Algebra  System settings

     * property threads:

         Maximal number of allowed theads in giac

     * property digits:

         Default digits number used for approximations:

         ::

             sage: from giacpy import giacsettings,libgiac
             sage: giacsettings.digits=20;giacsettings.digits
             20
             sage: libgiac.approx('1/7')
             0.14285714285714285714
             sage: giacsettings.digits=12;

     * property sqrtflag:

         Flag to allow sqrt extractions during solve and factorizations:

         ::

             sage: from giacpy import libgiac,giacsettings
             sage: giacsettings.sqrtflag=False;giacsettings.sqrtflag
             False
             sage: libgiac('x**2-2').factor()
             x^2-2
             sage: giacsettings.sqrtflag=True;
             sage: libgiac('x**2-2').factor()
             (x-sqrt(2))*(x+sqrt(2))


     property complexflag:

         Flag to allow complex number in solving equations or factorizations:

         ::

            sage: from giacpy import libgiac,giacsettings
            sage: giacsettings.complexflag=False;giacsettings.complexflag
            False
            sage: libgiac('x**2+4').factor()
            x^2+4
            sage: giacsettings.complexflag=True;
            sage: libgiac('x**2+4').factor()
            (x+2*i)*(x-2*i)


     * property eval_level:

         Recursive level of substitution of variables during an evaluation:

         ::

             sage: from giacpy import giacsettings,libgiac
             sage: giacsettings.eval_level=1
             sage: libgiac("purge(a):;b:=a;a:=1;b")
             "Done",a,1,a
             sage: giacsettings.eval_level=25; giacsettings.eval_level
             25
             sage: libgiac("purge(a):;b:=a;a:=1;b")
             "Done",a,1,1

     * property proba_epsilon:

         Maximum probability of a wrong answer with a probabilist algorithm.
         Set this number to 0 to disable probabilist algorithms (slower):

         ::

             sage: from giacpy import giacsettings,libgiac
             sage: giacsettings.proba_epsilon=0;libgiac('proba_epsilon')
             0.0
             sage: giacsettings.proba_epsilon=10^(-13)
             sage: libgiac('proba_epsilon')<10^(-14)
             False

     * property epsilon:

         Value used by the epsilon2zero function:

         ::

             sage: from giacpy import giacsettings,libgiac
             sage: giacsettings.epsilon=1e-10
             sage: P=libgiac('1e-11+x+5')
             sage: P==x+5
             False
             sage: (P.epsilon2zero()).simplify()
             x+5

     """

     def __repr__(self):
        return "Giac Settings"

     def _sage_doc_(self):
        return GiacSetting.__doc__


     property digits:
         r"""
         Default digits number used for approximations:

         ::

            sage: from giacpy import giacsettings,libgiac
            sage: giacsettings.digits=20;giacsettings.digits
            20
            sage: libgiac.approx('1/7')
            0.14285714285714285714
            sage: giacsettings.digits=12;

         """
         def __get__(self):
           return (self.cas_setup()[6])._val


         def __set__(self,value):
           l=Pygen('cas_setup()').eval()
           pl=[ i for i in l ]
           pl[6]=value
           Pygen('cas_setup(%s)'%(pl)).eval()


     property sqrtflag:
         r"""
         Flag to allow square roots in solving equations or factorizations:
         """
         def __get__(self):
           return (self.cas_setup()[9])._val == 1


         def __set__(self,value):
           l=Pygen('cas_setup()').eval()
           pl=[ i for i in l ]
           if value:
              pl[9]=1
           else:
              pl[9]=0
           Pygen('cas_setup(%s)'%(pl)).eval()


     property complexflag:
         r"""
         Flag to allow complex number in solving equations or factorizations:

         ::

            sage: from giacpy import libgiac,giacsettings
            sage: giacsettings.complexflag=False;giacsettings.complexflag
            False
            sage: libgiac('x**2+4').factor()
            x^2+4
            sage: giacsettings.complexflag=True;
            sage: libgiac('x**2+4').factor()
            (x+2*i)*(x-2*i)
         """
         def __get__(self):
           return (self.cas_setup()[2])._val == 1


         def __set__(self,value):
           l=Pygen('cas_setup()').eval()
           pl=[ i for i in l ]
           if value:
              pl[2]=1
           else:
              pl[2]=0
           Pygen('cas_setup(%s)'%(pl)).eval()


     property eval_level:
         r"""
         Recursive level of substitution of variables during an evaluation:

         ::

            sage: from giacpy import giacsettings,libgiac
            sage: giacsettings.eval_level=1
            sage: libgiac("purge(a):;b:=a;a:=1;b")
            "Done",a,1,a
            sage: giacsettings.eval_level=25; giacsettings.eval_level
            25
            sage: libgiac("purge(a):;b:=a;a:=1;b")
            "Done",a,1,1

         """
         def __get__(self):
           return (self.cas_setup()[7][3])._val


         def __set__(self,value):
           l=Pygen('cas_setup()').eval()
           pl=[ i for i in l ]
           pl[7]=[l[7][0],l[7][1],l[7][2],value]
           Pygen('cas_setup(%s)'%(pl)).eval()



     property proba_epsilon:
         r"""
         Maximum probability of a wrong answer with a probabilist algorithm.
         Set this number to 0 to disable probabilist algorithms (slower):

         ::

            sage: from giacpy import giacsettings,libgiac
            sage: giacsettings.proba_epsilon=0;libgiac('proba_epsilon')
            0.0
            sage: giacsettings.proba_epsilon=10^(-13)
            sage: libgiac('proba_epsilon')<10^(-14)
            False

         """
         def __get__(self):
           return (self.cas_setup()[5][1])._double


         def __set__(self,value):
           l=Pygen('cas_setup()').eval()
           pl=[ i for i in l ]
           pl[5]=[l[5][0],value]
           Pygen('cas_setup(%s)'%(pl)).eval()


     property epsilon:
         r"""
         Value used by the epsilon2zero function:

         ::

            sage: from giacpy import giacsettings,libgiac
            sage: giacsettings.epsilon=1e-10
            sage: P=libgiac('1e-11+x+5')
            sage: P==x+5
            False
            sage: (P.epsilon2zero()).simplify()
            x+5

         """
         def __get__(self):
           return (self.cas_setup()[5][0])._double


         def __set__(self,value):
           l=Pygen('cas_setup()').eval()
           pl=[ i for i in l ]
           pl[5]=[value,l[5][1]]
           Pygen('cas_setup(%s)'%(pl)).eval()

     property threads:
         """
         Maximal number of allowed theads in giac
         """
         def __get__(self):
           return (self.cas_setup()[7][0])._val

         def __set__(self,value):
           Pygen('threads:=%s'%(str(value))).eval()


########################################################
#                                                      #
#    The python class that points to a cpp giac gen    #
#                                                      #
########################################################
cdef class Pygen:

     cdef gen * gptr   #pointer to the corresponding C++ element of type giac::gen

     def __cinit__(self, s=None):

         #NB: the  != here gives problems with  the __richcmp__ function
         #if (s!=None):
         # so it's better to use isinstance
         if (isinstance(s,None.__class__)):
             # Do NOT replace with: self=GIACNULL  (cf the doctest in __repr__
             sig_on()
             self.gptr = new gen ((<Pygen>GIACNULL).gptr[0])
             sig_off()
             return

         if isinstance(s,int):       # This looks 100 faster than the str initialisation
           if PythonVersion3:
             #in python3 int and long are int
             if s.bit_length()< Pymaxint.bit_length():
                sig_on()
                self.gptr = new gen(<long long>s)
                sig_off()
             else:
                sig_on()
                self.gptr = new gen(pylongtogen(s))
                sig_off()
           else:
                sig_on()
                self.gptr = new gen(<long long>s)
                sig_off()

         #
         elif isinstance(s,Integer):       # for sage int (gmp)
             sig_on()
             if (abs(s)>Pymaxint):
                 self.gptr = new gen((<Integer>s).value)
             else:
                self.gptr = new gen(<long long>s)   #important for pow to have a int
             sig_off()


         elif isinstance(s,Rational):       # for sage rat (gmp)
             #self.gptr = new gen((<Pygen>(Pygen(s.numerator())/Pygen(s.denominator()))).gptr[0])
             # FIXME: it's slow
             sig_on()
             self.gptr = new gen(GIAC_rdiv(gen((<Integer>(s.numerator())).value),gen((<Integer>(s.denominator())).value)))
             sig_off()


         elif is_Matrix(s):
             s = Pygen(s.list()).list2mat(s.ncols())
             sig_on()
             self.gptr = new gen((<Pygen>s).gptr[0])
             sig_off()


         elif isinstance(s,float):
             sig_on()
             self.gptr = new gen(<double>s)
             sig_off()


         elif isinstance(s,long):
             #s=str(s)
             #self.gptr = new gen(<string>s,context_ptr)
             sig_on()
             self.gptr = new gen(pylongtogen(s))
             sig_off()


         elif isinstance(s,Pygen):
             #in the test: x,y=Pygen('x,y');((x+2*y).sin()).texpand()
             # the y are lost without this case.
             sig_on()
             self.gptr = new gen((<Pygen>s).gptr[0])
             sig_off()

         elif isinstance(s,listrange):
             sig_on()
             self.gptr = new gen(_wrap_pylist(<list>s),<short int>0)
             sig_off()

         elif isinstance(s,tuple):
             sig_on()
             self.gptr = new gen(_wrap_pylist(<tuple>s),<short int>1)
             sig_off()

         # Other types are converted with strings.
         else:
           sig_on()
           if not(isinstance(s,str)):  #modif python3
             s=s.__str__()
           self.gptr = new gen(<string>encstring23(s),context_ptr)
           sig_off()


     def __dealloc__(self):
         del self.gptr


     def __repr__(self):
        """
         sage: from giacpy import *
         sage: vide=Pygen()
         sage: vide
         NULL

        """
        try:
          # fast evaluation of the complexity of the gen. (it's not the number of char )
          sig_on()
          t=GIAC_taille(self.gptr[0], 6000)
          sig_off()
        except:
          raise RuntimeError
        if (t<6000) :
           sig_on()
           result=decstring23(GIAC_print(self.gptr[0], context_ptr).c_str()) #python3
           sig_off()
           return result
        else:
           sig_on()
           result=str(self.type)+"\nResult is too big for Display. If you really want to see it use print"
           sig_off()
           return result

     def __str__(self):
        #if self.gptr == NULL:
        #  return ''
        sig_on()
        result=decstring23(GIAC_print(self.gptr[0], context_ptr).c_str()) #python3
        sig_off()
        return result

     def __len__(self):
       #if (self.gptr != NULL):
         try:
           sig_on()
           rep=GIAC_size(self.gptr[0],context_ptr).val
           sig_off()
           #GIAC_size return a gen. we take the int: val
           return rep
         except:
           raise RuntimeError
       #else:
       #  raise MemoryError,"This Pygen is not associated to a giac gen"


     def __getitem__(self,i):  #TODO?: add gen support for indexes
        """
        Lists of 10^6 integers should be translated to giac easily

           sage: from giacpy import libgiac
           sage: l=libgiac(list(range(10^6)));l[5]   #python3
           5
           sage: l[35:50:7]
           [35,42,49]
           sage: l[-10^6]
           0
           sage: t=libgiac(tuple(range(10)))
           sage: t[:4:-1]
           9,8,7,6,5
           sage: x=libgiac('x'); sum([ x[i] for i in range(5)])^3
           (x[0]+x[1]+x[2]+x[3]+x[4])^3
           sage: A=libgiac.ranm(5,10); A[3,7]-A[3][7]
           0
           sage: A.transpose()[8,2]-A[2][8]
           0

        Crash test:

           sage: from giacpy import *
           sage: l=Pygen()
           sage: l[0]
           Traceback (most recent call last):
             File "<stdin>", line 1, in <module>
           ...
           IndexError: list index 0 out of range

        """
        cdef gen result

        if(self._type == 7) or (self._type == 12):   #if self is a list or a string
          if isinstance(i,int) or isinstance(i,Integer):
            n=len(self)
            if(i<n)and(-i<=n):
                if(i<0):
                  i=i+n
                try:
                  sig_on()
                  result=self.gptr[0][<int>i]
                  sig_off()
                  return _wrap_gen(result)
                except:
                  raise RuntimeError
            else:
                raise IndexError,'list index %s out of range'%(i)
          else:
            if isinstance(i,slice):
              sig_on()
              result=gen(_getgiacslice(self,i),<short int>self._subtype)
              sig_off()
              return _wrap_gen(result)
            # add support for multi indexes
            elif isinstance(i,tuple):
               if(len(i)==2):
                  return self[i[0]][i[1]]
               elif(len(i)==1):
                  # in case of a tuple like this: (3,)
                  return self[i[0]]
               else:
                  return self[i[0],i[1]][tuple(i[2:])]
            else:
              raise TypeError,'gen indexes are not yet implemented'
        # Here we add support to formal variable indexes:
        else:
          cmd='%s[%s]'%(self,i)
          ans=Pygen(cmd).eval()
          # if the answer is a string, it must be an error message because self is not a list or a string
          if (ans._type == 12):
            raise TypeError, "Error executing code in Giac\nCODE:\n\t%s\nGiac ERROR:\n\t%s"%(cmd, ans)
          return ans




     def __iter__(self):
       """
         Pygen lists of 10^6 elements should be yield

           sage: from time import clock
           sage: from giacpy import libgiac
           sage: l=libgiac(range(10^6))
           sage: t=clock();l1=[ i for i in l ];t1=clock()-t;'time for a list of 10^6  i ',t1        # doctest: +ELLIPSIS
           ('time for a list of 10^6  i ', ...)
           sage: t=clock();l1=[ i+i for i in l ];t2=clock()-t;'time for a list of 10^6  i+i ',t2    # doctest: +ELLIPSIS
           ('time for a list of 10^6  i+i ', ...)
           sage: t=clock();l1=[ 1+i for i in l ];t3=clock()-t;"time for a list of 10^6  1+i ",t3    # doctest: +ELLIPSIS
           ('time for a list of 10^6  1+i ', ...)
           sage: t=clock();l1=[ i+1 for i in l ];t4=clock()-t;"time for a list of 10^6  i+1 ",t4    # doctest: +ELLIPSIS
           ('time for a list of 10^6  i+1 ', ...)
           sage: L=libgiac(range(10))
           sage: iter(L).next() # test  trac 18841
           0

       """
       cdef int i
       for i in range(len(self)):
          yield self[i]


     def eval(self):
        cdef gen result
        try:
           sig_on()
           result=GIAC_protecteval(self.gptr[0],giacsettings.eval_level,context_ptr)
           sig_off()
           return _wrap_gen(result)
        except:
           raise RuntimeError


     def __add__(self, right):
         cdef gen result
         if isinstance(right,Pygen)==False:
            right=Pygen(right)
         # Curiously this case is important:
         # otherwise: f=1/(2+sin(5*x)) crash
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             sig_on()
             result= (<Pygen>self).gptr[0] + (<Pygen>right).gptr[0]
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError


     def __call__(self, *args):
         """
         sage: from giacpy import libgiac
         """
         cdef gen result
         n=len(args)
         if (n>1):
           #FIXME? improve with a vector, or improve Pygen(list)
           right=Pygen(args).eval()
         elif (n==1):
           right=Pygen(args[0])
         else:
           right=GIACNULL
         if isinstance(self,Pygen)==False:
           self=Pygen(self)
#        Some giac errors such as pari_locked are caught by the try
#        so we can't put the sig_on() in the try.
#        But now a keyboard interrupt fall back to this sig_on so
#        it may have left the giac pari locked.
         sig_on()
         try:
             result= ((<Pygen>self).gptr[0])((<Pygen>right).gptr[0],context_ptr)
             sig_off()
             return _wrap_gen(result)
         except:
            # The previous computation might have failed due to a pari_lock
            # So we will not raise an exception yet.
            tmp=Pygen('pari_unlock()').eval()
            # if pari was not locked in giac, we have locked it, so unlock it.
            if(tmp==0):
                 Pygen('pari_unlock()').eval()
                 sig_off()
                 raise
            else:
               try:
                  result= GIAC_eval((<Pygen>right).gptr[0],<int>1,context_ptr)
                  result= ((<Pygen>self).gptr[0])(result,context_ptr)
                  sig_off()
                  return _wrap_gen(result)
               except:
                  sig_off()
                  raise


     def __sub__(self, right):
         cdef gen result
         if isinstance(right,Pygen)==False:
            right=Pygen(right)
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             sig_on()
             result= (<Pygen>self).gptr[0] - (<Pygen>right).gptr[0]
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError

     def __mul__(self, right):
         """
         sage: from giacpy import libgiac
         sage: (sqrt(5)*libgiac('x')).factor() # BUG test could give 0
         sqrt(5)*x
         sage: (libgiac('x')*sqrt(5)).factor()
         sqrt(5)*x
         """
         cdef gen result
         if isinstance(right,Pygen)==False:
            right=Pygen(right)
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             #result= (<Pygen>self).gptr[0] * (<Pygen>right).gptr[0]
             #NB: with the natural previous method, the following error generated by
             #giac causes python to quit instead of an error message.
             #l=Pygen([1,2]);l.transpose()*l;
             sig_on()
             result= GIAC_giacmul((<Pygen>self).gptr[0] , (<Pygen>right).gptr[0],context_ptr)
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError

#PB / in python3 is truediv
     def __div__(self, right):
         """
         sage: from giacpy import libgiac
         sage: (sqrt(3)/libgiac('x')).factor()   # BUG test could give 0
         sqrt(3)/x
         sage: (libgiac('x')/sqrt(3)).factor()
         sqrt(3)*x/3
         """
         cdef gen result
         if isinstance(right,Pygen)==False:
            right=Pygen(right)
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             sig_on()
             result= GIAC_giacdiv((<Pygen>self).gptr[0] , (<Pygen>right).gptr[0],context_ptr)
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError

     def __truediv__(self, right):
         cdef gen result
         if isinstance(right,Pygen)==False:
            right=Pygen(right)
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             sig_on()
             result= (<Pygen>self).gptr[0] / (<Pygen>right).gptr[0]
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError


     def __pow__(self, right ,ignored):
         cdef gen result
         if isinstance(right,Pygen)==False:
            right=Pygen(right)
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             sig_on()
             result= GIAC_pow((<Pygen>self).gptr[0],(<Pygen>right).gptr[0], context_ptr )
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError

     def __mod__(self, right):
         cdef gen result
         if not isinstance(right,Pygen):
            right=Pygen(right)
         if not isinstance(self,Pygen):
            self=Pygen(self)

         try:
             #result= gen(GIAC_makenewvecteur((<Pygen>self).gptr[0],(<Pygen>right).gptr[0]),<short int>1)
             #to have an integer output:
             #result= GIAC_smod(result,context_ptr)
             #we give a modular output:
             sig_on()
             result= GIAC_giacmod((<Pygen>self).gptr[0],(<Pygen>right).gptr[0],context_ptr)
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError

     def __neg__(self):
         cdef gen result
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             sig_on()
             result= GIAC_neg((<Pygen>self).gptr[0])
             sig_off()
             return _wrap_gen(result)
         except:
             raise RuntimeError

     def __pos__(self):
         return self

     # To be able to use the eval function before the GiacMethods initialisation
     def cas_setup(self,*args):
        return Pygen('cas_setup')(self,*args)


     def savegen(self, str filename):
        """
          Archive a Pygen element to a file in giac compressed format.

          Use the loadgiacgen command to get back the Pygen from the file.

          In C++ these files can be opened with giac::unarchive

            sage: from giacpy import *
            sage: f=libgiac('(x+y+z+2)**10'); g=f.normal()
            sage: g.savegen("fichiertest")           #  doctest: +SKIP
            sage: a=loadgiacgen("fichiertest")    #  doctest: +SKIP
            sage: from tempfile import NamedTemporaryFile
            sage: F=NamedTemporaryFile()   # chose a temporary file for a test
            sage: g.savegen(F.name)
            sage: a=loadgiacgen(F.name)
            sage: a.factor()
            (x+y+z+2)^10
            sage: F.close()

        """

        try:
            sig_on()
            GIAC_archive( <string>encstring23(filename), (<Pygen>self).gptr[0], context_ptr)
            sig_off()
        except:
             raise RuntimeError




     # def htmlhelp(self, str lang='en'):
     #     """
     #     Open the giac  html  detailled help about self in an external  browser

     #     There are currently 3 supported languages: 'en', 'fr', 'el'

     #      sage: from giacpy import gcd
     #      sage: gcd.htmlhelp('fr')           # doctest: +SKIP

     #     """
     #     l={'fr':1 , 'en':2, 'el':4}
     #     if (not lang in ['en', 'fr', 'el']):
     #       lang='en'
     #     try:
     #       url=decstring23(browser_help(self.gptr[0],l[lang])) #python3
     #       giacbasedir=decstring23(GIAC_giac_aide_dir())  # python3
     #     except:
     #       raise RuntimeError,'giac docs dir not found'
     #     print(url)
     #     if os.access(giacbasedir,os.F_OK):
     #        url='file:'+url
     #        wwwbrowseropen(url)



     def _help(self):
        return self.findhelp().__str__()

     def help(self):
        return self._help()

     def _sage_doc_(self):
        return self._help()


     def __doc__(self):
        return self._help()

     # # # # # # # # # # # # # # # # #
     # sage addons
     # # # # # # # # # # # # # # # # #
     def _latex_(self):
        r"""
        You can output Giac expressions in latex.

        ::

            sage: from giacpy import libgiac
            sage: print latex(libgiac('(x^4 - y)/(y^2-3*x)'))
            \frac{(x^{4}-y)}{(y^{2}-3\cdot x)}

        """
        sig_on()
        result=decstring23(GIAC_gen2tex(self.gptr[0], context_ptr).c_str()) #python3
        sig_off()
        return result


     def _integer_(self,Z=None):
        """
        Convert giac integers or modular integers to sage Integers (via gmp)

        ::

            sage: from giacpy import *
            sage: a=libgiac('10'); b=libgiac('2**300')
            sage: a;type(ZZ(a))
            10
            <type 'sage.rings.integer.Integer'>
            sage: next_prime(b)
            2037035976334486086268445688409378161051468393665936250636140449354381299763336706183397533
           sage: c=libgiac('2 % nextprime(2**40)')
           sage: ZZ(c^1000)
           -233775163595
           sage: Mod(2,next_prime(2^40))^1000 - ZZ(c^1000)
           0
           sage: 2^320-(c^320).sage()
           0

        """
        cdef Integer n = PY_NEW(Integer)
        typ = self._type

        if(typ == 0):
            # giac _INT_  i.e int
            return Integer(self._val)

        elif(typ == 2):
            # giac _ZINT  i.e mpz_t
            sig_on()
            mpz_set(n.value,(self.gptr.ref_ZINTptr())[0])
            sig_off()
            return n

        elif(typ == 15):
            # self is a giac modulo
            sig_on()
            a = _wrap_gen( (self.gptr.ref_MODptr())[0])
            # It is useless to get the modulus here
            # because the result will be lift to ZZ.
            result = ZZ(a)
            sig_off()
            return result

        else:
           raise TypeError, "Cannot convert non giac integers to Integer"


     def _rational_(self,Z=None):
        """
        Convert giac rationals to sage rationals

        ::

           sage: from giacpy import *
           sage: a=libgiac('103993/33102')
           sage: b=QQ(a); b
           103993/33102
           sage: b == a.sage()
           True

        """

        typ = self._type

        # _INT_ or _ZINT
        if(typ == 0 or typ == 2):
           return QQ(ZZ(self))
        # _FRAC_
        elif(typ == 10):
            # giac _RAT_
            try:
               result=ZZ(self.numer())/ZZ(self.denom())
               return result
            except:
               RuntimeError, "Failed to convert to QQ"
        else:
           raise TypeError, "Cannot convert non giac _FRAC_ to QQ"


     def sage(self):
        r"""
        Convert a libgiac expression back to a Sage expression. (could be slow)

        This currently does not implement a parser for the Giac output language,
        therefore only very simple expressions will convert successfully.

        Lists are converted recursively to sage.

        CURRENT STATUS:

           ZZ, QQ, ZZ/nZZ, strings, are supported, other type are sent to the symbolic ring
           via strings. In particular symbolic expressions modulo n should be lift to ZZ
           before ( with % 0 ).

        EXAMPLES::

           sage: from giacpy import libgiac
           sage: m = libgiac('x^2 + 5*y')
           sage: m.sage()
           x^2 + 5*y

        ::

           sage: m = libgiac('sin(2*sqrt(1-x^2)) * (1 - cos(1/x))^2')
           sage: m.trigexpand().sage()
           2*cos(sqrt(-x^2 + 1))*cos(1/x)^2*sin(sqrt(-x^2 + 1)) - 4*cos(sqrt(-x^2 + 1))*cos(1/x)*sin(sqrt(-x^2 + 1)) + 2*cos(sqrt(-x^2 + 1))*sin(sqrt(-x^2 + 1))

        ::

           sage: a=libgiac(' 2 % 7')
           sage: (a.sage())^6
           1
           sage: a=libgiac('"une chaine"')
           sage: b=a.sage(); b + b
           'une chaineune chaine'
           sage: isinstance(b,str)
           True

        """
        typ = self._type

        if (typ != 7) :
            # self is not a list
            if ( typ == 0 or typ == 2):
               return ZZ(self)

            elif (typ == 10):
               return QQ(self)

            elif (typ == 15):
               # modular integer
               sig_on()
               a = _wrap_gen( (self.gptr.ref_MODptr())[0])
               b = _wrap_gen( (self.gptr.ref_MODptr())[1])
               result = IntegerModRing(ZZ(b))(ZZ(a))
               sig_off()
               return result

            elif (typ == 12):
               # string
               sig_on()
               result=eval(self.__str__())
               sig_off()
               return result

            else:
               return SR(self)

        else:
            # self is a list
            sig_on()
            result=[entry.sage() for entry in self]
            sig_off()
            return result


     def _symbolic_(self, R):
        r"""
        Convert self object to the ring R via a basic string evaluation. (slow)


        EXAMPLES::

            sage: from giacpy import *
            sage: u,v=var('u,v');a=libgiac('cos(u+v)').texpand()
            sage: simplify(SR(a)+sin(u)*sin(v))
            cos(u)*cos(v)
        """
        try:
            result=R(self.__str__())
            return result

        except Exception:
            raise NotImplementedError("Unable to parse Giac output: %s" % self.__repr__())



     def _matrix_(self, R=ZZ):
        r"""
        Return matrix over the (Sage) ring R  where self
        should be a  Giac matrix. The default ring is ZZ.


        EXAMPLES::

            sage: from giacpy import *
            sage: R.<x,y>=QQ[]
            sage: M=libgiac('matrix(4,4,(k,l)->(x^k-y^l))'); M
            //...
            matrix[[0,1-y,1-y^2,1-y^3],[x-1,x-y,x-y^2,x-y^3],[x^2-1,x^2-y,x^2-y^2,x^2-y^3],[x^3-1,x^3-y,x^3-y^2,x^3-y^3]]
            sage: M.eigenvals()       # random
            0,0,(x^3+x^2+x-y^3-y^2-y+sqrt(x^6+2*x^5+3*x^4-14*x^3*y^3+2*x^3*y^2+2*x^3*y+6*x^3+2*x^2*y^3-14*x^2*y^2+2*x^2*y+5*x^2+2*x*y^3+2*x*y^2-14*x*y+4*x+y^6+2*y^5+3*y^4+6*y^3+5*y^2+4*y-12))/2,(x^3+x^2+x-y^3-y^2-y-sqrt(x^6+2*x^5+3*x^4-14*x^3*y^3+2*x^3*y^2+2*x^3*y+6*x^3+2*x^2*y^3-14*x^2*y^2+2*x^2*y+5*x^2+2*x*y^3+2*x*y^2-14*x*y+4*x+y^6+2*y^5+3*y^4+6*y^3+5*y^2+4*y-12))/2
            sage: Z=matrix(M,R);Z
            [         0     -y + 1   -y^2 + 1   -y^3 + 1]
            [     x - 1      x - y   -y^2 + x   -y^3 + x]
            [   x^2 - 1    x^2 - y  x^2 - y^2 -y^3 + x^2]
            [   x^3 - 1    x^3 - y  x^3 - y^2  x^3 - y^3]
            sage: parent(Z)
            Full MatrixSpace of 4 by 4 dense matrices over Multivariate Polynomial Ring in x, y over Rational Field
        """
        cdef int c
        cdef int r
        v = self.dim()
        n = (v[0])._val
        m = (v[1])._val
        from sage.matrix.matrix_space import MatrixSpace
        M = MatrixSpace(R, n, m)
        sig_on()
        entries = [[R((self[r])[c]) for c in range(m)] for r in range(n)]
        sig_off()
        return M(entries)

     def _vector_(self, R=None):
        r"""
        Return vector over the (Sage) ring R where self
        should be a  Giac matrix. The default ring is ZZ.


        EXAMPLES::

            sage: from giacpy import *
            sage: v=libgiac(range(10))
            sage: vector(v+v)
            (0, 2, 4, 6, 8, 10, 12, 14, 16, 18)
            sage: vector(v+v/3,QQ)
            (0, 4/3, 8/3, 4, 16/3, 20/3, 8, 28/3, 32/3, 12)


        """
        if(isinstance(R, None.__class__)):
            R=ZZ

        v = self.dim()
        try:
           n = v._val
        except:
           raise TypeError, "Entry is not a giac vector"
        from sage.modules.free_module_element import vector
        sig_on()
        entries = [R(self[c]) for c in range(n)]
        sig_off()
        return vector(R,entries)

     # # # # # # # # # # # # # # #


     def mplot(self):
        """
        Basic export of some 2D plots to sage. Only generic plots are supported.
        lines, circles, ... are not implemented
        """
        xyscat=[]
        xyplot=[]
        plotdata=self
        if not plotdata.type()=='DOM_LIST':
           plotdata=[plotdata]

        sig_on()
        for G in plotdata:
             if G.dim()>2:  # it is not a pnt. Ex: scatterplot
                for g in G:
                    xyscat=xyscat+[[(g.real())._double,(g.im())._double]]


             else:
                if G[1].type()=='DOM_LIST':
                   l=G[1].op()
                else:
                   l=G[1][2].op()
                xyplot=[[(u.real())._double,(u.im())._double] for u in l]


        if (xyscat != []):
           result=scatter_plot(xyscat)

        else:
           result=line(xyplot)
        sig_off()

        return result




     # # # # # # # # # # # # # # # # # # # # # # # # #
     #           WARNING:
     #
     # Do not use things like != in  Pygen's __cinit__
     # with this __richcmp__ enabled
     # The methods will bug: a=Pygen(2); a.sin()
     #
     # # # # # # # # # # # # # # # # # # # # # # # # #

     def __richcmp__( self, other,op):
         if isinstance(other,Pygen)==False:
            other=Pygen(other)
         if isinstance(self,Pygen)==False:
            self=Pygen(self)

         try:
             sig_on()
             result= giacgenrichcmp((<Pygen>self).gptr[0],(<Pygen>other).gptr[0], op, context_ptr )
             sig_off()
         except:
             raise RuntimeError
         if result==1 :
             return True
         else:
             return False

     #
     # Some attributes of the gen class:
     #
     property _type:
         def __get__(self):
            sig_on()
            result=self.gptr.type
            sig_off()
            return result


     property _subtype:
         def __get__(self):
            sig_on()
            result=self.gptr.subtype
            sig_off()
            return result



     property _val:  # immediate int (type _INT_)
         """
         immediate int value of an _INT_ type gen.
         """
         def __get__(self):
            if(self._type == 0):
               sig_on()
               result=self.gptr.val
               sig_off()
               return result
            else:
               raise TypeError,"Cannot convert non _INT_ giac gen"


     property _double:  # immediate double (type _DOUBLE_)
         """
         immediate conversion to float for a gen of _DOUBLE_ type.
         """
         def __get__(self):
            if(self._type == 1):
               sig_on()
               result=self.gptr._DOUBLE_val
               sig_off()
               return result
            else:
               raise TypeError,"Cannot convert non _DOUBLE_ giac gen"

     property help:
         def __get__(self):
           return self.help()



     ###################################################
     # Add the others methods
     ###################################################
     #
     #  NB: with __getattr__ this is about 10 times slower: [a.normal() for i in range(10**4)]
     #      than [GiacMethods["normal"](a) for i in range(10**4)]
     #
     #     def __getattr__(self, name):
     #       return GiacMethods[str(name)](self)
     ##
     def Airy_Ai(self,*args):
        return GiacMethods['Airy_Ai'](self,*args)

     def Airy_Bi(self,*args):
        return GiacMethods['Airy_Bi'](self,*args)

     def Archive(self,*args):
        return GiacMethods['Archive'](self,*args)

     def BesselJ(self,*args):
        return GiacMethods['BesselJ'](self,*args)

     def BesselY(self,*args):
        return GiacMethods['BesselY'](self,*args)

     def Beta(self,*args):
        return GiacMethods['Beta'](self,*args)

     def BlockDiagonal(self,*args):
        return GiacMethods['BlockDiagonal'](self,*args)

     def Ci(self,*args):
        return GiacMethods['Ci'](self,*args)

     def Circle(self,*args):
        return GiacMethods['Circle'](self,*args)

     def Col(self,*args):
        return GiacMethods['Col'](self,*args)

     def CopyVar(self,*args):
        return GiacMethods['CopyVar'](self,*args)

     def Dirac(self,*args):
        return GiacMethods['Dirac'](self,*args)

     def Ei(self,*args):
        return GiacMethods['Ei'](self,*args)

     def Factor(self,*args):
        return GiacMethods['Factor'](self,*args)

     def GF(self,*args):
        return GiacMethods['GF'](self,*args)

     def Gamma(self,*args):
        return GiacMethods['Gamma'](self,*args)

     def Heaviside(self,*args):
        return GiacMethods['Heaviside'](self,*args)

     def JordanBlock(self,*args):
        return GiacMethods['JordanBlock'](self,*args)

     def LU(self,*args):
        return GiacMethods['LU'](self,*args)

     def Li(self,*args):
        return GiacMethods['Li'](self,*args)

     def Line(self,*args):
        return GiacMethods['Line'](self,*args)

     def LineHorz(self,*args):
        return GiacMethods['LineHorz'](self,*args)

     def LineTan(self,*args):
        return GiacMethods['LineTan'](self,*args)

     def LineVert(self,*args):
        return GiacMethods['LineVert'](self,*args)

     def Phi(self,*args):
        return GiacMethods['Phi'](self,*args)

     def Pi(self,*args):
        return GiacMethods['Pi'](self,*args)

     def Psi(self,*args):
        return GiacMethods['Psi'](self,*args)

     def QR(self,*args):
        return GiacMethods['QR'](self,*args)

     def RandSeed(self,*args):
        return GiacMethods['RandSeed'](self,*args)

     def Row(self,*args):
        return GiacMethods['Row'](self,*args)

     def SortA(self,*args):
        return GiacMethods['SortA'](self,*args)

     def SortD(self,*args):
        return GiacMethods['SortD'](self,*args)

     def UTPC(self,*args):
        return GiacMethods['UTPC'](self,*args)

     def UTPF(self,*args):
        return GiacMethods['UTPF'](self,*args)

     def UTPN(self,*args):
        return GiacMethods['UTPN'](self,*args)

     def UTPT(self,*args):
        return GiacMethods['UTPT'](self,*args)

     def VARS(self,*args):
        return GiacMethods['VARS'](self,*args)

     def VAS(self,*args):
        return GiacMethods['VAS'](self,*args)

     def VAS_positive(self,*args):
        return GiacMethods['VAS_positive'](self,*args)

     def Zeta(self,*args):
        return GiacMethods['Zeta'](self,*args)

     def a2q(self,*args):
        return GiacMethods['a2q'](self,*args)

     def abcuv(self,*args):
        return GiacMethods['abcuv'](self,*args)

     def about(self,*args):
        return GiacMethods['about'](self,*args)

     def abs(self,*args):
        return GiacMethods['abs'](self,*args)

     def abscissa(self,*args):
        return GiacMethods['abscissa'](self,*args)

     def accumulate_head_tail(self,*args):
        return GiacMethods['accumulate_head_tail'](self,*args)

     def acos(self,*args):
        return GiacMethods['acos'](self,*args)

     def acos2asin(self,*args):
        return GiacMethods['acos2asin'](self,*args)

     def acos2atan(self,*args):
        return GiacMethods['acos2atan'](self,*args)

     def acosh(self,*args):
        return GiacMethods['acosh'](self,*args)

     def acot(self,*args):
        return GiacMethods['acot'](self,*args)

     def acsc(self,*args):
        return GiacMethods['acsc'](self,*args)

     def add(self,*args):
        return GiacMethods['add'](self,*args)

     def additionally(self,*args):
        return GiacMethods['additionally'](self,*args)

     def adjoint_matrix(self,*args):
        return GiacMethods['adjoint_matrix'](self,*args)

     def affix(self,*args):
        return GiacMethods['affix'](self,*args)

     def algsubs(self,*args):
        return GiacMethods['algsubs'](self,*args)

     def algvar(self,*args):
        return GiacMethods['algvar'](self,*args)

     def alog10(self,*args):
        return GiacMethods['alog10'](self,*args)

     def altitude(self,*args):
        return GiacMethods['altitude'](self,*args)

     def angle(self,*args):
        return GiacMethods['angle'](self,*args)

     def angle_radian(self,*args):
        return GiacMethods['angle_radian'](self,*args)

     def angleat(self,*args):
        return GiacMethods['angleat'](self,*args)

     def angleatraw(self,*args):
        return GiacMethods['angleatraw'](self,*args)

     def ans(self,*args):
        return GiacMethods['ans'](self,*args)

     def apply(self,*args):
        return GiacMethods['apply'](self,*args)

     def approx(self,*args):
        return GiacMethods['approx'](self,*args)

     def arc(self,*args):
        return GiacMethods['arc'](self,*args)

     def arcLen(self,*args):
        return GiacMethods['arcLen'](self,*args)

     def arccos(self,*args):
        return GiacMethods['arccos'](self,*args)

     def arccosh(self,*args):
        return GiacMethods['arccosh'](self,*args)

     def arclen(self,*args):
        return GiacMethods['arclen'](self,*args)

     def arcsin(self,*args):
        return GiacMethods['arcsin'](self,*args)

     def arcsinh(self,*args):
        return GiacMethods['arcsinh'](self,*args)

     def arctan(self,*args):
        return GiacMethods['arctan'](self,*args)

     def arctanh(self,*args):
        return GiacMethods['arctanh'](self,*args)

     def area(self,*args):
        return GiacMethods['area'](self,*args)

     def areaat(self,*args):
        return GiacMethods['areaat'](self,*args)

     def areaatraw(self,*args):
        return GiacMethods['areaatraw'](self,*args)

     def areaplot(self,*args):
        return GiacMethods['areaplot'](self,*args)

     def arg(self,*args):
        return GiacMethods['arg'](self,*args)

     def array(self,*args):
        return GiacMethods['array'](self,*args)

     def asin(self,*args):
        return GiacMethods['asin'](self,*args)

     def asin2acos(self,*args):
        return GiacMethods['asin2acos'](self,*args)

     def asin2atan(self,*args):
        return GiacMethods['asin2atan'](self,*args)

     def asinh(self,*args):
        return GiacMethods['asinh'](self,*args)

     def assume(self,*args):
        return GiacMethods['assume'](self,*args)

     def at(self,*args):
        return GiacMethods['at'](self,*args)

     def atan(self,*args):
        return GiacMethods['atan'](self,*args)

     def atan2acos(self,*args):
        return GiacMethods['atan2acos'](self,*args)

     def atan2asin(self,*args):
        return GiacMethods['atan2asin'](self,*args)

     def atanh(self,*args):
        return GiacMethods['atanh'](self,*args)

     def atrig2ln(self,*args):
        return GiacMethods['atrig2ln'](self,*args)

     def augment(self,*args):
        return GiacMethods['augment'](self,*args)

     def autosimplify(self,*args):
        return GiacMethods['autosimplify'](self,*args)

     def avance(self,*args):
        return GiacMethods['avance'](self,*args)

     def avgRC(self,*args):
        return GiacMethods['avgRC'](self,*args)

     def axes(self,*args):
        return GiacMethods['axes'](self,*args)

     def back(self,*args):
        return GiacMethods['back'](self,*args)

     def baisse_crayon(self,*args):
        return GiacMethods['baisse_crayon'](self,*args)

     def bar_plot(self,*args):
        return GiacMethods['bar_plot'](self,*args)

     def barycenter(self,*args):
        return GiacMethods['barycenter'](self,*args)

     def base(self,*args):
        return GiacMethods['base'](self,*args)

     def basis(self,*args):
        return GiacMethods['basis'](self,*args)

     def batons(self,*args):
        return GiacMethods['batons'](self,*args)

     def bernoulli(self,*args):
        return GiacMethods['bernoulli'](self,*args)

     def besselJ(self,*args):
        return GiacMethods['besselJ'](self,*args)

     def besselY(self,*args):
        return GiacMethods['besselY'](self,*args)

     def betad(self,*args):
        return GiacMethods['betad'](self,*args)

     def betad_cdf(self,*args):
        return GiacMethods['betad_cdf'](self,*args)

     def betad_icdf(self,*args):
        return GiacMethods['betad_icdf'](self,*args)

     def bezier(self,*args):
        return GiacMethods['bezier'](self,*args)

     def bezout_entiers(self,*args):
        return GiacMethods['bezout_entiers'](self,*args)

     def binomial(self,*args):
        return GiacMethods['binomial'](self,*args)

     def binomial_cdf(self,*args):
        return GiacMethods['binomial_cdf'](self,*args)

     def binomial_icdf(self,*args):
        return GiacMethods['binomial_icdf'](self,*args)

     def bisection_solver(self,*args):
        return GiacMethods['bisection_solver'](self,*args)

     def bisector(self,*args):
        return GiacMethods['bisector'](self,*args)

     def bitand(self,*args):
        return GiacMethods['bitand'](self,*args)

     def bitor(self,*args):
        return GiacMethods['bitor'](self,*args)

     def bitxor(self,*args):
        return GiacMethods['bitxor'](self,*args)

     def blockmatrix(self,*args):
        return GiacMethods['blockmatrix'](self,*args)

     def border(self,*args):
        return GiacMethods['border'](self,*args)

     def boxwhisker(self,*args):
        return GiacMethods['boxwhisker'](self,*args)

     def brent_solver(self,*args):
        return GiacMethods['brent_solver'](self,*args)

     def cFactor(self,*args):
        return GiacMethods['cFactor'](self,*args)

     def cSolve(self,*args):
        return GiacMethods['cSolve'](self,*args)

     def cZeros(self,*args):
        return GiacMethods['cZeros'](self,*args)

     def camembert(self,*args):
        return GiacMethods['camembert'](self,*args)

     def canonical_form(self,*args):
        return GiacMethods['canonical_form'](self,*args)

     def cauchy(self,*args):
        return GiacMethods['cauchy'](self,*args)

     def cauchy_cdf(self,*args):
        return GiacMethods['cauchy_cdf'](self,*args)

     def cauchy_icdf(self,*args):
        return GiacMethods['cauchy_icdf'](self,*args)

     def cauchyd(self,*args):
        return GiacMethods['cauchyd'](self,*args)

     def cauchyd_cdf(self,*args):
        return GiacMethods['cauchyd_cdf'](self,*args)

     def cauchyd_icdf(self,*args):
        return GiacMethods['cauchyd_icdf'](self,*args)

     def cdf(self,*args):
        return GiacMethods['cdf'](self,*args)

     def ceil(self,*args):
        return GiacMethods['ceil'](self,*args)

     def ceiling(self,*args):
        return GiacMethods['ceiling'](self,*args)

     def center(self,*args):
        return GiacMethods['center'](self,*args)

     def center2interval(self,*args):
        return GiacMethods['center2interval'](self,*args)

     def centered_cube(self,*args):
        return GiacMethods['centered_cube'](self,*args)

     def centered_tetrahedron(self,*args):
        return GiacMethods['centered_tetrahedron'](self,*args)

     def cfactor(self,*args):
        return GiacMethods['cfactor'](self,*args)

     def cfsolve(self,*args):
        return GiacMethods['cfsolve'](self,*args)

     def changebase(self,*args):
        return GiacMethods['changebase'](self,*args)

     def char(self,*args):
        return GiacMethods['char'](self,*args)

     def charpoly(self,*args):
        return GiacMethods['charpoly'](self,*args)

     def chinrem(self,*args):
        return GiacMethods['chinrem'](self,*args)

     def chisquare(self,*args):
        return GiacMethods['chisquare'](self,*args)

     def chisquare_cdf(self,*args):
        return GiacMethods['chisquare_cdf'](self,*args)

     def chisquare_icdf(self,*args):
        return GiacMethods['chisquare_icdf'](self,*args)

     def chisquared(self,*args):
        return GiacMethods['chisquared'](self,*args)

     def chisquared_cdf(self,*args):
        return GiacMethods['chisquared_cdf'](self,*args)

     def chisquared_icdf(self,*args):
        return GiacMethods['chisquared_icdf'](self,*args)

     def chisquaret(self,*args):
        return GiacMethods['chisquaret'](self,*args)

     def cholesky(self,*args):
        return GiacMethods['cholesky'](self,*args)

     def chrem(self,*args):
        return GiacMethods['chrem'](self,*args)

     def circle(self,*args):
        return GiacMethods['circle'](self,*args)

     def circumcircle(self,*args):
        return GiacMethods['circumcircle'](self,*args)

     def classes(self,*args):
        return GiacMethods['classes'](self,*args)

     def coeff(self,*args):
        return GiacMethods['coeff'](self,*args)

     def coeffs(self,*args):
        return GiacMethods['coeffs'](self,*args)

     def col(self,*args):
        return GiacMethods['col'](self,*args)

     def colDim(self,*args):
        return GiacMethods['colDim'](self,*args)

     def colNorm(self,*args):
        return GiacMethods['colNorm'](self,*args)

     def colSwap(self,*args):
        return GiacMethods['colSwap'](self,*args)

     def coldim(self,*args):
        return GiacMethods['coldim'](self,*args)

     def collect(self,*args):
        return GiacMethods['collect'](self,*args)

     def colnorm(self,*args):
        return GiacMethods['colnorm'](self,*args)

     def color(self,*args):
        return GiacMethods['color'](self,*args)

     def colspace(self,*args):
        return GiacMethods['colspace'](self,*args)

     def colswap(self,*args):
        return GiacMethods['colswap'](self,*args)

     def comDenom(self,*args):
        return GiacMethods['comDenom'](self,*args)

     def comb(self,*args):
        return GiacMethods['comb'](self,*args)

     def combine(self,*args):
        return GiacMethods['combine'](self,*args)

     def comment(self,*args):
        return GiacMethods['comment'](self,*args)

     def common_perpendicular(self,*args):
        return GiacMethods['common_perpendicular'](self,*args)

     def companion(self,*args):
        return GiacMethods['companion'](self,*args)

     def compare(self,*args):
        return GiacMethods['compare'](self,*args)

     def complex(self,*args):
        return GiacMethods['complex'](self,*args)

     def complex_variables(self,*args):
        return GiacMethods['complex_variables'](self,*args)

     def complexroot(self,*args):
        return GiacMethods['complexroot'](self,*args)

     def concat(self,*args):
        return GiacMethods['concat'](self,*args)

     def cond(self,*args):
        return GiacMethods['cond'](self,*args)

     def cone(self,*args):
        return GiacMethods['cone'](self,*args)

     def confrac(self,*args):
        return GiacMethods['confrac'](self,*args)

     def conic(self,*args):
        return GiacMethods['conic'](self,*args)

     def conj(self,*args):
        return GiacMethods['conj'](self,*args)

     def conjugate_gradient(self,*args):
        return GiacMethods['conjugate_gradient'](self,*args)

     def cont(self,*args):
        return GiacMethods['cont'](self,*args)

     def contains(self,*args):
        return GiacMethods['contains'](self,*args)

     def content(self,*args):
        return GiacMethods['content'](self,*args)

     def contourplot(self,*args):
        return GiacMethods['contourplot'](self,*args)

     def convert(self,*args):
        return GiacMethods['convert'](self,*args)

     def convertir(self,*args):
        return GiacMethods['convertir'](self,*args)

     def convexhull(self,*args):
        return GiacMethods['convexhull'](self,*args)

     def coordinates(self,*args):
        return GiacMethods['coordinates'](self,*args)

     def copy(self,*args):
        return GiacMethods['copy'](self,*args)

     def correlation(self,*args):
        return GiacMethods['correlation'](self,*args)

     def cos(self,*args):
        return GiacMethods['cos'](self,*args)

     def cos2sintan(self,*args):
        return GiacMethods['cos2sintan'](self,*args)

     def cosh(self,*args):
        return GiacMethods['cosh'](self,*args)

     def cot(self,*args):
        return GiacMethods['cot'](self,*args)

     def cote(self,*args):
        return GiacMethods['cote'](self,*args)

     def count(self,*args):
        return GiacMethods['count'](self,*args)

     def count_eq(self,*args):
        return GiacMethods['count_eq'](self,*args)

     def count_inf(self,*args):
        return GiacMethods['count_inf'](self,*args)

     def count_sup(self,*args):
        return GiacMethods['count_sup'](self,*args)

     def courbe_parametrique(self,*args):
        return GiacMethods['courbe_parametrique'](self,*args)

     def courbe_polaire(self,*args):
        return GiacMethods['courbe_polaire'](self,*args)

     def covariance(self,*args):
        return GiacMethods['covariance'](self,*args)

     def covariance_correlation(self,*args):
        return GiacMethods['covariance_correlation'](self,*args)

     def cpartfrac(self,*args):
        return GiacMethods['cpartfrac'](self,*args)

     def crationalroot(self,*args):
        return GiacMethods['crationalroot'](self,*args)

     def crayon(self,*args):
        return GiacMethods['crayon'](self,*args)

     def cross(self,*args):
        return GiacMethods['cross'](self,*args)

     def crossP(self,*args):
        return GiacMethods['crossP'](self,*args)

     def cross_point(self,*args):
        return GiacMethods['cross_point'](self,*args)

     def cross_ratio(self,*args):
        return GiacMethods['cross_ratio'](self,*args)

     def crossproduct(self,*args):
        return GiacMethods['crossproduct'](self,*args)

     def csc(self,*args):
        return GiacMethods['csc'](self,*args)

     def csolve(self,*args):
        return GiacMethods['csolve'](self,*args)

     def cube(self,*args):
        return GiacMethods['cube'](self,*args)

     def cumSum(self,*args):
        return GiacMethods['cumSum'](self,*args)

     def cumsum(self,*args):
        return GiacMethods['cumsum'](self,*args)

     def cumulated_frequencies(self,*args):
        return GiacMethods['cumulated_frequencies'](self,*args)

     def curl(self,*args):
        return GiacMethods['curl'](self,*args)

     def current_sheet(self,*args):
        return GiacMethods['current_sheet'](self,*args)

     def curvature(self,*args):
        return GiacMethods['curvature'](self,*args)

     def curve(self,*args):
        return GiacMethods['curve'](self,*args)

     def cyan(self,*args):
        return GiacMethods['cyan'](self,*args)

     def cycle2perm(self,*args):
        return GiacMethods['cycle2perm'](self,*args)

     def cycleinv(self,*args):
        return GiacMethods['cycleinv'](self,*args)

     def cycles2permu(self,*args):
        return GiacMethods['cycles2permu'](self,*args)

     def cyclotomic(self,*args):
        return GiacMethods['cyclotomic'](self,*args)

     def cylinder(self,*args):
        return GiacMethods['cylinder'](self,*args)

     def dash_line(self,*args):
        return GiacMethods['dash_line'](self,*args)

     def dashdot_line(self,*args):
        return GiacMethods['dashdot_line'](self,*args)

     def dashdotdot_line(self,*args):
        return GiacMethods['dashdotdot_line'](self,*args)

     def dayofweek(self,*args):
        return GiacMethods['dayofweek'](self,*args)

     def deSolve(self,*args):
        return GiacMethods['deSolve'](self,*args)

     def debut_enregistrement(self,*args):
        return GiacMethods['debut_enregistrement'](self,*args)

     def degree(self,*args):
        return GiacMethods['degree'](self,*args)

     def delcols(self,*args):
        return GiacMethods['delcols'](self,*args)

     def delrows(self,*args):
        return GiacMethods['delrows'](self,*args)

     def deltalist(self,*args):
        return GiacMethods['deltalist'](self,*args)

     def denom(self,*args):
        return GiacMethods['denom'](self,*args)

     def densityplot(self,*args):
        return GiacMethods['densityplot'](self,*args)

     def derive(self,*args):
        return GiacMethods['derive'](self,*args)

     def deriver(self,*args):
        return GiacMethods['deriver'](self,*args)

     def desolve(self,*args):
        return GiacMethods['desolve'](self,*args)

     def dessine_tortue(self,*args):
        return GiacMethods['dessine_tortue'](self,*args)

     def det(self,*args):
        return GiacMethods['det'](self,*args)

     def det_minor(self,*args):
        return GiacMethods['det_minor'](self,*args)

     def developper(self,*args):
        return GiacMethods['developper'](self,*args)

     def developper_transcendant(self,*args):
        return GiacMethods['developper_transcendant'](self,*args)

     def dfc(self,*args):
        return GiacMethods['dfc'](self,*args)

     def dfc2f(self,*args):
        return GiacMethods['dfc2f'](self,*args)

     def diag(self,*args):
        return GiacMethods['diag'](self,*args)

     def diag(self,*args):
        return GiacMethods['diag'](self,*args)

     def diff(self,*args):
        return GiacMethods['diff'](self,*args)

     def dim(self,*args):
        return GiacMethods['dim'](self,*args)

     def display(self,*args):
        return GiacMethods['display'](self,*args)

     def disque(self,*args):
        return GiacMethods['disque'](self,*args)

     def disque_centre(self,*args):
        return GiacMethods['disque_centre'](self,*args)

     def distance(self,*args):
        return GiacMethods['distance'](self,*args)

     def distance2(self,*args):
        return GiacMethods['distance2'](self,*args)

     def distanceat(self,*args):
        return GiacMethods['distanceat'](self,*args)

     def distanceatraw(self,*args):
        return GiacMethods['distanceatraw'](self,*args)

     def divergence(self,*args):
        return GiacMethods['divergence'](self,*args)

     def divide(self,*args):
        return GiacMethods['divide'](self,*args)

     def divis(self,*args):
        return GiacMethods['divis'](self,*args)

     def division_point(self,*args):
        return GiacMethods['division_point'](self,*args)

     def divisors(self,*args):
        return GiacMethods['divisors'](self,*args)

     def divpc(self,*args):
        return GiacMethods['divpc'](self,*args)

     def dnewton_solver(self,*args):
        return GiacMethods['dnewton_solver'](self,*args)

     def dodecahedron(self,*args):
        return GiacMethods['dodecahedron'](self,*args)

     def dot(self,*args):
        return GiacMethods['dot'](self,*args)

     def dotP(self,*args):
        return GiacMethods['dotP'](self,*args)

     def dot_paper(self,*args):
        return GiacMethods['dot_paper'](self,*args)

     def dotprod(self,*args):
        return GiacMethods['dotprod'](self,*args)

     def droit(self,*args):
        return GiacMethods['droit'](self,*args)

     def droite_tangente(self,*args):
        return GiacMethods['droite_tangente'](self,*args)

     def dsolve(self,*args):
        return GiacMethods['dsolve'](self,*args)

     def e(self,*args):
        return GiacMethods['e'](self,*args)

     def e2r(self,*args):
        return GiacMethods['e2r'](self,*args)

     def ecart_type(self,*args):
        return GiacMethods['ecart_type'](self,*args)

     def ecart_type_population(self,*args):
        return GiacMethods['ecart_type_population'](self,*args)

     def egcd(self,*args):
        return GiacMethods['egcd'](self,*args)

     def egv(self,*args):
        return GiacMethods['egv'](self,*args)

     def egvl(self,*args):
        return GiacMethods['egvl'](self,*args)

     def eigVc(self,*args):
        return GiacMethods['eigVc'](self,*args)

     def eigVl(self,*args):
        return GiacMethods['eigVl'](self,*args)

     def eigenvals(self,*args):
        return GiacMethods['eigenvals'](self,*args)

     def eigenvalues(self,*args):
        return GiacMethods['eigenvalues'](self,*args)

     def eigenvectors(self,*args):
        return GiacMethods['eigenvectors'](self,*args)

     def eigenvects(self,*args):
        return GiacMethods['eigenvects'](self,*args)

     def element(self,*args):
        return GiacMethods['element'](self,*args)

     def eliminate(self,*args):
        return GiacMethods['eliminate'](self,*args)

     def ellipse(self,*args):
        return GiacMethods['ellipse'](self,*args)

     def entry(self,*args):
        return GiacMethods['entry'](self,*args)

     def envelope(self,*args):
        return GiacMethods['envelope'](self,*args)

     def epsilon(self,*args):
        return GiacMethods['epsilon'](self,*args)

     def epsilon2zero(self,*args):
        return GiacMethods['epsilon2zero'](self,*args)

     def equal(self,*args):
        return GiacMethods['equal'](self,*args)

     def equal2diff(self,*args):
        return GiacMethods['equal2diff'](self,*args)

     def equal2list(self,*args):
        return GiacMethods['equal2list'](self,*args)

     def equation(self,*args):
        return GiacMethods['equation'](self,*args)

     def equilateral_triangle(self,*args):
        return GiacMethods['equilateral_triangle'](self,*args)

     def erf(self,*args):
        return GiacMethods['erf'](self,*args)

     def erfc(self,*args):
        return GiacMethods['erfc'](self,*args)

     def error(self,*args):
        return GiacMethods['error'](self,*args)

     def euler(self,*args):
        return GiacMethods['euler'](self,*args)

     def euler_gamma(self,*args):
        return GiacMethods['euler_gamma'](self,*args)

     def eval_level(self,*args):
        return GiacMethods['eval_level'](self,*args)

     def evala(self,*args):
        return GiacMethods['evala'](self,*args)

     def evalb(self,*args):
        return GiacMethods['evalb'](self,*args)

     def evalc(self,*args):
        return GiacMethods['evalc'](self,*args)

     def evalf(self,*args):
        return GiacMethods['evalf'](self,*args)

     def evalm(self,*args):
        return GiacMethods['evalm'](self,*args)

     def even(self,*args):
        return GiacMethods['even'](self,*args)

     def evolute(self,*args):
        return GiacMethods['evolute'](self,*args)

     def exact(self,*args):
        return GiacMethods['exact'](self,*args)

     def exbisector(self,*args):
        return GiacMethods['exbisector'](self,*args)

     def excircle(self,*args):
        return GiacMethods['excircle'](self,*args)

     def execute(self,*args):
        return GiacMethods['execute'](self,*args)

     def exp(self,*args):
        return GiacMethods['exp'](self,*args)

     def exp2list(self,*args):
        return GiacMethods['exp2list'](self,*args)

     def exp2pow(self,*args):
        return GiacMethods['exp2pow'](self,*args)

     def exp2trig(self,*args):
        return GiacMethods['exp2trig'](self,*args)

     def expand(self,*args):
        return GiacMethods['expand'](self,*args)

     def expexpand(self,*args):
        return GiacMethods['expexpand'](self,*args)

     def expln(self,*args):
        return GiacMethods['expln'](self,*args)

     def exponential(self,*args):
        return GiacMethods['exponential'](self,*args)

     def exponential_cdf(self,*args):
        return GiacMethods['exponential_cdf'](self,*args)

     def exponential_icdf(self,*args):
        return GiacMethods['exponential_icdf'](self,*args)

     def exponential_regression(self,*args):
        return GiacMethods['exponential_regression'](self,*args)

     def exponential_regression_plot(self,*args):
        return GiacMethods['exponential_regression_plot'](self,*args)

     def exponentiald(self,*args):
        return GiacMethods['exponentiald'](self,*args)

     def exponentiald_cdf(self,*args):
        return GiacMethods['exponentiald_cdf'](self,*args)

     def exponentiald_icdf(self,*args):
        return GiacMethods['exponentiald_icdf'](self,*args)

     def expr(self,*args):
        return GiacMethods['expr'](self,*args)

     def extract_measure(self,*args):
        return GiacMethods['extract_measure'](self,*args)

     def ezgcd(self,*args):
        return GiacMethods['ezgcd'](self,*args)

     def f2nd(self,*args):
        return GiacMethods['f2nd'](self,*args)

     def fMax(self,*args):
        return GiacMethods['fMax'](self,*args)

     def fMin(self,*args):
        return GiacMethods['fMin'](self,*args)

     def fPart(self,*args):
        return GiacMethods['fPart'](self,*args)

     def faces(self,*args):
        return GiacMethods['faces'](self,*args)

     def facteurs_premiers(self,*args):
        return GiacMethods['facteurs_premiers'](self,*args)

     def factor(self,*args):
        return GiacMethods['factor'](self,*args)

     def factor_xn(self,*args):
        return GiacMethods['factor_xn'](self,*args)

     def factorial(self,*args):
        return GiacMethods['factorial'](self,*args)

     def factoriser(self,*args):
        return GiacMethods['factoriser'](self,*args)

     def factoriser_entier(self,*args):
        return GiacMethods['factoriser_entier'](self,*args)

     def factoriser_sur_C(self,*args):
        return GiacMethods['factoriser_sur_C'](self,*args)

     def factors(self,*args):
        return GiacMethods['factors'](self,*args)

     def fadeev(self,*args):
        return GiacMethods['fadeev'](self,*args)

     def false(self,*args):
        return GiacMethods['false'](self,*args)

     def falsepos_solver(self,*args):
        return GiacMethods['falsepos_solver'](self,*args)

     def fclose(self,*args):
        return GiacMethods['fclose'](self,*args)

     def fcoeff(self,*args):
        return GiacMethods['fcoeff'](self,*args)

     def fdistrib(self,*args):
        return GiacMethods['fdistrib'](self,*args)

     def fft(self,*args):
        return GiacMethods['fft'](self,*args)

     def fieldplot(self,*args):
        return GiacMethods['fieldplot'](self,*args)

     def find(self,*args):
        return GiacMethods['find'](self,*args)

     def findhelp(self,*args):
        return GiacMethods['findhelp'](self,*args)

     def fisher(self,*args):
        return GiacMethods['fisher'](self,*args)

     def fisher_cdf(self,*args):
        return GiacMethods['fisher_cdf'](self,*args)

     def fisher_icdf(self,*args):
        return GiacMethods['fisher_icdf'](self,*args)

     def fisherd(self,*args):
        return GiacMethods['fisherd'](self,*args)

     def fisherd_cdf(self,*args):
        return GiacMethods['fisherd_cdf'](self,*args)

     def fisherd_icdf(self,*args):
        return GiacMethods['fisherd_icdf'](self,*args)

     def flatten(self,*args):
        return GiacMethods['flatten'](self,*args)

     def float2rational(self,*args):
        return GiacMethods['float2rational'](self,*args)

     def floor(self,*args):
        return GiacMethods['floor'](self,*args)

     def fonction_derivee(self,*args):
        return GiacMethods['fonction_derivee'](self,*args)

     def fourier_an(self,*args):
        return GiacMethods['fourier_an'](self,*args)

     def fourier_bn(self,*args):
        return GiacMethods['fourier_bn'](self,*args)

     def fourier_cn(self,*args):
        return GiacMethods['fourier_cn'](self,*args)

     def fprint(self,*args):
        return GiacMethods['fprint'](self,*args)

     def frac(self,*args):
        return GiacMethods['frac'](self,*args)

     def fracmod(self,*args):
        return GiacMethods['fracmod'](self,*args)

     def frame_2d(self,*args):
        return GiacMethods['frame_2d'](self,*args)

     def frequencies(self,*args):
        return GiacMethods['frequencies'](self,*args)

     def frobenius_norm(self,*args):
        return GiacMethods['frobenius_norm'](self,*args)

     def froot(self,*args):
        return GiacMethods['froot'](self,*args)

     def fsolve(self,*args):
        return GiacMethods['fsolve'](self,*args)

     def fullparfrac(self,*args):
        return GiacMethods['fullparfrac'](self,*args)

     def funcplot(self,*args):
        return GiacMethods['funcplot'](self,*args)

     def function_diff(self,*args):
        return GiacMethods['function_diff'](self,*args)

     def fxnd(self,*args):
        return GiacMethods['fxnd'](self,*args)

     def gammad(self,*args):
        return GiacMethods['gammad'](self,*args)

     def gammad_cdf(self,*args):
        return GiacMethods['gammad_cdf'](self,*args)

     def gammad_icdf(self,*args):
        return GiacMethods['gammad_icdf'](self,*args)

     def gauss(self,*args):
        return GiacMethods['gauss'](self,*args)

     def gauss15(self,*args):
        return GiacMethods['gauss15'](self,*args)

     def gauss_seidel_linsolve(self,*args):
        return GiacMethods['gauss_seidel_linsolve'](self,*args)

     def gaussjord(self,*args):
        return GiacMethods['gaussjord'](self,*args)

     def gaussquad(self,*args):
        return GiacMethods['gaussquad'](self,*args)

     def gbasis(self,*args):
        return GiacMethods['gbasis'](self,*args)

     def gcd(self,*args):
        return GiacMethods['gcd'](self,*args)

     def gcdex(self,*args):
        return GiacMethods['gcdex'](self,*args)

     def genpoly(self,*args):
        return GiacMethods['genpoly'](self,*args)

     def geometric(self,*args):
        return GiacMethods['geometric'](self,*args)

     def geometric_cdf(self,*args):
        return GiacMethods['geometric_cdf'](self,*args)

     def geometric_icdf(self,*args):
        return GiacMethods['geometric_icdf'](self,*args)

     def getDenom(self,*args):
        return GiacMethods['getDenom'](self,*args)

     def getKey(self,*args):
        return GiacMethods['getKey'](self,*args)

     def getNum(self,*args):
        return GiacMethods['getNum'](self,*args)

     def getType(self,*args):
        return GiacMethods['getType'](self,*args)

     def greduce(self,*args):
        return GiacMethods['greduce'](self,*args)

     def groupermu(self,*args):
        return GiacMethods['groupermu'](self,*args)

     def hadamard(self,*args):
        return GiacMethods['hadamard'](self,*args)

     def half_cone(self,*args):
        return GiacMethods['half_cone'](self,*args)

     def half_line(self,*args):
        return GiacMethods['half_line'](self,*args)

     def halftan(self,*args):
        return GiacMethods['halftan'](self,*args)

     def halftan_hyp2exp(self,*args):
        return GiacMethods['halftan_hyp2exp'](self,*args)

     def halt(self,*args):
        return GiacMethods['halt'](self,*args)

     def hamdist(self,*args):
        return GiacMethods['hamdist'](self,*args)

     def harmonic_conjugate(self,*args):
        return GiacMethods['harmonic_conjugate'](self,*args)

     def harmonic_division(self,*args):
        return GiacMethods['harmonic_division'](self,*args)

     def has(self,*args):
        return GiacMethods['has'](self,*args)

     def hasard(self,*args):
        return GiacMethods['hasard'](self,*args)

     def head(self,*args):
        return GiacMethods['head'](self,*args)

     def hermite(self,*args):
        return GiacMethods['hermite'](self,*args)

     def hessenberg(self,*args):
        return GiacMethods['hessenberg'](self,*args)

     def hessian(self,*args):
        return GiacMethods['hessian'](self,*args)

     def heugcd(self,*args):
        return GiacMethods['heugcd'](self,*args)

     def hexagon(self,*args):
        return GiacMethods['hexagon'](self,*args)

     def hilbert(self,*args):
        return GiacMethods['hilbert'](self,*args)

     def histogram(self,*args):
        return GiacMethods['histogram'](self,*args)

     def hold(self,*args):
        return GiacMethods['hold'](self,*args)

     def homothety(self,*args):
        return GiacMethods['homothety'](self,*args)

     def horner(self,*args):
        return GiacMethods['horner'](self,*args)

     def hybrid_solver(self,*args):
        return GiacMethods['hybrid_solver'](self,*args)

     def hybridj_solver(self,*args):
        return GiacMethods['hybridj_solver'](self,*args)

     def hybrids_solver(self,*args):
        return GiacMethods['hybrids_solver'](self,*args)

     def hybridsj_solver(self,*args):
        return GiacMethods['hybridsj_solver'](self,*args)

     def hyp2exp(self,*args):
        return GiacMethods['hyp2exp'](self,*args)

     def hyperbola(self,*args):
        return GiacMethods['hyperbola'](self,*args)

     def iPart(self,*args):
        return GiacMethods['iPart'](self,*args)

     def iabcuv(self,*args):
        return GiacMethods['iabcuv'](self,*args)

     def ibasis(self,*args):
        return GiacMethods['ibasis'](self,*args)

     def ibpdv(self,*args):
        return GiacMethods['ibpdv'](self,*args)

     def ibpu(self,*args):
        return GiacMethods['ibpu'](self,*args)

     def icdf(self,*args):
        return GiacMethods['icdf'](self,*args)

     def ichinrem(self,*args):
        return GiacMethods['ichinrem'](self,*args)

     def ichrem(self,*args):
        return GiacMethods['ichrem'](self,*args)

     def icontent(self,*args):
        return GiacMethods['icontent'](self,*args)

     def icosahedron(self,*args):
        return GiacMethods['icosahedron'](self,*args)

     def id(self,*args):
        return GiacMethods['id'](self,*args)

     def identity(self,*args):
        return GiacMethods['identity'](self,*args)

     def idivis(self,*args):
        return GiacMethods['idivis'](self,*args)

     def idn(self,*args):
        return GiacMethods['idn'](self,*args)

     def iegcd(self,*args):
        return GiacMethods['iegcd'](self,*args)

     def ifactor(self,*args):
        return GiacMethods['ifactor'](self,*args)

     def ifactors(self,*args):
        return GiacMethods['ifactors'](self,*args)

     def igamma(self,*args):
        return GiacMethods['igamma'](self,*args)

     def igcd(self,*args):
        return GiacMethods['igcd'](self,*args)

     def igcdex(self,*args):
        return GiacMethods['igcdex'](self,*args)

     def ihermite(self,*args):
        return GiacMethods['ihermite'](self,*args)

     def ilaplace(self,*args):
        return GiacMethods['ilaplace'](self,*args)

     def im(self,*args):
        return GiacMethods['im'](self,*args)

     def imag(self,*args):
        return GiacMethods['imag'](self,*args)

     def image(self,*args):
        return GiacMethods['image'](self,*args)

     def implicitplot(self,*args):
        return GiacMethods['implicitplot'](self,*args)

     def inString(self,*args):
        return GiacMethods['inString'](self,*args)

     def in_ideal(self,*args):
        return GiacMethods['in_ideal'](self,*args)

     def incircle(self,*args):
        return GiacMethods['incircle'](self,*args)

     def indets(self,*args):
        return GiacMethods['indets'](self,*args)

     def inequationplot(self,*args):
        return GiacMethods['inequationplot'](self,*args)

     def inf(self,*args):
        return GiacMethods['inf'](self,*args)

     def infinity(self,*args):
        return GiacMethods['infinity'](self,*args)

     def infnorm(self,*args):
        return GiacMethods['infnorm'](self,*args)

     def insmod(self,*args):
        return GiacMethods['insmod'](self,*args)

     def int(self,*args):
        return GiacMethods['int'](self,*args)

     def intDiv(self,*args):
        return GiacMethods['intDiv'](self,*args)

     def integer(self,*args):
        return GiacMethods['integer'](self,*args)

     def integrate(self,*args):
        return GiacMethods['integrate'](self,*args)

     def integrer(self,*args):
        return GiacMethods['integrer'](self,*args)

     def inter(self,*args):
        return GiacMethods['inter'](self,*args)

     def interactive_odeplot(self,*args):
        return GiacMethods['interactive_odeplot'](self,*args)

     def interactive_plotode(self,*args):
        return GiacMethods['interactive_plotode'](self,*args)

     def interp(self,*args):
        return GiacMethods['interp'](self,*args)

     def interval(self,*args):
        return GiacMethods['interval'](self,*args)

     def interval2center(self,*args):
        return GiacMethods['interval2center'](self,*args)

     def inv(self,*args):
        return GiacMethods['inv'](self,*args)

     def inverse(self,*args):
        return GiacMethods['inverse'](self,*args)

     def inversion(self,*args):
        return GiacMethods['inversion'](self,*args)

     def invisible_point(self,*args):
        return GiacMethods['invisible_point'](self,*args)

     def invlaplace(self,*args):
        return GiacMethods['invlaplace'](self,*args)

     def invztrans(self,*args):
        return GiacMethods['invztrans'](self,*args)

     def iquo(self,*args):
        return GiacMethods['iquo'](self,*args)

     def iquorem(self,*args):
        return GiacMethods['iquorem'](self,*args)

     def iratrecon(self,*args):
        return GiacMethods['iratrecon'](self,*args)

     def irem(self,*args):
        return GiacMethods['irem'](self,*args)

     def isPrime(self,*args):
        return GiacMethods['isPrime'](self,*args)

     def is_collinear(self,*args):
        return GiacMethods['is_collinear'](self,*args)

     def is_concyclic(self,*args):
        return GiacMethods['is_concyclic'](self,*args)

     def is_conjugate(self,*args):
        return GiacMethods['is_conjugate'](self,*args)

     def is_coplanar(self,*args):
        return GiacMethods['is_coplanar'](self,*args)

     def is_cospheric(self,*args):
        return GiacMethods['is_cospheric'](self,*args)

     def is_cycle(self,*args):
        return GiacMethods['is_cycle'](self,*args)

     def is_element(self,*args):
        return GiacMethods['is_element'](self,*args)

     def is_equilateral(self,*args):
        return GiacMethods['is_equilateral'](self,*args)

     def is_harmonic(self,*args):
        return GiacMethods['is_harmonic'](self,*args)

     def is_harmonic_circle_bundle(self,*args):
        return GiacMethods['is_harmonic_circle_bundle'](self,*args)

     def is_harmonic_line_bundle(self,*args):
        return GiacMethods['is_harmonic_line_bundle'](self,*args)

     def is_inside(self,*args):
        return GiacMethods['is_inside'](self,*args)

     def is_isosceles(self,*args):
        return GiacMethods['is_isosceles'](self,*args)

     def is_orthogonal(self,*args):
        return GiacMethods['is_orthogonal'](self,*args)

     def is_parallel(self,*args):
        return GiacMethods['is_parallel'](self,*args)

     def is_parallelogram(self,*args):
        return GiacMethods['is_parallelogram'](self,*args)

     def is_permu(self,*args):
        return GiacMethods['is_permu'](self,*args)

     def is_perpendicular(self,*args):
        return GiacMethods['is_perpendicular'](self,*args)

     def is_prime(self,*args):
        return GiacMethods['is_prime'](self,*args)

     def is_pseudoprime(self,*args):
        return GiacMethods['is_pseudoprime'](self,*args)

     def is_rectangle(self,*args):
        return GiacMethods['is_rectangle'](self,*args)

     def is_rhombus(self,*args):
        return GiacMethods['is_rhombus'](self,*args)

     def is_square(self,*args):
        return GiacMethods['is_square'](self,*args)

     def ismith(self,*args):
        return GiacMethods['ismith'](self,*args)

     def isobarycenter(self,*args):
        return GiacMethods['isobarycenter'](self,*args)

     def isom(self,*args):
        return GiacMethods['isom'](self,*args)

     def isopolygon(self,*args):
        return GiacMethods['isopolygon'](self,*args)

     def isosceles_triangle(self,*args):
        return GiacMethods['isosceles_triangle'](self,*args)

     def isprime(self,*args):
        return GiacMethods['isprime'](self,*args)

     def ithprime(self,*args):
        return GiacMethods['ithprime'](self,*args)

     def jacobi_linsolve(self,*args):
        return GiacMethods['jacobi_linsolve'](self,*args)

     def jacobi_symbol(self,*args):
        return GiacMethods['jacobi_symbol'](self,*args)

     def jordan(self,*args):
        return GiacMethods['jordan'](self,*args)

     def keep_pivot(self,*args):
        return GiacMethods['keep_pivot'](self,*args)

     def ker(self,*args):
        return GiacMethods['ker'](self,*args)

     def kernel(self,*args):
        return GiacMethods['kernel'](self,*args)

     def kolmogorovd(self,*args):
        return GiacMethods['kolmogorovd'](self,*args)

     def kolmogorovt(self,*args):
        return GiacMethods['kolmogorovt'](self,*args)

     def l1norm(self,*args):
        return GiacMethods['l1norm'](self,*args)

     def l2norm(self,*args):
        return GiacMethods['l2norm'](self,*args)

     def lagrange(self,*args):
        return GiacMethods['lagrange'](self,*args)

     def laguerre(self,*args):
        return GiacMethods['laguerre'](self,*args)

     def laplace(self,*args):
        return GiacMethods['laplace'](self,*args)

     def laplacian(self,*args):
        return GiacMethods['laplacian'](self,*args)

     def latex(self,*args):
        return GiacMethods['latex'](self,*args)

     def lcm(self,*args):
        return GiacMethods['lcm'](self,*args)

     def lcoeff(self,*args):
        return GiacMethods['lcoeff'](self,*args)

     def ldegree(self,*args):
        return GiacMethods['ldegree'](self,*args)

     def left(self,*args):
        return GiacMethods['left'](self,*args)

     def left_rectangle(self,*args):
        return GiacMethods['left_rectangle'](self,*args)

     def legend(self,*args):
        return GiacMethods['legend'](self,*args)

     def legendre(self,*args):
        return GiacMethods['legendre'](self,*args)

     def legendre_symbol(self,*args):
        return GiacMethods['legendre_symbol'](self,*args)

     def length(self,*args):
        return GiacMethods['length'](self,*args)

     def lgcd(self,*args):
        return GiacMethods['lgcd'](self,*args)

     def lhs(self,*args):
        return GiacMethods['lhs'](self,*args)

     def ligne_chapeau_carre(self,*args):
        return GiacMethods['ligne_chapeau_carre'](self,*args)

     def ligne_chapeau_plat(self,*args):
        return GiacMethods['ligne_chapeau_plat'](self,*args)

     def ligne_chapeau_rond(self,*args):
        return GiacMethods['ligne_chapeau_rond'](self,*args)

     def ligne_polygonale(self,*args):
        return GiacMethods['ligne_polygonale'](self,*args)

     def ligne_polygonale_pointee(self,*args):
        return GiacMethods['ligne_polygonale_pointee'](self,*args)

     def ligne_tiret(self,*args):
        return GiacMethods['ligne_tiret'](self,*args)

     def ligne_tiret_point(self,*args):
        return GiacMethods['ligne_tiret_point'](self,*args)

     def ligne_tiret_pointpoint(self,*args):
        return GiacMethods['ligne_tiret_pointpoint'](self,*args)

     def ligne_trait_plein(self,*args):
        return GiacMethods['ligne_trait_plein'](self,*args)

     def limit(self,*args):
        return GiacMethods['limit'](self,*args)

     def limite(self,*args):
        return GiacMethods['limite'](self,*args)

     def lin(self,*args):
        return GiacMethods['lin'](self,*args)

     def line(self,*args):
        return GiacMethods['line'](self,*args)

     def line_inter(self,*args):
        return GiacMethods['line_inter'](self,*args)

     def line_paper(self,*args):
        return GiacMethods['line_paper'](self,*args)

     def line_segments(self,*args):
        return GiacMethods['line_segments'](self,*args)

     def linear_interpolate(self,*args):
        return GiacMethods['linear_interpolate'](self,*args)

     def linear_regression(self,*args):
        return GiacMethods['linear_regression'](self,*args)

     def linear_regression_plot(self,*args):
        return GiacMethods['linear_regression_plot'](self,*args)

     def lineariser(self,*args):
        return GiacMethods['lineariser'](self,*args)

     def lineariser_trigo(self,*args):
        return GiacMethods['lineariser_trigo'](self,*args)

     def linsolve(self,*args):
        return GiacMethods['linsolve'](self,*args)

     def linspace(self,*args):
        return GiacMethods['linspace'](self,*args)

     def lis_phrase(self,*args):
        return GiacMethods['lis_phrase'](self,*args)

     def list2mat(self,*args):
        return GiacMethods['list2mat'](self,*args)

     def listplot(self,*args):
        return GiacMethods['listplot'](self,*args)

     def lll(self,*args):
        return GiacMethods['lll'](self,*args)

     def ln(self,*args):
        return GiacMethods['ln'](self,*args)

     def lname(self,*args):
        return GiacMethods['lname'](self,*args)

     def lncollect(self,*args):
        return GiacMethods['lncollect'](self,*args)

     def lnexpand(self,*args):
        return GiacMethods['lnexpand'](self,*args)

     def locus(self,*args):
        return GiacMethods['locus'](self,*args)

     def log(self,*args):
        return GiacMethods['log'](self,*args)

     def log10(self,*args):
        return GiacMethods['log10'](self,*args)

     def logarithmic_regression(self,*args):
        return GiacMethods['logarithmic_regression'](self,*args)

     def logarithmic_regression_plot(self,*args):
        return GiacMethods['logarithmic_regression_plot'](self,*args)

     def logb(self,*args):
        return GiacMethods['logb'](self,*args)

     def logistic_regression(self,*args):
        return GiacMethods['logistic_regression'](self,*args)

     def logistic_regression_plot(self,*args):
        return GiacMethods['logistic_regression_plot'](self,*args)

     def lsmod(self,*args):
        return GiacMethods['lsmod'](self,*args)

     def lsq(self,*args):
        return GiacMethods['lsq'](self,*args)

     def lu(self,*args):
        return GiacMethods['lu'](self,*args)

     def lvar(self,*args):
        return GiacMethods['lvar'](self,*args)

     def mRow(self,*args):
        return GiacMethods['mRow'](self,*args)

     def mRowAdd(self,*args):
        return GiacMethods['mRowAdd'](self,*args)

     def magenta(self,*args):
        return GiacMethods['magenta'](self,*args)

     def makelist(self,*args):
        return GiacMethods['makelist'](self,*args)

     def makemat(self,*args):
        return GiacMethods['makemat'](self,*args)

     def makesuite(self,*args):
        return GiacMethods['makesuite'](self,*args)

     def makevector(self,*args):
        return GiacMethods['makevector'](self,*args)

     def map(self,*args):
        return GiacMethods['map'](self,*args)

     def maple2mupad(self,*args):
        return GiacMethods['maple2mupad'](self,*args)

     def maple2xcas(self,*args):
        return GiacMethods['maple2xcas'](self,*args)

     def maple_ifactors(self,*args):
        return GiacMethods['maple_ifactors'](self,*args)

     def maple_mode(self,*args):
        return GiacMethods['maple_mode'](self,*args)

     def markov(self,*args):
        return GiacMethods['markov'](self,*args)

     def mat2list(self,*args):
        return GiacMethods['mat2list'](self,*args)

     def mathml(self,*args):
        return GiacMethods['mathml'](self,*args)

     def matpow(self,*args):
        return GiacMethods['matpow'](self,*args)

     def matrix(self,*args):
        return GiacMethods['matrix'](self,*args)

     def matrix_norm(self,*args):
        return GiacMethods['matrix_norm'](self,*args)

     def max(self,*args):
        return GiacMethods['max'](self,*args)

     def maxnorm(self,*args):
        return GiacMethods['maxnorm'](self,*args)

     def mean(self,*args):
        return GiacMethods['mean'](self,*args)

     def median(self,*args):
        return GiacMethods['median'](self,*args)

     def median_line(self,*args):
        return GiacMethods['median_line'](self,*args)

     def member(self,*args):
        return GiacMethods['member'](self,*args)

     def mgf(self,*args):
        return GiacMethods['mgf'](self,*args)

     def mid(self,*args):
        return GiacMethods['mid'](self,*args)

     def middle_point(self,*args):
        return GiacMethods['middle_point'](self,*args)

     def midpoint(self,*args):
        return GiacMethods['midpoint'](self,*args)

     def min(self,*args):
        return GiacMethods['min'](self,*args)

     def mkisom(self,*args):
        return GiacMethods['mkisom'](self,*args)

     def mksa(self,*args):
        return GiacMethods['mksa'](self,*args)

     def modgcd(self,*args):
        return GiacMethods['modgcd'](self,*args)

     def mods(self,*args):
        return GiacMethods['mods'](self,*args)

     def montre_tortue(self,*args):
        return GiacMethods['montre_tortue'](self,*args)

     def moustache(self,*args):
        return GiacMethods['moustache'](self,*args)

     def moyal(self,*args):
        return GiacMethods['moyal'](self,*args)

     def moyenne(self,*args):
        return GiacMethods['moyenne'](self,*args)

     def mul(self,*args):
        return GiacMethods['mul'](self,*args)

     def mult_c_conjugate(self,*args):
        return GiacMethods['mult_c_conjugate'](self,*args)

     def mult_conjugate(self,*args):
        return GiacMethods['mult_conjugate'](self,*args)

     def multinomial(self,*args):
        return GiacMethods['multinomial'](self,*args)

     def multiplier_conjugue(self,*args):
        return GiacMethods['multiplier_conjugue'](self,*args)

     def multiplier_conjugue_complexe(self,*args):
        return GiacMethods['multiplier_conjugue_complexe'](self,*args)

     def multiply(self,*args):
        return GiacMethods['multiply'](self,*args)

     def mupad2maple(self,*args):
        return GiacMethods['mupad2maple'](self,*args)

     def mupad2xcas(self,*args):
        return GiacMethods['mupad2xcas'](self,*args)

     def nCr(self,*args):
        return GiacMethods['nCr'](self,*args)

     def nDeriv(self,*args):
        return GiacMethods['nDeriv'](self,*args)

     def nInt(self,*args):
        return GiacMethods['nInt'](self,*args)

     def nPr(self,*args):
        return GiacMethods['nPr'](self,*args)

     def nSolve(self,*args):
        return GiacMethods['nSolve'](self,*args)

     def ncols(self,*args):
        return GiacMethods['ncols'](self,*args)

     def negbinomial(self,*args):
        return GiacMethods['negbinomial'](self,*args)

     def negbinomial_cdf(self,*args):
        return GiacMethods['negbinomial_cdf'](self,*args)

     def negbinomial_icdf(self,*args):
        return GiacMethods['negbinomial_icdf'](self,*args)

     def newList(self,*args):
        return GiacMethods['newList'](self,*args)

     def newMat(self,*args):
        return GiacMethods['newMat'](self,*args)

     def newton(self,*args):
        return GiacMethods['newton'](self,*args)

     def newton_solver(self,*args):
        return GiacMethods['newton_solver'](self,*args)

     def newtonj_solver(self,*args):
        return GiacMethods['newtonj_solver'](self,*args)

     def nextperm(self,*args):
        return GiacMethods['nextperm'](self,*args)

     def nextprime(self,*args):
        return GiacMethods['nextprime'](self,*args)

     def nodisp(self,*args):
        return GiacMethods['nodisp'](self,*args)

     def non_recursive_normal(self,*args):
        return GiacMethods['non_recursive_normal'](self,*args)

     def nop(self,*args):
        return GiacMethods['nop'](self,*args)

     def nops(self,*args):
        return GiacMethods['nops'](self,*args)

     def norm(self,*args):
        return GiacMethods['norm'](self,*args)

     def normal(self,*args):
        return GiacMethods['normal'](self,*args)

     def normal_cdf(self,*args):
        return GiacMethods['normal_cdf'](self,*args)

     def normal_icdf(self,*args):
        return GiacMethods['normal_icdf'](self,*args)

     def normald(self,*args):
        return GiacMethods['normald'](self,*args)

     def normald_cdf(self,*args):
        return GiacMethods['normald_cdf'](self,*args)

     def normald_icdf(self,*args):
        return GiacMethods['normald_icdf'](self,*args)

     def normalize(self,*args):
        return GiacMethods['normalize'](self,*args)

     def normalt(self,*args):
        return GiacMethods['normalt'](self,*args)

     def nprimes(self,*args):
        return GiacMethods['nprimes'](self,*args)

     def nrows(self,*args):
        return GiacMethods['nrows'](self,*args)

     def nuage_points(self,*args):
        return GiacMethods['nuage_points'](self,*args)

     def nullspace(self,*args):
        return GiacMethods['nullspace'](self,*args)

     def numer(self,*args):
        return GiacMethods['numer'](self,*args)

     def octahedron(self,*args):
        return GiacMethods['octahedron'](self,*args)

     def odd(self,*args):
        return GiacMethods['odd'](self,*args)

     def odeplot(self,*args):
        return GiacMethods['odeplot'](self,*args)

     def odesolve(self,*args):
        return GiacMethods['odesolve'](self,*args)

     def op(self,*args):
        return GiacMethods['op'](self,*args)

     def open_polygon(self,*args):
        return GiacMethods['open_polygon'](self,*args)

     def ord(self,*args):
        return GiacMethods['ord'](self,*args)

     def order_size(self,*args):
        return GiacMethods['order_size'](self,*args)

     def ordinate(self,*args):
        return GiacMethods['ordinate'](self,*args)

     def orthocenter(self,*args):
        return GiacMethods['orthocenter'](self,*args)

     def orthogonal(self,*args):
        return GiacMethods['orthogonal'](self,*args)

     def osculating_circle(self,*args):
        return GiacMethods['osculating_circle'](self,*args)

     def p1oc2(self,*args):
        return GiacMethods['p1oc2'](self,*args)

     def p1op2(self,*args):
        return GiacMethods['p1op2'](self,*args)

     def pa2b2(self,*args):
        return GiacMethods['pa2b2'](self,*args)

     def pade(self,*args):
        return GiacMethods['pade'](self,*args)

     def parabola(self,*args):
        return GiacMethods['parabola'](self,*args)

     def parallel(self,*args):
        return GiacMethods['parallel'](self,*args)

     def parallelepiped(self,*args):
        return GiacMethods['parallelepiped'](self,*args)

     def parallelogram(self,*args):
        return GiacMethods['parallelogram'](self,*args)

     def parameq(self,*args):
        return GiacMethods['parameq'](self,*args)

     def parameter(self,*args):
        return GiacMethods['parameter'](self,*args)

     def paramplot(self,*args):
        return GiacMethods['paramplot'](self,*args)

     def parfrac(self,*args):
        return GiacMethods['parfrac'](self,*args)

     def pari(self,*args):
        return GiacMethods['pari'](self,*args)

     def part(self,*args):
        return GiacMethods['part'](self,*args)

     def partfrac(self,*args):
        return GiacMethods['partfrac'](self,*args)

     def pas_de_cote(self,*args):
        return GiacMethods['pas_de_cote'](self,*args)

     def pcar(self,*args):
        return GiacMethods['pcar'](self,*args)

     def pcar_hessenberg(self,*args):
        return GiacMethods['pcar_hessenberg'](self,*args)

     def pcoef(self,*args):
        return GiacMethods['pcoef'](self,*args)

     def pcoeff(self,*args):
        return GiacMethods['pcoeff'](self,*args)

     def perimeter(self,*args):
        return GiacMethods['perimeter'](self,*args)

     def perimeterat(self,*args):
        return GiacMethods['perimeterat'](self,*args)

     def perimeteratraw(self,*args):
        return GiacMethods['perimeteratraw'](self,*args)

     def perm(self,*args):
        return GiacMethods['perm'](self,*args)

     def perminv(self,*args):
        return GiacMethods['perminv'](self,*args)

     def permu2cycles(self,*args):
        return GiacMethods['permu2cycles'](self,*args)

     def permu2mat(self,*args):
        return GiacMethods['permu2mat'](self,*args)

     def permuorder(self,*args):
        return GiacMethods['permuorder'](self,*args)

     def perpen_bisector(self,*args):
        return GiacMethods['perpen_bisector'](self,*args)

     def perpendicular(self,*args):
        return GiacMethods['perpendicular'](self,*args)

     def peval(self,*args):
        return GiacMethods['peval'](self,*args)

     def pi(self,*args):
        return GiacMethods['pi'](self,*args)

     def piecewise(self,*args):
        return GiacMethods['piecewise'](self,*args)

     def pivot(self,*args):
        return GiacMethods['pivot'](self,*args)

     def pixoff(self,*args):
        return GiacMethods['pixoff'](self,*args)

     def pixon(self,*args):
        return GiacMethods['pixon'](self,*args)

     def plane(self,*args):
        return GiacMethods['plane'](self,*args)

     def playsnd(self,*args):
        return GiacMethods['playsnd'](self,*args)

     def plex(self,*args):
        return GiacMethods['plex'](self,*args)

     def plot(self,*args):
        return GiacMethods['plot'](self,*args)

     def plot3d(self,*args):
        return GiacMethods['plot3d'](self,*args)

     def plotarea(self,*args):
        return GiacMethods['plotarea'](self,*args)

     def plotcdf(self,*args):
        return GiacMethods['plotcdf'](self,*args)

     def plotcontour(self,*args):
        return GiacMethods['plotcontour'](self,*args)

     def plotdensity(self,*args):
        return GiacMethods['plotdensity'](self,*args)

     def plotfield(self,*args):
        return GiacMethods['plotfield'](self,*args)

     def plotfunc(self,*args):
        return GiacMethods['plotfunc'](self,*args)

     def plotimplicit(self,*args):
        return GiacMethods['plotimplicit'](self,*args)

     def plotinequation(self,*args):
        return GiacMethods['plotinequation'](self,*args)

     def plotlist(self,*args):
        return GiacMethods['plotlist'](self,*args)

     def plotode(self,*args):
        return GiacMethods['plotode'](self,*args)

     def plotparam(self,*args):
        return GiacMethods['plotparam'](self,*args)

     def plotpolar(self,*args):
        return GiacMethods['plotpolar'](self,*args)

     def plotseq(self,*args):
        return GiacMethods['plotseq'](self,*args)

     def plus_point(self,*args):
        return GiacMethods['plus_point'](self,*args)

     def pmin(self,*args):
        return GiacMethods['pmin'](self,*args)

     def point(self,*args):
        return GiacMethods['point'](self,*args)

     def point2d(self,*args):
        return GiacMethods['point2d'](self,*args)

     def point3d(self,*args):
        return GiacMethods['point3d'](self,*args)

     def poisson(self,*args):
        return GiacMethods['poisson'](self,*args)

     def poisson_cdf(self,*args):
        return GiacMethods['poisson_cdf'](self,*args)

     def poisson_icdf(self,*args):
        return GiacMethods['poisson_icdf'](self,*args)

     def polar(self,*args):
        return GiacMethods['polar'](self,*args)

     def polar_coordinates(self,*args):
        return GiacMethods['polar_coordinates'](self,*args)

     def polar_point(self,*args):
        return GiacMethods['polar_point'](self,*args)

     def polarplot(self,*args):
        return GiacMethods['polarplot'](self,*args)

     def pole(self,*args):
        return GiacMethods['pole'](self,*args)

     def poly2symb(self,*args):
        return GiacMethods['poly2symb'](self,*args)

     def polyEval(self,*args):
        return GiacMethods['polyEval'](self,*args)

     def polygon(self,*args):
        return GiacMethods['polygon'](self,*args)

     def polygone_rempli(self,*args):
        return GiacMethods['polygone_rempli'](self,*args)

     def polygonplot(self,*args):
        return GiacMethods['polygonplot'](self,*args)

     def polygonscatterplot(self,*args):
        return GiacMethods['polygonscatterplot'](self,*args)

     def polyhedron(self,*args):
        return GiacMethods['polyhedron'](self,*args)

     def polynom(self,*args):
        return GiacMethods['polynom'](self,*args)

     def polynomial_regression(self,*args):
        return GiacMethods['polynomial_regression'](self,*args)

     def polynomial_regression_plot(self,*args):
        return GiacMethods['polynomial_regression_plot'](self,*args)

     def position(self,*args):
        return GiacMethods['position'](self,*args)

     def poslbdLMQ(self,*args):
        return GiacMethods['poslbdLMQ'](self,*args)

     def posubLMQ(self,*args):
        return GiacMethods['posubLMQ'](self,*args)

     def potential(self,*args):
        return GiacMethods['potential'](self,*args)

     def pow2exp(self,*args):
        return GiacMethods['pow2exp'](self,*args)

     def power_regression(self,*args):
        return GiacMethods['power_regression'](self,*args)

     def power_regression_plot(self,*args):
        return GiacMethods['power_regression_plot'](self,*args)

     def powermod(self,*args):
        return GiacMethods['powermod'](self,*args)

     def powerpc(self,*args):
        return GiacMethods['powerpc'](self,*args)

     def powexpand(self,*args):
        return GiacMethods['powexpand'](self,*args)

     def powmod(self,*args):
        return GiacMethods['powmod'](self,*args)

     def prepend(self,*args):
        return GiacMethods['prepend'](self,*args)

     def preval(self,*args):
        return GiacMethods['preval'](self,*args)

     def prevperm(self,*args):
        return GiacMethods['prevperm'](self,*args)

     def prevprime(self,*args):
        return GiacMethods['prevprime'](self,*args)

     def primpart(self,*args):
        return GiacMethods['primpart'](self,*args)

     def prism(self,*args):
        return GiacMethods['prism'](self,*args)

     def product(self,*args):
        return GiacMethods['product'](self,*args)

     def projection(self,*args):
        return GiacMethods['projection'](self,*args)

     def proot(self,*args):
        return GiacMethods['proot'](self,*args)

     def propFrac(self,*args):
        return GiacMethods['propFrac'](self,*args)

     def propfrac(self,*args):
        return GiacMethods['propfrac'](self,*args)

     def psrgcd(self,*args):
        return GiacMethods['psrgcd'](self,*args)

     def ptayl(self,*args):
        return GiacMethods['ptayl'](self,*args)

     def purge(self,*args):
        return GiacMethods['purge'](self,*args)

     def pwd(self,*args):
        return GiacMethods['pwd'](self,*args)

     def pyramid(self,*args):
        return GiacMethods['pyramid'](self,*args)

     def q2a(self,*args):
        return GiacMethods['q2a'](self,*args)

     def qr(self,*args):
        return GiacMethods['qr'](self,*args)

     def quadric(self,*args):
        return GiacMethods['quadric'](self,*args)

     def quadrilateral(self,*args):
        return GiacMethods['quadrilateral'](self,*args)

     def quantile(self,*args):
        return GiacMethods['quantile'](self,*args)

     def quartile1(self,*args):
        return GiacMethods['quartile1'](self,*args)

     def quartile3(self,*args):
        return GiacMethods['quartile3'](self,*args)

     def quartiles(self,*args):
        return GiacMethods['quartiles'](self,*args)

     def quest(self,*args):
        return GiacMethods['quest'](self,*args)

     def quo(self,*args):
        return GiacMethods['quo'](self,*args)

     def quorem(self,*args):
        return GiacMethods['quorem'](self,*args)

     def quote(self,*args):
        return GiacMethods['quote'](self,*args)

     def r2e(self,*args):
        return GiacMethods['r2e'](self,*args)

     def radical_axis(self,*args):
        return GiacMethods['radical_axis'](self,*args)

     def radius(self,*args):
        return GiacMethods['radius'](self,*args)

     def ramene(self,*args):
        return GiacMethods['ramene'](self,*args)

     def rand(self,*args):
        return GiacMethods['rand'](self,*args)

     def randMat(self,*args):
        return GiacMethods['randMat'](self,*args)

     def randNorm(self,*args):
        return GiacMethods['randNorm'](self,*args)

     def randPoly(self,*args):
        return GiacMethods['randPoly'](self,*args)

     def randbinomial(self,*args):
        return GiacMethods['randbinomial'](self,*args)

     def randchisquare(self,*args):
        return GiacMethods['randchisquare'](self,*args)

     def randexp(self,*args):
        return GiacMethods['randexp'](self,*args)

     def randfisher(self,*args):
        return GiacMethods['randfisher'](self,*args)

     def randgeometric(self,*args):
        return GiacMethods['randgeometric'](self,*args)

     def randmarkov(self,*args):
        return GiacMethods['randmarkov'](self,*args)

     def randmatrix(self,*args):
        return GiacMethods['randmatrix'](self,*args)

     def randmultinomial(self,*args):
        return GiacMethods['randmultinomial'](self,*args)

     def randnorm(self,*args):
        return GiacMethods['randnorm'](self,*args)

     def randperm(self,*args):
        return GiacMethods['randperm'](self,*args)

     def randpoisson(self,*args):
        return GiacMethods['randpoisson'](self,*args)

     def randpoly(self,*args):
        return GiacMethods['randpoly'](self,*args)

     def randseed(self,*args):
        return GiacMethods['randseed'](self,*args)

     def randstudent(self,*args):
        return GiacMethods['randstudent'](self,*args)

     def randvector(self,*args):
        return GiacMethods['randvector'](self,*args)

     def rank(self,*args):
        return GiacMethods['rank'](self,*args)

     def ranm(self,*args):
        return GiacMethods['ranm'](self,*args)

     def ranv(self,*args):
        return GiacMethods['ranv'](self,*args)

     def rassembler_trigo(self,*args):
        return GiacMethods['rassembler_trigo'](self,*args)

     def rat_jordan(self,*args):
        return GiacMethods['rat_jordan'](self,*args)

     def rational(self,*args):
        return GiacMethods['rational'](self,*args)

     def rationalroot(self,*args):
        return GiacMethods['rationalroot'](self,*args)

     def ratnormal(self,*args):
        return GiacMethods['ratnormal'](self,*args)

     def rcl(self,*args):
        return GiacMethods['rcl'](self,*args)

     def rdiv(self,*args):
        return GiacMethods['rdiv'](self,*args)

     def re(self,*args):
        return GiacMethods['re'](self,*args)

     def read(self,*args):
        return GiacMethods['read'](self,*args)

     def readrgb(self,*args):
        return GiacMethods['readrgb'](self,*args)

     def readwav(self,*args):
        return GiacMethods['readwav'](self,*args)

     def real(self,*args):
        return GiacMethods['real'](self,*args)

     def realroot(self,*args):
        return GiacMethods['realroot'](self,*args)

     def reciprocation(self,*args):
        return GiacMethods['reciprocation'](self,*args)

     def rectangle(self,*args):
        return GiacMethods['rectangle'](self,*args)

     def rectangle_droit(self,*args):
        return GiacMethods['rectangle_droit'](self,*args)

     def rectangle_gauche(self,*args):
        return GiacMethods['rectangle_gauche'](self,*args)

     def rectangle_plein(self,*args):
        return GiacMethods['rectangle_plein'](self,*args)

     def rectangular_coordinates(self,*args):
        return GiacMethods['rectangular_coordinates'](self,*args)

     def recule(self,*args):
        return GiacMethods['recule'](self,*args)

     def red(self,*args):
        return GiacMethods['red'](self,*args)

     def reduced_conic(self,*args):
        return GiacMethods['reduced_conic'](self,*args)

     def reduced_quadric(self,*args):
        return GiacMethods['reduced_quadric'](self,*args)

     def ref(self,*args):
        return GiacMethods['ref'](self,*args)

     def reflection(self,*args):
        return GiacMethods['reflection'](self,*args)

     def regroup(self,*args):
        return GiacMethods['regroup'](self,*args)

     def rem(self,*args):
        return GiacMethods['rem'](self,*args)

     def remain(self,*args):
        return GiacMethods['remain'](self,*args)

     def remove(self,*args):
        return GiacMethods['remove'](self,*args)

     def reorder(self,*args):
        return GiacMethods['reorder'](self,*args)

     def residue(self,*args):
        return GiacMethods['residue'](self,*args)

     def resoudre(self,*args):
        return GiacMethods['resoudre'](self,*args)

     def resoudre_dans_C(self,*args):
        return GiacMethods['resoudre_dans_C'](self,*args)

     def resoudre_systeme_lineaire(self,*args):
        return GiacMethods['resoudre_systeme_lineaire'](self,*args)

     def resultant(self,*args):
        return GiacMethods['resultant'](self,*args)

     def reverse_rsolve(self,*args):
        return GiacMethods['reverse_rsolve'](self,*args)

     def revert(self,*args):
        return GiacMethods['revert'](self,*args)

     def revlex(self,*args):
        return GiacMethods['revlex'](self,*args)

     def revlist(self,*args):
        return GiacMethods['revlist'](self,*args)

     def rhombus(self,*args):
        return GiacMethods['rhombus'](self,*args)

     def rhombus_point(self,*args):
        return GiacMethods['rhombus_point'](self,*args)

     def rhs(self,*args):
        return GiacMethods['rhs'](self,*args)

     def right(self,*args):
        return GiacMethods['right'](self,*args)

     def right_rectangle(self,*args):
        return GiacMethods['right_rectangle'](self,*args)

     def right_triangle(self,*args):
        return GiacMethods['right_triangle'](self,*args)

     def risch(self,*args):
        return GiacMethods['risch'](self,*args)

     def rm_a_z(self,*args):
        return GiacMethods['rm_a_z'](self,*args)

     def rm_all_vars(self,*args):
        return GiacMethods['rm_all_vars'](self,*args)

     def rmbreakpoint(self,*args):
        return GiacMethods['rmbreakpoint'](self,*args)

     def rmmod(self,*args):
        return GiacMethods['rmmod'](self,*args)

     def rmwatch(self,*args):
        return GiacMethods['rmwatch'](self,*args)

     def romberg(self,*args):
        return GiacMethods['romberg'](self,*args)

     def rombergm(self,*args):
        return GiacMethods['rombergm'](self,*args)

     def rombergt(self,*args):
        return GiacMethods['rombergt'](self,*args)

     def rond(self,*args):
        return GiacMethods['rond'](self,*args)

     def root(self,*args):
        return GiacMethods['root'](self,*args)

     def rootof(self,*args):
        return GiacMethods['rootof'](self,*args)

     def roots(self,*args):
        return GiacMethods['roots'](self,*args)

     def rotate(self,*args):
        return GiacMethods['rotate'](self,*args)

     def rotation(self,*args):
        return GiacMethods['rotation'](self,*args)

     def round(self,*args):
        return GiacMethods['round'](self,*args)

     def row(self,*args):
        return GiacMethods['row'](self,*args)

     def rowAdd(self,*args):
        return GiacMethods['rowAdd'](self,*args)

     def rowDim(self,*args):
        return GiacMethods['rowDim'](self,*args)

     def rowNorm(self,*args):
        return GiacMethods['rowNorm'](self,*args)

     def rowSwap(self,*args):
        return GiacMethods['rowSwap'](self,*args)

     def rowdim(self,*args):
        return GiacMethods['rowdim'](self,*args)

     def rownorm(self,*args):
        return GiacMethods['rownorm'](self,*args)

     def rowspace(self,*args):
        return GiacMethods['rowspace'](self,*args)

     def rowswap(self,*args):
        return GiacMethods['rowswap'](self,*args)

     def rref(self,*args):
        return GiacMethods['rref'](self,*args)

     def rsolve(self,*args):
        return GiacMethods['rsolve'](self,*args)

     def same(self,*args):
        return GiacMethods['same'](self,*args)

     def sans_factoriser(self,*args):
        return GiacMethods['sans_factoriser'](self,*args)

     def saute(self,*args):
        return GiacMethods['saute'](self,*args)

     def scalarProduct(self,*args):
        return GiacMethods['scalarProduct'](self,*args)

     def scalar_product(self,*args):
        return GiacMethods['scalar_product'](self,*args)

     def scatterplot(self,*args):
        return GiacMethods['scatterplot'](self,*args)

     def schur(self,*args):
        return GiacMethods['schur'](self,*args)

     def sec(self,*args):
        return GiacMethods['sec'](self,*args)

     def secant_solver(self,*args):
        return GiacMethods['secant_solver'](self,*args)

     def segment(self,*args):
        return GiacMethods['segment'](self,*args)

     def select(self,*args):
        return GiacMethods['select'](self,*args)

     def semi_augment(self,*args):
        return GiacMethods['semi_augment'](self,*args)

     def seq(self,*args):
        return GiacMethods['seq'](self,*args)

     def seqplot(self,*args):
        return GiacMethods['seqplot'](self,*args)

     def seqsolve(self,*args):
        return GiacMethods['seqsolve'](self,*args)

     def series(self,*args):
        return GiacMethods['series'](self,*args)

     def shift(self,*args):
        return GiacMethods['shift'](self,*args)

     def shift_phase(self,*args):
        return GiacMethods['shift_phase'](self,*args)

     def sign(self,*args):
        return GiacMethods['sign'](self,*args)

     def signature(self,*args):
        return GiacMethods['signature'](self,*args)

     def signe(self,*args):
        return GiacMethods['signe'](self,*args)

     def similarity(self,*args):
        return GiacMethods['similarity'](self,*args)

     def simp2(self,*args):
        return GiacMethods['simp2'](self,*args)

     def simplex_reduce(self,*args):
        return GiacMethods['simplex_reduce'](self,*args)

     def simplifier(self,*args):
        return GiacMethods['simplifier'](self,*args)

     def simplify(self,*args):
        return GiacMethods['simplify'](self,*args)

     def simpson(self,*args):
        return GiacMethods['simpson'](self,*args)

     def simult(self,*args):
        return GiacMethods['simult'](self,*args)

     def sin(self,*args):
        return GiacMethods['sin'](self,*args)

     def sin2costan(self,*args):
        return GiacMethods['sin2costan'](self,*args)

     def sincos(self,*args):
        return GiacMethods['sincos'](self,*args)

     def single_inter(self,*args):
        return GiacMethods['single_inter'](self,*args)

     def sinh(self,*args):
        return GiacMethods['sinh'](self,*args)

     def sizes(self,*args):
        return GiacMethods['sizes'](self,*args)

     def slope(self,*args):
        return GiacMethods['slope'](self,*args)

     def slopeat(self,*args):
        return GiacMethods['slopeat'](self,*args)

     def slopeatraw(self,*args):
        return GiacMethods['slopeatraw'](self,*args)

     def smod(self,*args):
        return GiacMethods['smod'](self,*args)

     def snedecor(self,*args):
        return GiacMethods['snedecor'](self,*args)

     def snedecor_cdf(self,*args):
        return GiacMethods['snedecor_cdf'](self,*args)

     def snedecor_icdf(self,*args):
        return GiacMethods['snedecor_icdf'](self,*args)

     def snedecord(self,*args):
        return GiacMethods['snedecord'](self,*args)

     def snedecord_cdf(self,*args):
        return GiacMethods['snedecord_cdf'](self,*args)

     def snedecord_icdf(self,*args):
        return GiacMethods['snedecord_icdf'](self,*args)

     def solid_line(self,*args):
        return GiacMethods['solid_line'](self,*args)

     def solve(self,*args):
        return GiacMethods['solve'](self,*args)

     def somme(self,*args):
        return GiacMethods['somme'](self,*args)

     def sommet(self,*args):
        return GiacMethods['sommet'](self,*args)

     def sort(self,*args):
        return GiacMethods['sort'](self,*args)

     def sorta(self,*args):
        return GiacMethods['sorta'](self,*args)

     def sortd(self,*args):
        return GiacMethods['sortd'](self,*args)

     def soundsec(self,*args):
        return GiacMethods['soundsec'](self,*args)

     def sphere(self,*args):
        return GiacMethods['sphere'](self,*args)

     def spline(self,*args):
        return GiacMethods['spline'](self,*args)

     def split(self,*args):
        return GiacMethods['split'](self,*args)

     def sq(self,*args):
        return GiacMethods['sq'](self,*args)

     def sqrfree(self,*args):
        return GiacMethods['sqrfree'](self,*args)

     def sqrt(self,*args):
        return GiacMethods['sqrt'](self,*args)

     def square(self,*args):
        return GiacMethods['square'](self,*args)

     def square_point(self,*args):
        return GiacMethods['square_point'](self,*args)

     def srand(self,*args):
        return GiacMethods['srand'](self,*args)

     def sst(self,*args):
        return GiacMethods['sst'](self,*args)

     def sst_in(self,*args):
        return GiacMethods['sst_in'](self,*args)

     def star_point(self,*args):
        return GiacMethods['star_point'](self,*args)

     def start(self,*args):
        return GiacMethods['start'](self,*args)

     def stdDev(self,*args):
        return GiacMethods['stdDev'](self,*args)

     def stddev(self,*args):
        return GiacMethods['stddev'](self,*args)

     def stddevp(self,*args):
        return GiacMethods['stddevp'](self,*args)

     def steffenson_solver(self,*args):
        return GiacMethods['steffenson_solver'](self,*args)

     def student(self,*args):
        return GiacMethods['student'](self,*args)

     def student_cdf(self,*args):
        return GiacMethods['student_cdf'](self,*args)

     def student_icdf(self,*args):
        return GiacMethods['student_icdf'](self,*args)

     def studentd(self,*args):
        return GiacMethods['studentd'](self,*args)

     def studentt(self,*args):
        return GiacMethods['studentt'](self,*args)

     def sturm(self,*args):
        return GiacMethods['sturm'](self,*args)

     def sturmab(self,*args):
        return GiacMethods['sturmab'](self,*args)

     def sturmseq(self,*args):
        return GiacMethods['sturmseq'](self,*args)

     def style(self,*args):
        return GiacMethods['style'](self,*args)

     def subMat(self,*args):
        return GiacMethods['subMat'](self,*args)

     def subs(self,*args):
        return GiacMethods['subs'](self,*args)

     def subsop(self,*args):
        return GiacMethods['subsop'](self,*args)

     def subst(self,*args):
        return GiacMethods['subst'](self,*args)

     def substituer(self,*args):
        return GiacMethods['substituer'](self,*args)

     def subtype(self,*args):
        return GiacMethods['subtype'](self,*args)

     def sum(self,*args):
        return GiacMethods['sum'](self,*args)

     def sum_riemann(self,*args):
        return GiacMethods['sum_riemann'](self,*args)

     def suppress(self,*args):
        return GiacMethods['suppress'](self,*args)

     def surd(self,*args):
        return GiacMethods['surd'](self,*args)

     def svd(self,*args):
        return GiacMethods['svd'](self,*args)

     def swapcol(self,*args):
        return GiacMethods['swapcol'](self,*args)

     def swaprow(self,*args):
        return GiacMethods['swaprow'](self,*args)

     def switch_axes(self,*args):
        return GiacMethods['switch_axes'](self,*args)

     def sylvester(self,*args):
        return GiacMethods['sylvester'](self,*args)

     def symb2poly(self,*args):
        return GiacMethods['symb2poly'](self,*args)

     def symbol(self,*args):
        return GiacMethods['symbol'](self,*args)

     def syst2mat(self,*args):
        return GiacMethods['syst2mat'](self,*args)

     def tCollect(self,*args):
        return GiacMethods['tCollect'](self,*args)

     def tExpand(self,*args):
        return GiacMethods['tExpand'](self,*args)

     def table(self,*args):
        return GiacMethods['table'](self,*args)

     def tablefunc(self,*args):
        return GiacMethods['tablefunc'](self,*args)

     def tableseq(self,*args):
        return GiacMethods['tableseq'](self,*args)

     def tail(self,*args):
        return GiacMethods['tail'](self,*args)

     def tan(self,*args):
        return GiacMethods['tan'](self,*args)

     def tan2cossin2(self,*args):
        return GiacMethods['tan2cossin2'](self,*args)

     def tan2sincos(self,*args):
        return GiacMethods['tan2sincos'](self,*args)

     def tan2sincos2(self,*args):
        return GiacMethods['tan2sincos2'](self,*args)

     def tangent(self,*args):
        return GiacMethods['tangent'](self,*args)

     def tangente(self,*args):
        return GiacMethods['tangente'](self,*args)

     def tanh(self,*args):
        return GiacMethods['tanh'](self,*args)

     def taux_accroissement(self,*args):
        return GiacMethods['taux_accroissement'](self,*args)

     def taylor(self,*args):
        return GiacMethods['taylor'](self,*args)

     def tchebyshev1(self,*args):
        return GiacMethods['tchebyshev1'](self,*args)

     def tchebyshev2(self,*args):
        return GiacMethods['tchebyshev2'](self,*args)

     def tcoeff(self,*args):
        return GiacMethods['tcoeff'](self,*args)

     def tcollect(self,*args):
        return GiacMethods['tcollect'](self,*args)

     def tdeg(self,*args):
        return GiacMethods['tdeg'](self,*args)

     def tetrahedron(self,*args):
        return GiacMethods['tetrahedron'](self,*args)

     def texpand(self,*args):
        return GiacMethods['texpand'](self,*args)

     def thickness(self,*args):
        return GiacMethods['thickness'](self,*args)

     def throw(self,*args):
        return GiacMethods['throw'](self,*args)

     def title(self,*args):
        return GiacMethods['title'](self,*args)

     def titre(self,*args):
        return GiacMethods['titre'](self,*args)

     def tlin(self,*args):
        return GiacMethods['tlin'](self,*args)

     def tourne_droite(self,*args):
        return GiacMethods['tourne_droite'](self,*args)

     def tourne_gauche(self,*args):
        return GiacMethods['tourne_gauche'](self,*args)

     def trace(self,*args):
        return GiacMethods['trace'](self,*args)

     def trames(self,*args):
        return GiacMethods['trames'](self,*args)

     def tran(self,*args):
        return GiacMethods['tran'](self,*args)

     def translation(self,*args):
        return GiacMethods['translation'](self,*args)

     def transpose(self,*args):
        return GiacMethods['transpose'](self,*args)

     def trapeze(self,*args):
        return GiacMethods['trapeze'](self,*args)

     def trapezoid(self,*args):
        return GiacMethods['trapezoid'](self,*args)

     def triangle(self,*args):
        return GiacMethods['triangle'](self,*args)

     def triangle_paper(self,*args):
        return GiacMethods['triangle_paper'](self,*args)

     def triangle_plein(self,*args):
        return GiacMethods['triangle_plein'](self,*args)

     def triangle_point(self,*args):
        return GiacMethods['triangle_point'](self,*args)

     def trig2exp(self,*args):
        return GiacMethods['trig2exp'](self,*args)

     def trigcos(self,*args):
        return GiacMethods['trigcos'](self,*args)

     def trigexpand(self,*args):
        return GiacMethods['trigexpand'](self,*args)

     def trigsin(self,*args):
        return GiacMethods['trigsin'](self,*args)

     def trigtan(self,*args):
        return GiacMethods['trigtan'](self,*args)

     def trn(self,*args):
        return GiacMethods['trn'](self,*args)

     def true(self,*args):
        return GiacMethods['true'](self,*args)

     def trunc(self,*args):
        return GiacMethods['trunc'](self,*args)

     def truncate(self,*args):
        return GiacMethods['truncate'](self,*args)

     def tsimplify(self,*args):
        return GiacMethods['tsimplify'](self,*args)

     def ufactor(self,*args):
        return GiacMethods['ufactor'](self,*args)

     def ugamma(self,*args):
        return GiacMethods['ugamma'](self,*args)

     def unapply(self,*args):
        return GiacMethods['unapply'](self,*args)

     def unarchive(self,*args):
        return GiacMethods['unarchive'](self,*args)

     def unfactored(self,*args):
        return GiacMethods['unfactored'](self,*args)

     def uniform(self,*args):
        return GiacMethods['uniform'](self,*args)

     def uniform_cdf(self,*args):
        return GiacMethods['uniform_cdf'](self,*args)

     def uniform_icdf(self,*args):
        return GiacMethods['uniform_icdf'](self,*args)

     def uniformd(self,*args):
        return GiacMethods['uniformd'](self,*args)

     def uniformd_cdf(self,*args):
        return GiacMethods['uniformd_cdf'](self,*args)

     def uniformd_icdf(self,*args):
        return GiacMethods['uniformd_icdf'](self,*args)

     def unitV(self,*args):
        return GiacMethods['unitV'](self,*args)

     def unquote(self,*args):
        return GiacMethods['unquote'](self,*args)

     def user_operator(self,*args):
        return GiacMethods['user_operator'](self,*args)

     def usimplify(self,*args):
        return GiacMethods['usimplify'](self,*args)

     def valuation(self,*args):
        return GiacMethods['valuation'](self,*args)

     def vandermonde(self,*args):
        return GiacMethods['vandermonde'](self,*args)

     def variables_are_files(self,*args):
        return GiacMethods['variables_are_files'](self,*args)

     def variance(self,*args):
        return GiacMethods['variance'](self,*args)

     def version(self,*args):
        return GiacMethods['version'](self,*args)

     def vertices(self,*args):
        return GiacMethods['vertices'](self,*args)

     def vertices_abc(self,*args):
        return GiacMethods['vertices_abc'](self,*args)

     def vertices_abca(self,*args):
        return GiacMethods['vertices_abca'](self,*args)

     def vpotential(self,*args):
        return GiacMethods['vpotential'](self,*args)

     def weibull(self,*args):
        return GiacMethods['weibull'](self,*args)

     def weibull_cdf(self,*args):
        return GiacMethods['weibull_cdf'](self,*args)

     def weibull_icdf(self,*args):
        return GiacMethods['weibull_icdf'](self,*args)

     def weibulld(self,*args):
        return GiacMethods['weibulld'](self,*args)

     def weibulld_cdf(self,*args):
        return GiacMethods['weibulld_cdf'](self,*args)

     def weibulld_icdf(self,*args):
        return GiacMethods['weibulld_icdf'](self,*args)

     def widget_size(self,*args):
        return GiacMethods['widget_size'](self,*args)

     def wilcoxonp(self,*args):
        return GiacMethods['wilcoxonp'](self,*args)

     def wilcoxons(self,*args):
        return GiacMethods['wilcoxons'](self,*args)

     def wilcoxont(self,*args):
        return GiacMethods['wilcoxont'](self,*args)

     def writergb(self,*args):
        return GiacMethods['writergb'](self,*args)

     def writewav(self,*args):
        return GiacMethods['writewav'](self,*args)

     def xyztrange(self,*args):
        return GiacMethods['xyztrange'](self,*args)

     def zeros(self,*args):
        return GiacMethods['zeros'](self,*args)

     def ztrans(self,*args):
        return GiacMethods['ztrans'](self,*args)

     def type(self,*args):
        return GiacMethods['type'](self,*args)

     def zip(self,*args):
        return GiacMethods['zip'](self,*args)




##
################################################################
#   A wrapper from a cpp element of type giac gen to create    #
#   the Python object                                          #
################################################################
cdef inline _wrap_gen(gen  g)except +:

#   cdef Pygen pyg=Pygen('NULL')
# It is much faster with ''
#      [x-x for i in range(10**4)]
#      ('clock: ', 0.010000000000000009,
# than with 'NULL':
#      [x-x for i in range(10**4)]
#      ('clock: ', 1.1999999999999997,
#    #    #    #    #    #
# This is faster than with:
#    cdef Pygen pyg=Pygen('')
# ll=giac(range(10**6))
# ('clock: ', 0.40000000000000036, ' time: ', 0.40346789360046387)
# gg=[1 for i in ll]
# ('clock: ', 0.669999999999999, ' time: ', 0.650738000869751)
#
# But beware when changing the None case in  Pygen init.
#
    sig_on()
    cdef Pygen pyg=Pygen()
    del pyg.gptr # Pygen.__cinit__() always creates a gen. So we have to delete it here.
    pyg.gptr=new gen(g)
    sig_off()
    return pyg
#    if(pyg.gptr !=NULL):
#      return pyg
#    else:
#      raise MemoryError,"empty gen"



################################################################
#    A wrapper from a python list to a vector of gen           #
################################################################

cdef  vecteur _wrap_pylist(L) except +:
   cdef vecteur  * V
   cdef int i

   if (isinstance(L,tuple) or isinstance(L,listrange)):
      n=len(L)
      V=new vecteur()

      sig_on()
      for i in range(n):
        V.push_back((<Pygen>Pygen(L[i])).gptr[0])
      sig_off()
      return V[0]
   else:
     raise TypeError,"argument must be a tuple or a list"


#################################
#  slice wrapper for a giac list
#################################
cdef  vecteur _getgiacslice(Pygen L,slice sl) except +:
   cdef vecteur  * V
   cdef int u

   if (L.type()=="DOM_LIST"):
      n=len(L)
      V=new vecteur()

      sig_on()
#      for u in range(n)[sl]:   #pb python3
      (b,e,st)=sl.indices(n)
      for u in range(b,e,st):
        V.push_back((L.gptr[0])[u])
      sig_off()
      return V[0]
   else:
     raise TypeError,"argument must be a Pygen list and a slice"





cdef  gen pylongtogen(a) except +:
   #                                                                     #
   # basic conversion of Python long integers to gen via Horner's Method #
   #                                                                     #

   aneg=False
   cdef gen g=gen(<int>0)
   cdef gen M

   if (a<0):
     aneg=True
     a=-a
   if Pyversioninfo >= (2,7):
      size=a.bit_length()  # bit_length python >= 2.7 required.
      shift=Pymaxint.bit_length()-1
   else:
      size=math.trunc(math.log(a,2))+1
      shift=math.trunc(math.log(Pymaxint))
   M=gen(<long long>(1<<shift))

   while (size>=shift):
     size=size-shift
     i=int(a>>size)
     g=(g*M+gen(<long long>i))
     a=a-(i<<size)

   g=g*gen(<long long>(1<<size))+gen(<long long> a)
   if aneg:
     g=-g
   return g;



#############################################################
# Examples of python functions directly implemented from giac
#############################################################
#def giaceval(Pygen self):
#    cdef gen result
#    try:
#      result=GIAC_protecteval(self.gptr[0],1,context_ptr)
#      return _wrap_gen(result)
#    except:
#      raise
#
#
#def giacfactor(Pygen self):
#
#    cdef gen result
#    try:
#      result=GIAC_factor(self.gptr[0],context_ptr)
#      return _wrap_gen(result)
#    except:
#      raise
#
#
#
#def giacfactors(Pygen self):
#    cdef gen result
#    try:
#      result=GIAC_factors(self.gptr[0],context_ptr)
#      return _wrap_gen(result)
#    except:
#      raise
#
#
#
#
#def giacnormal(Pygen self):
#    cdef gen result
#    try:
#      result=GIAC_normal(self.gptr[0],context_ptr)
#      return _wrap_gen(result)
#    except:
#      raise
#
#
#def giacgcd(Pygen a, Pygen b):
#    cdef gen result
#    try:
#      result=gen( GIAC_makenewvecteur(a.gptr[0],b.gptr[0]) ,<short int>1)
#      result=GIAC_gcd(result,context_ptr)
#      return _wrap_gen(result)
#    except:
#      raise





#############################################################
#  Most giac keywords
#############################################################
include 'keywords.pxi'
GiacMethods={}





#############################################################
# Some convenient settings
#############################################################
Pygen('printpow(1)').eval() ; # default power is ^
Pygen('add_language(1)').eval(); # Add the french keywords in the giac library language.
# FIXME: print I for sqrt(-1) instead of i
# GIAC_try_parse_i(False,context_ptr); (does not work??)


for i in mostkeywords:
   tmp=Pygen(i)
   # in the sage version we remove:    globals()[i]=tmp
   GiacMethods[i]=tmp

# We put the giac names that should not be exported to Python in moremethods.
for i in moremethods:
   tmp=Pygen(i)
   GiacMethods[i]=tmp


# To avoid conflicts we export only these few ones.  Most giac keywords will be
# avaible through: libgiac.keywordname
__all__=['Pygen','giacsettings','libgiac','loadgiacgen']


##
def moretests():
   """
     sage: from giacpy import libgiac
     sage: libgiac(3^100)
     515377520732011331036461129765621272702107522001
     sage: libgiac(-3^100)
     -515377520732011331036461129765621272702107522001
     sage: libgiac(-11^1000)
     -2469932918005826334124088385085221477709733385238396234869182951830739390375433175367866116456946191973803561189036523363533798726571008961243792655536655282201820357872673322901148243453211756020067624545609411212063417307681204817377763465511222635167942816318177424600927358163388910854695041070577642045540560963004207926938348086979035423732739933235077042750354729095729602516751896320598857608367865475244863114521391548985943858154775884418927768284663678512441565517194156946312753546771163991252528017732162399536497445066348868438762510366191040118080751580689254476068034620047646422315123643119627205531371694188794408120267120500325775293645416335230014278578281272863450085145349124727476223298887655183167465713337723258182649072572861625150703747030550736347589416285606367521524529665763903537989935510874657420361426804068643262800901916285076966174176854351055183740078763891951775452021781225066361670593917001215032839838911476044840388663443684517735022039957481918726697789827894303408292584258328090724141496484460001

   """
   return


def loadgiacgen(str filename):
        """
          Open a file in giac compressed format to create a Pygen element.

          Use the save method from Pygen elements to create such files.

          In C++ these files can be opened with giac::unarchive and created with
          giac::archive

            sage: from giacpy import *
            sage: g=libgiac.texpand('cos(10*a+5*b)')
            sage: g.save("fichiertest")           #  doctest: +SKIP
            sage: a=loadgiacgen("fichiertest")    #  doctest: +SKIP
            sage: from tempfile import NamedTemporaryFile
            sage: F=NamedTemporaryFile()   # chose a temporary file for a test
            sage: g.savegen(F.name)
            sage: a=loadgiacgen(F.name)
            sage: a.tcollect()
            cos(10*a+5*b)
            sage: F.close()

        """
        cdef gen result
        try:
            sig_on()
            result=GIAC_unarchive( <string>encstring23(filename), context_ptr)
            sig_off()
            return _wrap_gen(result)
        except:
             raise RuntimeError



class GiacInstance:

    def __init__(self):
       self.__dict__.update(GiacMethods)


    def __call__(self,s):
       return _giac(s)


    def _sage_doc_(self):
       return _giac.__doc__


    def eval(self, code, strip=True, **kwds):

        if strip:
             code = code.replace("\n","").strip()
        return self(code)


    __doc__ = _giac.__doc__



libgiac=GiacInstance()

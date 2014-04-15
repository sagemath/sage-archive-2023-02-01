r"""
Solving ordinary differential equations

This file contains functions useful for solving differential equations
which occur commonly in a 1st semester differential equations
course. For another numerical solver see :meth:`ode_solver` function
and optional package Octave.

Commands:

- ``desolve`` - Computes the "general solution" to a 1st or 2nd order
  ODE via Maxima.

- ``desolve_laplace`` - Solves an ODE using laplace transforms via
  Maxima. Initials conditions are optional.

- ``desolve_system`` - Solves any size system of 1st order odes using
  Maxima. Initials conditions are optional.

- ``desolve_rk4`` - Solves numerically IVP for one first order
  equation, returns list of points or plot

- ``desolve_system_rk4`` - Solves numerically IVP for system of first
  order equations, returns list of points

- ``desolve_odeint`` - Solves numerically a system of first-order ordinary
  differential equations using ``odeint`` from scipy.integrate module.

- ``eulers_method`` - Approximate solution to a 1st order DE,
  presented as a table.

- ``eulers_method_2x2`` - Approximate solution to a 1st order system
  of DEs, presented as a table.

- ``eulers_method_2x2_plot`` - Plots the sequence of points obtained
  from Euler's method.

AUTHORS:

- David Joyner (3-2006) - Initial version of functions

- Marshall Hampton (7-2007) - Creation of Python module and testing

- Robert Bradshaw (10-2008) - Some interface cleanup.

- Robert Marik (10-2009) - Some bugfixes and enhancements

"""

##########################################################################
#  Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>, Marshall Hampton,
#  Robert Marik <marik@mendelu.cz>
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  http://www.gnu.org/licenses/
##########################################################################

from sage.interfaces.maxima import Maxima
from sage.plot.all import line
from sage.symbolic.expression import is_SymbolicEquation
from sage.symbolic.ring import is_SymbolicVariable
from sage.calculus.functional import diff
from sage.misc.decorators import rename_keyword

maxima = Maxima()

def desolve(de, dvar, ics=None, ivar=None, show_method=False, contrib_ode=False):
    r"""
    Solves a 1st or 2nd order linear ODE via maxima. Including IVP and BVP.

    *Use* ``desolve? <tab>`` *if the output in truncated in notebook.*

    INPUT:

    - ``de`` - an expression or equation representing the ODE

    - ``dvar`` - the dependent variable (hereafter called ``y``)

    - ``ics`` - (optional) the initial or boundary conditions

      - for a first-order equation, specify the initial ``x`` and ``y``

      - for a second-order equation, specify the initial ``x``, ``y``,
        and ``dy/dx``, i.e. write `[x_0, y(x_0), y'(x_0)]`

      - for a second-order boundary solution, specify initial and
        final ``x`` and ``y`` boundary conditions, i.e. write `[x_0, y(x_0), x_1, y(x_1)]`.

      - gives an error if the solution is not SymbolicEquation (as happens for
        example for Clairaut equation)

    - ``ivar`` - (optional) the independent variable (hereafter called
      x), which must be specified if there is more than one
      independent variable in the equation.

    - ``show_method`` - (optional) if true, then Sage returns pair
      ``[solution, method]``, where method is the string describing
      method which has been used to get solution (Maxima uses the
      following order for first order equations: linear, separable,
      exact (including exact with integrating factor), homogeneous,
      bernoulli, generalized homogeneous) - use carefully in class,
      see below for the example of the equation which is separable but
      this property is not recognized by Maxima and equation is solved
      as exact.

    - ``contrib_ode`` - (optional) if true, desolve allows to solve
      clairaut, lagrange, riccati and some other equations. May take
      a long time and thus turned off by default.  Initial conditions
      can be used only if the result is one SymbolicEquation (does not
      contain singular solution, for example)

    OUTPUT:

    In most cases returns SymbolicEquation which defines the solution
    implicitly.  If the result is in the form y(x)=... (happens for
    linear eqs.), returns the right-hand side only.  The possible
    constant solutions of separable ODE's are omitted.


    EXAMPLES::

        sage: x = var('x')
        sage: y = function('y', x)
        sage: desolve(diff(y,x) + y - 1, y)
        (c + e^x)*e^(-x)

    ::

        sage: f = desolve(diff(y,x) + y - 1, y, ics=[10,2]); f
        (e^10 + e^x)*e^(-x)

    ::

        sage: plot(f)

    We can also solve second-order differential equations.::

        sage: x = var('x')
        sage: y = function('y', x)
        sage: de = diff(y,x,2) - y == x
        sage: desolve(de, y)
        k2*e^(-x) + k1*e^x - x


    ::

        sage: f = desolve(de, y, [10,2,1]); f
        -x + 7*e^(x - 10) + 5*e^(-x + 10)

    ::

        sage: f(x=10)
        2

    ::

        sage: diff(f,x)(x=10)
        1

    ::

        sage: de = diff(y,x,2) + y == 0
        sage: desolve(de, y)
        k2*cos(x) + k1*sin(x)

    ::

        sage: desolve(de, y, [0,1,pi/2,4])
        cos(x) + 4*sin(x)

    ::

        sage: desolve(y*diff(y,x)+sin(x)==0,y)
        -1/2*y(x)^2 == c - cos(x)

    Clairot equation: general and singular solutions::

        sage: desolve(diff(y,x)^2+x*diff(y,x)-y==0,y,contrib_ode=True,show_method=True)
        [[y(x) == c^2 + c*x, y(x) == -1/4*x^2], 'clairault']

    For equations involving more variables we specify independent variable::

        sage: a,b,c,n=var('a b c n')
        sage: desolve(x^2*diff(y,x)==a+b*x^n+c*x^2*y^2,y,ivar=x,contrib_ode=True)
        [[y(x) == 0, (b*x^(n - 2) + a/x^2)*c^2*u == 0]]

    ::

        sage: desolve(x^2*diff(y,x)==a+b*x^n+c*x^2*y^2,y,ivar=x,contrib_ode=True,show_method=True)
        [[[y(x) == 0, (b*x^(n - 2) + a/x^2)*c^2*u == 0]], 'riccati']


    Higher orded, not involving independent variable::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y).expand()
        1/6*y(x)^3 + k1*y(x) == k2 + x

    ::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,1,3]).expand()
        1/6*y(x)^3 - 5/3*y(x) == x - 3/2

    ::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,1,3],show_method=True)
        [1/6*y(x)^3 - 5/3*y(x) == x - 3/2, 'freeofx']

    Separable equations - Sage returns solution in implicit form::

        sage: desolve(diff(y,x)*sin(y) == cos(x),y)
        -cos(y(x)) == c + sin(x)

    ::

        sage: desolve(diff(y,x)*sin(y) == cos(x),y,show_method=True)
        [-cos(y(x)) == c + sin(x), 'separable']

    ::

        sage: desolve(diff(y,x)*sin(y) == cos(x),y,[pi/2,1])
        -cos(y(x)) == -cos(1) + sin(x) - 1

    Linear equation - Sage returns the expression on the right hand side only::

        sage: desolve(diff(y,x)+(y) == cos(x),y)
        1/2*((cos(x) + sin(x))*e^x + 2*c)*e^(-x)

    ::

        sage: desolve(diff(y,x)+(y) == cos(x),y,show_method=True)
        [1/2*((cos(x) + sin(x))*e^x + 2*c)*e^(-x), 'linear']

    ::

        sage: desolve(diff(y,x)+(y) == cos(x),y,[0,1])
        1/2*(cos(x)*e^x + e^x*sin(x) + 1)*e^(-x)

    This ODE with separated variables is solved as
    exact. Explanation - factor does not split `e^{x-y}` in Maxima
    into `e^{x}e^{y}`::

        sage: desolve(diff(y,x)==exp(x-y),y,show_method=True)
        [-e^x + e^y(x) == c, 'exact']

    You can solve Bessel equations. You can also use initial
    conditions, but you cannot put (sometimes desired) initial
    condition at x=0, since this point is singlar point of the
    equation. Anyway, if the solution should be bounded at x=0, then
    k2=0.::

        sage: desolve(x^2*diff(y,x,x)+x*diff(y,x)+(x^2-4)*y==0,y)
        k1*bessel_J(2, x) + k2*bessel_Y(2, x)

    Difficult ODE produces error::

        sage: desolve(sqrt(y)*diff(y,x)+e^(y)+cos(x)-sin(x+y)==0,y) # not tested
        Traceback (click to the left for traceback)
        ...
        NotImplementedError, "Maxima was unable to solve this ODE. Consider to set option contrib_ode to True."

    Difficult ODE produces error - moreover, takes a long time ::

        sage: desolve(sqrt(y)*diff(y,x)+e^(y)+cos(x)-sin(x+y)==0,y,contrib_ode=True) # not tested

    Some more types od ODE's::

        sage: desolve(x*diff(y,x)^2-(1+x*y)*diff(y,x)+y==0,y,contrib_ode=True,show_method=True)
        [[y(x) == c*e^x, y(x) == c + log(x)], 'factor']

    ::

        sage: desolve(diff(y,x)==(x+y)^2,y,contrib_ode=True,show_method=True)
        [[[x == c - arctan(sqrt(t)), y(x) == -x - sqrt(t)], [x == c + arctan(sqrt(t)), y(x) == -x + sqrt(t)]], 'lagrange']

    These two examples produce error (as expected, Maxima 5.18 cannot
    solve equations from initial conditions). Current Maxima 5.18
    returns false answer in this case!::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,2]).expand() # not tested
        Traceback (click to the left for traceback)
        ...
        NotImplementedError, "Maxima was unable to solve this ODE. Consider to set option contrib_ode to True."

    ::

        sage: desolve(diff(y,x,2)+y*(diff(y,x,1))^3==0,y,[0,1,2],show_method=True) # not tested
        Traceback (click to the left for traceback)
        ...
        NotImplementedError, "Maxima was unable to solve this ODE. Consider to set option contrib_ode to True."

    Second order linear ODE::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y)
        (k2*x + k1)*e^(-x) + 1/2*sin(x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,show_method=True)
        [(k2*x + k1)*e^(-x) + 1/2*sin(x), 'variationofparameters']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,1])
        1/2*(7*x + 6)*e^(-x) + 1/2*sin(x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,1],show_method=True)
        [1/2*(7*x + 6)*e^(-x) + 1/2*sin(x), 'variationofparameters']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,pi/2,2])
        3*(x*(e^(1/2*pi) - 2)/pi + 1)*e^(-x) + 1/2*sin(x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == cos(x),y,[0,3,pi/2,2],show_method=True)
        [3*(x*(e^(1/2*pi) - 2)/pi + 1)*e^(-x) + 1/2*sin(x), 'variationofparameters']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y)
        (k2*x + k1)*e^(-x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,show_method=True)
        [(k2*x + k1)*e^(-x), 'constcoeff']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,1])
        (4*x + 3)*e^(-x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,1],show_method=True)
        [(4*x + 3)*e^(-x), 'constcoeff']

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,pi/2,2])
        (2*x*(2*e^(1/2*pi) - 3)/pi + 3)*e^(-x)

    ::

        sage: desolve(diff(y,x,2)+2*diff(y,x)+y == 0,y,[0,3,pi/2,2],show_method=True)
        [(2*x*(2*e^(1/2*pi) - 3)/pi + 3)*e^(-x), 'constcoeff']

    TESTS:

    Trac #9961 fixed (allow assumptions on the dependent variable in desolve)::

        sage: y=function('y',x); assume(x>0); assume(y>0)
        sage: sage.calculus.calculus.maxima('domain:real')  # needed since Maxima 5.26.0 to get the answer as below
        real
        sage: desolve(x*diff(y,x)-x*sqrt(y^2+x^2)-y == 0, y, contrib_ode=True)
        [x - arcsinh(y(x)/x) == c]

    Trac #10682 updated Maxima to 5.26, and it started to show a different
    solution in the complex domain for the ODE above::

        sage: sage.calculus.calculus.maxima('domain:complex')  # back to the default complex domain
        complex
        sage: desolve(x*diff(y,x)-x*sqrt(y^2+x^2)-y == 0, y, contrib_ode=True)
        [1/2*(2*x^2*sqrt(x^(-2)) - 2*x*sqrt(x^(-2))*arcsinh(y(x)/sqrt(x^2)) -
            2*x*sqrt(x^(-2))*arcsinh(y(x)^2/(x*sqrt(y(x)^2))) +
            log(4*(2*x^2*sqrt((x^2*y(x)^2 + y(x)^4)/x^2)*sqrt(x^(-2)) + x^2 +
            2*y(x)^2)/x^2))/(x*sqrt(x^(-2))) == c]

    Trac #6479 fixed::

        sage: x = var('x')
        sage: y = function('y', x)
        sage: desolve( diff(y,x,x) == 0, y, [0,0,1])
        x

    ::

        sage: desolve( diff(y,x,x) == 0, y, [0,1,1])
        x + 1

    Trac #9835 fixed::

        sage: x = var('x')
        sage: y = function('y', x)
        sage: desolve(diff(y,x,2)+y*(1-y^2)==0,y,[0,-1,1,1])
        Traceback (most recent call last):
        ...
        NotImplementedError: Unable to use initial condition for this equation (freeofx).

    Trac #8931 fixed::

        sage: x=var('x'); f=function('f',x); k=var('k'); assume(k>0)
        sage: desolve(diff(f,x,2)/f==k,f,ivar=x)
        k1*e^(sqrt(k)*x) + k2*e^(-sqrt(k)*x)


    AUTHORS:

    - David Joyner (1-2006)

    - Robert Bradshaw (10-2008)

    - Robert Marik (10-2009)

    """
    if is_SymbolicEquation(de):
        de = de.lhs() - de.rhs()
    if is_SymbolicVariable(dvar):
        raise ValueError("You have to declare dependent variable as a function, eg. y=function('y',x)")
    # for backwards compatibility
    if isinstance(dvar, list):
        dvar, ivar = dvar
    elif ivar is None:
        ivars = de.variables()
        ivars = [t for t in ivars if t is not dvar]
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = ivars[0]
    de00 = de._maxima_()
    P = de00.parent()
    dvar_str=P(dvar.operator()).str()
    ivar_str=P(ivar).str()
    de00 = de00.str()
    def sanitize_var(exprs):
        t = exprs.replace("'"+dvar_str+"(_SAGE_VAR_"+ivar_str+")",dvar_str)
        return t.replace("'"+dvar_str+"("+ivar_str+")",dvar_str)
    de0 = sanitize_var(de00)
    ode_solver="ode2"
    cmd="(TEMP:%s(%s,%s,%s), if TEMP=false then TEMP else substitute(%s=%s(%s),TEMP))"%(ode_solver,de0,dvar_str,ivar_str,dvar_str,dvar_str,ivar_str)
    # we produce string like this
    # ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y(x),x)
    soln = P(cmd)

    if str(soln).strip() == 'false':
        if contrib_ode:
            ode_solver="contrib_ode"
            P("load('contrib_ode)")
            cmd="(TEMP:%s(%s,%s,%s), if TEMP=false then TEMP else substitute(%s=%s(%s),TEMP))"%(ode_solver,de0,dvar_str,ivar_str,dvar_str,dvar_str,ivar_str)
            # we produce string like this
            # (TEMP:contrib_ode(x*('diff(y,x,1))^2-(x*y+1)*'diff(y,x,1)+y,y,x), if TEMP=false then TEMP else substitute(y=y(x),TEMP))
            soln = P(cmd)
            if str(soln).strip() == 'false':
                raise NotImplementedError("Maxima was unable to solve this ODE.")
        else:
            raise NotImplementedError("Maxima was unable to solve this ODE. Consider to set option contrib_ode to True.")

    if show_method:
        maxima_method=P("method")

    if (ics is not None):
        if not is_SymbolicEquation(soln.sage()):
            if not show_method:
                maxima_method=P("method")
            raise NotImplementedError("Unable to use initial condition for this equation (%s)."%(str(maxima_method).strip()))
        if len(ics) == 2:
            tempic=(ivar==ics[0])._maxima_().str()
            tempic=tempic+","+(dvar==ics[1])._maxima_().str()
            cmd="(TEMP:ic1(%s(%s,%s,%s),%s),substitute(%s=%s(%s),TEMP))"%(ode_solver,de00,dvar_str,ivar_str,tempic,dvar_str,dvar_str,ivar_str)
            cmd=sanitize_var(cmd)
            # we produce string like this
            # (TEMP:ic2(ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y,x),x=0,y=3,'diff(y,x)=1),substitute(y=y(x),TEMP))
            soln=P(cmd)
        if len(ics) == 3:
            #fixed ic2 command from Maxima - we have to ensure that %k1, %k2 do not depend on variables, should be removed when fixed in Maxima
            P("ic2_sage(soln,xa,ya,dya):=block([programmode:true,backsubst:true,singsolve:true,temp,%k2,%k1,TEMP_k], \
                noteqn(xa), noteqn(ya), noteqn(dya), boundtest('%k1,%k1), boundtest('%k2,%k2), \
                temp: lhs(soln) - rhs(soln), \
                TEMP_k:solve([subst([xa,ya],soln), subst([dya,xa], lhs(dya)=-subst(0,lhs(dya),diff(temp,lhs(xa)))/diff(temp,lhs(ya)))],[%k1,%k2]), \
                if not freeof(lhs(ya),TEMP_k) or not freeof(lhs(xa),TEMP_k) then return (false), \
                temp: maplist(lambda([zz], subst(zz,soln)), TEMP_k), \
                if length(temp)=1 then return(first(temp)) else return(temp))")
            tempic=P(ivar==ics[0]).str()
            tempic=tempic+","+P(dvar==ics[1]).str()
            tempic=tempic+",'diff("+dvar_str+","+ivar_str+")="+P(ics[2]).str()
            cmd="(TEMP:ic2_sage(%s(%s,%s,%s),%s),substitute(%s=%s(%s),TEMP))"%(ode_solver,de00,dvar_str,ivar_str,tempic,dvar_str,dvar_str,ivar_str)
            cmd=sanitize_var(cmd)
            # we produce string like this
            # (TEMP:ic2(ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y,x),x=0,y=3,'diff(y,x)=1),substitute(y=y(x),TEMP))
            soln=P(cmd)
            if str(soln).strip() == 'false':
                raise NotImplementedError("Maxima was unable to solve this IVP. Remove the initial condition to get the general solution.")
        if len(ics) == 4:
            #fixed bc2 command from Maxima - we have to ensure that %k1, %k2 do not depend on variables, should be removed when fixed in Maxima
            P("bc2_sage(soln,xa,ya,xb,yb):=block([programmode:true,backsubst:true,singsolve:true,temp,%k1,%k2,TEMP_k], \
                noteqn(xa), noteqn(ya), noteqn(xb), noteqn(yb), boundtest('%k1,%k1), boundtest('%k2,%k2), \
                TEMP_k:solve([subst([xa,ya],soln), subst([xb,yb],soln)], [%k1,%k2]), \
                if not freeof(lhs(ya),TEMP_k) or not freeof(lhs(xa),TEMP_k) then return (false), \
                temp: maplist(lambda([zz], subst(zz,soln)),TEMP_k), \
                if length(temp)=1 then return(first(temp)) else return(temp))")
            cmd="bc2_sage(%s(%s,%s,%s),%s,%s=%s,%s,%s=%s)"%(ode_solver,de00,dvar_str,ivar_str,P(ivar==ics[0]).str(),dvar_str,P(ics[1]).str(),P(ivar==ics[2]).str(),dvar_str,P(ics[3]).str())
            cmd="(TEMP:%s,substitute(%s=%s(%s),TEMP))"%(cmd,dvar_str,dvar_str,ivar_str)
            cmd=sanitize_var(cmd)
            # we produce string like this
            # (TEMP:bc2(ode2('diff(y,x,2)+2*'diff(y,x,1)+y-cos(x),y,x),x=0,y=3,x=%pi/2,y=2),substitute(y=y(x),TEMP))
            soln=P(cmd)
            if str(soln).strip() == 'false':
                raise NotImplementedError("Maxima was unable to solve this BVP. Remove the initial condition to get the general solution.")

    soln=soln.sage()
    if is_SymbolicEquation(soln) and soln.lhs() == dvar:
        # Remark: Here we do not check that the right hand side does not depend on dvar.
        # This probably will not hapen for soutions obtained via ode2, anyway.
        soln = soln.rhs()
    if show_method:
        return [soln,maxima_method.str()]
    else:
        return soln


#def desolve_laplace2(de,vars,ics=None):
##     """
##     Solves an ODE using laplace transforms via maxima. Initial conditions
##     are optional.

##     INPUT:
##         de    -- a lambda expression representing the ODE
##                  (eg, de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")
##         vars  -- a list of strings representing the variables
##                  (eg, vars = ["x","f"], if x is the independent
##                   variable and f is the dependent variable)
##         ics   -- a list of numbers representing initial conditions,
##                  with symbols allowed which are represented by strings
##                  (eg, f(0)=1, f'(0)=2 is ics = [0,1,2])

##     EXAMPLES:
##         sage: from sage.calculus.desolvers import desolve_laplace
##         sage: x = var('x')
##         sage: f = function('f', x)
##         sage: de = lambda y: diff(y,x,x) - 2*diff(y,x) + y
##         sage: desolve_laplace(de(f(x)),[f,x])
##          #x*%e^x*(?%at('diff('f(x),x,1),x=0))-'f(0)*x*%e^x+'f(0)*%e^x
##         sage: desolve_laplace(de(f(x)),[f,x],[0,1,2])  ## IC option does not work
##          #x*%e^x*(?%at('diff('f(x),x,1),x=0))-'f(0)*x*%e^x+'f(0)*%e^x

##     AUTHOR: David Joyner (1st version 1-2006, 8-2007)
##     """
#    ######## this method seems reasonable but doesn't work for some reason
#    name0 = vars[0]._repr_()[0:(len(vars[0]._repr_())-2-len(str(vars[1])))]
#    name1 = str(vars[1])
#    #maxima("de:"+de+";")
#    if ics!=None:
#        ic0 = maxima("ic:"+str(vars[1])+"="+str(ics[0]))
#        d = len(ics)
#        for i in range(d-1):
#            maxima(vars[0](vars[1])).diff(vars[1],i).atvalue(ic0,ics[i+1])
#    de0 = de._maxima_()
#    #cmd = "desolve("+de+","+vars[1]+"("+vars[0]+"));"
#    #return maxima.eval(cmd)
#    return de0.desolve(vars[0]).rhs()


def desolve_laplace(de, dvar, ics=None, ivar=None):
    """
    Solves an ODE using laplace transforms. Initials conditions are optional.

    INPUT:

    - ``de`` - a lambda expression representing the ODE (eg, de =
      diff(y,x,2) == diff(y,x)+sin(x))

    - ``dvar`` - the dependent variable (eg y)

    - ``ivar`` - (optional) the independent variable (hereafter called
      x), which must be specified if there is more than one
      independent variable in the equation.

    - ``ics`` - a list of numbers representing initial conditions, (eg,
      f(0)=1, f'(0)=2 is ics = [0,1,2])

    OUTPUT:

    Solution of the ODE as symbolic expression

    EXAMPLES::

        sage: u=function('u',x)
        sage: eq = diff(u,x) - exp(-x) - u == 0
        sage: desolve_laplace(eq,u)
        1/2*(2*u(0) + 1)*e^x - 1/2*e^(-x)

    We can use initial conditions::

        sage: desolve_laplace(eq,u,ics=[0,3])
        -1/2*e^(-x) + 7/2*e^x

    The initial conditions do not persist in the system (as they persisted
    in previous versions)::

        sage: desolve_laplace(eq,u)
        1/2*(2*u(0) + 1)*e^x - 1/2*e^(-x)

    ::

        sage: f=function('f', x)
        sage: eq = diff(f,x) + f == 0
        sage: desolve_laplace(eq,f,[0,1])
        e^(-x)

    ::

        sage: x = var('x')
        sage: f = function('f', x)
        sage: de = diff(f,x,x) - 2*diff(f,x) + f
        sage: desolve_laplace(de,f)
        -x*e^x*f(0) + x*e^x*D[0](f)(0) + e^x*f(0)

    ::

        sage: desolve_laplace(de,f,ics=[0,1,2])
        x*e^x + e^x

    TESTS:

    Trac #4839 fixed::

        sage: t=var('t')
        sage: x=function('x', t)
        sage: soln=desolve_laplace(diff(x,t)+x==1, x, ics=[0,2])
        sage: soln
        e^(-t) + 1

    ::

        sage: soln(t=3)
        e^(-3) + 1

    AUTHORS:

    - David Joyner (1-2006,8-2007)

    - Robert Marik (10-2009)
    """
    #This is the original code from David Joyner (inputs and outputs strings)
    #maxima("de:"+de._repr_()+"=0;")
    #if ics!=None:
    #    d = len(ics)
    #    for i in range(0,d-1):
    #        ic = "atvalue(diff("+vars[1]+"("+vars[0]+"),"+str(vars[0])+","+str(i)+"),"+str(vars[0])+"="+str(ics[0])+","+str(ics[1+i])+")"
    #        maxima(ic)
    #
    #cmd = "desolve("+de._repr_()+","+vars[1]+"("+vars[0]+"));"
    #return maxima(cmd).rhs()._maxima_init_()

    ## verbatim copy from desolve - begin
    if is_SymbolicEquation(de):
        de = de.lhs() - de.rhs()
    if is_SymbolicVariable(dvar):
        raise ValueError("You have to declare dependent variable as a function, eg. y=function('y',x)")
    # for backwards compatibility
    if isinstance(dvar, list):
        dvar, ivar = dvar
    elif ivar is None:
        ivars = de.variables()
        ivars = [t for t in ivars if t != dvar]
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = ivars[0]
    ## verbatim copy from desolve - end

    dvar_str = str(dvar)
    def sanitize_var(exprs):  # 'y(x) -> y(x)
        t = exprs.replace("'"+dvar_str,dvar_str)
        return t.replace("'_SAGE_VAR_"+dvar_str,dvar_str)
    de0=de._maxima_()
    P = de0.parent()
    i = dvar_str.find('(')
    dvar_str = dvar_str[:i+1] + '_SAGE_VAR_' + dvar_str[i+1:]
    cmd = sanitize_var("desolve("+de0.str()+","+dvar_str+")")
    soln=P(cmd).rhs()
    if str(soln).strip() == 'false':
        raise NotImplementedError("Maxima was unable to solve this ODE.")
    soln=soln.sage()
    if ics!=None:
        d = len(ics)
        for i in range(0,d-1):
            soln=eval('soln.substitute(diff(dvar,ivar,i)('+str(ivar)+'=ics[0])==ics[i+1])')
    return soln


def desolve_system(des, vars, ics=None, ivar=None):
    """
    Solves any size system of 1st order ODE's. Initials conditions are optional.

    Onedimensional systems are passed to :meth:`desolve_laplace`.

    INPUT:

    - ``des`` - list of ODEs

    - ``vars`` - list of dependent variables

    - ``ics`` - (optional) list of initial values for ivar and vars

    - ``ivar`` - (optional) the independent variable, which must be
      specified if there is more than one independent variable in the
      equation.

    EXAMPLES::

        sage: t = var('t')
        sage: x = function('x', t)
        sage: y = function('y', t)
        sage: de1 = diff(x,t) + y - 1 == 0
        sage: de2 = diff(y,t) - x + 1 == 0
        sage: desolve_system([de1, de2], [x,y])
        [x(t) == (x(0) - 1)*cos(t) - (y(0) - 1)*sin(t) + 1,
         y(t) == (y(0) - 1)*cos(t) + (x(0) - 1)*sin(t) + 1]

    Now we give some initial conditions::

        sage: sol = desolve_system([de1, de2], [x,y], ics=[0,1,2]); sol
        [x(t) == -sin(t) + 1, y(t) == cos(t) + 1]

    ::

        sage: solnx, solny = sol[0].rhs(), sol[1].rhs()
        sage: plot([solnx,solny],(0,1))  # not tested
        sage: parametric_plot((solnx,solny),(0,1))  # not tested

    TESTS:

    Trac #9823 fixed::

        sage: t = var('t')
        sage: x = function('x', t)
        sage: de1 = diff(x,t) + 1 == 0
        sage: desolve_system([de1], [x])
        -t + x(0)

    AUTHORS:

    - Robert Bradshaw (10-2008)
    """
    if len(des)==1:
        return desolve_laplace(des[0], vars[0], ics=ics, ivar=ivar)
    ivars = set([])
    for i, de in enumerate(des):
        if not is_SymbolicEquation(de):
            des[i] = de == 0
        ivars = ivars.union(set(de.variables()))
    if ivar is None:
        ivars = ivars - set(vars)
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = list(ivars)[0]
    dvars = [v._maxima_() for v in vars]
    if ics is not None:
        ivar_ic = ics[0]
        for dvar, ic in zip(dvars, ics[1:]):
            dvar.atvalue(ivar==ivar_ic, ic)
    soln = dvars[0].parent().desolve(des, dvars)
    if str(soln).strip() == 'false':
        raise NotImplementedError("Maxima was unable to solve this system.")
    soln = list(soln)
    for i, sol in enumerate(soln):
        soln[i] = sol.sage()
    if ics is not None:
        ivar_ic = ics[0]
        for dvar, ic in zip(dvars, ics[:1]):
            dvar.atvalue(ivar==ivar_ic, dvar)
    return soln


def desolve_system_strings(des,vars,ics=None):
    r"""
    Solves any size system of 1st order ODE's. Initials conditions are optional.

    This function is obsolete, use desolve_system.

    INPUT:

    - ``de`` - a list of strings representing the ODEs in maxima
      notation (eg, de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")

    - ``vars`` - a list of strings representing the variables (eg,
      vars = ["s","x","y"], where s is the independent variable and
      x,y the dependent variables)

    - ``ics`` - a list of numbers representing initial conditions
      (eg, x(0)=1, y(0)=2 is ics = [0,1,2])

    WARNING:

        The given ics sets the initial values of the dependent vars in
        maxima, so subsequent ODEs involving these variables will have
        these initial conditions automatically imposed.

    EXAMPLES::

        sage: from sage.calculus.desolvers import desolve_system_strings
        sage: s = var('s')
        sage: function('x', s)
        x(s)

    ::

        sage: function('y', s)
        y(s)

    ::

        sage: de1 = lambda z: diff(z[0],s) + z[1] - 1
        sage: de2 = lambda z: diff(z[1],s) - z[0] + 1
        sage: des = [de1([x(s),y(s)]),de2([x(s),y(s)])]
        sage: vars = ["s","x","y"]
        sage: desolve_system_strings(des,vars)
        ["(1-'y(0))*sin(s)+('x(0)-1)*cos(s)+1", "('x(0)-1)*sin(s)+('y(0)-1)*cos(s)+1"]

    ::

        sage: ics = [0,1,-1]
        sage: soln = desolve_system_strings(des,vars,ics); soln
        ['2*sin(s)+1', '1-2*cos(s)']

    ::

        sage: solnx, solny = map(SR, soln)
        sage: RR(solnx(s=3))
        1.28224001611973

    ::

        sage: P1 = plot([solnx,solny],(0,1))
        sage: P2 = parametric_plot((solnx,solny),(0,1))

    Now type show(P1), show(P2) to view these.


    AUTHORS:

    - David Joyner (3-2006, 8-2007)
    """
    d = len(des)
    dess = [de._maxima_init_() + "=0" for de in des]
    for i in range(d):
        cmd="de:" + dess[int(i)] + ";"
        maxima.eval(cmd)
    desstr = "[" + ",".join(dess) + "]"
    d = len(vars)
    varss = list("'" + vars[i] + "(_SAGE_VAR_" + vars[0] + ")" for i in range(1,d))
    varstr = "[" + ",".join(varss) + "]"
    if ics is not None:
        #d = len(ics) ## must be same as len(des)
        for i in range(1,d):
            ic = "atvalue('" + vars[i] + "(_SAGE_VAR_"+vars[0] + ")," + "_SAGE_VAR_"\
             + str(vars[0]) + "=" + str(ics[0]) + "," + str(ics[i]) + ")"
            maxima.eval(ic)
    cmd = "desolve(" + desstr + "," + varstr + ");"
    soln = maxima(cmd)
    return [f.rhs()._maxima_init_().replace("_SAGE_VAR_"+vars[0],vars[0]) for f in soln]

@rename_keyword(deprecation=6094, method="algorithm")
def eulers_method(f,x0,y0,h,x1,algorithm="table"):
    r"""
    This implements Euler's method for finding numerically the
    solution of the 1st order ODE ``y' = f(x,y)``, ``y(a)=c``. The "x"
    column of the table increments from ``x0`` to ``x1`` by ``h`` (so
    ``(x1-x0)/h`` must be an integer). In the "y" column, the new
    y-value equals the old y-value plus the corresponding entry in the
    last column.

    *For pedagogical purposes only.*

    EXAMPLES::

        sage: from sage.calculus.desolvers import eulers_method
        sage: x,y = PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1)
             x                    y                  h*f(x,y)
             0                    1                   -2
           1/2                   -1                 -7/4
             1                -11/4                -11/8

    ::

        sage: x,y = PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1,algorithm="none")
        [[0, 1], [1/2, -1], [1, -11/4], [3/2, -33/8]]

    ::

        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y = PolynomialRing(RR,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1,algorithm="None")
        [[0, 1], [1/2, -1.0], [1, -2.7], [3/2, -4.0]]

    ::

        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y=PolynomialRing(RR,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1)
             x                    y                  h*f(x,y)
             0                    1                 -2.0
           1/2                 -1.0                 -1.7
             1                 -2.7                 -1.3

    ::

        sage: x,y=PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,1,1,1/3,2)
                 x                    y                  h*f(x,y)
                 1                    1                  1/3
               4/3                  4/3                    1
               5/3                  7/3                 17/9
                 2                 38/9                83/27

    ::

        sage: eulers_method(5*x+y-5,0,1,1/2,1,algorithm="none")
        [[0, 1], [1/2, -1], [1, -11/4], [3/2, -33/8]]

    ::

        sage: pts = eulers_method(5*x+y-5,0,1,1/2,1,algorithm="none")
        sage: P1 = list_plot(pts)
        sage: P2 = line(pts)
        sage: (P1+P2).show()

    AUTHORS:

    - David Joyner
    """
    if algorithm=="table":
        print("%10s %20s %25s"%("x","y","h*f(x,y)"))
    n=int((1.0)*(x1-x0)/h)
    x00=x0; y00=y0
    soln = [[x00,y00]]
    for i in range(n+1):
        if algorithm=="table":
            print("%10r %20r %20r"%(x00,y00,h*f(x00,y00)))
        y00 = y00+h*f(x00,y00)
        x00=x00+h
        soln.append([x00,y00])
    if algorithm!="table":
        return soln

@rename_keyword(deprecation=6094, method="algorithm")
def eulers_method_2x2(f,g, t0, x0, y0, h, t1,algorithm="table"):
    r"""
    This implements Euler's method for finding numerically the
    solution of the 1st order system of two ODEs

    ``x' = f(t, x, y), x(t0)=x0.``

    ``y' = g(t, x, y), y(t0)=y0.``

    The "t" column of the table increments from `t_0` to `t_1` by `h`
    (so `\\frac{t_1-t_0}{h}` must be an integer). In the "x" column,
    the new x-value equals the old x-value plus the corresponding
    entry in the next (third) column.  In the "y" column, the new
    y-value equals the old y-value plus the corresponding entry in the
    next (last) column.

    *For pedagogical purposes only.*

    EXAMPLES::

        sage: from sage.calculus.desolvers import eulers_method_2x2
        sage: t, x, y = PolynomialRing(QQ,3,"txy").gens()
        sage: f = x+y+t; g = x-y
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1,algorithm="none")
        [[0, 0, 0], [1/3, 0, 0], [2/3, 1/9, 0], [1, 10/27, 1/27], [4/3, 68/81, 4/27]]

    ::

        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    0                         0                    0                    0
           1/3                    0                       1/9                    0                    0
           2/3                  1/9                      7/27                    0                 1/27
             1                10/27                     38/81                 1/27                  1/9

    ::

        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: t,x,y=PolynomialRing(RR,3,"txy").gens()
        sage: f = x+y+t; g = x-y
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    0                      0.00                    0                 0.00
           1/3                 0.00                      0.13                 0.00                 0.00
           2/3                 0.13                      0.29                 0.00                0.043
             1                 0.41                      0.57                0.043                 0.15

    To numerically approximate `y(1)`, where `(1+t^2)y''+y'-y=0`,
    `y(0)=1`, `y'(0)=-1`, using 4 steps of Euler's method, first
    convert to a system: `y_1' = y_2`, `y_1(0)=1`; `y_2' =
    \\frac{y_1-y_2}{1+t^2}`, `y_2(0)=-1`.::

         sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
         sage: t, x, y=PolynomialRing(RR,3,"txy").gens()
         sage: f = y; g = (x-y)/(1+t^2)
         sage: eulers_method_2x2(f,g, 0, 1, -1, 1/4, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    1                     -0.25                   -1                 0.50
           1/4                 0.75                     -0.12                -0.50                 0.29
           1/2                 0.63                    -0.054                -0.21                 0.19
           3/4                 0.63                   -0.0078               -0.031                 0.11
             1                 0.63                     0.020                0.079                0.071

    To numerically approximate y(1), where `y''+ty'+y=0`, `y(0)=1`, `y'(0)=0`::

        sage: t,x,y=PolynomialRing(RR,3,"txy").gens()
        sage: f = y; g = -x-y*t
        sage: eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    1                      0.00                    0                -0.25
           1/4                  1.0                    -0.062                -0.25                -0.23
           1/2                 0.94                     -0.11                -0.46                -0.17
           3/4                 0.88                     -0.15                -0.62                -0.10
             1                 0.75                     -0.17                -0.68               -0.015

    AUTHORS:

    - David Joyner
    """
    if algorithm=="table":
        print("%10s %20s %25s %20s %20s"%("t", "x","h*f(t,x,y)","y", "h*g(t,x,y)"))
    n=int((1.0)*(t1-t0)/h)
    t00 = t0; x00 = x0; y00 = y0
    soln = [[t00,x00,y00]]
    for i in range(n+1):
        if algorithm=="table":
            print("%10r %20r %25r %20r %20r"%(t00,x00,h*f(t00,x00,y00),y00,h*g(t00,x00,y00)))
        x01 = x00 + h*f(t00,x00,y00)
        y00 = y00 + h*g(t00,x00,y00)
        x00 = x01
        t00 = t00 + h
        soln.append([t00,x00,y00])
    if algorithm!="table":
        return soln

def eulers_method_2x2_plot(f,g, t0, x0, y0, h, t1):
    r"""
    Plots solution of ODE

    This plots the soln in the rectangle ``(xrange[0],xrange[1])
    x (yrange[0],yrange[1])`` and plots using Euler's method the
    numerical solution of the 1st order ODEs `x' = f(t,x,y)`,
    `x(a)=x_0`, `y' = g(t,x,y)`, `y(a) = y_0`.

    *For pedagogical purposes only.*

    EXAMPLES::

        sage: from sage.calculus.desolvers import eulers_method_2x2_plot

    The following example plots the solution to
    `\theta''+\sin(\theta)=0`, `\theta(0)=\frac 34`, `\theta'(0) =
    0`.  Type ``P[0].show()`` to plot the solution,
    ``(P[0]+P[1]).show()`` to plot `(t,\theta(t))` and
    `(t,\theta'(t))`::

        sage: f = lambda z : z[2]; g = lambda z : -sin(z[1])
        sage: P = eulers_method_2x2_plot(f,g, 0.0, 0.75, 0.0, 0.1, 1.0)
    """
    n=int((1.0)*(t1-t0)/h)
    t00 = t0; x00 = x0; y00 = y0
    soln = [[t00,x00,y00]]
    for i in range(n+1):
        x01 = x00 + h*f([t00,x00,y00])
        y00 = y00 + h*g([t00,x00,y00])
        x00 = x01
        t00 = t00 + h
        soln.append([t00,x00,y00])
    Q1 = line([[x[0],x[1]] for x in soln], rgbcolor=(1/4,1/8,3/4))
    Q2 = line([[x[0],x[2]] for x in soln], rgbcolor=(1/2,1/8,1/4))
    return [Q1,Q2]

def desolve_rk4_determine_bounds(ics,end_points=None):
    """
    Used to determine bounds for numerical integration.

    - If end_points is None, the interval for integration is from ics[0]
      to ics[0]+10

    - If end_points is a or [a], the interval for integration is from min(ics[0],a)
      to max(ics[0],a)

    - If end_points is [a,b], the interval for integration is from min(ics[0],a)
      to max(ics[0],b)

    EXAMPLES::

        sage: from sage.calculus.desolvers import desolve_rk4_determine_bounds
        sage: desolve_rk4_determine_bounds([0,2],1)
        (0, 1)

    ::

        sage: desolve_rk4_determine_bounds([0,2])
        (0, 10)

    ::

        sage: desolve_rk4_determine_bounds([0,2],[-2])
        (-2, 0)

    ::

        sage: desolve_rk4_determine_bounds([0,2],[-2,4])
        (-2, 4)

    """
    if end_points is None:
        return((ics[0],ics[0]+10))
    if not isinstance(end_points,list):
        end_points=[end_points]
    if len(end_points)==1:
        return (min(ics[0],end_points[0]),max(ics[0],end_points[0]))
    else:
        return (min(ics[0],end_points[0]),max(ics[0],end_points[1]))


def desolve_rk4(de, dvar, ics=None, ivar=None, end_points=None, step=0.1, output='list', **kwds):
    """
    Solves numerically one first-order ordinary differential
    equation. See also ``ode_solver``.

    INPUT:

    input is similar to ``desolve`` command. The differential equation can be
    written in a form close to the plot_slope_field or desolve command

    - Variant 1 (function in two variables)

      - ``de`` - right hand side, i.e. the function `f(x,y)` from ODE `y'=f(x,y)`

      - ``dvar`` - dependent variable (symbolic variable declared by var)

    - Variant 2 (symbolic equation)

      - ``de`` - equation, including term with ``diff(y,x)``

      - ``dvar``` - dependent variable (declared as funciton of independent variable)

    - Other parameters

      - ``ivar`` - should be specified, if there are more variables or if the equation is autonomous

      - ``ics`` - initial conditions in the form [x0,y0]

      - ``end_points`` - the end points of the interval

        - if end_points is a or [a], we integrate on between min(ics[0],a) and max(ics[0],a)
        - if end_points is None, we use end_points=ics[0]+10

        - if end_points is [a,b] we integrate on between min(ics[0],a) and max(ics[0],b)

      - ``step`` - (optional, default:0.1) the length of the step (positive number)

      - ``output`` - (optional, default: 'list') one of 'list',
        'plot', 'slope_field' (graph of the solution with slope field)

    OUTPUT:

    Returns a list of points, or plot produced by list_plot,
    optionally with slope field.


    EXAMPLES::

        sage: from sage.calculus.desolvers import desolve_rk4

    Variant 2 for input - more common in numerics::

        sage: x,y=var('x y')
        sage: desolve_rk4(x*y*(2-y),y,ics=[0,1],end_points=1,step=0.5)
        [[0, 1], [0.5, 1.12419127425], [1.0, 1.46159016229]]

    Variant 1 for input - we can pass ODE in the form used by
    desolve function In this example we integrate bakwards, since
    ``end_points < ics[0]``::

        sage: y=function('y',x)
        sage: desolve_rk4(diff(y,x)+y*(y-1) == x-2,y,ics=[1,1],step=0.5, end_points=0)
        [[0.0, 8.90425710896], [0.5, 1.90932794536], [1, 1]]

    Here we show how to plot simple pictures. For more advanced
    aplications use list_plot instead. To see the resulting picture
    use ``show(P)`` in Sage notebook. ::

        sage: x,y=var('x y')
        sage: P=desolve_rk4(y*(2-y),y,ics=[0,.1],ivar=x,output='slope_field',end_points=[-4,6],thickness=3)

    ALGORITHM:

    4th order Runge-Kutta method. Wrapper for command ``rk`` in
    Maxima's dynamics package.  Perhaps could be faster by using
    fast_float instead.

    AUTHORS:

    - Robert Marik (10-2009)
    """
    if ics is None:
        raise ValueError("No initial conditions, specify with ics=[x0,y0].")

    if ivar is None:
        ivars = de.variables()
        ivars = [t for t in ivars if t != dvar]
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = ivars[0]

    if not is_SymbolicVariable(dvar):
        from sage.calculus.var import var
        from sage.calculus.all import diff
        from sage.symbolic.relation import solve
        if is_SymbolicEquation(de):
            de = de.lhs() - de.rhs()
        dummy_dvar=var('dummy_dvar')
        # consider to add warning if the solution is not unique
        de=solve(de,diff(dvar,ivar),solution_dict=True)
        if len(de) != 1:
            raise NotImplementedError("Sorry, cannot find explicit formula for right-hand side of the ODE.")
        de=de[0][diff(dvar,ivar)].subs(dvar==dummy_dvar)
    else:
        dummy_dvar=dvar

    step=abs(step)
    de0=de._maxima_()
    maxima("load('dynamics)")
    lower_bound,upper_bound=desolve_rk4_determine_bounds(ics,end_points)
    sol_1, sol_2 = [],[]
    if lower_bound<ics[0]:
        cmd="rk(%s,%s,%s,[%s,%s,%s,%s])\
        "%(de0.str(),'_SAGE_VAR_'+str(dummy_dvar),str(ics[1]),'_SAGE_VAR_'+str(ivar),str(ics[0]),lower_bound,-step)
        sol_1=maxima(cmd).sage()
        sol_1.pop(0)
        sol_1.reverse()
    if upper_bound>ics[0]:
        cmd="rk(%s,%s,%s,[%s,%s,%s,%s])\
        "%(de0.str(),'_SAGE_VAR_'+str(dummy_dvar),str(ics[1]),'_SAGE_VAR_'+str(ivar),str(ics[0]),upper_bound,step)
        sol_2=maxima(cmd).sage()
        sol_2.pop(0)
    sol=sol_1
    sol.extend([[ics[0],ics[1]]])
    sol.extend(sol_2)

    if output=='list':
        return sol
    from sage.plot.plot import list_plot
    from sage.plot.plot_field import plot_slope_field
    R = list_plot(sol,plotjoined=True,**kwds)
    if output=='plot':
        return R
    if output=='slope_field':
        XMIN=sol[0][0]
        YMIN=sol[0][1]
        XMAX=XMIN
        YMAX=YMIN
        for s,t in sol:
            if s>XMAX:XMAX=s
            if s<XMIN:XMIN=s
            if t>YMAX:YMAX=t
            if t<YMIN:YMIN=t
        return plot_slope_field(de,(ivar,XMIN,XMAX),(dummy_dvar,YMIN,YMAX))+R

    raise ValueError("Option output should be 'list', 'plot' or 'slope_field'.")

def desolve_system_rk4(des, vars, ics=None, ivar=None, end_points=None, step=0.1):
    r"""
    Solves numerically system of first-order ordinary differential
    equations using the 4th order Runge-Kutta method. Wrapper for
    Maxima command ``rk``. See also ``ode_solver``.

    INPUT:

    input is similar to desolve_system and desolve_rk4 commands

    - ``des`` - right hand sides of the system

    - ``vars`` - dependent variables

    - ``ivar`` - (optional) should be specified, if there are more variables or
      if the equation is autonomous and the independent variable is
      missing

    - ``ics`` - initial conditions in the form [x0,y01,y02,y03,....]

    - ``end_points`` - the end points of the interval

      - if end_points is a or [a], we integrate on between min(ics[0],a) and max(ics[0],a)
      - if end_points is None, we use end_points=ics[0]+10

      - if end_points is [a,b] we integrate on between min(ics[0],a) and max(ics[0],b)

    - ``step`` -- (optional, default: 0.1) the length of the step

    OUTPUT:

    Returns a list of points.

    EXAMPLES::

        sage: from sage.calculus.desolvers import desolve_system_rk4

    Lotka Volterra system::

        sage: from sage.calculus.desolvers import desolve_system_rk4
        sage: x,y,t=var('x y t')
        sage: P=desolve_system_rk4([x*(1-y),-y*(1-x)],[x,y],ics=[0,0.5,2],ivar=t,end_points=20)
        sage: Q=[ [i,j] for i,j,k in P]
        sage: LP=list_plot(Q)

        sage: Q=[ [j,k] for i,j,k in P]
        sage: LP=list_plot(Q)

    ALGORITHM:

    4th order Runge-Kutta method. Wrapper for command ``rk`` in Maxima's
    dynamics package.  Perhaps could be faster by using ``fast_float``
    instead.

    AUTHOR:

    - Robert Marik (10-2009)
    """

    if ics is None:
        raise ValueError("No initial conditions, specify with ics=[x0,y01,y02,...].")

    ivars = set([])

    for de in des:
        ivars = ivars.union(set(de.variables()))
    if ivar is None:
        ivars = ivars - set(vars)
        if len(ivars) != 1:
            raise ValueError("Unable to determine independent variable, please specify.")
        ivar = list(ivars)[0]

    dess = [de._maxima_().str() for de in des]
    desstr = "[" + ",".join(dess) + "]"
    varss = [varsi._maxima_().str() for varsi in vars]
    varstr = "[" + ",".join(varss) + "]"
    x0=ics[0]
    icss = [ics[i]._maxima_().str() for i in range(1,len(ics))]
    icstr = "[" + ",".join(icss) + "]"
    step=abs(step)

    maxima("load('dynamics)")
    lower_bound,upper_bound=desolve_rk4_determine_bounds(ics,end_points)
    sol_1, sol_2 = [],[]
    if lower_bound<ics[0]:
        cmd="rk(%s,%s,%s,[%s,%s,%s,%s])\
        "%(desstr,varstr,icstr,'_SAGE_VAR_'+str(ivar),str(x0),lower_bound,-step)
        sol_1=maxima(cmd).sage()
        sol_1.pop(0)
        sol_1.reverse()
    if upper_bound>ics[0]:
        cmd="rk(%s,%s,%s,[%s,%s,%s,%s])\
        "%(desstr,varstr,icstr,'_SAGE_VAR_'+str(ivar),str(x0),upper_bound,step)
        sol_2=maxima(cmd).sage()
        sol_2.pop(0)
    sol=sol_1
    sol.append(ics)
    sol.extend(sol_2)

    return sol

def desolve_odeint(des, ics, times, dvars, ivar=None, compute_jac=False, args=()
, rtol=None, atol=None, tcrit=None, h0=0.0, hmax=0.0, hmin=0.0, ixpr=0
, mxstep=0, mxhnil=0, mxordn=12, mxords=5, printmessg=0):
    r"""
    Solves numerically a system of first-order ordinary differential equations
    using ``odeint`` from scipy.integrate module.

    INPUT:

    - ``des``  -- right hand sides of the system

    - ``ics``  -- initial conditions

    - ``times`` -- a sequence of time points in which the solution must be found

    - ``dvars`` -- dependent variables. ATTENTION: the order must be the same as
      in des, that means: d(dvars[i])/dt=des[i]

    - ``ivar`` -- independent variable, optional.

    - ``compute_jac`` -- boolean. If True, the Jacobian of des is computed and
      used during the integration of Stiff Systems. Default value is False.

    Other Parameters (taken from the documentation of odeint function from
      scipy.integrate module)

    - ``rtol``, ``atol`` : float
      The input parameters rtol and atol determine the error
      control performed by the solver.  The solver will control the
      vector, e, of estimated local errors in y, according to an
      inequality of the form:

        max-norm of (e / ewt) <= 1

      where ewt is a vector of positive error weights computed as:

        ewt = rtol * abs(y) + atol

      rtol and atol can be either vectors the same length as y or scalars.

    - ``tcrit`` : array
      Vector of critical points (e.g. singularities) where integration
      care should be taken.

    - ``h0`` : float, (0: solver-determined)
      The step size to be attempted on the first step.

    - ``hmax`` : float, (0: solver-determined)
      The maximum absolute step size allowed.

    - ``hmin`` : float, (0: solver-determined)
      The minimum absolute step size allowed.

    - ``ixpr`` : boolean.
      Whether to generate extra printing at method switches.

    - ``mxstep`` : integer, (0: solver-determined)
      Maximum number of (internally defined) steps allowed for each
      integration point in t.

    - ``mxhnil`` : integer, (0: solver-determined)
      Maximum number of messages printed.

    - ``mxordn`` : integer, (0: solver-determined)
      Maximum order to be allowed for the nonstiff (Adams) method.

    - ``mxords`` : integer, (0: solver-determined)
      Maximum order to be allowed for the stiff (BDF) method.

    OUTPUT:

    Returns a list with the solution of the system at each time in times.

    EXAMPLES:

    Lotka Volterra Equations::

        sage: from sage.calculus.desolvers import desolve_odeint
        sage: x,y=var('x,y')
        sage: f=[x*(1-y),-y*(1-x)]
        sage: sol=desolve_odeint(f,[0.5,2],srange(0,10,0.1),[x,y])
        sage: p=line(zip(sol[:,0],sol[:,1]))
        sage: p.show()

    Lorenz Equations::

        sage: x,y,z=var('x,y,z')
        sage: # Next we define the parameters
        sage: sigma=10
        sage: rho=28
        sage: beta=8/3
        sage: # The Lorenz equations
        sage: lorenz=[sigma*(y-x),x*(rho-z)-y,x*y-beta*z]
        sage: # Time and initial conditions
        sage: times=srange(0,50.05,0.05)
        sage: ics=[0,1,1]
        sage: sol=desolve_odeint(lorenz,ics,times,[x,y,z],rtol=1e-13,atol=1e-14)

    One-dimensional Stiff system::

        sage: y= var('y')
        sage: epsilon=0.01
        sage: f=y^2*(1-y)
        sage: ic=epsilon
        sage: t=srange(0,2/epsilon,1)
        sage: sol=desolve_odeint(f,ic,t,y,rtol=1e-9,atol=1e-10,compute_jac=True)
        sage: p=points(zip(t,sol))
        sage: p.show()

    Another Stiff system with some optional parameters with no
    default value::

        sage: y1,y2,y3=var('y1,y2,y3')
        sage: f1=77.27*(y2+y1*(1-8.375*1e-6*y1-y2))
        sage: f2=1/77.27*(y3-(1+y1)*y2)
        sage: f3=0.16*(y1-y3)
        sage: f=[f1,f2,f3]
        sage: ci=[0.2,0.4,0.7]
        sage: t=srange(0,10,0.01)
        sage: v=[y1,y2,y3]
        sage: sol=desolve_odeint(f,ci,t,v,rtol=1e-3,atol=1e-4,h0=0.1,hmax=1,hmin=1e-4,mxstep=1000,mxords=17)

    AUTHOR:

    - Oriol Castejon (05-2010)
    """

    from scipy.integrate import odeint
    from sage.ext.fast_eval import fast_float
    from sage.calculus.functions import jacobian

    if ivar==None:
        if len(dvars)==0 or len(dvars)==1:
            if len(dvars)==1:
                des=des[0]
                dvars=dvars[0]
            all_vars = set(des.variables())
        else:
            all_vars = set([])
            for de in des:
                all_vars.update(set(de.variables()))
        if is_SymbolicVariable(dvars):
            ivars = all_vars - set([dvars])
        else:
            ivars = all_vars - set(dvars)

        if len(ivars)==1:
            ivar = ivars.pop()
        elif not ivars:
            from sage.symbolic.ring import var
            try:
                safe_names = [ 't_' + str(dvar) for dvar in dvars ]
            except TypeError:  # not iterable
                safe_names = [ 't_' + str(dvars) ]
            ivar = map(var, safe_names)
        else:
            raise ValueError("Unable to determine independent variable, please specify.")

    # one-dimensional systems:
    if is_SymbolicVariable(dvars):
        func = fast_float(des,dvars,ivar)
        if not compute_jac:
            Dfun=None
        else:
            J = diff(des,dvars)
            J = fast_float(J,dvars,ivar)
            Dfun = lambda y,t: [J(y,t)]

    # n-dimensional systems:
    else:
        desc = []
        variabs = dvars[:]
        variabs.append(ivar)
        for de in des:
            desc.append(fast_float(de,*variabs))

        def func(y,t):
            v = list(y[:])
            v.append(t)
            return [dec(*v) for dec in desc]

        if not compute_jac:
            Dfun=None
        else:
            J = jacobian(des,dvars)
            J = [list(v) for v in J]
            J = fast_float(J,*variabs)
            def Dfun(y,t):
                v = list(y[:])
                v.append(t)
                return [[element(*v) for element in row] for row in J]


    sol=odeint(func, ics, times, args=args, Dfun=Dfun, rtol=rtol, atol=atol,
        tcrit=tcrit, h0=h0, hmax=hmax, hmin=hmin, ixpr=ixpr, mxstep=mxstep,
        mxhnil=mxhnil, mxordn=mxordn, mxords=mxords, printmessg=printmessg)

    return sol

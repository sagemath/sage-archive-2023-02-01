"""
This file contains functions useful for solving differential equations
which occur commonly in a 1st semester differential equations course.

* desolve -- Computes the "general solution" to a 1st or 2nd order
             ODE via maxima.

* desolve_laplace -- Solves an ODE using laplace transforms via maxima.
                     Initials conditions are optional.

* desolve_system -- Solves any size system of 1st order odes using maxima.
                    Initials conditions are optional.

* eulers_method -- Approximate solution to a 1st order DE, presented as a table.

* eulers_method_2x2 -- Approximate solution to a 1st order system of DEs,
                       presented as a table.

* eulers_method_2x2_plot -- Plots the sequence of points obtained from
                    Euler's method.

AUTHORS: David Joyner (3-2006)     -- Initial version of functions
         Marshall Hampton (7-2007) -- Creation of Python module and testing
         Robert Bradshaw (10-2008) -- Some interface cleanup.

Other functions for solving DEs are given in functions/elementary.py.

"""

##########################################################################
#  Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>, Marshall Hampton
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  http://www.gnu.org/licenses/
##########################################################################

from sage.interfaces.maxima import MaximaElement, Maxima
from sage.plot.all import line
from sage.symbolic.expression import is_SymbolicEquation

maxima = Maxima()

def desolve(de, dvar, ics=None, ivar=None):
    """
    Solves a 1st or 2nd order linear ODE via maxima. Initials conditions
    are not given.

    INPUT:
        de    -- an expression or equation representing the ODE
        dvar  -- the dependent variable (hereafter called y)
        ics   -- (optional) the initial or boundary conditions
                    for a first-order equation, specify the initial x and y
                    for a second-order equation, specify the initial x, y, and dy/dx
                    for a second-order boundary solution, specify initial and final x and y
        ivar  -- (optional) the independent variable  (hereafter called x),
                    which must be specified if there is more than one
                    independent variable in the equation.

    EXAMPLES:
        sage: x = var('x')
        sage: y = function('y', x)
        sage: desolve(diff(y,x) + y - 1, y)
        (c + e^x)*e^(-x)
        sage: f = desolve(diff(y,x) + y - 1, y, ics=[10,2]); f
        (e^10 + e^x)*e^(-x)
        sage: plot(f)

    We can also solve second-order differential equations.
        sage: x = var('x')
        sage: y = function('y', x)
        sage: de = diff(y,x,2) - y == x
        sage: desolve(de, y)
        k1*e^x + k2*e^(-x) - x
        sage: f = desolve(de, y, [10,2,1]); f
        1/2*(y(10) + 12)*e^(x - 10) + 1/2*(e^10*y(10) + 8*e^10)*e^(-x) - x
        sage: f(x=10).expand().simplify()
        y(10)
        sage: diff(f,x)(x=10).expand().simplify()
        1

    AUTHOR: David Joyner (1-2006)
            Robert Bradshaw (10-2008)
    """
    if is_SymbolicEquation(de):
        de = de.lhs() - de.rhs()
    # for backwards compatibility
    if isinstance(dvar, list):
        dvar, ivar = dvar
    elif ivar is None:
        ivars = de.variables()
        ivars = [t for t in ivars if t != dvar]
        if len(ivars) != 1:
            raise ValueError, "Unable to determine independent variable, please specify."
        ivar = ivars[0]
    de0 = de._maxima_()
    soln = de0.ode2(dvar, ivar)
    if str(soln).strip() == 'false':
        raise NotImplementedError, "Maxima was unable to solve this system."
    def to_eqns(lhs, exprs):
        eqns = []
        for lhs, expr in zip(lhs, exprs):
            if is_SymbolicEquation(expr):
                eqns.append(expr)
            else:
                if lhs == dvar and len(exprs) == 2:
                    ivar_ic = exprs[0] # figure this out...
                    lhs = lhs.subs(lhs.default_variable()==ivar_ic)
                eqns.append(lhs == expr)
        return eqns
    if ics is not None:
        if len(ics) == 2:
            ics = to_eqns([ivar, dvar], ics)
            soln = soln.ic1(ics[0], ics[1])
        if len(ics) == 3:
            ics = to_eqns([ivar, dvar, dvar.derivative(ivar)], ics)
            soln = soln.ic2(*ics)
        if len(ics) == 4:
            ics = to_eqns([ivar, dvar, ivar, dvar], ics)
            soln = soln.bc(*ics)
    if soln.lhs() == dvar:
        soln = soln.rhs()
    return soln.sage()

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

def desolve_laplace(de,vars,ics=None):
    """
    Solves an ODE using laplace transforms via maxima. Initials conditions
    are optional.

    INPUT:
        de    -- a lambda expression representing the ODE
                 (eg, de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")
        vars  -- a list of strings representing the variables
                 (eg, vars = ["x","f"], if x is the independent
                  variable and f is the dependent variable)
        ics   -- a list of numbers representing initial conditions,
                 with symbols allowed which are represented by strings
                 (eg, f(0)=1, f'(0)=2 is ics = [0,1,2])

    EXAMPLES:
        sage: from sage.calculus.desolvers import desolve_laplace
        sage: x = var('x')
        sage: function('f', x)
        f(x)
        sage: de = lambda y: diff(y,x,x) - 2*diff(y,x) + y
        sage: desolve_laplace(de(f(x)),["x","f"])
        "x*%e^x*('at('diff('f(x),x,1),x=0))-'f(0)*x*%e^x+'f(0)*%e^x"
        sage: desolve_laplace(de(f(x)),["x","f"],[0,1,2])
         'x*%e^x+%e^x'

    WARNING:
        The second Sage command in the above example sets the values of f(0) and f'(0)
        in Maxima, so subsequent ODEs involving these variables will have these initial conditions
        automatically imposed.

    AUTHOR: David Joyner (1-2006,8-2007)
    """
    if ics!=None:
        d = len(ics)
        for i in range(0,d-1):
            ic = "atvalue(diff('"+vars[1]+"("+vars[0]+"),"+str(vars[0])+","+str(i)+"),"+str(vars[0])+"="+str(ics[0])+","+str(ics[1+i])+")"
            maxima(ic)
            #print i,ic
    cmd = "desolve("+de._maxima_init_()+",'"+vars[1]+"("+vars[0]+"));"
    return maxima(cmd).rhs()._maxima_init_()

def desolve_system(des, vars, ics=None, ivar=None):
    """
    Solves any size system of 1st order odes using maxima. Initials conditions
    are optional.

    INPUT:
        des   -- list of ODEs
        vars  -- list of dependent variables
        ics   -- (optional) list of initial values for ivar and vars
        ivar  -- (optional) the independent variable, which must be specified
                    if there is more than one independent variable in the equation.

    EXAMPLES:
        sage: t = var('t')
        sage: x = function('x', t)
        sage: y = function('y', t)
        sage: de1 = diff(x,t) + y - 1 == 0
        sage: de2 = diff(y,t) - x + 1 == 0
        sage: desolve_system([de1, de2], [x,y])
        [x(t) == (x(0) - 1)*cos(t) - (y(0) - 1)*sin(t) + 1,
         y(t) == (x(0) - 1)*sin(t) + (y(0) - 1)*cos(t) + 1]

    Now we give some initial conditions:
        sage: sol = desolve_system([de1, de2], [x,y], ics=[0,1,2]); sol
        [x(t) == -sin(t) + 1, y(t) == cos(t) + 1]
        sage: solnx, solny = sol[0].rhs(), sol[1].rhs()
        sage: plot([solnx,solny],0,1)
        sage: parametric_plot((solnx,solny),0,1)

    AUTHOR: Robert Bradshaw (10-2008)
    """
    ivars = set([])
    for i, de in enumerate(des):
        if not is_SymbolicEquation(de):
            des[i] = de == 0
        ivars = ivars.union(set(de.variables()))
    if ivar is None:
        ivars = ivars - set(vars)
        if len(ivars) != 1:
            raise ValueError, "Unable to determine independent variable, please specify."
        ivar = list(ivars)[0]
    dvars = [v._maxima_() for v in vars]
    if ics is not None:
        ivar_ic = ics[0]
        for dvar, ic in zip(dvars, ics[1:]):
            dvar.atvalue(ivar==ivar_ic, ic)
    soln = dvars[0].parent().desolve(des, dvars)
    if str(soln).strip() == 'false':
        raise NotImplementedError, "Maxima was unable to solve this system."
    soln = list(soln)
    for i, sol in enumerate(soln):
        soln[i] = sol.sage()
    if ics is not None:
        ivar_ic = ics[0]
        for dvar, ic in zip(dvars, ics[:1]):
            dvar.atvalue(ivar==ivar_ic, dvar)
    return soln


def desolve_system_strings(des,vars,ics=None):
    """
    Solves any size system of 1st order odes using maxima. Initials conditions
    are optional.


    INPUT:
        de    -- a list of strings representing the ODEs in maxima notation
                 (eg, de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")
        vars  -- a list of strings representing the variables
                 (eg, vars = ["s","x","y"], where s is the independent variable
                  and x,y the dependent variables)
        ics   -- a list of numbers representing initial conditions
                 (eg, x(0)=1, y(0)=2 is ics = [0,1,2])

    WARNING:
        The given ics sets the initial values of the dependent vars in maxima, so
        subsequent ODEs involving these variables will have these initial conditions
        automatically imposed.

    EXAMPLES:
        sage: from sage.calculus.desolvers import desolve_system_strings
        sage: s = var('s')
        sage: function('x', s)
        x(s)
        sage: function('y', s)
        y(s)
        sage: de1 = lambda z: diff(z[0],s) + z[1] - 1
        sage: de2 = lambda z: diff(z[1],s) - z[0] + 1
        sage: des = [de1([x(s),y(s)]),de2([x(s),y(s)])]
        sage: vars = ["s","x","y"]
        sage: desolve_system_strings(des,vars)
        ["(1-'y(0))*sin(s)+('x(0)-1)*cos(s)+1", "('x(0)-1)*sin(s)+('y(0)-1)*cos(s)+1"]
        sage: ics = [0,1,-1]
        sage: soln = desolve_system_strings(des,vars,ics); soln
        ['2*sin(s)+1', '1-2*cos(s)']
        sage: solnx, solny = map(SR, soln)
        sage: RR(solnx(s=3))
        1.28224001611973
        sage: P1 = plot([solnx,solny],0,1)
        sage: P2 = parametric_plot((solnx,solny),0,1)

        Now type show(P1), show(P2) to view these.


    AUTHOR: David Joyner (3-2006, 8-2007)
    """
    d = len(des)
    dess = [de._maxima_init_() + "=0" for de in des]
    for i in range(d):
        cmd="de:" + dess[int(i)] + ";"
        maxima.eval(cmd)
    desstr = "[" + ",".join(dess) + "]"
    d = len(vars)
    varss = list("'" + vars[i] + "(" + vars[0] + ")" for i in range(1,d))
    varstr = "[" + ",".join(varss) + "]"
    if ics is not None:
        #d = len(ics) ## must be same as len(des)
        for i in range(1,d):
            ic = "atvalue('" + vars[i] + "("+vars[0] + ")," + str(vars[0]) + "=" + str(ics[0]) + "," + str(ics[i]) + ")"
            maxima.eval(ic)
    cmd = "desolve(" + desstr + "," + varstr + ");"
    soln = maxima(cmd)
    return [f.rhs()._maxima_init_() for f in soln]

def eulers_method(f,x0,y0,h,x1,method="table"):
    """
    This implements Euler's method for finding numerically the solution of the 1st order
    ODE y' = f(x,y), y(a)=c. The "x" column of the table increments from
    x0 to x1 by h (so (x1-x0)/h must be an integer). In the "y" column,
    the new y-value equals the old y-value plus the corresponding entry in the
    last column.

    *For pedagogical purposes only.*

    EXAMPLES:
        sage: from sage.calculus.desolvers import eulers_method
        sage: x,y = PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1)
             x                    y                  h*f(x,y)
             0                    1                   -2
           1/2                   -1                 -7/4
             1                -11/4                -11/8
        sage: x,y = PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1,method="none")
        [[0, 1], [1/2, -1], [1, -11/4], [3/2, -33/8]]
        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y = PolynomialRing(RR,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1,method="None")
        [[0, 1], [1/2, -1.0], [1, -2.7], [3/2, -4.0]]
        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y=PolynomialRing(RR,2,"xy").gens()
        sage: eulers_method(5*x+y-5,0,1,1/2,1)
             x                    y                  h*f(x,y)
             0                    1                 -2.0
           1/2                 -1.0                 -1.7
             1                 -2.7                 -1.3
        sage: x,y=PolynomialRing(QQ,2,"xy").gens()
        sage: eulers_method(5*x+y-5,1,1,1/3,2)
                 x                    y                  h*f(x,y)
                 1                    1                  1/3
               4/3                  4/3                    1
               5/3                  7/3                 17/9
                 2                 38/9                83/27
        sage: eulers_method(5*x+y-5,0,1,1/2,1,method="none")
        [[0, 1], [1/2, -1], [1, -11/4], [3/2, -33/8]]
        sage: pts = eulers_method(5*x+y-5,0,1,1/2,1,method="none")
        sage: P1 = list_plot(pts)
        sage: P2 = line(pts)
        sage: (P1+P2).show()

    AUTHOR: David Joyner
    """
    if method=="table":
        print "%10s %20s %25s"%("x","y","h*f(x,y)")
    n=int((1.0)*(x1-x0)/h)
    x00=x0; y00=y0
    soln = [[x00,y00]]
    for i in range(n+1):
        if method=="table":
            print "%10r %20r %20r"%(x00,y00,h*f(x00,y00))
        y00 = y00+h*f(x00,y00)
        x00=x00+h
        soln.append([x00,y00])
    if method!="table":
        return soln

def eulers_method_2x2(f,g, t0, x0, y0, h, t1,method="table"):
    """
    This implements Euler's method for finding numerically the solution of the 1st order
    system of two ODEs
        x' = f(t, x, y), x(t0)=x0.
        y' = g(t, x, y), y(t0)=y0.
    The "t" column of the table increments from t0 to t1 by
    h (so (t1-t0)/h must be an integer). In the "x" column, the new x-value equals the
    old x-value plus the corresponding entry in the next (third) column.
    In the "y" column, the new y-value equals the old y-value plus the
    corresponding entry in the next (last) column.

    *For pedagogical purposes only.*

    EXAMPLES:
        sage: from sage.calculus.desolvers import eulers_method_2x2
        sage: t, x, y = PolynomialRing(QQ,3,"txy").gens()
        sage: f = x+y+t; g = x-y
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1,method="none")
        [[0, 0, 0], [1/3, 0, 0], [2/3, 1/9, 0], [1, 10/27, 1/27], [4/3, 68/81, 4/27]]
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    0                         0                    0                    0
           1/3                    0                       1/9                    0                    0
           2/3                  1/9                      7/27                    0                 1/27
             1                10/27                     38/81                 1/27                  1/9
        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: t,x,y=PolynomialRing(RR,3,"txy").gens()
        sage: f = x+y+t; g = x-y
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    0                      0.00                    0                 0.00
           1/3                 0.00                      0.13                 0.00                 0.00
           2/3                 0.13                      0.29                 0.00                0.043
             1                 0.41                      0.57                0.043                 0.15

    To numerically approximate y(1), where (1+t^2)y''+y'-y=0, y(0)=1,y'(0)=-1,
    using 4 steps of Euler's method, first convert to a system:
    y1' = y2, y1(0)=1; y2' = (y1-y2)/(1+t^2), y2(0)=-1.

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

    To numerically approximate y(1), where y''+ty'+y=0, y(0)=1,y'(0)=0:

        sage: t,x,y=PolynomialRing(RR,3,"txy").gens()
        sage: f = y; g = -x-y*t
        sage: eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
             t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
             0                    1                      0.00                    0                -0.25
           1/4                  1.0                    -0.062                -0.25                -0.23
           1/2                 0.94                     -0.11                -0.46                -0.17
           3/4                 0.88                     -0.15                -0.62                -0.10
             1                 0.75                     -0.17                -0.68               -0.015

    AUTHOR: David Joyner
    """
    if method=="table":
        print "%10s %20s %25s %20s %20s"%("t", "x","h*f(t,x,y)","y", "h*g(t,x,y)")
    n=int((1.0)*(t1-t0)/h)
    t00 = t0; x00 = x0; y00 = y0
    soln = [[t00,x00,y00]]
    for i in range(n+1):
        if method=="table":
            print "%10r %20r %25r %20r %20r"%(t00,x00,h*f(t00,x00,y00),y00,h*g(t00,x00,y00))
        x01 = x00 + h*f(t00,x00,y00)
        y00 = y00 + h*g(t00,x00,y00)
        x00 = x01
        t00 = t00 + h
        soln.append([t00,x00,y00])
    if method!="table":
        return soln

def eulers_method_2x2_plot(f,g, t0, x0, y0, h, t1):
    """
    This plots the soln in the rectangle
    (xrange[0],xrange[1]) x (yrange[0],yrange[1])
    and plots using Euler's method the numerical solution of
    the 1st order ODEs x' = f(t,x,y), x(a)=x0, y' = g(t,x,y), y(a) = y0.

    *For pedagogical purposes only.*

    ===
    The following example plots the solution to theta''+sin(theta)=0,
    theta(0)=3/4, theta'(0) = 0.  Type P[0].show() to plot the solution,
    (P[0]+P[1]).show() to plot (t,theta(t)) and (t,theta'(t)).
    ===
    EXAMPLES:
        sage: from sage.calculus.desolvers import eulers_method_2x2_plot
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

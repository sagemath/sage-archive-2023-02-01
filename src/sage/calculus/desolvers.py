"""
This file contains functions useful for solving differential equations
which occur commonly in a 1st semster differential equations course.

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
from sage.plot.plot import LineFactory
line = LineFactory()

maxima = Maxima()

def desolve(de,vars):
    """
    Solves a 1st or 2nd order ODE via maxima. Initials conditions
    are not given.

    INPUT:
        de    -- a string representing the ODE
                 (eg, de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")
        vars  -- a list of strings representing the variables
                 (eg, vars = ["x","y"], if x is the independent
                  variable and y is the dependent variable)

    EXAMPLES:
        sage: desolvers.desolve("x**2*diff(y,x) + 3*y*x = sin(x)/x",["x","y"])
        y=(%c-cos(x))/x^3

    AUTHOR: David Joyner (1-2006)
    """
    maxima("depends("+vars[1]+","+vars[0]+");")
    maxima("de:"+de+";")
    cmd = "ode2("+de+","+vars[1]+","+vars[0]+");"
    print maxima.eval(cmd)
    #Maxima.eval(cmd)

def desolve_laplace(de,vars,ics=None):
    """
    Solves an ODE using laplace transforms via maxima. Initials conditions
    are optional.

    INPUT:
        de    -- a string representing the ODE
                 (eg, de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")
        vars  -- a list of strings representing the variables
                 (eg, vars = ["x","f"], if x is the independent
                  variable and f is the dependent variable)
        ics   -- a list of numbers representing initial conditions,
                 with symbols allowed which are represented by strings
                 (eg, f(0)=1, f'(0)=2 is ics = [0,1,2])

    EXAMPLES:
        sage: desolvers.desolve_laplace("diff(f(x),x,2)=2*diff(f(x),x)-f(x)",["x","f"])
        f(x)=x*%e^x*(?%at('diff(f(x),x,1),x=0))-f(0)*x*%e^x+f(0)*%e^x
        sage: desolvers.desolve_laplace("diff(f(x),x,2)=2*diff(f(x),x)-f(x)",["x","f"],[0,1,2])
        f(x)=x*%e^x+%e^x

    WARNING:
        The second SAGE command in the above example sets the values of f(0) and f'(0) in maxima, so
        subsequent ODEs involving these variables will have these initial conditions
        automatically imposed.

    AUTHOR: David Joyner (1-2006)
    """
    maxima("de:"+de+";")
    if ics!=None:
        d = len(ics)
        for i in range(0,d-1):
            ic = "atvalue(diff("+vars[1]+"("+vars[0]+"),"+str(vars[0])+","+str(i)+"),"+str(vars[0])+"="+str(ics[0])+","+str(ics[1+i])+")"
            maxima(ic)
            #print i,ic
    cmd = "desolve("+de+","+vars[1]+"("+vars[0]+"));"
    print maxima.eval(cmd)
    #maxima.eval(cmd)

def desolve_system(des,vars,ics=None):
    """
    Solves any size system of 1st order odes using maxima. Initials conditions
    are optional.


    INPUT:
        de    -- a list of strings representing the ODEs in maxima notation
                 (eg, de = "diff(f(x),x,2)=diff(f(x),x)+sin(x)")
        vars  -- a list of strings representing the variables
                 (eg, vars = ["t","x","y"], where t is the indeendent variable
                  and x,y the dependent variables)
        ics   -- a list of numbers representing initial conditions
                 (eg, x(0)=1, y(0)=2 is ics = [0,1,2])

    WARNING:
        The given ics sets the initial values of the dependent vars in maxima, so
        subsequent ODEs involving these variables will have these initial conditions
        automatically imposed.

    EXAMPLES:
        sage: des=["'diff(x(t),t)=-4*y(t)","'diff(y(t),t)=-x(t)"]
        sage: vars = ["t","x","y"]
        sage: desolvers.desolve_system(des,vars)
        [x(t)=(2*y(0)+x(0))*%e^-(2*t)/2-(2*y(0)-x(0))*%e^(2*t)/2,y(t)=(2*y(0)-x(0))*%e^(2*t)/4+(2*y(0)+x(0))*%e^-(2*t)/4]
	sage: des=["'diff(x(t),t)=-4*y(t)","'diff(y(t),t)=-x(t)"]
    	sage: vars = ["t","x","y"]
    	sage: desolvers.desolve_system(des,vars,[0,270,90])
        [x(t)=45*%e^(2*t)+225*%e^-(2*t),y(t)=225*%e^-(2*t)/2-45*%e^(2*t)/2]
    	sage: des=["'diff(x(t),t)=-4*y(t)","'diff(y(t),t)=-z(t)","'diff(z(t),t)=-2*x(t)"]
    	sage: vars = ["t","x","y","z"]
    	sage: ics=[0,270,90,100]
    	sage: desolvers.desolve_system(des,vars,ics)
    	[x(t)=%e^t*(260*cos(sqrt(3)*t)/3-80*sin(sqrt(3)*t)/sqrt(3))+550*%e^-(2*t)/3,y(t)=%e^t*(85*sin(sqrt(3)*t)/sqrt(3)-5*cos(sqrt(3)*t)/3)+275*%e^-(2*t)/3,z(t)=%e^t*(-90*sin(sqrt(3)*t)/sqrt(3)-250*cos(sqrt(3)*t)/3)+550*%e^-(2*t)/3]

    AUTHOR: David Joyner (3-2006)
    """
    d = len(des)
    for i in range(d):
        cmd="de:"+des[int(i)]+";"
        maxima(cmd)
    desstr = "["+",".join(des)+"]"
    #print desstr
    d = len(vars)
    varss = list(vars[i]+"("+vars[0]+")" for i in range(1,d))
    varstr = "["+",".join(varss)+"]"
    #print varstr
    if ics!=None:
        d = len(ics)
        for i in range(1,d):
            ic = "atvalue("+vars[int(i)]+"("+vars[0]+"),"+str(vars[0])+"="+str(ics[0])+","+str(ics[int(i)])+")"
            maxima(ic)
    cmd = "desolve("+desstr+","+varstr+");"
    #print cmd
    print maxima.eval(cmd)
    #maxima.eval(cmd)

def eulers_method(f,x0,y0,h,x1,method="table"):
    """
    This implements Euler's method for
    finding numerically the solution of the 1st order
    ODE y' = f(x,y), y(a)=c. The "x" column of the table
    increments from x0 to x1 by h (so (x1-x0)/h must be
    an integer). In the "y" column, the new y-value equals the
    old y-value plus the corresponding entry in the
    last column.

    EXAMPLES:
        sage: x,y = PolynomialRing(QQ,2,"xy").gens()
        sage: desolvers.eulers_method(5*x+y-5,0,1,1/2,1,method="none")
        [[0, 1], [1/2, -1], [1, -11/4], [3/2, -33/8]]
        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y = PolynomialRing(RR,2,"xy").gens()
        sage: desolvers.eulers_method(5*x+y-5,0,1,1/2,1,method="None")
        [[0, 1], [1/2, -1.0], [1, -2.7], [3/2, -4.0]]
        sage: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
        sage: x,y=PolynomialRing(RR,2,"xy").gens()
        sage: desolvers.eulers_method(5*x+y-5,0,1,1/2,1)
         x                    y                  h*f(x,y)
         0                    1                -2.00
       1/2                -1.00                -1.75
         1                -2.75                -1.37

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
    h (so (t1-t0)/h must be an integer).
    In the "x" column, the new x-value equals the
    old x-value plus the corresponding entry in the next
    (third) column. In the "y" column, the new y-value equals the
    old y-value plus the corresponding entry in the next
    (last) column.

    EXAMPLES:
        sage: t, x, y = PolynomialRing(QQ,3,"txy").gens()
        sage: f = x+y+t; g = x-y
        sage: eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1,method="none")
        [[0, 0, 0], [1/3, 0, 0], [2/3, 1/9, 0], [1, 10/27, 1/27], [4/3, 68/81, 4/27]]
        sage: desolvers.eulers_method_2x2(f,g, 0, 0, 0, 1/3, 1)
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
         0                    0                     0.000                    0                0.000
       1/3                0.000                     0.125                0.000                0.000
       2/3                0.125                     0.282                0.000               0.0430
         1                0.407                     0.563               0.0430                0.141

    To numerically approximate y(1), where (1+t^2)y''+y'-y=0, y(0)=1,y'(0)=-1,
    using 4 steps of Euler's method, first convert to a system:
    y1' = y2, y1(0)=1; y2' = (y1-y2)/(1+t^2), y2(0)=-1.

         sage.: RR = RealField(sci_not=0, prec=4, rnd='RNDU')
         sage.: t, x, y=PolynomialRing(RR,3,"txy").gens()
         sage.: f = y; g = (x-y)/(1+t^2)
         sage.:. eulers_method_2x2(f,g, 0, 1, -1, 1/4, 1)
         t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
         0                    1                    -0.250                   -1                0.500
       1/4                0.750                    -0.125               -0.500                0.282
       1/2                0.625                   -0.0546               -0.218                0.188
       3/4                0.625                  -0.00781              -0.0312                0.110
         1                0.625                    0.0196               0.0782               0.0704

    To numerically approximate y(1), where y''+ty'+y=0, y(0)=1,y'(0)=0:

        sage.: t,x,y=PolynomialRing(RR,3,"txy").gens()
        sage.: f = y; g = -x-y*t
        sage.:. eulers_method_2x2(f,g, 0, 1, 0, 1/4, 1)
         t                    x                h*f(t,x,y)                    y           h*g(t,x,y)
         0                    1                     0.000                    0               -0.250
       1/4                 1.00                   -0.0625               -0.250               -0.234
       1/2                0.938                    -0.117               -0.468               -0.171
       3/4                0.875                    -0.156               -0.625               -0.101
         1                0.750                    -0.171               -0.687              -0.0156

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
    ===
    The following example plots the solution to theta''+sin(theta)=0,
    theta(0)=3/4, theta'(0) = 0.  Type P[0].show() to plot the solution,
    (P[0]+P[1]).show() to plot (t,theta(t)) and (t,theta'(t)).
    ===
    EXAMPLES:
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

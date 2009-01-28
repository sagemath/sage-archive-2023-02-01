"""
Numerical Root Finding and Optimization

AUTHOR:
    -- William Stein (2007)
"""
from sage.modules.free_module_element import vector
from sage.rings.real_double import RDF

def find_root(f, a, b, xtol=10e-13, rtol=4.5e-16, maxiter=100, full_output=False):
    """
    Numerically find a root of f on the closed interval [a,b]
    (or [b,a]) if possible, where f is a function in the one variable.

    INPUT:
        f -- a function of one variable or symbolic equality
        a, b -- endpoints of the interval
        xtol, rtol -- the routine converges when a root is known
                to lie within xtol of the value return. Should be
                >= 0.  The routine modifies this to take into
                account the relative precision of doubles.
        maxiter -- integer; if convergence is not achieved in
                maxiter iterations, an error is raised. Must be >= 0.
        full_output -- bool (default: False), if True, also return
                object that contains information about convergence.

    EXAMPLES:
    An example involving an algebraic polynomial function.
        sage: R.<x> = QQ[]
        sage: f = (x+17)*(x-3)*(x-1/8)^3
        sage: find_root(f, 0,4)
        2.9999999999999951
        sage: find_root(f, 0,1)  # note -- precision of answer isn't very good on some machines.
        0.124999...
        sage: find_root(f, -20,-10)
        -17.0

    In Pomerance book on primes he asserts that the famous Riemann
    Hypothesis is equivalent to the statement that the function f(x)
    defined below is positive for all $x\geq 2.01$
        sage: def f(x):
        ...       return sqrt(x) * log(x) - abs(Li(x) - prime_pi(x))

    We find where $f$ equals, i.e., what value that is slightly smaller
    than $2.01$ that could have been used in the formulation of the Riemann
    Hypothesis:
        sage: find_root(f, 2, 4, rtol=0.0001)
        2.0082590205656166

    This agrees with the plot:
        sage: plot(f,2,2.01)
    """
    try:
        return f.find_root(a=a,b=b,xtol=xtol,rtol=rtol,maxiter=maxiter,full_output=full_output)
    except AttributeError:
        pass
    a = float(a); b = float(b)
    if a > b:
        a, b = b, a
    left = f(a)
    right = f(b)
    if left > 0 and right > 0:
        # Refine further -- try to find a point where this
        # function is negative in the interval
        val, s = find_minimum_on_interval(f, a, b)
        if val > 0:
            if val < rtol:
                if full_output:
                    return s, "No extra data"
                else:
                    return s
            raise RuntimeError, "f appears to have no zero on the interval"
        # If we found such an s, then we just instead find
        # a root between left and s or s and right.
        a = s   # arbitrary choice -- maybe should try both and take one that works?

    elif left < 0 and right < 0:
        # Refine further
        val, s = find_maximum_on_interval(f, a, b)
        if val < 0:
            if abs(val) < rtol:
                if full_output:
                    return s, "No extra data"
                else:
                    return s
            raise RuntimeError, "f appears to have no zero on the interval"
        a = s

    import scipy.optimize
    return scipy.optimize.brentq(f, a, b,
                                 full_output=full_output, xtol=xtol, rtol=rtol, maxiter=maxiter)

def find_maximum_on_interval(f, a, b, tol=1.48e-08, maxfun=500):
    """
    Numerically find the maximum of the expression f on the interval
    [a,b] (or [b,a]) along with the point at which the maximum is attained.

    See the documentation for \code{find_minimum_on_interval}
    for more details.

    EXAMPLES:
        sage: f = lambda x: x*cos(x)
        sage: find_maximum_on_interval(f, 0,5)
        (0.561096338191..., 0.8603335890...)
        sage: find_maximum_on_interval(f, 0,5, tol=0.1, maxfun=10)
        (0.561090323458..., 0.857926501456...)
    """
    def g(z):
        r"""
        Returns the negative of the input function f. Finding the maximum
        of f(z) on [a,b] is equivalent to finding th minimum of -f(z) on
        [a,b].

        EXAMPLES:

        """
        return -f(z)
    minval, x = find_minimum_on_interval(g, a=a, b=b, tol=tol, maxfun=maxfun)
    return -minval, x

def find_minimum_on_interval(f, a, b, tol=1.48e-08, maxfun=500):
    """
    Numerically find the minimum of the expression self on the
    interval [a,b] (or [b,a]) and the point at which it attains that
    minimum.  Note that self must be a function of (at most) one
    variable.

    INPUT:
        a,b -- endpoints of interval on which to minimize self.
        tol -- the convergence tolerance
        maxfun -- maximum function evaluations

    OUTPUT:
        minval -- (float) the minimum value that self takes on in the interval [a,b]
        x -- (float) the point at which self takes on the minimum value

    EXAMPLES:
        sage: f = lambda x: x*cos(x)
        sage: find_minimum_on_interval(f, 1, 5)
        (-3.28837139559..., 3.4256184695...)
        sage: find_minimum_on_interval(f, 1, 5, tol=1e-3)
        (-3.28837136189098..., 3.42575079030572...)
        sage: find_minimum_on_interval(f, 1, 5, tol=1e-2, maxfun=10)
        (-3.28837084598..., 3.4250840220...)
        sage: show(plot(f, 0, 20))
        sage: find_minimum_on_interval(f, 1, 15)
        (-9.4772942594..., 9.5293344109...)

    ALGORITHM: Uses scipy.optimize.fminbound which uses Brent's method.

    AUTHOR:
         -- William Stein (2007-12-07)
    """
    try:
        return f.find_minimum_on_interval(a=a, b=b, tol=tol,maxfun=maxfun)
    except AttributeError:
        pass
    a = float(a); b = float(b)
    import scipy.optimize
    xmin, fval, iter, funcalls = scipy.optimize.fminbound(f, a, b, full_output=1, xtol=tol, maxfun=maxfun)
    return fval, xmin


def minimize(func,x0,gradient=None,hessian=None,algorithm="default",**args):
    r"""
    This function is an interface to a variety of algorithms for computing the minimum of
    a function of several variables.

    INPUT:
       func - Either a symbolic function or a Python function whose argument is a tuple with n components

       x0 - Initial point for finding minimum.

       gradient - Optional gradient function. This will be computed automatically for symbolic functions.

                  For Python functions, it allows the use of algorithms requiring derivatives.  It should
                  accept a tuple of arguments and return a numpy array containing the partial derivatives
                  at that point.

       hessian -  Optional hessian function. This will be computed automatically for symbolic functions.

                  For Python functions, it allows the use of algorithms requiring derivatives. It should accept
                  a tuple of arguments and return a numpy array containing the 2nd partial derivatives
                  of the function.

       algorithm - String specifying algorithm to use. Options are
                   "default" (for Python functions, the simplex method is the default)
                             (for symbolic functions bfgs is the default)
                   "simplex"
                   "powell"
                   "bfgs" - (broyden-fletcher-goldfarb-shannon) - requires gradient
                   "cg" - (conjugate-gradient) -requires gradient
                   "ncg" - (newton-conjugate gradient) - requires gradient and hessian


    EXAMPLES:
        sage: vars=var('x y z')
        sage: f=100*(y-x^2)^2+(1-x)^2+100*(z-y^2)^2+(1-y)^2
        sage: minimize(f,[.1,.3,.4],disp=0)
        (1.00..., 1.00..., 1.00...)

        sage: minimize(f,[.1,.3,.4],algorithm="ncg",disp=0)
        (0.9999999..., 0.999999..., 0.999999...)

      Same example with just Python functions:

        sage: def rosen(x): # The Rosenbrock function
        ...      return sum(100.0r*(x[1r:]-x[:-1r]**2.0r)**2.0r + (1r-x[:-1r])**2.0r)
        sage: minimize(rosen,[.1,.3,.4],disp=0)
        (1.00..., 1.00..., 1.00...)

      Same example with a pure Python function and a Python function to compute the gradient

        sage: def rosen(x): # The Rosenbrock function
        ...      return sum(100.0r*(x[1r:]-x[:-1r]**2.0r)**2.0r + (1r-x[:-1r])**2.0r)
        sage: import numpy
        sage: from numpy import zeros
        sage: def rosen_der(x):
        ...      xm = x[1r:-1r]
        ...      xm_m1 = x[:-2r]
        ...      xm_p1 = x[2r:]
        ...      der = zeros(x.shape,dtype=float)
        ...      der[1r:-1r] = 200r*(xm-xm_m1**2r) - 400r*(xm_p1 - xm**2r)*xm - 2r*(1r-xm)
        ...      der[0] = -400r*x[0r]*(x[1r]-x[0r]**2r) - 2r*(1r-x[0])
        ...      der[-1] = 200r*(x[-1r]-x[-2r]**2r)
        ...      return der
        sage: minimize(rosen,[.1,.3,.4],gradient=rosen_der,algorithm="bfgs",disp=0)
        (1.00...,  1.00..., 1.00...)

    """
    from sage.calculus.calculus import SymbolicExpression
    import scipy
    from scipy import optimize
    if isinstance(func,SymbolicExpression):
        var_list=func.variables()
        var_names=map(str,var_list)
        fast_f=func._fast_float_(*var_names)
        f=lambda p: fast_f(*p)
        gradient_list=func.gradient()
        fast_gradient_functions=[gradient_list[i]._fast_float_(*var_names)  for i in xrange(len(gradient_list))]
        gradient=lambda p: scipy.array([ a(*p) for a in fast_gradient_functions])
    else:
        f=func

    if algorithm=="default":
        if gradient==None:
            min=optimize.fmin(f,map(float,x0),**args)
        else:
            min= optimize.fmin_bfgs(f,map(float,x0),fprime=gradient,**args)
    else:
        if algorithm=="simplex":
            min= optimize.fmin(f,map(float,x0),**args)
        elif algorithm=="bfgs":
            min= optimize.fmin_bfgs(f,map(float,x0),fprime=gradient,**args)
        elif algorithm=="cg":
            min= optimize.fmin_cg(f,map(float,x0),fprime=gradient,**args)
        elif algorithm=="powell":
            min= optimize.fmin_powell(f,map(float,x0),**args)
        elif algorithm=="ncg":
            if isinstance(func,SymbolicExpression):
                hess=func.hessian()
                hess_fast= [ [a._fast_float_(*var_names) for a in row] for row in hess]
                hessian=lambda p: [[a(*p) for a in row] for row in hess_fast]
                hessian_p=lambda p,v: scipy.dot(scipy.array(hessian(p)),v)
                min= optimize.fmin_ncg(f,map(float,x0),fprime=gradient,fhess=hessian,fhess_p=hessian_p,**args)
    return vector(RDF,min)

def minimize_constrained(func,cons,x0,gradient=None,algorithm='default', **args):
    r"""
    Minimize a function with constraints.

    INPUT:
       func - Either a symbolic function, or a Python function whose argument is a tuple with n components

       cons - constraints. This should be either a function or list of functions that must be positive.
              Alternatively, the constraints can be specified as a list of intervals that define the
              region we are minimizing in.

              If the constraints are specified as functions, the functions should be functions of a tuple
              with n components (assuming n variables).

              If the constraints are specifed as a list of intervals and there are no constraints for a given
              variable, that component can be (None,None).

       x0 - Initial point for finding minimum

       algorithm - Optional, specify the algorithm to use:
                   "default"  - default choices
                   "l-bfgs-b" - only effective, if you specify bound constraints
                                Ref.: C. Zhu, R. H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B,
                                      FORTRAN routines for large scale bound constrained optimization (1997),
                                      ACM Transactions on Mathematical Software, Vol 23, Num. 4, pp. 550 - 560.

       gradient - Optional gradient function. This will be computed automatically for symbolic functions.
                  This is only used when the constraints are specified as a list of intervals.

      EXAMPLES:
        Let us maximize $x+y-50$  subject to the following constraints: $50*x+24*y<=2400$,
        $30*x+33*y<=2100$, $x>=45$, and $y>=5$.

          sage: y = var('y')
          sage: f = lambda p: -p[0]-p[1]+50
          sage: c_1 = lambda p: p[0]-45
          sage: c_2 = lambda p: p[1]-5
          sage: c_3 = lambda p: -50*p[0]-24*p[1]+2400
          sage: c_4 = lambda p: -30*p[0]-33*p[1]+2100
          sage: a = minimize_constrained(f,[c_1,c_2,c_3,c_4],[2,3])
          sage: a
          (45.0, 6.25)

        Let's find a minimum of sin(x*y):

          sage: x,y = var('x y')
          sage: f = sin(x*y)
          sage: minimize_constrained(f, [(None,None),(4,10)],[5,5])
          (4.8..., 4.8...)

        Check, if L-BFGS-B finds the same minimum:
          sage: minimize_constrained(f, [(None,None),(4,10)],[5,5], algorithm='l-bfgs-b')
          (4.7..., 4.9...)

        Rosenbrock funcion, [http://en.wikipedia.org/wiki/Rosenbrock_function]:
          sage: from scipy.optimize import rosen, rosen_der
          sage: minimize_constrained(rosen, [(-50,-10),(5,10)],[1,1],gradient=rosen_der,algorithm='l-bfgs-b')
          (-10.0, 10.0)
          sage: minimize_constrained(rosen, [(-50,-10),(5,10)],[1,1],algorithm='l-bfgs-b')
          (-10.0, 10.0)

    """
    from sage.calculus.calculus import SymbolicExpression
    import scipy
    from scipy import optimize
    function_type=type(lambda x,y: x+y)

    if isinstance(func,SymbolicExpression):
        var_list=func.variables()
        var_names=map(str,var_list)
        fast_f=func._fast_float_(*var_names)
        f=lambda p: fast_f(*p)
        gradient_list=func.gradient()
        fast_gradient_functions=[gradient_list[i]._fast_float_(*var_names)  for i in xrange(len(gradient_list))]
        gradient=lambda p: scipy.array([ a(*p) for a in fast_gradient_functions])
    else:
        f=func

    if isinstance(cons,list):
        if isinstance(cons[0],tuple) or isinstance(cons[0],list) or cons[0]==None:
            if gradient!=None:
                if algorithm=='l-bfgs-b':
                    min= optimize.fmin_l_bfgs_b(f,x0,gradient,bounds=cons, iprint=-1, **args)[0]
                else:
                    min= optimize.fmin_tnc(f,x0,gradient,bounds=cons,messages=0,**args)[0]
            else:
                if algorithm=='l-bfgs-b':
                    min= optimize.fmin_l_bfgs_b(f,x0,approx_grad=True,bounds=cons,iprint=-1, **args)[0]
                else:
                    min= optimize.fmin_tnc(f,x0,approx_grad=True,bounds=cons,messages=0,**args)[0]

        elif isinstance(cons[0],function_type):
            min= optimize.fmin_cobyla(f,x0,cons,iprint=0,**args)
    elif isinstance(cons, function_type):
        min= optimize.fmin_cobyla(f,x0,cons,iprint=0,**args)
    return vector(RDF,min)


def linear_program(c,G,h,A=None,b=None):
    """
     Solves the dual linear programs:

     \begin{itemize}
       \item Minimize  $c'*x$ subject to $G*x + s = h$, $A*x=b$, and $s>=0$ where
       ' denotes transpose.
       \item Maximize  $-h'*z - b'*y$ subject to $G'*z + A'*y + c = 0$ and $z>=0$.
     \end{itemize}

    INPUT:
        c - a vector
        G - a matrix
        h - a vector
        A - a matrix
        b - a vector

        These can be over any field that can be turned into a floating point number.

    OUTPUT:
        A dictionary sol with keys x,s,y,z corresponding to the variables above

        sol['x'] - the solution to the linear program
        sol['s'] - the slack variables for the solution
        sol['z'] , sol['y'] - solutions to the dual program

    EXAMPLES:
      First, we minimize $-4x_1 - 5x_2$ subject to $2x_1 + x_2 <=3$,  $x_1 +  2x_2 <=3$,
      $x_1 >= 0$, and $x_2 >= 0$.
        sage: c=vector(RDF,[-4,-5])
        sage: G=matrix(RDF,[[2,1],[1,2],[-1,0],[0,-1]])
        sage: h=vector(RDF,[3,3,0,0])
        sage: sol=linear_program(c,G,h)
        sage: sol['x']
        (0.999..., 1.000...)

      Next, we maximize $x+y-50$ subject to $50*x+24*y<=2400$, $30*x+33*y<=2100$, $x>=45$,
      and $y>=5$.
        sage: v=vector([-1.0,-1.0,-1.0])
        sage: m=matrix([[50.0,24.0,0.0],[30.0,33.0,0.0],[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0],[0.0,0.0,-1.0]])
        sage: h=vector([2400.0,2100.0,-45.0,-5.0,1.0,-1.0])
        sage: sol=linear_program(v,m,h)
        sage: sol['x']
        (45.000000..., 6.2499999...3, 1.00000000...)

    """
    from cvxopt.base import matrix as m
    from cvxopt import solvers
    solvers.options['show_progress']=False
    c_=m(c.base_extend(RDF).numpy())
    G_=m(G.base_extend(RDF).numpy())
    h_=m(h.base_extend(RDF).numpy())
    if A!=None and b!=None:
        A_=m(A.base_extend(RDF).numpy())
        b_=m(b.base_extend(RDF).numpy())
        sol=solvers.lp(c_,G_,h_,A_,b_)
    else:
        sol=solvers.lp(c_,G_,h_)
    x=vector(RDF,list(sol['x']))
    s=vector(RDF,list(sol['s']))
    y=vector(RDF,list(sol['y']))
    z=vector(RDF,list(sol['z']))
    return  {'x':x,'s':s,'y':y,'z':z}


def find_fit(data, model, initial_guess = None, parameters = None, variables = None, solution_dict = False):
    r"""
    Finds numerical estimates for the parameters of the function model to give a best fit to data.

    INPUT:
       data - A two dimensional table of floating point numbers of the form
              [[x_{1,1}, x_{1,2}, \ldots, x_{1,k}, f_1],
               [x_{2,1}, x_{2,2}, \ldots, x_{2,k}, f_2],
               \ldots,
               [x_{n,1}, x_{n,2}, \ldots, x_{n,k}, f_n]]
              given as either a list of lists, matrix, or numpy array.

       model - Either a symbolic expression, symbolic function, or a Python function.
               model has to be a function of the variables (x_1, x_2, \ldots, x_k) and
               free parameters (a_1, a_2, \ldots, a_l).

       initial_guess - (default: None) Initial estimate for the parameters (a_1, a_2, \ldots, a_l),
                       given as either a list, tuple, vector or numpy array.
                       If None the default estimate for each parameter is 1.

       parameters - (default: None) A list of the parameters (a_1, a_2, \ldots, a_l),
                    If model is a symbolic function it is ignored, and the free parameters
                    of the symbolic function are used.

       variables - (default: None) A list of the variables (x_1, x_2, \ldots, x_k).
                   If model is a symbolic function it is ignored, and the variables of
                   the symbolic function are used.

       solution_dict - (default: False) if True, return the solution
                       as a dictionary rather than an equation.

      EXAMPLES:
        First we create some datapoints of a sine function with some random perturbations:
          sage: data = [(i, 1.2 * sin(0.5*i-0.2) + 0.1 * normalvariate(0, 1)) for i in xsrange(0, 4*pi, 0.2)]
          sage: var('a, b, c, x')
          (a, b, c, x)

        We define a function with free parameters $a$, $b$ and $c$:
          sage: model(x) = a * sin(b * x - c)

        We search for the parameters that give the best fit to the data:
          sage: find_fit(data, model)
          [a == 1.21..., b == 0.49..., c == 0.19...]

        We can also use a python function for the model:
          sage: def f(x, a, b, c): return a * sin(b * x - c)
          sage: find_fit(data, f, parameters = [a, b, c], variables = [x], solution_dict = True)
          {a: 1.21..., c: 0.19..., b: 0.49...}

        We search for a formula for the n-th prime number:
          sage: dataprime = [(i, nth_prime(i)) for i in xrange(1, 5000, 100)]
          sage: find_fit(dataprime, a * x * log(b * x), parameters = [a, b], variables = [x])
          [a == 1.11..., b == 1.24...]

      ALGORITHM: Uses scipy.optimize.leastsq which in turn uses MINPACK's lmdif and lmder algorithms.
    """
    import numpy

    if not isinstance(data, numpy.ndarray):
        try:
            data = numpy.array(data, dtype = float)
        except (ValueError, TypeError):
            raise TypeError, "data has to be a list of lists, a matrix, or a numpy array"
    elif data.dtype == object:
        raise ValueError, "the entries of data have to be of type float"

    if data.ndim != 2:
        raise ValueError, "data has to be a two dimensional table of floating point numbers"

    from sage.calculus.calculus import SymbolicExpression, CallableSymbolicExpression
    if isinstance(model, CallableSymbolicExpression):
        parameters = list(model.variables())
        variables = list(model.arguments())
        for v in variables:
            parameters.remove(v)

    if data.shape[1] != len(variables) + 1:
        raise ValueError, "each row of data needs %d entries, only %d entries given" % (len(variables) + 1, data.shape[1])

    if parameters == None or len(parameters) == 0 or \
       variables == None or len(variables) == 0:
        raise ValueError, "no variables given"

    if initial_guess == None:
        initial_guess = len(parameters) * [1]

    if not isinstance(initial_guess, numpy.ndarray):
        try:
            initial_guess = numpy.array(initial_guess, dtype = float)
        except (ValueError, TypeError):
            raise TypeError, "initial_guess has to be a list, tuple, or numpy array"
    elif initial_guess.dtype == object:
        raise ValueError, "the entries of initial_guess have to be of type float"

    if len(initial_guess) != len(parameters):
        raise ValueError, "length of initial_guess does not coincide with the number of parameters"

    if isinstance(model, SymbolicExpression):
        var_list = variables + parameters
        var_names = map(str, var_list)
        func = model._fast_float_(*var_names)
    else:
        func = model

    def function(x_data, params):
        result = numpy.zeros(len(x_data))
        for row in xrange(len(x_data)):
            fparams = numpy.hstack((x_data[row], params)).tolist()
            result[row] = func(*fparams)
        return result

    def error_function(params, x_data, y_data):
        result = numpy.zeros(len(x_data))
        for row in xrange(len(x_data)):
            fparams = x_data[row].tolist() + params.tolist()
            result[row] = func(*fparams)
        return result - y_data

    x_data = data[:, 0:len(variables)]
    y_data = data[:, -1]

    from scipy.optimize import leastsq
    estimated_params, d = leastsq(error_function, initial_guess, args = (x_data, y_data))

    if isinstance(estimated_params, float):
        estimated_params = [estimated_params]
    else:
        estimated_params = estimated_params.tolist()

    if solution_dict:
       dict = {}
       for item in zip(parameters, estimated_params):
           dict[item[0]] = item[1]
       return dict

    return [item[0] == item[1] for item in zip(parameters, estimated_params)]

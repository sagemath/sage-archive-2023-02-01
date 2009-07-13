********
Calculus
********

Here are some examples of calculus symbolic computations using
Sage. They use the Maxima interface.

Work is being done to make the commands for the symbolic
calculations given below more intuitive and natural. At the moment,
we use the maxima class interface.


.. index::
   pair: calculus; differentiation

Differentiation
===============

Differentiation:

::

    sage: var('x k w')
    (x, k, w)
    sage: f = x^3 * e^(k*x) * sin(w*x); f
    x^3*e^(k*x)*sin(w*x)
    sage: f.diff(x)
    k*x^3*e^(k*x)*sin(w*x) + w*x^3*e^(k*x)*cos(w*x) + 3*x^2*e^(k*x)*sin(w*x)
    sage: latex(f.diff(x))
    k x^{3} e^{k x} \sin\left(w x\right) + w x^{3} e^{k x} \cos\left(w x\right) + 3 \, x^{2} e^{k x} \sin\left(w x\right)

If you type ``view(f.diff('x'))`` another window will open up
displaying the compiled output. In the notebook, you can enter

::

    f = maxima('x^3 * %e^(k*x) * sin(w*x)')
    show(f)
    show(f.diff('x'))

into a cell and press ``shift-enter`` for a similar result. You can
also call Maxima indirectly using the commands

::

    R = PolynomialRing(QQ,"x")
    x = R.gen()
    p = x^2 + 1
    show(p.derivative())
    show(p.integral())

in a notebook cell, or

::

    sage: R = PolynomialRing(QQ,"x")
    sage: x = R.gen()
    sage: p = x^2 + 1
    sage: p.derivative()
    2*x
    sage: p.integral()
    1/3*x^3 + x

on the command line.  At this point you can also type
``view(p.derivative())`` or ``view(p.integral())`` to open a new
window with output typeset by LaTeX.

.. index::
   pair: calculus; critical points

Critical points
---------------

You can find critical points of a piecewise defined function:

::

    sage: x = PolynomialRing(RationalField(), 'x').gen()
    sage: f1 = x^0
    sage: f2 = 1-x
    sage: f3 = 2*x
    sage: f4 = 10*x-x^2
    sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
    sage: f.critical_points()
    [5.0]

.. index:: Taylor series, power series

Power series
------------

Taylor series:

::

    sage: var('f0 k x')
    (f0, k, x)
    sage: g = f0/sinh(k*x)^4
    sage: g.taylor(x, 0, 3)
    -62/945*f0*k^2*x^2 + 11/45*f0 - 2/3*f0/(k^2*x^2) + f0/(k^4*x^4)
    sage: maxima(g).powerseries('x',0)
    16*f0*('sum((2^(2*i1-1)-1)*bern(2*i1)*k^(2*i1-1)*x^(2*i1-1)/(2*i1)!,i1,0,inf))^4

Of course, you can view the LaTeX-ed version of this using
``view(g.powerseries('x',0))``.

The Maclaurin and power series of
:math:`\log({\frac{\sin(x)}{x}})`:

::

    sage: f = log(sin(x)/x)
    sage: f.taylor(x, 0, 10)
    -1/467775*x^10 - 1/37800*x^8 - 1/2835*x^6 - 1/180*x^4 - 1/6*x^2
    sage: [bernoulli(2*i) for i in range(1,7)]
    [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730]
    sage: maxima(f).powerseries(x,0)
    ('sum((-1)^i2*2^(2*i2)*bern(2*i2)*x^(2*i2)/(i2*(2*i2)!),i2,1,inf))/2

.. index::
   pair: calculus; integration

Integration
===========

Numerical integration is discussed in  :ref:`section-riemannsums` below.

Sage can integrate some simple functions on its own:

::

    sage: f = x^3
    sage: f.integral()
    1/4*x^4
    sage: integral(x^3,x)
    1/4*x^4
    sage: f = x*sin(x^2)
    sage: integral(f,x)
    -1/2*cos(x^2)

Sage can also compute symbolic definite integrals involving limits.

::

    sage: var('x, k, w')
    (x, k, w)
    sage: f = x^3 * e^(k*x) * sin(w*x)
    sage: f.integrate(x)
    -(((k^6*w + 3*k^4*w^3 + 3*k^2*w^5 + w^7)*x^3 - 24*k^3*w + 24*k*w^3 - 6*(k^5*w + 2*k^3*w^3 + k*w^5)*x^2 + 6*(3*k^4*w + 2*k^2*w^3 - w^5)*x)*e^(k*x)*cos(w*x) - ((k^7 + 3*k^5*w^2 + 3*k^3*w^4 + k*w^6)*x^3 - 6*k^4 + 36*k^2*w^2 - 6*w^4 - 3*(k^6 + k^4*w^2 - k^2*w^4 - w^6)*x^2 + 6*(k^5 - 2*k^3*w^2 - 3*k*w^4)*x)*e^(k*x)*sin(w*x))/(k^8 + 4*k^6*w^2 + 6*k^4*w^4 + 4*k^2*w^6 + w^8)
    sage: integrate(1/x^2, x, 1, infinity)
    1


.. index: convolution

Convolution
-----------

You can find the convolution of any piecewise defined function with
another (off the domain of definition, they are assumed to be
zero). Here is :math:`f`, :math:`f*f`, and :math:`f*f*f`,
where :math:`f(x)=1`, :math:`0<x<1`:

::

    sage: x = PolynomialRing(QQ, 'x').gen()
    sage: f = Piecewise([[(0,1),1*x^0]])
    sage: g = f.convolution(f)
    sage: h = f.convolution(g)
    sage: P = f.plot(); Q = g.plot(rgbcolor=(1,1,0)); R = h.plot(rgbcolor=(0,1,1))

To view this, type ``show(P+Q+R)``.


.. _section-riemannsums:

Riemann and trapezoid sums for integrals
----------------------------------------

Regarding numerical approximation of :math:`\int_a^bf(x)\, dx`,
where :math:`f` is a piecewise defined function, can


-  compute (for plotting purposes) the piecewise linear function
   defined by the trapezoid rule for numerical integration based on a
   subdivision into :math:`N` subintervals

-  the approximation given by the trapezoid rule,

-  compute (for plotting purposes) the piecewise constant function
   defined by the Riemann sums (left-hand, right-hand, or midpoint) in
   numerical integration based on a subdivision into :math:`N`
   subintervals,

-  the approximation given by the Riemann sum approximation.


::

    sage: f1(x) = x^2
    sage: f2(x) = 5-x^2
    sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
    sage: f.trapezoid(4)
    Piecewise defined function with 4 parts, [[(0, 1/2), 1/2*x],
    [(1/2, 1), 9/2*x - 2], [(1, 3/2), 1/2*x + 2],
    [(3/2, 2), -7/2*x + 8]]
    sage: f.riemann_sum_integral_approximation(6,mode="right")
    19/6
    sage: f.integral()
    Piecewise defined function with 2 parts,
    [[(0, 1), x |--> 1/3*x^3], [(1, 2), x |--> -1/3*x^3 + 5*x - 13/3]]
    sage: f.integral(definite=True)
    3

.. index: Laplace transform

Laplace transforms
------------------

If you have a piecewise-defined polynomial function then there is a
"native" command for computing Laplace transforms. This calls
Maxima but it's worth noting that Maxima cannot handle (using the
direct interface illustrated in the last few examples) this type of
computation.

::

    sage: var('x s')
    (x, s)
    sage: f1(x) = 1
    sage: f2(x) = 1-x
    sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
    sage: f.laplace(x, s)
    (s + 1)*e^(-2*s)/s^2 - e^(-s)/s + 1/s - e^(-s)/s^2

For other "reasonable" functions, Laplace transforms can be
computed using the Maxima interface:

::

    sage: var('k, s, t')
    (k, s, t)
    sage: f = 1/exp(k*t)
    sage: f.laplace(t,s)
    1/(k + s)

is one way to compute LT's and

::

    sage: var('s, t')
    (s, t)
    sage: f = t^5*exp(t)*sin(t)
    sage: L = laplace(f, t, s); L
    3840*(s - 1)^5/(s^2 - 2*s + 2)^6 - 3840*(s - 1)^3/(s^2 - 2*s + 2)^5 +
    720*(s - 1)/(s^2 - 2*s + 2)^4

is another way.

.. index:
   pair: differential equations; solve

Ordinary differential equations
===============================

Symbolically solving ODEs can be done using 's interface with
Maxima. Numerical solution of ODEs can be done using 's interface
with Octave (an experimental package), or routines in the GSL (Gnu
Scientific Library).

You can solve ODE's symbolically in Sage using the Maxima interface
(do not type the ``...``):

::

    sage: maxima.de_solve('diff(y,x,2) + 3*x = y', ['x','y'], [1,1,1])
    y=3*x-2*%e^(x-1)
    sage: maxima.de_solve('diff(y,x,2) + 3*x = y', ['x','y'])
    y=%k1*%e^x+%k2*%e^-x+3*x
    sage: maxima.de_solve('diff(y,x) + 3*x = y', ['x','y'])
    y=(%c-3*(-x-1)*%e^-x)*%e^x
    sage: maxima.de_solve('diff(y,x) + 3*x = y', ['x','y'],[1,1])
    y=-%e^-1*(5*%e^x-3*%e*x-3*%e)

    sage: maxima.clear('x'); maxima.clear('f')
    sage: maxima.de_solve_laplace("diff(f(x),x,2) = 2*diff(f(x),x)-f(x)",\
    ...   ["x","f"], [0,1,2])
    f(x)=x*%e^x+%e^x

    sage: maxima.clear('x'); maxima.clear('f')
    sage: f = maxima.de_solve_laplace("diff(f(x),x,2) = 2*diff(f(x),x)-f(x)",\
    ...   ["x","f"])
    sage: f
    f(x)=x*%e^x*('at('diff(f(x),x,1),x=0))-f(0)*x*%e^x+f(0)*%e^x
    sage: f.display2d()
                                                   !
                                       x  d        !                  x          x
                            f(x) = x %e  (-- (f(x))!     ) - f(0) x %e  + f(0) %e
                                          dx       !
                                                   !x = 0

.. index:
   pair: differential equations; plot

If you have ``Octave`` and ``gnuplot`` installed,

::

    sage: octave.de_system_plot(['x+y','x-y'], [1,-1], [0,2]) # optional octave required

yields the two plots :math:`(t,x(t)), (t,y(t))` on the same graph
(the :math:`t`-axis is the horizonal axis) of the system of ODEs

.. math::
    x' = x+y, x(0) = 1; y' = x-y, y(0) = -1,

for :math:`0 <= t <= 2`. Another way this system can be solved is to use the
command ``desolve_system`` in ``calculus/examples``.

.. skip

::

    sage: attach os.environ['SAGE_ROOT'] + '/examples/calculus/desolvers.sage'
    sage: des = ["'diff(x(t),t)=x(t)+y(t)","'diff(y(t),t)=x(t)-y(t)"]
    sage: vars = ["t","x","y"]
    sage: ics = [0,1,-1]
    sage: desolve_system(des,vars,ics)
    [x(t)=cosh(2^(1/2)*t),y(t)=2*sinh(2^(1/2)*t)/sqrt(2)-cosh(2^(1/2)*t)]

The output of this command is *not* a pair of functions.

Finally, can solve linear DEs using power series:

::

    sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)
    sage: a = 2 - 3*t + 4*t^2 + O(t^10)
    sage: b = 3 - 4*t^2 + O(t^7)
    sage: f = a.solve_linear_de(prec=5, b=b, f0=3/5)
    sage: f
    3/5 + 21/5*t + 33/10*t^2 - 38/15*t^3 + 11/24*t^4 + O(t^5)
    sage: f.derivative() - a*f - b
    O(t^4)

Fourier series of periodic functions
====================================

If :math:`f(x)` is a piecewise-defined polynomial function on
:math:`-L<x<L` then the Fourier series

.. math::
   f(x) \sim \frac{a_0}{2} + \sum_{n=1}^\infty \left[a_n\cos\left(\frac{n\pi x}{L}\right) +
   b_n\sin\left(\frac{n\pi x}{L}\right)\right]


converges. In addition to computing the coefficients
:math:`a_n,b_n`, it will also compute the partial sums (as a
string), plot the partial sums (as a function of :math:`x` over
:math:`(-L,L)`, for comparison with the plot of :math:`f(x)`
itself), compute the value of the FS at a point, and similar
computations for the cosine series (if :math:`f(x)` is even) and
the sine series (if :math:`f(x)` is odd). Also, it will plot the
partial F.S. Cesaro mean sums (a "smoother" partial sum
illustrating how the Gibbs phenomenon is mollified).

::

    sage: f1 = lambda x: -1
    sage: f2 = lambda x: 2
    sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
    sage: f.fourier_series_cosine_coefficient(5,pi)
    -3/5/pi
    sage: f.fourier_series_sine_coefficient(2,pi)
    -3/pi
    sage: f.fourier_series_partial_sum(3,pi)
    -3*sin(2*x)/pi + sin(x)/pi - 3*cos(x)/pi + 1/4

Type ``show(f.plot_fourier_series_partial_sum(15,pi,-5,5))`` and
``show(f.plot_fourier_series_partial_sum_cesaro(15,pi,-5,5))``
(and be patient) to view the partial sums.

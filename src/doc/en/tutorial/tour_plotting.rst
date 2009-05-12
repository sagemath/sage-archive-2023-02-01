.. _section-plot:

Plotting
========

Sage can produce two-dimensional and three-dimensional plots.

Two-dimensional Plots
---------------------

In two dimensions, Sage can draw circles, lines, and polygons;
plots of functions in rectangular coordinates; and also polar
plots, contour plots and vector field plots. We present examples of
some of these here. For more examples of plotting with Sage, see
:ref:`section-systems` and :ref:`section-maxima`, and also the
"Sage Constructions" documentation.

This command produces a yellow circle of radius 1, centered at the
origin:

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0))

You can also produce a filled circle:

::

    sage: circle((0,0), 1, rgbcolor=(1,1,0), fill=True)

You can also create a circle by assigning it to a variable; this
does not plot it:

::

    sage: c = circle((0,0), 1, rgbcolor=(1,1,0))

To plot it, use ``c.show()`` or ``show(c)``, as follows:

.. link

::

    sage: c.show()

Alternatively, evaluating ``c.save('filename.png')`` will save the
plot to the given file.

Now, these 'circles' look more like ellipses because the axes are
scaled differently.  You can fix this:

.. link

::

    sage: c.show(aspect_ratio=1)

The command ``show(c, aspect_ratio=1)`` accomplishes the same
thing, or you can save the picture using
``c.save('filename.png', aspect_ratio=1)``.

It's easy to plot basic functions:

::

    sage: plot(cos, (-5,5))

Once you specify a variable name, you can create parametric plots
also:

::

    sage: x = var('x')
    sage: parametric_plot((cos(x),sin(x)^3),0,2*pi,rgbcolor=hue(0.6))

You can combine several plots by adding them:

::

    sage: x = var('x')
    sage: p1 = parametric_plot((cos(x),sin(x)),0,2*pi,rgbcolor=hue(0.2))
    sage: p2 = parametric_plot((cos(x),sin(x)^2),0,2*pi,rgbcolor=hue(0.4))
    sage: p3 = parametric_plot((cos(x),sin(x)^3),0,2*pi,rgbcolor=hue(0.6))
    sage: show(p1+p2+p3, axes=false)

A good way to produce filled-in shapes is to produce a list of
points (``L`` in the example below) and then use the ``polygon``
command to plot the shape with boundary formed by those points. For
example, here is a green deltoid:

::

    sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),\
    ...   2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,3/4,1/2))
    sage: p

Type ``show(p, axes=false)`` to see this without any axes.

You can add text to a plot:

::

    sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),\
    ...   6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
    sage: p = polygon(L, rgbcolor=(1/8,1/4,1/2))
    sage: t = text("hypotrochoid", (5,4), rgbcolor=(1,0,0))
    sage: show(p+t)

Calculus teachers draw the following plot frequently on the board:
not just one branch of arcsin but rather several of them: i.e., the
plot of :math:`y=\sin(x)` for :math:`x` between :math:`-2\pi`
and :math:`2\pi`, flipped about the 45 degree line. The following
Sage commands construct this:

::

    sage: v = [(sin(x),x) for x in srange(-2*float(pi),2*float(pi),0.1)]
    sage: line(v)

Since the tangent function has a larger range than sine, if you use
the same trick to plot the inverse tangent, you should change the
minimum and maximum coordinates for the *x*-axis:

::

    sage: v = [(tan(x),x) for x in srange(-2*float(pi),2*float(pi),0.01)]
    sage: show(line(v), xmin=-20, xmax=20)

Sage also computes polar plots, contour plots and vector field plots
(for special types of functions). Here is an example of a contour
plot:

::

    sage: f = lambda x,y: cos(x*y)
    sage: contour_plot(f, (-4, 4), (-4, 4))

Three-Dimensional Plots
-----------------------

Sage produces three-dimensional plots using an open source package
called [Jmol]_. Here are a few examples:

Yellow Whitney's umbrella
http://en.wikipedia.org/wiki/Whitney_umbrella:

::

    sage: u, v = var('u,v')
    sage: fx = u*v
    sage: fy = u
    sage: fz = v^2
    sage: parametric_plot3d([fx, fy, fz], (u, -1, 1), (v, -1, 1),
    ...   frame=False, color="yellow")


Once you have evaluated ``parametric_plot3d``, so that the plot is visible,
you can click and drag on it to rotate the figure.

Cross cap http://en.wikipedia.org/wiki/Cross-cap:

::

    sage: u, v = var('u,v')
    sage: fx = (1+cos(v))*cos(u)
    sage: fy = (1+cos(v))*sin(u)
    sage: fz = -tanh((2/3)*(u-pi))*sin(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ...   frame=False, color="red")

Twisted torus:

::

    sage: u, v = var('u,v')
    sage: fx = (3+sin(v)+cos(u))*cos(2*v)
    sage: fy = (3+sin(v)+cos(u))*sin(2*v)
    sage: fz = sin(u)+2*cos(v)
    sage: parametric_plot3d([fx, fy, fz], (u, 0, 2*pi), (v, 0, 2*pi),
    ...   frame=False, color="red")


.. [Jmol] Jmol: an open-source Java viewer for chemical structures in 3D http://www.jmol.org/
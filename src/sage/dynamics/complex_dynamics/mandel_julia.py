r"""
Mandelbrot and Julia sets

Plots the Mandelbrot set for the map $f(z)=z^2+c$ in the complex plane.

AUTHORS:

- Ben Barros

"""

#*****************************************************************************
#       Copyright (C) 2017 BEN BARROS <bbarros@slu.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def mandelbrot_plot(**kwds):
    r"""
    Interactive plot of the Mandelbrot set for the map $f(z) = z^2 + c$.

    kwds:

    - ``x_center`` -- float (optional - default: ``-1.0``), Real part of center point.

    - ``y_center`` -- float (optional - default: ``0.0``), Imaginary part of center point.

    - ``image_width`` -- float (optional - default: ``4.0``), width of image in the complex plane.

    - ``max_iteration`` -- int (optional - default: ``500``), maximum number of iterations the map $f(z)$.

    - ``pixel_count`` -- int (optional - default: ``500``), side length of image in number of pixels.

    - ``base_color`` -- RGB color (optional - default: ``[70,40,240]``) color used to determine the coloring of set.

    - ``interacts`` -- boolean (optional - default: ``True``), controls whether plot will have interactive functionality.
        For a static plot of the Mandelbrot set, set ``interacts`` to ``False``.

    OUTPUT:

    Plot of the Mandelbrot set.

    NOTEBOOK EXAMPLES:

    Interactive plot of the Mandelbrot set::

        mandelbrot_plot()

    If we would like to adjust the size of the image, change the value of ``pixel_count``::

        mandelbrot_plot(pixel_count= 1000)

    To change the coloring of the set, change the value of ``base_color``::

        mandelbrot_plot(base_color= [200,90,50])

    We can also change the initial position of our interactive plot using key words::

        mandelbrot_plot(x_center= -0.75, y_center= 0.25, image_width= 1/2)

    To plot a static version of the Mandelbrot set, set ``interacts`` to ``False``::

        mandelbrot_plot(interacts= False)

    We can also increase the number of iterations by changing ``max_iteration``. This is helpful when zooming
    into small portions of the set::

        mandelbrot_plot(interacts= False, x_center= -1.1, y_center= 0.23, image_width= 1/4096, max_iteration= 10000)
    """

    from sage.dynamics.complex_dynamics.mandel_julia_helper import _fast_mandel_plot
    from sagenb.notebook.interact import interact
    from sagenb.notebook.interact import slider
    from sagenb.notebook.interact import input_box


    x_center = kwds.pop("x_center",-1.0)
    y_center = kwds.pop("y_center",0.0)
    image_width = kwds.pop("image_width",4.0)
    max_iteration = kwds.pop("max_iteration",500)
    pixel_count = kwds.pop("pixel_count",500)
    base_color = kwds.pop("base_color",[70,40,240])
    interacts = kwds.pop("interacts",True)

    if interacts:
        @interact(layout={'bottom':[['real_center'],['im_center'],['width']],'top':[['iterations']]})
        def _(real_center = input_box(x_center, 'Real'),
            im_center = input_box(y_center,'Imaginary'),
            width = slider([2**(-i) for i in range(-2,15)],label= 'Width of Image',default= image_width),
            iterations = input_box(max_iteration,'Max Number of Iterations')):
            print "Center: %s + %s*i" % (real_center,im_center)
            _fast_mandel_plot(real_center, im_center, width, iterations, pixel_count, base_color).show()

    else:
        _fast_mandel_plot(x_center, y_center, image_width, max_iteration, pixel_count, base_color).show()

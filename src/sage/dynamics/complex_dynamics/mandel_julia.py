r"""
Mandelbrot and Julia sets

Plots the Mandelbrot and Julia sets for general polynomial maps in the complex
plane.

The Mandelbrot set is the set of complex numbers `c` for which the map
`f_c(z)` does not diverge when iterated from `z = 0`. This set of complex
numbers can be visualized by plotting each value for `c` in the complex plane.
The Mandelbrot set is often an example of a fractal when plotted in the complex
plane.

The Julia set for a given parameter `c` is the set of complex numbers for which
the function `f_c(z)` is bounded under iteration.

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

from __future__ import absolute_import, division
from sage.dynamics.complex_dynamics.mandel_julia_helper import (fast_mandelbrot_plot,
                                                                fast_external_ray,
                                                                convert_to_pixels,
                                                                get_line,
                                                                fast_julia_plot,
                                                                polynomial_mandelbrot,
                                                                julia_helper)
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from sagenb.notebook.interact import (interact,
                                      slider,
                                      input_box,
                                      color_selector,
                                      checkbox)
from sage.plot.colors import Color
from sage.repl.image import Image
from sage.functions.log import logb
from sage.rings.all import QQ, CC
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.misc.prandom import randint
from sage.calculus.var import var

EPS = 0.00001


def mandelbrot_plot(f=None, **kwds):
    r"""
    Plot of the Mandelbrot set for a general polynomial map `f_c(z)`.

    REFERENCE:

    [Dev2005]_

    INPUT:


    - ``f`` -- map (optional - default: ``z^2 + c``), polynomial map used to
    plot the Mandelbrot set.

    - ``parameter`` -- variable (optional - default: ``c``), parameter variable
    used to plot the Mandelbrot set.

    - ``x_center`` -- double (optional - default: ``-1.0``), Real part of center
    point.

    - ``y_center`` -- double (optional - default: ``0.0``), Imaginary part of
    center point.

    - ``image_width`` -- double (optional - default: ``4.0``), width of image
    in the complex plane.

    - ``max_iteration`` -- long (optional - default: ``500``), maximum number of
    iterations the map ``f_c(z)``.

    - ``pixel_count`` -- long (optional - default: ``500``), side length of
    image in number of pixels.

    - ``base_color`` -- Hex color (optional - default: ``tomato``) color
    used to determine the coloring of set.

    - ``level_sep`` -- long (optional - default: 1) number of iterations
    between each color level.

    - ``number_of_colors`` -- long (optional - default: 30) number of colors
    used to plot image.

    - ``interact`` -- boolean (optional - default: ``False``), controls whether
    plot will have interactive functionality.

    OUTPUT:

    24-bit RGB image of the Mandelbrot set in the complex plane.

    EXAMPLES:

    ::

        sage: mandelbrot_plot()
        500x500px 24-bit RGB image

    ::

        sage: mandelbrot_plot(pixel_count=1000)
        1000x1000px 24-bit RGB image

    ::

        sage: mandelbrot_plot(x_center=-1.11, y_center=0.2283, image_width=1/128, # long time
        ....: max_iteration=2000, number_of_colors=500, base_color=[40, 100, 100])
        500x500px 24-bit RGB image

    To display an interactive plot of the Mandelbrot in the Notebook, set
    ``interact`` to ``True``. (This is only implemented for ``z^2 + c``)::

        sage: mandelbrot_plot(interact=True)
        <html>...</html>

    ::

        sage: mandelbrot_plot(interact=True, x_center=-0.75, y_center=0.25,
        ....: image_width=1/2, number_of_colors=75)
        <html>...</html>

    Polynomial maps can be defined over a multivariate polynomial ring or a
    univariate polynomial ring tower::

        sage: R.<z,c> = CC[]
        sage: f = z^2 + c
        sage: mandelbrot_plot(f)  # not tested
        500x500px 24-bit RGB image

    ::

        sage: B.<z> = CC[]
        sage: R.<c> = B[]
        sage: f = z^5 + c
        sage: mandelbrot_plot(f) # not tested
        500x500px 24-bit RGB image

    When the polynomial is defined over a multivariate polynomial ring it is
    necessary to specify the parameter variable (default parameter is ``c``)::

        sage: R.<a,b> = CC[]
        sage: f = a^2 + b^3
        sage: mandelbrot_plot(f, parameter=b) # not tested
        500x500px 24-bit RGB image

    Interact functionality is not implemented for general polynomial maps::

        sage: R.<z,c> = CC[]
        sage: f = z^3 + c
        sage: mandelbrot_plot(f, interact=True) # not tested
        NotImplementedError: Interact only implemented for z^2 + c
    """
    parameter = kwds.pop("parameter", None)
    x_center = kwds.pop("x_center", 0.0)
    y_center = kwds.pop("y_center", 0.0)
    image_width = kwds.pop("image_width", 4.0)
    max_iteration = kwds.pop("max_iteration", None)
    pixel_count = kwds.pop("pixel_count", 500)
    level_sep = kwds.pop("level_sep", 1)
    number_of_colors = kwds.pop("number_of_colors", 30)
    interacts = kwds.pop("interact", False)
    base_color = kwds.pop("base_color", Color('tomato'))
    # Check if user specified maximum number of iterations
    given_iterations = True
    if max_iteration is None:
        # Set default to 500 for z^2 + c map
        max_iteration = 500
        given_iterations = False

    from ipywidgets.widgets import FloatSlider, IntSlider, ColorPicker, interact
    widgets = dict(
                   x_center = FloatSlider(min=-1.0, max=1.0, step=EPS,
                                          value=x_center, description="Real center"),
                   y_center = FloatSlider(min=-1.0, max=1.0, step=EPS,
                                          value=y_center, description="Imag center"),
                   image_width = FloatSlider(min=EPS, max=4.0, step=EPS,
                                             value=image_width, description="Zoom"),
                   max_iteration = IntSlider(min=0, max=600,
                                             value=max_iteration, description="Iterations"),
                   pixel_count = IntSlider(min=10, max=600,
                                           value=pixel_count, description="Pixels"),
                   level_sep = IntSlider(min=1, max=20,
                                         value=level_sep, description="Color sep"),
                   color_num = IntSlider(min=1, max=100,
                                         value=number_of_colors, description="# Colors"),
                   base_color = ColorPicker(value=Color(base_color).html_color(),
                                            description="Base color"),
                   )
            
    if f is None:
        # Quadratic map f = z^2 + c
        
        if interacts:
            return interact(**widgets).widget(fast_mandelbrot_plot)
        
        else:
            return fast_mandelbrot_plot(x_center, y_center, image_width,
             max_iteration, pixel_count, level_sep, number_of_colors,
             base_color)

    else:
        if parameter is None:
            c = var('c')
            parameter = c

        P = f.parent()

        if P.base_ring() is CC or P.base_ring() is CDF:
            gen_list = list(P.gens())
            parameter = gen_list.pop(gen_list.index(parameter))
            variable = gen_list.pop()

        elif P.base_ring().base_ring() is CC or P.base_ring().base_ring() is CDF:
            parameter = P.gen()
            variable = P.base().gen()

        else:
            raise ValueError("Base ring must be a complex field")

        if f == variable**2 + parameter:
            # Quadratic map f = z^2 + c
            if interacts:
                return interact(**widgets).widget(fast_mandelbrot_plot)

            else:
                return fast_mandelbrot_plot(x_center, y_center, image_width,
                 max_iteration, pixel_count, level_sep, number_of_colors,
                 base_color)
        else:
            if interacts:
                raise NotImplementedError("Interact only implemented for z^2 + c")
            else:
                # Set default of max_iteration to 50 for general polynomial maps
                # This prevents the function from being very slow by default
                if not given_iterations:
                    max_iteration = 50

                # Mandelbrot of General Polynomial Map
                return polynomial_mandelbrot(f, parameter, x_center, y_center, \
                 image_width, max_iteration, pixel_count, level_sep, \
                 number_of_colors, base_color)

def external_ray(theta, **kwds):
    r"""
    Draws the external ray(s) of a given angle (or list of angles)
    by connecting a finite number of points that were approximated using
    Newton's method. The algorithm used is described in a paper by
    Tomoki Kawahira.

    REFERENCE:

    [Kaw2009]_

    INPUT:

    - ``theta`` -- double or list of doubles, angles between 0 and 1 inclusive.

    kwds:

    - ``image`` -- 24-bit RGB image (optional - default: None) user specified
     image of Mandelbrot set.

    - ``D`` -- long (optional - default: ``25``) depth of the approximation.
     As ``D`` increases, the external ray gets closer to the boundary of the
     Mandelbrot set. If the ray doesn't reach the boundary of the Mandelbrot
     set, increase ``D``.

    - ``S`` -- long (optional - default: ``10``) sharpness of the approximation.
     Adjusts the number of points used to approximate the external ray (number
     of points is equal to ``S*D``). If ray looks jagged, increase ``S``.

    - ``R`` -- long (optional - default: ``100``) radial parameter. If ``R`` is
     large, the external ray reaches sufficiently close to infinity. If ``R`` is
     too small, Newton's method may not converge to the correct ray.

    - ``prec`` -- long (optional - default: ``300``) specifies the bits of
     precision used by the Complex Field when using Newton's method to compute
     points on the external ray.

    - ``ray_color`` -- RGB color (optional - default: ``[255, 255, 255]``) color
     of the external ray(s).

    OUTPUT:

    24-bit RGB image of external ray(s) on the Mandelbrot set.

    EXAMPLES::

        sage: external_ray(1/3)
        500x500px 24-bit RGB image

    ::

        sage: external_ray(0.6, ray_color=[255, 0, 0])
        500x500px 24-bit RGB image

    ::

        sage: external_ray([0, 0.2, 0.4, 0.7])
        500x500px 24-bit RGB image

    ::

        sage: external_ray([i/5 for i in range(1,5)])
        500x500px 24-bit RGB image

    WARNING:

    If you are passing in an image, make sure you specify
    which parameters to use when drawing the external ray.
    For example, the following is incorrect::

        sage: M = mandelbrot_plot(x_center=0)  # not tested
        sage: external_ray(5/7, image=M)       # not tested
        500x500px 24-bit RGB image

    To get the correct external ray, we adjust our parameters::

        sage: M = mandelbrot_plot(x_center=0)
        sage: external_ray(5/7, x_center=0, image=M)
        500x500px 24-bit RGB image

    .. TODO::

        The ``copy()`` function for bitmap images needs to be implemented
        in Sage.
    """
    x_0 = kwds.get("x_center", -1)
    y_0 = kwds.get("y_center", 0)
    plot_width = kwds.get("image_width", 4)
    pixel_width = kwds.get("pixel_count", 500)
    depth = kwds.get("D", 25)
    sharpness = kwds.get("S", 10)
    radial_parameter = kwds.get("R", 100)
    precision = kwds.get("prec", 300)
    precision = max(precision, -logb(pixel_width * 0.001, 2).round() + 10)
    ray_color = kwds.get("ray_color", [255] * 3)
    image = kwds.get("image", None)
    if image is None:
        image = mandelbrot_plot(x_center=x_0, **kwds)

    # Make a copy of the bitmap image.
    old_pixel = image.pixels()
    M = Image('RGB', (pixel_width, pixel_width))
    pixel = M.pixels()
    for i in range(pixel_width):
        for j in range(pixel_width):
            pixel[i, j] = old_pixel[i, j]

    # Make sure that theta is a list so loop below works
    if type(theta) != list:
        theta = [theta]

    # Check if theta is in the interval [0,1]
    for angle in theta:
        if angle < 0 or angle > 1:
            raise ValueError("values for theta must be in "
                             "the closed interval [0,1].")

    # Loop through each value for theta in list and plot the external ray.
    for angle in theta:
        E = fast_external_ray(angle, D=depth, S=sharpness, R=radial_parameter,
                              prec=precision, image_width=plot_width,
                              pixel_count=pixel_width)

        # Convert points to pixel coordinates.
        pixel_list = convert_to_pixels(E, x_0, y_0, plot_width, pixel_width)

        # Find the pixels between points in pixel_list.
        extra_points = []
        for i in range(len(pixel_list) - 1):
            if min(pixel_list[i + 1]) >= 0 and max(pixel_list[i + 1]) < pixel_width:
                for j in get_line(pixel_list[i], pixel_list[i + 1]):
                    extra_points.append(j)

        # Add these points to pixel_list to fill in gaps in the ray.
        pixel_list += extra_points

        # Remove duplicates from list.
        pixel_list = list(set(pixel_list))

        # Check if point is in window and if it is, plot it on the image to
        # create an external ray.
        for k in pixel_list:
            if max(k) < pixel_width and min(k) >= 0:
                pixel[int(k[0]), int(k[1])] = tuple(ray_color)
    return M


def julia_plot(c=-1,
               x_center=0.0,
               y_center=0.0,
               image_width=4.0,
               max_iteration=500,
               pixel_count=500,
               base_color='steelblue',
               level_sep=1,
               number_of_colors=50,
               point_color='yellow',
               interact=False,
               mandelbrot=True,
               period=None):
    r"""
    Plots the Julia set of a given complex `c` value. Users can specify whether
    they would like to display the Mandelbrot side by side with the Julia set.

    The Julia set of a given `c` value is the set of complex numbers for which
    the function `Q_c(z)=z^2+c` is bounded under iteration. The Julia set can
    be visualized by plotting each point in the set in the complex plane.
    Julia sets are examples of fractals when plotted in the complex plane.

    ALGORITHM:

    Define the map `Q_c(z) = z^2 + c` for some `c \in \mathbb{C}`. For every
    `p \in \mathbb{C}`, if `|Q_{c}^{k}(p)| > 2` for some `k \geq 0`,
    then `Q_{c}^{n}(p) \to \infty`. Let `N` be the maximum number of iterations.
    Compute the first `N` points on the orbit of `p` under `Q_c`. If for
    any `k < N`, `|Q_{c}^{k}(p)| > 2`, we stop the iteration and assign a color
    to the point `p` based on how quickly `p` escaped to infinity under
    iteration of `Q_c`. If `|Q_{c}^{i}(p)| \leq 2` for all `i \leq N`, we assume
    `p` is in the Julia set and assign the point `p` the color black.

    INPUT:

    - ``c`` -- complex (optional - default: ``-1``), complex point `c` that
      determines the Julia set.

    - ``period`` -- list (optional - default: ``None``), returns the Julia set
      for a random `c` value with the given (formal) cycle structure.

    - ``mandelbrot`` -- boolean (optional - default: ``True``), when set to
      ``True``, an image of the Mandelbrot set is appended to the right of the
      Julia set.

    - ``point_color`` -- RGB color (optional - default: ``'tomato'``),
      color of the point `c` in the Mandelbrot set (any valid input for Color).

    - ``x_center`` -- double (optional - default: ``-1.0``), Real part
      of center point.

    - ``y_center`` -- double (optional - default: ``0.0``), Imaginary part
      of center point.

    - ``image_width`` -- double (optional - default: ``4.0``), width of image
      in the complex plane.

    - ``max_iteration`` -- long (optional - default: ``500``), maximum number
      of iterations the map `Q_c(z)`.

    - ``pixel_count`` -- long (optional - default: ``500``), side length of
      image in number of pixels.

    - ``base_color`` -- hex color (optional - default: ``'steelblue'``), color
      used to determine the coloring of set (any valid input for Color).

    - ``level_sep`` -- long (optional - default: 1), number of iterations
      between each color level.

    - ``number_of_colors`` -- long (optional - default: 30), number of colors
      used to plot image.

    - ``interact`` -- boolean (optional - default: ``False``), controls whether
      plot will have interactive functionality.

    OUTPUT:

    24-bit RGB image of the Julia set in the complex plane.

    EXAMPLES::

        sage: julia_plot()
        1001x500px 24-bit RGB image

    To display only the Julia set, set ``mandelbrot`` to ``False``::

        sage: julia_plot(mandelbrot=False)
        500x500px 24-bit RGB image

    To display an interactive plot of the Julia set in the Notebook,
    set ``interact`` to ``True``::

        sage: julia_plot(interact=True)
        interactive(children=(FloatSlider(value=-1.0, description=u'Real c'...

    To return the Julia set of a random `c` value with (formal) cycle structure
    `(2,3)`, set ``period = [2,3]``::

        sage: julia_plot(period=[2,3])
        1001x500px 24-bit RGB image

    To return all of the Julia sets of `c` values with (formal) cycle structure
    `(2,3)`::

        sage: period = [2,3] # not tested
        ....: R.<c> = QQ[]
        ....: P.<x,y> = ProjectiveSpace(R,1)
        ....: f = DynamicalSystem([x^2+c*y^2, y^2])
        ....: L = f.dynatomic_polynomial(period).subs({x:0,y:1}).roots(ring=CC)
        ....: c_values = [k[0] for k in L]
        ....: for c in c_values:
        ....:     julia_plot(c)
    """
    if period is not None:
        R = PolynomialRing(QQ, 'c')
        c = R.gen()
        x, y = ProjectiveSpace(R, 1, 'x,y').gens()
        f = DynamicalSystem([x**2 + c * y**2, y**2])
        L = f.dynatomic_polynomial(period).subs({x: 0, y: 1}).roots(ring=CC)
        c = L[randint(0, len(L) - 1)][0]

    c = CC(c)
    c_real = c.real()
    c_imag = c.imag()

    base_color = Color(base_color)
    point_color = Color(point_color)

    if interact:
        from ipywidgets.widgets import FloatSlider, IntSlider, ColorPicker, interact
        widgets = dict(
            c_real = FloatSlider(min=-2.0, max=2.0, step=EPS,
                                 value=c_real, description="Real c"),
            c_imag = FloatSlider(min=-2.0, max=2.0, step=EPS,
                                 value=c_imag, description="Imag c"),
            x_center = FloatSlider(min=-1.0, max=1.0, step=EPS,
                                   value=x_center, description="Real center"),
            y_center = FloatSlider(min=-1.0, max=1.0, step=EPS,
                                   value=y_center, description="Imag center"),
            image_width = FloatSlider(min=EPS, max=4.0, step=EPS,
                                      value=image_width, description="Zoom"),
            max_iteration = IntSlider(min=0, max=600,
                                      value=max_iteration, description="Iterations"),
            pixel_count = IntSlider(min=10, max=600,
                                    value=pixel_count, description="Pixels"),
            level_sep = IntSlider(min=1, max=20,
                                  value=level_sep, description="Color sep"),
            color_num = IntSlider(min=1, max=100,
                                  value=number_of_colors, description="# Colors"),
            base_color = ColorPicker(value=base_color.html_color(),
                                     description="Base color"),
        )
        if mandelbrot:
            widgets["point_color"] = ColorPicker(value=point_color.html_color(),
                                                 description="Point color")
            return interact(**widgets).widget(julia_helper)
        else:
            return interact(**widgets).widget(fast_julia_plot)

    if mandelbrot:
        return julia_helper(c_real, c_imag, x_center, y_center,
                            image_width, max_iteration, pixel_count,
                            level_sep,
                            number_of_colors, base_color, point_color)

    else:
        return fast_julia_plot(c_real, c_imag, x_center, y_center,
                               image_width, max_iteration, pixel_count,
                               level_sep,
                               number_of_colors, base_color)

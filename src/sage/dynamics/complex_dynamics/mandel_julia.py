r"""
Mandelbrot and Julia sets

Plots the Mandelbrot and Julia sets for the map `Q_c(z)=z^2+c` in the complex
plane.

The Mandelbrot set is the set of complex numbers `c` for which the function
`Q_c(z)=z^2+c` does not diverge when iterated from `z = 0`. This set of complex
numbers can be visualized by plotting each value for `c` in the complex plane.
The Mandelbrot set is an example of a fractal when plotted in the complex plane.

The Julia set for a given `c` is the set of complex numbers for which the
function `Q_c(z)=z^2+c` is bounded under iteration.

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
                                                                fast_julia_plot,
                                                                julia_helper)
from sagenb.notebook.interact import interact
from sagenb.notebook.interact import slider
from sagenb.notebook.interact import input_box
from sagenb.notebook.interact import color_selector
from sagenb.notebook.interact import checkbox
from sage.plot.colors import Color
from sage.repl.image import Image
from sage.functions.log import function_log as log
from sage.rings.rational_field import QQ
from sage.rings.all import CC
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.categories.homset import End
from sage.misc.prandom import randint

def mandelbrot_plot(**kwds):
    r"""
    Interactive plot of the Mandelbrot set for the map `Q_c(z) = z^2 + c`.

    ALGORITHM:

    Let each pixel in the image be a point `c \in \mathbb{C}` and define the
    map `Q_c(z) = z^2 + c`. If `|Q_{c}^{k}(c)| > 2` for some `k \geq 0`, it
    follows that `Q_{c}^{n}(c) \to \infty`. Let `N` be the maximum number of
    iterations. Compute the first `N` points on the orbit of `0` under `Q_c`.
    If for any `k < N`, `|Q_{c}^{k}(0)| > 2`, we stop the iteration and assign
    a color to the point `c` based on how quickly `0` escaped to infinity under
    iteration of `Q_c`. If `|Q_{c}^{i}(0)| \leq 2` for all `i \leq N`, we assume
    `c` is in the Mandelbrot set and assign the point `c` the color black.

    REFERENCE:

    [Devaney]_

    kwds:

    - ``x_center`` -- double (optional - default: ``-1.0``), Real part of center point.

    - ``y_center`` -- double (optional - default: ``0.0``), Imaginary part of center point.

    - ``image_width`` -- double (optional - default: ``4.0``), width of image in the complex plane.

    - ``max_iteration`` -- long (optional - default: ``500``), maximum number of iterations the map ``Q_c(z)``.

    - ``pixel_count`` -- long (optional - default: ``500``), side length of image in number of pixels.

    - ``base_color`` -- RGB color (optional - default: ``[40, 40, 40]``) color used to determine the coloring of set.

    - ``iteration_level`` -- long (optional - default: 1) number of iterations between each color level

    - ``number_of_colors`` -- long (optional - default: 30) number of colors used to plot image

    - ``interact`` -- boolean (optional - default: ``False``), controls whether plot will have interactive functionality.

    OUTPUT:

    24-bit RGB image of the Mandelbrot set in the complex plane.

    EXAMPLES:

    ::

        sage: mandelbrot_plot() # long time
        500x500px 24-bit RGB image

    ::

        sage: mandelbrot_plot(pixel_count=1000) # long time
        1000x1000px 24-bit RGB image

    ::

        sage: mandelbrot_plot(x_center=-1.11, y_center=0.2283, image_width=1/128, # long time
        ....: max_iteration=2000, number_of_colors=500, base_color=[40, 100, 100])
        500x500px 24-bit RGB image

    To display an interactive plot of the Mandelbrot set in the Notebook, set ``interact`` to ``True``::

        sage: mandelbrot_plot(interact=True)
        <html>...</html>

    ::

        sage: mandelbrot_plot(interact=True, x_center=-0.75, y_center=0.25,
        ....: image_width=1/2, number_of_colors=75)
        <html>...</html>
    """

    x_center = kwds.pop("x_center", -1.0)
    y_center = kwds.pop("y_center", 0.0)
    image_width = kwds.pop("image_width", 4.0)
    max_iteration = kwds.pop("max_iteration", 500)
    pixel_count = kwds.pop("pixel_count", 500)
    base_color = kwds.pop("base_color", [40, 40, 40])
    iteration_level = kwds.pop("iteration_level", 1)
    number_of_colors = kwds.pop("number_of_colors", 30)
    interacts = kwds.pop("interact", False)

    if interacts:
        @interact(layout={'bottom':[['real_center'], ['im_center'], ['width']],
         'top':[['iterations'], ['level_sep'], ['color_num'], ['image_color']]})
        def _(real_center=input_box(x_center, 'Real'),
            im_center=input_box(y_center, 'Imaginary'),
            width=input_box(image_width, 'Width of Image'),
            iterations=input_box(max_iteration, 'Max Number of Iterations'),
            level_sep=input_box(iteration_level, 'Iterations between Colors'),
            color_num=input_box(number_of_colors, 'Number of Colors'),
            image_color=color_selector(default=Color([j/255 for j in base_color]),
             label="Image Color", hide_box=True)):
            return fast_mandelbrot_plot(real_center, im_center, width,
             iterations, pixel_count, level_sep, color_num, image_color).show()

    else:
        return fast_mandelbrot_plot(x_center, y_center, image_width, max_iteration,
         pixel_count, iteration_level, number_of_colors, base_color)

def julia_plot(c=-1, **kwds):
    r"""
    Plots the Julia set of a given complex `c` value. Users can specify whether
    they would like to diplay the Mandelbrot side by side with the Julia set.

    The Julia set of a given `c` value is the set of complex numbers for which
    the function `Q_c(z)=z^2+c` is bounded under iteration. The Julia set can
    be visualized by plotting each point in the set in the complex plane.
    Julia sets are examples of fractals when plotted in the complex plane.

    ALGORITHM:

    Define the map `Q_c(z) = z^2 + c` for some `c \in \mathbb{C}`. For every
    `p \in \mathbb{C}`, if `|Q_{p}^{k}(c)| > 2` for some `k \geq 0`,
    then `Q_{p}^{n}(c) \to \infty`. Let `N` be the maximum number of iterations.
    Compute the first `N` points on the orbit of `p` under `Q_c`. If for
    any `k < N`, `|Q_{p}^{k}(c)| > 2`, we stop the iteration and assign a color
    to the point `c` based on how quickly `c` escaped to infinity under
    iteration of `Q_c`. If `|Q_{p}^{i}(c)| \leq 2` for all `i \leq N`, we assume
    `p` is in the Julia set and assign the point `p` the color black.

    INPUT:

    - ``c`` -- complex (optional - default: ``-1``), complex point `c` that determines the Julia set.

    kwds:

    - ``period`` -- list (optional - default: ``None``), returns the Julia set for a random `c` value with the given cycle structure.

    - ``mandelbrot`` -- boolean (optional - default: ``True``), when set to ``True``, an image of the Mandelbrot set is appended to the right of the Julia set.

    - ``point_color`` -- RGB color (optional - default: ``[255, 0, 0]``), color of the point `c` in the Mandelbrot set.

    - ``x_center`` -- double (optional - default: ``-1.0``), Real part of center point.

    - ``y_center`` -- double (optional - default: ``0.0``), Imaginary part of center point.

    - ``image_width`` -- double (optional - default: ``4.0``), width of image in the complex plane.

    - ``max_iteration`` -- long (optional - default: ``500``), maximum number of iterations the map `Q_c(z)`.

    - ``pixel_count`` -- long (optional - default: ``500``), side length of image in number of pixels.

    - ``base_color`` -- RGB color (optional - default: ``[40, 40, 40]``), color used to determine the coloring of set.

    - ``iteration_level`` -- long (optional - default: 1), number of iterations between each color level.

    - ``number_of_colors`` -- long (optional - default: 30), number of colors used to plot image.

    - ``interact`` -- boolean (optional - default: ``False``), controls whether plot will have interactive functionality.

    OUTPUT:

    24-bit RGB image of the Julia set in the complex plane.

    EXAMPLES::

        sage: julia_plot()
        1001x500px 24-bit RGB image

    To display only the Julia set, set ``mandelbrot`` to ``False``::

        sage: julia_plot(mandelbrot=False)
        500x500px 24-bit RGB image

    To display an interactive plot of the Julia set in the Notebook, set ``interact`` to ``True``::

        sage: julia_plot(interact=True)
        <html>...</html>

    To return the Julia set of a random `c` value with cycle structure `(2,3)`, set ``period = [2,3]``::

        sage: julia_plot(period=[2,3])
        1001x500px 24-bit RGB image

    To return all of the Julia sets of `c` values with cycle structure `(2,3)`::

        sage: period = [2,3] # not tested
        ....: R.<c> = QQ[]
        ....: P.<x,y> = ProjectiveSpace(R,1)
        ....: R = P.coordinate_ring()
        ....: H = End(P)
        ....: f = H([x^2+c*y^2,y^2])
        ....: L = f.dynatomic_polynomial(period).subs({x:0,y:1}).roots(ring=CC)
        ....: c_values = [k[0] for k in L]
        ....: for c in c_values:
        ....:     julia_plot(c)
    """

    x_center = kwds.pop("x_center", 0.0)
    y_center = kwds.pop("y_center", 0.0)
    image_width = kwds.pop("image_width", 4.0)
    max_iteration = kwds.pop("max_iteration", 500)
    pixel_count = kwds.pop("pixel_count", 500)
    base_color = kwds.pop("base_color", [50, 50, 50])
    iteration_level = kwds.pop("iteration_level", 1)
    number_of_colors = kwds.pop("number_of_colors", 50)
    point_color = kwds.pop("point_color", [255, 0, 0])
    interacts = kwds.pop("interact", False)
    mandelbrot = kwds.pop("mandelbrot", True)
    period = kwds.pop("period", None)

    if period is not None:
        R = PolynomialRing(QQ, 'c')
        c = R.gen()
        P = ProjectiveSpace(R, 1, 'x,y')
        x,y = P.gens()
        H = End(P)
        f = H([x**2+c*y**2, y**2])
        L = f.dynatomic_polynomial(period).subs({x:0,y:1}).roots(ring=CC)
        c = L[randint(0,len(L)-1)][0]

    c_real = CC(c).real()
    c_imag = CC(c).imag()

    if interacts:
        @interact(layout={'bottom':[['real_center'], ['im_center'], ['width']],
         'top':[['iterations'], ['level_sep'], ['color_num'], ['mandel'],
         ['cx'], ['cy']], 'right':[['image_color'], ['pt_color']]})
        def _(cx = input_box(c_real, '$Re(c)$'),
            cy = input_box(c_imag, '$Im(c)$'),
            real_center=input_box(x_center, 'Real Center'),
            im_center=input_box(y_center, 'Imaginary Center'),
            width=input_box(image_width, 'Width of Image'),
            iterations=input_box(max_iteration, 'Max Number of Iterations'),
            level_sep=input_box(iteration_level, 'Iterations between Colors'),
            color_num=input_box(number_of_colors, 'Number of Colors'),
            image_color=color_selector(default=Color([j/255 for j in base_color]),
             label="Image Color", hide_box=True),
            pt_color=color_selector(default=Color([j/255 for j in point_color]),
             label="Point Color", hide_box=True),
            mandel=checkbox(mandelbrot, label='Mandelbrot set')):

            if mandel:
                return julia_helper(cx, cy, real_center, im_center,
                 width, iterations, pixel_count, level_sep, color_num,
                 image_color, pt_color).show()

            else:
                return fast_julia_plot(cx, cy, real_center, im_center,
                 width, iterations, pixel_count, level_sep, color_num,
                 image_color).show()

    else:
        if mandelbrot:
            return julia_helper(c_real, c_imag, x_center, y_center,
             image_width, max_iteration, pixel_count, iteration_level,
             number_of_colors, base_color, point_color)

        else:
            return fast_julia_plot(c_real, c_imag, x_center, y_center,
             image_width, max_iteration, pixel_count, iteration_level,
             number_of_colors, base_color)

r"""
Mandelbrot and Julia sets

Plots the Mandelbrot set for the map `Q_c(z)=z^2+c` in the complex plane.

The Mandelbrot set is the set of complex numbers `c` for which the function
`Q_c(z)=z^2+c` does not diverge when iterated from `z = 0`. This set of complex
numbers can be visualized by plotting each value for `c` in the complex plane.
The Mandelbrot set is an example of a fractal when plotted in the complex plane.

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
from sage.dynamics.complex_dynamics.mandel_julia_helper import \
 fast_mandelbrot_plot, fast_external_ray, convert_to_pixels, get_line
from sagenb.notebook.interact import interact, slider, input_box, color_selector
from sage.plot.colors import Color
from sage.repl.image import Image
from sage.functions.log import function_log as log


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

    [Dev2005]_

    kwds:

    - ``x_center`` -- double (optional - default: ``-1.0``), Real part of center point.

    - ``y_center`` -- double (optional - default: ``0.0``), Imaginary part of center point.

    - ``image_width`` -- double (optional - default: ``4.0``), width of image in the complex plane.

    - ``max_iteration`` -- long (optional - default: ``500``), maximum number of iterations the map ``Q_c(z)``.

    - ``pixel_count`` -- long (optional - default: ``500``), side length of image in number of pixels.

    - ``base_color`` -- RGB color (optional - default: ``[40, 40, 40]``) color used to determine the coloring of set.

    - ``iteration_level`` -- long (optional - default: 1) number of iterations between each color level.

    - ``number_of_colors`` -- long (optional - default: 30) number of colors used to plot image.

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
            return fast_mandelbrot_plot(real_center, im_center, width, iterations,
             pixel_count, level_sep, color_num, image_color).show()

    else:
        return fast_mandelbrot_plot(x_center, y_center, image_width, max_iteration,
         pixel_count, iteration_level, number_of_colors, base_color)

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

    - ``image`` -- 24-bit RGB image (optional - default: None) user specified image of Mandelbrot set.

    - ``D`` -- long (optional - default: ``25``) depth of the approximation. As ``D`` increases, the external ray gets closer to the boundary of the Mandelbrot set. If the ray doesn't reach the boundary of the Mandelbrot set, increase ``D``.

    - ``S`` -- long (optional - default: ``10``) sharpness of the approximation. Adjusts the number of points used to approximate the external ray (number of points is equal to ``S*D``). If ray looks jagged, increase ``S``.

    - ``R`` -- long (optional - default: ``100``) radial parameter. If ``R`` is large, the external ray reaches sufficiently close to infinity. If ``R`` is too small, Newton's method may not converge to the correct ray.

    - ``prec`` -- long (optional - default: ``300``) specifies the bits of precision used by the Complex Field when using Newton's method to compute points on the external ray.

    - ``ray_color`` -- RGB color (optional - default: ``[255, 255, 255]``) color of the external ray(s).

    OUTPUT:

    24-bit RGB image of external ray(s) on the Mandelbrot set.

    EXAMPLES::

        sage: external_ray(1/3)
        500x500px 24-bit RGB image

    ::

        sage: external_ray(0.6, ray_color=[255, 0, 0])
        500x500px 24-bit RGB image

    ::

        sage: external_ray([0, 0.2, 0.4, 0.7]) # long time
        500x500px 24-bit RGB image

    ::

        sage: external_ray([i/5 for i in range(1,5)]) # long time
        500x500px 24-bit RGB image

    WARNING:

    If you are passing in an image, make sure you specify
    which parameters to use when drawing the external ray.
    For example, the following is incorrect::

        sage: M = mandelbrot_plot(x_center=0) # not tested
        sage: external_ray(5/7, image=M) # not tested
        500x500px 24-bit RGB image

    To get the correct external ray, we adjust our parameters::

        sage: M = mandelbrot_plot(x_center=0) # not tested
        sage: external_ray(5/7, x_center=0, image=M) # not tested
        500x500px 24-bit RGB image

    TODO:

    The ``copy()`` function for bitmap images needs to be implemented in Sage.
    """

    x_0 = kwds.get("x_center", -1)
    y_0 = kwds.get("y_center", 0)
    plot_width = kwds.get("image_width", 4)
    pixel_width = kwds.get("pixel_count", 500)
    depth = kwds.get("D", 25)
    sharpness = kwds.get("S", 10)
    radial_parameter = kwds.get("R", 100)
    precision = kwds.get("prec", 300)
    precision = max(precision, -log(pixel_width * 0.001, 2).round() + 10)
    ray_color = kwds.get("ray_color", [255]*3)
    image = kwds.get("image", None)
    if image is None:
        image = mandelbrot_plot(**kwds)

    # Make a copy of the bitmap image.
    # M = copy(image)
    old_pixel = image.pixels()
    M = Image('RGB', (pixel_width, pixel_width))
    pixel = M.pixels()
    for i in range(pixel_width):
        for j in range(pixel_width):
            pixel[i,j] = old_pixel[i,j]

    # Make sure that theta is a list so loop below works
    if type(theta) != list:
        theta = [theta]

    # Check if theta is in the invterval [0,1]
    for angle in theta:
        if angle < 0 or angle > 1:
            raise \
            ValueError("values for theta must be in the closed interval [0,1].")

    # Loop through each value for theta in list and plot the external ray.
    for angle in theta:
        E = fast_external_ray(angle, D=depth, S=sharpness, R=radial_parameter,
         prec=precision, image_width=plot_width, pixel_count=pixel_width)

        # Convert points to pixel coordinates.
        pixel_list = convert_to_pixels(E, x_0, y_0, plot_width, pixel_width)

        # Find the pixels between points in pixel_list.
        extra_points = []
        for i in range(len(pixel_list) - 1):
            if min(pixel_list[i+1]) >= 0 and max(pixel_list[i+1]) < pixel_width:
                for j in get_line(pixel_list[i], pixel_list[i+1]):
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

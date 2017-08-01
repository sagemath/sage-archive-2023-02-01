r"""
Mandelbrot and Julia sets (Cython helper)

This is the helper file providing functionality for mandel_julia.py.

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
from sage.plot.colors import Color
from sage.repl.image import Image
from copy import copy
from cysignals.signals cimport sig_check

def fast_mandelbrot_plot(double x_center, double y_center, double image_width,
 long max_iteration, long pixel_count, long level_sep, long color_num, base_color):

    r"""
    Plots the Mandelbrot set in the complex plane for the map `Q_c(z) = z^2 + c`.

    INPUT:

    - ``x_center`` -- double, real part of the center point in the complex plane.

    - ``y_center`` -- double, imaginary part of the center point in the complex plane.

    - ``image_width`` -- double, width of the image in the complex plane.

    - ``max_iteration`` -- long, maximum number of iterations the map `Q_c(z)` considered.

    - ``pixel_count`` -- long, side length of image in number of pixels.

    - ``level_sep`` -- long, number of iterations between each color level

    - ``color_num`` -- long, number of colors used to plot image

    - ``base_color`` -- list, RGB color used to determine the coloring of set.

    OUTPUT:

    - A 24-bit RGB image of the Mandelbrot set in the complex plane

    EXAMPLES:

    Plot the Mandelbrot set with the center point `-1 + 0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandelbrot_plot
        sage: fast_mandelbrot_plot(-1, 0, 4, 500, 600, 1, 20, [40, 40, 40]) # long time
        600x600px 24-bit RGB image

    We can focus on smaller parts of the set by adjusting image_width::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandelbrot_plot
        sage: fast_mandelbrot_plot(-1.11, 0.2283, 1/128, 2000, 500, 1, 500, [40, 100, 100]) # long time
        500x500px 24-bit RGB image
    """

    cdef long i, j, col, row, level, color_value, iteration
    cdef double k, x_corner, y_corner, step_size, x_coor, y_coor, new_x, new_y
    cdef M, pixel, color_list

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Initialize an image to the color black and access the pixels
    M = Image("RGB", (pixel_count,pixel_count), 'black')
    pixel = M.pixels()

    # Take the given base color and create a list of evenly spaced
    # colors between the given base color and white. The number of
    # colors in the list depends on the variable color_num.
    if type(base_color) == Color:
        # Convert Color to RGB list
        base_color = [int(k*255) for k in base_color]
    color_list = []
    for i in range(color_num):
        sig_check()
        color_list.append(copy(base_color))
        for j in range(3):
            color_list[i][j] += i * (255 - color_list[i][j]) // color_num
        color_list[i] = tuple(color_list[i])

    # First, we determine the complex coordinates of the point in the top left
    # corner of the image. Then, we loop through each pixel in the image and
    # assign it complex coordinates relative to the image's top left corner.
    x_corner = x_center - image_width/2
    y_corner = y_center + image_width/2
    step_size = image_width / pixel_count
    for col in range(pixel_count):
        x_coor = x_corner + col*step_size
        for row in range(pixel_count):
            sig_check()
            y_coor = y_corner - row*step_size

            # We compute the orbit of 0 under the map Q(z) = z^2 + c
            # until we either reach the maximum number of iterations
            # or find a point in the orbit with modulus greater than 2
            new_x, new_y = 0.0, 0.0
            iteration = 0
            while (new_x**2 + new_y**2 <= 4.0 and iteration < max_iteration):
                sig_check()
                new_x, new_y = new_x**2 - new_y**2 + x_coor, \
                 2*new_x*new_y + y_coor
                iteration += 1

            # If the point escapes to infinity, assign the point a color
            # based on how fast it escapes. The more iterations it takes for
            # a point to escape to infinity, the lighter its color will be.
            # Otherwise, assume the point is in the Mandelbrot set and leave
            # it black.
            if iteration != max_iteration:
                # Assign each point a level based on its number of iterations.
                level = iteration // level_sep
                # Assign the pixel a color based on it's level. If we run out
                # of colors, assign it the last color in the list.
                if level < color_num:
                    pixel[col,row] = color_list[level]
                else:
                    pixel[col,row] = color_list[-1]
    return M

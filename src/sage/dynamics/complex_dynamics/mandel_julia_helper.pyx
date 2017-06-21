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

def fast_mandel_plot(float x_center, float y_center, float image_width,
 int max_iteration, int pixel_count, int level_sep, int color_num, base_color):

    r"""
    Plots the Mandelbrot set in the complex plane for the map `Q_c(z) = z^2 + c`.

    INPUT:

    - ``x_center`` -- float, real part of the center point in the complex plane.

    - ``y_center`` -- float, imaginary part of the center point in the complex plane.

    - ``image_width`` -- float, width of the image in the complex plane.

    - ``max_iteration`` -- int, maximum number of iterations the map `f(z)` considered.

    - ``pixel_count`` -- int, side length of image in number of pixels.

    - ``level_sep`` -- int, number of iterations between each color level

    - ``color_num`` -- int, number of colors used to plot image

    - ``base_color`` -- list, RGB color used to determine the coloring of set.

    OUTPUT:

    - A 24-bit RGB image of the Mandelbrot set in the complex plane

    EXAMPLES:

    Plot the Mandelbrot set with the center point `-1 + 0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandel_plot
        sage: fast_mandel_plot(-1, 0, 4, 500, 600, 1, 20, [40, 40, 40])
        Launched png viewer for 500x500px 24-bit RGB image

    We can focus on smaller parts of the set by adjusting image_width::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandel_plot
        sage: fast_mandel_plot(-0.75, 0.10, 1/4, 500, 600, 10, 25, [40, 40, 40])
        Launched png viewer for 500x500px 24-bit RGB image
    """

    from sage.plot.colors import Color
    from sage.repl.image import Image
    from copy import copy

    cdef int color_value, row, col, iteration, i, j
    cdef float new_x, new_y, x_coor, y_coor
    cdef color_list, M, pixel

    # reflect image about x-axis
    y_center *= -1

    M = Image("RGB", (pixel_count,pixel_count), 'black') # create image
    pixel = M.pixels() # get pixels

    image_width = abs(image_width)

    color_list = []
    for i in range(color_num):
        color_list.append(copy(base_color))
        for j in range(3):
            color_list[i][j] += i*(255-color_list[i][j])/color_num
        color_list[i] = tuple(color_list[i])

    for row in range(pixel_count):
        x_coor = x_center + image_width*(row-pixel_count/2)/pixel_count # width of image in cartesian coordinates
        for col in range(pixel_count): # loop through pixels
            y_coor = y_center + image_width*(col-pixel_count/2)/pixel_count

            # compute the orbit of 0 under the map f(z) = z^2 + c
            new_x,new_y = (0.0, 0.0)
            iteration = 0

            while (new_x**2 + new_y**2 <= 4.0 and iteration < max_iteration): # escape condition
                new_x,new_y = new_x**2 - new_y**2 + x_coor, 2*new_x*new_y + y_coor
                iteration += 1

            if iteration != max_iteration:
                level = iteration/level_sep

            if level < color_num:
                pixel[row,col] = color_list[level]
            else:
                pixel[row,col] = color_list[-1]

            if iteration == max_iteration:
                pixel[row,col] = (0,0,0)
    return M

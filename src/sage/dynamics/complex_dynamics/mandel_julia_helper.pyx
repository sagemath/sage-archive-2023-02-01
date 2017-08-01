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

cpdef convert_to_pixels(point_list, double x_0, double y_0, double width,
 long number_of_pixels):
    r"""
    Converts cartesian coordinates to pixels within a specified window.

    INPUT:

    - ``point_list`` -- list of tuples, points in cartesian coordinates.

    - ``x_0`` -- double, x-coordinate of the center of the image.

    - ``y_0`` -- double, y-coordinate of the center of the image.

    - ``width`` -- double, width of visible window in caresian coordinates.

    - ``number_of_pixels`` -- long, width of image in pixels.

    OUTPUT:

    List of tuples of integers representing the pixels.

    EXAMPLES::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import convert_to_pixels
        sage: convert_to_pixels([(-1,3),(0,-4),(5,0)], 0, 0, 12, 100)
        [(42, 25), (50, 83), (92, 50)]
    """
    cdef:
        k, pixel_list, x_corner, y_corner, step_size
        long x_pixel, y_pixel
    pixel_list = []

    # Compute top left corner of window and step size
    x_corner = x_0 - width/2
    y_corner = y_0 + width/2
    step_size = number_of_pixels / width

    # Convert each point in list to pixel coordinates
    for k in point_list:
        sig_check()
        x_pixel = round((k[0] - x_corner) * step_size)
        y_pixel = round((y_corner - k[1]) * step_size)
        pixel_list.append((x_pixel, y_pixel))
    return pixel_list

cpdef fast_julia_plot(double c_real, double c_imag,
  double x_center=0, double y_center=0, double image_width=4,
  long max_iteration=500, long pixel_count=500, long level_sep=2,
  long color_num=40, base_color=[50, 50, 50]):
    r"""
    Plots the Julia set for a given `c` value in the complex plane for the map `Q_c(z) = z^2 + c`.

    INPUT:

    - ``c_real`` -- double, Real part of `c` value that determines Julia set.

    - ``c_imag`` -- double, Imaginary part of `c` value that determines Julia set.

    - ``x_center`` -- double (optional - default: ``0.0``), Real part of center point.

    - ``y_center`` -- double (optional - default: ``0.0``), Imaginary part of center point.

    - ``image_width`` -- double (optional - default: ``4.0``), width of image in the complex plane.

    - ``max_iteration`` -- long (optional - default: ``500``), maximum number of iterations the map ``Q_c(z)``.

    - ``pixel_count`` -- long (optional - default: ``500``), side length of image in number of pixels.

    - ``level_sep`` -- long (optional - default: ``2``), number of iterations between each color level.

    - ``color_num`` -- long (optional - default: ``40``), number of colors used to plot image.

    - ``base_color`` -- RGB color (optional - default: ``[50, 50, 50]``), color used to determine the coloring of set.

    OUTPUT:

    24-bit RGB image of the Julia set in the complex plane.

    EXAMPLES:

    Plot the Julia set for `c=-1+0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_julia_plot
        sage: fast_julia_plot(-1, 0)
        500x500px 24-bit RGB image

    We can focus on smaller parts of the set by adjusting ``image_width``::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_julia_plot
        sage: fast_julia_plot(-0.7, 0.3, x_center=.75, image_width=0.01, color_num=100)
        500x500px 24-bit RGB image
    """

    cdef long i, j, col, row, level, color_value, iteration
    cdef double k, x_corner, y_corner, step_size, x_coor, y_coor, new_x, new_y
    cdef M, pixel, color_list

    # Make sure image_width is positive
    image_width = abs(image_width)

    # Initialize an image to the color black and access the pixels
    J = Image("RGB", (pixel_count,pixel_count), 'black')
    Jp = J.pixels()

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

            # We compute the orbit of c under the map Q(z) = z^2 + c
            # until we either reach the maximum number of iterations
            # or find a point in the orbit with modulus greater than 2
            new_x, new_y = x_coor, y_coor
            iteration = 0
            while (new_x**2 + new_y**2 <= 4.0 and iteration < max_iteration):
                sig_check()
                new_x, new_y = new_x**2 - new_y**2 + c_real, \
                 2*new_x*new_y + c_imag
                iteration += 1

            # If the point escapes to infinity, assign the point a color
            # based on how fast it escapes. The more iterations it takes for
            # a point to escape to infinity, the lighter its color will be.
            # Otherwise, assume the point is in the Julia set and leave
            # it black.
            if iteration != max_iteration:
                # Assign each point a level based on its number of iterations.
                level = iteration // level_sep
                # Assign the pixel a color based on it's level. If we run out
                # of colors, assign it the last color in the list.
                if level < color_num:
                    Jp[col,row] = color_list[level]
                else:
                    Jp[col,row] = color_list[-1]

    return J

cpdef julia_helper(double c_real, double c_imag, double x_center=0,
 double y_center=0, double image_width=4, long max_iteration=500,
 long pixel_count=500, long level_sep=2, long color_num=40,
 base_color=[50, 50, 50], point_color=[255, 0, 0]):

    r"""
    Helper function that returns the image of a Julia set for a given
    `c` value side by side with the Mandelbrot set with a point denoting
    the `c` value.

    INPUT:

    - ``c_real`` -- double, Real part of `c` value that determines Julia set.

    - ``c_imag`` -- double, Imaginary part of `c` value that determines Julia set.

    - ``x_center`` -- double (optional - default: ``0.0``), Real part of center point.

    - ``y_center`` -- double (optional - default: ``0.0``), Imaginary part of center point.

    - ``image_width`` -- double (optional - default: ``4.0``), width of image in the complex plane.

    - ``max_iteration`` -- long (optional - default: ``500``), maximum number of iterations the map ``Q_c(z)``.

    - ``pixel_count`` -- long (optional - default: ``500``), side length of image in number of pixels.

    - ``level_sep`` -- long (optional - default: ``2``), number of iterations between each color level.

    - ``color_num`` -- long (optional - default: ``40``), number of colors used to plot image.

    - ``base_color`` -- RGB color (optional - default: ``[50, 50, 50]``), color used to determine the coloring of set.

    - ``point_color`` -- RGB color (optional - default: ``[255, 0, 0]``), color of the point `c` in the Mandelbrot set.

    OUTPUT:

    24-bit RGB image of the Julia and Mandelbrot sets in the complex plane.

    EXAMPLES:

    Plot the Julia set for `c=-1+0i`::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import julia_helper
        sage: julia_helper(-1,0)
        1001x500px 24-bit RGB image

    We can focus on smaller parts of the set by adjusting ``image_width``::

        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import julia_helper
        sage: julia_helper(-.5, .6, y_center=0.178, image_width=0.01) # long time
        1001x500px 24-bit RGB image
    """

    cdef:
        int i, j
        M, Mp, G, Gp, J, Jp, CP

    # Initialize the Julia set
    J = fast_julia_plot(c_real, c_imag, x_center, y_center, image_width,
     max_iteration, pixel_count, level_sep, color_num, base_color)
    Jp = J.pixels()

    # Initialize the image with Julia set on left side
    # Add white border between images
    G = Image("RGB", (2*pixel_count+1,pixel_count), 'white')
    Gp = G.pixels()
    for i in range(pixel_count):
        for j in range(pixel_count):
            Gp[i,j] = Jp[i,j]

    # Plot the Mandelbrot set on the right side
    M = fast_mandelbrot_plot(-1, 0, 4, 500, pixel_count, 1, 30, base_color)
    Mp = M.pixels()
    for i in range(pixel_count+1,2*pixel_count):
        for j in range(pixel_count):
            Gp[i,j] = Mp[int(i-pixel_count),j]

    # Convert Color to RGB list if necessary
    if type(point_color) == Color:
        point_color = [int(k*255) for k in point_color]

    # Add a cross representing c-value to the Mandelbrot set.
    CP = convert_to_pixels([(c_real, c_imag)], -1, 0, 4, pixel_count)
    for i in range(-3,4):
        # Loop through x and y coordinates and check if they are in image
        if min(CP[0][0]+i, CP[0][1]) >= 0 and \
         max(CP[0][0]+i, CP[0][1]) < pixel_count:
            Gp[CP[0][0]+i+pixel_count+1, CP[0][1]] = tuple(point_color)
        if min(CP[0][0], CP[0][1]+i) >= 0 and \
         max(CP[0][0], CP[0][1]+i) < pixel_count:
            Gp[CP[0][0]+pixel_count+1, CP[0][1]+i] = tuple(point_color)

    return G

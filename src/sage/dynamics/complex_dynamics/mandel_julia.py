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

from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandel_plot
from sagenb.notebook.interact import interact
from sagenb.notebook.interact import slider
from sagenb.notebook.interact import input_box

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

    - ``x_center`` -- float (optional - default: ``-1.0``), Real part of center point.

    - ``y_center`` -- float (optional - default: ``0.0``), Imaginary part of center point.

    - ``image_width`` -- float (optional - default: ``4.0``), width of image in the complex plane.

    - ``max_iteration`` -- int (optional - default: ``500``), maximum number of iterations the map ``f(z)``.

    - ``pixel_count`` -- int (optional - default: ``500``), side length of image in number of pixels.

    - ``base_color`` -- RGB color (optional - default: ``[40, 40, 40]``) color used to determine the coloring of set.

    - ``iteration_level`` -- int (optional - default: 1) number of iterations between each color level

    - ``number_of_colors`` -- int (optional - default: 30) number of colors used to plot image

    - ``interact`` -- boolean (optional - default: ``True``), controls whether plot will have interactive functionality.

    OUTPUT:

    Interactive 24-bit RGB image of the Mandelbrot set in the complex plane.

    EXAMPLES:

    ::

        sage: mandelbrot_plot() # not tested

    ::

        sage: mandelbrot_plot(pixel_count=1000) # not tested

    ::

        sage: mandelbrot_plot(base_color=[70, 40, 240]) # not tested

    ::

        sage: mandelbrot_plot(x_center=-0.75, y_center=0.25, image_width=1/2, number_of_colors=75) # not tested

    To use the function outside of the notebook, you must set ``interact`` to False::

        sage: mandelbrot_plot(interact=False) # not tested
        Launched png viewer for 500x500px 24-bit RGB image

    ::

        sage: mandelbrot_plot(interact=False, x_center=-1.11, y_center=0.2283, image_width=1/128, # not tested
        ....: max_iteration=2000, number_of_colors=500, base_color=[40, 100, 100])
        Launched png viewer for 500x500px 24-bit RGB image

    TESTS:

    sage: mandelbrot_plot()
    <html><!--notruncate-->
            <div padding=6 id="div-interact-0">
              <table width=800px height=20px bgcolor="#c5c5c5" cellpadding=15>
                <tr>
                  <td bgcolor="#f9f9f9" valign=top align=left>
                <table>
                  <tr><td colspan=3><table><tr><td align=right><font color="black">Max Number of Iterations&nbsp;</font></td><td><input type="text" value="500" size=80 onchange="interact(0, {variable: 'iterations', adapt_number: 4, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Iterations between Colors&nbsp;</font></td><td><input type="text" value="1" size=80 onchange="interact(0, {variable: 'level_sep', adapt_number: 5, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Number of Colors&nbsp;</font></td><td><input type="text" value="30" size=80 onchange="interact(0, {variable: 'color_num', adapt_number: 6, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">RGB Color&nbsp;</font></td><td><input type="text" value="[40, 40, 40]" size=80 onchange="interact(0, {variable: 'image_color', adapt_number: 7, value: encode64(this.value)}, 1)"></input></td>
    </tr></table></td></tr>
                  <tr><td></td><td style='width: 100%;'>
            <div id="cell-interact-0"><?__SAGE__START>
              <table border=0 bgcolor="white" width=100%>
                <tr>
                  <td bgcolor="white" align=left valign=top>
                    <pre><?__SAGE__TEXT></pre>
                  </td>
                </tr>
                <tr>
                  <td align=left valign=top><?__SAGE__HTML></td>
                </tr>
              </table><?__SAGE__END>
            </div></td><td></td></tr>
                  <tr><td colspan=3><table><tr><td align=right><font color="black">Real&nbsp;</font></td><td><input type="text" value="-1.0" size=80 onchange="interact(0, {variable: 'real_center', adapt_number: 1, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Imaginary&nbsp;</font></td><td><input type="text" value="0.0" size=80 onchange="interact(0, {variable: 'im_center', adapt_number: 2, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Width of Image&nbsp;</font></td><td><input type="text" value="4.0" size=80 onchange="interact(0, {variable: 'width', adapt_number: 3, value: encode64(this.value)}, 1)"></input></td>
    </tr></table></td></tr>
                </table></td>
                </tr>
              </table>
            </div></html>

    sage: mandelbrot_plot(x_center=-0.75, y_center=0.25, image_width=1/2, number_of_colors=75, base_color=[70, 40, 240], pixel_count=1000)
    <html><!--notruncate-->
            <div padding=6 id="div-interact-0">
              <table width=800px height=20px bgcolor="#c5c5c5" cellpadding=15>
                <tr>
                  <td bgcolor="#f9f9f9" valign=top align=left>
                <table>
                  <tr><td colspan=3><table><tr><td align=right><font color="black">Max Number of Iterations&nbsp;</font></td><td><input type="text" value="500" size=80 onchange="interact(0, {variable: 'iterations', adapt_number: 11, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Iterations between Colors&nbsp;</font></td><td><input type="text" value="1" size=80 onchange="interact(0, {variable: 'level_sep', adapt_number: 12, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Number of Colors&nbsp;</font></td><td><input type="text" value="75" size=80 onchange="interact(0, {variable: 'color_num', adapt_number: 13, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">RGB Color&nbsp;</font></td><td><input type="text" value="[70, 40, 240]" size=80 onchange="interact(0, {variable: 'image_color', adapt_number: 14, value: encode64(this.value)}, 1)"></input></td>
    </tr></table></td></tr>
                  <tr><td></td><td style='width: 100%;'>
            <div id="cell-interact-0"><?__SAGE__START>
              <table border=0 bgcolor="white" width=100%>
                <tr>
                  <td bgcolor="white" align=left valign=top>
                    <pre><?__SAGE__TEXT></pre>
                  </td>
                </tr>
                <tr>
                  <td align=left valign=top><?__SAGE__HTML></td>
                </tr>
              </table><?__SAGE__END>
            </div></td><td></td></tr>
                  <tr><td colspan=3><table><tr><td align=right><font color="black">Real&nbsp;</font></td><td><input type="text" value="-0.750000000000000" size=80 onchange="interact(0, {variable: 'real_center', adapt_number: 8, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Imaginary&nbsp;</font></td><td><input type="text" value="0.250000000000000" size=80 onchange="interact(0, {variable: 'im_center', adapt_number: 9, value: encode64(this.value)}, 1)"></input></td>
    </tr><tr><td align=right><font color="black">Width of Image&nbsp;</font></td><td><input type="text" value="1/2" size=80 onchange="interact(0, {variable: 'width', adapt_number: 10, value: encode64(this.value)}, 1)"></input></td>
    </tr></table></td></tr>
                </table></td>
                </tr>
              </table>
            </div></html>

    sage: mandelbrot_plot(interact=False)

    sage: mandelbrot_plot(interact=False, x_center=-1.11, y_center=0.2283, image_width=1/128, # long time
    ....: max_iteration=2000, number_of_colors=500, base_color=[40, 100, 100])
    """

    x_center = kwds.pop("x_center", -1.0)
    y_center = kwds.pop("y_center", 0.0)
    image_width = kwds.pop("image_width", 4.0)
    max_iteration = kwds.pop("max_iteration", 500)
    pixel_count = kwds.pop("pixel_count", 500)
    base_color = kwds.pop("base_color", [40, 40, 40])
    iteration_level = kwds.pop("iteration_level", 1)
    number_of_colors = kwds.pop("number_of_colors", 30)
    interacts = kwds.pop("interact", True)

    if interacts:
        @interact(layout={'bottom':[['real_center'], ['im_center'], ['width']],
         'top':[['iterations'], ['level_sep'], ['color_num'], ['image_color']]})
        def _(real_center=input_box(x_center, 'Real'),
            im_center=input_box(y_center, 'Imaginary'),
            width=input_box(image_width, 'Width of Image'),
            iterations=input_box(max_iteration, 'Max Number of Iterations'),
            level_sep=input_box(iteration_level, 'Iterations between Colors'),
            color_num=input_box(number_of_colors, 'Number of Colors'),
            image_color=input_box(base_color, 'RGB Color')):
            fast_mandel_plot(real_center, im_center, width, iterations,
             pixel_count, level_sep, color_num, image_color).show()

    else:
        fast_mandel_plot(x_center, y_center, image_width, max_iteration,
         pixel_count, iteration_level, number_of_colors, base_color).show()

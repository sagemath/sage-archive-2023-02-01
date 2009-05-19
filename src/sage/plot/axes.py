r"""
2D Axes

SAGE provides several different axes for its 2D plotting functionality.
The 'look' of the axes are similar to what Mathematica provides for its
2D plotting.

To create the axes we uses matplotlib lines to draw all axes and ticks.
For the axes labels we use matplotlib text.  It should be noted that
we also explicitly turn off the matplotlib axes.

The following axes types are supported:
\begin{itemize}
    \item add_xy__axes -- plain axes in the center of a graphic
    \item add_xy_frame_axes -- axes draw around the perimeter of a graphic
    \item add_xy_matrix_frame_axes -- axes around a \code{matrix_plot}
\end{itemize}

"""
#*****************************************************************************
#  Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com> and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from math import floor, log
from sage.structure.sage_object import SageObject
import sage.misc.misc
from copy import copy

class Axes(SageObject):
    """
    Axes for SAGE 2D Graphics.

    Set all axis properties and then add one of
    the following axes to the current (matplotlib) subplot:
    add_xy__axes
    add_xy_frame_axes
    add_xy_matrix_frame_axes

    """
    def __init__(self, color=(0,0,0), fontsize=8, linewidth=0.6,axes_labels=None,
                 axes_label_color=(0,0,0), tick_label_color=(0,0,0)):
        self.__color = color
        self.__tick_label_color = tick_label_color
        self.__axes_labels = axes_labels
        self.__axes_label_color = axes_label_color
        self.__fontsize = fontsize
        self.__linewidth = linewidth
        self.__draw_x_axis = True
        self.__draw_y_axis = True

    def _tasteful_ticks(self, minval, maxval):
        """
        This function finds spacing for axes tick marks that are well spaced.
        Were 'well spaced' means for any given min and max values
        the tick spacing should look even and visually nice (tasteful).

        """
        minval, maxval = float(minval), float(maxval)
        absmin, absmax = abs(minval), abs(maxval)
        # Initialize the domain flags:
        onlyneg, domneg, onlypos, dompos = False, False, False, False

        # Is the range: *only or dominantly* and  *negative or positive*?
        if absmin > absmax:
            if maxval < 0:
                onlyneg = True
            else:
                domneg  = True

            # Is the stepsize going to be < 1?
            if absmin < 1:
                n = 0
                s = str(absmin).split('.')[1]
                for c in s:
                    n+=1
                    if c != '0':
                        break
                p = -(n-1)
                d0 = eval(s[n-1])
                #string may be of length 1
                try:
                    d1 = eval(s[n])
                except IndexError:
                    d1 = 0

            #the stepsize will be 1 or greater:
            else:
                if absmin >= 10:
                    sl = [s for s in str(int(absmin))]
                    d0 = eval(sl[0])
                    d1 = eval(sl[1])
                    p = len(sl)
                else:
                    sl = str(absmin).split('.')
                    d0 = eval(sl[0])
                    d1 = eval(sl[1][0])
                    p = 1

        else: #this means: abs(minval) < abs(maxval)
            if minval > 0:
                 onlypos = True
            else:
                 dompos = True
            #is the stepsize going to be < 1?
            if absmax < 1:
                n = 0
                s = str(absmax).split('.')[1]
                for c in s:
                    n+=1
                    if c != '0':
                        break
                p = -(n-1)
                d0 = eval(s[n-1])
                try:
                    sn = s[n]
                except IndexError:
                    sn = "0"
                d1 = eval(sn)
            #the stepsize will be 1 or greater:
            else:
                if maxval >= 10:
                    sl = [s for s in str(int(absmax))]
                    d0 = eval(sl[0])
                    d1 = eval(sl[1])
                    p = len(sl)
                else:
                    sl = str(absmax).split('.')
                    d0 = eval(sl[0])
                    d1 = eval(sl[1][0])
                    p = 1

        #choose a step size depending on either
        #1st or 2nd digits, d0 and d1, of given maxval or minval
        o0 = 10**(p-1)
        o1 = 10**(p-2)
        #fundamental step sizes: [1,2,2.5,5,10]
        # 'fundamental' means that we split the x and y
        # ranges into widths from the above step sizes
        if d0 == 1:
            if d1 > 5 and p!=1:
                funda = 2.5
                step = funda*o1
            elif d1 < 5 and p==1:
                funda = 2.5
                step = funda*o1
            elif d1 < 5 and p>1:
                funda = 2.5
                step = funda*o1
            else:
                funda = 5
                step = funda*o1
        elif d0 == 2 or (d0 == 3 and p > 1):
            funda = 5
            step = funda*o1
        elif d0 in [3, 4, 5, 6]:
            funda = 1
            step = funda*o0
        else:
            funda = 2
            step = funda*o0

        #the 'fundamental' range
        fundrange = sage.misc.misc.srange(0, d0*o0 + d1*o1 + step, step)

        #Now find the tick step list for major ticks (tslmajor)
        #'oppaxis' is the positioning value of the other axis.
        if onlyneg:
            tslmajor = self._in_range([-x for x in fundrange], minval, maxval)
            tslmajor.sort()
            oppaxis = tslmajor[-1]

        elif domneg or dompos:
            tslmajor = self._in_range([-x for x in fundrange] + fundrange, minval, maxval)
            tslmajor.sort()
            oppaxis = 0

        else: #onlypos
            tslmajor = self._in_range(fundrange, minval, maxval)
            tslmajor.sort()
            oppaxis = tslmajor[0]

        return tslmajor, oppaxis, step

    def _in_range(self, v, minval, maxval):
        "find axis values set in a given range"
        return list(set([x for x in v if x >= minval and x <= maxval]))

    def _trunc(self, x, digits_before_the_decimal):
        s = '%f'%float(x)
        i = s.find('.')
        t = s[:i - digits_before_the_decimal]
        if digits_before_the_decimal > 0:
            t += '0'* digits_before_the_decimal
        return float(eval(t))

    def _format_tick_string(self, s):
        s = str(s)
        if s[-2:] == '.0':
            return s[:-2]
        return s

    def _tasteless_ticks(self, minval, maxval, num_pieces):
        minval0 = minval
        maxval0 = maxval
        rnd   = int(floor(log(maxval - minval)/log(10)))
        if rnd < 0:
            rnd -= 1

        step  = (maxval - minval)/float(num_pieces)
        minval = self._trunc(minval, rnd)
        maxval = self._trunc(maxval + step, rnd)

        step  = (maxval - minval)/float(num_pieces)
        tslmajor = sage.misc.misc.srange(minval, minval+(num_pieces+1)*step, step)
        tslmajor = self._in_range(tslmajor, minval0, maxval0)

        oppaxis = 0
        if maxval <= 0:  # only negative
            oppaxis = tslmajor[-1]
        elif minval >= 0:
            oppaxis = tslmajor[0]

        return tslmajor, oppaxis, step

    def _find_axes(self, minval, maxval):
        """
        Try to find axis tick positions that are well spaced
        """
        if minval >= maxval:
            raise ValueError, "maxval >= minval is required"

        # If there is a small differences between max and min values
        # compared to the size of the largest of (abs(maxval), abs(minval))
        # the function 'tasteless_ticks' is used, which in common usage is rare.
        if (abs((maxval - minval)/float(max(abs(maxval),abs(minval)))) < 0.2):
            tslmajor, oppaxis, step = self._tasteless_ticks(minval, maxval, 10)
        else:
            tslmajor, oppaxis, step = self._tasteful_ticks(minval, maxval)
        min = tslmajor[0] - step
        tslminor = sage.misc.misc.srange(min, maxval + 0.2*step, 0.2*step)
        tslminor = self._in_range(tslminor, minval, maxval)
        return oppaxis, step, tslminor, tslmajor

    def _draw_axes(self, subplot, axes, xmin, xmax, ymin, ymax, x_axis_ypos, y_axis_xpos):
        from matplotlib import lines
        if isinstance(axes, (list, tuple)) and len(axes) == 2 and \
        (axes[0] in [True, False]) and (axes[1] in [True, False]):
            self.__draw_x_axis = axes[0]
            self.__draw_y_axis = axes[1]
            #draw the x-axes?
            if self.__draw_x_axis:
                subplot.add_line(lines.Line2D([xmin, xmax], [x_axis_ypos, x_axis_ypos],
                                 color=self.__color, linewidth=float(self.__linewidth)))
            #draw y axis line?
            if self.__draw_y_axis:
                subplot.add_line(lines.Line2D([y_axis_xpos, y_axis_xpos],[ymin, ymax],
                                  color=self.__color, linewidth=float(self.__linewidth)))
        else: #draw them both
            subplot.add_line(lines.Line2D([xmin, xmax], [x_axis_ypos, x_axis_ypos],
                             color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(lines.Line2D([y_axis_xpos, y_axis_xpos],[ymin, ymax],
                              color=self.__color, linewidth=float(self.__linewidth)))

    def _draw_axes_labels(self, subplot, axes_labels, xmin, xmax, ymin, ymax, xstep, ystep, x_axis_ypos, y_axis_xpos, pad=0.2):
        al = axes_labels
        if not isinstance(al, (list,tuple)) or len(al) != 2:
            raise TypeError, "axes_labels must be a list of two strings."
        #draw x-axis label if there is a x-axis:
        fontsize = int(self.__fontsize)
        if self.__draw_x_axis:
            s = str(al[0])
            subplot.text(xmax + pad*xstep, x_axis_ypos, s,
                         fontsize = fontsize,
                         color = self.__axes_label_color, horizontalalignment="left",
                         verticalalignment="center", family="monospace")
            xmax += 0.0025*(xmax-xmin)*len(s) * fontsize

        #draw y-axis label if there is a y-axis
        if self.__draw_y_axis:
            subplot.text(y_axis_xpos, ymax + 2*pad*ystep, str(al[1]),
                         fontsize = fontsize,
                         color = self.__axes_label_color,
                         horizontalalignment="left", verticalalignment="center", family="monospace")
            ymax += 0.075*(ymax-ymin)
        return xmin, xmax, ymin, ymax


    def add_xy_axes(self, subplot, xmin, xmax, ymin, ymax, axes=True,
                    ticks="automatic", axesstyle="automatic", axes_labels=None):
        r"""
        \code{_add_xy_axes} is used when the 'save' method
        of any Graphics object is called.

        Additionally this function uses the function '_find_axes'
        from axis.py which attempts to find aesthetically pleasing
        tick and label spacing values.

        Some definitons of the parameters:

        y_axis_xpos : "where on the x-axis to draw the y-axis"
        xstep : "the spacing between major tick marks"
        xtslminor : "x-axis minor tick step list"
        xtslmajor : "x-axis major tick step list"
        yltheight : "where the top of the major ticks go"
        ystheight : "where the top of the minor ticks go"
        ylabel : "where the ylabel is drawn"
        xlabel : "where the xlabel is drawn"

        """
        from matplotlib import lines
        xmin = float(xmin); xmax=float(xmax); ymin=float(ymin); ymax=float(ymax)
        yspan = ymax - ymin
        xspan = xmax - xmin

        #evalute find_axes for x values and y ticks
        y_axis_xpos, xstep, xtslminor, xtslmajor = self._find_axes(xmin, xmax)
        yltheight = 0.015*xspan
        ystheight = 0.25*yltheight
        ylabel = y_axis_xpos - 2*ystheight

        #evalute find_axes for y values and x ticks
        x_axis_ypos, ystep, ytslminor, ytslmajor = self._find_axes(ymin, ymax)
        xltheight = 0.015*yspan
        xstheight = 0.25*xltheight
        xlabel = x_axis_ypos - xltheight

        if axes:
            self._draw_axes(subplot, axes, xmin, xmax, ymin, ymax, x_axis_ypos, y_axis_xpos)

        #the x-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for x in xtslmajor:
            if x == y_axis_xpos:
                continue
            s = self._format_tick_string(x)
            subplot.text(x, xlabel, s, fontsize=int(self.__fontsize), horizontalalignment="center",
                        color=self.__tick_label_color, verticalalignment="top")
            subplot.add_line(lines.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xltheight],
                        color=self.__color, linewidth=float(self.__linewidth)))

        #now draw the x-axis minor tick marks
        for x in xtslminor:
            subplot.add_line(lines.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xstheight],
                        color=self.__color, linewidth=float(self.__linewidth)))

        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for y in ytslmajor:
            if y == x_axis_ypos:
                continue
            s = self._format_tick_string(y)
            subplot.text(ylabel, y, s, fontsize=int(self.__fontsize), verticalalignment="center",
                        color=self.__tick_label_color, horizontalalignment="right")
            subplot.add_line(lines.Line2D([y_axis_xpos, y_axis_xpos + yltheight], [y, y],
                    color=self.__color, linewidth=float(self.__linewidth)))

        #now draw the x-axis minor tick marks
        for y in ytslminor:
            subplot.add_line(lines.Line2D([y_axis_xpos, y_axis_xpos + ystheight], [y, y],
                color=self.__color, linewidth=float(self.__linewidth)))

        # now draw the x and y axis labels
        if self.__axes_labels:
            xmin, xmax, ymin, ymax = self._draw_axes_labels(subplot, self.__axes_labels, xmin, xmax, ymin, ymax, xstep, ystep,
                                   x_axis_ypos, y_axis_xpos, pad=0.2)

        return xmin, xmax, ymin, ymax


    def _draw_frame(self, subplot, xmins, xmaxs, ymins, ymaxs):
        """
        Draw a frame around a graphic at the given
        (scaled out) x and y min and max values.
        """
        from matplotlib import lines
        #border horizontal axis:
        #bottom:
        subplot.add_line(lines.Line2D([xmins, xmaxs], [ymins, ymins],
                                        color=self.__color, linewidth=float(self.__linewidth)))
        #top:
        subplot.add_line(lines.Line2D([xmins, xmaxs], [ymaxs, ymaxs],
                                        color=self.__color, linewidth=float(self.__linewidth)))
        #border vertical axis:
        #left:
        subplot.add_line(lines.Line2D([xmins, xmins], [ymins, ymaxs],
                                        color=self.__color, linewidth=float(self.__linewidth)))
        #right:
        subplot.add_line(lines.Line2D([xmaxs, xmaxs], [ymins, ymaxs],
                                        color=self.__color, linewidth=float(self.__linewidth)))


    def add_xy_frame_axes(self, subplot, xmin, xmax, ymin, ymax,
                          axes_with_no_ticks=False, axes_labels=None):
        r"""
        Draw a frame around the perimeter of a graphic.

        Only major tick marks are drawn on a frame axes.

        If \code{axes_with_no_ticks} is true, then also draw
        centered axes with no tick marks.

        """
        from matplotlib import lines
        xmin = float(xmin); xmax=float(xmax); ymin=float(ymin); ymax=float(ymax)
        yspan = ymax - ymin
        xspan = xmax - xmin

        #evalute find_axes for x values and y ticks
        y_axis_xpos, xstep, xtslminor, xtslmajor = self._find_axes(xmin, xmax)
        yltheight = 0.015 * xspan
        ystheight = 0.25  * yltheight
        #ylabel    = y_axis_xpos - 2*ystheight
        ylabel    = -2*ystheight

        #evalute find_axes for y values and x ticks
        x_axis_ypos, ystep, ytslminor, ytslmajor = self._find_axes(ymin, ymax)
        xltheight = 0.015 * yspan
        xstheight = 0.25  * xltheight
        #xlabel    = x_axis_ypos - xltheight
        xlabel    = -xltheight

        #scale the axes out from the actual plot
        xmins, xmaxs, ymins, ymaxs = self._adjustments_for_frame(xmin, xmax, ymin, ymax)

        #now draw the frame border:
        self._draw_frame(subplot, xmins, xmaxs, ymins, ymaxs)

        #these are the centered axes, like in regular plot, but with no ticks
        if axes_with_no_ticks:
            #the x axis line
            subplot.add_line(lines.Line2D([xmins, xmaxs], [x_axis_ypos, x_axis_ypos],
                                            color=self.__color, linewidth=float(self.__linewidth)))

            #the y axis line
            subplot.add_line(lines.Line2D([y_axis_xpos, y_axis_xpos],[ymins, ymaxs],
                                            color=self.__color, linewidth=float(self.__linewidth)))

        #the x-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for x in xtslmajor:
            s = self._format_tick_string(x)
            subplot.text(x, xlabel + ymins, s, fontsize=int(self.__fontsize),
                         horizontalalignment="center", verticalalignment="top")

        #now draw the x-axis minor tick marks
        for x in xtslminor:
            subplot.add_line(lines.Line2D([x, x], [ymins, xstheight + ymins],
                        color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(lines.Line2D([x, x], [ymaxs, ymaxs - xstheight],
                        color=self.__color, linewidth=float(self.__linewidth)))


        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for y in ytslmajor:
            s = self._format_tick_string(y)
            subplot.text(ylabel + xmins, y, s, fontsize=int(self.__fontsize),
                         verticalalignment="center", horizontalalignment="right")

        #now draw the x-axis minor tick marks
        for y in ytslminor:
            subplot.add_line(lines.Line2D([xmins, ystheight + xmins], [y, y],
                color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(lines.Line2D([xmaxs, xmaxs - ystheight], [y, y],
                color=self.__color, linewidth=float(self.__linewidth)))

    def _adjustments_for_frame(self, xmin, xmax, ymin, ymax):
        r"""
        Scale the axes out from the actual plot to accommodate a frame.

        INPUT:
            xmin, xmax, ymin, ymax -- numbers

        OUTPUT:
            xmin, xmax, ymin, ymax -- numbers

        TESTS:
            sage: from sage.plot.axes import Axes
            sage: Axes()._adjustments_for_frame(-10,40,10,35)
            (-11.0, 41.0, 9.5, 35.5)
        """
        xmin = float(xmin); xmax=float(xmax); ymin=float(ymin); ymax=float(ymax)
        yspan = ymax - ymin
        xspan = xmax - xmin
        ys = 0.02*yspan
        xs = 0.02*xspan
        ymin -= ys
        ymax += ys
        xmin -= xs
        xmax += xs
        return xmin, xmax, ymin, ymax

    def add_xy_matrix_frame_axes(self, subplot, xmin, xmax, ymin, ymax):
        """
        Draw a frame around a \code{matrix_plot}.

        The tick marks drawn on the frame correspond to
        the ith row and jth column of the matrix.

        """
        from matplotlib import lines
        xmax = int(xmax)
        ymax = int(ymax)

        #evalute find_axes for x values and y ticks
        # the > 14 is just for appearance improvement
        if xmax > 19:
            #get nice tick spacing and assure we get largest point
            tl, opax, step = self._tasteful_ticks(0, xmax)
            #print tl, opax, step
            xtl = self._tasteful_ticks(0, xmax)[0] + [xmax-1]
            xrm = [float(x+0.5) for x in xtl]
            xtslmajor = [int(n) for n in xtl]
        else:
            xtl = sage.misc.misc.srange(0, xmax)
            xrm = [float(x+0.5) for x in xtl]
            xtslmajor = [int(n) for n in xtl]
        yltheight = 0.015*xmax
        ystheight = 0.25*yltheight
        ylabel = -2*ystheight

        #evalute find_axes for y values and x ticks
        # the > 14 is just for appearance improvement
        if ymax > 19:
            #get nice tick spacing and assure we get largest point
            tl, opax, step = self._tasteful_ticks(0, ymax)
            yrm = [0.5]+[float(y+0.5+1) for y in tl[1:-1]]+[ymax-0.5]
            yrm.reverse()
            ytslmajor = [0] + tl[1:-1] + [ymax-1]
        else:
            ytl = sage.misc.misc.srange(0, ymax)
            ytslmajor = [int(n) for n in ytl]
            yrm = [float(y+0.5) for y in ytslmajor]
            yrm.reverse()
        xltheight = 0.015*ymax
        xstheight = 0.25*xltheight
        xlabel = -xltheight

        #scale the axes out from the actual plot
        xs = 0.02*xmax
        ys = 0.02*ymax
        xmins = -xs
        xmaxs = xmax + xs
        ymins = -ys
        ymaxs = ymax + ys

        #now draw the frame border:
        self._draw_frame(subplot, xmins, xmaxs, ymins, ymaxs)

        #the x-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for xr, x in zip(xrm, xtslmajor):
            s = self._format_tick_string(x)
            subplot.text(xr, xlabel + ymins, s, fontsize=int(self.__fontsize),
                         horizontalalignment="center", verticalalignment="top")
            subplot.text(xr, -2*xlabel + ymaxs, s, fontsize=int(self.__fontsize),
                         horizontalalignment="center", verticalalignment="top")
            subplot.add_line(lines.Line2D([xr, xr], [ymins, xstheight + ymins],
                        color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(lines.Line2D([xr, xr], [ymaxs, ymaxs - xstheight],
                        color=self.__color, linewidth=float(self.__linewidth)))

        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for yr, y in zip(yrm, ytslmajor):
            s = self._format_tick_string(y)
            subplot.text(ylabel + xmins, yr, s, fontsize=int(self.__fontsize),
                         verticalalignment="center", horizontalalignment="right")
            subplot.text(-2*ylabel + xmaxs, yr, s, fontsize=int(self.__fontsize),
                         verticalalignment="center", horizontalalignment="left")
            subplot.add_line(lines.Line2D([xmins, ystheight + xmins], [yr, yr],
                color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(lines.Line2D([xmaxs, xmaxs - ystheight], [yr, yr],
                color=self.__color, linewidth=float(self.__linewidth)))

class GridLines(SageObject):
    """
    Grid lines for SAGE 2D Graphics.

    See the docstring for Graphics.show for examples.
    """
    def __init__(self, gridlines=None, gridlinesstyle=None,
            vgridlinesstyle=None, hgridlinesstyle=None):
        r"""
        Add horizontal and vertical grid lines to a Graphics object.

        INPUT:
            gridlines    -- (default: None) can be any of the following:
                            1. None, False: do not add grid lines.
                            2. True, "automatic", "major": add grid lines
                               at major ticks of the axes.
                            3. "minor": add grid at major and minor ticks.
                            4. [xlist,ylist]: a tuple or list containing
                               two elements, where xlist (or ylist) can be
                               any of the following.
                               4a. None, False: don't add horizontal (or
                                   vertical) lines.
                               4b. True, "automatic", "major": add
                                   horizontal (or vertical) grid lines at
                                   the major ticks of the axes.
                               4c. "minor": add horizontal (or vertical)
                                   grid lines at major and minor ticks of
                                   axes.
                               4d. an iterable yielding numbers n or pairs
                                   (n,opts), where n is the coordinate of
                                   the line and opt is a dictionary of
                                   MATPLOTLIB options for rendering the
                                   line.
            gridlinesstyle,
            hgridlinesstyle,
            vgridlinesstyle
                         -- (default: None) a dictionary of MATPLOTLIB
                            options for the rendering of the grid lines,
                            the horizontal grid lines or the vertical grid
                            lines, respectively.

        TESTS:
            sage: from sage.plot.axes import GridLines
            sage: GridLines()
            <class 'sage.plot.axes.GridLines'>
            sage: gl = GridLines(False)
            sage: gl = GridLines(True)
            sage: gl = GridLines("automatic")
            sage: gl = GridLines("major")
            sage: gl = GridLines("minor")
            sage: gl = GridLines([True,False])
            sage: gl = GridLines(["minor","major"])
            sage: gl = GridLines(["automatic",None])
            sage: gl = GridLines([range(-10,10,2), lambda x,y:srange(x,y,0.5)])
            sage: gl = GridLines(None, dict(color="red"),
            ...     dict(linestyle=":"), dict(color="blue"))
            sage: gl = GridLines(None, dict(rgbcolor="red"),
            ...     dict(linestyle=":"), dict(color="blue"))
        """
        self.__gridlines = gridlines

        defaultstyle = dict(color=(0.3,0.3,0.3),linewidth=0.4)
        if gridlinesstyle is not None:
            rgbcolor_keyword_support(gridlinesstyle)
            defaultstyle.update(gridlinesstyle)
        self.__gridlinesstyle = [copy(defaultstyle),copy(defaultstyle)]
        if vgridlinesstyle is not None:
            rgbcolor_keyword_support(vgridlinesstyle)
            self.__gridlinesstyle[0].update(vgridlinesstyle)
        if hgridlinesstyle is not None:
            rgbcolor_keyword_support(hgridlinesstyle)
            self.__gridlinesstyle[1].update(hgridlinesstyle)

    def add_gridlines(self, subplot, xmin, xmax, ymin, ymax, frame=False):
        # Process the input to get valid gridline data.
        r"""
        Add the grid lines to a subplot object.

        INPUT:
            subplot     -- an instance of matplotlib.axes.Subplot
            xmin, xmax  -- $x$ range of the Graphics object containing subplot
            ymin, ymax  -- $y$ range of the Graphics object containing subplot
            frame       -- (default: False) if True, then adjust the lengths of
                           the grid lines to touch connect to the frame.

        OUTPUT:
            None (modifies subplot)

        TESTS:
            sage: from sage.plot.axes import GridLines
            sage: from matplotlib.figure import Figure
            sage: subplot = Figure().add_subplot(111)
            sage: lims = [-10,20,10,35]

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines().add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            0

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines(False).add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            0

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines(True).add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            13

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines("automatic").add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            13

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines("major").add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            13

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines("minor").add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            57

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines([True,False]).add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            7

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines(["minor","major"]).add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            37

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines(["automatic",None]).add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            7

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines([range(-10,10,2), lambda x,y:srange(x,y,0.5)]).add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            60

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines("automatic", dict(color="red"),
            ...     dict(linestyle=":"),
            ...     dict(color="blue")).add_gridlines(subplot,*lims)
            sage: len(subplot.lines)
            13

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines([1,2,3]).add_gridlines(subplot,*lims)
            Traceback (most recent call last):
            ...
            TypeError: gridlines is not a list or tuple of length 2

            sage: subplot = Figure().add_subplot(111)
            sage: GridLines([1,2]).add_gridlines(subplot,*lims)
            Traceback (most recent call last):
            ...
            TypeError: elements of gridlines need to be iterable: [1, 2]
        """
        points = [[xmin, xmax], [ymin, ymax]]
        if self.__gridlines is None or self.__gridlines is False:
            return
        elif self.__gridlines in ["major", "automatic"] or self.__gridlines is True:
            self.__gridlines = [
                    self._get_ticks_locations(points[0]),
                    self._get_ticks_locations(points[1])
                    ]
        elif self.__gridlines == "minor":
            self.__gridlines = [
                    self._get_ticks_locations(points[0], self.__gridlines),
                    self._get_ticks_locations(points[1], self.__gridlines)
                    ]
        else:
            try:
                gridlines = [None]*2
                gridlines[0], gridlines[1] = self.__gridlines
            except ValueError:
                raise TypeError, "gridlines is not a list or tuple of length 2"

            for i in range(2):
                if gridlines[i] is None or gridlines[i] is False:
                    gridlines[i] = []
                elif gridlines[i] in ["major", "automatic"] or gridlines[i] is True:
                    gridlines[i] = self._get_ticks_locations(points[i])
                elif gridlines[i] == "minor":
                    gridlines[i] = self._get_ticks_locations(points[i],gridlines[i])
                elif callable(gridlines[i]):
                    gridlines[i] = gridlines[i](*points[i])

            if not (hasattr(gridlines[0],'__iter__') and
                    hasattr(gridlines[1],'__iter__')):
                raise TypeError, "elements of gridlines need to be iterable: %s" \
                    % gridlines
            self.__gridlines = gridlines

        # add the gridlines to subplot.
        if frame is True:
            xmin, xmax, ymin, ymax = \
                    self._get_adjustments_for_frame(*points)
            points = [[xmin, xmax], [ymin, ymax]]

        new_gridlines = []
        for i in range(2):
            new_list = []
            for entry in self.__gridlines[i]:
                kwds = copy(self.__gridlinesstyle[i])
                if hasattr(entry,'__len__'):
                    if isinstance(entry, (list, tuple)) and len(entry) == 2:
                        val = float(entry[0])
                        rgbcolor_keyword_support(entry[1])
                        kwds.update(entry[1])
                    else:
                        val = float(entry)
                        kwds = copy(self.__gridlinesstyle[i])
                else:
                    val = float(entry)
                    kwds = copy(self.__gridlinesstyle[i])
                new_list.append([val,kwds])
            new_gridlines.append(new_list)
        xlines, ylines = new_gridlines

        # draw the grid lines
        from matplotlib import lines
        # horizontal lines
        for (yval, ykwds) in ylines:
            subplot.add_line(
                    lines.Line2D(points[0],[yval,yval],**ykwds)
                    )
        # vertical lines
        for (xval, xkwds) in xlines:
            subplot.add_line(
                    lines.Line2D([xval,xval],points[1],**xkwds)
                    )

    def _get_ticks_locations(self, interval, ticks="major"):
        r"""
        Find the locations of the major and/or minor ticks of the axes
        in the interval.

        INPUT:
            interval -- an interval as a pair of numbers
            ticks -- "major" or "minor". If "minor", then return also
                the locations of the minor ticks.

        OUTPUT:
            list -- the locations of the ticks on the axes

        TESTS:
            sage: from sage.plot.axes import GridLines
            sage: GridLines()._get_ticks_locations([-10,20])
            [-10, -5, 0, 5, 10, 15, 20]
            sage: GridLines()._get_ticks_locations([10,35],"minor")
            [10.0, 11.0, 12.0, 13.0, ..., 32.0, 33.0, 34.0, 35.0]

        """
        # Axes._find_axes[2] returns locations of minor ticks
        # Axes._find_axes[3] returns locations of major ticks
        if ticks == "minor":
            minorticks = True
        else:
            minorticks = False
        return Axes()._find_axes(*interval)[2 if minorticks else 3]

    def _get_adjustments_for_frame(self, xinterval, yinterval):
        r"""
        Returns new limits for axes to accommodate a frame drawn around the
        plot.

        INPUT:
            xinterval -- x-axis interval as pairs of numbers
            yinterval -- y-axis interval as pairs of numbers

        OUTPUT:
            xmin, xmax, ymin, ymax -- numbers

        TESTS:
            sage: from sage.plot.axes import GridLines
            sage: GridLines()._get_adjustments_for_frame([-10,40],[10,35])
            (-11.0, 41.0, 9.5, 35.5)
        """
        return Axes()._adjustments_for_frame(*(xinterval+yinterval))

def rgbcolor_keyword_support(d):
    r"""
    Change the rgbcolor key to color.

    NOTE: The rgbcolor keyword is a synonym for color and not a matplotlib
    option.

    INPUT:
        d -- a dictionary

    OUTPUT:
        None -- modifies d

    TESTS:
        sage: from sage.plot.axes import rgbcolor_keyword_support
        sage: d = dict(rgbcolor="red",other="blue")
        sage: rgbcolor_keyword_support(d)
        sage: d
        {'color': 'red', 'other': 'blue'}

        sage: d = dict(color="red",other="blue")
        sage: rgbcolor_keyword_support(d)
        sage: d
        {'color': 'red', 'other': 'blue'}

        sage: d = dict(color="red", rgbcolor="red")
        sage: rgbcolor_keyword_support(d)
        sage: d
        {'color': 'red'}

        sage: d = dict(color="red", rgbcolor="black")
        sage: rgbcolor_keyword_support(d)
        Traceback (most recent call last):
        ...
        TypeError: specify only one of color or rgbcolor
    """
    if d.has_key("rgbcolor"):
        if d.has_key("color"):
            if d["rgbcolor"] == d["color"]:
                del d["rgbcolor"]
            else:
                raise TypeError, "specify only one of color or rgbcolor"
        else:
            d["color"] = d["rgbcolor"]
            del d["rgbcolor"]

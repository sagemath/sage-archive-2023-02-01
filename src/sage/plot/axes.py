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

class Axes(SageObject):
    """
    Axes for SAGE 2D Graphics.

    Set all axis properties and then add one of
    the following axes to the current (matplotlib) subplot:
    add_xy__axes
    add_xy_frame_axes
    add_xy_matrix_frame_axes

    """
    def __init__(self, color=(0,0,0), fontsize=6, linewidth=0.6,axes_label=None,
                 axes_label_color=(0,0,0), tick_color=(0,0,0), tick_label_color=(0,0,0)):
        self.__color = color
        self.__tick_color = tick_color
        self.__tick_label_color = tick_label_color
        self.__axes_label = axes_label
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
        absmax, absmin = abs(minval), abs(maxval)
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
        import matplotlib.patches as patches
        if isinstance(axes, (list, tuple)) and len(axes) == 2 and \
        (axes[0] in [True, False]) and (axes[1] in [True, False]):
            self.__draw_x_axis = axes[0]
            self.__draw_y_axis = axes[1]
            #draw the x-axes?
            if self.__draw_x_axis:
                subplot.add_line(patches.Line2D([xmin, xmax], [x_axis_ypos, x_axis_ypos],
                                 color=self.__color, linewidth=float(self.__linewidth)))
            #draw y axis line?
            if self.__draw_y_axis:
                subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos],[ymin, ymax],
                                  color=self.__color, linewidth=float(self.__linewidth)))
        else: #draw them both
            subplot.add_line(patches.Line2D([xmin, xmax], [x_axis_ypos, x_axis_ypos],
                             color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos],[ymin, ymax],
                              color=self.__color, linewidth=float(self.__linewidth)))

    def _draw_axes_labels(self, subplot, axes_label, xmax, ymax, xstep, ystep, x_axis_ypos, y_axis_xpos, pad=0.2):
        al = axes_label
        if not isinstance(al, (list,tuple)) or len(al) != 2:
            raise TypeError, "axes_label must be a list of two strings."
        #draw x-axis label if there is a x-axis:
        if self.__draw_x_axis:
            subplot.text(xmax + pad*xstep, x_axis_ypos, str(al[0]), fontsize=int(self.__fontsize),
                       color=self.__axes_label_color, horizontalalignment="center", verticalalignment="center")
        #draw y-axis label if there is a y-axis
        if self.__draw_x_axis:
            subplot.text(y_axis_xpos, ymax + pad*ystep, str(al[1]), fontsize=int(self.__fontsize),
                        color=self.__axes_label_color, horizontalalignment="center", verticalalignment="center")


    def add_xy_axes(self, subplot, xmin, xmax, ymin, ymax, axes=True,
                    ticks="automatic", axesstyle="automatic", axes_label=None):
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
        import matplotlib.patches as patches
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
                        color=self.__tick_color, verticalalignment="top")
            subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xltheight],
                        color=self.__color, linewidth=float(self.__linewidth)))

        #now draw the x-axis minor tick marks
        for x in xtslminor:
            subplot.add_line(patches.Line2D([x, x], [x_axis_ypos, x_axis_ypos + xstheight],
                        color=self.__color, linewidth=float(self.__linewidth)))

        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for y in ytslmajor:
            if y == x_axis_ypos:
                continue
            s = self._format_tick_string(y)
            subplot.text(ylabel, y, s, fontsize=int(self.__fontsize), verticalalignment="center",
                        color=self.__tick_color, horizontalalignment="right")
            subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + yltheight], [y, y],
                    color=self.__color, linewidth=float(self.__linewidth)))

        #now draw the x-axis minor tick marks
        for y in ytslminor:
            subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos + ystheight], [y, y],
                color=self.__color, linewidth=float(self.__linewidth)))

        #now draw the x and y axis labels
        if self.__axes_label:
            self._draw_axes_labels(subplot, self.__axes_label, xmax, ymax, xstep, ystep,
                                   x_axis_ypos, y_axis_xpos, pad=0.2)


    def _draw_frame(self, subplot, xmins, xmaxs, ymins, ymaxs):
        """
        Draw a frame around a graphic at the given
        (scaled out) x and y min and max values.
        """
        #border horizontal axis:
        #bottom:
        import matplotlib.patches as patches
        subplot.add_line(patches.Line2D([xmins, xmaxs], [ymins, ymins],
                                        color=self.__color, linewidth=float(self.__linewidth)))
        #top:
        subplot.add_line(patches.Line2D([xmins, xmaxs], [ymaxs, ymaxs],
                                        color=self.__color, linewidth=float(self.__linewidth)))
        #border vertical axis:
        #left:
        subplot.add_line(patches.Line2D([xmins, xmins], [ymins, ymaxs],
                                        color=self.__color, linewidth=float(self.__linewidth)))
        #right:
        subplot.add_line(patches.Line2D([xmaxs, xmaxs], [ymins, ymaxs],
                                        color=self.__color, linewidth=float(self.__linewidth)))


    def add_xy_frame_axes(self, subplot, xmin, xmax, ymin, ymax,
                          axes_with_no_ticks=False, axes_label=None):
        r"""
        Draw a frame around the perimeter of a graphic.

        Only major tick marks are drawn on a frame axes.

        If \code{axes_with_no_ticks} is true, then also draw
        centered axes with no tick marks.

        """
        import matplotlib.patches as patches
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
        ys = 0.02*yspan
        xs = 0.02*xspan
        ymins = ymin - ys
        ymaxs = ymax + ys
        xmins = xmin - xs
        xmaxs = xmax + xs

        #now draw the frame border:
        self._draw_frame(subplot, xmins, xmaxs, ymins, ymaxs)

        #these are the centered axes, like in regular plot, but with no ticks
        if axes_with_no_ticks:
            #the x axis line
            subplot.add_line(patches.Line2D([xmins, xmaxs], [x_axis_ypos, x_axis_ypos],
                                            color=self.__color, linewidth=float(self.__linewidth)))

            #the y axis line
            subplot.add_line(patches.Line2D([y_axis_xpos, y_axis_xpos],[ymins, ymaxs],
                                            color=self.__color, linewidth=float(self.__linewidth)))

        #the x-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for x in xtslmajor:
            s = self._format_tick_string(x)
            subplot.text(x, xlabel + ymins, s, fontsize=int(self.__fontsize),
                         horizontalalignment="center", verticalalignment="top")

        #now draw the x-axis minor tick marks
        for x in xtslminor:
            subplot.add_line(patches.Line2D([x, x], [ymins, xstheight + ymins],
                        color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(patches.Line2D([x, x], [ymaxs, ymaxs - xstheight],
                        color=self.__color, linewidth=float(self.__linewidth)))


        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for y in ytslmajor:
            s = self._format_tick_string(y)
            subplot.text(ylabel + xmins, y, s, fontsize=int(self.__fontsize),
                         verticalalignment="center", horizontalalignment="right")

        #now draw the x-axis minor tick marks
        for y in ytslminor:
            subplot.add_line(patches.Line2D([xmins, ystheight + xmins], [y, y],
                color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(patches.Line2D([xmaxs, xmaxs - ystheight], [y, y],
                color=self.__color, linewidth=float(self.__linewidth)))


    def add_xy_matrix_frame_axes(self, subplot, xmin, xmax, ymin, ymax):
        """
        Draw a frame around a \code{matrix_plot}.

        The tick marks drawn on the frame correspond to
        the ith row and jth column of the matrix.

        """
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
            subplot.add_line(patches.Line2D([xr, xr], [ymins, xstheight + ymins],
                        color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(patches.Line2D([xr, xr], [ymaxs, ymaxs - xstheight],
                        color=self.__color, linewidth=float(self.__linewidth)))

        #the y-axis ticks and labels
        #first draw major tick marks and their corresponding values
        for yr, y in zip(yrm, ytslmajor):
            s = self._format_tick_string(y)
            subplot.text(ylabel + xmins, yr, s, fontsize=int(self.__fontsize),
                         verticalalignment="center", horizontalalignment="right")
            subplot.text(-2*ylabel + xmaxs, yr, s, fontsize=int(self.__fontsize),
                         verticalalignment="center", horizontalalignment="left")
            subplot.add_line(patches.Line2D([xmins, ystheight + xmins], [yr, yr],
                color=self.__color, linewidth=float(self.__linewidth)))
            subplot.add_line(patches.Line2D([xmaxs, xmaxs - ystheight], [yr, yr],
                color=self.__color, linewidth=float(self.__linewidth)))


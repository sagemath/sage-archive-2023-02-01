r"""nodoctest
Use Matplotlib to draw graphics.

Matplotlib is included with \SAGE, but the version built with \sage
doesn't have any GUI built by default.  Instead, we take care to make
sure the backends (for producing png, eps, svg, etc.)  are built by
default.  This is much easier -- getting matplotlib's GUI's to build
in any configuration is nontrivial and greatly enhances the dependency
requirements for \SAGE.  Supporting the backends though is easy, and is
all that we need for, e.g., producing graphs that will be embedded in
the web browser interface to \SAGE.

We spent a lot of time deciding between matplotlib and gnuplot as the
canonical ``included with SAGE'' graphing system.  Gnuplot lost
primarily because of its non-GPL compatible license, which basically
disqualified it from further consideration for being in the core of
SAGE.  Matplotlib produces very beautiful graphics, so I'm quite happy
to go with it though, and it's nice that the code is all Python so
very readable.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

_pylab = None
def pylab():
    global _pylab
    if _pylab is None:
        try:
            import pylab
        except ImportError:
            raise RuntimeError, "You must have the optional matplotlib package installed to use the matplotlib functions.  See http://modular.ucsd.edu/SAGE/packages/optional/."
        _pylab = pylab
    return _pylab

def _adapt(args, kwds):
    from sage.structure.element import RingElement
    args = list(args)
    for i in range(len(args)):
        if isinstance(args[i], rings.RingElement):
            args[i] = float(args[i])
    args = tuple(args)
    for k in kwds:
        if isinstance(kwds[k], rings.RingElement):
            kwds[k] = float(kwds[k])
    return args, kwds

def plot(xvals, yvals, *args, **kwds):
    args, kwds = _adapt(args, kwds)
    return pylab().plot(xvals, yvals, *args, **kwds)

def plot_function(f, xvals, *args, **kwds):
    """
    Draw a plot of the function f evaluated at the values in xvals.

    To use this command you should run \sage with the command
    \code{sage -mpl}.  Also the optional matplotlib library (and
    dependencies) must be installed.

    INPUT:
        f -- a function of a singular argument
        xvals -- the values that will appear along the horizontal axis,
                 and at which f will be evaluated.
        *args, **kwds -- all additional parameters are passed along
                 to the plot command.
    OUTPUT:
        a plot of f is returned and displayed on-screen.

    EXAMPLES:
        sage: p = mpl.plot_function(lambda x : x*x, range(10))

    The arange command is similar to range but allows for a third
    step argument.
        sage: p = mpl.plot_function(sin, arange(0,7,0.1))
    """
    args, kwds = _adapt(args, kwds)
    xvals = [float(x) for x in xvals]
    yvals = [float(f(x)) for x in xvals]
    return plot(xvals, yvals, *args, **kwds)

def plot_array(array, *args, **kwds):
    args, kwds = _adapt(args, kwds)
    return pylab().imshow(array, *args, **kwds)



###########################

def arange(start, stop=None, step=1, typecode=None):
    """
    Just like range() except it returns an array whose type can be
    specified by the keyword argument typecode.
    """
    return pylab().arange(start, stop, step, typecode)


def subplot(nrows=1, ncols=1, number=0):
    """
    EXAMPLES:
        sage: mpl.subplot(2,2,0)
        sage: mpl.plot(sin, arange(0,10,0.1))
        sage: mpl.subplot(2,2,1)
        sage: mpl.plot(cos, arange(0,30,0.1))
        sage: mpl.subplot(2,2,2)
        sage: mpl.plot(tan, arange(0,30,0.1))
        sage: mpl.subplot(2,2,3)
        sage: mpl.plot(lambda x: x, arange(0,30,0.01))
    """
    return pylab().subplot(int(nrows), int(ncols), int(number+1))

def set_ranges(xmin, xmax, ymin=None, ymax=None):
    if ymin == None:
        _, _, ymin, ymax = pylab().axis()
    return pylab().axis([float(xmin), float(xmax), float(ymin), float(ymax)])

def get_ranges():
    return pylab().axis()

def axis():
    return pylab().gca()

def new(number=0, *args, **kwds):
    args, kwds = _adapt(args, kwds)
    return pylab().figure(int(number+1), *args, **kwds)

def select(number=0, *args, **kwds):
    args, kwds = _adapt(args, kwds)
    return pylab().figure(int(number+1), *args, **kwds)

def clear(n=None):
    if not (n is None):
        select(n)
    pylab().clf()

def clear_axis(n=None):
    if not (n is None):
        select_figure(n)
    pylab().cla()

def save(filename, *args, **kwds):
    """
    Format determined by filename include
        ps, eps, jpg, png
    """
    args, kwds = _adapt(args, kwds)
    pylab().savefig(str(filename), *args, **kwds)


###########################

def update():
    pylab().draw()

def text(x, y, string, *args, **kwds):
    args, kwds = _adapt(args, kwds)
    pylab().text(x, y, string, *args, **kwds)
    update()

def set_tick_labels(axes='x', labels=[], *args, **kwds):
    args, kwds = _adapt(args, kwds)
    if axes == 'x':
        axis().set_xticklabels(labels)
    elif axes == 'y':
        axis().set_yticklabels(labels)
    else:
        raise ValueError, "Unknown axes %s"%axes
    update()

def axes_label(axes='x', label='x', *args, **kwds):
    args, kwds = _adapt(args, kwds)
    if axes == 'x':
        pylab().xlabel(label, *args, **kwds)
    elif axes == 'y':
        pylab().ylabel(label, *args, **kwds)
    else:
        raise ValueError, "Unknown axes %s"%axes
    update()

def title(string, *args, **kwds):
    args, kwds = _adapt(args, kwds)
    pylab().title(string, *args, **kwds)
    update()

def _axis(axes='x'):
    if axes == 'x':
        return axis().xaxis
    elif axes == 'y':
        return axis().yaxis
    else:
        raise ValueError, "axis %s not defined"%axis

def set_tick_locations(axes='x', locator='auto', *args, **kwds):
    """
    INPUT:
        axes -- either 'x' or 'y'
        locator -- 'null' -- no ticks
                   'index' -- locator for index plots
                   'linear' -- evenly spaced ticks from min to max
                   'log', base=10.0, subs=[1.0] -- logarithmic ticks from min to max
                                   at base^i * subs[j]
                   'multiple' -- ticks and range are a multiple of base
                   'auto' (default) -- choose a multiple locator and dynamically reassign
    """
    args, kwds = _adapt(args, kwds)
    A = _axis(axes)
    A.set_major_locator(pylab().matplotlib.ticker.__dict__[\
                     '%sLocator'%locator.capitalize()](*args))
    update()




def set_tick_formatter(axes='x', format='scalar', *args, **kwds):
    """
    INPUT:
        axes -- either 'x' or 'y'
        format -- 'null' -- no labels on the ticks
                  'fixed' -- set the strings manually for the labels
                  'func'  -- user defined function sets the labels:
                                argument is f(x, pos=0), which returns
                                the string for tick val x at position pos
                  'formatStr' -- (note capital S), use a sprintf format string
                  'index' -- cycle through fixed strings by tick position
                  'scalar' -- (default); autopick
                  'log' -- formatter for log axes
                  'date' -- use a strftime string to format the date
    """
    args, kwds = _adapt(args, kwds)
    A = _axis(axes)
    A.set_major_formatter(pylab().matplotlib.ticker.__dict__[\
                     '%sFormatter'%format.capitalize()](*args))
    update()

def set_tick_format(axes='x', format_str='%d'):
    """
    INPUT:
        axes -- either 'x' or 'y'
        format_str -- string, the format string
    """
    A = _axis(axes)
    A.set_major_formatter(pylab().matplotlib.ticker.FormatStrFormatter(format_str))
    update()


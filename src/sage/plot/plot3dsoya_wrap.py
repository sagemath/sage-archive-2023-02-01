"""


AUTHORS:
    -- William Stein
    -- Josh Kantor
"""

##############################################################################
#       Copyright (C) 2006 Josh Kantor and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
##############################################################################


def plot3dsoya(f, p=(0,0), side=1.0, res = 8, fineness=25,
               ignore_bad_values=True):
    r"""
    Return a 3d plot of the values of the function f(x,y) on
    the \emph{square} centered at p with given side length.
    Use show on the output to view it.

    WARNING:
        1. You \emph{must} have the optional soya3d package.
        2. Sometimes if soya crashes you will have to explicitly
           kill the soya3d window.

    KEYBOARD CONTROLS:
       ...
       (and printout when plot appears)

    INPUT:
        f -- object that can be called with two arguments as input, i.e.,
                so that f(x,y) makes sense and returns a value that is
                coercible to a float.
        p -- (default: (0,0)) 2-tuple (x,y) that defined the point
                the plot is centered at.
        side -- (default: 1.0) float, length of the sides of the square whose
                image we are plotting
        res --  (default: 8), integer, linear resolution, i.e., a lower bound
                on the number of sample points in the x or y direction.  Note
                that internally res is rounded up to the nearest power of 2.
        fineness -- (default: 25), integer, specifies the smoothness of the
                color map, i.e., the number of shades as the color map varies
                from blue (=small) to red (=big).
        ignore_bad_values -- (default: True), bool; if True any values of f(x,y)
                are replaced by 0.

    OUTPUT:
        A 3d plot object with a show method.

    EXAMPLES:

    A periodic surface:
        sage: p = plot3dsoya(lambda x,y : sin(x)+cos(y), (0,0), 10)   # optional: requires soya3d

    Now use the show command to view the surface:
        sage.: p.show()

    The Riemann Zeta function
        p = plot3dsoya(lambda x,y : abs(zeta(x + y*I)), (0.5,5), 10, res=4)

    """
    try:
        p = (float(p[0]),float(p[1]))
    except (IndexError, KeyError, TypeError):
        raise TypeError, "p (=%s) must be a point, i.e., a 2-tuple or list of length 2"%p

    try:
        side = float(side)
    except TypeError:
        raise TypeError, "side (=%s) must be coercible to a float"%side

    try:
        res = float(res)
    except TypeError:
        raise TypeError, "res (=%s) must be coercible to a float"%res

    try:
        fineness = int(fineness)
    except TypeError:
        raise TypeError, "fineness (=%s) must be an integer."%fineness




    import plot3dsoya
    return plot3dsoya.plot3dsoya(f, p, side/2, res, fineness, ignore_bad_values)

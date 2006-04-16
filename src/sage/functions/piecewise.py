import sage.plot.plot
import sage.interfaces.all

class Piecewise:
    """
    Implements a very simple class of piecewise defined functions.

    TODO: Lots!

    AUTHOR: David Joyner (2006-04)

    """
    def __init__(self, list_of_pairs):
        r"""
        \code{list_of_pairs} is a list of pairs (fcn,I), where fcn is
        a SAGE function (such as a polynomial over RR, or functions
        using the lambda notation), and I is an interval such as I =
        interval(1,3).

        We assume that these definitions are consistent. We assume
        that x is the variable used.
        """
        self._length = len(list_of_pairs)
        self._intervals = [x[0] for x in list_of_pairs]
        self._functions = [x[1] for x in list_of_pairs]
        self._list = list_of_pairs

    def list(self):
        return self._list

    def length(self):
        return self._length

    def _repr_(self):
        return 'Piecewise defined functions with %s parts, %s'%(
            self.length(),self.list())

    def intervals(self):
        return self._intervals

    def functions(self):
        return self._functions


## FIX THIS -- it isn't right as is, with
##             the p[1](x), which doesn't make any sense (what is x?)
##     def integral(self):
##         """
##         Returns the definite integral (as computed by maxima)
##         $\sum_I \int_I self|_I$, as I runs over the intervals
##         belonging to self.

##         EXAMPLES:
##             sage: f1 = lambda x:1
##             sage: f2 = lambda x:1-x
##             sage: f = Piecewise([[interval(0,1),f1],[interval(1,2),f2]])
##             sage: f.integral()
##             1/2
##         """
##         maxima = sage.interfaces.all.maxima
##         ints = [maxima('%s'%p[1](x)).integral('x', p[0][0], p[0][1]) \
##                 for p in self.list()]
##         return add(ints)

    def plot(self):
        """
        Returns the plot of self.

        EXAMPLES:
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f = Piecewise([[interval(0,1),f1],[interval(1,2),f2]])
            sage: P = f.plot()
        """
        plot = sage.plot.plot.plot
        plots = [plot(p[1], p[0][0], p[0][1], rgbcolor=(1,0,0) ) \
                 for p in self.list()]
        return plots


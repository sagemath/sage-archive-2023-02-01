"""
Generic plotting interface

(Not done -- in progress)


TODO:
   -- create a list of primitive, probably based on graph.py, which is based on
      pstricks, which is based on postscript.
   -- make classes for all primitives, like Path below.
   -- make matplotlib, gnuplot, and latex functions
   -- add plot methods that use these to primitives to many sage objects, e.g.,
      elliptic curves, transcedental functions, riemann zeta (on complex values),
      modular forms, l-series, matrices
   -- worry about where output appears (e.g., on screen with matplotlib, in web browser,
      in html/dvi/pdf log, etc.)
   -- investigate viewer sage workspace with pygame.
"""


from sage.ext.sage_object import SageObject

#import matplotlib
#matplotlib.use('Agg')

#import pylab
#from sage.interfaces.all import gnuplot

import sage.misc.all as misc

class Plot(SageObject):
    def __init__(self, desc, png=None):
        self.__desc = desc
        self.__png = png

    def _repr_(self):
        return str(self.__desc)

    def png(self, filename):
        if self.__png is None:
            raise NotImplementedError
        import os
        os.system('cp %s %s'%(self.__png, filename))

class PlotPrimitive(SageObject):
    def gnuplot(self):
        raise NotImplementedError

    def latex(self):    # pstricks
        raise NotImplementedError

    def matplotlib(self):
        raise NotImplementedError


class Path(PlotPrimitive):
    def __init__(self, points):
        self.__points = points

    def _repr_(self):
        return "A path"

    def matplotlib(self):
        points = self.__points
        x_list = [float(P[0]) for P in points]
        y_list = [float(P[1]) for P in points]
        pylab.plot(x_list, y_list)

    def gnuplot(self):
        fname = misc.SAGE_TMP + '/gnuplot.data'
        O = open(fname,'w')
        for x,y in self.__points:
            O.write('%f\t%f\n'%(x,y))
        O.close()
        cmd = 'plot "%s"'%fname
        gnuplot.plot(cmd)

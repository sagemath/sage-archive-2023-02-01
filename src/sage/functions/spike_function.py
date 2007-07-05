

import math

from sage.plot.all import line
from sage.modules.free_module_element import vector
from sage.rings.all import RDF

class SpikeFunction:
    """
    INPUT:
        v -- list of pairs (x, height)
        eps -- parameter that determines approximation to a true spike

    OUTPUT:
        a function with spikes at each point in v with the given height.
    """
    def __init__(self, v, eps=0.0000001):
        if len(v) == 0:
           v = [(0,0)]
        v = [(float(x[0]), float(x[1])) for x in v]
        notify = False

        for i in reversed(range(len(v)-1)):
            if v[i+1][0] - v[i][0] <= eps:
                notify = True
                del v[i+1]

        if notify:
            print "Some overlapping spikes have been deleted."
            print "You might want to use a smaller value for eps."

        v.sort()
        self.v = v
        self.eps = eps
        self.support = [float(x[0]) for x in self.v]
        self.height = [float(x[1]) for x in self.v]

    def __repr__(self):
        return "A spike function with spikes at %s"%self.support

    def _eval(self, x):
        eps = self.eps
        x = float(x)
        for i in range(len(self.support)):
            z = self.support[i]
            if z - eps <= x and x <= z + eps:
                return self.height[i], i
        return float(0), -1

    def __call__(self, x):
        return self._eval(x)[0]

    def plot_fft(self, samples=2**12, xmin=None, xmax=None,  **kwds):
        w = self.vector(samples = samples, xmin=xmin, xmax=xmax)
        xmin, xmax = self._ranges(xmin, xmax)
        z = w.fft()
        k = vector(RDF, [abs(z[i]) for i in range(int(len(z)/2))])
        return k.plot(xmin=0, xmax=xmax, **kwds)

    def vector(self, samples=2**16, xmin=None, xmax=None):
        v = vector(RDF, samples)
        xmin, xmax = self._ranges(xmin, xmax)
        delta = (xmax - xmin)/samples
        w = int(math.ceil(self.eps/delta))
        for i in range(len(self.support)):
            x = self.support[i]
            if x > xmax:
                break
            h = self.height[i]
            j = int((x - xmin)/delta)
            for k in range(j, min(samples, j+w)):
                v[k] = h
        return v

    def _ranges(self, xmin, xmax):
        width = (self.support[-1] + self.support[0])/float(2)
        if xmin is None:
           xmin = self.support[0] - width/float(5)
        if xmax is None:
           xmax = self.support[-1] + width/float(5)
        if xmax <= xmin:
           xmax = xmin + 1
        return xmin, xmax

    def plot(self, xmin=None, xmax=None, **kwds):
        v = []
        xmin, xmax = self._ranges(xmin, xmax)
        x = xmin
        eps = self.eps
        while x < xmax:
            y, i = self._eval(x)
            v.append( (x, y) )
            if i != -1:
                x0 = self.support[i] + eps
                v.extend([(x0,y), (x0,0)])
                if i+1 < len(self.support):
                    x = self.support[i+1] - eps
                    v.append( (x, 0) )
                else:
                    x = xmax
                    v.append( (xmax, 0) )
            else:
                new_x = None
                for j in range(len(self.support)):
                    if self.support[j] - eps > x:
                        new_x = self.support[j] - eps
                        break
                if new_x is None:
                    new_x = xmax
                v.append( (new_x, 0) )
                x = new_x
        L = line(v, **kwds)
        L.xmin(xmin-1); L.xmax(xmax)
        return L



spike_function = SpikeFunction

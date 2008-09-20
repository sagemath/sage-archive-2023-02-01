############################################################################
#  Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

"""
Animated plots

EXAMPLES:
We plot a circle shooting up to the right:

    sage: a = animate([circle((i,i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
    ...               xmin=0,ymin=0,xmax=2,ymax=2,figsize=[2,2]) # optional -- requires convert command
    sage: a.show() # optional -- requires convert command

"""

import os
import shutil

from sage.structure.sage_object import SageObject

import plot
import sage.misc.misc
import sage.misc.viewer

class Animation(SageObject):
    """
    Return an animation of a sequence of plots of objects.

    INPUT:
        v -- list of SAGE objects. These should preferably be graphics objects,
             but if they aren't then plot is called on them.
        xmin, xmax, ymin, ymax -- the ranges of the x and y axis.
        **kwds -- all additional inputs are passed onto the rendering
              command.  E.g., use figsize to adjust the resolution and
              aspect ratio.

    EXAMPLES:
        sage: a = animate([sin(x + float(k)) for k in srange(0,4,0.3)],
        ...                xmin=0, xmax=2*pi, figsize=[2,1]) # optional -- requires convert command
        sage: a # optional
        Animation with 14 frames
        sage: a[:5] # optional
        Animation with 5 frames
        sage: a.show()          # optional
        sage: a[:5].show()      # optional

    We draw an animation of drawing a parabola:
        sage: step = 0.1
        sage: L = Graphics()
        sage: v = []
        sage: for i in srange(0,1,step):
        ...       L += line([(i,i^2),(i+step,(i+step)^2)], rgbcolor=(1,0,0), thickness=2)
        ...       v.append(L)
        ...
        sage: a = animate(v, xmin=0, ymin=0)
        sage: a.show()
        sage: show(L)

    TESTS:
    This illustrates ticket \#2066 is fixed (setting axes ranges when an endpoint is 0):
        sage: animate([plot(sin, -1,1)], xmin=0, ymin=0)._Animation__xmin
        0
    """
    def __init__(self, v,
                 xmin=None, xmax=None, ymin=None, ymax=None,
                 **kwds):
        w = []
        for x in v:
            if not isinstance(x, plot.Graphics):
                if xmin is None:
                    xmin = 0
                if xmax is None:
                    xmax = 1
                x = plot.plot(x, (xmin, xmax))
            w.append(x)
        if len(w) == 0:
            w = [plot.Graphics()]
        self.__frames = w
        G = w[0]
        if xmin is None:
            xmin = G.xmin()
        if xmax is None:
            xmax = G.xmax()
        if ymin is None:
            ymin = G.ymin()
        if ymax is None:
            ymax = G.ymax()
        self.__xmin = xmin
        self.__xmax = xmax
        self.__ymin = ymin
        self.__ymax = ymax
        self.__kwds = kwds

        self._set_axes()

    def __getitem__(self, i):
        """
        Get a frame from an animation.

        EXAMPLES:
            sage: a = animate([x, x^2, x^3, x^4]) # optional -- requires convert command
            sage: a[2].show()       # optional
        """
        return self.__frames[i]

    def __getslice__(self, *args):
        """
        Slice this animation returning a subanimation.

        EXAMPLES:
            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ...               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2]) # optional -- requires convert command
            ...
            sage: a # optional
            Animation with 10 frames
            sage: a.show() # optional
            sage: a[3:7] # optional
            Animation with 4 frames
            sage: a[3:7].show() # optional
        """
        return Animation(self.__frames.__getslice__(*args), xmin=self.__xmin,
                       xmax = self.__xmax, ymin = self.__ymin,
                       ymax = self.__ymax, **self.__kwds)

    def _set_axes(self):
        for F in self.__frames:
            F.xmin(self.__xmin)
            F.xmax(self.__xmax)
            F.ymin(self.__ymin)
            F.ymax(self.__ymax)

    def _repr_(self):
        return "Animation with %s frames"%(len(self.__frames))

    def __add__(self, other):
        """
        Add two animations.  This has the effect of
        superimposing the two animinations frame-by-frame.

        EXAMPLES:
        We add and multiply two animations.

            sage: a = animate([circle((i,0),1) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1]) # optional -- requires convert command
            sage: a.show()   # optional
            sage: b = animate([circle((0,i),1,hue=0) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2]) # optional
            sage: b.show() # optional
            sage: (a*b).show()    # optional
            sage: (a+b).show()    # optional
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = dict(self.__kwds)
        for k, v in other.__kwds.iteritems():
            if not kwds.has_key(k):
                kwds[k] = v

        if len(self.__frames) > len(other.__frames):
            frames = self.__frames
            for i in range(len(other.__frames)):
                frames[i] += other.__frames[i]
        else:
            frames = other.__frames
            for i in range(len(self.__frames)):
                frames[i] += self.__frames[i]

        return Animation(frames,
                       xmin = min(self.__xmin, other.__xmin),
                       xmax = max(self.__xmax, other.__xmax),
                       ymin = min(self.__ymin, other.__ymin),
                       ymax = max(self.__ymax, other.__ymax),
                       **kwds)

    def __mul__(self, other):
        """
        Multiply two animations.  This has the effect of
        appending the two animinations (the second comes
        after the first).

        EXAMPLES:
        We add and multiply two animations.
            sage: a = animate([circle((i,0),1,thickness=20*i) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1], axes=False) # optional -- requires convert command
            sage: a.show()     # optional
            sage: b = animate([circle((0,i),1,hue=0,thickness=20*i) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2], axes=False) # optional
            sage: b.show()             # optional
            sage: (a*b).show()         # optional
            sage: (a+b).show()         # optional
        """
        if not isinstance(other, Animation):
            other = Animation(other)
        kwds = dict(self.__kwds)
        for k, v in other.__kwds.iteritems():
            if not kwds.has_key(k):
                kwds[k] = v
        return Animation(self.__frames + other.__frames,
                       xmin = min(self.__xmin, other.__xmin),
                       xmax = max(self.__xmax, other.__xmax),
                       ymin = min(self.__ymin, other.__ymin),
                       ymax = max(self.__ymax, other.__ymax),
                       **kwds)

    def png(self, dir=None):
        """
        Return the absolute path to a temp directory that contains the rendered
        png's of all the images in this animation.

        EXAMPLES:
            sage: a = animate([plot(x^2 + n) for n in range(4)]) # optional -- requires convert command
            sage: d = a.png() # optional  -- directory where the files are
            sage: v = os.listdir(d); v.sort(); v # optional
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
        """
        try:
            return self.__png_dir
        except AttributeError:
            pass
        d = sage.misc.misc.tmp_dir()
        xmin = self.__xmin
        xmax = self.__xmax
        ymin = self.__ymin
        ymax = self.__ymax
        G = self.__frames
        for i in range(len(G)):
            filename = '%s/%s'%(d,sage.misc.misc.pad_zeros(i,8))
            G[i].save(filename + '.png',
                      xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, **self.__kwds)
        self.__png_dir = d
        return d

    def graphics_array(self, ncols=3):
        """
        Return a graphics array with the given number of columns
        with plots of the frames of this animation.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: v = [E.change_ring(GF(p)).plot(pointsize=30) for p in [97, 101, 103, 107]]
            sage: a = animate(v, xmin=0, ymin=0) # optional -- requires convert command
            sage: a # optional
            Animation with 4 frames
            sage: a.show()        # optional

            sage: g = a.graphics_array() # optional
            sage: print g # optional
            Graphics Array of size 1 x 3
            sage: g.show(figsize=[4,1]) # optional

            sage: g = a.graphics_array(ncols=2) # optional
            sage: print g # optional
            Graphics Array of size 2 x 2
            sage: g.show('sage.png')         # optional
        """
        n = len(self.__frames)
        ncols = int(ncols)
        return plot.graphics_array(self.__frames, int(n/ncols),  ncols)

##     def html(self):
##         d = self.png()

##         for filename in os.path.listdir(d):
##             shutil.copyfile(d + '/' + filename, filename)

    def gif(self, delay=20, outfile=None, iterations=0):
        """
        Returns an animated gif composed from rendering the
        graphics objects in self.

        This function will only work if the Imagemagick command line
        tools package is installed, i.e., you have the ``\code{convert}'' command.

        INPUT:
            delay -- (default: 20) delay in hundredths of a second between frames
            outfile -- file that the animated gif gets saved to
            iterations -- integer (default: 0); number of iterations of
                          animation.  If 0, loop forever.

        AUTHOR:
            -- William Stein
        """
        if not outfile:
            outfile = sage.misc.misc.graphics_filename(ext='gif')
        if not outfile.endswith('.gif'):
            outfile += '.gif'
        outfile = os.path.abspath(outfile)
        d = self.png()
        cmd = 'cd "%s"; convert -delay %s -loop %s *.png "%s"'%(d, int(delay), int(iterations), outfile)
        os.system(cmd)

    def show(self, delay=20, iterations=0):
        """
        Show this animation.

        Currently this is done by default using an animated gif, though this
        could change in the future.
        """
        if plot.DOCTEST_MODE:
            self.gif(delay = delay, iterations = iterations)
            return

        if plot.EMBEDDED_MODE:
            self.gif(delay = delay, iterations = iterations)
        else:
            filename = sage.misc.misc.tmp_filename() + '.gif'
            self.gif(delay=delay, outfile=filename, iterations=iterations)
            os.system('%s %s 2>/dev/null 1>/dev/null &'%(
                sage.misc.viewer.browser(), filename))

    def save(self, filename=None):
        if filename is None or filename.endswith('.gif'):
            self.gif(outfile=filename)
            return
        elif filename.endswith('.sobj'):
            SageObject.save(self, filename)
            return
        else:
            raise ValueError, "Unable to save to a file with the extension '%s'"%(
                os.path.splitext(filename)[1][1:])






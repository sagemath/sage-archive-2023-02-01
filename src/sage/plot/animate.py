############################################################################
#  Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

"""
Animated plots

EXAMPLES:
We plot a circle shooting up to the right::

    sage: a = animate([circle((i,i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
    ...               xmin=0,ymin=0,xmax=2,ymax=2,figsize=[2,2])
    sage: a.show() # optional -- requires convert command

"""

import os
import shutil

from sage.structure.sage_object import SageObject

import plot
import sage.misc.misc
import sage.misc.viewer

class Animation(SageObject):
    r"""
    Return an animation of a sequence of plots of objects.

    INPUT:


    -  ``v`` - list of Sage objects. These should
       preferably be graphics objects, but if they aren't then plot is
       called on them.

    -  ``xmin, xmax, ymin, ymax`` - the ranges of the x and
       y axes.

    -  ``**kwds`` - all additional inputs are passed onto
       the rendering command. E.g., use figsize to adjust the resolution
       and aspect ratio.


    EXAMPLES::

        sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.3)],
        ...                xmin=0, xmax=2*pi, figsize=[2,1])
        sage: a
        Animation with 21 frames
        sage: a[:5]
        Animation with 5 frames
        sage: a.show()          # optional -- requires convert command
        sage: a[:5].show()      # optional

    The ``show`` function takes arguments to specify the
    delay between frames (measured in hundredths of a second, default
    value 20) and the number of iterations (default value 0, which
    means to iterate forever). To iterate 4 times with half a second
    between each frame::

        sage: a.show(delay=50, iterations=4) # optional

    An animation of drawing a parabola::

        sage: step = 0.1
        sage: L = Graphics()
        sage: v = []
        sage: for i in srange(0,1,step):
        ...       L += line([(i,i^2),(i+step,(i+step)^2)], rgbcolor=(1,0,0), thickness=2)
        ...       v.append(L)
        sage: a = animate(v, xmin=0, ymin=0)
        sage: a.show() # optional -- requires convert command
        sage: show(L)

    TESTS: This illustrates that ticket #2066 is fixed (setting axes
    ranges when an endpoint is 0)::

        sage: animate([plot(sin, -1,1)], xmin=0, ymin=0)._Animation__kwds['xmin']
        0
    """
    def __init__(self, v, **kwds):
        w = []
        for x in v:
            if not isinstance(x, plot.Graphics):
                x = plot.plot(x, (kwds.get('xmin',0), kwds.get('xmax', 1)))
            w.append(x)
        if len(w) == 0:
            w = [plot.Graphics()]
        self.__frames = w
        self.__kwds = kwds

    def _combine_kwds(self, *kwds_tuple):
        """
        Returns a dictionary which is a combination of the all the
        dictionaries in kwds_tuple. This also does the appropriate thing
        for taking the mins and maxes of all of the x/y mins/maxes.

        EXAMPLES::

            sage: a = animate([plot(sin, -1,1)], xmin=0, ymin=0)
            sage: kwds1 = {'a':1, 'b':2, 'xmin':2, 'xmax':5}
            sage: kwds2 = {'b':3, 'xmin':0, 'xmax':4}
            sage: kwds = a._combine_kwds(kwds1, kwds2)
            sage: list(sorted(kwds.items()))
            [('a', 1), ('b', 3), ('xmax', 5), ('xmin', 0)]
        """
        new_kwds = {}

        for kwds in kwds_tuple:
            new_kwds.update(kwds)

        import __builtin__
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            values = [v for v in [kwds.get(name, None) for kwds in kwds_tuple] if v is not None]
            if values:
                new_kwds[name] = getattr(__builtin__, name[1:])(*values)
        return new_kwds


    def __getitem__(self, i):
        """
        Get a frame from an animation.

        EXAMPLES::

            sage: a = animate([x, x^2, x^3, x^4])
            sage: a[2].show()       # optional -- requires convert command
        """
        return self.__frames[i]

    def __getslice__(self, *args):
        """
        Slice this animation returning a subanimation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ...               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: a
            Animation with 10 frames
            sage: a.show() # optional -- requires convert command
            sage: a[3:7]
            Animation with 4 frames
            sage: a[3:7].show() # optional
        """
        return Animation(self.__frames.__getslice__(*args), **self.__kwds)

    def _repr_(self):
        """
        Print representation for an animation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ...               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: a
            Animation with 10 frames
            sage: a._repr_()
            'Animation with 10 frames'
        """
        return "Animation with %s frames"%(len(self.__frames))

    def __add__(self, other):
        """
        Add two animations. This has the effect of superimposing the two
        animations frame-by-frame.

        EXAMPLES: We add and multiply two animations.

        ::

            sage: a = animate([circle((i,0),1) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1])
            sage: a.show()   # optional -- requires convert command
            sage: b = animate([circle((0,i),1,hue=0) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2])
            sage: b.show() # optional
            sage: (a*b).show()    # optional
            sage: (a+b).show()    # optional
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self.__kwds, other.__kwds)

        #Combine the frames
        m = max(len(self.__frames), len(other.__frames))
        frames = [a+b for a,b in zip(self.__frames, other.__frames)]
        frames += self.__frames[m:] + other.__frames[m:]

        return Animation(frames, **kwds)

    def __mul__(self, other):
        """
        Multiply two animations. This has the effect of appending the two
        animations (the second comes after the first).

        EXAMPLES: We add and multiply two animations.

        ::

            sage: a = animate([circle((i,0),1,thickness=20*i) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1], axes=False)
            sage: a.show()     # optional -- requires convert command
            sage: b = animate([circle((0,i),1,hue=0,thickness=20*i) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2], axes=False)
            sage: b.show()             # optional
            sage: (a*b).show()         # optional
            sage: (a+b).show()         # optional
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self.__kwds, other.__kwds)

        return Animation(self.__frames + other.__frames, **kwds)

    def png(self, dir=None):
        """
        Return the absolute path to a temp directory that contains the
        rendered png's of all the images in this animation.

        EXAMPLES::

            sage: a = animate([plot(x^2 + n) for n in range(4)])
            sage: d = a.png()
            sage: v = os.listdir(d); v.sort(); v
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
        """
        try:
            return self.__png_dir
        except AttributeError:
            pass
        d = sage.misc.misc.tmp_dir()
        G = self.__frames
        for i, frame in enumerate(self.__frames):
            filename = '%s/%s'%(d,sage.misc.misc.pad_zeros(i,8))
            frame.save(filename + '.png', **self.__kwds)
        self.__png_dir = d
        return d

    def graphics_array(self, ncols=3):
        """
        Return a graphics array with the given number of columns with plots
        of the frames of this animation.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: v = [E.change_ring(GF(p)).plot(pointsize=30) for p in [97, 101, 103, 107]]
            sage: a = animate(v, xmin=0, ymin=0)
            sage: a
            Animation with 4 frames
            sage: a.show()        # optional -- requires convert command

        ::

            sage: g = a.graphics_array()
            sage: print g
            Graphics Array of size 1 x 3
            sage: g.show(figsize=[4,1]) # optional

        ::

            sage: g = a.graphics_array(ncols=2)
            sage: print g
            Graphics Array of size 2 x 2
            sage: g.show('sage.png')         # optional
        """
        n = len(self.__frames)
        ncols = int(ncols)
        return plot.graphics_array(self.__frames, int(n/ncols),  ncols)

    def gif(self, delay=20, savefile=None, iterations=0, show_path=False):
        r"""
        Returns an animated gif composed from rendering the graphics
        objects in self.

        This function will only work if the ImageMagick command line tools
        package is installed, i.e., you have the "``convert``"
        command.

        INPUT:


        -  ``delay`` - (default: 20) delay in hundredths of a
           second between frames

        -  ``savefile`` - file that the animated gif gets saved
           to

        -  ``iterations`` - integer (default: 0); number of
           iterations of animation. If 0, loop forever.

        -  ``show_path`` - boolean (default: False); if True,
           print the path to the saved file


        If savefile is not specified: in notebook mode, display the
        animation; otherwise, save it to a default file name.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: a.gif()              # optional -- requires convert command
            sage: a.gif(delay=35, iterations=3)       # optional
            sage: a.gif(savefile='my_animation.gif')  # optional
            sage: a.gif(savefile='my_animation.gif', show_path=True) # optional
            Animation saved to .../my_animation.gif.

        .. note::

           If ImageMagick is not installed, you will get an error
           message like this::

              /usr/local/share/sage/local/bin/sage-native-execute: 8: convert:
              not found

              Error: ImageMagick does not appear to be installed. Saving an
              animation to a GIF file or displaying an animation requires
              ImageMagick, so please install it and try again.

           See www.imagemagick.org, for example.

        AUTHORS:

        - William Stein
        """
        if not savefile:
            savefile = sage.misc.misc.graphics_filename(ext='gif')
        if not savefile.endswith('.gif'):
            savefile += '.gif'
        savefile = os.path.abspath(savefile)
        d = self.png()
        cmd = 'cd "%s"; sage-native-execute convert -delay %s -loop %s *.png "%s"'%(d, int(delay), int(iterations), savefile)
        from subprocess import check_call, CalledProcessError
        try:
            check_call(cmd, shell=True)
            if show_path:
                print "Animation saved to file %s." % savefile
        except (CalledProcessError, OSError):
            print ""
            print "Error: ImageMagick does not appear to be installed. Saving an"
            print "animation to a GIF file or displaying an animation requires"
            print "ImageMagick, so please install it and try again."
            print ""
            print "See www.imagemagick.org, for example."

    def show(self, delay=20, iterations=0):
        r"""
        Show this animation.

        INPUT:


        -  ``delay`` - (default: 20) delay in hundredths of a
           second between frames

        -  ``iterations`` - integer (default: 0); number of
           iterations of animation. If 0, loop forever.


        .. note::

           Currently this is done using an animated gif, though this
           could change in the future. This requires that the ImageMagick
           command line tools package be installed, i.e., that you have the
           ``convert`` command.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: a.show()       # optional -- requires convert command

        The preceding will loop the animation forever. If you want to show
        only three iterations instead::

            sage: a.show(iterations=3)    # optional

        To put a half-second delay between frames::

            sage: a.show(delay=50)        # optional

        .. note::

           If ImageMagick is not installed, you will get an error
           message like this::

              /usr/local/share/sage/local/bin/sage-native-execute: 8: convert:
              not found

              Error: ImageMagick does not appear to be installed. Saving an
              animation to a GIF file or displaying an animation requires
              ImageMagick, so please install it and try again.

            See www.imagemagick.org, for example.
        """
        if plot.DOCTEST_MODE:
            self.gif(delay = delay, iterations = iterations)
            return

        if plot.EMBEDDED_MODE:
            self.gif(delay = delay, iterations = iterations)
        else:
            filename = sage.misc.misc.tmp_filename() + '.gif'
            self.gif(delay=delay, savefile=filename, iterations=iterations)
            os.system('%s %s 2>/dev/null 1>/dev/null &'%(
                sage.misc.viewer.browser(), filename))

    def save(self, filename=None, show_path=False):
        """
        Save this animation into a gif or sobj file.

        INPUT:


        -  ``filename`` - (default: None) name of save file

        -  ``show_path`` - boolean (default: False); if True,
           print the path to the saved file


        If filename is None, then in notebook mode, display the animation;
        othewise, save the animation to a gif file. If filename ends in
        '.gif', save to a gif file. If filename ends in '.sobj', save to an
        sobj file.

        In all other cases, print an error.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: a.save()         # optional -- requires convert command
            sage: a.save('wave.gif')   # optional
            sage: a.save('wave.gif', show_path=True)   # optional
            Animation saved to file .../wave.gif.
            sage: a.save('wave0.sobj')  # optional
            sage: a.save('wave1.sobj', show_path=True)  # optional
            Animation saved to file wave1.sobj.
        """
        if filename is None or filename.endswith('.gif'):
            self.gif(savefile=filename, show_path=show_path)
            return
        elif filename.endswith('.sobj'):
            SageObject.save(self, filename)
            if show_path:
                print "Animation saved to file %s." % filename
            return
        else:
            raise ValueError, "Unable to save to a file with the extension '%s'"%(
                os.path.splitext(filename)[1][1:])

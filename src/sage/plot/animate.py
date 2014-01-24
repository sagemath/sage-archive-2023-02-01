r"""
Animated plots

Animations are generated from a list (or other iterable) of graphics
objects.  Images are produced by calling the ``save_image`` method on
each input object, and using ImageMagick's ``convert`` program [IM] or
``ffmpeg`` [FF] to generate an animation.  The output format is GIF by
default, but can be any of the formats supported by ``convert`` or
``ffmpeg``.

.. Warning::

    Note that ImageMagick and FFmpeg are not included with Sage, and
    must be installed by the user.  On unix systems, type ``which
    convert`` at a command prompt to see if ``convert`` (part of the
    ImageMagick suite) is installed.  If it is, you will be given its
    location.  Similarly, you can check for ``ffmpeg`` with ``which
    ffmpeg``.  See [IM] or [FF] for installation instructions.

EXAMPLES:

The sine function::

    sage: sines = [plot(c*sin(x), (-2*pi,2*pi), color=Color(c,0,0), ymin=-1,ymax=1) for c in sxrange(0,1,.2)]
    sage: a = animate(sines)
    sage: a
    Animation with 5 frames
    sage: a.show()  # optional -- ImageMagick

Animate using ffmpeg instead of ImageMagick::

    sage: f = sage.misc.temporary_file.tmp_filename(ext='.gif')
    sage: a.save(filename=f,use_ffmpeg=True) # optional -- ffmpeg

An animated :class:`sage.plot.graphics.GraphicsArray` of rotating ellipses::

    sage: E = animate((graphics_array([[ellipse((0,0),a,b,angle=t,xmin=-3,xmax=3)+circle((0,0),3,color='blue') for a in range(1,3)] for b in range(2,4)]) for t in sxrange(0,pi/4,.15)))
    sage: E         # animations produced from a generator do not have a known length
    Animation with unknown number of frames
    sage: E.show()  # optional -- ImageMagick

A simple animation of a circle shooting up to the right::

    sage: c = animate([circle((i,i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
    ....:               xmin=0,ymin=0,xmax=2,ymax=2,figsize=[2,2])
    sage: c.show() # optional -- ImageMagick


Animations of 3d objects::

    sage: var('s,t')
    (s, t)
    sage: def sphere_and_plane(x):
    ....:     return sphere((0,0,0),1,color='red',opacity=.5)+parametric_plot3d([t,x,s],(s,-1,1),(t,-1,1),color='green',opacity=.7)
    sage: sp = animate([sphere_and_plane(x) for x in sxrange(-1,1,.3)])
    sage: sp[0]      # first frame
    sage: sp[-1]     # last frame
    sage: sp.show()  # optional -- ImageMagick

    sage: (x,y,z) = var('x,y,z')
    sage: def frame(t):
    ....:     return implicit_plot3d((x^2 + y^2 + z^2), (x, -2, 2), (y, -2, 2), (z, -2, 2), plot_points=60, contour=[1,3,5], region=lambda x,y,z: x<=t or y>=t or z<=t)
    sage: a = animate([frame(t) for t in srange(.01,1.5,.2)])
    sage: a[0]       # first frame
    sage: a.show()   # optional -- ImageMagick

If the input objects do not have a ``save_image`` method, then the
animation object attempts to make an image by calling its internal
method :meth:`sage.plot.animate.Animation.make_image`.  This is
illustrated by the following example::

    sage: t = var('t')
    sage: a = animate((sin(c*pi*t) for c in sxrange(1,2,.2)))
    sage: a.show()  # optional -- ImageMagick


AUTHORS:

- William Stein
- John Palmieri
- Niles Johnson (2013-12): Expand to animate more graphics objects

REFERENCES:

.. [IM] http://www.imagemagick.org

.. [FF]  http://www.ffmpeg.org

"""

############################################################################
#  Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

import os

from sage.structure.sage_object import SageObject
from sage.misc.temporary_file import tmp_filename, tmp_dir
import plot
import sage.misc.misc
import sage.misc.viewer


def animate(frames, **kwds):
    r"""
    Animate a list of frames by creating a
    :class:`sage.plot.animate.Animation` object.

    EXAMPLES::

        sage: t = var('t')
        sage: a = animate((cos(c*pi*t) for c in sxrange(1,2,.2)))
        sage: a.show()  # optional -- ImageMagick
    """
    return Animation(frames, **kwds)

class Animation(SageObject):
    r"""
    Return an animation of a sequence of plots of objects.

    INPUT:


    - ``v`` - iterable of Sage objects. These should preferably be
      graphics objects, but if they aren't then :meth:`make_image` is
      called on them.

    - ``xmin, xmax, ymin, ymax`` - the ranges of the x and y axes.

    - ``**kwds`` - all additional inputs are passed onto the rendering
      command. E.g., use figsize to adjust the resolution and aspect
      ratio.


    EXAMPLES::

        sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.3)],
        ....:                xmin=0, xmax=2*pi, figsize=[2,1])
        sage: a
        Animation with 21 frames
        sage: a[:5]
        Animation with 5 frames
        sage: a.show()          # optional -- ImageMagick
        sage: a[:5].show()      # optional -- ImageMagick

    The :meth:`show` method takes arguments to specify the
    delay between frames (measured in hundredths of a second, default
    value 20) and the number of iterations (default value 0, which
    means to iterate forever). To iterate 4 times with half a second
    between each frame::

        sage: a.show(delay=50, iterations=4) # optional -- ImageMagick

    An animation of drawing a parabola::

        sage: step = 0.1
        sage: L = Graphics()
        sage: v = []
        sage: for i in srange(0,1,step):
        ....:       L += line([(i,i^2),(i+step,(i+step)^2)], rgbcolor=(1,0,0), thickness=2)
        ....:       v.append(L)
        sage: a = animate(v, xmin=0, ymin=0)
        sage: a.show() # optional -- ImageMagick
        sage: show(L)

    TESTS:

    This illustrates that :trac:`2066` is fixed (setting axes
    ranges when an endpoint is 0)::

        sage: animate([plot(sin, -1,1)], xmin=0, ymin=0)._kwds['xmin']
        0

    We check that :trac:`7981` is fixed::

        sage: a = animate([plot(sin(x + float(k)), (0, 2*pi), ymin=-5, ymax=5)
        ....:            for k in srange(0,2*pi,0.3)])
        sage: a.show() # optional -- ImageMagick

    Do not convert input iterator to a list::

        sage: a = animate((x^p for p in sxrange(1,2,.1))); a
        Animation with unknown number of frames
        sage: a._frames
        <generator object ...

    """
    def __init__(self, v=None, **kwds):
        r"""
        Return an animation of a sequence of plots of objects.  See
        documentation of :func:`animate` for more details and
        examples.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.3)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1]) # indirect doctest
            sage: a
            Animation with 21 frames
        """
        self._frames = v
        self._kwds = kwds

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

        Test that the bug reported in :trac:`12107` has been fixed::

            sage: kwds3 = {}
            sage: kwds4 = {'b':3, 'xmin':0, 'xmax':4}
            sage: a._combine_kwds(kwds3, kwds4)['xmin']
            0
        """
        new_kwds = {}

        for kwds in kwds_tuple:
            new_kwds.update(kwds)

        import __builtin__
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            values = [v for v in [kwds.get(name, None) for kwds in kwds_tuple] if v is not None]
            if values:
                new_kwds[name] = getattr(__builtin__, name[1:])(values)
        return new_kwds


    def __getitem__(self, i):
        """
        Get a frame from an animation or
        slice this animation returning a subanimation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ....:               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: a
            Animation with 10 frames
            sage: frame2 = a[2]  # indirect doctest
            sage: frame2.show()
            sage: a.show() # optional -- ImageMagick
            sage: a[3:7]   # indirect doctest
            Animation with 4 frames
            sage: a[3:7].show() # optional -- ImageMagick
        """
        if isinstance(i, slice):
            return Animation(self._frames[i], **self._kwds)
        else:
            return self._frames[i]

    def _repr_(self):
        """
        Print representation for an animation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ....:               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: a
            Animation with 10 frames
            sage: a._repr_()
            'Animation with 10 frames'
        """
        try:
            num = len(self)
        except TypeError:
            num = "unknown number of"
        return "Animation with %s frames"%num

    def __add__(self, other):
        """
        Add two animations. This has the effect of superimposing the two
        animations frame-by-frame.

        EXAMPLES::

            sage: a = animate([circle((i,0),1) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1])
            sage: a.show()        # optional -- ImageMagick
            sage: b = animate([circle((0,i),1,hue=0) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2])
            sage: b.show()        # optional -- ImageMagick
            sage: s = a+b         # indirect doctest
            sage: len(a), len(b)
            (5, 5)
            sage: len(s)
            5
            sage: s.show()        # optional -- ImageMagick
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self._kwds, other._kwds)

        #Combine the frames
        m = max(len(self._frames), len(other._frames))
        frames = [a+b for a,b in zip(self._frames, other._frames)]
        frames += self._frames[m:] + other._frames[m:]

        return Animation(frames, **kwds)

    def __mul__(self, other):
        """
        Multiply two animations. This has the effect of appending the two
        animations (the second comes after the first).

        EXAMPLES::

            sage: a = animate([circle((i,0),1,thickness=20*i) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1], axes=False)
            sage: a.show()             # optional -- ImageMagick
            sage: b = animate([circle((0,i),1,hue=0,thickness=20*i) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2], axes=False)
            sage: b.show()             # optional -- ImageMagick
            sage: p = a*b              # indirect doctest
            sage: len(a), len(b)
            (5, 5)
            sage: len(p)
            10
            sage: (a*b).show()         # optional -- ImageMagick
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self._kwds, other._kwds)

        return Animation(self._frames + other._frames, **kwds)

    def __len__(self):
        """
        Length of self

        EXAMPLES::

            sage: a = animate([circle((i,0),1,thickness=20*i) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1], axes=False)
            sage: len(a)
            5
        """
        return len(self._frames)

    def make_image(self, frame, filename, **kwds):
        r"""
        Given a frame which has no ``save_image()`` method, make a graphics
        object and save it as an image with the given filename.  By default, this is
        :meth:`sage.plot.plot.plot`.  To make animations of other objects,
        override this method in a subclass.

        EXAMPLES::

            sage: from sage.plot.animate import Animation
            sage: class MyAnimation(Animation):
            ....:    def make_image(self, frame, filename, **kwds):
            ....:        P = parametric_plot(frame[0], frame[1], **frame[2])
            ....:        P.save_image(filename,**kwds)

            sage: t = var('t')
            sage: x = lambda t: cos(t)
            sage: y = lambda n,t: sin(t)/n
            sage: B = MyAnimation([([x(t), y(i+1,t)],(t,0,1), {'color':Color((1,0,i/4)), 'aspect_ratio':1, 'ymax':1}) for i in range(4)])

            sage: d = B.png()
            sage: v = os.listdir(d); v.sort(); v
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
            sage: B.show()  # not tested

            sage: class MyAnimation(Animation):
            ....:    def make_image(self, frame, filename, **kwds):
            ....:        G = frame.plot()
            ....:        G.set_axes_range(floor(G.xmin()),ceil(G.xmax()),floor(G.ymin()),ceil(G.ymax()))
            ....:        G.save_image(filename, **kwds)

            sage: B = MyAnimation([graphs.CompleteGraph(n) for n in range(7,11)], figsize=5)
            sage: d = B.png()
            sage: v = os.listdir(d); v.sort(); v
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
            sage: B.show()  # not tested

        """
        p = plot.plot(frame)
        p.save_image(filename, **kwds)

    def png(self, dir=None):
        r"""
        Render PNG images of the frames in this animation, saving them
        in ``dir``.  Return the absolute path to that directory.  If
        the frames have been previously rendered and ``dir`` is
        ``None``, just return the directory in which they are stored.

        When ``dir`` is other than ``None``, force re-rendering of
        frames.

        INPUT:

        - ``dir`` -- Directory in which to store frames.  Default
          ``None``; in this case, a temporary directory will be
          created for storing the frames.

        EXAMPLES::

            sage: a = animate([plot(x^2 + n) for n in range(4)])
            sage: d = a.png()
            sage: v = os.listdir(d); v.sort(); v
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
        """
        if dir is None:
            try:
                return self._png_dir
            except AttributeError:
                pass
            d = tmp_dir()
        else:
            d = dir
        G = self._frames
        for i, frame in enumerate(self._frames):
            filename = '%s/%s'%(d,sage.misc.misc.pad_zeros(i,8))
            try:
                frame.save_image(filename + '.png', **self._kwds)
            except AttributeError:
                self.make_image(frame, filename + '.png', **self._kwds)
        self._png_dir = d
        return d

    def graphics_array(self, ncols=3):
        r"""
        Return a :class:`sage.plot.graphics.GraphicsArray` with plots of the
        frames of this animation, using the given number of columns.
        The frames must be acceptable inputs for
        :class:`sage.plot.graphics.GraphicsArray`.


        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: v = [E.change_ring(GF(p)).plot(pointsize=30) for p in [97, 101, 103, 107]]
            sage: a = animate(v, xmin=0, ymin=0)
            sage: a
            Animation with 4 frames
            sage: a.show() # optional -- ImageMagick

        Modify the default arrangement of array::

            sage: g = a.graphics_array(); print g
            Graphics Array of size 2 x 3
            sage: g.show(figsize=[4,1]) # optional

        Specify different arrangement of array and save with different file name::

            sage: g = a.graphics_array(ncols=2); print g
            Graphics Array of size 2 x 2
            sage: g.show('sage.png') # optional

        Frames can be specified as a generator too; it is internally converted to a list::

            sage: t = var('t')
            sage: b = animate((plot(sin(c*pi*t)) for c in sxrange(1,2,.2)))
            sage: g = b.graphics_array(); print g
            Graphics Array of size 2 x 3
            sage: g.show() # optional
        """
        ncols = int(ncols)
        frame_list = list(self._frames)
        n = len(frame_list)
        nrows, rem = divmod(n,ncols)
        if rem > 0:
            nrows += 1
        return plot.graphics_array(frame_list, nrows,  ncols)

    def gif(self, delay=20, savefile=None, iterations=0, show_path=False,
            use_ffmpeg=False):
        r"""
        Returns an animated gif composed from rendering the graphics
        objects in self.

        This method will only work if either (a) the ImageMagick
        software suite is installed, i.e., you have the ``convert``
        command or (b) ``ffmpeg`` is installed.  See
        [IM] for more about ImageMagick, and see
        [FF] for more about ``ffmpeg``.  By default, this
        produces the gif using ``convert`` if it is present.  If this
        can't find ``convert`` or if ``use_ffmpeg`` is True, then it
        uses ``ffmpeg`` instead.

        INPUT:

        -  ``delay`` - (default: 20) delay in hundredths of a
           second between frames

        -  ``savefile`` - file that the animated gif gets saved
           to

        -  ``iterations`` - integer (default: 0); number of
           iterations of animation. If 0, loop forever.

        -  ``show_path`` - boolean (default: False); if True,
           print the path to the saved file

        - ``use_ffmpeg`` - boolean (default: False); if True, use
          'ffmpeg' by default instead of 'convert'.

        If ``savefile`` is not specified: in notebook mode, display the
        animation; otherwise, save it to a default file name.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir()
            sage: a.gif()              # not tested
            sage: a.gif(savefile=dir + 'my_animation.gif', delay=35, iterations=3)  # optional -- ImageMagick
            sage: a.gif(savefile=dir + 'my_animation.gif', show_path=True) # optional -- ImageMagick
            Animation saved to .../my_animation.gif.
            sage: a.gif(savefile=dir + 'my_animation_2.gif', show_path=True, use_ffmpeg=True) # optional -- ffmpeg
            Animation saved to .../my_animation_2.gif.

        .. note::

           If neither ffmpeg nor ImageMagick is installed, you will
           get an error message like this::

              Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
              animation to a GIF file or displaying an animation requires one of these
              packages, so please install one of them and try again.

              See www.imagemagick.org and www.ffmpeg.org for more information.
        """
        from sage.misc.sage_ostools import have_program
        have_convert = have_program('convert')
        have_ffmpeg = self._have_ffmpeg()
        if use_ffmpeg or not have_convert:
            if have_ffmpeg:
                self.ffmpeg(savefile=savefile, show_path=show_path,
                            output_format='.gif', delay=delay,
                            iterations=iterations)
            else:
                if not have_convert:
                    msg = """
Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
animation to a GIF file or displaying an animation requires one of these
packages, so please install one of them and try again.

See www.imagemagick.org and www.ffmpeg.org for more information."""
                else:
                    msg = """
Error: ffmpeg does not appear to be installed.  Download it from
www.ffmpeg.org, or use 'convert' to produce gifs instead."""
                raise OSError, msg
        else:
            if not savefile:
                savefile = tmp_filename(ext='.gif')
            if not savefile.endswith('.gif'):
                savefile += '.gif'
            savefile = os.path.abspath(savefile)
            d = self.png()
            cmd = ( 'cd "%s"; sage-native-execute convert -dispose Background '
                    '-delay %s -loop %s *.png "%s"' ) % ( d, int(delay),
                        int(iterations), savefile )
            from subprocess import check_call, CalledProcessError
            try:
                check_call(cmd, shell=True)
                if show_path:
                    print "Animation saved to file %s." % savefile
            except (CalledProcessError, OSError):
                msg = """
Error: Cannot generate GIF animation.  Verify that convert
(ImageMagick) or ffmpeg is installed, and that the objects passed to
the animate command can be saved in PNG image format.

See www.imagemagick.org and www.ffmpeg.org for more information."""
                raise OSError, msg

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
           could change in the future. This requires that either
           ffmpeg or the ImageMagick suite (in particular, the
           ``convert`` command) is installed.

        See also the :meth:`ffmpeg` method.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: a.show()       # optional -- ImageMagick

        The preceding will loop the animation forever. If you want to show
        only three iterations instead::

            sage: a.show(iterations=3)    # optional -- ImageMagick

        To put a half-second delay between frames::

            sage: a.show(delay=50)        # optional -- ImageMagick

        .. note::

           If you don't have ffmpeg or ImageMagick installed, you will
           get an error message like this::

              Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
              animation to a GIF file or displaying an animation requires one of these
              packages, so please install one of them and try again.

              See www.imagemagick.org and www.ffmpeg.org for more information.
        """
        if sage.doctest.DOCTEST_MODE:
            filename = tmp_filename(ext='.gif')
            self.gif(savefile=filename, delay=delay, iterations=iterations)
            return

        if plot.EMBEDDED_MODE:
            self.gif(delay = delay, iterations = iterations)
        else:
            filename = tmp_filename(ext='.gif')
            self.gif(delay=delay, savefile=filename, iterations=iterations)
            os.system('%s %s 2>/dev/null 1>/dev/null &'%(
                sage.misc.viewer.browser(), filename))

    def _have_ffmpeg(self):
        """
        Return True if the program 'ffmpeg' is installed.  See
        www.ffmpeg.org to download ffmpeg.

        EXAMPLES::

            sage: a = animate([plot(sin, -1,1)], xmin=0, ymin=0)
            sage: a._have_ffmpeg() # random: depends on whether ffmpeg is installed
            False
        """
        from sage.misc.sage_ostools import have_program
        return have_program('ffmpeg')

    def ffmpeg(self, savefile=None, show_path=False, output_format=None,
               ffmpeg_options='', delay=None, iterations=0, pix_fmt='rgb24'):
        r"""
        Returns a movie showing an animation composed from rendering
        the frames in self.

        This method will only work if ffmpeg is installed.  See
        http://www.ffmpeg.org for information about ffmpeg.

        INPUT:

        -  ``savefile`` - file that the mpeg gets saved to.

        .. warning:

            This will overwrite ``savefile`` if it already exists.

        - ``show_path`` - boolean (default: False); if True, print the
          path to the saved file

        - ``output_format`` - string (default: None); format and
          suffix to use for the video.  This may be 'mpg', 'mpeg',
          'avi', 'gif', or any other format that ffmpeg can handle.
          If this is None and the user specifies ``savefile`` with a
          suffix, say ``savefile='animation.avi'``, try to determine the
          format ('avi' in this case) from that file name.  If no file
          is specified or if the suffix cannot be determined, 'mpg' is
          used.

        - ``ffmpeg_options`` - string (default: ''); this string is
          passed directly to ffmpeg.

        - ``delay`` - integer (default: None); delay in hundredths of a
          second between frames.  The framerate is 100/delay.
          This is not supported for mpeg files: for mpegs, the frame
          rate is always 25 fps.

        - ``iterations`` - integer (default: 0); number of iterations
          of animation. If 0, loop forever.  This is only supported
          for animated gif output and requires ffmpeg version 0.9 or
          later.  For older versions, set ``iterations=None``.

        - ``pix_fmt`` - string (default: 'rgb24'); used only for gif
          output.  Different values such as 'rgb8' or 'pal8' may be
          necessary depending on how ffmpeg was installed.  Set
          ``pix_fmt=None`` to disable this option.

        If ``savefile`` is not specified: in notebook mode, display
        the animation; otherwise, save it to a default file name.  Use
        :func:`sage.misc.misc.set_verbose` with ``level=1`` to see
        additional output.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir()
            sage: a.ffmpeg(savefile=dir + 'new.mpg')       # optional -- ffmpeg
            sage: a.ffmpeg(savefile=dir + 'new.avi')       # optional -- ffmpeg
            sage: a.ffmpeg(savefile=dir + 'new.gif')       # optional -- ffmpeg
            sage: a.ffmpeg(savefile=dir + 'new.mpg', show_path=True) # optional -- ffmpeg
            Animation saved to .../new.mpg.

        .. note::

           If ffmpeg is not installed, you will get an error message
           like this::

              Error: ffmpeg does not appear to be installed. Saving an animation to
              a movie file in any format other than GIF requires this software, so
              please install it and try again.

              See www.ffmpeg.org for more information.


        TESTS::

            sage: a.ffmpeg(output_format='gif',delay=30,iterations=5)     # optional -- ffmpeg
        """
        if not self._have_ffmpeg():
            msg = """Error: ffmpeg does not appear to be installed. Saving an animation to
a movie file in any format other than GIF requires this software, so
please install it and try again."""
            raise OSError, msg
        else:
            if savefile is None:
                if output_format is None:
                    output_format = '.mpg'
                else:
                    if output_format[0] != '.':
                        output_format = '.'+output_format
                savefile = tmp_filename(ext=output_format)
            else:
                if output_format is None:
                    suffix = os.path.splitext(savefile)[1]
                    if len(suffix) > 0:
                        output_format = suffix
                    else:
                        output_format = '.mpg'
            if not savefile.endswith(output_format):
                savefile += output_format
            early_options = ''
            if output_format == '.gif':
                # We try to set reasonable options for gif output.
                #
                # Older versions of ffmpeg (before 0.9, summer 2011)
                # use the option -loop_output instead of -loop.
                # Setting iterations=None is a way of preventing sage
                # from adding the -loop option.  A separate
                # -loop_output option can be added with the
                # ffmpeg_options argument.
                if iterations is not None:
                    loop_cmd = '-loop {0} '.format(iterations)
                else:
                    loop_cmd = ''
                # A pix_fmt value is required for some but not all
                # ffmpeg installations.  Setting pix_fmt=None will
                # prevent sage from adding this option, and it may be
                # controlled separately through ffmpeg_options.
                if pix_fmt is not None:
                    pix_fmt_cmd = '-pix_fmt {0} '.format(pix_fmt)
                else:
                    pix_fmt_cmd = ''
                ffmpeg_options += ' {0}{1}'.format(pix_fmt_cmd,loop_cmd)
            if delay is not None and output_format != '.mpeg' and output_format != '.mpg':
                early_options += ' -r %s ' % int(100/delay)
            savefile = os.path.abspath(savefile)
            pngdir = self.png()
            pngs = os.path.join(pngdir, "%08d.png")
            # For ffmpeg, it seems that some options, like '-g ... -r
            # ...', need to come before the input file names, while
            # some options, like '-pix_fmt rgb24', need to come
            # afterwards.  Hence 'early_options' and 'ffmpeg_options'
            cmd = 'cd "%s"; sage-native-execute ffmpeg -y -f image2 %s -i %s %s %s' % (pngdir, early_options, pngs, ffmpeg_options, savefile)
            from subprocess import check_call, CalledProcessError, PIPE
            try:
                if sage.misc.misc.get_verbose() > 0:
                    set_stderr = None
                else:
                    set_stderr = PIPE
                sage.misc.misc.verbose("Executing '%s'" % cmd,level=1)
                sage.misc.misc.verbose("\n---- ffmpeg output below ----\n")
                check_call(cmd, shell=True, stderr=set_stderr)
                if show_path:
                    print "Animation saved to file %s." % savefile
            except (CalledProcessError, OSError):
                print "Error running ffmpeg."
                raise

    def save(self, filename=None, show_path=False, use_ffmpeg=False):
        """
        Save this animation.

        INPUT:

        -  ``filename`` - (default: None) name of save file

        -  ``show_path`` - boolean (default: False); if True,
           print the path to the saved file

        - ``use_ffmpeg`` - boolean (default: False); if True, use
          'ffmpeg' by default instead of 'convert' when creating GIF
          files.

        If filename is None, then in notebook mode, display the
        animation; otherwise, save the animation to a GIF file.  If
        filename ends in '.sobj', save to an sobj file.  Otherwise,
        try to determine the format from the filename extension
        ('.mpg', '.gif', '.avi', etc.).  If the format cannot be
        determined, default to GIF.

        For GIF files, either ffmpeg or the ImageMagick suite must be
        installed.  For other movie formats, ffmpeg must be installed.
        An sobj file can be saved with no extra software installed.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir()
            sage: a.save()         # not tested
            sage: a.save(dir + 'wave.gif')   # optional -- ImageMagick
            sage: a.save(dir + 'wave.gif', show_path=True)   # optional -- ImageMagick
            Animation saved to file .../wave.gif.
            sage: a.save(dir + 'wave.avi', show_path=True)   # optional -- ffmpeg
            Animation saved to file .../wave.avi.
            sage: a.save(dir + 'wave0.sobj')
            sage: a.save(dir + 'wave1.sobj', show_path=True)
            Animation saved to file .../wave1.sobj.
        """
        if filename is None:
            suffix = '.gif'
        else:
            suffix = os.path.splitext(filename)[1]
            if len(suffix) == 0:
                suffix = '.gif'

        if filename is None or suffix == '.gif':
            self.gif(savefile=filename, show_path=show_path,
                     use_ffmpeg=use_ffmpeg)
            return
        elif suffix == '.sobj':
            SageObject.save(self, filename)
            if show_path:
                print "Animation saved to file %s." % filename
            return
        else:
            self.ffmpeg(savefile=filename, show_path=show_path)
            return

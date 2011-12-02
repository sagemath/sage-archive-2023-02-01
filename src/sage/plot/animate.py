"""
Animated plots

EXAMPLES:
We plot a circle shooting up to the right::

    sage: a = animate([circle((i,i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
    ...               xmin=0,ymin=0,xmax=2,ymax=2,figsize=[2,2])
    sage: a.show() # optional -- ImageMagick
"""

############################################################################
#  Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

import os

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
        sage: a.show()          # optional -- ImageMagick
        sage: a[:5].show()      # optional -- ImageMagick

    The ``show`` function takes arguments to specify the
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
        ...       L += line([(i,i^2),(i+step,(i+step)^2)], rgbcolor=(1,0,0), thickness=2)
        ...       v.append(L)
        sage: a = animate(v, xmin=0, ymin=0)
        sage: a.show() # optional -- ImageMagick
        sage: show(L)

    TESTS: This illustrates that ticket #2066 is fixed (setting axes
    ranges when an endpoint is 0)::

        sage: animate([plot(sin, -1,1)], xmin=0, ymin=0)._Animation__kwds['xmin']
        0

    We check that Trac #7981 is fixed::

        sage: a = animate([plot(sin(x + float(k)), (0, 2*pi), ymin=-5, ymax=5)
        ...            for k in srange(0,2*pi,0.3)])
        sage: a.show() # optional -- ImageMagick
    """
    def __init__(self, v, **kwds):
        r"""
        Return an animation of a sequence of plots of objects.

        See documentation of ``animate`` for more details and examples.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.3)],
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: a
            Animation with 21 frames
        """
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

        Test that the bug reported in ticket #12107 has been fixed::

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
        Get a frame from an animation.

        EXAMPLES::

            sage: a = animate([x, x^2, x^3, x^4])
            sage: a[2].show()       # optional -- ImageMagick
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
            sage: a.show() # optional -- ImageMagick
            sage: a[3:7]
            Animation with 4 frames
            sage: a[3:7].show() # optional -- ImageMagick
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
            sage: a.show()   # optional -- ImageMagick
            sage: b = animate([circle((0,i),1,hue=0) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2])
            sage: b.show() # optional
            sage: (a*b).show()    # optional -- ImageMagick
            sage: (a+b).show()    # optional -- ImageMagick
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
            sage: a.show()     # optional -- ImageMagick
            sage: b = animate([circle((0,i),1,hue=0,thickness=20*i) for i in srange(0,2,0.4)],
            ...                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2], axes=False)
            sage: b.show()             # optional -- ImageMagick
            sage: (a*b).show()         # optional -- ImageMagick
            sage: (a+b).show()         # optional -- ImageMagick
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self.__kwds, other.__kwds)

        return Animation(self.__frames + other.__frames, **kwds)

    def png(self, dir=None):
        """
        Return the absolute path to a temp directory that contains the
        rendered PNG's of all the images in this animation.

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
            sage: a.show()        # optional -- ImageMagick

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

    def gif(self, delay=20, savefile=None, iterations=0, show_path=False,
            use_ffmpeg=False):
        r"""
        Returns an animated gif composed from rendering the graphics
        objects in self.

        This function will only work if either (a) the ImageMagick
        software suite is installed, i.e., you have the ``convert``
        command or (b) ``ffmpeg`` is installed.  See
        www.imagemagick.org for more about ImageMagic, and see
        www.ffmpeg.org for more about ``ffmpeg``.  By default, this
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

        If savefile is not specified: in notebook mode, display the
        animation; otherwise, save it to a default file name.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir() + '/'
            sage: a.gif()              # not tested
            sage: a.gif(savefile=dir + 'my_animation.gif', delay=35, iterations=3)  # optional -- ImageMagick
            sage: a.gif(savefile=dir + 'my_animation.gif', show_path=True) # optional -- ImageMagick
            Animation saved to .../my_animation.gif.
            sage: a.gif(savefile=dir + 'my_animation.gif', show_path=True, use_ffmpeg=True) # optional -- ffmpeg
            Animation saved to .../my_animation.gif.

        .. note::

           If neither ffmpeg nor ImageMagick is installed, you will
           get an error message like this::

              Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
              animation to a GIF file or displaying an animation requires one of these
              packages, so please install one of them and try again.

              See www.imagemagick.org and www.ffmpeg.org for more information.

        AUTHORS:

        - William Stein
        """
        from sage.misc.sage_ostools import have_program
        have_convert = have_program('convert')
        have_ffmpeg = self._have_ffmpeg()
        if use_ffmpeg or not have_convert:
            if have_ffmpeg:
                self.ffmpeg(savefile=savefile, show_path=show_path,
                            output_format='gif', delay=delay,
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
                savefile = sage.misc.misc.graphics_filename(ext='gif')
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
Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
animation to a GIF file or displaying an animation requires one of these
packages, so please install one of them and try again.

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
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
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
        if plot.DOCTEST_MODE:
            filename = sage.misc.misc.tmp_filename() + '.gif'
            self.gif(savefile=filename, delay=delay, iterations=iterations)
            return

        if plot.EMBEDDED_MODE:
            self.gif(delay = delay, iterations = iterations)
        else:
            filename = sage.misc.misc.tmp_filename() + '.gif'
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
               ffmpeg_options='', delay=None, iterations=0, verbose=False):
        """
        Returns a movie showing an animation composed from rendering
        the graphics objects in self.

        This function will only work if ffmpeg is installed.  See
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
          suffix, say 'savefile=animation.avi', try to determine the
          format ('avi' in this case) from that file name.  If no file
          is specified or if the suffix cannot be determined, 'mpg' is
          used.

        - ``ffmpeg_options`` - string (default: ''); this string is
          passed directly to ffmpeg.

        - ``delay`` - integer (default: None) delay in hundredths of a
          second between frames; i.e., the framerate is 100/delay.
          This is not supported for mpeg files: for mpegs, the frame
          rate is always 25 fps.

        - ``iterations`` - integer (default: 0); number of iterations
          of animation. If 0, loop forever.  This is only supported
          for animated gif output.

        - ``verbose`` - boolean (default: False); if True, print
          messages produced by the ffmpeg command.

        If savefile is not specified: in notebook mode, display the
        animation; otherwise, save it to a default file name.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir() + '/'
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
        """
        if not self._have_ffmpeg():
            msg = """Error: ffmpeg does not appear to be installed. Saving an animation to
a movie file in any format other than GIF requires this software, so
please install it and try again."""
            raise OSError, msg
        else:
            if not savefile:
                if output_format is None:
                    output_format = 'mpg'
                savefile = sage.misc.misc.graphics_filename(ext=output_format)
            else:
                if output_format is None:
                    suffix = os.path.splitext(savefile)[1]
                    if len(suffix) > 0:
                        suffix = suffix.lstrip('.')
                        output_format = suffix
                    else:
                        output_format = 'mpg'
            if not savefile.endswith('.' + output_format):
                savefile += '.' + output_format
            early_options = ''
            if output_format == 'gif':
                ffmpeg_options += ' -pix_fmt rgb24 -loop_output %s ' % iterations
            if delay is not None and output_format != 'mpeg' and output_format != 'mpg':
                early_options += ' -r %s -g 3 ' % int(100/delay)
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
                if verbose:
                    print "Executing %s " % cmd
                    check_call(cmd, shell=True)
                else:
                    check_call(cmd, shell=True, stderr=PIPE)
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
        determined, default to gif.

        For GIF files, either ffmpeg or the ImageMagick suite must be
        installed.  For other movie formats, ffmpeg must be installed.
        An sobj file can be saved with no extra software installed.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ...                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir() + '/'
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

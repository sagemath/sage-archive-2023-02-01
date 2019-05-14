# -*- coding: utf-8 -*-
r"""
Graphics arrays and insets

"""
from __future__ import print_function, absolute_import
import os
from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject
from sage.misc.temporary_file import tmp_filename
from .graphics import Graphics, ALLOWED_EXTENSIONS, _parse_figsize


class MultiGraphics(WithEqualityById, SageObject):
    r"""
    Base class for objects composed of :class:`~sage.plot.graphics.Graphics`
    objects.

    """
    def __init__(self):
        r"""
        Initialize the attributes common to all MultiGraphics objects.

        TESTS::

            sage: from sage.plot.multigraphics import MultiGraphics
            sage: a = MultiGraphics()
            sage: a._glist
            []

        """
        self._glist = []

    def _repr_(self):
        r"""
        Representation of ``self``.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n, (x,0,1), color=R[n]) for n in range(6)]
            sage: graphics_array(L, 2, 3)
            Graphics Array of size 2 x 3

        """
        return str(self)

    def _rich_repr_(self, display_manager, **kwds):
        r"""
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        .. TODO::

           This method is identical to Graphics._rich_repr_ so it could be
           inherited from a common base class

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: g = graphics_array([Graphics(), Graphics()], 1, 2)
            sage: g._rich_repr_(dm)
            OutputImagePng container

        """
        types = display_manager.types
        prefer_raster = (
            ('.png', types.OutputImagePng),
            ('.jpg', types.OutputImageJpg),
            ('.gif', types.OutputImageGif),
        )
        prefer_vector = (
            ('.svg', types.OutputImageSvg),
            ('.pdf', types.OutputImagePdf),
        )
        graphics = display_manager.preferences.graphics
        if graphics == 'disable':
            return
        elif graphics == 'raster' or graphics is None:
            preferred = prefer_raster + prefer_vector
        elif graphics == 'vector':
            preferred = prefer_vector + prefer_raster
        else:
            raise ValueError('unknown graphics output preference')
        for file_ext, output_container in preferred:
            if output_container in display_manager.supported_output():
                return display_manager.graphics_from_save(
                    self.save, kwds, file_ext, output_container)

    def __getitem__(self, i):
        """
        Return the ``i``th element of the list of graphics composing ``self``.

        EXAMPLES:

        We can access and view individual plots::

            sage: M = [[plot(x^2)], [plot(x^3)]]
            sage: H = graphics_array(M)
            sage: H[1]
            Graphics object consisting of 1 graphics primitive

        Another example::

            sage: L = [plot(sin(k*x), (x,-pi,pi)) + circle((k,k), 1, color='red')
            ....:      for k in range(10)]
            sage: G = graphics_array(L, 5, 2)
            sage: G[3]
            Graphics object consisting of 2 graphics primitives

        """
        return self._glist[i]

    def __setitem__(self, i, g):
        """
        Set the ``i``th element of the list of graphics composing ``self``.

        EXAMPLES::

            sage: M = [[plot(x^2)], [plot(x^3)]]
            sage: H = graphics_array(M)
            sage: H[1] # the plot of x^3
            Graphics object consisting of 1 graphics primitive

        Now we change it::

            sage: H[1] = circle((1,1), 2) + points([(1,2), (3,2), (5,5)], color='purple')
            sage: H[1] # a circle and some purple points
            Graphics object consisting of 2 graphics primitives

        """
        self._glist[i] = g

    def __len__(self):
        """
        Total number of elements graphics object composing ``self``.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n, (x,0,1), color=R[n]) for n in range(6)]
            sage: G = graphics_array(L, 2, 3)
            sage: len(G)
            6

        """
        return len(self._glist)

    def matplotlib(self, figure=None, figsize=None, **kwds):
        r"""
        Create a matplotlib ``Figure`` object from ``self``.

        INPUT:

        - ``figure`` -- (default: ``None``) matplotlib figure (class
          ``matplotlib.figure.Figure``) on which ``self`` is to be displayed;
          if none is provided, the figure will be created from the parameter
          ``figsize``

        - ``figsize`` -- (default: ``None``) width or [width, height] in inches
          of the matplotlib figure in case ``figure`` is ``None``; if none
          is provided matplotlib's default is used

        - ``kwds`` -- options passed to the
          :meth:`~sage.plot.graphics.Graphics.matplotlib` method of
          each graphics object constituting ``self``

        OUTPUT:

        - a ``matplotlib.figure.Figure`` object

        EXAMPLES:

        """
        from matplotlib.figure import Figure
        glist = self._glist
        dims = self._dims
        if dims == 0:  # empty MultiGraphics
            glist = [Graphics()]
            dims = 1
        # If no matplotlib figure is provided, it is created here:
        if figure is None:
            if figsize is not None:
                figsize = _parse_figsize(figsize)
            figure = Figure(figsize=figsize)
        global do_verify
        do_verify = True
        for i, g in zip(range(1, dims + 1), glist):
            # Creation of the matplotlib Axes object "subplot" for g
            subplot = self._add_subplot(figure, i)
            # Setting the options for g.matplotlib:
            options = {}
            options.update(Graphics.SHOW_OPTIONS)  # default options for show()
            options['legend_options'] = Graphics.LEGEND_OPTIONS  # default legend options
            options.update(g._extra_kwds)  # options set in g
            options.update(kwds)
            # We get rid of options that are not relevant for g.matplotlib():
            options.pop('dpi', None)
            options.pop('transparent', None)
            options.pop('fig_tight', None)
            # Adding g to the figure via subplot:
            g.matplotlib(figure=figure, sub=subplot, verify=do_verify,
                         **options)
        return figure

    def save(self, filename, figsize=None, **kwds):
        r"""
        Save ``self``.

        INPUT:

        - ``filename`` -- string. The filename and the image format
          given by the extension, which can be one of the following:

            * ``.eps``,

            * ``.pdf``,

            * ``.png``,

            * ``.ps``,

            * ``.sobj`` (for a Sage object you can load later),

            * ``.svg``,

            * empty extension will be treated as ``.sobj``.

        - ``figsize`` -- (default: ``None``) width or [width, height] in inches
          of the matplotlib figure in case ``figure`` is ``None``; if none
          is provided matplotlib's default is used

        EXAMPLES::

            sage: F = tmp_filename(ext='.png')
            sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in [1..3]]
            sage: G = graphics_array(L)
            sage: G.save(F, dpi=500, axes=False)  # long time (6s on sage.math, 2012)

        TESTS::

            sage: graphics_array([]).save(F)
            sage: graphics_array([[]]).save(F)
        """
        from matplotlib import rcParams
        figure = self.matplotlib(figsize=figsize, **kwds)
        transparent = kwds.get('transparent', Graphics.SHOW_OPTIONS['transparent'])
        fig_tight = kwds.get('fig_tight', Graphics.SHOW_OPTIONS['fig_tight'])
        dpi = kwds.get('dpi', Graphics.SHOW_OPTIONS['dpi'])
        ext = os.path.splitext(filename)[1].lower()
        if ext not in ALLOWED_EXTENSIONS:
            raise ValueError("allowed file extensions for images are '" +
                             "', '".join(ALLOWED_EXTENSIONS) + "'!")
        rc_backup = (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'],
                     rcParams['text.usetex']) # save the rcParams
        # You can output in PNG, PS, EPS, PDF, PGF, or SVG format, depending
        # on the file extension.
        # PGF is handled by a different backend
        if ext == '.pgf':
            from sage.misc.sage_ostools import have_program
            latex_implementations = [i for i in ["xelatex", "pdflatex",
                                                 "lualatex"]
                                     if have_program(i)]
            if not latex_implementations:
                raise ValueError("Matplotlib requires either xelatex, "
                                 "lualatex, or pdflatex.")
            if latex_implementations[0] == "pdflatex":
                # use pdflatex and set font encoding as per
                # matplotlib documentation:
                # https://matplotlib.org/users/pgf.html#pgf-tutorial
                pgf_options= {"pgf.texsystem": "pdflatex",
                              "pgf.preamble": [
                                  r"\usepackage[utf8x]{inputenc}",
                                  r"\usepackage[T1]{fontenc}"
                              ]
                }
            else:
                pgf_options = {
                    "pgf.texsystem": latex_implementations[0],
                }
            rcParams.update(pgf_options)
            from matplotlib.backends.backend_pgf import FigureCanvasPgf
            figure.set_canvas(FigureCanvasPgf(figure))
        # matplotlib looks at the file extension to see what the renderer should be.
        # The default is FigureCanvasAgg for PNG's because this is by far the most
        # common type of files rendered, like in the notebook, for example.
        # if the file extension is not '.png', then matplotlib will handle it.
        else:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            figure.set_canvas(FigureCanvasAgg(figure))
        # tight_layout adjusts the *subplot* parameters so ticks aren't cut off, etc.
        figure.tight_layout()
        opts = dict(dpi=dpi, transparent=transparent)
        if fig_tight is True:
            opts['bbox_inches'] = 'tight'
        figure.savefig(filename, **opts)
        # Restore the rcParams to the original, possibly user-set values
        (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'],
                                       rcParams['text.usetex']) = rc_backup

    def save_image(self, filename=None, *args, **kwds):
        r"""
        Save an image representation of ``self``.  The image type is
        determined by the extension of the filename.  For example,
        this could be ``.png``, ``.jpg``, ``.gif``, ``.pdf``,
        ``.svg``.  Currently this is implemented by calling the
        :meth:`save` method of self, passing along all arguments and
        keywords.

        .. NOTE::

            Not all image types are necessarily implemented for all
            graphics types.  See :meth:`save` for more details.

        EXAMPLES::

            sage: plots = [[plot(m*cos(x + n*pi/4), (x,0, 2*pi)) for n in range(3)] for m in range(1,3)]
            sage: G = graphics_array(plots)
            sage: G.save_image(tmp_filename(ext='.png'))

        """
        self.save(filename, *args, **kwds)

    def _latex_(self, **kwds):
        """
        Return a string plotting ``self`` with PGF.

        INPUT:

        All keyword arguments will be passed to the plotter.

        OUTPUT:

        A string of PGF commands to plot ``self``

        EXAMPLES::

            sage: A = graphics_array([plot(sin), plot(cos)])
            sage: A._latex_()[:40]
            '%% Creator: Matplotlib, PGF backend\n%%\n%'

        """
        tmpfilename = tmp_filename(ext='.pgf')
        self.save(filename=tmpfilename, **kwds)
        with open(tmpfilename, "r") as tmpfile:
            latex_list = tmpfile.readlines()
        return ''.join(latex_list)

    def show(self, **kwds):
        r"""
        Show ``self`` immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        OPTIONAL INPUT:

        -  ``dpi`` - dots per inch

        - ``figsize`` -- (default: ``None``) width or [width, height] in inches
          of the matplotlib figure in case ``figure`` is ``None``; if none
          is provided matplotlib's default is used

        -  ``axes`` - (default: True)

        -  ``fontsize`` - positive integer

        -  ``frame`` - (default: False) draw a frame around the
           image

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        EXAMPLES:

        This draws a graphics array with four trig plots and no
        axes in any of the plots::

            sage: G = graphics_array([[plot(sin), plot(cos)], [plot(tan), plot(sec)]])
            sage: G.show(axes=False)

        .. PLOT::

            G = graphics_array([[plot(sin), plot(cos)], [plot(tan), plot(sec)]])
            sphinx_plot(G, axes=False)

        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, **kwds)

    def plot(self):
        """
        Draw a 2D plot of this graphics object, which just returns this
        object since this is already a 2D graphics object.

        EXAMPLES::

            sage: g1 = plot(cos(20*x)*exp(-2*x), 0, 1)
            sage: g2 = plot(2*exp(-30*x) - exp(-3*x), 0, 1)
            sage: S = graphics_array([g1, g2], 2, 1)
            sage: S.plot() is S
            True
        """
        return self




# ****************************************************************************


class GraphicsArray(MultiGraphics):
    """
    GraphicsArray takes a (`m` x `n`) list of lists of
    graphics objects and plots them all on one canvas.

    """
    def __init__(self, array):
        """
        Constructor for ``GraphicsArray`` class.  Normally used only
        via :func:`graphics_array` function.

        INPUT: a list or list of lists/tuples, all of which are graphics objects

        EXAMPLES::

            sage: L = [plot(sin(k*x),(x,-pi,pi)) for k in range(10)]
            sage: G = graphics_array(L)
            sage: G.ncols()
            10
            sage: M = [[plot(x^2)],[plot(x^3)]]
            sage: H = graphics_array(M)
            sage: str(H[1])
            'Graphics object consisting of 1 graphics primitive'

        TESTS::

            sage: L = [[plot(sin),plot(cos)],[plot(tan)]]
            sage: graphics_array(L)
            Traceback (most recent call last):
            ...
            TypeError: array (=[[Graphics object consisting of 1 graphics primitive, Graphics object consisting of 1 graphics primitive], [Graphics object consisting of 1 graphics primitive]]) must be a list of lists of Graphics objects
            sage: G = plot(x,(x,0,1))
            sage: graphics_array(G)
            Traceback (most recent call last):
            ...
            TypeError: array (=Graphics object consisting of 1 graphics primitive) must be a list of lists of Graphics objects
            sage: G = [[plot(x,(x,0,1)),x]]
            sage: graphics_array(G)
            Traceback (most recent call last):
            ...
            TypeError: every element of array must be a Graphics object

            sage: hash(graphics_array([])) # random
            42
        """
        MultiGraphics.__init__(self)
        if not isinstance(array, (list, tuple)):
            raise TypeError("array (=%s) must be a list of lists of Graphics objects"%(array))
        array = list(array)
        self._rows = len(array)
        if self._rows > 0:
            if not isinstance(array[0], (list, tuple)):
                array = [array]
                self._rows = 1
            self._cols = len(array[0])
        else:
            self._cols = 0
        self._dims = self._rows*self._cols
        for row in array: # basically flatten the list
            if not isinstance(row, (list, tuple)) or len(row) != self._cols:
                raise TypeError("array (=%s) must be a list of lists of Graphics objects"%(array))
            for g in row:
                if not isinstance(g, Graphics):
                    raise TypeError("every element of array must be a Graphics object")
                self._glist.append(g)

    def __str__(self):
        """
        String representation of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: G.__str__()
            'Graphics Array of size 2 x 3'
            sage: str(G)
            'Graphics Array of size 2 x 3'
        """
        return "Graphics Array of size %s x %s"%(self._rows, self._cols)

    def _add_subplot(self, figure, index, **options):
        r"""
        Add a subplot to a given matplotlib ``Figure``, the position of
        which is governed by a given element of ``self``.

        This method encapsulates the matplotlib method ``Figure.add_subplot``
        and is intended to be called by :meth:`MultiGraphics.save`

        INPUT:

        - ``figure`` -- a matplotlib Figure
        - ``index `` -- integer specifiying the element of ``self``
        - ``options`` -- extra options to be passed to ``Figure.add_subplot``

        OUTPUT:

        - a matplotlib ``Axes`` object

        """
        if self._dims == 0:
            rows = 1
            cols = 1
        else:
            rows = self._rows
            cols = self._cols
        return figure.add_subplot(rows, cols, index, **options)


    def nrows(self):
        """
        Number of rows of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: G.nrows()
            2
            sage: graphics_array(L).nrows()
            1
        """
        return self._rows

    def ncols(self):
        """
        Number of columns of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n,(x,0,1),color=R[n]) for n in range(6)]
            sage: G = graphics_array(L,2,3)
            sage: G.ncols()
            3
            sage: graphics_array(L).ncols()
            6
        """
        return self._cols

    def append(self, g):
        """
        Appends a graphic to the array.  Currently
        not implemented.

        TESTS::

            sage: from sage.plot.multigraphics import GraphicsArray
            sage: G = GraphicsArray([plot(sin),plot(cos)])
            sage: G.append(plot(tan))
            Traceback (most recent call last):
            ...
            NotImplementedError: Appending to a graphics array is not yet implemented
        """
        # Not clear if there is a way to do this
        raise NotImplementedError('Appending to a graphics array is not yet implemented')


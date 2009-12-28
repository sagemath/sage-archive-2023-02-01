Make SageTeX known to TeX
===========================

Sage is largely self-contained, but some parts do need some intervention
to work properly. SageTeX is one such part.

The SageTeX package allows one to embed computations and plots from Sage
into a LaTeX document. SageTeX is installed in Sage by default, but to
use SageTeX with your LaTeX documents, you need to make your TeX
installation aware of it before it will work.

The key to this is that TeX needs to be able to find ``sagetex.sty``,
which can be found in
``SAGE_ROOT/local/share/texmf/tex/generic/sagetex/``, where
``SAGE_ROOT`` is the directory where you built or installed Sage. If
TeX can find ``sagetex.sty``, then SageTeX will work. There are several
ways to accomplish this.

- The first and simplest way is simply to copy ``sagetex.sty`` into the
  same directory as your LaTeX document. Since the current directory is
  always searched when typesetting a document, this will always work.

  There are a couple small problems with this, however: the first is
  that you will end up with many unnecessary copies of ``sagetex.sty``
  scattered around your computer. The second and more serious problem is
  that if you upgrade Sage and get a new version of SageTeX, the Python
  code and LaTeX code for SageTeX may no longer match, causing errors.

- The second way is to use the ``TEXINPUTS`` environment variable. If
  you are using the bash shell, you can do

  ::

      export TEXINPUTS="SAGE_ROOT/local/share/texmf//:"

  where ``SAGE_ROOT`` is the location of your Sage installation. Note
  that the double slash and colon at the end of that line are important.
  Thereafter, TeX and friends will find the SageTeX style file. If you
  want to make this change permanent, you can add the above line to your
  ``.bashrc`` file. If you are using a different shell, you may have to
  modify the above command to make the environment variable known; see
  your shell's documentation for how to do that.

  One flaw with this method is that if you use applications like
  TeXShop, Kile, or Emacs/AucTeX, they will not necessarily pick up the
  environment variable, since when they run LaTeX, they may do so
  outside your usual shell environment.

  If you ever move your Sage installation, or install a new version into
  a new directory, you'll need to update the above command to reflect
  the new value of ``SAGE_ROOT``.

- The third (and best) way to make TeX aware of ``sagetex.sty`` is to
  copy that file into a convenient place in your home directory. In most
  TeX distributions, the ``texmf`` directory in your home directory is
  automatically searched for packages. To find out exactly what this
  directory is, do the following on the command line::

      kpsewhich -var-value=TEXMFHOME

  which will print out a directory, such as ``/home/drake/texmf`` or
  ``/Users/drake/Library/texmf``. Copy the ``tex/`` directory from
  ``SAGE_ROOT/local/share/texmf/`` into your home ``texmf`` directory
  with a command like

  ::

      cp -R SAGE_ROOT/local/share/texmf/tex TEXMFHOME

  where ``SAGE_ROOT`` is, as usual, replaced with the location of your
  Sage installation and ``TEXMFHOME`` is the result of the
  ``kpsewhich`` command above.

  If you upgrade Sage and discover that SageTeX no longer works, you can
  simply repeat these steps and the Sage and TeX parts of SageTeX will
  again be synchronized.

.. _sagetex_installation_multiuser:

- For installation on a multiuser system, you just modify the above
  instructions appropriately to copy ``sagetex.sty`` into a systemwide
  TeX directory. Instead of the directory ``TEXMFHOME``, probably the
  best choice is to use the result of

  ::

      kpsewhich -var-value=TEXMFLOCAL

  which will likely produce something like ``/usr/local/share/texmf``.
  Copy the ``tex`` directory as above into the ``TEXMFLOCAL``
  directory. Now you need to update TeX's database of packages, which
  you can do simply by running

  ::

      texhash TEXMFLOCAL

  as root, replacing ``TEXMFLOCAL`` appropriately. Now all users of your
  system will have access to the LaTeX package, and if they can also run
  Sage, they will be able to use SageTeX.

.. warning::

  it's very important that the file ``sagetex.sty`` that LaTeX uses when
  typesetting your document match the version of SageTeX that Sage is
  using. If you upgrade your Sage installation, you really should delete
  all the old versions of ``sagetex.sty`` floating around.

  Because of this problem, we recommend copying the SageTeX files into
  your home directory's texmf directory (the third method above). Then
  there is only one thing you need to do (copy a directory) when you
  upgrade Sage to insure that SageTeX will work properly.

SageTeX documentation
---------------------

While not strictly part of installation, it bears mentioning here that
the documentation for SageTeX is maintained in
``SAGE_ROOT/local/share/texmf/tex/generic/sagetex/sagetexpackage.pdf``.
There is also an example file in the same directory -- see
``example.tex`` and ``example.pdf``, the pre-built result of typesetting
that file with LaTeX and Sage. You can also get those files from
:abbr:`CTAN (Comprehensive TeX Archive Network)`:
http://www.ctan.org/tex-archive/macros/latex/contrib/sagetex/.

However, be advised that SageTeX is "officially" distributed as part of
Sage, and the version available on CTAN may not be the most current
version.

SageTeX and TeXLive
-------------------

One potentially confusing issue is that the popular TeX distribution
`TeXLive 2009 <http://www.tug.org/texlive/>`_ includes SageTeX. This may
seem nice, but with SageTeX, it's important that the Sage bits and LaTeX
bits be synchronized -- which is a problem in this case, since both Sage
and SageTeX are updated frequently, and TeXLive is not. In fact, at the
time of this writing (January 2010), many Linux distributions are still
shipping TeXLive 2007.

Because of this, it is *strongly recommended* that you always install
the LaTeX part of SageTeX from Sage, as described above. The
instructions above will insure that both halves of SageTeX are
compatible and will work properly. Using TeXLive to provide the LaTeX
side of SageTeX is not supported.


The Documentation
=================

You do not need to install the documentation separately. This
includes the tutorial, the manual, the programming guide, and this
installation guide. The tutorial is an overview about the
philosophy of Sage, and some examples of how to use Sage. The
reference manual explains what is available in Sage and gives
examples of how to use it.

For Sage to build the HTML, DVI, and PDF versions of the
documentation, it requires a working LaTeX system with ``latex2html``.
The source to create the Sage documentation is included with the
source tarball. This documentation is in the subdirectory
``devel/doc-*``. To build it, change into this directory and type
``./rebuild``. All the HTML documentation is in the ``doc``
subdirectory (which is actually a link to ``devel/doc-*/html``).

    Some systems exhibit a LaTeX bug where the ``url`` command gets defined
    twice, and this prevents building the documentation which use LINKS
    links. A hack to get around this, at least on an Ubuntu 5.10
    system, is to comment out line 103 of
    ``/usr/share/texmf/tex/latex/misc/url.sty``. If you can figure out a
    much better fix, please post it on the ``sage-devel`` Google group.

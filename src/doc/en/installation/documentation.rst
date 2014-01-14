
The Documentation
=================

You do not need to install the documentation separately: it is
included with Sage. The Sage standard documentation includes a guided
tour of Sage, a tutorial, a reference manual, a developer's guide, an
installation guide, and other documents. The tutorial is a good
starting point to learn how to use Sage. The reference manual
describes what is available in Sage along with examples on how to use
specific commands.

When you build Sage from source, by default the HTML version of the
standard documentation is built as well. But note that in this way the
HTML version does not have links to the PDF version. To build the HTML
or PDF version of the documentation yourself, use the general command ::

    sage -docbuild {document} {format}

For example, the command ::

    sage -docbuild reference html

builds the HTML version of the reference manual.

You can choose to have the built HTML version of the documentation
link to the PDF version. To do so, you need to build both the HTML and
PDF versions. To have the HTML version link to the PDF version, do ::

    sage -docbuild all html
    sage -docbuild all pdf

Type ``sage -docbuild -H`` to see a list of available options for
building the documentation or any part of it. See the file
``SAGE_ROOT/Makefile`` for further information on how the
documentation is built by default.

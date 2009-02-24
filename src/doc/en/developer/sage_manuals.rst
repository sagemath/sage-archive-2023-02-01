.. _chapter-sage_manuals:

=================
The Sage Manuals
=================

This chapter describes how to modify the Sage manuals. Sage's
manuals are written in ReST, and to edit them, you just need to
edit the appropriate file.

The Sage manuals and the corresponding files to edit:

-  The Sage tutorial: ``SAGE_ROOT/devel/sage/doc/en/tutorial``

-  The Sage developer's guide:
   ``SAGE_ROOT/devel/sage/doc/en/developer``

-  Constructions in Sage:
   ``SAGE_ROOT/devel/sage/doc/en/constructions``

-  The Sage installation guide:
   ``SAGE_ROOT/devel/sage/doc/en/installation``

-  The Sage reference manual: some of this is contained in the file
   ``SAGE_ROOT/devel/sage/doc/en/reference``, but most of it is
   automatically generated from the Sage source code.

-  Additional, more specialied  manuals can be found under
   ``SAGE_ROOT/devel/sage/doc/en`` as well.

.. note::

   You can edit manuals that have been translated into another language
   by replacing the ``en/`` above with the appropriate two letter
   language code.  For example, the French tutorial is located in
   ``SAGE_ROOT/devel/sage/doc/fr/tutorial``

Editing the Manuals
-------------------

If, for example, you want to change the Sage tutorial, then you should
start by modifying the files in
``SAGE_ROOT/devel/sage/doc/en/tutorial/``. Then to build a PDF file
with your changes, type::

    sage -docbuild tutorial pdf

You'll get a file ``tutorial.pdf`` in
``SAGE_ROOT/devel/sage/doc/output/pdf/en/tutorial`` which you should
inspect.  You can build the HTML version
of the tutorial by typing::

    sage -docbuild tutorial html

Once you've done this, you can access the new HTML version from the
notebook interface to Sage by clicking the ``Help`` link.  For more
detailed information about building the documentation, see
:ref:`building`.

You should also run

::

    sage -tp 1 SAGE_ROOT/devel/sage/doc/en/tutorial/

to test all of the examples in the tutorial - see
:ref:`chapter-testing`.

Finally, you might want to share your changes with the Sage
community. To do this, use Mercurial (see :ref:`chapter-mercurial`) to
produce patch files, and submit them to the Sage trac server.

.. _building:

Building the Manuals
--------------------

All of the Sage manuals are built using the ``sage -docbuild``
script.  The content of the ``sage -docbuild`` script is defined in
``SAGE_ROOT/devel/sage/doc/common/builder.py``.  It is a thin wrapper
around the ``sphinx-build`` script which does all of the real work.
It is designed to be a replacement for the default Makefiles generated
by the ``sphinx-quickstart`` script.

Document Names
~~~~~~~~~~~~~~
A name of a document is of the form::

    lang/name

where ``lang`` is a two-letter language code, and name is the
descriptive name of the document.  If the language isn't specified,
then it defaults to English (``en``).  The following two commands do
the exact same thing::

    sage -docbuild tutorial html
    sage -docbuild en/tutorial html

To specify the French version of the tutorial, you would simply run::

    sage -docbuild fr/tutorial html

Output Formats
~~~~~~~~~~~~~~

The Sage documentation build system currently supports all of the
output formats that Sphinx does.

For more detailed information, see the documentation on builders at
http://sphinx.pocoo.org/builders.html .
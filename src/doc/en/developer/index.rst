
======================================
Welcome to the Sage Developer's Guide!
======================================

Everybody who uses Sage is encouraged to contribute something back to
Sage at some point. You could:

* Add examples to the documentation
* Find bugs or typos
* Fix a bug
* Implement a new function
* Create a new class, create a fast new C library, etc.

This document describes how to write programs using Sage, how to modify
and extend the core Sage libraries, and how to modify Sage's
documentation. We also discuss how to share your new and modified code
with other Sage users around the globe.

All development takes place on or via `the Sage Trac server
<http://trac.sagemath.org>`_,
including bug reports, fixes, new functionality, and discussions about
approaches to particular tickets.  Once you start writing code for Sage,
you will want to carefully read
:ref:`conventions and guidelines <section-writing-code-for-sage>` we use.

Depending on your previous knowledge, there are several places you can
start learning about the source code revision control process.

-  First, although it is possible to try out bugfixes and explore the
   code without having a developer account, it is best to :ref:`acquire a
   Trac account <section-trac-account>` first, then :ref:`configure git
   <section-git-setup-name>` for use with Trac.
-  An overview of the Sage development process, assuming you have ``git``
   installed and know the basics of how to use it, is in the :ref:`concise
   development guide <chapter-walkthrough>`.
-  More advanced :ref:`tricks and tips <section-git-tricks-and-tips>` for
   ``git`` are linked below.
-  For those unfamiliar with revision control, please start by reading
   about :ref:`collaborative development with Git-Trac <chapter-git_trac>`,
   which provides some easier interface with git and Trac, both for newbies
   and power users.
-  Alternately, one can do certain amounts of Sage development without
   having git installed, by using Sage's own internal installation of git
   and the :ref:`Sage dev scripts <chapter-devscript>`.  This is mainly
   intended as a bridge to full use of git once one becomes more comfortable
   with the system.
-  Finally, if you've never worked on software before, don't forget you
   will need the `prerequisites to compile
   <http://www.sagemath.org/doc/installation/source.html#prerequisites>`_
   in order to make your changes in the source code work.

No matter where you start, good luck and welcome to Sage development!


Walk-Through and First Steps
============================

.. toctree::
   :maxdepth: 3

   walk_through
   git_trac
   dev_script

Git and Trac Reference
======================

.. toctree::
   :maxdepth: 3

   git_setup
   trac


.. _section-git-tricks-and-tips:

Git Tricks & Tips
=================

.. toctree::
   :maxdepth: 3

   manual_git
   git_background
   advanced_git
   workflows

.. _section-writing-code-for-sage:

Writing Code for Sage
=====================

.. toctree::
   :maxdepth: 3

   coding_basics
   coding_in_python
   coding_in_cython
   coding_in_other
   packaging
   packaging_old_spkgs
   doctesting
   sage_manuals


Sage Notebook Developer Guide
=============================

.. toctree::
   :maxdepth: 3

   sagenb/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License <http://creativecommons.org/licenses/by-sa/3.0/>`_.




======================================
Welcome to the Sage Developer's Guide!
======================================

Everybody who uses Sage is encouraged to contribute something back to
Sage at some point. You could:

* Add examples to the documentation
* Find bugs or typos
* Fix a bug
* Implement a new function
* Contribute a useful tutorial for a mathematical topic
* Translate an existing document to a new language
* Create a new class, create a fast new C library, etc.

This document describes how to write programs using Sage, how to modify
and extend the core Sage libraries, and how to modify Sage's
documentation. We also discuss how to share your new and modified code
with other Sage users around the globe.

It is not necessary to memorize this entire guide to begin working on
Sage, but careful reading of different sections now will be well worth
the effort in seamless contributions later.  There are four main things
to be aware of.

- All development takes place on or via `the Sage Trac server
  <http://trac.sagemath.org>`_, including bug reports, fixes,
  new functionality, and discussions about approaches to particular tickets.
  If you don't have an account on it, read about how to :ref:`acquire a
  Trac account <section-trac-account>`.  This is recommended even if you
  only want to report bugs or request new functionality, not necessarily
  to help make changes to Sage.

- Next, if you've never worked on software before, you will want to read
  about the `prerequisites to compile
  <http://www.sagemath.org/doc/installation/source.html#prerequisites>`_
  from the installation guide.  This will allow you to
  make your changes in the source code work.  Pay close attention
  to any system-specific requirements.

- Once you start writing code for Sage, you will want to carefully read the
  :ref:`conventions and guidelines <section-writing-code-for-sage>` we use.
  (Looking at newer files and functionality within Sage is another way to
  get a sense for the the general style, but refer here for details.)

  - There is an entire section on how to modify or add to the various
    :ref:`manuals and tutorials <chapter-sage_manuals>`,
    including localizing to other languages.

- Finally, in order to share changes with the Sage community, you will
  need to learn some basics of the source code revision control process.
  There are several places to start, depending upon your previous knowledge.

  - Don't forget to :ref:`acquire a Trac account <section-trac-account>`.
  - First, you will need to :ref:`install the 'git' revision control
    software <section-git-install>` if you don't already have it.
  - Then you will need to go through a short process to :ref:`configure git
    <section-git-setup-name>` for use with Trac.
  - Assuming you have ``git`` installed and know the basics of how to use it,
    the next step is the overview of the Sage development flow in the
    :ref:`concise development guide <chapter-walkthrough>`.

    - (More advanced :ref:`tricks and tips <section-git-tricks-and-tips>` for
      ``git`` are linked below.)

  - For those unfamiliar with ``git`` or revision control, please start by
    reading about :ref:`collaborative development with Git-Trac <chapter-git_trac>`,
    which provides some easier interface with git and Trac, both for newbies
    and power users.
  - Alternately, one can do certain amounts of Sage development without
    having ``git`` installed, by using Sage's own internal installation of ``git``
    and the :ref:`Sage dev scripts <chapter-devscript>`.  This is mainly
    intended as a bridge to full use of git once one becomes more comfortable
    with the system.

No matter where you start, good luck and welcome to Sage development!


Git and Trac
============

First Steps with Git
--------------------

Sage uses git for version control.

.. toctree::
   :maxdepth: 3

   git_setup
   walk_through

Sage Trac and tickets
---------------------

All changes to Sage source code require a ticket on the
`Sage trac server <http://trac.sagemath.org>`_.

.. toctree::
   :maxdepth: 3

   trac

Git and Trac integration
------------------------

Putting your local changes on a Trac ticket.

.. toctree::
   :maxdepth: 3

   git_trac
   dev_script

.. _section-git-tricks-and-tips:

Git Tricks & Tips
-----------------

When ``git trac`` is not enough.

.. toctree::
   :maxdepth: 3

   manual_git
   git_background
   advanced_git
   workflows


.. _section-writing-code-for-sage:

Writing Code for Sage
=====================

Basics of Writing and Testing Sage Code
---------------------------------------

.. toctree::
   :maxdepth: 3

   coding_basics
   doctesting

Contributing to Manuals and Tutorials
-------------------------------------

.. toctree::
   :maxdepth: 3

   sage_manuals

Sage Coding Details
-------------------

.. toctree::
   :maxdepth: 3

   coding_in_python
   coding_in_cython
   coding_in_other

Packaging Third-Party Code
--------------------------

.. toctree::
   :maxdepth: 3

   packaging
   packaging_old_spkgs

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

.. _using-git:

Working with *Sage Notebook* source code
================================================

Getting Started
---------------

Development of the Sage notebook currently occurs on Github using the
Git revision control system.  However, since Sage ships with
Mercurial, a Mercurial repository is provided in the spkg which
mirrors the Git repository.  You can use Mercurial to investigate the
history and to create small patches.

Howver, in order to do development in the Sage notebook, including
submitting change requests to upstream, we highly recommend
installing Git (see below).

To update to the latest development source, run the commands below,
where ``SAGE_ROOT`` is the root directory of the Sage installation,
and where ``hackdir`` is a directory you create for working on code
changes (it need not have the name or location given below).

.. warning:: This will create a new sagenb repository ignoring any
   changes you have made to the files.

::

    mkdir ~/hackdir
    cd ~/hackdir
    git clone git://github.com/sagemath/sagenb.git sagenb-git
    cd SAGE_ROOT/devel
    rm sagenb
    ln -s ~/hackdir/sagenb sagenb
    cd sagenb
    ../../sage setup.py develop

What this has done is to create a new directory, move to that
directory, and create a clone of the most up-to-date version of
the upstream notebook sources there.  Then we remove a symbolic
link ``sagenb`` in the Sage folder and replace it with a link
to your clone of upstream, finally making sure that the notebook
has the correct dependencies.

An advantage of having the separate directory for sagenb is that
you would later be able to keep it and do development work in it
even when you upgrade Sage, or even if you accidentally destroy your
Sage installation somehow.

The rest of these instructions is some very generic documentation,
slightly adapted to help develop the notebook using Git and Github.

The most important section involves how to update your new sagenb
source repository and create a "fork" of the master copy, so that you
will be able to request your changes to be merged in the Sage notebook,
called a "pull request"; see :ref:`git-development`.  There are also
a number of helpful guides to Git and Github at :ref:`git-resources`.


.. toctree::
   :maxdepth: 2

   git_intro
   git_install
   following_latest
   patching
   git_development
   git_resources


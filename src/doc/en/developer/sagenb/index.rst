.. _using-git:

Working with *Sage Notebook* source code
================================================

Getting Started
---------------

Development of the Sage notebook currently occurs on Github using the
Git revision control system.  However, since Sage ships with
Mercurial, a Mercurial repository is provided in the spkg which
mirrors the Git repository.  You can use Mercurial to investigate the
history and to create and submit small patches.

In order to do serious development in the Sage notebook, we highly
recommend installing Git (see below).

To update to the latest development source, run the commands below,
where ``SAGE_ROOT`` is the root directory of the Sage installation.

.. warning:: This will create a new sagenb repository ignoring any
   changes you have made to the files.

::

    cd SAGE_ROOT/devel
    git clone git://github.com/sagemath/sagenb.git sagenb-git
    rm sagenb
    ln -s sagenb-git sagenb
    cd sagenb
    ../../sage setup.py develop

Now you can use the following instructions to develop the Sage
notebook using Git and Github.

.. toctree::
   :maxdepth: 2

   git_intro
   git_install
   following_latest
   patching
   git_development
   git_resources


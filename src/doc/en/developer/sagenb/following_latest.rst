.. _following-latest:

=============================
 Following the latest source
=============================

These are the instructions if you just want to follow the latest
*Sage Notebook* source, but you don't need to do any development for now.

The steps are:

* :ref:`install-git`
* get local copy of the `Sage Notebook github`_ git repository
* update local copy from time to time

Get the local copy of the code
==============================

From the command line::

   git clone git://github.com/sagemath/sagenb.git

You now have a copy of the code tree in the new ``sagenb`` directory.

Updating the code
=================

From time to time you may want to pull down the latest code.  Do this with::

   cd sagenb
   git pull

The tree in ``sagenb`` will now have the latest changes from the initial
repository.

.. include:: links.inc

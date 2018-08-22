.. _installation-guide:

Welcome to the Sage Installation Guide
======================================

You can install Sage either from the package manager of your distribution, a
pre-built binary tarball or from its sources.

Installing Sage from your distribution package manager is the preferred and
faster solution (dependencies will be automatically taken care of and Sage will
be using your system Python). It is the case at least for the following
GNU/Linux distributions: Debian version >= 9, Ubuntu version >= 18.04 and
archlinux. In this case, install the packages ``sagemath`` and
``sagemath-jupyter``. For GNU/Linux Gentoo, you might want to give a try to
`sage-on-gentoo <https://github.com/cschwan/sage-on-gentoo>`_.

If your operating system does not provide SageMath, then using a pre-built
binary is the fastest method. See the section
:ref:`sec-installation-from-binaries`.

By compiling Sage from its sources you run a slightly more up-to-date
version. You can also modify it and contribute back to the project. Compiling
Sage should be simpler than you're used to with most software, since much
testing is done on a wide range of computers. Though, it might take up to
4 hours on a recent computer. For that option, go to the section
:ref:`sec-installation-from-sources`.

Note that there are other alternatives to use Sage that avoids installing it:

- the `Sage Debian Live USB key <https://sagedebianlive.metelu.net/>`_: a full
  featured USB key that contains a whole Linux distribution including SageMath.
  This might be an option if you fail installing SageMath on your operating
  system.

- `CoCalc <https://cocalc.com/>`_: an online service that provides sagemath and
  many other tools.

- `Sage Cell Server <https://sagecell.sagemath.org/>`_: an online service for
  elementary Sage computations.

The rest of this document describes how to install Sage from pre-built
binaries and from sources.

.. toctree::
   :maxdepth: 2

   binary
   source
   launching
   standard_packages
   troubles

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License`__.

__ http://creativecommons.org/licenses/by-sa/3.0/


Pre-built Binary Install
========================

Linux and OS X
--------------

Installation from a pre-built binary tarball should in the long run
be the easiest and fastest way to install Sage. This is not
necessarily the case right now. Note that Sage is itself a
programming environment, so building it from source guarantees you
maximum flexibility in the long run. Nonetheless, we provide
pre-built binaries.

Assumptions: You have a computer with at least 550 megabytes free
disk space and the operating system is Linux (32-bit or 64-bit) or
OS X.

Highly Recommended: It is highly recommended that you have LaTeX
installed.

Download the latest tarball from
http://www.sagemath.org/download.html . For example, it might be
called ``sage-x.y.z-x86_64-Linux.tgz``. Unpack it on your computer
in a directory which you have permissions:

::

        tar zxvf sage-x.y.z-x86_64-Linux.tgz

Change into the directory just created, e.g.,
``sage-x.y.z-x86_64-Linux`` and type ``./sage`` to run Sage. You can
move the directory ``sage-x.y.z-x86_64-Linux`` anywhere, and still
run ``sage`` from it. You can also copy ``sage`` and put it anywhere,
e.g., ``/usr/local/bin/``, but you'll have likely have to edit the
``ROOT="....."`` line at the top.

We currently distribute .dmg files for OSX. But we would like to
make Sage more of a native application. Work for that is ongoing,
but help is always welcome.


Microsoft Windows
-----------------

The best way to install Sage on Windows is to get the free VMware
player and use the VMware Sage appliance, which is available at
http://www.sagemath.org/bin/microsoft_windows/ . Be sure to read README.txt
in that directory.

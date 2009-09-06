
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

Assumptions: You have a computer with at least 2 GB of free
disk space and the operating system is Linux (32-bit or 64-bit) or
OS X (10.5.x).

Highly Recommended: It is highly recommended that you have LaTeX
installed.

Download the latest binary tarball from
http://www.sagemath.org/download.html . For example, it might be
called ``sage-x.y.z-x86_64-Linux.tgz``. Unpack it on your computer
in a directory which you have permission to read and write:

::

        tar zxvf sage-x.y.z-x86_64-Linux.tgz

Change into the directory just created, e.g.,
``sage-x.y.z-x86_64-Linux`` and type ``./sage`` to run Sage. The first
time you run that command, the necessary environment variables are set
for your system. If this command doesn't result in any errors, you
should then quit Sage using ``exit`` and then start up Sage again
using the command ``./sage -br main``. This would regenerate the
necessary files for your local Sage installation.

You can move the directory ``sage-x.y.z-x86_64-Linux`` anywhere and still
run ``sage`` from it. You can also copy ``sage`` and put it anywhere,
e.g., ``/usr/local/bin/``, but you would have to edit the
``ROOT="....."`` line at the top. If you have moved the directory
``sage-x.y.z-x86_64-Linux`` to somewhere else in your system, it is
recommended that you use a terminal program to cd to the new Sage top
level directory and type ``./sage`` to reset all necessary Sage
environment variables. Then exit Sage with ``exit`` and load Sage
again with ``./sage -br main`` to regenerate necessary files for your
local Sage installation.

We currently distribute ``.dmg`` files for OS X 10.4.x and 10.5.x. But
we would like to make Sage more of a native application. Work for that
is ongoing, but help is always welcome.


Microsoft Windows
-----------------

The best way to install Sage on Windows is to get the free VMware
player and use the VMware Sage appliance, which is available at
http://www.sagemath.org/download-windows.html . Be sure to read the
file README.txt at
http://www.sagemath.org/bin/microsoft_windows/README.txt .

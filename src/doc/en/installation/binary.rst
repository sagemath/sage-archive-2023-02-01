
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

You can move the resulting directory ``sage-x.y.z-x86_64-Linux``
anywhere and still run ``sage`` from it. You can also copy the file
``sage`` from that directory and put it anywhere, e.g.,
``/usr/local/bin/``, but then you have to edit the
``SAGE_ROOT="....."`` line at the top of that file, replacing the dots
with the path to the Sage directory ``sage-x.y.z-x86_64-Linux``.  As
long as ``/usr/local/bin`` is in your ``$PATH``, you can then type
``sage`` from the command line to run Sage.  Another approach is to
create a symbolic link from ``sage-x.y.z-x86_64-Linux`` to, say,
``/usr/local/share/sage``::

    ln -s /.../path_to/.../sage-x.y.z-x86_64-Linux /usr/local/share/sage

Then put ``/usr/local/share/sage`` in your ``$PATH``.  If you do this,
you can type ``sage`` from the command line to run Sage.  Also, if you
install a different version of Sage, you just have to delete the old
link and create one from the new directory to
``/usr/local/share/sage``.

Finally, you can also combine these two approaches, copying ``sage``
to ``/usr/local/bin/``, creating a link from
``sage-x.y.z-x86_64-Linux`` to ``/usr/local/share/sage``, and editing
the file ``/usr/local/bin/sage``: change the line ::

  SAGE_ROOT="....."

to ::

  SAGE_ROOT="/usr/local/share/sage"

When you want to install a new version of Sage, just delete the old
link and create a new one; you shouldn't have to replace or modify the
file ``/usr/local/bin/sage``.

The first time you run Sage, and any time you move the Sage directory
or create a link as above, you may see a message saying

::

   The Sage install tree may have moved.
   Regenerating Python.pyo and .pyc files that hardcode the install PATH
   (please wait at most a few minutes)...
   Do not interrupt this.

We currently distribute ``.dmg`` files for OS X 10.4.x and 10.5.x. But
we would like to make Sage more of a native application. Work for that
is ongoing, but help is always welcome.


Microsoft Windows
-----------------

The best way to install Sage on Windows is to install
`VirtualBox for Windows <http://www.virtualbox.org/wiki/Downloads>`_
and then download and install the VirtualBox distribution of Sage. See
`this URL <http://www.sagemath.org/download-windows.html>`_ for
further instructions on installing Sage on Windows. Be sure to read the
file `README.txt <http://www.sagemath.org/mirror/win/README.txt>`_.

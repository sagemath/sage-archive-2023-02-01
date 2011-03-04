Sage .spkg files
================

The directory SAGE_ROOT/spkg/standard contains spkg's. In a source
install, these are all Sage spkg files (actually .tar or .tar.bz2
files), which are the source code that defines Sage. In a binary
install some of these may be small placeholder files to save space.

Sage packages are distributed as .spkg files, which are .tar.bz2 files
(or tar files) but have the extension .spkg to discourage
confusion. Although Sage packages are packed using tar and/or bzip2,
please note that .spkg files contain control information (installation
scripts and metadata) that are necessary for building and installing
them. For source distributions, when you compile Sage the file
SAGE_ROOT/makefile takes care of the unpacking, compilation, and
installation of Sage packages for you. For more information on the
structure of .spkg files, please refer to the Sage Developer's Guide
in your local installation of Sage at

SAGE_ROOT/sage/doc/output/html/en/developer/index.html

If you cannot locate that file in your local installation of Sage, you
might want to consider (re)building the standard Sage documentation
using this command:

SAGE_ROOT/sage -docbuild all html

or visit the URL

http://www.sagemath.org/doc/developer/

Additional Sage packages can be found at

http://www.sagemath.org/download-packages.html

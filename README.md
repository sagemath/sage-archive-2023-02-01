#  Sage: Open Source Mathematical Software

>   "Creating a Viable Open Source Alternative to
>    Magma, Maple, Mathematica, and MATLAB"

>   Copyright (C) 2005-2016 The Sage Development Team

   http://www.sagemath.org

The Sage Library is GPLv2+, and included packages have compatible OSS
licenses (see COPYING.txt). Over 400 people have contributed code to
Sage. Please see the following web page for a list:

   http://www.sagemath.org/development-map.html

In many cases, documentation for modules and functions list the
authors.


Getting Started
---------------

This README.txt contains build instructions for Sage. If you downloaded  
a binary, you do not need to do anything; just execute:  

    ./sage

from the command line. If you downloaded the sources, please read below  
on how to build Sage and work around common issues.  

If you have questions or encounter problems, please do not hesitate  
to email the sage-support mailing list:  

   http://groups.google.com/group/sage-support


Supported Platforms
-------------------

Sage fully supports several Linux distributions, recent versions of  
Mac OS X, Windows (using virtualization), as well as a number of  
Solaris and OpenSolaris releases.  

There should be no serious bugs in an officially released version of  
Sage on any of the fully supported platforms, but any major issues with  
a particular release will be documented on an errata page:  

   http://wiki.sagemath.org/errata

Ports are in progress to some other, less common platforms. The list of  
supported platforms and their current statuses are given at the  
following web page:  

   http://wiki.sagemath.org/SupportedPlatforms

If you are interested in helping port Sage to a new platform, please let  
us know at the sage-devel mailing list:  

   http://groups.google.com/group/sage-devel


Quick Instructions to Build from Source
---------------------------------------

The following steps briefly outline the process of building Sage from  
source. More detailed instructions, including how to build faster on  
multicore machines, are contained later in this README and in the  
Installation Guide:  

   http://doc.sagemath.org/html/en/installation

__1. Make sure you have the dependencies and 5 GB of free disk space.__

   >* __All Linux versions:__ gcc, make, m4, perl, ranlib, and tar.  

   >* __Fedora or RedHat systems:__ the perl-ExtUtils-MakeMaker package.  
   (install these using your package manager)  

   >* __OS X:__ Xcode. Make sure you have installed the most recent version  
   of Xcode. With recent versions of OS X (OS X Lion or later), you  
   can install Xcode for free from the App Store. For pre-Lion  
   versions of OS X, you can download Xcode from  
   http://developer.apple.com/downloads/.  

   >* With OS X, you also need to install the "command line tools". When  
   using OS X Mavericks, after installing Xcode, run this command from  
   a terminal window:  

         xcode-select --install

   >* Then click "Install" in the pop-up window.

   >* When using OS X Mountain Lion or earlier, you need to install the  
   command line tools from Xcode: run Xcode; then from the File  
   menu, choose "Preferences", then the "Downloads" tab, and then  
   "Install" the Command Line Tools.  

   >* Other platforms: See detailed instructions below.

__2. Extract the tarball:__

         tar xvf sage-*.tar

__3. cd into the Sage directory and type make:__

         cd sage-*/
         make

   That's it! Everything is automatic and non-interactive. The build  
   should work fine on all fully supported platforms. If it does not, we  
   want to know!  

   If you'd like to contribute to Sage, be sure to read the
   Developer's Guide:

   http://doc.sagemath.org/html/en/developer/index.html


Environment Variables
---------------------

There are a lot of environment variables which control the install  
process of Sage, see:  

   http://doc.sagemath.org/html/en/installation/source.html#environment-variables


Implementation
--------------

Sage has significant components written in the following languages:  
   > C/C++, Python, Cython, Lisp, and Fortran. Lisp (ECL), Python, and Cython  
   are built as part of Sage and a GNU Fortran (gfortran) binary is  
   included (OS X only), so you do not need them in order to build Sage.  


More Detailed Instructions to Build from Source
-----------------------------------------------

__1.__ Make sure you have about __5 GB of free disk space__.

__2. Install build dependencies.__

   >* __Linux:__ See quick instructions above.

   >* __OS X:__ Make sure you have XCode version >= 2.4, i.e. "gcc -v" should  
   output build >= 5363. If you don't, go to:  
      http://developer.apple.com/

   > sign up, and download the free XCode package. Only OS X >= 10.4 is  
   supported.

   >* __Solaris and OpenSolaris:__ Building Sage on these platforms is more  
   tricky than on Linux or OS X. For details on how to build Sage on  
   these platforms, see:  
      http://wiki.sagemath.org/solaris

   >* __Windows:__ Download and install VirtualBox, and then download the  
   Sage virtual appliance. For details, see:  
      http://wiki.sagemath.org/SageAppliance

   >* __NOTE:__ On some operating systems, it might be necessary to install  
   gas/as, gld/ld, gnm/nm. On most platforms, these are automatically  
   installed when you install the programs listed above.  

__3. Extract the Sage source tarball__ and cd into a directory with no  
   spaces in it. If you have a machine with 4 processors, say, type  
   the following to configure the build script to perform a parallel  
   compilation of Sage using 4 jobs:

         export MAKE="make -j4"

   >* (With 4 processors, you might also consider "-j5" or "-j6" --  
   building with more jobs than CPU cores can speed things up.)  
   You might in addition pass a "-l" flag to "make": this  
   sets a load limit, so for example if you execute  

         export MAKE="make -j4 -l5.5"

   > then "make" won't start more than one job at a time if the system  
   load average is above 5.5. See:  
   http://www.gnu.org/software/make/manual/make.html#Options-Summary  
   http://www.gnu.org/software/make/manual/make.html#Parallel.  

   > *If you want to run the test suite for each individual spkg as it is  
   installed, type:  

         export SAGE_CHECK="yes"

   > before starting the Sage build. This will run each test suite and  
   will raise an error if any failures occur. Python's test suite has  
   been disabled by default, because it causes failures on most  
   systems. To renable the Python testsuite, set the environment  
   variable SAGE_CHECK_PACKAGES to "python".

   >* To start the build, type:  

         make

__4. Wait about 20 minutes to 14 days, depending on your computer__ (it took  
   about 2 weeks to build Sage on the T-Mobile G1 Android cell phone).  

__5. Type "./sage" to try it out.__  

__6. Optional__: Type "make ptest" to test all examples in the documentation  
   (over 200,000 lines of input!) -- this takes from 10 minutes to  
   several hours. Don't get too disturbed if there are 2 to 3 failures,  
   but always feel free to email the section of logs/ptest.log that  
   contains errors to the sage-support mailing list. If there are  
   numerous failures, there was a serious problem with your build.  

__7.__ The HTML version of the documentation is built during the compilation  
   process of Sage and resides in the directory:

         $SAGE_ROOT/local/share/doc/sage/html/

   >* __Optional__: If you want to build the PDF version (requires LaTeX)  
   of the documentation, run:  

         make doc-pdf

__8. Optional__: It is recommended that you install the optional GAP
   database by typing:

         ./sage -i database_gap

   This will download the package and install it. While you're at it,
   you might install other optional packages of interest to you: type
   "./sage --optional" to get a list.

__9. Optional__: It is recommended that you have both LaTeX and the  
   ImageMagick tools (e.g. the "convert" command) installed since some  
   plotting functionality benefits from it.  

__10. Optional__: Read this if you are intending to run a Sage notebook  
    server for multiple users. For security (i.e., to run  
    "notebook(secure=True)") you want to access the server using the  
    HTTPS protocol. First, install OpenSSL and the OpenSSL development  
    headers on your system if they are not already installed. Then  
    install pyOpenSSL by building Sage and then typing  

         ./sage -i pyopenssl

   >* Note that this command requires internet access. Alternatively,  
   "make ssl" builds Sage and installs pyOpenSSL.

Installation from Docker
------------------------
__1. Installing Docker__

If you don't have docker install on your system, the following steps would guide you through.

    sudo apt-get update

    sudo apt-get install curl

    curl -fsSL https://get.docker.com/ | sh

    curl -fsSL https://get.docker.com/gpg | sudo apt-key add

If you're not willing to run a random shell script, please see the installation instructions for your distribution.

__2. Docker Build for Sage__  

There is a prebuilt docker image at https://registry.hub.docker.com/u/sagemath/sage/. Docker will automatically download and use it if you ask it to run the sagemath/sage image. So just running the following command will bring up a Sage command line without having installed anything before (except docker of course):

    docker run --rm -i -t sagemath/sage

All arguments after the sagemath/sage are passed to sage. If you want to run a sage notebook server you probably want to run it detached -d instead of interactively -i, which you can do like this:

    docker run --rm -p 127.0.0.1:8080:8080 -d -t sagemath/sage --notebook=sagenb
    docker run --rm -p 127.0.0.1:8080:8080 -d -t sagemath/sage --notebook=ipython --ip='*' --port=8080

Note that the sage notebook server runs on port 8080 inside the container. The above command maps it to 8080 outside of the container, so you can connect with your usual web browser at http://127.0.0.1:8080. But only from the computer running docker. If you want everybody on the internet to be able to connect, use -p 0.0.0.0:8080:8080 instead. You can change the port to any other port number.

     docker run --rm -i -t sagemath/sage -sh

Each docker run commant creates a new container, independent of all other containers. In particular, there is no storage shared: Each container is a completely independent file system. To keep files files after you quit the dockerized Sage, you should:

 * Not delete the container when exiting Sage (leave out the --rm command line option above).

 * Give containers a name so you know which one is which (add the --name=... commandline option).

 If you then start a container again you will get back the same files that you had when you quit that particular container. For example, say you run the container for the first time using:

        docker run --name=sage-cmdline -i -t sagemath/sage

 Then interact with Sage, create files, etc. Then quit Sage, the containers exits but is still stored (does not appear in docker ps but does appear in docker ps -a). To restart the container named sage-cmdline again, you would then just call

        docker start -i sage-cmdline

Only when you delete the container using docker rm sage-cmdline any file modifications that you made in that container are permanently lost, unless you explicitly save them with docker commit.

__Or__

You can see the containers in the docker with:

      sudo docker ps -a

And then start the sage container with :

      sudo docker start -ai ID-NAME

ID-NAME is the container ID which contains the sage. This will launch the jupyter notebook in the localhost with the port you provided before.


Following the steps would install Sage 6.7. Similar steps could be followed for the latest build, just you have to get the docker 7.5 image and build it.






Problems
--------

* If you have problems building Sage, check the Sage Installation Guide,  
and also note the following. Each separate component of Sage is  
contained in an spkg; these are stored in spkg/standard/. As each one  
is built, a build log is stored in logs/pkgs/, so you can browse these  
to find error messages. If an spkg fails to build, the whole build  
process will stop soon after, so check the most recent log files  
first, or run  

         grep -li "^Error" logs/pkgs/*

from the top-level Sage directory to find log files with error  
messages in them.  Send (a small part of) the relevant log file to the  
sage-devel mailing list, making sure to include at least some of the  
error messages; probably someone there will have some helpful  
suggestions.  


Supported Compilers
-------------------

* Sage includes a GCC (_GNU Compiler Collection_) package. In order to  
build Sage, you need a C compiler which can build GCC and its  
prerequisites. gcc version 4.0.1 or later should probably work. On  
Solaris or OpenSolaris, building with the Sun compiler should also work.  

* The GCC package in Sage is not always installed. It is determined  
automatically whether it needs to be installed. You can override this  
by setting the environment variable **SAGE_INSTALL_GCC=yes** (to force  
installation of GCC) or **SAGE_INSTALL_GCC=no** (to disable installation of  
GCC). If you don't want to install GCC, you need to have recent  
versions of gcc, g++ and gfortran; moreover, the versions must be equal.  

* There are some known problems with old assemblers, in particular when  
building the ECM package. You should ensure that your assembler  
understands all instructions for your processor. On Linux, this means  
you need a recent version of binutils; on OS X you need a recent version  
of XCode.  


Directory Layout
----------------

__Simplified directory layout (only essential files/directories):__  
```
SAGE_ROOT                 Root directory (sage-x.y.z in Sage tarball)
├── build
|   ├── deps              Dependency information of packages
|   └── pkgs              Every package is a subdirectory here
|       ├── atlas
|       ...
|       └── zn_poly
├── COPYING.txt           Copyright information
├── local                 Compiled packages are installed here
|   ├── bin               Executables
|   ├── include           C/C++ headers
|   ├── lib               Shared libraries
|   ├── share             Databases, architecture-independent data
|   └── var
|       ├── sage          List of installed packages
|       └── tmp           Temporary files when building Sage
├── logs
|   ├── dochtml.log       Log of the documentation build
|   ├── install.log       Full install log
|   └── pkgs              Build logs of individual packages
|       ├── atlas-3.10.1.p7.log
|       ...
|       └── zn_poly-0.9.p11.log
├── Makefile              Running "make" uses this file
├── README.txt            This file
├── sage                  Script to start Sage
├── src                   All of Sage source (not third-party packages)
|   ├── bin               Scripts that Sage uses internally
|   ├── doc               Sage documentation
|   └── sage              The Sage library source code
├── upstream              Source tarballs of packages
|   ├── atlas-3.10.1.tar.bz2
|   ...
|   └── zn_poly-0.9.tar.bz2
└── VERSION.txt
```
For more details, see:

    http://doc.sagemath.org/html/en/developer/coding_basics.html#files-and-directory-structure


Relocation
----------

You *should* be able to move the sage-x.y.z/ directory anywhere you  
want. If you copy the sage script or make a symbolic link to it, you  
should modify the script to reflect this (as instructed at the top of  
the script). It is best if the path to Sage does not have any spaces in  
it.  

For a system-wide installation, as root you can move the sage-x.y.z/  
directory to a system-wide directory. Afterwards, you need to start up  
Sage as root at least once prior to using the system-wide Sage as a  
normal user. See the Installation Guide for further information on  
performing a system-wide installation:

   http://doc.sagemath.org/html/en/installation/source.html#installation-in-a-multiuser-environment

If you find anything that doesn't work correctly after you moved the  
directory, please email the sage-support mailing list.  


Redistribution
--------------

Your local Sage install is almost exactly the same as any 'developer'  
install. You can make changes to documentation, source, etc., and very  
easily package the complete results up for redistribution just like we  
do.  

1. To make your own source tarball of Sage, type:  

         sage --sdist

2. To make a binary distribution with your currently installed packages,  
   type:

         sage --bdist

3. To make a binary that will run on the widest range of target  
   machines, set the SAGE_FAT_BINARY environment variable to "yes"  
   before building Sage:  

         export SAGE_FAT_BINARY="yes"
         make distclean && make
         ./sage --bdist

In all cases, the result is placed in the directory "$SAGE_ROOT/dist/".  


Changes to Included Software
----------------------------

All software included with Sage is copyrighted by the respective authors  
and released under an open source license that is __GPL version 3 or  
later__ compatible. See the file COPYING.txt for more details.  

Sources are in unmodified (as far as possible) tarballs in the  
__"$SAGE_ROOT/upstream"__ directory. The remaining description, version  
information, patches, and build scripts are in the accompanying  
__"$SAGE_ROOT/build/pkgs/<packagename>"__ directory. This directory is  
part of the Sage git repository.  

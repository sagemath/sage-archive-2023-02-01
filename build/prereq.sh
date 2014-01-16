#!/bin/sh
###########################################
## Check for prerequisite programs
###########################################

UNAME=`uname`
RELEASE=`uname -r`

# Run this from SAGE_ROOT
cd ..

echo "Starting prerequisite check."
echo "Machine: `uname -a`"
echo ""

if [ "x$SAGE_PORT" = x ]; then
    case "${UNAME}-${RELEASE}" in
    # Unsupported platforms.
    SunOS-5.[0-9])
        echo "Sage is not supported on any version of Solaris earlier than 10."
        echo "Sage has been tested on the first release of Solaris 10"
        echo "(03/2005) and works on that. Sage may or may not work with"
        echo "your version of Solaris."
        echo ""
        echo "More information can be found about Sage on Solaris"
        echo "on the Wiki at http://wiki.sagemath.org/solaris"
        echo ""
        fail=yes;;
    HP-UX*)
        echo "You are attempting to build Sage on HP's HP-UX operating system,"
        echo "which is not a supported platform for Sage yet though"
        echo "some work has been done on HP-UX. A port does not look to"
        echo "be particularly difficult. Some information can be"
        echo "found on the Sage Wiki at http://wiki.sagemath.org/HP-UX"
        echo ""
        echo "If you would like to help port Sage to HP-UX,"
        echo "please join the sage-devel discussion list - see"
        echo "http://groups.google.com/group/sage-devel"
        echo "The Sage community would also appreciate any patches you submit."
        echo ""
        fail=yes;;
    AIX*)
        echo "You are attempting to build Sage on IBM's AIX operating system,"
        echo "which is not a supported platform for Sage yet. Things may or"
        echo "may not work. If you would like to help port Sage to AIX,"
        echo "please join the sage-devel discussion list - see"
        echo "http://groups.google.com/group/sage-devel"
        echo "The Sage community would also appreciate any patches you submit."
        echo ""
        fail=yes;;
    IRIX*)
        echo "You are attempting to build Sage on SGI's IRIX operating system,"
        echo "which is not a supported platform for Sage yet. Things may or"
        echo "may not work. If you would like to help port Sage to IRIX,"
        echo "please join the sage-devel discussion list - see"
        echo "http://groups.google.com/group/sage-devel"
        echo "The Sage community would also appreciate any patches you submit."
        echo ""
        fail=yes;;
    Tru64*)
        echo "You are attempting to build Sage on HP's Tru64 operating system,"
        echo "which is not a supported platform for Sage yet. Things may or"
        echo "may not work. If you would like to help port Sage to Tru64,"
        echo "please join the sage-devel discussion list - see"
        echo "http://groups.google.com/group/sage-devel"
        echo "The Sage community would also appreciate any patches you submit."
        echo ""
        fail=yes;;
    FreeBSD*)
        echo "You are attempting to build Sage on the FreeBSD operating system,"
        echo "which is not a supported platform for Sage yet, though"
        echo "developers are working on adding FreeBSD support. Things may or"
        echo "may not work. If you would like to help port Sage to FreeBSD,"
        echo "please join the sage-devel discussion list - see"
        echo "http://groups.google.com/group/sage-devel"
        echo "The Sage community would also appreciate any patches you submit."
        echo ""
        fail=yes;;
    # Supported platforms.
    Linux*|Darwin*|SunOS*|CYGWIN*)
        # For SunOS, not all versions are supported, but we checked for
        # unsupported versions above.
        #
        # On Cygwin the result of `uname` depends on the version of
        # Windows it is running on, e.g. CYGWIN_NT-5.1 on 32 bits
        # Windows XP or CYGWIN_NT-6.1-WOW64 on 64 bits Windows 7 which
        # currently makes no difference for us.
        fail=no;;
    # Wildcard for other unsupported platforms.
    *)
        echo "You are attempting to build Sage on $UNAME,"
        echo "which is not a supported platform for Sage yet. Things may or"
        echo "may not work. If you would like to help port Sage to $UNAME,"
        echo "please join the sage-devel discussion list - see"
        echo "http://groups.google.com/group/sage-devel"
        echo "The Sage community would also appreciate any patches you submit."
        echo ""
        fail=yes;;
    esac

    # Test versions of "tar" and "make".
    have_gnu_tar=`tar --version 2>&1 | grep GNU`
    have_bsd_tar=`tar --version 2>&1 | grep bsdtar`
    have_gnu_make=`make --version 2>&1 | grep GNU`

    if [ -z "${have_gnu_tar}${have_bsd_tar}" ] ; then
         echo "You MUST use either the GNU or BSD version of tar, as some parts of the"
         echo "Sage source code are contained in tar files, which were produced with GNU"
         echo "tar. These are not POSIX compliant tar files (unless they were created with"
         echo "the --posix option), so can not always be extracted with other 'tar'"
         echo "programs."
         echo ""
         echo "Please ensure that a GNU or BSD version of tar is first in your path."
         echo ""
         fail=yes

         if [ -f "/usr/sfw/bin/gtar" ] ; then
            echo "There is a copy of GNU tar, but renamed to 'gtar' in"
            echo "the directory /usr/sfw/bin. It is suggested that you copy"
            echo "/usr/sfw/bin/gtar to a different directory, renaming it"
            echo "'tar' instead of gtar. Then ensure that directory is in"
            echo "you path before /usr/bin, so whenever you type 'tar' you"
            echo "run the GNU version."
            echo ""
            echo "Executing something like"
            echo ""
            echo '$ mkdir $HOME/bins-for-sage'
            echo '$ cp /usr/sfw/bin/gtar $HOME/bins-for-sage/tar'
            echo '$ export PATH=$HOME/bins-for-sage:$PATH'
            echo ""
            echo "will put the GNU version of 'tar' in your path before /usr/bin/."
            echo ""
            echo "If you type 'tar --version', it should be obvious"
            echo "whether the GNU version of 'tar' is first in your path or"
            echo "not."
            echo ""
        fi
    fi

    if [ -z "${have_gnu_make}" ] ; then
         echo "Due to numerous GNUisms, you MUST use the GNU version of 'make'"
         echo "to build Sage. Other versions will not work."
         echo ""
         fail=yes

         if [ -f "/usr/sfw/bin/gmake" ] ; then
            echo "There is a copy of GNU make, but renamed to 'gmake' in"
            echo "the directory /usr/sfw/bin. It is suggested that you copy"
            echo "/usr/sfw/bin/gmake to a different directory, renaming it"
            echo "'make' instead of gmake. Then ensure that directory is in"
            echo "you path before /usr/ccs/bin, so whenever you type 'make' you"
            echo "run the GNU version."
            echo ""
            echo "Executing something like"
            echo ""
            echo '$ mkdir $HOME/bins-for-sage'
            echo '$ cp /usr/sfw/bin/gmake $HOME/bins-for-sage/make'
            echo '$ export PATH=$HOME/bins-for-sage:$PATH'
            echo ""
            echo "will put the GNU version of 'make' in your path before /usr/ccs/bin/."
            echo ""
            echo "If you type 'make --version', it should be obvious"
            echo "whether the GNU version of 'make' is first in your path or"
            echo "not."
            echo ""
        fi
    fi

    if [ x$fail = xyes ]; then
        echo "To get past this message and try building Sage anyway,"
        echo "export the variable SAGE_PORT to something non-empty."
        exit 1
    fi
fi

# Make sure the configure script exists
$MAKE bootstrap || exit $?

# A reasonably sophisticated test is performed in a configure
# script, which checks compilers exist, their version numbers,
# the fact GNU and non-GNU compilers are not mixed etc.
./configure --disable-maintainer-mode $PREREQ_OPTIONS
if [ $? -ne 0 ]; then
    echo >&2 "You do not have all of the prerequisites needed to build Sage"
    echo >&2 "from source. See the errors above."
    if [ "x$SAGE_PORT" = x ]; then
        echo >&2 "If you would like to try the build anyway (to help porting)"
        echo >&2 "export the variable 'SAGE_PORT' to something non-empty."
        exit 1
    else
        echo >&2 "However, since 'SAGE_PORT' is set, we will try to build anyway."
    fi
fi

if pwd | grep >/dev/null " "; then
    echo ""
    echo "*********************************************************"
    echo ""
    echo " ERROR: Sage will (probably) not build correctly if there is a space"
    echo " in the path to the current directory."
    echo ""
    echo " Path = `pwd`"
    echo ""
    echo "*********************************************************"
    exit 1
fi


###########################################################################
# (OS X only)
# Sage will probably not build at all if either Fink or MacPorts can be
# found, and the error messages can be extremely confusing.  Even if it does
# build, the product will probably be wrong.  This runs a basic check to
# find them. Once the Sage build process is perfected, this won't be necessary.
# dphilp 15/9/2008
###########################################################################

if [ "$UNAME" = "Darwin" ]; then
    # Warning: xcodebuild does not seem to be maintained in Xcode 4.3
    # or later, so do not rely on the variable XCODE_VERS with OS X
    # 10.7 or later.
    XCODE_VERS=`xcodebuild -version 2> /dev/null | grep Xcode | sed -e 's/[A-Za-z ]//g'`
    if [ -z $XCODE_VERS ]; then
        XCODE_VERS="2"
    fi
    XCODE_VERS_MAJOR=`echo $XCODE_VERS | cut '-d.' -f1`
    DARWIN_VERSION=`uname -r | cut '-d.' -f1`
    echo "***************************************************"
    echo "***************************************************"
    if [ $DARWIN_VERSION -gt 10 ]; then
        echo "You are using OS X Lion (or later)."
        echo "You are strongly advised to install Apple's latest Xcode"
        echo "unless you already have it. You can install this using"
        echo "the App Store. Also, make sure you install Xcode's"
        echo "Command Line Tools -- see Sage's README.txt."
    elif [ $XCODE_VERS_MAJOR -gt 2 ]; then
        echo "You are using Xcode version $XCODE_VERS."
        echo "You are strongly advised to install Apple's latest Xcode"
        echo "unless you already have it. You can download this from"
        echo "http://developer.apple.com/downloads/."
        echo "If using Xcode 4.3 or later, make sure you install Xcode's"
        echo "Command Line Tools -- see Sage's README.txt."
    else
        echo "You are using Xcode version 1 or 2"
        echo "WARNING: You are strongly advised to install the"
        echo "latest version of Apple's Xcode for your platform,"
        echo "unless you already have it."
        if [ $DARWIN_VERSION -eq 10 ]; then
           echo "Probably you need Xcode 3.2.6"
        elif [ $DARWIN_VERSION -eq 9 ]; then
           echo "Probably you need Xcode 3.1.4"
        elif [ $DARWIN_VERSION -lt 9 ]; then
           echo "Probably you need Xcode 2.5"
        fi
    fi
    echo "***************************************************"
    echo "***************************************************"
    SAGE_ABORT=no
    # Try to find ports automatically.
    PORTS_PATH=`which port`
    if [ -f "`which port`" ]; then
        echo "Found MacPorts in " $PORTS_PATH
        SAGE_ABORT="yes"
    fi

    # Try to find fink automatically.
    FINK_PATH=`which fink`
    if [ -f "`which fink`" ]; then
        echo "Found Fink in " $FINK_PATH
        SAGE_ABORT="yes"
    fi

    if [ "$SAGE_ABORT" = "yes" ]; then
        echo ""
        echo "*********************************************************"
        echo ""
        echo "Found either MacPorts or Fink in your PATH, which potentially wrecks the Sage build process."
        if [ "$SAGE_COMPILE_DESPITE_PORTS_AND_FINK" ]; then
            echo "Continuing because SAGE_COMPILE_DESPITE_PORTS_AND_FINK is set."
            echo ""
            echo "*********************************************************"
            echo ""
        else
            echo "You should make sure MacPorts and Fink cannot be found.  Either:"
            echo "(1) rename /opt/local and /sw, or"
            echo "(2) change PATH and DYLD_LIBRARY_PATH"
            echo "(Once Sage is built, you can restore them.)"
            echo ""
            echo "*********************************************************"
            echo ""
            exit 1
        fi
    fi
fi

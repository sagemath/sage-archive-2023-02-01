# Very loosely based on configure.ac in prereq-0.3 written by William Stein.
# Version 0.7 written by David Kirkby, released under the GPL version 2.
# in January 2010

# If you are going to update this, please stick the recommended layout
# in the autoconf manual - i.e.

# First check for programs
# Next check for libraries
# Next check for header files.
# Next check for types.
# Next check for structures.
# Next check compiler characteristics
# Next check for library functions
# Next check for system services

# We do not want to update this configure.ac each time a new Sage is released
# so I am putting 'x.y.z'

AC_INIT([Sage],[x.y.z], sage-devel@googlegroups.com)
AC_PREREQ([2.60])

AC_COPYRIGHT([GPL version 2 or (at your option) any later version])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([config])

# We don't really use automake, but configure needs aclocal and the
# automake library files (like config.guess).
AM_INIT_AUTOMAKE([1.9.6 foreign])

# Allow "configure --disable-maintainer-mode" to disable timestamp checking
AM_MAINTAINER_MODE([enable])

#---------------------------------------------------------
# We need to run this configure script with bash
if test -z "$BASH_VERSION$CONFIG_SHELL"
then
    CONFIG_SHELL=bash
    export CONFIG_SHELL
    if $CONFIG_SHELL -c "exit 0"
    then
        exec $CONFIG_SHELL $0 "$@"
    else
        AC_MSG_NOTICE([The 'bash' shell is needed to build AC_PACKAGE_NAME])
        AC_MSG_NOTICE([All modern systems will have the 'bash' shell installed somewhere])
        if test -d /opt/OpenSource/bin
        then
           AC_MSG_NOTICE([On HP-UX you may try adding /opt/OpenSource/bin to your path])
        fi
        if test -d /opt/pware/bin
        then
           AC_MSG_NOTICE([On AIX you may try adding /opt/pware/bin to your path])
        fi
        AC_MSG_ERROR(['bash' not found])
    fi
fi

#---------------------------------------------------------
buggy_gcc_version1="4.0.0"
minimum_gcc_version_for_no_hassle="4.0.1"
minimum_gcc_version_for_debugging_purposes="3.4.0"

#---------------------------------------------------------
#   Process options and environment variables

# Check --enable-compiler-checks (yes by default)
AC_ARG_ENABLE([compiler-checks], [Check versions and presence of C, C++ and Fortran compilers (default: yes)],
    [enable_compiler_checks=$enableval], [enable_compiler_checks=yes])

# Import environment variables.
source src/bin/sage-env || AC_MSG_ERROR([failed to source sage-env])

#---------------------------------------------------------
#   Hardware architecture
#
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

#---------------------------------------------------------
# Check some programs needed actually exist.
AC_CHECK_PROG(found_ar, ar, yes, no)
if test x$found_ar != xyes
then
    AC_MSG_NOTICE([Sorry, the 'ar' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems it can be found in /usr/ccs/bin ])
    AC_MSG_NOTICE(['ar' is also part of the GNU 'binutils' package.])
    AC_MSG_ERROR([Exiting, as the archiver 'ar' can not be found.])
fi

AC_CHECK_PROG(found_m4, m4, yes, no)
if test x$found_m4 != xyes
then
    AC_MSG_NOTICE([Sorry, the 'm4' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems it can be found in /usr/ccs/bin])
    AC_MSG_NOTICE([See also http://www.gnu.org/software/m4/])
    AC_MSG_ERROR([Exiting, as the macro processor 'm4' can not be found.])
fi 

AC_CHECK_PROG(found_ranlib, ranlib, yes, no)
if test x$found_ranlib != xyes
then
    AC_MSG_NOTICE([Sorry, the 'ranlib' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems it can be found in /usr/ccs/bin ])
    AC_MSG_NOTICE(['ranlib' is also part of the GNU 'binutils' package.])
    AC_MSG_ERROR([Exiting, as 'ranlib' can not be found.])
fi

AC_CHECK_PROG(found_strip, strip, yes, no)
if test x$found_strip != xyes
then
    AC_MSG_NOTICE([Sorry, the 'strip' command must be in the path to build AC_PACKAGE_NAME])
    AC_MSG_NOTICE([On some systems 'strip' can be found in /usr/ccs/bin ])
    AC_MSG_NOTICE(['strip' is also part of the GNU 'binutils' package.])
    AC_MSG_ERROR([Exiting, as 'strip' can not be found.])
fi

# Check for Latex, the use of which is less important in Sage than
# it used to be, as it was at one time required to build any documentation
# but this is no longer so.
AC_CHECK_PROG(found_latex, latex, yes, no)
if test x$found_latex != xyes
then
    AC_MSG_WARN([You do not have 'latex', which is recommended, but not])
    AC_MSG_WARN([required. Latex is only really used for building pdf])
    AC_MSG_WARN([documents and for %latex mode in the AC_PACKAGE_NAME notebook.])
fi

# Check that perl is available, with version 5.8.0 or later.
# Some packages need perl, however it is not clear whether Sage really
# requires version >= 5.8.0.  The R package *used* to require it, but
# not anymore. -- Jeroen Demeyer
AC_PATH_PROG([PERL],[perl])
AX_PROG_PERL_VERSION([5.8.0],[],[
    AC_MSG_ERROR([Exiting, since AC_PACKAGE_NAME requires perl-5.8.0 or later])
])

# To build Python on multi-arch Debian-based systems, we need
# dpkg-architecture. Since we need dpkg-architecture to determine
# whether we're on a multi-arch system and require dpkg-architecture,
# we simply require it always on Debian-based systems.
AC_CHECK_PROG(found_dpkg, dpkg, yes, no)
AC_CHECK_PROG(found_dpkg_arch, dpkg-architecture, yes, no)
if test x$found_dpkg = xyes && test x$found_dpkg_arch = xno
then
    AC_MSG_NOTICE([You do not have 'dpkg-architecture', which is required to build])
    AC_MSG_NOTICE([Python on multi-arch Debian-based systems. This includes all recent])
    AC_MSG_NOTICE([Debian and Ubuntu systems. You can install this with:])
    AC_MSG_NOTICE([  sudo apt-get install dpkg-dev])
    AC_MSG_ERROR([Exiting, since AC_PACKAGE_NAME requires dpkg-architecture on Debian])
fi

#---------------------------------------------------------
#       C/C++/Fortran compilers

# First check for programs we need.
AC_LANG(C)
AC_PROG_CC()
if test "x$CC" = x
then
    AC_MSG_ERROR([Exiting, since a C compiler was not found.])
fi
AC_PROG_CPP()
AC_PROG_CXX()
AC_PROG_FC()
if test $enable_compiler_checks = yes
then
    if test "x$CXX" = x
    then
        AC_MSG_ERROR([Exiting, since a C++ compiler was not found.])
    fi
    if test "x$FC" = x
    then
        AC_MSG_ERROR([Exiting, since a Fortran compiler was not found.])
    fi
fi

# As suggested in the autoconf manual, the characteristics of the compilers
# are checked later.


# Next one should check for header files.
# complex.h is one that might not exist on older systems.
AC_CHECK_HEADER([complex.h],[],[
    AC_MSG_ERROR([Exiting, since you do not have the 'complex.h' header file.])
])

# Next one should check for types.
# None needed

# Next one should check for structures.
# None needed

# Next one should check for compiler characterists.

# Check that we can compile C99 code
AC_PROG_CC_C99()
if test $enable_compiler_checks = yes
then
    if test "x$ac_cv_prog_cc_c99" = xno
    then
        AC_MSG_ERROR([Exiting, as your C compiler cannot compile C99 code])
    fi
fi

# Check the Fortran compiler accepts free-format source code
# (as opposed to the older fixed-format style from Fortran 77).
# This helps verify the compiler works too, so if some idiot
# sets SAGE_FORTRAN to /usr/bin/ls, we will at least know it's
# not a working Fortran compiler.
if test $enable_compiler_checks = yes
then
   # see http://www.gnu.org/software/hello/manual/autoconf/Fortran-Compiler.html
   AC_FC_FREEFORM([],
   [
   AC_MSG_NOTICE([Your Fortran compiler does not accept free-format source code])
   AC_MSG_NOTICE([which means the compiler is either seriously broken, or])
   AC_MSG_NOTICE([is too old to build Sage.])
   AC_MSG_ERROR([Exiting, as the Fortran compiler is not suitable])
   ])
fi

if test $enable_compiler_checks = yes
then
    # Check that all compilers (C, C++, Fortan) are either all GNU
    # compiler or all non-GNU compilers.  If not, there is a problem, as
    # mixing GNU and non-GNU compilers is likely to cause problems.
    if test x$GCC = xyes && test x$GXX != xyes
    then
        AC_MSG_NOTICE([You are trying to use gcc but not g++])
        AC_MSG_ERROR([The mixing of GNU and non-GNU compilers is not permitted])
    fi
    if test x$GXX = xyes && test x$ac_cv_fc_compiler_gnu != xyes
    then
        AC_MSG_NOTICE([You are trying to use g++ but not gfortran])
        AC_MSG_ERROR([The mixing of GNU and non-GNU compilers is not permitted])
    fi
    if test x$ac_cv_fc_compiler_gnu = xyes && test x$GCC != xyes
    then
        AC_MSG_NOTICE([You are trying to use gfortran but not gcc])
        AC_MSG_ERROR([The mixing of GNU and non-GNU compilers is not permitted])
    fi
fi

# The following tests check the version of the compilers (if GNU)
# are all the same. If non-GNU compilers are used, then no such
# checks are performed.
if test $enable_compiler_checks = yes
then
if test x$GCC = xyes
then
    # Thank you to Andrew W. Nosenko andrew.w.nosenko@gmail.com
    # who answered my query about testing of gcc versions on
    # the autoconf@gnu.org mailing list.
    # AS_VERSION_COMPARE(ver-1, ver-2, [action-if-less], [action-if-eq], [action-if-greater])
    AX_GCC_VERSION
    AX_GXX_VERSION
    AS_VERSION_COMPARE([$GCC_VERSION], [$GXX_VERSION],
    [
    AC_MSG_NOTICE([gcc ($GCC_VERSION) and g++ ($GXX_VERSION) are not the same version])
    AC_MSG_NOTICE([which they must be. Check your setting of CC and CXX])
    AC_MSG_ERROR([Exiting, since the C and C++ compilers have different versions])
    ],[],[
    AC_MSG_NOTICE([gcc ($GCC_VERSION) and g++ ($GXX_VERSION) are not the same version])
    AC_MSG_NOTICE([which they must be. Check your setting of CC and CXX])
    AC_MSG_ERROR([Exiting, since the C and C++ compilers have different versions])])

    # In the paragraph below, 'gfortran' is used to indicate the GNU Fortran
    # compiler, though it might be called something else.

    # It's not easily possible to determine the Fortran version, as
    # gfortran -dumpversion did not until GCC 4.5 return just the
    # the version number, but the same as gfortran --version
    # for example:

    # drkirkby@hawk:~$ gcc -dumpversion
    # 4.3.4

    # drkirkby@hawk:~$ g++ -dumpversion
    # 4.3.4

    # drkirkby@hawk:~$ gfortran -dumpversion
    # GNU Fortran (GCC) 4.3.4
    # Copyright (C) 2008 Free Software Foundation, Inc.
    # GNU Fortran comes with NO WARRANTY, to the extent permitted by law.
    # You may redistribute copies of GNU Fortran
    # under the terms of the GNU General Public License.
    # For more information about these matters, see the file named COPYING

    # This behaviour is fixed in the gcc 4.5 branch. Since we need to
    # support older versions of the compiler, we can't change this.

    # But I would expect that the version will be on the output
    # of the compiler followed by -dumpversion (e.g. fortran -dumpversion)

    # So we grep for the known gcc version on the output of gfortran -dumpversion.
    # and hope we find the same string. If so, they are almost certainly
    # the same version.
    fortran_version_string="`$FC -dumpversion | grep $GCC_VERSION 2>&1`"

    if test "x$fortran_version_string" = x
    then
        AC_MSG_NOTICE([Although gcc and g++ are both version $GCC_VERSION])
        AC_MSG_NOTICE([the Fortran compiler $FC is some other version.])
        AC_MSG_NOTICE([The output from $FC --version is below.])
        echo ""
        $FC --version 2>&1
        echo ""
        AC_MSG_ERROR([Exiting, since the Fortran compiler is not the same version as the C and C++ compilers])
    else
        AC_MSG_NOTICE([Excellent, the C, C++ and Fortran compilers are all GCC $GCC_VERSION])
    fi

    # Exit if the version of GCC is known to be too old, and so old
    # we have no intension whatsoever of trying to make Sage work with it.
    AS_VERSION_COMPARE([$GCC_VERSION], [$minimum_gcc_version_for_debugging_purposes],[
        AC_MSG_NOTICE([GCC $GCC_VERSION is too old and can not build AC_PACKAGE_NAME. ])
        AC_MSG_NOTICE([Please use a GCC of at least $minimum_gcc_version_for_no_hassle ])
        AC_MSG_NOTICE([There are no plans whatsoever to support GCC $GCC_VERSION])
        AC_MSG_ERROR([Exiting, due to the use of a version of GCC that is too old])
    ], [],[])

    # Exit if Sage is *precisely* version $buggy_gcc_version1, as that is very buggy.
    # If any later versions of GCC are found to be buggy (which would never surprise me)
    # We can easily add a buggy_gcc_version2 and repeat the next 5 lines.
    # At the time of writing (28th September 2009) that is version 4.0.0,
    # but rather than hard-code that, it is set as a variable.
    AS_VERSION_COMPARE([$GCC_VERSION], [$buggy_gcc_version1],[],[
        AC_MSG_NOTICE([GCC $buggy_gcc_version1  is very buggy and can not build AC_PACKAGE_NAME.])
        AC_MSG_NOTICE([Please use a gcc of at least $minimum_gcc_version_for_no_hassle.])
        AC_MSG_ERROR([Exiting, due to the use of a version of GCC that is too buggy])
     ],[])

    # Issue a warning if gcc is too old to work with all packages, but
    # someone wants to try getting one or more packages work with
    # an earlier gcc. At the time of writing, (28th Sept 2009), ratpoints
    # has such an issue, requiring version 4.0.1, but we would like to
    # get it to work with version 3.4.x
    AS_VERSION_COMPARE([$GCC_VERSION], [$minimum_gcc_version_for_no_hassle],[
        AC_MSG_NOTICE([******************************************************])
        AC_MSG_NOTICE([******************************************************])
        AC_MSG_NOTICE([******************************************************])
        AC_MSG_NOTICE([GCC $GCC_VERSION is too old and can not build AC_PACKAGE_NAME.])
        AC_MSG_NOTICE([Please use a gcc of at least $minimum_gcc_version_for_no_hassle])
        AC_MSG_NOTICE([if you just want AC_PACKAGE_NAME to build without problems.])
        AC_MSG_NOTICE([])
        if test "${SAGE_USE_OLD_GCC+set}" = set; then
           AC_MSG_NOTICE([Since the variable SAGE_USE_OLD_GCC was set, the])
           AC_MSG_NOTICE([build will continue, but it will fail without changes])
           AC_MSG_NOTICE([to the Sage source code. You can be 100% sure of that.])
        else
           AC_MSG_NOTICE([If you want to try building Sage with a GCC 3.4.x     ])
           AC_MSG_NOTICE([with a view to debugging the problems which stop it   ])
           AC_MSG_NOTICE([working on a gcc 3.4 series compiler, set the         ])
           AC_MSG_NOTICE([environment variable SAGE_USE_OLD_GCC to something non])
           AC_MSG_NOTICE([empty. But if you just want Sage to work, upgrade GCC ])
           AC_MSG_NOTICE([******************************************************])
           AC_MSG_NOTICE([******************************************************])
           AC_MSG_ERROR([Exiting, due to the use of a version of GCC that is too old])
        fi
    ],
[
AC_MSG_NOTICE([Good, gcc and g++ are are just new enough as GCC $minimum_gcc_version_for_no_hassle])
AC_MSG_NOTICE([is the minimum version needed needed to build AC_PACKAGE_NAME])
]
,
[
AC_MSG_NOTICE([Excellent, GCC $GCC_VERSION is later than the minimum])
AC_MSG_NOTICE([needed to build AC_PACKAGE_NAME, which is GCC version $minimum_gcc_version_for_no_hassle])
])

    # AS_VERSION_COMPARE(ver-1, ver-2, [action-if-less], [action-if-eq], [action-if-greater])
else
    AC_MSG_WARN([You have a non-GNU compiler, but AC_PACKAGE_NAME has never been built])
    AC_MSG_WARN([successfully with any compiler other than GCC, despite])
    AC_MSG_WARN([some attempts made on Solaris to use Sun's compiler, which])
    AC_MSG_WARN([produces faster code than GCC. However, the AC_PACKAGE_NAME developers])
    AC_MSG_WARN([want AC_PACKAGE_NAME to work with other compilers, so please try.])
    AC_MSG_WARN([The AC_PACKAGE_NAME developers would welcome any feedback you can give.])
    AC_MSG_WARN([Please visit http://groups.google.com/group/sage-devel])
    AC_MSG_WARN([If you just want to use AC_PACKAGE_NAME, we suggest the use of])
    AC_MSG_WARN([GCC of at least version $minimum_gcc_version_for_no_hassle])
fi
fi  # test $enable_compiler_checks = yes

# Testing for library functions
# First check for something that should be in any maths library (sqrt).
AC_CHECK_LIB(m,sqrt,[],[
                       AC_MSG_NOTICE([This system has no maths library installed.])
                       # On AIX this is not installed by default - strange as that might seem.
                       # but is in a fileset bos.adt.libm. However, the fileset bos.adt
                       # includes other things that are probably useful.
                       if test "x`uname`" = 'xAIX'
                       then
                          AC_MSG_NOTICE([On AIX, libm is contained in the bos.adt.libm fileset.])
                          AC_MSG_NOTICE([Actually, we recommend to install the complete bos.adt fileset.])
                          AC_MSG_NOTICE([This needs to be performed by a system administrator.])
                       fi
                       AC_MSG_ERROR([Exiting, since a maths library was not found.])
                       ])

# Next check for system services

# Test for old versions of Solaris, and issue a warning if too old.
# Advise people to upgrade their Solaris version, but issue a warning
# in the unlikely event their hardware is too old to run Solaris 10
# or later.

case "$build_os" in
    solaris2.?|solaris2.?.*)
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([You are using an older version of Solaris.               ])
    AC_MSG_WARN([AC_PACKAGE_NAME is only supported on Solaris 10 or later, though it ])
    AC_MSG_WARN([may work on older releases.                              ])
    if test "x`uname -m`" = 'xsun4m'  ||  test "x`uname -m`" = 'xsun4c'
    then
        AC_MSG_WARN([Unfortunately, your hardware is too old to update to     ])
        AC_MSG_WARN([Solaris 10, so you may experience problems building AC_PACKAGE_NAME.])
        AC_MSG_WARN([In any case, it would be very slow on such old hardware. ])
    else
        AC_MSG_WARN([You can download the latest Solaris version free from])
        AC_MSG_WARN([http://www.sun.com/software/solaris/ and install it])
        AC_MSG_WARN([Unless you have other reasons for particularly wishing])
        AC_MSG_WARN([to an older release of Solaris, upgrading would be  ])
        AC_MSG_WARN([the easiest way to run AC_PACKAGE_NAME on your hardware.       ])
    fi
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([If you would like to help port AC_PACKAGE_NAME to this older release])
    AC_MSG_WARN([of Solaris, then please join the sage-devel discussion   ])
    AC_MSG_WARN([list, where development issues are discussed. See: ])
    AC_MSG_WARN([http://groups.google.com/group/sage-devel])
    AC_MSG_WARN([The AC_PACKAGE_NAME community would also appreciate any patches you ])
    AC_MSG_WARN([submit, especially for Solaris 8 and 9.                  ])
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([*********************************************************])
    ;;
esac

# Test for old versions of Darwin. It's hard to know exactly what each version
# will show for build_os, but a search with Google for config.log files show
# some values that can not be found on
# http://en.wikipedia.org/wiki/Darwin_(operating_system) (retrieved 10/24/2009)
# Hence the following are educated guesses, based on a search of config.log files
# posted by people

case "$build_os" in
    darwin5.?|darwin5.?.?|darwin6.?|darwin6.?.?|darwin7.?|darwin7.?.?)
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([AC_PACKAGE_NAME has never been built on OS X 10.3 (Panther)])
    AC_MSG_WARN([or earlier. The oldest version of OS X successfully used ])
    AC_MSG_WARN([is OS X version 10.4 (Tiger). You might consider updating])
    AC_MSG_WARN([your version of OS X if your hardware permits this, but  ])
    AC_MSG_WARN([Apple charges for upgrades of OS X.])
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([If you would like to help port AC_PACKAGE_NAME to this older release])
    AC_MSG_WARN([of OS X, then please join the sage-devel discussion   ])
    AC_MSG_WARN([list, where development issues are discussed. See: ])
    AC_MSG_WARN([http://groups.google.com/group/sage-devel])
    AC_MSG_WARN([*********************************************************])
    AC_MSG_WARN([*********************************************************])
    ;;
esac

dnl AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([build/Makefile-auto])
AC_CONFIG_MACRO_DIR([m4])

AC_OUTPUT


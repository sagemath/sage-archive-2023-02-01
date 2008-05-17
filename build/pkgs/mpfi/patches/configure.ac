# initialisation autoconf
AC_INIT(mpfi,1.3.4-RC4,[Nathalie.Revol@ens-lyon.fr,Fabrice.Rouillier@inria.fr])

# initialisation automake
AM_INIT_AUTOMAKE(mpfi,1.3.4-RC4)

# fichier contenant les definitions de toutes
# les macros definies par Autoconf
AC_CONFIG_HEADERS([mpfi_config.h])

# pour la construction d'une librairie dynamique
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# traitement des sous-repertoires
MPFI_SUBDIRS="src tests doc"
AC_SUBST(MPFI_SUBDIRS)

# ajout des options sur la ligne de commande du configure
# Initialisation des repertoires GMP
AC_ARG_WITH(gmp-dir,
	    AC_HELP_STRING([--with-gmp-dir=DIR : GMP installation directory]),
            with_gmp_include=$withval/include
	    with_gmp_lib=$withval/lib
            with_mpfr_include=$withval/include
	    with_mpfr_lib=$withval/lib)
AC_ARG_WITH(gmp-libpath,
	    AC_HELP_STRING([--with-gmp-libpath=DIR : GMP libraries directory]),
            with_gmp_lib=$withval)
AC_ARG_WITH(gmp-incpath,
	    AC_HELP_STRING([--with-gmp-incpath=DIR : GMP include directory]),
	    with_gmp_include=$withval)
# Initialisation des repertoires de MPFR
AC_ARG_WITH(mpfr-dir,
	    AC_HELP_STRING([--with-mpfr-dir=DIR : MPFR installation directory]),
            with_mpfr_include=$withval/include
	    with_mpfr_lib=$withval/lib)
AC_ARG_WITH(mpfr-libpath,
            AC_HELP_STRING([--with-mpfr-libpath=DIR : MPFR library directory]),
            with_mpfr_lib=$withval)
AC_ARG_WITH(mpfr-incpath,
            AC_HELP_STRING([--with-mpfr-incpath=DIR : MPFR include directory]),
            with_mpfr_include=$withval)

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([string.h])

# Add include directory for GMP and MPFR
if test -d "$with_gmp_include"; then
  CPPFLAGS="$CPPFLAGS -I$with_gmp_include"
else
  gmp_include=
fi
if test -d "$with_mpfr_include"; then
  CPPFLAGS="$CPPFLAGS -I$with_mpfr_include"
else
  mpfr_include=
fi

# Checks for GMP
AC_MSG_CHECKING(for gmp.h)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
]], [[]])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)
     AC_MSG_ERROR([gmp.h may be missing ${with_gmp_include:+in $with_gmp_include}])
])

AC_MSG_CHECKING(for valid GMP)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#if (__GNU_MP_VERSION < 4)
# error "min: GMP 4.0.0"
#endif
]], [[]])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)
     AC_MSG_ERROR([GMP 4.0.0 min required])
])

# Checks for MPFR
AC_MSG_CHECKING(for mpfr.h)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#include "mpfr.h"
]], [[]])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)
     AC_MSG_ERROR([mpfr.h may be missing ${with_mpfr_include:+in $with_mpfr_include}])
])

# Checks for MPFR lib (Before GMP!)
if ` test "$with_mpfr_lib" `
then
  AC_MSG_CHECKING(MPFR library)
        if test -r "$with_mpfr_lib/libmpfr.so"
        then
          LDADD="$LDADD -L$with_gmp_lib -lmpfr"
        else
           if test -r "$with_mpfr_lib/libmpfr.dylib"
           then
             LDADD="$LDADD -L$with_gmp_lib -lmpfr"
           else
             if test -r "$with_mpfr_lib/libmpfr.dll.a"
             then
               LDADD="$LDADD -L$with_gmp_lib -lmpfr"
             else
               AC_MSG_ERROR([$with_mpfr_lib/libmpfr.so or libmpfr.dylib not found])
             fi
           fi
        fi
  AC_MSG_RESULT(yes)
else
  AC_CHECK_LIB(mpfr, main, , AC_MSG_ERROR([Library MPFR not found]))
fi

# Checks for GMP lib
if ` test "$with_gmp_lib" `
then
  AC_MSG_CHECKING(GMP library)
        if test -r "$with_gmp_lib/libgmp.so"
        then
          LDADD="$LDADD -L$with_gmp_lib -lgmp"
        else
          if test -r "$with_gmp_lib/libgmp.dylib"
          then
            LDADD="$LDADD -L$with_gmp_lib -lgmp"
          else
            if test -r "$with_gmp_lib/libgmp.dll.a"
            then
               LDADD="$LDADD -L$with_gmp_lib -lgmp"
            else
               AC_MSG_ERROR([$with_gmp_lib/libgmp.so or libgmp.dylib not found])
            fi
          fi
        fi
  AC_MSG_RESULT(yes)
else
  AC_CHECK_LIB(gmp, main, , AC_MSG_ERROR([Library GMP not found]))
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([memset])
AC_CANONICAL_HOST

AC_SUBST(INCLUDES)
AC_SUBST(LDADD)
AC_SUBST(LDFLAGS)
AC_SUBST(CFLAGS)

# generation des fichiers suivants par autoconf
AC_OUTPUT([Makefile src/Makefile tests/Makefile doc/Makefile])

# pour la configurer la distrib
# configure --with-gmp-dir=/users/rouillie/DEVELOP/i686_Linux_libc6/GMP --prefix=/tmp/mpfi

# pour fabriquer la distrib
# aclocal
# autoconf
# autoheader
# automake -a

# pour faire un tar de la distrib (apres l'avoir fabriquee et configuree)
# make dist

# pour ajouter des fichiers source : voir les Makefile.am
# dans les sous-repertoires

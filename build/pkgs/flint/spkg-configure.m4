SAGE_SPKG_CONFIGURE([flint], [
    AC_ARG_WITH([flint],
    [AS_HELP_STRING([--with-flint=system],
        [use the system FLINT, if possible (default)])]
    [AS_HELP_STRING([--with-flint=install],
        [use the Sage SPKG for FLINT])])

dnl Just part the options here
    case "$with_flint" in
        system) ;;
        install) ;;
        "") with_flint=system;;
        *)
            AC_MSG_ERROR([allowed values for --with-flint are system and install]);;
    esac
    
    case "$with_flint" in
        system)
               AC_CHECK_HEADER(flint/flint.h, [], [sage_spkg_install_flint=yes])
          dnl fmpz_mat_is_hadamard appears in Flint 2.5.0
               AC_SEARCH_LIBS([fmpz_mat_is_hadamard], [flint], [], [sage_spkg_install_flint=yes])
          dnl check that NTL is linked in
               AC_SEARCH_LIBS([fmpz_poly_get_ZZX], [flint], [break], [sage_spkg_install_flint=yes])
          dnl 
               AC_MSG_CHECKING([that GC is not enabled in Flint... ])
               AC_RUN_IFELSE([
                  AC_LANG_PROGRAM([[#include <flint/flint.h>]], [[return HAVE_GC;]])],
                  [AC_MSG_RESULT([GC not enabled. Good.])],
		  [AC_MSG_RESULT([GC enabled. Incompatible with Sage.])
		   sage_spkg_install_flint=yes])

            if test x$sage_spkg_install_flint = xyes; then
               AC_SUBST(SAGE_FLINT_PREFIX, ['$SAGE_LOCAL'])
               AC_MSG_RESULT([using Sage's flint SPKG])
            else
               AC_SUBST(SAGE_FLINT_PREFIX, [''])
               AC_MSG_RESULT([using flint library from the system])
            fi
            ;;
        install)
            sage_spkg_install_flint=yes
            AC_SUBST(SAGE_FLINT_PREFIX, ['$SAGE_LOCAL'])
            AC_MSG_RESULT([using Sage's flint SPKG])
            ;;
    esac

])

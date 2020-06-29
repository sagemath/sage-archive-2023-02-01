SAGE_SPKG_CONFIGURE([freetype], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_LIBPNG])
    AC_MSG_CHECKING([Installing libpng? ])
    if test x$sage_spkg_install_libpng = xyes; then
      AC_MSG_RESULT([yes; install freetype as well])
      sage_spkg_install_freetype=yes
    else
      AC_MSG_RESULT([no])
      dnl freetype versions are libtool's ones, cf trac #30014
      PKG_CHECK_MODULES([FREETYPE], [freetype2 >= 16.1], [], [sage_spkg_install_freetype=yes])
    fi
    if test x$sage_spkg_install_freetype = xyes; then
      AC_SUBST(SAGE_FREETYPE_PREFIX, ['$SAGE_LOCAL'])
    else
      AC_SUBST(SAGE_FREETYPE_PREFIX, [''])
    fi
])



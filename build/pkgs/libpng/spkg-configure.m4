SAGE_SPKG_CONFIGURE([libpng], [
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_ZLIB])
    AC_MSG_CHECKING([Installing zlib? ])
    if test x$sage_spkg_install_zlib = xyes; then
      AC_MSG_RESULT([Yes. Install libpng as well.])
      sage_spkg_install_libpng=yes
    else
      AC_MSG_RESULT([No.])
      dnl First try checking for libpng with pkg-config
      PKG_CHECK_MODULES([LIBPNG], [libpng >= 1.2], [], [
        dnl Fallback to manually grubbing around for headers and libs
        AC_CHECK_HEADERS([png.h], [break], [sage_spkg_install_libpng=yes])
        AC_SEARCH_LIBS([png_get_io_ptr], [png], [], [sage_spkg_install_libpng=yes])
      ])
    fi
])

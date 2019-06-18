SAGE_SPKG_CONFIGURE([zlib], [
    AX_CHECK_ZLIB([ 
      PKG_CHECK_MODULES([LIBPNG], [libpng >= 1.2], [], [
         dnl inflateValidate is needed for Sage's libpng, newer than 1.2; this ensures
         dnl we have the minimum required for building zlib version
         AC_CHECK_LIB([z], [inflateValidate], [], [sage_spkg_install_zlib=yes])
      ])], [sage_spkg_install_zlib=yes])
])

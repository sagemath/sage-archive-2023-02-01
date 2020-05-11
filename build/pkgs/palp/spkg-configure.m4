SAGE_SPKG_CONFIGURE([palp], [
   AC_PATH_PROG([PALP], [poly.x])
   AS_IF([test -z "$ac_cv_path_PALP"], [sage_spkg_install_palp=yes])
])

SAGE_SPKG_CONFIGURE([eclib], [
    SAGE_SPKG_DEPCHECK([ntl pari flint], [
    dnl Trac #31443: use existing eclib only if the version reported by pkg-config is correct
    PKG_CHECK_MODULES([ECLIB],
                      [eclib = v20210317],
                      [sage_spkg_install_eclib=no],
                      [sage_spkg_install_eclib=yes])
    AC_PATH_PROG([MWRANK], [mwrank])
    AS_IF([test -z "$ac_cv_path_MWRANK"], [sage_spkg_install_eclib=yes])
  ])
])

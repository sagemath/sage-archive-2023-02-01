SAGE_SPKG_CONFIGURE([linbox], [
  SAGE_SPKG_DEPCHECK([fflas-ffpack flint fplll givaro gmp iml m4ri m4rie mpfr ntl], [
    PKG_CHECK_MODULES([LINBOX],
                      [linbox >= 1.6.3],
                      [],
                      [sage_spkg_install_linbox=yes])
  ])
])

SAGE_SPKG_CONFIGURE([gsl], [
    m4_pushdef([SAGE_GSL_MINVER],["2.5"])
    SAGE_SPKG_DEPCHECK([atlas openblas], [
      PKG_CHECK_MODULES([GSL], [gsl >= $SAGE_GSL_MINVER], [], [sage_spkg-install_gsl=yes])
    ])
    m4_popdef([SAGE_GSL_MINVER])
])

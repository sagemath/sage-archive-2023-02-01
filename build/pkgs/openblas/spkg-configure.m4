SAGE_SPKG_CONFIGURE([openblas], [
    PKG_CHECK_MODULES([OPENBLAS], [openblas >= 0.3.5 blas lapack], [], [sage_spkg_install_openblas=yes])
])

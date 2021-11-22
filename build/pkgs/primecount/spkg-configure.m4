SAGE_SPKG_CONFIGURE([primecount], [
    m4_pushdef([SAGE_PRIMECOUNT_MINVER],["7.1"])
    SAGE_SPKG_DEPCHECK([primesieve], [
      dnl Checking for primecount with pkg-config
      PKG_CHECK_MODULES([PRIMECOUNT], [primecount >= $SAGE_PRIMECOUNT_MINVER], [ ], [
          AC_CHECK_HEADER([primecount.h], [
           AC_SEARCH_LIBS([primecount_pi], [primecount], [
             dnl rely on primesieve being new enough
             ],
              [sage_spkg_install_primecount=yes])
          ], [sage_spkg_install_primecount=yes])
          sage_spkg_install_primecount=yes])
    ])
    m4_popdef([SAGE_PRIMECOUNT_MINVER])
])


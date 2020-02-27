SAGE_SPKG_CONFIGURE([fplll], [
  SAGE_SPKG_DEPCHECK([mpfr], [
    dnl If we're using the system mpfr, use pkgconfig to determine
    dnl if there's a usable system copy of fplll. Unless there's
    dnl a system that ships fplll without fplll.pc file, falling
    dnl back to a manual header/library search is pointless.
    PKG_CHECK_MODULES([FPLLL],
                      [fplll >= 5.3],
                      [],
                      [sage_spkg_install_fplll=yes])
  ],
  [ dnl If we're installing sage's mpfr, then we have to install
    dnl its fplll, too.
    sage_spkg_install_fplll=yes
  ])
])

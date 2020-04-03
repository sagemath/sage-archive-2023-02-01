SAGE_SPKG_CONFIGURE([ppl], [
  SAGE_SPKG_DEPCHECK([glpk gmp mpir], [
    # If our dependencies come from the system, then we can use
    # the system ppl, too. This macro works sort-of like the
    # PKG_CHECK_MODULES macro, defining e.g. PPL_CFLAGS when a
    # suitable version of PPL is detected. But notably, it doesn't
    # define PPL_LIBS.
    AM_PATH_PPL([1.2],
                [sage_spkg_install_ppl=no],
                [sage_spkg_install_ppl=yes])
  ],
  [ # Some of its dependencies are installed as SPKGs, so install the
    # ppl SPKG as well.
    sage_spkg_install_ppl=yes
  ])
])

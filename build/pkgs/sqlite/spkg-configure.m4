SAGE_SPKG_CONFIGURE([sqlite], [
  SAGE_SPKG_DEPCHECK([readline], [
    PKG_CHECK_MODULES([SQLITE],
                      [sqlite3 >= 3.29.0],
                      [],
                      [sage_spkg_install_sqlite=yes])
  ])
])

SAGE_SPKG_CONFIGURE([cddlib], [
    SAGE_SPKG_DEPCHECK([gmp mpir], [
        AC_PATH_PROG([CDDEXEC], [cddexec])
        AS_IF([test x$CDDEXEC = x], [
           AC_MSG_NOTICE([cddexec not found. Installing cddlib])
           sage_spkg_install_cddlib=yes
           ], [
           AC_PATH_PROG([CDDEXECGMP], [cddexec_gmp])
           AS_IF([test x$CDDEXECGMP = x], [
              AC_MSG_NOTICE([cddexec_gmp not found. Installing cddlib])
              sage_spkg_install_cddlib=yes])
           ])
    ])
])

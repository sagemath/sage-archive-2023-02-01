SAGE_SPKG_CONFIGURE([ecl], [
    SAGE_SPKG_DEPCHECK([gc gmp mpir], [
	        AC_PATH_PROG([ECL_CONFIG], [ecl-config])
	        AS_IF([test x$ECL_CONFIG = x], [
	           AC_MSG_NOTICE([ecl-config not found. Installing ecl])
	           sage_spkg_install_ecl=yes])
    ])
])

SAGE_SPKG_CONFIGURE([ecl], [
    SAGE_SPKG_DEPCHECK([gc gmp mpir], [
	        AC_PATH_PROG([ECL], [ecl])
	        AS_IF([test x$ECL = x], [
	           AC_MSG_NOTICE([ecl not found. Installing ecl])
	           sage_spkg_install_ecl=yes])
    ])
])

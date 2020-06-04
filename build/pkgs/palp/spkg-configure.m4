AC_DEFUN([SAGE_TEST_PALP_PROGS], [
    AC_PATH_PROG(PALP[$1], [[$1].x])
    AS_IF([test x$PALP[$1] = x], [sage_spkg_install_palp=yes])
    m4_foreach([suff], [4, 5, 6, 11], [
         AC_PATH_PROG(PALP[$1]suff, [[$1][-]suff[d.x]])
         AS_IF([test x$PALP[$1]suff = x], [sage_spkg_install_palp=yes])
    ])
])

SAGE_SPKG_CONFIGURE([palp], [
   dnl m4_foreach([palpprog], [[poly], [class], [nef], [cws]], [
     SAGE_TEST_PALP_PROGS(poly)
   dnl ])
])

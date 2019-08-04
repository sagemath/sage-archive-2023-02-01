SAGE_SPKG_CONFIGURE([pari], [
    dnl See gp_version below on how the version is computed from MAJV.MINV.PATCHV
    m4_pushdef([SAGE_PARI_MINVER],["133889"])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_GMP])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_READLINE])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI_ELLDATA])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI_GALDATA])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI_GALPOL])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI_SEADATA])
    AC_REQUIRE([SAGE_SPKG_CONFIGURE_PARI_SEADATA_SMALL])
    AC_MSG_CHECKING([installing gmp/mpir? ])
    if test x$sage_spkg_install_mpir = xyes -o x$sage_spkg_install_gmp = xyes; then
        AC_MSG_RESULT([yes; install pari as well])
        sage_spkg_install_pari=yes
    else
      AC_MSG_RESULT([no])
      AC_MSG_CHECKING([installing readline or PARI/GP packages? ])
      if test x$sage_spkg_install_readline = xyes -o x$sage_spkg_install_pari_galdata = xyes -o x$sage_spkg_install_pari_galpol = xyes \
         -o x$sage_spkg_install_pari_seadata_small = xyes -o x$sage_spkg_install_pari_seadata = xyes -o x$sage_spkg_install_pari_elldata = xyes; then
          AC_MSG_RESULT([yes; install pari as well])
          sage_spkg_install_pari=yes
      else
        AC_MSG_RESULT([no])
        AC_CHECK_HEADER([pari/pari.h], [], [sage_spkg_install_pari=yes])
        dnl matpermanent appears in pari 2.11
        AC_SEARCH_LIBS([matpermanent], [pari], [
          if test x$GP != x; then
              AC_MSG_CHECKING([getting GP's version ])
              gp_version=`echo "v=version(); v[[1]]<<16 + v[[2]]<<8 + v[[3]]" | $GP -qf`
              AC_MSG_RESULT([$gp_version])
              AC_MSG_CHECKING([comparing GP and libpari versions])
              AC_LANG_PUSH(C)
              AC_RUN_IFELSE([AC_LANG_PROGRAM([
                  [#include <pari/pari.h>]],
	         [[long vers;]
                  [pari_init(5000000, 2);]
                  [vers=paricfg_version_code;]
                  [pari_close()];
                  [return vers!=$gp_version;]])],
                [AC_MSG_RESULT([libpari's and GP's versions match. Good])],
	        [AC_MSG_RESULT([libpari's version does not match GP's version. Not good])
		         sage_spkg_install_pari=yes])
              AC_MSG_CHECKING([is GP's version good enough? ])
              AX_COMPARE_VERSION([$gp_version], [ge], [$SAGE_PARI_MINVER], [
                  AC_MSG_RESULT([yes])
                  AC_MSG_CHECKING([getting GP's datadir])
                  gp_datadir=`echo "default(datadir)" | $GP -qf`
                  AC_MSG_RESULT([$gp_datadir])
                  AC_MSG_CHECKING([comparing GP's and libpari's datadirs])
                  AC_RUN_IFELSE([AC_LANG_PROGRAM([
                      [#include <pari/pari.h>]],
	             [[int t;]
                      [pari_init(5000000, 2);]
                      [t=strcmp(paricfg_datadir,$gp_datadir);]
                      [pari_close()];
                      [return t;]])],
                    [AC_MSG_RESULT([libpari's and GP's datadirs match. Good])],
	            [AC_MSG_RESULT([libpari's datadir does not match GP's datadir. Not good])
		         sage_spkg_install_pari=yes])
                 ], [
                  AC_MSG_RESULT([no])
                  sage_spkg_install_pari=yes])
              AC_LANG_POP()
          else
              sage_spkg_install_pari=yes
          fi
        ], [sage_spkg_install_pari=yes])
      fi
    fi
    m4_popdef([SAGE_PARI_MINVER])
], [], [], [
    if test x$sage_spkg_install_pari = xyes; then
        AC_SUBST(SAGE_PARI_PREFIX, ['$SAGE_LOCAL'])
        AC_MSG_RESULT([using Sage's pari SPKG])
    else
       AC_SUBST(SAGE_PARI_PREFIX, [''])
       AC_MSG_RESULT([using pari/gp from the system])
    fi
])

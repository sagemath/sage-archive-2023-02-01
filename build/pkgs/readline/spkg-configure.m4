SAGE_SPKG_CONFIGURE([readline], [
  AC_REQUIRE([SAGE_SPKG_CONFIGURE_NCURSES])
  AC_MSG_CHECKING([Installing ncurses? ])
  if test x$sage_spkg_install_ncurses = xyes; then
     AC_MSG_RESULT([Yes. Install readline as well.])
     sage_spkg_install_readline=yes
  else
    AC_MSG_RESULT([No.])
  dnl First try checking for readline with pkg-config
    PKG_CHECK_MODULES([READLINE], [readline >= 6.0], [],
      [AC_CHECK_HEADERS([readline/readline.h],
  dnl rl_replace_line is not present in macos's readline
        [AC_SEARCH_LIBS([rl_replace_line], [readline], [break],
                        [sage_spkg_install_readline=yes])], 
        [sage_spkg_install_readline=yes])],
      [sage_spkg_install_readline=yes])
  fi
])

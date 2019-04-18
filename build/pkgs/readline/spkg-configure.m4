SAGE_SPKG_CONFIGURE([readline], [
    dnl First try checking for readline with pkg-config
    PKG_CHECK_MODULES([READLINE], [readline >= 6.0], [],
      [AC_CHECK_HEADERS([readline/readline.h],
        [dnl AC_SEARCH_LIBS([wresize], [readline], [break],
         dnl              [sage_spkg_install_readline=yes])
        ], 
        [sage_spkg_install_readline=yes])],
      [sage_spkg_install_readline=yes])
])

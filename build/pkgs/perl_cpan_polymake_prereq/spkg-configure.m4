SAGE_SPKG_CONFIGURE(
    [perl_cpan_polymake_prereq], [
    m4_pushdef([MODULES], m4_include(build/pkgs/perl_cpan_polymake_prereq/distros/cpan.txt))
    AX_PROG_PERL_MODULES(MODULES,
      [],
      [sage_spkg_install_perl_cpan_polymake_prereq=yes
       AS_CASE([SAGE_OPTIONAL_INSTALLED_PACKAGES],
         [*polymake*], [
           dnl We do not exit with AC_MSG_ERROR so that the system package notice at
           dnl the end of configure will be displayed.
           AC_MSG_WARN([Optional package polymake needs a working installation of Perl and modules ]MODULES)
         ]
       )
    ])
])

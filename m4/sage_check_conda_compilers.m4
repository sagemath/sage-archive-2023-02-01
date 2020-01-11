AC_DEFUN([SAGE_CHECK_CONDA_COMPILERS], [
    AC_MSG_CHECKING([whether a conda environment is active])
    AS_IF([test "x$CONDA_PREFIX" != x], [have_conda_active=yes], [have_conda_active=no])
    AC_MSG_RESULT($have_conda_active)
    AS_IF([test $have_conda_active = yes], [
        dnl A conda environment is active.
        dnl #27699: Conda compiler packages must be installed
        need_pkgs="c-compiler cxx-compiler fortran-compiler"
        for pkg in $need_pkgs; do
            AC_MSG_CHECKING([whether conda package "$pkg" is installed in the current conda environment])
            AS_IF([$CONDA_EXE list -c -f $pkg | grep -q $pkg], [have_pkg=yes], [have_pkg=no])
            AC_MSG_RESULT($have_pkg)
            AS_IF([test $have_pkg = no],
                  [AC_MSG_ERROR([A conda environment ($CONDA_DEFAULT_ENV) is active, but it is missing
the following conda packages required for building Sage:
    $need_pkgs
For building Sage, either:
- activate a conda environment that has these packages, using:
    conda activate ENVIRONMENT
- or install these conda packages, using
    conda install $need_pkgs
- or deactivate conda by
    conda deactivate
  (this command may need to be repeated).])
                  ])
        done
    ])
])

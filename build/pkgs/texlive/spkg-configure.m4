m4_pushdef([required_latex_packages], [fontspec,xunicode,xltxtra,amssymb,amsfonts,graphicx,mathrsfs,textcomp,tikz,tikz-qtree,iftex,tkz-berge,tkz-graph,xy,babel,subfigure,hyperref,hypcap,xr])
SAGE_SPKG_CONFIGURE([texlive], [
    sage_spkg_install_texlive=no
    dnl pdflatex --version
    m4_foreach([latex_package], [required_latex_packages], [
        AC_MSG_CHECKING([for latex package ]latex_package)
        AS_IF([kpsewhich ]latex_package[.sty >& AS_MESSAGE_LOG_FD 2>&1], [
            AC_MSG_RESULT([yes])
            sage_spkg_install_texlive=yes
        ], [
            AC_MSG_RESULT([no])
        ])
    ])
])
m4_popdef([required_latex_packages])

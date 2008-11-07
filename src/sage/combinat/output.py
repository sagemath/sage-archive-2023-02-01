from sage.misc.latex import latex

def tex_from_array(a):
    """
    EXAMPLES:
        sage: from sage.combinat.output import tex_from_array
        sage: print tex_from_array([[1,2,3],[4,5]])
        {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{ccc}
        \cline{1-1}\cline{2-2}\cline{3-3}
        \lr{1}&\lr{2}&\lr{3}\\
        \cline{1-1}\cline{2-2}\cline{3-3}
        \lr{4}&\lr{5}\\
        \cline{1-1}\cline{2-2}
        \end{array}$}
        }
    """
    rows = len(a)
    cols = len(a[0])
    s = ""
    s += "{\\def\\lr#1{\\multicolumn{1}{|@{\\hspace{.6ex}}c@{\\hspace{.6ex}}|}{\\raisebox{-.3ex}{$#1$}}}\n"
    s += "\\raisebox{-.6ex}{$"
    s += "\\begin{array}[b]{"+"c"*cols+"}\n"
    for j in range(cols):
        if a[0][j] is not None:
            s += "\\cline{" + str(j+1) + "-" + str(j+1) + "}"
    s += "\n"
    for i in range(rows):
        cols = len(a[i])
        if len(a[i])>0 and a[i][0] is not None:
            for j in range(cols):
                if a[i][j] is not None:
                    s += "\\lr{" + str(a[i][j]) + "}"
                    if j < cols-1:
                        s += "&"
            s += "\\\\\n"
            for j in range(cols):
                if a[i][j] is not None:
                    s += "\\cline{" + str(j+1) + "-" + str(j+1) + "}"
            s += "\n"
    s += "\\end{array}$}\n}"
    return s



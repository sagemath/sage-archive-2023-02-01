from sage.misc.latex import latex

def tex_from_array(a):
    """
    EXAMPLES:
        sage: from sage.combinat.output import tex_from_array
        sage: print tex_from_array([[1,2],[3,4]])
        {\def\lr#1#2#3{\multicolumn{1}{#1@{\hspace{.6ex}}c@{\hspace{.6ex}}#2}{\raisebox{-.3ex}{$#3$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{cc}
        \cline{1-1}\cline{2-2}%
        \lr{|}{|}{1}&\lr{|}{|}{2}\\ %
        \cline{1-1}\cline{2-2}%
        \lr{|}{|}{3}&\lr{|}{|}{4}\\ %
        \cline{1-1}\cline{2-2}%
        \end{array}$}
        }
    """
    rows = len(a)
    cols = len(a[0])
    s = ""
    s += "{\\def\\lr#1#2#3{\\multicolumn{1}{#1@{\\hspace{.6ex}}c@{\\hspace{.6ex}}#2}{\\raisebox{-.3ex}{$#3$}}}\n"
    s += "\\raisebox{-.6ex}{$"
    s += "\\begin{array}[b]{"+"c"*cols+"}\n"
    for i in range(rows):
        for j in range(cols):
            s += "\\cline{" + str(j+1) + "-" + str(j+1) + "}"
            if a[i][j] is None:
                break

        s += "%\n"
        for j in range(cols):
            if a[i][j] is not None:
                s += "\\lr"
                s += "{|}"
                s += "{|}"
                s += "{" + str(a[i][j]) + "}"
            if j != cols-1:
                s += "&"
        s += "\\\\ %\n"

    for j in range(cols):
        if a[i][j] is None:
            break
        s += "\\cline{" + str(j+1) + "-" + str(j+1) + "}"

    s += "%\n"
    s += "\\end{array}$}\n}"
    return s



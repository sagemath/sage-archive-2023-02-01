def version():
    """
    Return the version of SAGE.

    INPUT:
       nothing
    OUTPUT:
       str

    EXAMPLES:
       sage.: print version()
       SAGE Version x.y.z, Build Date: xxxx-xx-xx-xxxx
    """
    import sage.version
    return 'SAGE Version %s, Build Date: %s'%(sage.version.version, sage.version.date)

def banner_text():
    bars = "-"*56
    s = bars + '\n'
    s += "| %-52s |\n"%version()
    s += "| %-52s |\n"%'Distributed under the GNU General Public License V2.'
    s += "| %-52s |\n"%'For help type <anything>? or <anything>??.'
    s += bars
    return s


def banner():
    print banner_text() + '\n'



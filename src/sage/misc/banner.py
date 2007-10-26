def version():
    """
    Return the version of SAGE.

    INPUT:
       nothing
    OUTPUT:
       str

    EXAMPLES:
       sage: version()
       'SAGE Version ..., Release Date: ...'
    """
    import sage.version
    return 'SAGE Version %s, Release Date: %s'%(sage.version.version, sage.version.date)

def banner_text():
    bars = "-"*70
    s = bars
    s += "\n| %-66s |\n"%version()
    s += "| %-66s |\n"%'Type notebook() for the GUI, and license() for information.'
    #s += "| %-66s |\n"%'Distributed under the GNU General Public License V2.'
    s += bars
    return s


def banner():
    print banner_text()



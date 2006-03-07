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

def banner():
    bars = "-"*56
    print bars
    print "| %-52s |"%version()
    print "| %-52s |"%'Distributed under the GNU General Public License V2.'
    print "| %-52s |"%'For help type <anything>? or <anything>??.'
    print bars
    print ""


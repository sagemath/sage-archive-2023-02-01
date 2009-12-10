# This file is part of the OLD Sage notebook and is NOT actively developed,
# maintained, or supported.  As of Sage v4.1.2, all notebook development has
# moved to the separate Sage Notebook project:
#
# http://nb.sagemath.org/
#
# The new notebook is installed in Sage as an spkg (e.g., sagenb-0.3.spkg).
#
# Please visit the project's home page for more information, including directions on
# obtaining the latest source code.  For notebook-related development and support,
# please consult the sage-notebook discussion group:
#
# http://groups.google.com/group/sage-notebook

import interact

def javascript(s):
    print '<html><script>%s</script></html>'%s

def cell_id():
    return interact.SAGE_CELL_ID

def draggable(option=''):
    """
    INPUT:
        option -- string
            "" (default) -- Creates new draggables on the current cell.
            "enable" -- Enable the draggable functionality.
            "disable" -- Temporarily disable the draggable functionality.
    """
    s = '$("#cell_outer_%s").draggable("%s")'%(cell_id(), option)
    javascript(s)

def resizable():
    s = '$("#cell_outer_%s").resizable()'%cell_id()
    javascript(s)


import manipulate

def javascript(s):
    print '<html><script>%s</script></html>'%s

def cell_id():
    return manipulate.SAGE_CELL_ID

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

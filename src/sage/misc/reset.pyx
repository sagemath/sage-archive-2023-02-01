import sys

def reset():
    """
    Delete all user defined variables, reset all globals variables
    back to their default state, and reset all interfaces to other
    computer algebra systems.
    """
    G = globals()  # this is the reason the code must be in SageX.
    T = type(sys)
    for k in G.keys():
        if k[0] != '_' and type(k) != T:
            del G[k]
    restore()
    reset_interfaces()


def restore():
    """
    Restore all predefined global variables to their default values.
    """
    G = globals()  # this is the reason the code must be in SageX.
    import sage.all
    for k, v in sage.all.__dict__.iteritems():
        G[k] = v


def reset_interfaces():
    from sage.interfaces.quit import expect_quitall
    expect_quitall()

##     import sys
##     M = sys.modules
##     for k in M.keys():
##         if 'sage.interfaces' in k:
##             if not M[k] is None:
##                 reload(M[k])

##     import sage.interfaces.all
##     G = globals()
##     for k, v in sage.interfaces.all.__dict__.iteritems():
##         G[k] = v



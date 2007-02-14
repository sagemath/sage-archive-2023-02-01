in_dsage_mode = False

def start_dsage_console():
    if in_dsage_mode:
        return
    import dsage_console
    shell = dsage_console.IPShellTwisted(
        argv=[],
        user_ns=globals())
    print 'Starting Distributed SAGE console...'
    global in_dsage_mode
    in_dsage_mode = True
    shell.mainloop()


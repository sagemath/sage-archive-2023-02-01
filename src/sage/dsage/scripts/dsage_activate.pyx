in_console_mode = False

def start_dsage_console():
    if in_console_mode:
        return
    import dsage_console
    shell = dsage_console.IPShellTwisted(
        argv=[],
        user_ns=globals())
    print 'Starting Distributed SAGE...'
    global in_console_mode
    in_console_mode = True
    shell.mainloop()


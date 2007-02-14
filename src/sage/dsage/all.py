from dsage import dsage
from sage.dsage.interface.dsage_interface import BlockingDSage

def DSage(server=None, port=8081, username=None,
          pubkey_file=None, privkey_file=None):
    """
    This object represents a connection to the distributed SAGE server.

    INPUT:
        server -- hostname of the DSAGE server (str)
        port -- port of the server (int)
        username -- username stored on the server (str)
        pubkey_file -- file that stores the users public key
        privkey_file -- file that stores the users private key
    """
    import sage.dsage.scripts.dsage_activate as activate
    if not activate.in_dsage_mode:
        raise ValueError, "You must first turn on distributed SAGE using dsage.console()."

    from sage.dsage.interface.dsage_interface import DSage
    return DSage(server = server, port = port, username = username,
                 pubkey_file = pubkey_file, privkey_file = privkey_file)


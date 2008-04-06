import socket
import gnutls
import gnutls.connection
import gnutls.errors

class GnuTLSSocketSSL:
    """
    This class provides a bare-bones replacement to socket.ssl in the case
    that openssl is not installed/detected when Python is built.

    socket.ssl MUST be set before urllib is imported.
    """
    def __init__(self, sock, key_file=None, cert_file=None):
        self.creds = gnutls.connection.X509Credentials()
        self.session = gnutls.connection.ClientSession(sock, self.creds)
        self.session.handshake()

    def send(self, *args, **kwds):
        self.session.send(*args, **kwds)

    write = send

    def recv(self, *args, **kwds):
        try:
            s = self.session.recv(*args, **kwds)
        except gnutls.errors.GNUTLSError:
            raise socket.sslerror(socket.SSL_ERROR_EOF)
        return s

    read = recv

    def close(self):
        self.session.bye()
        self.session.shutdown()
        self.session.close()

def require_SSL():
    """
    If ssl does not already exist in the socket module, supply our gnutls
    version.
    """
    if not hasattr(socket, "ssl"):
        socket.ssl = GnuTLSSocketSSL

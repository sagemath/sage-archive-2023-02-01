import sage.dev.trac_interface

class DoctestTracInterface(sage.dev.trac_interface.TracInterface):
    def __init__(self, config, UI, server):
        sage.dev.trac_interface.TracInterface.__init__(self, config, UI)

        self._server = server

    @property
    def _anonymous_server_proxy(self):
        from server_proxy import DoctestServerProxy
        return DoctestServerProxy(self._server)

    @property
    def _authenticated_server_proxy(self):
        from server_proxy import AuthenticatedDoctestServerProxy
        return AuthenticatedDoctestServerProxy(self._server, self._username, self._password)
